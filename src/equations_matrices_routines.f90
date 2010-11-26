!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all equations matrix and rhs routines.
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

!> This module handles all equations matrix and rhs routines.
MODULE EQUATIONS_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES
  USE LINKEDLIST_ROUTINES
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes EQUATIONS_MATRICES_ROUTINES::EquationsMatrixStructureTypes
  !> \brief Equations matrices structure (sparsity) types
  !> \see EQUATIONS_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_DIAGONAL_STRUCTURE=3 !<Diagonal matrix structure. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
  !>@}


  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_LumpingTypes EQUATIONS_MATRICES_ROUTINES::LumpingTypes
  !> \brief Equations matrix lumping types
  !> \see EQUATIONS_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_UNLUMPED=1 !<The matrix is not lumped \see EQUATIONS_MATRICES_ROUTINES_LumpingTypes,EQUATIONS_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_LUMPED=2 !<The matrix is "mass" lumped \see EQUATIONS_MATRICES_ROUTINES_LumpingTypes,EQUATIONS_MATRICES_ROUTINES
 !>@}
 
  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_EquationsMatricesSparsityTypes EQUATIONS_MATRICES_ROUTINES::EquationsMatricesSparsityTypes
  !> \brief Equations matrices sparsity types
  !> \see EQUATIONS_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_SPARSE_MATRICES=1 !<Use sparse equations matrices \see EQUATIONS_MATRICES_ROUTINES_EquationsMatricesSparsityTypes,EQUATIONS_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_FULL_MATRICES=2 !<Use fully populated equation matrices \see EQUATIONS_MATRICES_ROUTINES_EquationsMatricesSparsityTypes,EQUATIONS_MATRICES_ROUTINES
  !>@}

  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_SelectMatricesTypes EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes
  !> \brief The types of selection available for the equations matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_ALL=1 !<Select all the equations matrices and vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic equations matrices and vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
   INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_LINEAR_ONLY=3 !<Select only the linear equations matrices and vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear equations matrices and vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian equations matrix \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual equations vector \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RHS_ONLY=7 !<Select only the RHS equations vector \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_SOURCE_ONLY=8 !<Select only the RHS equations vector \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY=9 !<Select only the RHS and residual equations vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RHS_SOURCE_ONLY=10 !<Assemble only the RHS and source equations vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY=11 !<Assemble only the residual and source equations vectors\see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES 
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_VECTORS_ONLY=12 !<Assemble only the equations vectors \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATIONS_MATRIX_NO_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE

  PUBLIC EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED

  PUBLIC EQUATIONS_MATRICES_SPARSE_MATRICES,EQUATIONS_MATRICES_FULL_MATRICES

  PUBLIC EQUATIONS_MATRICES_CREATE_FINISH,EQUATIONS_MATRICES_CREATE_START,EQUATIONS_MATRICES_DESTROY

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC EQUATIONS_MATRICES_ELEMENT_ADD,EQUATIONS_MATRICES_ELEMENT_CALCULATE,EQUATIONS_MATRICES_ELEMENT_INITIALISE, &
    & EQUATIONS_MATRICES_ELEMENT_FINALISE,EQUATIONS_MATRICES_VALUES_INITIALISE

  PUBLIC EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE,EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE, &
    & EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE,EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP

  PUBLIC EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE,EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE, &
    & EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE,EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP

  PUBLIC EQUATIONS_MATRICES_ALL,EQUATIONS_MATRICES_LINEAR_ONLY,EQUATIONS_MATRICES_NONLINEAR_ONLY,EQUATIONS_MATRICES_JACOBIAN_ONLY, &
    & EQUATIONS_MATRICES_RESIDUAL_ONLY,EQUATIONS_MATRICES_RHS_ONLY,EQUATIONS_MATRICES_SOURCE_ONLY, &
    & EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY,EQUATIONS_MATRICES_RHS_SOURCE_ONLY,EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY, &
    & EQUATIONS_MATRICES_VECTORS_ONLY
  
  PUBLIC EQUATIONS_MATRICES_OUTPUT

  PUBLIC EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET,EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET, &
    & EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET

  PUBLIC EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET,EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET

  PUBLIC EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET,EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalise the equations Jacobian and deallocate all memory
  SUBROUTINE EQUATIONS_JACOBIAN_FINALISE(EQUATIONS_JACOBIAN,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: EQUATIONS_JACOBIAN !<A pointer to the equations Jacobian to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("EQUATIONS_JACOBIAN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_JACOBIAN)) THEN
      IF(ASSOCIATED(EQUATIONS_JACOBIAN%JACOBIAN)) CALL DISTRIBUTED_MATRIX_DESTROY(EQUATIONS_JACOBIAN%JACOBIAN,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_JACOBIAN%ELEMENT_JACOBIAN,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(EQUATIONS_JACOBIAN%ELEMENT_RESIDUAL,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_JACOBIAN_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_JACOBIAN_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_JACOBIAN_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_JACOBIAN_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the equations Jacobian.
  SUBROUTINE EQUATIONS_JACOBIAN_INITIALISE(NONLINEAR_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer to the equations matrices nonlinear matrices to initialise the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_JACOBIAN_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
      EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
      IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
        EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            IF(ASSOCIATED(NONLINEAR_MATRICES%JACOBIAN)) THEN
              CALL FLAG_ERROR("Nonlinear matrices Jacobian is already associated.",ERR,ERROR,*998)
            ELSE
              ALLOCATE(NONLINEAR_MATRICES%JACOBIAN,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations Jacobian.",ERR,ERROR,*999)
              NONLINEAR_MATRICES%JACOBIAN%NONLINEAR_MATRICES=>NONLINEAR_MATRICES
              NONLINEAR_MATRICES%JACOBIAN%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
              NONLINEAR_MATRICES%JACOBIAN%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
              NONLINEAR_MATRICES%JACOBIAN%NUMBER_OF_COLUMNS=NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%NUMBER_OF_COLUMNS
              NONLINEAR_MATRICES%JACOBIAN%UPDATE_JACOBIAN=.TRUE.
              NONLINEAR_MATRICES%JACOBIAN%FIRST_ASSEMBLY=.TRUE.
              NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%JACOBIAN=>NONLINEAR_MATRICES%JACOBIAN
              NULLIFY(NONLINEAR_MATRICES%JACOBIAN%JACOBIAN)
              CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(NONLINEAR_MATRICES%JACOBIAN%ELEMENT_JACOBIAN,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(NONLINEAR_MATRICES%JACOBIAN%ELEMENT_RESIDUAL,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Nonlinear matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_JACOBIAN_INITIALISE")
    RETURN
999 CALL EQUATIONS_JACOBIAN_FINALISE(NONLINEAR_MATRICES%JACOBIAN,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_JACOBIAN_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_JACOBIAN_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_JACOBIAN_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the equations matrices and RHS for the the equations
  SUBROUTINE EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<The pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    type(LinkedList),pointer :: list(:) 
    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)

    CALL ENTERS("EQUATIONS_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have already been finished.",ERR,ERROR,*998)
      ELSE
        EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          ROW_DOMAIN_MAP=>EQUATIONS_MAPPING%ROW_DOFS_MAPPING
          IF(ASSOCIATED(ROW_DOMAIN_MAP)) THEN
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
              !Dynamic matrices
              DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
              IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
                !Now create the individual dynamic equations matrices
                DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
                  EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
                  IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                    COLUMN_DOMAIN_MAP=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
                      !Create the distributed equations matrix
                      CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,EQUATIONS_MATRICES% &
                           & DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(EQUATIONS_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(EQUATIONS_MATRIX%MATRIX,EQUATIONS_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                      !Calculate and set the matrix structure/sparsity pattern
                      IF(EQUATIONS_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                        & EQUATIONS_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                        CALL EQUATIONS_MATRIX_STRUCTURE_CALCULATE(EQUATIONS_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                          & list,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_LINKLIST_SET(EQUATIONS_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(EQUATIONS_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(EQUATIONS_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                        IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                      ENDIF
                      CALL DISTRIBUTED_MATRIX_CREATE_FINISH(EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Column domain map for dynamic matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                      & " is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx                
              ELSE
                CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)                
              ENDIF
            ENDIF
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            IF(ASSOCIATED(LINEAR_MATRICES)) THEN
              !Linear matrices
              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
              IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                !Now create the individual linear equations matrices
                DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                  EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
                  IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                    COLUMN_DOMAIN_MAP=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
                      !Create the distributed equations matrix
                      CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,EQUATIONS_MATRICES% &
                           & LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(EQUATIONS_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(EQUATIONS_MATRIX%MATRIX,EQUATIONS_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                      !Calculate and set the matrix structure/sparsity pattern
                      IF(EQUATIONS_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                        & EQUATIONS_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                        CALL EQUATIONS_MATRIX_STRUCTURE_CALCULATE(EQUATIONS_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                          & list,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_LINKLIST_SET(EQUATIONS_MATRIX%MATRIX,LIST,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(EQUATIONS_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(EQUATIONS_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                        IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                      ENDIF
                      CALL DISTRIBUTED_MATRIX_CREATE_FINISH(EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Column domain map for linear matrix number "// &
                        & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                      & " is not associated."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
              ELSE
                CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)                
              ENDIF
            ENDIF
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
              !Nonlinear matrices
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                !Set up the Jacobian matrix
                JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
                IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                  COLUMN_DOMAIN_MAP=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%COLUMN_DOFS_MAPPING
                  IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
!!TODO: Set the distributed matrix not to allocate the data if the Jacobian is not calculated.
                    !Create the distributed Jacobian matrix
                    CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,EQUATIONS_MATRICES%NONLINEAR_MATRICES% &
                         & JACOBIAN%JACOBIAN,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(JACOBIAN_MATRIX%JACOBIAN,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(JACOBIAN_MATRIX%JACOBIAN,JACOBIAN_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                    !Calculate and set the matrix structure/sparsity pattern
                    IF(JACOBIAN_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                      & JACOBIAN_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                      CALL JACOBIAN_MATRIX_STRUCTURE_CALCULATE(JACOBIAN_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                        & ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(JACOBIAN_MATRIX%JACOBIAN,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(JACOBIAN_MATRIX%JACOBIAN,ROW_INDICES,COLUMN_INDICES, &
                        & ERR,ERROR,*999)
                      IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                      IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                    ENDIF
                    CALL DISTRIBUTED_MATRIX_CREATE_FINISH(JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Column domain map is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
                ENDIF
                !Set up the residual vector                
                CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,EQUATIONS_MATRICES%NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NONLINEAR_MATRICES%RESIDUAL,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_CREATE_FINISH(NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDIF
            RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              !Set up the equations RHS vector          
              CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,EQUATIONS_MATRICES%RHS_VECTOR%VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(RHS_VECTOR%VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(RHS_VECTOR%VECTOR,ERR,ERROR,*999)
            ENDIF
            SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
            IF(ASSOCIATED(SOURCE_VECTOR)) THEN
              !Set up the equations source vector          
              CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,EQUATIONS_MATRICES%SOURCE_VECTOR%VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOURCE_VECTOR%VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOURCE_VECTOR%VECTOR,ERR,ERROR,*999)
            ENDIF
            !Finish up
            EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Row domain map is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of the equations matrices and rhs for the the equations
  SUBROUTINE EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<The pointer to the equations to create the equations matrices for
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<On return, a pointer to the equations matrices being created.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR    

    CALL ENTERS("EQUATIONS_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN      
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
          CALL FLAG_ERROR("Equations matrices is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(EQUATIONS_MATRICES)
          !Initialise the equations matrices
          CALL EQUATIONS_MATRICES_INITIALISE(EQUATIONS,ERR,ERROR,*999)
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_CREATE_START")
    RETURN
999 CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS%EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the equations matrices
  SUBROUTINE EQUATIONS_MATRICES_DESTROY(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer the equations matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("EQUATIONS_MATRICES_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_DESTROY",ERR,ERROR)    
    CALL EXITS("EQUATIONS_MATRICES_DESTROY")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MATRICES_DESTROY

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the element matrix. Old CMISS name MELGE.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(ELEMENT_MATRIX,UPDATE_MATRIX,ROW_ELEMENT_NUMBER,COLUMN_ELEMENT_NUMBER, &
    & ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix to calculate
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if the element matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: ROW_ELEMENT_NUMBER !<The row element number to calculate
    INTEGER(INTG), INTENT(IN) :: COLUMN_ELEMENT_NUMBER !<The column element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLS_FIELD_VARIABLE !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,global_ny,local_ny,node,node_idx
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      IF(ASSOCIATED(COLS_FIELD_VARIABLE)) THEN
        ELEMENT_MATRIX%NUMBER_OF_ROWS=0
        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
        IF(UPDATE_MATRIX) THEN
          IF(ASSOCIATED(ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE)) THEN
            !Row and columns variable is the same.
            DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              IF(ROW_ELEMENT_NUMBER>=1.AND.ROW_ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                  ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                  ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                  ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                  ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP(ROW_ELEMENT_NUMBER)
                  global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                  ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                  ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                  ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                  ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ROW_ELEMENT_NUMBER)%BASIS
                  DO node_idx=1,BASIS%NUMBER_OF_NODES
                    node=ELEMENTS_TOPOLOGY%ELEMENTS(ROW_ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                    DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                      derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ROW_ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                      local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node)
                      global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO !derivative_idx
                  ENDDO !node_idx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The interpolation type of "// &
                    & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of rows field variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                END SELECT
              ELSE
                LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ROW_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                  & " of rows field variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & ". The element number must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !component_idx
          ELSE
            !Row and column variables are different
            !Row mapping
            DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              IF(ROW_ELEMENT_NUMBER>=1.AND.ROW_ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                  ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP(ROW_ELEMENT_NUMBER)
                  ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                  ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ROW_ELEMENT_NUMBER)%BASIS
                  DO node_idx=1,BASIS%NUMBER_OF_NODES
                    node=ELEMENTS_TOPOLOGY%ELEMENTS(ROW_ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                    DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                      derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ROW_ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                      local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node)
                      ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                      ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                    ENDDO !derivative_idx
                  ENDDO !node_idx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The interpolation type of "// &
                    & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of rows field variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                END SELECT
              ELSE
                LOCAL_ERROR="Row element number "//TRIM(NUMBER_TO_VSTRING(ROW_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                  & " of rows field variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & ". The element number must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !component_idx
            !Column mapping
            DO component_idx=1,COLS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              IF(COLUMN_ELEMENT_NUMBER>=1.AND.COLUMN_ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                SELECT CASE(COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                  ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                  ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                    & ELEMENT_PARAM2DOF_MAP(COLUMN_ELEMENT_NUMBER)
                  global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                  ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                  ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(COLUMN_ELEMENT_NUMBER)%BASIS
                  DO node_idx=1,BASIS%NUMBER_OF_NODES
                    node=ELEMENTS_TOPOLOGY%ELEMENTS(COLUMN_ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                    DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                      derivative=ELEMENTS_TOPOLOGY%ELEMENTS(COLUMN_ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                      local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node)
                      global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO !derivative_idx
                  ENDDO !node_idx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The interpolation type of "// &
                    & TRIM(NUMBER_TO_VSTRING(COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of column field variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(COLS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
                END SELECT
              ELSE
                LOCAL_ERROR="Column element number "//TRIM(NUMBER_TO_VSTRING(COLUMN_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                  & " of column field variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(COLS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & ". The element number must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !component_idx
          ENDIF
          ELEMENT_MATRIX%MATRIX=0.0_DP
        ENDIF
      ELSE
        CALL FLAG_ERROR("Columns field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise an element matrix and deallocate all memory
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE):: ELEMENT_MATRIX !<The element matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT_MATRIX%ROW_DOFS)) DEALLOCATE(ELEMENT_MATRIX%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%COLUMN_DOFS)) DEALLOCATE(ELEMENT_MATRIX%COLUMN_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%MATRIX)) DEALLOCATE(ELEMENT_MATRIX%MATRIX)
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element matrix.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !The element matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE",ERR,ERROR,*999)

    ELEMENT_MATRIX%EQUATIONS_MATRIX_NUMBER=0
    ELEMENT_MATRIX%NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
       
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets up the element matrix for the row and column field variables.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(ELEMENT_MATRIX,ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLS_FIELD_VARIABLE !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP",ERR,ERROR,*998)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      IF(ASSOCIATED(COLS_FIELD_VARIABLE)) THEN
        ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=ROWS_FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS* &
          & ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=COLS_FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS* &
          & COLS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        IF(ALLOCATED(ELEMENT_MATRIX%ROW_DOFS)) THEN
          CALL FLAG_ERROR("Element matrix row dofs already allocated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix row dofs.",ERR,ERROR,*999)
        ENDIF
        IF(ALLOCATED(ELEMENT_MATRIX%COLUMN_DOFS)) THEN
          CALL FLAG_ERROR("Element matrix column dofs already allocated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix column dofs.",ERR,ERROR,*999)
        ENDIF
        IF(ALLOCATED(ELEMENT_MATRIX%MATRIX)) THEN
          CALL FLAG_ERROR("Element matrix already allocated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(ELEMENT_MATRIX%MATRIX(ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS,ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Columns field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP")
    RETURN
999 CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(ELEMENT_MATRIX,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations rhs of the element rhs vector. Old CMISS name MELGE.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(ELEMENT_VECTOR,UPDATE_VECTOR,ELEMENT_NUMBER,ROWS_FIELD_VARIABLE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element vector to calculate.
    LOGICAL :: UPDATE_VECTOR !<Is .TRUE. if the element vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,local_ny,node,node_idx
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      !Calculate the rows for the element vector
      ELEMENT_VECTOR%NUMBER_OF_ROWS=0
      IF(UPDATE_VECTOR) THEN
        DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
          IF(ELEMENT_NUMBER>=1.AND.ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
            SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
              ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP(ELEMENT_NUMBER)
              ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
              ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO node_idx=1,BASIS%NUMBER_OF_NODES
                node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                  derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node)
                  ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                  ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
                ENDDO !derivative_idx
              ENDDO !node_idx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interpolation type of "// &
                & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                & " is invalid for component number "// &
                & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " of rows field variable type "// &
                & TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
            END SELECT
          ELSE
            LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
              & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " of rows field variable type "//TRIM(NUMBER_TO_VSTRING(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
              & ". The element number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !component_idx
        ELEMENT_VECTOR%VECTOR=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise an element vector and deallocate all memory
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE):: ELEMENT_VECTOR !<The element vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(ELEMENT_VECTOR%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_VECTOR%VECTOR)) DEALLOCATE(ELEMENT_VECTOR%VECTOR)
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element vector
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !The element vector to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE",ERR,ERROR,*999)

    ELEMENT_VECTOR%NUMBER_OF_ROWS=0
    ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
       
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets up the element vector for the row field variables.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(ELEMENT_VECTOR,ROWS_FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !<The element vector to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP",ERR,ERROR,*998)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=ROWS_FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS* &
        & ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      IF(ALLOCATED(ELEMENT_VECTOR%ROW_DOFS)) THEN
        CALL FLAG_ERROR("Element vector row dofs is already allocated.",ERR,ERROR,*999)        
      ELSE
        ALLOCATE(ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector row dofs.",ERR,ERROR,*999)
      ENDIF
      IF(ALLOCATED(ELEMENT_VECTOR%VECTOR)) THEN
        CALL FLAG_ERROR("Element vector vector already allocated.",ERR,ERROR,*999)        
      ELSE
        ALLOCATE(ELEMENT_VECTOR%VECTOR(ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector vector.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP")
    RETURN
999 CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(ELEMENT_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP

  !
  !================================================================================================================================
  !

  !>Adds the element matrices and rhs vector into the equations matrices and rhs vector.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,matrix_idx,row_idx
    REAL(DP) :: SUM
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_ELEMENT_ADD()")
#endif

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
      IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
        !Add the element matrices
        DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
          EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
            IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
              !Handle lumped matrices
              IF(EQUATIONS_MATRIX%LUMPED) THEN
                DO row_idx=1,EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS
                  SUM=0.0_DP
                  DO column_idx=1,EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS
                    SUM=SUM+EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx,column_idx)
                    EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx,column_idx)=0.0_DP
                  ENDDO !column_idx
                  EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx,row_idx)=SUM
                  !Add the element matrice into the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(EQUATIONS_MATRIX%MATRIX,EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS(row_idx), &
                    & EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(row_idx),EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx, &
                    & row_idx),ERR,ERROR,*999)
                ENDDO !row_idx
              ELSE
                !Add the element matrice into the distributed equations matrix
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(EQUATIONS_MATRIX%MATRIX,EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
                  & ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
      IF(ASSOCIATED(LINEAR_MATRICES)) THEN
        !Add the element matrices
        DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
            IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
              !Handle lumped matrices
              IF(EQUATIONS_MATRIX%LUMPED) THEN
                DO row_idx=1,EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS
                  SUM=0.0_DP
                  DO column_idx=1,EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS
                    SUM=SUM+EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx,column_idx)
                    EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx,column_idx)=0.0_DP
                  ENDDO !column_idx
                  EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx,row_idx)=SUM
                  !Add the element matrice into the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(EQUATIONS_MATRIX%MATRIX,EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS(row_idx), &
                    & EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(row_idx),EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(row_idx, &
                    & row_idx),ERR,ERROR,*999)
                ENDDO !row_idx
              ELSE
                !Add the element matrice into the distributed equations matrix
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(EQUATIONS_MATRIX%MATRIX,EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                  & EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:EQUATIONS_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
                  & ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
        IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
          IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
            !Add in Jacobian element matrices
            CALL DISTRIBUTED_MATRIX_VALUES_ADD(JACOBIAN_MATRIX%JACOBIAN,JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS(1: &
              & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS(1: &
              & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1: &
              & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS,1:JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS), &
              & ERR,ERROR,*999)
          ENDIF
!!TODO: work out what to do with jacobian and re2 with FD jacobian calc.
        ELSE
          CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
        ENDIF
        IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
          !Add the residual element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(NONLINEAR_MATRICES%RESIDUAL,NONLINEAR_MATRICES%ELEMENT_RESIDUAL%ROW_DOFS(1: &
            & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%NUMBER_OF_ROWS),NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(1:NONLINEAR_MATRICES% &
            & ELEMENT_RESIDUAL%NUMBER_OF_ROWS),ERR,ERROR,*999)
        ENDIF
      ENDIF
      RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        IF(RHS_VECTOR%UPDATE_VECTOR) THEN
          !Add the rhs element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(RHS_VECTOR%VECTOR,RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS(1: &
            & RHS_VECTOR%ELEMENT_VECTOR%NUMBER_OF_ROWS),RHS_VECTOR%ELEMENT_VECTOR%VECTOR(1:RHS_VECTOR% &
            & ELEMENT_VECTOR%NUMBER_OF_ROWS),ERR,ERROR,*999)
        ENDIF
      ENDIF
      SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
      IF(ASSOCIATED(SOURCE_VECTOR)) THEN
        IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
          !Add the rhs element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOURCE_VECTOR%VECTOR,SOURCE_VECTOR%ELEMENT_VECTOR%ROW_DOFS(1: &
            & SOURCE_VECTOR%ELEMENT_VECTOR%NUMBER_OF_ROWS),SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(1:SOURCE_VECTOR% &
            & ELEMENT_VECTOR%NUMBER_OF_ROWS),ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_ELEMENT_ADD()")
#endif
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_ADD")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_ADD",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_ADD")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices and rhs of the element matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          !Calculate the row and columns for the dynamic equations matrices
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(EQUATIONS_MATRIX%ELEMENT_MATRIX,EQUATIONS_MATRIX%UPDATE_MATRIX, &
                  & ELEMENT_NUMBER,ELEMENT_NUMBER,FIELD_VARIABLE,FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          !Calculate the row and columns for the linear equations matrices
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(EQUATIONS_MATRIX%ELEMENT_MATRIX,EQUATIONS_MATRIX%UPDATE_MATRIX, &
                  & ELEMENT_NUMBER,ELEMENT_NUMBER,FIELD_VARIABLE,FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          !Calculate the rows and columns of the Jacobian
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
              FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE
              CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN,JACOBIAN_MATRIX%UPDATE_JACOBIAN, &
                & ELEMENT_NUMBER,ELEMENT_NUMBER,FIELD_VARIABLE,FIELD_VARIABLE,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
            ENDIF
            !Calculate the rows the equations residual
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,NONLINEAR_MATRICES% &
              & UPDATE_RESIDUAL,ELEMENT_NUMBER,FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
          IF(ASSOCIATED(RHS_MAPPING)) THEN
            !Calculate the rows  for the equations RHS
            FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(RHS_VECTOR%ELEMENT_VECTOR,RHS_VECTOR%UPDATE_VECTOR,ELEMENT_NUMBER, &
              & FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations mapping rhs mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
        IF(ASSOCIATED(SOURCE_VECTOR)) THEN
          !Calculate the rows the equations source. The number of rows is not set by the source field so take the number of rows
          !from the RHS vector in the first instance.
          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
          IF(ASSOCIATED(RHS_MAPPING)) THEN
            FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(SOURCE_VECTOR%ELEMENT_VECTOR,SOURCE_VECTOR%UPDATE_VECTOR, &
              & ELEMENT_NUMBER,FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations mapping rhs mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_ELEMENT_CALCULATE()")
#endif
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information and deallocate all memory
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<The equations matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
      IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
        !Finalise the dynamic element matrices
        DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
          EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
            EQUATIONS_MATRIX%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
            EQUATIONS_MATRIX%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
            IF(ALLOCATED(EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS)) DEALLOCATE(EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS)
            IF(ALLOCATED(EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS)) DEALLOCATE(EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS)
            IF(ALLOCATED(EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX)) DEALLOCATE(EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX)
          ELSE
            LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
      IF(ASSOCIATED(LINEAR_MATRICES)) THEN
        !Finalise the linear element matrices
        DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
            EQUATIONS_MATRIX%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
            EQUATIONS_MATRIX%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
            IF(ALLOCATED(EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS)) DEALLOCATE(EQUATIONS_MATRIX%ELEMENT_MATRIX%ROW_DOFS)
            IF(ALLOCATED(EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS)) DEALLOCATE(EQUATIONS_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS)
            IF(ALLOCATED(EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX)) DEALLOCATE(EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX)
          ELSE
            LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
        IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
          JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MAX_NUMBER_OF_ROWS=0
          JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MAX_NUMBER_OF_COLUMNS=0
          IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS)
          IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS)
          IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX)
          JACOBIAN_MATRIX%ELEMENT_RESIDUAL%MAX_NUMBER_OF_ROWS=0
          IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_RESIDUAL%ROW_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_RESIDUAL%ROW_DOFS)
          IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_RESIDUAL%VECTOR)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_RESIDUAL%VECTOR)
        ELSE
          CALL FLAG_ERROR("Nonlinear matrices Jacobian is not associated.",ERR,ERROR,*999)
        ENDIF
        NONLINEAR_MATRICES%ELEMENT_RESIDUAL%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%ROW_DOFS)) DEALLOCATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%ROW_DOFS)
        IF(ALLOCATED(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR)) DEALLOCATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR)
      ENDIF
      RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        !Finalise the element vector
        RHS_VECTOR%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)
      ENDIF
      SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
      IF(ASSOCIATED(SOURCE_VECTOR)) THEN
        !Finalise the element source vector
        SOURCE_VECTOR%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(SOURCE_VECTOR%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(SOURCE_VECTOR%ELEMENT_VECTOR%ROW_DOFS)
        IF(ALLOCATED(SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the equations matrices
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !The equations matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          !Initialise the dynamic element matrices
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(EQUATIONS_MATRIX%ELEMENT_MATRIX,FIELD_VARIABLE,FIELD_VARIABLE, &
                  & ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          !Initialise the linear element matrices
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(matrix_idx)%VARIABLE
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(EQUATIONS_MATRIX%ELEMENT_MATRIX,FIELD_VARIABLE,FIELD_VARIABLE, &
                  & ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          !Initialise the Jacobian element matrices
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
              FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE
              CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(JACOBIAN_MATRIX%ELEMENT_JACOBIAN,FIELD_VARIABLE,FIELD_VARIABLE, &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(JACOBIAN_MATRIX%ELEMENT_RESIDUAL,FIELD_VARIABLE,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
            ENDIF
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          !Initialise the RHS element vector
          RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
          IF(ASSOCIATED(RHS_MAPPING)) THEN
            FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(RHS_VECTOR%ELEMENT_VECTOR,FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("RHS mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
        IF(ASSOCIATED(SOURCE_VECTOR)) THEN
          !Initialise the source element vector. Note that the number of rows in the source vector is taken, for now, from the RHS
          !vector
          IF(ASSOCIATED(RHS_VECTOR)) THEN
            !Initialise the RHS element vector
            RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
            IF(ASSOCIATED(RHS_MAPPING)) THEN
              FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
              CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(SOURCE_VECTOR%ELEMENT_VECTOR,FIELD_VARIABLE,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("RHS mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise a equations matrix and deallocate all memory
  SUBROUTINE EQUATIONS_MATRIX_FINALISE(EQUATIONS_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX !<A pointer to the equations matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("EQUATIONS_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
      IF(ASSOCIATED(EQUATIONS_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_MATRIX%TEMP_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(EQUATIONS_MATRIX%TEMP_VECTOR,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the dynamic equations matrix.
  SUBROUTINE EQUATIONS_MATRIX_DYNAMIC_INITIALISE(DYNAMIC_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES !<A pointer to the dynamic matrices to initialise the dynamic equations matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The dynamic matrix number in the dynamic equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
        EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
          EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
              IF(ASSOCIATED(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
                LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is already associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                ALLOCATE(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrix.",ERR,ERROR,*999)
                EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
                EQUATIONS_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
                EQUATIONS_MATRIX%DYNAMIC_MATRICES=>DYNAMIC_MATRICES
                NULLIFY(EQUATIONS_MATRIX%LINEAR_MATRICES)
                EQUATIONS_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
                EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
                EQUATIONS_MATRIX%LUMPED=.FALSE.
                EQUATIONS_MATRIX%UPDATE_MATRIX=.TRUE.
                EQUATIONS_MATRIX%FIRST_ASSEMBLY=.TRUE.
                EQUATIONS_MATRIX%NUMBER_OF_COLUMNS=DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%EQUATIONS_MATRIX=>EQUATIONS_MATRIX
                NULLIFY(EQUATIONS_MATRIX%MATRIX)
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
                NULLIFY(EQUATIONS_MATRIX%TEMP_VECTOR)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dynamic matrices equations matrices is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified dynamic matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Dynamic matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRIX_FINALISE(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRIX_DYNAMIC_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialise the linear equations matrix.
  SUBROUTINE EQUATIONS_MATRIX_LINEAR_INITIALISE(LINEAR_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES !<A pointer to the linear matrices to initialise the linear equations matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The linear matrix number in the linear equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRIX_LINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(LINEAR_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
        EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
          EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
              IF(ASSOCIATED(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
                LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is already associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                ALLOCATE(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrix.",ERR,ERROR,*999)
                EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
                EQUATIONS_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
                NULLIFY(EQUATIONS_MATRIX%DYNAMIC_MATRICES)
                EQUATIONS_MATRIX%LINEAR_MATRICES=>LINEAR_MATRICES
                EQUATIONS_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
                EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
                EQUATIONS_MATRIX%LUMPED=.FALSE.
                EQUATIONS_MATRIX%UPDATE_MATRIX=.TRUE.
                EQUATIONS_MATRIX%FIRST_ASSEMBLY=.TRUE.
                EQUATIONS_MATRIX%NUMBER_OF_COLUMNS=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%EQUATIONS_MATRIX=>EQUATIONS_MATRIX
                NULLIFY(EQUATIONS_MATRIX%MATRIX)
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
                NULLIFY(EQUATIONS_MATRIX%TEMP_VECTOR)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Linear matrices equations matrices is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified linear matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRIX_LINEAR_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRIX_FINALISE(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRIX_LINEAR_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRIX_LINEAR_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRIX_LINEAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices dynamic matrices and deallocates all memory
  SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_FINALISE(DYNAMIC_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES !<A pointer to the equation matrices dynamic matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
     
    CALL ENTERS("EQUATIONS_MATRICES_DYNAMIC_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
      IF(ALLOCATED(DYNAMIC_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(DYNAMIC_MATRICES%MATRICES,1)
          CALL EQUATIONS_MATRIX_FINALISE(DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(DYNAMIC_MATRICES%MATRICES)
      ENDIF
      IF(ASSOCIATED(DYNAMIC_MATRICES%TEMP_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(DYNAMIC_MATRICES%TEMP_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(DYNAMIC_MATRICES)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_DYNAMIC_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices dynamic matrices
  SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equation matrices to initialise the dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
     
    CALL ENTERS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE",ERR,ERROR,*998)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%DYNAMIC_MATRICES)) THEN
            CALL FLAG_ERROR("Equations matrices dynamic matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices dynamic matrices.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%DYNAMIC_MATRICES%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
            EQUATIONS_MATRICES%DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
            ALLOCATE(EQUATIONS_MATRICES%DYNAMIC_MATRICES%MATRICES(DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices dynamic matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              NULLIFY(EQUATIONS_MATRICES%DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL EQUATIONS_MATRIX_DYNAMIC_INITIALISE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
            NULLIFY(EQUATIONS_MATRICES%DYNAMIC_MATRICES%TEMP_VECTOR)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices equations mapping is not associated.",ERR,ERROR,*998)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_DYNAMIC_FINALISE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices linear matrices and deallocates all memory
  SUBROUTINE EQUATIONS_MATRICES_LINEAR_FINALISE(LINEAR_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES !<A pointer to the equation matrices linear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
     
    CALL ENTERS("EQUATIONS_MATRICES_LINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_MATRICES)) THEN
      IF(ALLOCATED(LINEAR_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(LINEAR_MATRICES%MATRICES,1)
          CALL EQUATIONS_MATRIX_FINALISE(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(LINEAR_MATRICES%MATRICES)
      ENDIF
      DEALLOCATE(LINEAR_MATRICES)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_LINEAR_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_LINEAR_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices linear matrices
  SUBROUTINE EQUATIONS_MATRICES_LINEAR_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equation matrices to initialise the linear matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
     
    CALL ENTERS("EQUATIONS_MATRICES_LINEAR_INITIALISE",ERR,ERROR,*998)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
        IF(ASSOCIATED(LINEAR_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%LINEAR_MATRICES)) THEN
            CALL FLAG_ERROR("Equations matrices linear matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%LINEAR_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices linear matrices.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%LINEAR_MATRICES%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
            EQUATIONS_MATRICES%LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES=LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
            ALLOCATE(EQUATIONS_MATRICES%LINEAR_MATRICES%MATRICES(LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices linear matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              NULLIFY(EQUATIONS_MATRICES%LINEAR_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL EQUATIONS_MATRIX_LINEAR_INITIALISE(EQUATIONS_MATRICES%LINEAR_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices equations mapping is not associated.",ERR,ERROR,*998)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_LINEAR_FINALISE(EQUATIONS_MATRICES%LINEAR_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_LINEAR_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_LINEAR_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices nonlinear matrices and deallocates all memory
  SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_FINALISE(NONLINEAR_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer to the equation matrices nonlinear matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("EQUATIONS_MATRICES_NONLINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
      CALL EQUATIONS_JACOBIAN_FINALISE(NONLINEAR_MATRICES%JACOBIAN,ERR,ERROR,*999)
      IF(ASSOCIATED(NONLINEAR_MATRICES%RESIDUAL)) CALL DISTRIBUTED_VECTOR_DESTROY(NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ERR,ERROR,*999)
      DEALLOCATE(NONLINEAR_MATRICES)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_NONLINEAR_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices nonlinear matrices
  SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equation matrices to initialise the nonlinear matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%NONLINEAR_MATRICES)) THEN
            CALL FLAG_ERROR("Equations matrices nonlinear matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices nonlinear matrices.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
            NULLIFY(EQUATIONS_MATRICES%NONLINEAR_MATRICES%JACOBIAN)
            CALL EQUATIONS_JACOBIAN_INITIALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,ERR,ERROR,*999)
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%UPDATE_RESIDUAL=.TRUE.
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(EQUATIONS_MATRICES%NONLINEAR_MATRICES%RESIDUAL)
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_NONLINEAR_FINALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Outputs the equations matrices
  SUBROUTINE EQUATIONS_MATRICES_OUTPUT(ID,EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    
    CALL ENTERS("EQUATIONS_MATRICES_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"Equations matrices:",ERR,ERROR,*999)
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Dynamic matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of dynamic matrices = ",DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,ERR,ERROR,*999)
          DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
            EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
              CALL WRITE_STRING_VALUE(ID,"Equations matrix : ",matrix_idx,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Linear matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of linear matrices = ",LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,ERR,ERROR,*999)
          DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
            EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
              CALL WRITE_STRING_VALUE(ID,"Equations matrix : ",matrix_idx,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Nonlinear matrices and vector:",ERR,ERROR,*999)
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN        
            CALL WRITE_STRING(ID,"Jacobian matrix:",ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_OUTPUT(ID,JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(NONLINEAR_MATRICES%RESIDUAL)) THEN
            CALL WRITE_STRING(ID,"Residual vector:",ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Nonlinear matrices residual is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          CALL WRITE_STRING(ID,"RHS vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,RHS_VECTOR%VECTOR,ERR,ERROR,*999)
        ENDIF
        SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
        IF(ASSOCIATED(SOURCE_VECTOR)) THEN
          CALL WRITE_STRING(ID,"Source vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,SOURCE_VECTOR%VECTOR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_OUTPUT")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_OUTPUT",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_OUTPUT")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_OUTPUT
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices RHS vector and deallocates all memory
  SUBROUTINE EQUATIONS_MATRICES_RHS_FINALISE(RHS_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR !<A pointer to the equation matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("EQUATIONS_MATRICES_RHS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_VECTOR)) THEN
      IF(ASSOCIATED(RHS_VECTOR%VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(RHS_VECTOR%VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(RHS_VECTOR)
    ENDIF      
     
    CALL EXITS("EQUATIONS_MATRICES_RHS_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_RHS_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_RHS_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_RHS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices RHS vector
  SUBROUTINE EQUATIONS_MATRICES_RHS_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equation matrices to initialise the rhs vector for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_RHS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
        IF(ASSOCIATED(RHS_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%RHS_VECTOR)) THEN
            CALL FLAG_ERROR("Equations matrices RHS vector is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%RHS_VECTOR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices RHS vector.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%RHS_VECTOR%UPDATE_VECTOR=.TRUE.
            EQUATIONS_MATRICES%RHS_VECTOR%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(EQUATIONS_MATRICES%RHS_VECTOR%VECTOR)
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(EQUATIONS_MATRICES%RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices equation mapping is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_RHS_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_RHS_FINALISE(EQUATIONS_MATRICES%RHS_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_RHS_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_RHS_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_RHS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the equations matrices source vector and deallocates all memory
  SUBROUTINE EQUATIONS_MATRICES_SOURCE_FINALISE(SOURCE_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR !<A pointer to the equation matrices source vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("EQUATIONS_MATRICES_SOURCE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOURCE_VECTOR)) THEN
      IF(ASSOCIATED(SOURCE_VECTOR%VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(SOURCE_VECTOR%VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(SOURCE_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(SOURCE_VECTOR)
    ENDIF      
     
    CALL EXITS("EQUATIONS_MATRICES_SOURCE_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_SOURCE_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_SOURCE_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_SOURCE_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations matrices source vector
  SUBROUTINE EQUATIONS_MATRICES_SOURCE_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equation matrices to initialise the source vector for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_SOURCE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%SOURCE_VECTOR)) THEN
            CALL FLAG_ERROR("Equations matrices source vector is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%SOURCE_VECTOR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices source vector.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%SOURCE_VECTOR%UPDATE_VECTOR=.TRUE.
            EQUATIONS_MATRICES%SOURCE_VECTOR%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(EQUATIONS_MATRICES%SOURCE_VECTOR%VECTOR)
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(EQUATIONS_MATRICES%SOURCE_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrices equation mapping is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_SOURCE_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_SOURCE_FINALISE(EQUATIONS_MATRICES%SOURCE_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_SOURCE_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_SOURCE_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_SOURCE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Sets the lumping of the linear equations matrices
  SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES,LUMPING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: LUMPING_TYPE(:) !<LUMPING_TYPE(matrix_idx). The lumping type for the matrix_idx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have already been finished.",ERR,ERROR,*999)
      ELSE
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          IF(SIZE(LUMPING_TYPE,1)==DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                SELECT CASE(LUMPING_TYPE(matrix_idx))
                CASE(EQUATIONS_MATRIX_UNLUMPED)
                  EQUATIONS_MATRIX%LUMPED=.FALSE.
                CASE(EQUATIONS_MATRIX_LUMPED)
                  EQUATIONS_MATRIX%LUMPED=.TRUE.        
                CASE DEFAULT
                  LOCAL_ERROR="The specified lumping type of "//TRIM(NUMBER_TO_VSTRING(LUMPING_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the lumping type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(LUMPING_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the dynamic equations matrices
  SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have already been finished.",ERR,ERROR,*999)
      ELSE
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          IF(SIZE(STORAGE_TYPE,1)==DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
            DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                SELECT CASE(STORAGE_TYPE(matrix_idx))
                CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
                CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
                CASE DEFAULT
                  LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the dynamic matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the linear equations matrices
  SUBROUTINE EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          IF(SIZE(STORAGE_TYPE,1)==LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
            DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                SELECT CASE(STORAGE_TYPE(matrix_idx))
                CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
                CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
                CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                  EQUATIONS_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
                CASE DEFAULT
                  LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the linear matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of linear matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET

 !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE !< The storage type for the Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            SELECT CASE(STORAGE_TYPE)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              JACOBIAN_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
            CASE DEFAULT
              LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE,"*",ERR,ERROR))// &
                & " for the Jacobian matrix is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the dynamic equations matrices
  SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The storage type for the matrix_idx'th dynamic equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
           DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
              EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                SELECT CASE(STRUCTURE_TYPE(matrix_idx))
                CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
                CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_FEM_STRUCTURE
                CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for dynamic matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the linear equations matrices
  SUBROUTINE EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The storage type for the matrix_idx'th linear equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
           DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
              EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                SELECT CASE(STRUCTURE_TYPE(matrix_idx))
                CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
                CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_FEM_STRUCTURE
                CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for linear matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of linear matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//")."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE !<STRUCTURE_TYPE. The storage type for the Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN          
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            SELECT CASE(STRUCTURE_TYPE)
            CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
              JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
            CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
              JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_FEM_STRUCTURE
            CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
              JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
            CASE DEFAULT
              LOCAL_ERROR="The specified strucutre type of "// &
                & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE,"*",ERR,ERROR))//" for the Jacobian matrix is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Finalise the equations matrices and deallocate all memory.
  SUBROUTINE EQUATIONS_MATRICES_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   
    CALL ENTERS("EQUATIONS_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      CALL EQUATIONS_MATRICES_DYNAMIC_FINALISE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_LINEAR_FINALISE(EQUATIONS_MATRICES%LINEAR_MATRICES,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_NONLINEAR_FINALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_RHS_FINALISE(EQUATIONS_MATRICES%RHS_VECTOR,ERR,ERROR,*999)      
      CALL EQUATIONS_MATRICES_SOURCE_FINALISE(EQUATIONS_MATRICES%SOURCE_VECTOR,ERR,ERROR,*999)      
      DEALLOCATE(EQUATIONS_MATRICES)
    ENDIF
       
    CALL EXITS("EQUATIONS_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the equations matrices for the equations.
  SUBROUTINE EQUATIONS_MATRICES_INITIALISE(EQUATIONS,ERR,ERROR,*)
    
     !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to initialise the equations matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MATRICES)) THEN
        CALL FLAG_ERROR("Equations matrices is already associated for this equations.",ERR,ERROR,*998)
      ELSE
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
            ALLOCATE(EQUATIONS%EQUATIONS_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations equations matrices.",ERR,ERROR,*999)
            EQUATIONS%EQUATIONS_MATRICES%EQUATIONS=>EQUATIONS
            EQUATIONS%EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED=.FALSE.
            EQUATIONS%EQUATIONS_MATRICES%EQUATIONS_MAPPING=>EQUATIONS_MAPPING
            NULLIFY(EQUATIONS%EQUATIONS_MATRICES%SOLVER_MAPPING)
            EQUATIONS%EQUATIONS_MATRICES%NUMBER_OF_ROWS=EQUATIONS_MAPPING%NUMBER_OF_ROWS
            EQUATIONS%EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS=EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
            EQUATIONS%EQUATIONS_MATRICES%NUMBER_OF_GLOBAL_ROWS=EQUATIONS_MAPPING%NUMBER_OF_GLOBAL_ROWS
            NULLIFY(EQUATIONS%EQUATIONS_MATRICES%DYNAMIC_MATRICES)
            NULLIFY(EQUATIONS%EQUATIONS_MATRICES%LINEAR_MATRICES)
            NULLIFY(EQUATIONS%EQUATIONS_MATRICES%NONLINEAR_MATRICES)
            NULLIFY(EQUATIONS%EQUATIONS_MATRICES%RHS_VECTOR)
            NULLIFY(EQUATIONS%EQUATIONS_MATRICES%SOURCE_VECTOR)            
            CALL EQUATIONS_MATRICES_DYNAMIC_INITIALISE(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)            
            CALL EQUATIONS_MATRICES_LINEAR_INITIALISE(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)            
            CALL EQUATIONS_MATRICES_NONLINEAR_INITIALISE(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)            
            CALL EQUATIONS_MATRICES_RHS_INITIALISE(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)            
            CALL EQUATIONS_MATRICES_SOURCE_INITIALISE(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)            
          ELSE
            CALL FLAG_ERROR("Equations mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_MATRICES_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS%EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialise the values of the equations matrices and vectors to the given value e.g., 0.0_DP
  SUBROUTINE EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,SELECTION_TYPE,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices to initialise the values for
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The selection of equations matrices to be initialised \see EQUATIONS_MATRICES_ROUTINES::SelectMatricesTypes,EQUATION_MATRICES_ROUTINES
    REAL(DP), INTENT(IN) :: VALUE !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    
    CALL ENTERS("EQUATIONS_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY) THEN
        DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
        IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
          DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
            EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
              IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(EQUATIONS_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY) THEN
        LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
        IF(ASSOCIATED(LINEAR_MATRICES)) THEN
          DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
            EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
              IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(EQUATIONS_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_JACOBIAN_ONLY) THEN
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIAN
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
              CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(JACOBIAN_MATRIX%JACOBIAN,VALUE,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RESIDUAL_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_VECTORS_ONLY) THEN
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
            CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(NONLINEAR_MATRICES%RESIDUAL,VALUE,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RHS_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RHS_SOURCE_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_VECTORS_ONLY) THEN    
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          IF(RHS_VECTOR%UPDATE_VECTOR) THEN
            CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(RHS_VECTOR%VECTOR,VALUE,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_DYNAMIC_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_LINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_SOURCE_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RHS_SOURCE_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_VECTORS_ONLY) THEN    
        SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
        IF(ASSOCIATED(SOURCE_VECTOR)) THEN
          IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
            CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(SOURCE_VECTOR%VECTOR,VALUE,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EQUATIONS_MATRICES_VALUES_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_VALUES_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_VALUES_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_VALUES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a equations matrix.
  SUBROUTINE EQUATIONS_MATRIX_STRUCTURE_CALCULATE(EQUATIONS_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES,list,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX !<A pointer to the equations matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    type(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::  column_idx,DUMMY_ERR,elem_idx,global_column,local_column,local_ny,MATRIX_NUMBER,mk,mp,ne,nh,nh2,nn,nnk,np, &
      & NUMBER_OF_COLUMNS,nyy,nyyg,npg,nhg,local_cols,local_dof
    INTEGER(INTG), ALLOCATABLE :: COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES,DOMAIN_NODES_GLOBAL
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: DEPENDENT_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    INTEGER(INTG), POINTER :: BOUNDARY_NODES_LIST(:)
    integer(INTG),allocatable:: row_array(:)
    CALL ENTERS("EQUATIONS_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*998)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
      IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
        IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
          MATRIX_NUMBER=EQUATIONS_MATRIX%MATRIX_NUMBER
          SELECT CASE(EQUATIONS_MATRIX%STRUCTURE_TYPE)
          CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
            CALL FLAG_ERROR("There is no structure to calculate for a matrix with no structure.",ERR,ERROR,*998)
          CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
            SELECT CASE(EQUATIONS_MATRIX%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
              DYNAMIC_MATRICES=>EQUATIONS_MATRIX%DYNAMIC_MATRICES
              IF(ASSOCIATED(DYNAMIC_MATRICES).OR.ASSOCIATED(LINEAR_MATRICES)) THEN
                IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                  EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
                ELSE
                  EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
                ENDIF
                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                  EQUATIONS=>EQUATIONS_MATRICES%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
                    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                      DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
                      LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                      IF(ASSOCIATED(DYNAMIC_MAPPING).OR.ASSOCIATED(LINEAR_MAPPING)) THEN
                        EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                        IF(ASSOCIATED(EQUATIONS_SET)) THEN
                          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                            IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                              FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%VARIABLE
                            ELSE
                              FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(MATRIX_NUMBER)%VARIABLE
                            ENDIF
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DEPENDENT_DOFS_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                              IF(ASSOCIATED(DEPENDENT_DOFS_DOMAIN_MAPPING)) THEN
                                DEPENDENT_DOFS_PARAM_MAPPING=>FIELD_VARIABLE%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(DEPENDENT_DOFS_PARAM_MAPPING)) THEN
                                  !Allocate lists
                                  ALLOCATE(COLUMN_INDICES_LISTS(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",ERR,ERROR,*999)
                                  !Allocate row indices
                                  ALLOCATE(ROW_INDICES(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",ERR,ERROR,*999)
                                  ROW_INDICES(1)=1
                                  
                                  !First, loop over the rows and calculate the number of non-zeros
                                  NUMBER_OF_NON_ZEROS=0
                                  DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                    IF(DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(1,local_ny)==FIELD_NODE_DOF_TYPE) THEN
                                      nyy=DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(2,local_ny)!value for a particular field dof (local_ny)
                                      np=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(2,nyy)!node number (np) of the field parameter
                                      nh=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(3,nyy)!component number (nh) of the field parameter
                                      DOMAIN_NODES=>FIELD_VARIABLE%COMPONENTS(nh)%DOMAIN%TOPOLOGY%NODES
                                      
                                      !Set up list
                                      NULLIFY(COLUMN_INDICES_LISTS(local_ny)%PTR)
                                      CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                                      CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,DOMAIN_NODES%NODES(np)% &
                                        & NUMBER_OF_SURROUNDING_ELEMENTS*FIELD_VARIABLE%COMPONENTS(nh)% &
                                        & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                                      CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      !Loop over all elements containing the dof
                                      DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                                        ne=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
                                        DO nh2=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                          DOMAIN_ELEMENTS=>FIELD_VARIABLE%COMPONENTS(nh2)%DOMAIN%TOPOLOGY%ELEMENTS
                                          BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                                          DO nn=1,BASIS%NUMBER_OF_NODES
                                            mp=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                            DO nnk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                                              mk=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                              !Find the local and global column and add the global column to the indices list
                                              local_column=FIELD_VARIABLE%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP(mk,mp)
                                              global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                          
                                              CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                                                
                                            ENDDO !mk
                                          ENDDO !nn
                                        ENDDO !nh2
                                      ENDDO !elem_idx
                                      CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS, &
                                        & ERR,ERROR,*999)
                                      NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                      ROW_INDICES(local_ny+1)=NUMBER_OF_NON_ZEROS+1
                                    ELSE
                                      LOCAL_ERROR="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                                        & " is not a node based dof."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ENDDO !local_ny
                                  
                                  
                                  !Allocate and setup the column locations
                                  ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)

                                  ALLOCATE(list(DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL))

                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",ERR,ERROR,*999)
                                  DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                   
                                    CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS,COLUMNS, &
                                      & ERR,ERROR,*999)        
                                    DO column_idx=1,NUMBER_OF_COLUMNS
                                      !COLUMNS store the list of nonzero column indices for each local row (local_ny)
                                      COLUMN_INDICES(ROW_INDICES(local_ny)+column_idx-1)=COLUMNS(column_idx) 

                                      ! global to local columns
                                       IF(ASSOCIATED(LINEAR_MAPPING).OR.ASSOCIATED(DYNAMIC_MAPPING)) THEN
                                         IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                                           local_cols=equations_matrices%equations_mapping%dynamic_mapping &
                                             & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                             & (COLUMNS(column_idx))%LOCAL_NUMBER(1)
                                           local_dof = local_cols
                                           ! Column to dof mapping?
                                           !local_dof=equations_matrices%equations_mapping%dynamic_mapping% &
                                            ! & equations_matrix_to_var_maps(1)%column_to_dof_map(local_cols)
                                         ELSE
                                           local_cols=equations_matrices%equations_mapping%linear_mapping &
                                             & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                             & (COLUMNS(column_idx))%LOCAL_NUMBER(1)
                                           local_dof = local_cols
                                         ENDIF
                                       ENDIF
                                       nyyg=DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(2,local_dof)
                                       npg=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(2,nyyg)
                                       nhg=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(3,nyyg)
                                       DOMAIN_NODES=>FIELD_VARIABLE%COMPONENTS(nhg)%DOMAIN%TOPOLOGY%NODES
                            
                                      ! Check whether boundary node    
                                      IF(DOMAIN_NODES%NODES(npg)%BOUNDARY_NODE)THEN
                                        CALL LinkedList_Add(list(COLUMNS(column_idx)),local_ny)
                                      ENDIF
                                    
                                    ENDDO !column_idx
                                    DEALLOCATE(COLUMNS)                                    
                                  ENDDO !local_ny

                                 
                                  IF(DIAGNOSTICS1) THEN
                                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix structure:",ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix number : ",MATRIX_NUMBER, &
                                      & ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                      & DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                      & DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                      & NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                                    IF(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                                      & DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL/=0) THEN
                                      SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(DEPENDENT_DOFS_DOMAIN_MAPPING% &
                                        & TOTAL_NUMBER_OF_LOCAL*DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,DP)*100.0_DP
                                      CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F6.2", &
                                        & ERR,ERROR,*999)
                                    ENDIF
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DEPENDENT_DOFS_DOMAIN_MAPPING% &
                                      & TOTAL_NUMBER_OF_LOCAL+1,8,8,ROW_INDICES,'("  Row indices    :",8(X,I13))', &
                                      & '(18X,8(X,I13))',ERR,ERROR,*999)
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8,COLUMN_INDICES, &
                                      & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Dependent dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Dependent dofs domain mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Either equations mapping dynamic mapping or linear mapping is not associated.", &
                          & ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dynamic or linear matrices equations matrices is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Either equations matrix dynamic or linear matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
            CALL FLAG_ERROR("There is not structure to calculate for a diagonal matrix.",ERR,ERROR,*998)
          CASE DEFAULT
            LOCAL_ERROR="The matrix structure type of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MATRIX%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Column indices is already associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indieces is already associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("EQUATIONS_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ALLOCATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(local_ny)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_ny
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
998 CALL ERRORS("EQUATIONS_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRIX_STRUCTURE_CALCULATE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRIX_STRUCTURE_CALCULATE

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a Jacobian matrix.
  SUBROUTINE JACOBIAN_MATRIX_STRUCTURE_CALCULATE(JACOBIAN_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX !<A pointer to the Jacobian matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::  column_idx,DUMMY_ERR,elem_idx,global_column,local_column,local_ny,mk,mp,ne,nh,nh2,nn,nnk,np, &
      & NUMBER_OF_COLUMNS,nyy
    INTEGER(INTG), ALLOCATABLE :: COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: DEPENDENT_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("JACOBIAN_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*998)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
      IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
        IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
          SELECT CASE(JACOBIAN_MATRIX%STRUCTURE_TYPE)
          CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*998)
          CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
            SELECT CASE(JACOBIAN_MATRIX%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
              IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                  EQUATIONS=>EQUATIONS_MATRICES%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
                    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                      NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                      IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                        EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                        IF(ASSOCIATED(EQUATIONS_SET)) THEN
                          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                            FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DEPENDENT_DOFS_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                              IF(ASSOCIATED(DEPENDENT_DOFS_DOMAIN_MAPPING)) THEN
                                DEPENDENT_DOFS_PARAM_MAPPING=>FIELD_VARIABLE%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(DEPENDENT_DOFS_PARAM_MAPPING)) THEN
                                  !Allocate lists
                                  ALLOCATE(COLUMN_INDICES_LISTS(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",ERR,ERROR,*999)
                                  !Allocate row indices
                                  ALLOCATE(ROW_INDICES(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",ERR,ERROR,*999)
                                  ROW_INDICES(1)=1
                                  !First, loop over the rows and calculate the number of non-zeros
                                  NUMBER_OF_NON_ZEROS=0
                                  DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                    SELECT CASE(DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(1,local_ny))
                                    CASE(FIELD_CONSTANT_INTERPOLATION)
                                      CALL FLAG_ERROR("Constant interpolation is not implemented yet.",ERR,ERROR,*999)
                                    CASE(FIELD_NODE_DOF_TYPE)
                                      nyy=DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(2,local_ny)
                                      np=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(2,nyy)
                                      nh=DEPENDENT_DOFS_PARAM_MAPPING%NODE_DOF2PARAM_MAP(3,nyy)
                                      DOMAIN_NODES=>FIELD_VARIABLE%COMPONENTS(nh)%DOMAIN%TOPOLOGY%NODES
                                      DOMAIN_ELEMENTS=>FIELD_VARIABLE%COMPONENTS(nh)%DOMAIN%TOPOLOGY%ELEMENTS
                                      !Set up list
                                      NULLIFY(COLUMN_INDICES_LISTS(local_ny)%PTR)
                                      CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                                      CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,DOMAIN_NODES%NODES(np)% &
                                        & NUMBER_OF_SURROUNDING_ELEMENTS*FIELD_VARIABLE%COMPONENTS(nh)% &
                                        & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                                      CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      !Loop over all elements containing the dof
                                      DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                                        ne=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
                                        DO nh2=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                          SELECT CASE(FIELD_VARIABLE%COMPONENTS(nh2)%INTERPOLATION_TYPE)
                                          CASE(FIELD_CONSTANT_INTERPOLATION)
                                            ! do nothing? this will probably never be encountered...?
                                          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                            local_column=FIELD_VARIABLE%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                              & ELEMENT_PARAM2DOF_MAP(ne)
                                            global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                            CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                                          CASE(FIELD_NODE_BASED_INTERPOLATION)
                                            DOMAIN_ELEMENTS=>FIELD_VARIABLE%COMPONENTS(nh2)%DOMAIN%TOPOLOGY%ELEMENTS
                                            BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                                            DO nn=1,BASIS%NUMBER_OF_NODES
                                              mp=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                              DO nnk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                                                mk=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                                !Find the local and global column and add the global column to the indices list
                                                local_column=FIELD_VARIABLE%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP(mk,mp)
                                                global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                                CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                                              ENDDO !mk
                                            ENDDO !nn
                                          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                            CALL FLAG_ERROR("Grid point based interpolation is not implemented yet.",& 
                                              & ERR,ERROR,*999)
                                          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                            CALL FLAG_ERROR("Gauss point based interpolation is not implemented yet.",&
                                              & ERR,ERROR,*999)
                                          CASE DEFAULT
                                            LOCAL_ERROR="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                                              & " has invalid interpolation type."
                                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                          END SELECT
                                        ENDDO !nh2
                                      ENDDO !elem_idx
                                      CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS, &
                                        & ERR,ERROR,*999)
                                      NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                      ROW_INDICES(local_ny+1)=NUMBER_OF_NON_ZEROS+1
                                    CASE(FIELD_ELEMENT_DOF_TYPE)
                                      ! row corresponds to a variable that's element-wisely interpolated
                                      nyy=DEPENDENT_DOFS_PARAM_MAPPING%DOF_TYPE(2,local_ny)          ! nyy = index in ELEMENT_DOF2PARAM_MAP
                                      nh=DEPENDENT_DOFS_PARAM_MAPPING%ELEMENT_DOF2PARAM_MAP(2,nyy)   ! current variable component
                                      ne=DEPENDENT_DOFS_PARAM_MAPPING%ELEMENT_DOF2PARAM_MAP(1,nyy)   ! current element (i.e. corresponds to current dof)
                                      DOMAIN_ELEMENTS=>FIELD_VARIABLE%COMPONENTS(nh)%DOMAIN%TOPOLOGY%ELEMENTS
                                      BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                                      !Set up list
                                      NULLIFY(COLUMN_INDICES_LISTS(local_ny)%PTR)
                                      CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                                      CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_ny)%PTR, &
                                        & FIELD_VARIABLE%COMPONENTS(nh)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS+1, &
                                        & ERR,ERROR,*999) ! size = all nodal dofs + itself
                                      CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      DO nh2=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                        SELECT CASE(FIELD_VARIABLE%COMPONENTS(nh2)%INTERPOLATION_TYPE)
                                        CASE(FIELD_CONSTANT_INTERPOLATION)
                                          CALL FLAG_ERROR("Constant interpolation is not implemented yet.",ERR,ERROR,*999)
                                        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                          ! it's assumed that element-based variables arne't directly coupled
                                          ! put a diagonal entry
                                          local_column=FIELD_VARIABLE%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                            & ELEMENT_PARAM2DOF_MAP(ne)
                                          global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                          CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                                        CASE(FIELD_NODE_BASED_INTERPOLATION)
                                          ! loop over all nodes in the element (and dofs belonging to them)
                                          DO nn=1,BASIS%NUMBER_OF_NODES
                                            mp=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                            DO nnk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                                              mk=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                              !Find the local and global column and add the global column to the indices list
                                              local_column=FIELD_VARIABLE%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP(mk,mp)
                                              global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                              CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                                            ENDDO !mk
                                          ENDDO !nn
                                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                          CALL FLAG_ERROR("Grid point based interpolation is not implemented yet.",ERR,ERROR,*999)
                                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                          CALL FLAG_ERROR("Gauss point based interpolation is not implemented yet.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                                            & " has invalid interpolation type."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                      ENDDO !nh2
                                      ! clean up the list
                                      CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(local_ny)%PTR,ERR,ERROR,*999)
                                      CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS, &
                                        & ERR,ERROR,*999)
                                      NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                      ROW_INDICES(local_ny+1)=NUMBER_OF_NON_ZEROS+1
                                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                      CALL FLAG_ERROR("Grid point based interpolation is not implemented yet.",ERR,ERROR,*999)
                                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                      CALL FLAG_ERROR("Gauss point based interpolation is not implemented yet.",ERR,ERROR,*999)
                                    CASE DEFAULT
                                      LOCAL_ERROR="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                                        & " has an invalid type."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                  ENDDO !local_ny
                                  !Allocate and setup the column locations
                                  ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",ERR,ERROR,*999)
                                  DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                    CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS,COLUMNS, &
                                      & ERR,ERROR,*999)
                                    DO column_idx=1,NUMBER_OF_COLUMNS
                                      COLUMN_INDICES(ROW_INDICES(local_ny)+column_idx-1)=COLUMNS(column_idx)
                                    ENDDO !column_idx
                                    DEALLOCATE(COLUMNS)
                                  ENDDO !local_ny
                                  IF(DIAGNOSTICS1) THEN
                                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Jacobian matrix structure:",ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                      & DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                      & DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                      & NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                                    IF(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                                      & DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL/=0) THEN
                                      SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(DEPENDENT_DOFS_DOMAIN_MAPPING% &
                                        & TOTAL_NUMBER_OF_LOCAL*DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,DP)*100.0_DP
                                      CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F5.2", &
                                        & ERR,ERROR,*999)
                                    ENDIF
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DEPENDENT_DOFS_DOMAIN_MAPPING% &
                                      & TOTAL_NUMBER_OF_LOCAL+1,8,8,ROW_INDICES,'("  Row indices    :",8(X,I13))', &
                                      & '(18X,8(X,I13))',ERR,ERROR,*999)
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8,COLUMN_INDICES, &
                                      & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Dependent dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Dependent dofs domain mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations matrix nonlinear matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(JACOBIAN_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The matrix structure type of "// &
              & TRIM(NUMBER_TO_VSTRING(JACOBIAN_MATRIX%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Column indices is already associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indices is already associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
    ENDIF

      
    CALL EXITS("JACOBIAN_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ALLOCATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(local_ny)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_ny
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
998 CALL ERRORS("JACOBIAN_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    CALL EXITS("JACOBIAN_MATRIX_STRUCTURE_CALCULATE")
    RETURN 1
  END SUBROUTINE JACOBIAN_MATRIX_STRUCTURE_CALCULATE

  !
  !================================================================================================================================
  !
 
END MODULE EQUATIONS_MATRICES_ROUTINES
