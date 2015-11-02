!> \file
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
!> Contributor(s): David Ladd
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

#include "macros.h"  

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
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_NODAL_STRUCTURE=4 !<Nodal matrix structure. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
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

  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_JacobianCalculationTypes EQUATIONS_MATRICES_ROUTINES:JacobianCalculationTypes
  !> \brief Jacobian calculation types
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED=1 !<Use finite differencing to calculate the Jacobian
  INTEGER(INTG), PARAMETER :: EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED=2 !<Use an analytic Jacobian evaluation
  !>@}
  !>@}

  !Module types

  !Module variables

  !Interfaces

  INTERFACE EquationsMatrices_CreateFinish
    MODULE PROCEDURE EQUATIONS_MATRICES_CREATE_FINISH
  END INTERFACE EquationsMatrices_CreateFinish

  INTERFACE EquationsMatrices_CreateStart
    MODULE PROCEDURE EQUATIONS_MATRICES_CREATE_START
  END INTERFACE EquationsMatrices_CreateStart

  INTERFACE EquationsMatrices_Destroy
    MODULE PROCEDURE EQUATIONS_MATRICES_DESTROY
  END INTERFACE EquationsMatrices_Destroy

  INTERFACE EquationsMatrices_DynamicLumpingTypeSet
    MODULE PROCEDURE EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET
  END INTERFACE EquationsMatrices_DynamicLumpingTypeSet

  INTERFACE EquationsMatrices_DynamicStorageTypeSet
    MODULE PROCEDURE EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET
  END INTERFACE EquationsMatrices_DynamicStorageTypeSet

  INTERFACE EquationsMatrices_ElementInitialise
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_INITIALISE
  END INTERFACE EquationsMatrices_ElementInitialise

  INTERFACE EquationsMatrices_ElementFinalise
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_FINALISE
  END INTERFACE EquationsMatrices_ElementFinalise

  INTERFACE EquationsMatrices_ElementAdd
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_ADD
  END INTERFACE EquationsMatrices_ElementAdd
  
  INTERFACE EquationsMatrices_JacobianElementAdd
    MODULE PROCEDURE EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD
  END INTERFACE EquationsMatrices_JacobianElementAdd
  
  INTERFACE EquationsMatrices_ElementCalculate
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_CALCULATE
  END INTERFACE EquationsMatrices_ElementCalculate

  INTERFACE EquationsMatrices_ValuesInitialise
    MODULE PROCEDURE EQUATIONS_MATRICES_VALUES_INITIALISE
  END INTERFACE EquationsMatrices_ValuesInitialise

  INTERFACE EquationsMatrices_ElementMatrixFinalise
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE
  END INTERFACE EquationsMatrices_ElementMatrixFinalise

  INTERFACE EquationsMatrices_ElementMatrixCalculate
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE
  END INTERFACE EquationsMatrices_ElementMatrixCalculate

  INTERFACE EquationsMatrices_ElementMatrixSetup
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP
  END INTERFACE EquationsMatrices_ElementMatrixSetup

  INTERFACE EquationsMatrices_ElementVectorFinalise
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE
  END INTERFACE EquationsMatrices_ElementVectorFinalise

  INTERFACE EquationsMatrices_ElementVectorCalculate
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE
  END INTERFACE EquationsMatrices_ElementVectorCalculate

  INTERFACE EquationsMatrices_ElementVectorSetup
    MODULE PROCEDURE EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP
  END INTERFACE EquationsMatrices_ElementVectorSetup

  INTERFACE EquationsMatrices_LinearStorageTypeSet
    MODULE PROCEDURE EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET
  END INTERFACE EquationsMatrices_LinearStorageTypeSet

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  INTERFACE EquationsMatrices_NonlinearStorageTypeSet
    MODULE PROCEDURE EquationsMatrices_NonlinearStorageTypeSet0
    MODULE PROCEDURE EquationsMatrices_NonlinearStorageTypeSet1
  END INTERFACE EquationsMatrices_NonlinearStorageTypeSet

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  INTERFACE EquationsMatrices_NonlinearStructureTypeSet
    MODULE PROCEDURE EquationsMatrices_NonlinearStructureTypeSet0
    MODULE PROCEDURE EquationsMatrices_NonlinearStructureTypeSet1
  END INTERFACE EquationsMatrices_NonlinearStructureTypeSet

  INTERFACE EquationsMatrices_Output
    MODULE PROCEDURE EQUATIONS_MATRICES_OUTPUT
  END INTERFACE EquationsMatrices_Output

  INTERFACE EquationsMatrices_JacobianOutput
    MODULE PROCEDURE EQUATIONS_MATRICES_JACOBIAN_OUTPUT
  END INTERFACE EquationsMatrices_JacobianOutput

  PUBLIC EQUATIONS_MATRIX_NO_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE, &
    & EQUATIONS_MATRIX_NODAL_STRUCTURE

  PUBLIC EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED

  PUBLIC EQUATIONS_MATRICES_SPARSE_MATRICES,EQUATIONS_MATRICES_FULL_MATRICES

  PUBLIC EQUATIONS_MATRICES_ALL,EQUATIONS_MATRICES_LINEAR_ONLY,EQUATIONS_MATRICES_NONLINEAR_ONLY,EQUATIONS_MATRICES_JACOBIAN_ONLY, &
    & EQUATIONS_MATRICES_RESIDUAL_ONLY,EQUATIONS_MATRICES_RHS_ONLY,EQUATIONS_MATRICES_SOURCE_ONLY, &
    & EQUATIONS_MATRICES_RHS_RESIDUAL_ONLY,EQUATIONS_MATRICES_RHS_SOURCE_ONLY,EQUATIONS_MATRICES_RESIDUAL_SOURCE_ONLY, &
    & EQUATIONS_MATRICES_VECTORS_ONLY
  
  PUBLIC EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED,EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED

  PUBLIC EQUATIONS_MATRICES_CREATE_FINISH,EQUATIONS_MATRICES_CREATE_START

  PUBLIC EquationsMatrices_CreateFinish,EquationsMatrices_CreateStart

  PUBLIC EQUATIONS_MATRICES_DESTROY

  PUBLIC EquationsMatrices_Destroy

  PUBLIC EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET,EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET

  PUBLIC EquationsMatrices_DynamicLumpingTypeSet,EquationsMatrices_DynamicStorageTypeSet,EquationsMatrices_DynamicStructureTypeSet

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC EQUATIONS_MATRICES_ELEMENT_INITIALISE,EQUATIONS_MATRICES_ELEMENT_FINALISE

  PUBLIC EquationsMatrices_ElementInitialise,EquationsMatrices_ElementFinalise
  
  PUBLIC EQUATIONS_MATRICES_ELEMENT_ADD,EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD,EQUATIONS_MATRICES_ELEMENT_CALCULATE, &
    & EQUATIONS_MATRICES_VALUES_INITIALISE

  PUBLIC EquationsMatrices_ElementAdd,EquationsMatrices_JacobianElementAdd,EquationsMatrices_ElementCalculate, &
    & EquationsMatrices_ValuesInitialise

  PUBLIC EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE

  PUBLIC EquationsMatrices_ElementMatrixFinalise,EquationsMatrices_ElementMatrixInitialise
  
  PUBLIC EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE,EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP

  PUBLIC EquationsMatrices_ElementMatrixCalculate,EquationsMatrices_ElementMatrixSetup

  PUBLIC EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE

  PUBLIC EquationsMatrices_ElementVectorFinalise,EquationsMatrices_ElementVectorInitialise

  PUBLIC EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE,EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP

  PUBLIC EquationsMatrices_ElementVectorCalculate,EquationsMatrices_ElementVectorSetup

  PUBLIC EquationsMatrices_JacobianTypesSet

  PUBLIC EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET

  PUBLIC EquationsMatrices_LinearStorageTypeSet,EquationsMatrices_LinearStructureTypeSet

  PUBLIC EquationsMatrices_NodalInitialise,EquationsMatrices_NodalFinalise

  PUBLIC EquationsMatrices_NodeAdd,EquationsMatrices_JacobianNodeAdd,EquationsMatrices_NodalCalculate

  PUBLIC EquationsMatrices_NodalMatrixFinalise,EquationsMatrices_NodalMatrixInitialise

  PUBLIC EquationsMatrices_NodalMatrixCalculate,EquationsMatrices_NodalMatrixSetup

  PUBLIC EquationsMatrices_NodalVectorFinalise,EquationsMatrices_NodalVectorInitialise

  PUBLIC EquationsMatrices_NodalVectorCalculate,EquationsMatrices_NodalVectorSetup

  PUBLIC EquationsMatrices_NonlinearStorageTypeSet,EquationsMatrices_NonlinearStructureTypeSet

  PUBLIC EQUATIONS_MATRICES_OUTPUT,EQUATIONS_MATRICES_JACOBIAN_OUTPUT

  PUBLIC EquationsMatrices_Output,EquationsMatrices_JacobianOutput

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
    
    ENTERS("EQUATIONS_JACOBIAN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_JACOBIAN)) THEN
      IF(ASSOCIATED(EQUATIONS_JACOBIAN%JACOBIAN)) CALL DISTRIBUTED_MATRIX_DESTROY(EQUATIONS_JACOBIAN%JACOBIAN,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_JACOBIAN%ELEMENT_JACOBIAN,ERR,ERROR,*999)
      CALL EquationsMatrices_NodalMatrixFinalise(EQUATIONS_JACOBIAN%NodalJacobian,ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_JACOBIAN_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_JACOBIAN_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_JACOBIAN_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the equations Jacobian.
  SUBROUTINE EQUATIONS_JACOBIAN_INITIALISE(NONLINEAR_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES !<A pointer to the equations matrices nonlinear matrices to initialise the Jacobian for
    INTEGER(INTG), INTENT(IN) :: MATRIX_NUMBER !<The index of the Jacobian matrix to initialise for the nonlinear matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("EQUATIONS_JACOBIAN_INITIALISE",ERR,ERROR,*998)
 
    IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
      EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
      IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
        EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            IF(ALLOCATED(NONLINEAR_MATRICES%JACOBIANS)) THEN
              IF(ASSOCIATED(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR)) THEN
                CALL FlagError("Nonlinear matrices Jacobian is already associated.",ERR,ERROR,*998)
              ELSE
                ALLOCATE(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate equations Jacobian.",ERR,ERROR,*999)
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN_NUMBER=MATRIX_NUMBER
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%NONLINEAR_MATRICES=>NONLINEAR_MATRICES
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%NUMBER_OF_COLUMNS= &
                    & NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%UPDATE_JACOBIAN=.TRUE.
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%FIRST_ASSEMBLY=.TRUE.
                NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(MATRIX_NUMBER)%JACOBIAN=>NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR
                NULLIFY(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN)
                CALL EquationsMatrices_ElementMatrixInitialise(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR% &
                    & ELEMENT_JACOBIAN,ERR,ERROR,*999)
                CALL EquationsMatrices_NodalMatrixInitialise(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR% &
                    & NodalJacobian,ERR,ERROR,*999)
                NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR%JACOBIAN_CALCULATION_TYPE= &
                  & EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED
              ENDIF
            ELSE
              CALL FlagError("Equations matrices nonlinear matrieces Jacobian is not allocated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FlagError("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Nonlinear matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_JACOBIAN_INITIALISE")
    RETURN
999 CALL EQUATIONS_JACOBIAN_FINALISE(NONLINEAR_MATRICES%JACOBIANS(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_JACOBIAN_INITIALISE",ERR,ERROR)
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

    ENTERS("EQUATIONS_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have already been finished.",ERR,ERROR,*998)
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
                        CALL EquationsMatrix_StructureCalculate(EQUATIONS_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
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
                        & TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is not associated."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                      & " is not associated."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx                
              ELSE
                CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)                
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
                        CALL EquationsMatrix_StructureCalculate(EQUATIONS_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
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
                        & TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is not associated."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                      & " is not associated."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO !matrix_idx
              ELSE
                CALL FlagError("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)                
              ENDIF
            ENDIF
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
              !Nonlinear matrices
              NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
              IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                !Set up the Jacobian matrices
                DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
                  JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
                  IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                    COLUMN_DOMAIN_MAP=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(matrix_idx)%COLUMN_DOFS_MAPPING
                    IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
  !!TODO: Set the distributed matrix not to allocate the data if the Jacobian is not calculated.
                      !Create the distributed Jacobian matrix
                      CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(JACOBIAN_MATRIX%JACOBIAN,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(JACOBIAN_MATRIX%JACOBIAN,JACOBIAN_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                      !Calculate and set the matrix structure/sparsity pattern
                      IF(JACOBIAN_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                        & JACOBIAN_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                        CALL JacobianMatrix_StructureCalculate(JACOBIAN_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(JACOBIAN_MATRIX%JACOBIAN,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                        CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(JACOBIAN_MATRIX%JACOBIAN,ROW_INDICES,COLUMN_INDICES, &
                          & ERR,ERROR,*999)
                        IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                        IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                      ENDIF
                      CALL DISTRIBUTED_MATRIX_CREATE_FINISH(JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
                    ELSE
                      CALL FlagError("Column domain map is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Jacobian matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is not associated."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDDO
                !Set up the residual vector                
                CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,EQUATIONS_MATRICES%NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NONLINEAR_MATRICES%RESIDUAL,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_CREATE_FINISH(NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
                !Initialise the residual vector to zero for time dependent problems so that the previous residual is set to zero
                CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(NONLINEAR_MATRICES%RESIDUAL,0.0_DP,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
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
            CALL FlagError("Row domain map is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("EQUATIONS_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_CREATE_FINISH",ERR,ERROR)
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

    ENTERS("EQUATIONS_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN      
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
          CALL FlagError("Equations matrices is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(EQUATIONS_MATRICES)
          !Initialise the equations matrices
          CALL EQUATIONS_MATRICES_INITIALISE(EQUATIONS,ERR,ERROR,*999)
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        ENDIF
      ELSE
        CALL FlagError("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_CREATE_START")
    RETURN
999 CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS%EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_CREATE_START",ERR,ERROR)
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

    ENTERS("EQUATIONS_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Equations matrices is not associated",ERR,ERROR,*999)
    ENDIF
        
    EXITS("EQUATIONS_MATRICES_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_DESTROY",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE EQUATIONS_MATRICES_DESTROY

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the element matrix. Old CMISS name MELGE.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(ELEMENT_MATRIX,UPDATE_MATRIX,ROW_ELEMENT_NUMBERS,COLUMN_ELEMENT_NUMBERS, &
    & ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !<The element matrix to calculate
    LOGICAL :: UPDATE_MATRIX !<Is .TRUE. if the element matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: ROW_ELEMENT_NUMBERS(:) !<The row element number to calculate
    INTEGER(INTG), INTENT(IN) :: COLUMN_ELEMENT_NUMBERS(:) !<The column element number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROWS_FIELD_VARIABLE !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLS_FIELD_VARIABLE !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,global_ny,local_ny,node,node_idx,version,dataPointIdx, &
      & localDataPointNumber,elementIdx,rowElementNumber,colElementNumber
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(ROWS_FIELD_VARIABLE)) THEN
      IF(ASSOCIATED(COLS_FIELD_VARIABLE)) THEN
        ELEMENT_MATRIX%NUMBER_OF_ROWS=0
        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
        IF(UPDATE_MATRIX) THEN
          IF(ASSOCIATED(ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE)) THEN
            !Row and columns variable is the same.
            DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              DO elementIdx=1,SIZE(ROW_ELEMENT_NUMBERS)
                rowElementNumber=ROW_ELEMENT_NUMBERS(elementIdx)
                IF(rowElementNumber>=1.AND.rowElementNumber<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                      & ELEMENTS(rowElementNumber)
                    global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%BASIS
                    DO node_idx=1,BASIS%NUMBER_OF_NODES
                      node=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_NODES(node_idx)
                      DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                        derivative=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                        version=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%elementVersions(derivative_idx,node_idx)
                        local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                          & DERIVATIVES(derivative)%VERSIONS(version)
                        global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                        ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                        ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                      ENDDO !derivative_idx
                    ENDDO !node_idx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                    decompositionData=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
                    DO dataPointIdx=1,decompositionData%elementDataPoint(rowElementNumber)%numberOfProjectedData
                      localDataPointNumber=decompositionData%elementDataPoint(rowElementNumber)% &
                        & dataIndices(dataPointIdx)%localNumber
                      local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                        & DATA_POINTS(localDataPointNumber)
                      global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "// &
                      & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                      & " of rows field variable type "// &
                      & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)          
                  END SELECT
                ELSE
                  LOCAL_ERROR="Element number "//TRIM(NumberToVString(rowElementNumber,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                    & " of rows field variable type "// &
                    & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NumberToVString(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elementIdx
            ENDDO !component_idx
          ELSE
            !Row and column variables are different
            !Row mapping
            DO component_idx=1,ROWS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              DO elementIdx=1,SIZE(ROW_ELEMENT_NUMBERS)
                rowElementNumber=ROW_ELEMENT_NUMBERS(elementIdx)
                IF(rowElementNumber>=1.AND.rowElementNumber<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  SELECT CASE(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                      & ELEMENTS(rowElementNumber)
                    ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%BASIS
                    DO node_idx=1,BASIS%NUMBER_OF_NODES
                      node=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_NODES(node_idx)
                      DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                        derivative=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                        version=ELEMENTS_TOPOLOGY%ELEMENTS(rowElementNumber)%elementVersions(derivative_idx,node_idx)
                        local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                          & DERIVATIVES(derivative)%VERSIONS(version)
                        ELEMENT_MATRIX%NUMBER_OF_ROWS=ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                        ELEMENT_MATRIX%ROW_DOFS(ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                      ENDDO !derivative_idx
                    ENDDO !node_idx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                    decompositionData=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
                    DO dataPointIdx=1,decompositionData%elementDataPoint(colElementNumber)%numberOfProjectedData
                      localDataPointNumber=decompositionData%elementDataPoint(colElementNumber)% &
                        & dataIndices(dataPointIdx)%localNumber
                      local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                        & DATA_POINTS(localDataPointNumber)
                      global_ny=ROWS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "// &
                      & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                      & " of rows field variable type "// &
                      & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)          
                  END SELECT
                ELSE
                  LOCAL_ERROR="Row element number "//TRIM(NumberToVString(rowElementNumber,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                    & " of rows field variable type "// &
                    & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NumberToVString(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elementIdx
            ENDDO !component_idx
            !Column mapping
            DO component_idx=1,COLS_FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              ELEMENTS_TOPOLOGY=>COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              DO elementIdx=1,SIZE(COLUMN_ELEMENT_NUMBERS)
                colElementNumber=COLUMN_ELEMENT_NUMBERS(elementIdx)
                IF(colElementNumber>=1.AND.colElementNumber<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  SELECT CASE(COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                      & ELEMENTS(colElementNumber)
                    global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%BASIS
                    DO node_idx=1,BASIS%NUMBER_OF_NODES
                      node=ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%ELEMENT_NODES(node_idx)
                      DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                        derivative=ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                        version=ELEMENTS_TOPOLOGY%ELEMENTS(colElementNumber)%elementVersions(derivative_idx,node_idx)
                        local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                          & DERIVATIVES(derivative)%VERSIONS(version)
                        global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                        ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                      ENDDO !derivative_idx
                    ENDDO !node_idx
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                    decompositionData=>COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
                    DO dataPointIdx=1,decompositionData%elementDataPoint(colElementNumber)%numberOfProjectedData
                      localDataPointNumber=decompositionData%elementDataPoint(colElementNumber)% &
                        & dataIndices(dataPointIdx)%localNumber
                      local_ny=COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                        & DATA_POINTS(localDataPointNumber)
                      global_ny=COLS_FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      ELEMENT_MATRIX%NUMBER_OF_COLUMNS=ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      ELEMENT_MATRIX%COLUMN_DOFS(ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "// &
                      & TRIM(NumberToVString(COLS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                      & " of column field variable type "// &
                      & TRIM(NumberToVString(COLS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)          
                  END SELECT
                ELSE
                  LOCAL_ERROR="Column element number "//TRIM(NumberToVString(colElementNumber,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                    & " of column field variable type "// &
                    & TRIM(NumberToVString(COLS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NumberToVString(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !elementIdx
            ENDDO !component_idx
          ENDIF
          ELEMENT_MATRIX%MATRIX=0.0_DP
        ENDIF
      ELSE
        CALL FlagError("Columns field variable is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR,*999)

    ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
    IF(ALLOCATED(ELEMENT_MATRIX%ROW_DOFS)) DEALLOCATE(ELEMENT_MATRIX%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%COLUMN_DOFS)) DEALLOCATE(ELEMENT_MATRIX%COLUMN_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%MATRIX)) DEALLOCATE(ELEMENT_MATRIX%MATRIX)
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element matrix.
  SUBROUTINE EquationsMatrices_ElementMatrixInitialise(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !The element matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ElementMatrixInitialise",ERR,ERROR,*999)

    ELEMENT_MATRIX%EQUATIONS_MATRIX_NUMBER=0
    ELEMENT_MATRIX%NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
       
    EXITS("EquationsMatrices_ElementMatrixInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementMatrixInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EquationsMatrices_ElementMatrixInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the element matrix for the row and column field variables.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(elementMatrix,rowsFieldVariable,columnsFieldVariable, &
    & rowsNumberOfElements,colsNumberOfElements,err,error,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: elementMatrix !<The element matrix to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: columnsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(IN)  :: rowsNumberOfElements !<Number of elements in the row variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(IN)  :: colsNumberOfElements !<Number of elements in the col variables whose dofs are present in this element matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr, componentIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP",err,error,*998)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      IF(ASSOCIATED(columnsFieldVariable)) THEN
        elementMatrix%MAX_NUMBER_OF_ROWS = 0
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          elementMatrix%MAX_NUMBER_OF_ROWS=elementMatrix%MAX_NUMBER_OF_ROWS+ &
            & rowsFieldVariable%COMPONENTS(componentIdx)%maxNumberElementInterpolationParameters
        ENDDO
        elementMatrix%MAX_NUMBER_OF_ROWS=elementMatrix%MAX_NUMBER_OF_ROWS*rowsNumberOfElements
        elementMatrix%MAX_NUMBER_OF_COLUMNS = 0
        DO componentIdx=1,columnsFieldVariable%NUMBER_OF_COMPONENTS
          elementMatrix%MAX_NUMBER_OF_COLUMNS=elementMatrix%MAX_NUMBER_OF_COLUMNS+ &
            & columnsFieldVariable%COMPONENTS(componentIdx)%maxNumberElementInterpolationParameters
        ENDDO
        elementMatrix%MAX_NUMBER_OF_COLUMNS=elementMatrix%MAX_NUMBER_OF_COLUMNS*colsNumberOfElements
        IF(ALLOCATED(elementMatrix%ROW_DOFS)) THEN
          CALL FlagError("Element matrix row dofs already allocated.",err,error,*999)
        ELSE
          ALLOCATE(elementMatrix%ROW_DOFS(elementMatrix%MAX_NUMBER_OF_ROWS),STAT=err)
          IF(ERR/=0) CALL FlagError("Could not allocate element matrix row dofs.",err,error,*999)
        ENDIF
        IF(ALLOCATED(elementMatrix%COLUMN_DOFS)) THEN
          CALL FlagError("Element matrix column dofs already allocated.",err,error,*999)
        ELSE
          ALLOCATE(elementMatrix%COLUMN_DOFS(elementMatrix%MAX_NUMBER_OF_COLUMNS),STAT=err)
          IF(ERR/=0) CALL FlagError("Could not allocate element matrix column dofs.",err,error,*999)
        ENDIF
        IF(ALLOCATED(elementMatrix%MATRIX)) THEN
          CALL FlagError("Element matrix already allocated.",err,error,*999)
        ELSE
          ALLOCATE(elementMatrix%MATRIX(elementMatrix%MAX_NUMBER_OF_ROWS,elementMatrix%MAX_NUMBER_OF_COLUMNS),STAT=err)
          IF(ERR/=0) CALL FlagError("Could not allocate element matrix.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Columns field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP")
    RETURN
999 CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(elementMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP",err,error)
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
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,local_ny,node,node_idx,version,dataPointIdx,localDataPointNumber
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DecompositionDataPointsType), POINTER :: decompositionData
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE",ERR,ERROR,*999)

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
              local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                & ELEMENTS(ELEMENT_NUMBER)
              ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
              ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO node_idx=1,BASIS%NUMBER_OF_NODES
                node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                  derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                  version=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%elementVersions(derivative_idx,node_idx)
                  local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & DERIVATIVES(derivative)%VERSIONS(version)
                  ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                  ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
                ENDDO !derivative_idx
              ENDDO !node_idx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
              decompositionData=>ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%DECOMPOSITION%TOPOLOGY%dataPoints
              DO dataPointIdx=1,decompositionData%elementDataPoint(ELEMENT_NUMBER)%numberOfProjectedData
                localDataPointNumber=decompositionData%elementDataPoint(ELEMENT_NUMBER)% &
                  & dataIndices(dataPointIdx)%localNumber
                local_ny=ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP% &
                  & DATA_POINTS(localDataPointNumber)
                ELEMENT_VECTOR%NUMBER_OF_ROWS=ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                ELEMENT_VECTOR%ROW_DOFS(ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
              ENDDO
            CASE DEFAULT
              LOCAL_ERROR="The interpolation type of "// &
                & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                & " is invalid for component number "// &
                & TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
                & " of rows field variable type "// &
                & TRIM(NumberToVString(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))//"."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)          
            END SELECT
          ELSE
            LOCAL_ERROR="Element number "//TRIM(NumberToVString(ELEMENT_NUMBER,"*",ERR,ERROR))// &
              & " is invalid for component number "//TRIM(NumberToVString(component_idx,"*",ERR,ERROR))// &
              & " of rows field variable type "//TRIM(NumberToVString(ROWS_FIELD_VARIABLE%VARIABLE_TYPE,"*",ERR,ERROR))// &
              & ". The element number must be between 1 and "// &
              & TRIM(NumberToVString(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !component_idx
        ELEMENT_VECTOR%VECTOR=0.0_DP
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(ELEMENT_VECTOR%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_VECTOR%VECTOR)) DEALLOCATE(ELEMENT_VECTOR%VECTOR)
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR)

    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element vector
  SUBROUTINE EquationsMatrices_ElementVectorInitialise(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !The element vector to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_ElementVectorInitialise",ERR,ERROR,*999)

    ELEMENT_VECTOR%NUMBER_OF_ROWS=0
    ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
       
    EXITS("EquationsMatrices_ElementVectorInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_ElementVectorInitialise",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EquationsMatrices_ElementVectorInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the element vector for the row field variables.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(elementVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: elementVector !<The element vector to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,componentIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP",err,error,*998)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      elementVector%MAX_NUMBER_OF_ROWS = 0
      DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
        elementVector%MAX_NUMBER_OF_ROWS=elementVector%MAX_NUMBER_OF_ROWS+ &
          & rowsFieldVariable%COMPONENTS(componentIdx)%maxNumberElementInterpolationParameters
      ENDDO
      IF(ALLOCATED(elementVector%ROW_DOFS)) THEN
        CALL FlagError("Element vector row dofs is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(elementVector%ROW_DOFS(elementVector%MAX_NUMBER_OF_ROWS),STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate element vector row dofs.",err,error,*999)
      ENDIF
      IF(ALLOCATED(elementVector%VECTOR)) THEN
        CALL FlagError("Element vector vector already allocated.",err,error,*999)
      ELSE
        ALLOCATE(elementVector%VECTOR(elementVector%MAX_NUMBER_OF_ROWS),STAT=err)
        IF(ERR/=0) CALL FlagError("Could not allocate element vector vector.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP")
    RETURN
999 CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(elementVector,DUMMY_ERR,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP",err,error)
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

    ENTERS("EQUATIONS_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

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
            LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
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
            LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
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
      CALL FlagError("Equations matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_ELEMENT_ADD()")
#endif
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_ADD")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_ADD",ERR,ERROR)
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,COL_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_ELEMENT_CALCULATE()")
#endif

    ENTERS("EQUATIONS_MATRICES_ELEMENT_CALCULATE",ERR,ERROR,*999)

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
                  & [ELEMENT_NUMBER],[ELEMENT_NUMBER],FIELD_VARIABLE,FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
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
                  & [ELEMENT_NUMBER],[ELEMENT_NUMBER],FIELD_VARIABLE,FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FlagError("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          !Calculate the rows and columns of the Jacobian
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(1)%VARIABLE !Row field variable
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                COL_FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(matrix_idx)%VARIABLE
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN,JACOBIAN_MATRIX%UPDATE_JACOBIAN, &
                  & [ELEMENT_NUMBER],[ELEMENT_NUMBER],FIELD_VARIABLE,COL_FIELD_VARIABLE,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Jacobian matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
            !Calculate the rows of the equations residual
            RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
            IF(ASSOCIATED(RHS_MAPPING)) THEN
              FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            ELSE
              FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
            ENDIF
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,NONLINEAR_MATRICES% &
              & UPDATE_RESIDUAL,ELEMENT_NUMBER,FIELD_VARIABLE,ERR,ERROR,*999)
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL_CALCULATED=0
          ELSE
            CALL FlagError("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
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
            CALL FlagError("Equations mapping rhs mapping is not associated.",ERR,ERROR,*999)
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
            CALL FlagError("Equations mapping rhs mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_ELEMENT_CALCULATE()")
#endif
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices and rhs of the nodal matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE EquationsMatrices_NodalCalculate(equationsMatrices,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The nodal number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: sourceVector
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,columnFieldVariable
    TYPE(VARYING_STRING) :: localError

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_NodalCalculate()")
#endif

    ENTERS("EquationsMatrices_NodalCalculate",err,error,*999)

    IF(ASSOCIATED(equationsMatrices)) THEN
      equationsMapping=>equationsMatrices%EQUATIONS_MAPPING
      IF(ASSOCIATED(equationsMapping)) THEN
        dynamicMatrices=>equationsMatrices%DYNAMIC_MATRICES
        IF(ASSOCIATED(dynamicMatrices)) THEN
          !Calculate the row and columns for the dynamic equations matrices
          dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
          IF(ASSOCIATED(dynamicMapping)) THEN
            DO matrixIdx=1,dynamicMatrices%NUMBER_OF_DYNAMIC_MATRICES
              equationsMatrix=>dynamicMatrices%MATRICES(matrixIdx)%PTR
              IF(ASSOCIATED(equationsMatrix)) THEN
                fieldVariable=>dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE
                CALL EquationsMatrices_NodalMatrixCalculate(equationsMatrix%NodalMatrix,equationsMatrix%UPDATE_MATRIX, &
                  & nodeNumber,nodeNumber,fieldVariable,fieldVariable,err,error,*999)
              ELSE
                localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
                  & " is not associated."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          ELSE
            CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        linearMatrices=>equationsMatrices%LINEAR_MATRICES
        IF(ASSOCIATED(linearMatrices)) THEN
          !Calculate the row and columns for the linear equations matrices
          linearMapping=>equationsMapping%LINEAR_MAPPING
          IF(ASSOCIATED(linearMapping)) THEN
            DO matrixIdx=1,linearMatrices%NUMBER_OF_LINEAR_MATRICES
              equationsMatrix=>linearMatrices%MATRICES(matrixIdx)%PTR
              IF(ASSOCIATED(equationsMatrix)) THEN
                fieldVariable=>linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE
                CALL EquationsMatrices_NodalMatrixCalculate(equationsMatrix%NodalMatrix,equationsMatrix%UPDATE_MATRIX, &
                  & nodeNumber,nodeNumber,fieldVariable,fieldVariable,err,error,*999)
              ELSE
                localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
                  & " is not associated."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          ELSE
            CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
        IF(ASSOCIATED(nonlinearMatrices)) THEN
          !Calculate the rows and columns of the Jacobian
          nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
          IF(ASSOCIATED(nonlinearMapping)) THEN
            fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE !Row field variable
            DO matrixIdx=1,nonlinearMatrices%NUMBER_OF_JACOBIANS
              jacobianMatrix=>nonlinearMatrices%JACOBIANS(matrixIdx)%PTR
              IF(ASSOCIATED(jacobianMatrix)) THEN
                columnFieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE
                CALL EquationsMatrices_NodalMatrixCalculate(jacobianMatrix%NodalJacobian,jacobianMatrix%UPDATE_JACOBIAN, &
                  & nodeNumber,nodeNumber,fieldVariable,columnFieldVariable,err,error,*999)
              ELSE
                CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
              ENDIF
            ENDDO
            !Calculate the rows of the equations residual
            rhsMapping=>equationsMapping%RHS_MAPPING
            IF(ASSOCIATED(rhsMapping)) THEN
              fieldVariable=>rhsMapping%RHS_VARIABLE
            ELSE
              fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
            ENDIF
            CALL EquationsMatrices_NodalVectorCalculate(nonlinearMatrices%NodalResidual,nonlinearMatrices% &
              & UPDATE_RESIDUAL,nodeNumber,fieldVariable,err,error,*999)
            nonlinearMatrices%NodalResidualCalculated=0
          ELSE
            CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        rhsVector=>equationsMatrices%RHS_VECTOR
        IF(ASSOCIATED(rhsVector)) THEN
          rhsMapping=>equationsMapping%RHS_MAPPING
          IF(ASSOCIATED(rhsMapping)) THEN
            !Calculate the rows  for the equations RHS
            fieldVariable=>rhsMapping%RHS_VARIABLE
            CALL EquationsMatrices_NodalVectorCalculate(rhsVector%NodalVector,rhsVector%UPDATE_VECTOR,nodeNumber, &
              & fieldVariable,err,error,*999)
          ELSE
            CALL FlagError("Equations mapping rhs mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        sourceVector=>equationsMatrices%SOURCE_VECTOR
        IF(ASSOCIATED(sourceVector)) THEN
          !Calculate the rows the equations source. The number of rows is not set by the source field so take the number of rows
          !from the RHS vector in the first instance.
          rhsMapping=>equationsMapping%RHS_MAPPING
          IF(ASSOCIATED(rhsMapping)) THEN
            fieldVariable=>rhsMapping%RHS_VARIABLE
            CALL EquationsMatrices_NodalVectorCalculate(sourceVector%NodalVector,sourceVector%UPDATE_VECTOR, &
              & nodeNumber,fieldVariable,err,error,*999)
          ELSE
            CALL FlagError("Equations mapping rhs mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations mapping is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not allocated",err,error,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_NodalCalculate()")
#endif
    
    EXITS("EquationsMatrices_NodalCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalCalculate",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations matrices of the nodal matrix.
  SUBROUTINE EquationsMatrices_NodalMatrixCalculate(nodalMatrix,updateMatrix,rowNodeNumber,columnNodeNumber, &
    & rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix to calculate
    LOGICAL :: updateMatrix !<Is .TRUE. if the nodal matrix is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowNodeNumber !<The row nodal number to calculate
    INTEGER(INTG), INTENT(IN) :: columnNodeNumber !<The column nodal number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx
    INTEGER(INTG) :: localRow,globalRow,localColumn,globalColumn
    INTEGER(INTG) :: numberOfDerivatives,numberOfVersions,versionIdx,derivativeIdx
    TYPE(DOMAIN_NODES_TYPE), POINTER :: nodesTopology
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_NodalMatrixCalculate",err,error,*999)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      IF(ASSOCIATED(colsFieldVariable)) THEN
        nodalMatrix%numberOfRows=0
        nodalMatrix%numberOfColumns=0
        IF(updateMatrix) THEN
          IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
            !Row and columns variable is the same.
            DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
              nodesTopology=>rowsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES
              IF(rowNodeNumber>=1.AND.rowNodeNumber<=nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
                SELECT CASE(rowsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  globalRow=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localRow)
                  nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                  nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                  nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                  nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalRow
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                    & ELEMENTS(rowNodeNumber)
                  globalRow=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localRow)
                  nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                  nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                  nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                  nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalRow
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  numberOfDerivatives=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
                    & NUMBER_OF_DERIVATIVES
                  DO derivativeIdx=1,numberOfDerivatives
                    numberOfVersions=rowsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES%NODES(rowNodeNumber)% &
                      & DERIVATIVES(derivativeIdx)%numberOfVersions
                    DO versionIdx=1,numberOfVersions
                      localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                        & NODES(rowNodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                      globalRow=rowsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localRow)
                      nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                      nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                      nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                      nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalRow
                    ENDDO !versionIdx
                  ENDDO !derivativeIdx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The interpolation type of "// &
                    & TRIM(NumberToVString(rowsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                    & " is invalid for component number "// &
                    & TRIM(NumberToVString(componentIdx,"*",err,error))// &
                    & " of rows field variable type "// &
                    & TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)          
                END SELECT
              ELSE
                localError="Nodal number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
                  & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                  & " of rows field variable type "// &
                  & TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
                  & ". The nodal number must be between 1 and "// &
                  & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !componentIdx
          ELSE
            !Row and column variables are different
            !Row mapping
            DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
              nodesTopology=>rowsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES
              IF(rowNodeNumber>=1.AND.rowNodeNumber<=nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
                SELECT CASE(rowsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                  nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                    & ELEMENTS(rowNodeNumber)
                  nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                  nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  numberOfDerivatives=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
                    & NUMBER_OF_DERIVATIVES
                  DO derivativeIdx=1,numberOfDerivatives
                    numberOfVersions=colsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES%NODES(rowNodeNumber)% &
                      & DERIVATIVES(derivativeIdx)%numberOfVersions
                    DO versionIdx=1,numberOfVersions
                      localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                        & NODES(rowNodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                      nodalMatrix%numberOfRows=nodalMatrix%numberOfRows+1
                      nodalMatrix%rowDofs(nodalMatrix%numberOfRows)=localRow
                    ENDDO !versionIdx
                  ENDDO !derivativeIdx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The interpolation type of "// &
                    & TRIM(NumberToVString(rowsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                    & " is invalid for component number "// &
                    & TRIM(NumberToVString(componentIdx,"*",err,error))// &
                    & " of rows field variable type "// &
                    & TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)          
                END SELECT
              ELSE
                localError="Row nodal number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
                  & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                  & " of rows field variable type "// &
                  & TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
                  & ". The nodal number must be between 1 and "// &
                  & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !componentIdx
            !Column mapping
            DO componentIdx=1,colsFieldVariable%NUMBER_OF_COMPONENTS
              nodesTopology=>colsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES
              IF(columnNodeNumber>=1.AND.columnNodeNumber<=nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
                SELECT CASE(colsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  localColumn=colsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  globalColumn=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                  nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                  nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  localColumn=colsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                    & ELEMENTS(columnNodeNumber)
                  globalColumn=colsFieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                  nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                  nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=globalColumn
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  numberOfDerivatives=colsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
                    & NUMBER_OF_DERIVATIVES
                  DO derivativeIdx=1,numberOfDerivatives
                    numberOfVersions=colsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES%NODES(rowNodeNumber)% &
                      & DERIVATIVES(derivativeIdx)%numberOfVersions
                    DO versionIdx=1,numberOfVersions
                      localRow=colsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                        & NODES(rowNodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                      nodalMatrix%numberOfColumns=nodalMatrix%numberOfColumns+1
                      nodalMatrix%columnDofs(nodalMatrix%numberOfColumns)=localRow
                    ENDDO !versionIdx
                  ENDDO !derivativeIdx
                CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The interpolation type of "// &
                    & TRIM(NumberToVString(colsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                    & " is invalid for component number "// &
                    & TRIM(NumberToVString(componentIdx,"*",err,error))// &
                    & " of column field variable type "// &
                    & TRIM(NumberToVString(colsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)          
                END SELECT
              ELSE
                localError="Column nodal number "//TRIM(NumberToVString(columnNodeNumber,"*",err,error))// &
                  & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
                  & " of column field variable type "// &
                  & TRIM(NumberToVString(colsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
                  & ". The nodal number must be between 1 and "// &
                  & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !componentIdx
          ENDIF
          nodalMatrix%matrix=0.0_DP
        ENDIF
      ELSE
        CALL FlagError("Columns field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NodalMatrixCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalMatrixCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the equations rhs of the nodal rhs vector.
  SUBROUTINE EquationsMatrices_NodalVectorCalculate(nodalVector,updateVector,rowNodeNumber,rowsFieldVariable, &
    & err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !<The nodal vector to calculate.
    LOGICAL :: updateVector !<Is .TRUE. if the nodal vector is to be updated, .FALSE. if not.
    INTEGER(INTG), INTENT(IN) :: rowNodeNumber !<The nodal number to calculate
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,localRow
    INTEGER(INTG) :: numberOfDerivatives,numberOfVersions,versionIdx,derivativeIdx
    TYPE(DOMAIN_NODES_TYPE), POINTER :: nodesTopology
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_NodalVectorCalculate",err,error,*999)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      !Calculate the rows for the nodal vector
      nodalVector%numberOfRows=0
      IF(updateVector) THEN
        DO componentIdx=1,rowsFieldVariable%NUMBER_OF_COMPONENTS
          nodesTopology=>rowsFieldVariable%components(componentIdx)%domain%topology%nodes
          IF(rowNodeNumber>=1.AND.rowNodeNumber<=nodesTopology%TOTAL_NUMBER_OF_NODES) THEN
            SELECT CASE(rowsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
              nodalVector%numberOfRows=nodalVector%numberOfRows+1
              nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP% &
                & ELEMENTS(rowNodeNumber)
              nodalVector%numberOfRows=nodalVector%numberOfRows+1
              nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              numberOfDerivatives=rowsFieldVariable%components(componentIdx)%domain%topology%nodes%nodes(rowNodeNumber)% &
                & NUMBER_OF_DERIVATIVES
              DO derivativeIdx=1,numberOfDerivatives
                numberOfVersions=rowsFieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES%NODES(rowNodeNumber)% &
                  & DERIVATIVES(derivativeIdx)%numberOfVersions
                DO versionIdx=1,numberOfVersions
                  localRow=rowsFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
                    & nodes(rowNodeNumber)%derivatives(derivativeIdx)%versions(versionIdx)
                  nodalVector%numberOfRows=nodalVector%numberOfRows+1
                  nodalVector%rowDofs(nodalVector%numberOfRows)=localRow
                ENDDO !versionIdx
              ENDDO !derivativeIdx
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The interpolation type of "// &
                & TRIM(NumberToVString(rowsFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE,"*",err,error))// &
                & " is invalid for component number "// &
                & TRIM(NumberToVString(componentIdx,"*",err,error))// &
                & " of rows field variable type "// &
                & TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)          
            END SELECT
          ELSE
            localError="Node number "//TRIM(NumberToVString(rowNodeNumber,"*",err,error))// &
              & " is invalid for component number "//TRIM(NumberToVString(componentIdx,"*",err,error))// &
              & " of rows field variable type "//TRIM(NumberToVString(rowsFieldVariable%VARIABLE_TYPE,"*",err,error))// &
              & ". The nodal number must be between 1 and "// &
              & TRIM(NumberToVString(nodesTopology%TOTAL_NUMBER_OF_NODES,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !componentIdx
        nodalVector%vector=0.0_DP
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NodalVectorCalculate")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NodalVectorCalculate

  !
  !================================================================================================================================
  !

  !>Adds the nodal matrices and rhs vector into the equations matrices and rhs vector.
  SUBROUTINE EquationsMatrices_NodeAdd(equationsMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnIdx,matrixIdx,rowIdx
    REAL(DP) :: sum
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: sourceVector
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_NodeAdd()")
#endif

    ENTERS("EquationsMatrices_NodeAdd",err,error,*999)

    IF(ASSOCIATED(equationsMatrices)) THEN
      dynamicMatrices=>equationsMatrices%DYNAMIC_MATRICES
      IF(ASSOCIATED(dynamicMatrices)) THEN
        !Add the nodal matrices
        DO matrixIdx=1,dynamicMatrices%NUMBER_OF_DYNAMIC_MATRICES
          equationsMatrix=>dynamicMatrices%MATRICES(matrixIdx)%PTR
          IF(ASSOCIATED(equationsMatrix)) THEN
            IF(equationsMatrix%UPDATE_MATRIX) THEN
              !Handle lumped matrices
              IF(equationsMatrix%LUMPED) THEN
                DO rowIdx=1,equationsMatrix%NodalMatrix%numberOfRows
                  sum=0.0_DP
                  DO columnIdx=1,equationsMatrix%NodalMatrix%numberOfColumns
                    sum=sum+equationsMatrix%NodalMatrix%matrix(rowIdx,columnIdx)
                    equationsMatrix%NodalMatrix%matrix(rowIdx,columnIdx)=0.0_DP
                  ENDDO !columnIdx
                  equationsMatrix%NodalMatrix%matrix(rowIdx,rowIdx)=sum
                  !Add the nodal matrice into the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(equationsMatrix%matrix,equationsMatrix%NodalMatrix%rowDofs(rowIdx), &
                    & equationsMatrix%NodalMatrix%columnDofs(rowIdx),equationsMatrix%NodalMatrix%matrix(rowIdx, &
                    & rowIdx),err,error,*999)
                ENDDO !rowIdx
              ELSE
                !Add the nodal matrice into the distributed equations matrix
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(equationsMatrix%matrix,equationsMatrix%NodalMatrix%rowDofs(1: &
                  & equationsMatrix%NodalMatrix%numberOfRows),equationsMatrix%NodalMatrix%columnDofs(1: &
                  & equationsMatrix%NodalMatrix%numberOfColumns),equationsMatrix%NodalMatrix%matrix(1: &
                  & equationsMatrix%NodalMatrix%numberOfRows,1:equationsMatrix%NodalMatrix%numberOfColumns), &
                  & err,error,*999)
              ENDIF
            ENDIF
          ELSE
            localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
              & " is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      linearMatrices=>equationsMatrices%LINEAR_MATRICES
      IF(ASSOCIATED(linearMatrices)) THEN
        !Add the nodal matrices
        DO matrixIdx=1,linearMatrices%NUMBER_OF_LINEAR_MATRICES
          equationsMatrix=>linearMatrices%MATRICES(matrixIdx)%PTR
          IF(ASSOCIATED(equationsMatrix)) THEN
            IF(equationsMatrix%UPDATE_MATRIX) THEN
              !Handle lumped matrices
              IF(equationsMatrix%LUMPED) THEN
                DO rowIdx=1,equationsMatrix%NodalMatrix%numberOfRows
                  sum=0.0_DP
                  DO columnIdx=1,equationsMatrix%NodalMatrix%numberOfColumns
                    sum=sum+equationsMatrix%NodalMatrix%matrix(rowIdx,columnIdx)
                    equationsMatrix%NodalMatrix%matrix(rowIdx,columnIdx)=0.0_DP
                  ENDDO !columnIdx
                  equationsMatrix%NodalMatrix%matrix(rowIdx,rowIdx)=sum
                  !Add the nodal matrice into the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_VALUES_ADD(equationsMatrix%matrix,equationsMatrix%NodalMatrix%rowDofs(rowIdx), &
                    & equationsMatrix%NodalMatrix%columnDofs(rowIdx),equationsMatrix%NodalMatrix%matrix(rowIdx, &
                    & rowIdx),err,error,*999)
                ENDDO !rowIdx
              ELSE
                !Add the nodal matrice into the distributed equations matrix
                CALL DISTRIBUTED_MATRIX_VALUES_ADD(equationsMatrix%matrix,equationsMatrix%NodalMatrix%rowDofs(1: &
                  & equationsMatrix%NodalMatrix%numberOfRows),equationsMatrix%NodalMatrix%columnDofs(1: &
                  & equationsMatrix%NodalMatrix%numberOfColumns),equationsMatrix%NodalMatrix%matrix(1: &
                  & equationsMatrix%NodalMatrix%numberOfRows,1:equationsMatrix%NodalMatrix%numberOfColumns), &
                  & err,error,*999)
              ENDIF
            ENDIF
          ELSE
            localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
              & " is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        IF(nonlinearMatrices%UPDATE_RESIDUAL) THEN
          !Add the residual nodal vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(nonlinearMatrices%RESIDUAL,nonlinearMatrices%NodalResidual%rowDofs(1: &
            & nonlinearMatrices%NodalResidual%numberOfRows),nonlinearMatrices%NodalResidual%vector(1:nonlinearMatrices% &
            & NodalResidual%numberOfRows),err,error,*999)
        ENDIF
      ENDIF
      rhsVector=>equationsMatrices%RHS_VECTOR
      IF(ASSOCIATED(rhsVector)) THEN
        IF(rhsVector%UPDATE_VECTOR) THEN
          !Add the rhs nodal vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(rhsVector%vector,rhsVector%NodalVector%rowDofs(1: &
            & rhsVector%NodalVector%numberOfRows),rhsVector%NodalVector%vector(1:rhsVector% &
            & NodalVector%numberOfRows),err,error,*999)
        ENDIF
      ENDIF
      sourceVector=>equationsMatrices%SOURCE_VECTOR
      IF(ASSOCIATED(sourceVector)) THEN
        IF(sourceVector%UPDATE_VECTOR) THEN
          !Add the rhs nodal vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(sourceVector%vector,sourceVector%NodalVector%rowDofs(1: &
            & sourceVector%NodalVector%numberOfRows),sourceVector%NodalVector%vector(1:sourceVector% &
            & NodalVector%numberOfRows),err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not allocated.",err,error,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_NodeAdd()")
#endif
    
    EXITS("EquationsMatrices_NodeAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodeAdd",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodeAdd

  !
  !================================================================================================================================
  !

  !>Initialise the nodal calculation information for the equations matrices
  SUBROUTINE EquationsMatrices_NodalInitialise(equationsMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices !The equations matrices to initialise the nodal information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: sourceVector
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,columnFieldVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_NodalInitialise",err,error,*999)

    IF(ASSOCIATED(equationsMatrices)) THEN
      equationsMapping=>equationsMatrices%EQUATIONS_MAPPING
      IF(ASSOCIATED(equationsMapping)) THEN
        dynamicMatrices=>equationsMatrices%DYNAMIC_MATRICES
        IF(ASSOCIATED(dynamicMatrices)) THEN
          !Initialise the dynamic nodal matrices
          dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
          IF(ASSOCIATED(dynamicMapping)) THEN
            DO matrixIdx=1,dynamicMatrices%NUMBER_OF_DYNAMIC_MATRICES
              equationsMatrix=>dynamicMatrices%MATRICES(matrixIdx)%PTR
              IF(ASSOCIATED(equationsMatrix)) THEN
                fieldVariable=>dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE
                CALL EquationsMatrices_NodalMatrixSetup(equationsMatrix%NodalMatrix,fieldVariable,fieldVariable, &
                  & err,error,*999)
              ELSE
                localError="Equations dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
                  & " is not associated."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          ELSE
            CALL FlagError("Equations mapping dynamic mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        linearMatrices=>equationsMatrices%LINEAR_MATRICES
        IF(ASSOCIATED(linearMatrices)) THEN
          !Initialise the linear nodal matrices
          linearMapping=>equationsMapping%LINEAR_MAPPING
          IF(ASSOCIATED(linearMapping)) THEN
            DO matrixIdx=1,linearMatrices%NUMBER_OF_LINEAR_MATRICES
              equationsMatrix=>linearMatrices%MATRICES(matrixIdx)%PTR
              IF(ASSOCIATED(equationsMatrix)) THEN
                fieldVariable=>linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixIdx)%VARIABLE
                CALL EquationsMatrices_NodalMatrixSetup(equationsMatrix%NodalMatrix,fieldVariable,fieldVariable, &
                  & err,error,*999)
              ELSE
                localError="Equations linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
                  & " is not associated."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !matrixIdx
          ELSE
            CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
        IF(ASSOCIATED(nonlinearMatrices)) THEN
          !Initialise the Jacobian nodal matrices
          nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
          IF(ASSOCIATED(nonlinearMapping)) THEN
            fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
            DO matrixIdx=1,nonlinearMatrices%NUMBER_OF_JACOBIANS
              jacobianMatrix=>nonlinearMatrices%JACOBIANS(matrixIdx)%PTR
              IF(ASSOCIATED(jacobianMatrix)) THEN
                columnFieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixIdx)%VARIABLE
                CALL EquationsMatrices_NodalMatrixSetup(jacobianMatrix%NodalJacobian,fieldVariable,columnFieldVariable, &
                  & err,error,*999)
              ELSE
                CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
              ENDIF
            ENDDO
            !Use RHS variable for residual vector, otherwise first nonlinear variable if no RHS
            rhsMapping=>equationsMapping%RHS_MAPPING
            IF(ASSOCIATED(rhsMapping)) THEN
              fieldVariable=>rhsMapping%RHS_VARIABLE
            ELSE
              fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
            ENDIF
            CALL EquationsMatrices_NodalVectorSetup(nonlinearMatrices%NodalResidual,fieldVariable,err,error,*999)
            nonlinearMatrices%NodalResidualCalculated=0
          ELSE
            CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        rhsVector=>equationsMatrices%RHS_VECTOR
        IF(ASSOCIATED(rhsVector)) THEN
          !Initialise the RHS nodal vector
          rhsMapping=>equationsMapping%RHS_MAPPING
          IF(ASSOCIATED(rhsMapping)) THEN
            fieldVariable=>rhsMapping%RHS_VARIABLE
            CALL EquationsMatrices_NodalVectorSetup(rhsVector%NodalVector,fieldVariable,err,error,*999)
          ELSE
            CALL FlagError("RHS mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
        sourceVector=>equationsMatrices%SOURCE_VECTOR
        IF(ASSOCIATED(sourceVector)) THEN
          !Initialise the source nodal vector. Note that the number of rows in the source vector is taken, for now, from the RHS
          !vector
          IF(ASSOCIATED(rhsVector)) THEN
            !Initialise the RHS nodal vector
            rhsMapping=>equationsMapping%RHS_MAPPING
            IF(ASSOCIATED(rhsMapping)) THEN
              fieldVariable=>rhsMapping%RHS_VARIABLE
              CALL EquationsMatrices_NodalVectorSetup(sourceVector%NodalVector,fieldVariable,err,error,*999)
            ELSE
              CALL FlagError("RHS mapping is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices mapping is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NodalInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalInitialise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the nodal matrix for the row and column field variables.
  SUBROUTINE EquationsMatrices_NodalMatrixSetup(nodalMatrix,rowsFieldVariable,colsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !<The nodal matrix to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable !<A pointer to the field variable associated with the columns
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_NodalMatrixSetup",err,error,*998)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      IF(ASSOCIATED(colsFieldVariable)) THEN
        nodalMatrix%maxNumberOfRows=rowsFieldVariable%maxNumberNodeInterpolationParameters* &
          & rowsFieldVariable%NUMBER_OF_COMPONENTS
        nodalMatrix%maxNumberOfColumns=colsFieldVariable%maxNumberNodeInterpolationParameters* &
          & colsFieldVariable%NUMBER_OF_COMPONENTS
        IF(ALLOCATED(nodalMatrix%rowDofs)) THEN
          CALL FlagError("Nodal matrix row dofs already allocated.",err,error,*999)
        ELSE
          ALLOCATE(nodalMatrix%rowDofs(nodalMatrix%maxNumberOfRows),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate nodal matrix row dofs.",err,error,*999)
        ENDIF
        IF(ALLOCATED(nodalMatrix%columnDofs)) THEN
          CALL FlagError("Nodal matrix column dofs already allocated.",err,error,*999)
        ELSE
          ALLOCATE(nodalMatrix%columnDofs(nodalMatrix%maxNumberOfColumns),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate nodal matrix column dofs.",err,error,*999)
        ENDIF
        IF(ALLOCATED(nodalMatrix%matrix)) THEN
          CALL FlagError("Nodal matrix already allocated.",err,error,*999)
        ELSE
          ALLOCATE(nodalMatrix%matrix(nodalMatrix%maxNumberOfRows,nodalMatrix%maxNumberOfColumns),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate nodal matrix.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Columns field variable is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NodalMatrixSetup")
    RETURN
999 CALL EquationsMatrices_NodalMatrixFinalise(nodalMatrix,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NodalMatrixSetup",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalMatrixSetup

  !
  !================================================================================================================================
  !

  !>Sets up the nodal vector for the row field variables.
  SUBROUTINE EquationsMatrices_NodalVectorSetup(nodalVector,rowsFieldVariable,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !<The nodal vector to setup
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowsFieldVariable !<A pointer to the field variable associated with the rows
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsMatrices_NodalVectorSetup",err,error,*998)

    IF(ASSOCIATED(rowsFieldVariable)) THEN
      nodalVector%maxNumberOfRows=rowsFieldVariable%maxNumberNodeInterpolationParameters* &
        & rowsFieldVariable%NUMBER_OF_COMPONENTS
      IF(ALLOCATED(nodalVector%rowDofs)) THEN
        CALL FlagError("Nodal vector row dofs is already allocated.",err,error,*999)        
      ELSE
        ALLOCATE(nodalVector%rowDofs(nodalVector%maxNumberOfRows),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate nodal vector row dofs.",err,error,*999)
      ENDIF
      IF(ALLOCATED(nodalVector%vector)) THEN
        CALL FlagError("Nodal vector vector already allocated.",err,error,*999)        
      ELSE
        ALLOCATE(nodalVector%vector(nodalVector%maxNumberOfRows),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate nodal vector vector.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Rows field variable is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NodalVectorSetup")
    RETURN
999 CALL EquationsMatrices_NodalVectorFinalise(nodalVector,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsMatrices_NodalVectorSetup",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalVectorSetup

  !
  !================================================================================================================================
  !

  !>Finalise the nodal calculation information and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalFinalise(equationsMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices !<The equations matrices for which to finalise the nodals
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: sourceVector
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsMatrices_NodalFinalise",err,error,*999)

    IF(ASSOCIATED(equationsMatrices)) THEN
      dynamicMatrices=>equationsMatrices%DYNAMIC_MATRICES
      IF(ASSOCIATED(dynamicMatrices)) THEN
        !Finalise the dynamic nodal matrices
        DO matrixIdx=1,dynamicMatrices%NUMBER_OF_DYNAMIC_MATRICES
          equationsMatrix=>dynamicMatrices%MATRICES(matrixIdx)%PTR
          IF(ASSOCIATED(equationsMatrix)) THEN
            equationsMatrix%NodalMatrix%maxNumberOfRows=0
            equationsMatrix%NodalMatrix%maxNumberOfColumns=0
            IF(ALLOCATED(equationsMatrix%NodalMatrix%rowDofs)) DEALLOCATE(equationsMatrix%NodalMatrix%rowDofs)
            IF(ALLOCATED(equationsMatrix%NodalMatrix%columnDofs)) DEALLOCATE(equationsMatrix%NodalMatrix%columnDofs)
            IF(ALLOCATED(equationsMatrix%NodalMatrix%matrix)) DEALLOCATE(equationsMatrix%NodalMatrix%matrix)
          ELSE
            localError="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
              & " is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      linearMatrices=>equationsMatrices%LINEAR_MATRICES
      IF(ASSOCIATED(linearMatrices)) THEN
        !Finalise the linear nodal matrices
        DO matrixIdx=1,linearMatrices%NUMBER_OF_LINEAR_MATRICES
          equationsMatrix=>linearMatrices%MATRICES(matrixIdx)%PTR
          IF(ASSOCIATED(equationsMatrix)) THEN
            equationsMatrix%NodalMatrix%maxNumberOfRows=0
            equationsMatrix%NodalMatrix%maxNumberOfColumns=0
            IF(ALLOCATED(equationsMatrix%NodalMatrix%rowDofs)) DEALLOCATE(equationsMatrix%NodalMatrix%rowDofs)
            IF(ALLOCATED(equationsMatrix%NodalMatrix%columnDofs)) DEALLOCATE(equationsMatrix%NodalMatrix%columnDofs)
            IF(ALLOCATED(equationsMatrix%NodalMatrix%matrix)) DEALLOCATE(equationsMatrix%NodalMatrix%matrix)
          ELSE
            localError="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
              & " is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        DO matrixIdx=1,nonlinearMatrices%NUMBER_OF_JACOBIANS
          jacobianMatrix=>nonlinearMatrices%JACOBIANS(matrixIdx)%PTR
          IF(ASSOCIATED(jacobianMatrix)) THEN
            jacobianMatrix%NodalJacobian%maxNumberOfRows=0
            jacobianMatrix%NodalJacobian%maxNumberOfColumns=0
            IF(ALLOCATED(jacobianMatrix%NodalJacobian%rowDofs)) DEALLOCATE(jacobianMatrix%NodalJacobian%rowDofs)
            IF(ALLOCATED(jacobianMatrix%NodalJacobian%columnDofs)) DEALLOCATE(jacobianMatrix%NodalJacobian%columnDofs)
            IF(ALLOCATED(jacobianMatrix%NodalJacobian%matrix)) DEALLOCATE(jacobianMatrix%NodalJacobian%matrix)
          ELSE
            CALL FlagError("Nonlinear matrices Jacobian number "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
                & " is not associated.",err,error,*999)
          ENDIF
        ENDDO
        nonlinearMatrices%NodalResidual%maxNumberOfRows=0
        IF(ALLOCATED(nonlinearMatrices%NodalResidual%rowDofs)) DEALLOCATE(nonlinearMatrices%NodalResidual%rowDofs)
        IF(ALLOCATED(nonlinearMatrices%NodalResidual%vector)) DEALLOCATE(nonlinearMatrices%NodalResidual%vector)
      ENDIF
      rhsVector=>equationsMatrices%RHS_VECTOR
      IF(ASSOCIATED(rhsVector)) THEN
        !Finalise the nodal vector
        rhsVector%NodalVector%maxNumberOfRows=0
        IF(ALLOCATED(rhsVector%NodalVector%rowDofs)) DEALLOCATE(rhsVector%NodalVector%rowDofs)
        IF(ALLOCATED(rhsVector%NodalVector%vector)) DEALLOCATE(rhsVector%NodalVector%vector)
      ENDIF
      sourceVector=>equationsMatrices%SOURCE_VECTOR
      IF(ASSOCIATED(sourceVector)) THEN
        !Finalise the nodal source vector
        sourceVector%NodalVector%maxNumberOfRows=0
        IF(ALLOCATED(sourceVector%NodalVector%rowDofs)) DEALLOCATE(sourceVector%NodalVector%rowDofs)
        IF(ALLOCATED(sourceVector%NodalVector%vector)) DEALLOCATE(sourceVector%NodalVector%vector)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NodalFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the nodal matrix.
  SUBROUTINE EquationsMatrices_NodalMatrixInitialise(nodalMatrix,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType) :: nodalMatrix !The nodal matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_NodalMatrixInitialise",err,error,*999)

    nodalMatrix%equationsMatrixNumber=0
    nodalMatrix%numberOfRows=0
    nodalMatrix%numberOfColumns=0
    nodalMatrix%maxNumberOfRows=0
    nodalMatrix%maxNumberOfColumns=0
       
    EXITS("EquationsMatrices_NodalMatrixInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixInitialise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalMatrixInitialise

  !
  !================================================================================================================================
  !

  !>Finalise an nodal matrix and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalMatrixFinalise(nodalMatrix,err,error,*)

    !Argument variables
    TYPE(NodalMatrixType):: nodalMatrix !<The nodal matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_NodalMatrixFinalise",err,error,*999)

    IF(ALLOCATED(nodalMatrix%rowDofs)) DEALLOCATE(nodalMatrix%rowDofs)
    IF(ALLOCATED(nodalMatrix%columnDofs)) DEALLOCATE(nodalMatrix%columnDofs)
    IF(ALLOCATED(nodalMatrix%matrix)) DEALLOCATE(nodalMatrix%matrix)
    
    EXITS("EquationsMatrices_NodalMatrixFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalMatrixFinalise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalMatrixFinalise

  !
  !================================================================================================================================
  !

  !>Initialise the nodal vector
  SUBROUTINE EquationsMatrices_NodalVectorInitialise(nodalVector,err,error,*)

    !Argument variables
    TYPE(NodalVectorType) :: nodalVector !The nodal vector to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsMatrices_NodalVectorInitialise",err,error,*999)

    nodalVector%numberOfRows=0
    nodalVector%maxNumberOfRows=0
       
    EXITS("EquationsMatrices_NodalVectorInitialise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorInitialise",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalVectorInitialise

  !
  !================================================================================================================================
  !

  !>Finalise an nodal vector and deallocate all memory
  SUBROUTINE EquationsMatrices_NodalVectorFinalise(nodalVector,err,error,*)

    !Argument variables
    TYPE(NodalVectorType):: nodalVector !<The nodal vector to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("EquationsMatrices_NodalVectorFinalise",err,error,*999)

    IF(ALLOCATED(nodalVector%rowDofs)) DEALLOCATE(nodalVector%rowDofs)
    IF(ALLOCATED(nodalVector%vector)) DEALLOCATE(nodalVector%vector)
    
    EXITS("EquationsMatrices_NodalVectorFinalise")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NodalVectorFinalise",err,error)

    RETURN 1
  END SUBROUTINE EquationsMatrices_NodalVectorFinalise

  !
  !================================================================================================================================
  !

  !>Adds the Jacobian matrices into the equations Jacobian.
  SUBROUTINE EquationsMatrices_JacobianNodeAdd(equationsMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianMatrixIdx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(VARYING_STRING) :: localError
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_JacobianNodeAdd()")
#endif

    ENTERS("EquationsMatrices_JacobianNodeAdd",err,error,*999)

    IF(ASSOCIATED(equationsMatrices)) THEN
      nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        DO jacobianMatrixIdx=1,nonlinearMatrices%NUMBER_OF_JACOBIANS
          jacobianMatrix=>nonlinearMatrices%JACOBIANS(jacobianMatrixIdx)%PTR
          IF(ASSOCIATED(jacobianMatrix)) THEN
            IF(jacobianMatrix%UPDATE_JACOBIAN) THEN
              !Add in Jacobian element matrices
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(jacobianMatrix%jacobian,jacobianMatrix%NodalJacobian%rowDofs(1: &
                & jacobianMatrix%NodalJacobian%numberOfRows),jacobianMatrix%NodalJacobian%columnDofs(1: &
                & jacobianMatrix%NodalJacobian%numberOfColumns),jacobianMatrix%NodalJacobian%matrix(1: &
                & jacobianMatrix%NodalJacobian%numberOfRows,1:jacobianMatrix%NodalJacobian%numberOfColumns), &
                & err,error,*999)
            ENDIF
          ELSE
            localError="Jacobian matrix for Jacobian matrix index "// &
              & TRIM(NumberToVString(jacobianMatrixIdx,"*",err,error))//" is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDDO !jacobianMatrixIdx
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not allocated.",err,error,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_JacobianNodeAdd()")
#endif
    
    EXITS("EquationsMatrices_JacobianNodeAdd")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianNodeAdd",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_JacobianNodeAdd

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
    
    ENTERS("EQUATIONS_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
      IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
        !Finalise the dynamic element matrices
        DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
          EQUATIONS_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
            CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
      IF(ASSOCIATED(LINEAR_MATRICES)) THEN
        !Finalise the linear element matrices
        DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
            CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ENDIF
      NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MAX_NUMBER_OF_ROWS=0
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MAX_NUMBER_OF_COLUMNS=0
            IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS)
            IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS)
            IF(ALLOCATED(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX)) DEALLOCATE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX)
          ELSE
            CALL FlagError("Nonlinear matrices Jacobian number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                & " is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO
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
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_FINALISE",ERR,ERROR)
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
    INTEGER(INTG) :: rowsNumberOfElements,colsNumberOfElements !Number of elements in the row and col variables whose dofs are present in the element matrix
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
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,COL_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("EQUATIONS_MATRICES_ELEMENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      rowsNumberOfElements=1
      colsNumberOfElements=1
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
                  & rowsNumberOfElements,colsNumberOfElements,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
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
                  & rowsNumberOfElements,colsNumberOfElements,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Equations linear matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            CALL FlagError("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          !Initialise the Jacobian element matrices
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
            FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                COL_FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(matrix_idx)%VARIABLE
                CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(JACOBIAN_MATRIX%ELEMENT_JACOBIAN,FIELD_VARIABLE,COL_FIELD_VARIABLE, &
                  & rowsNumberOfElements,colsNumberOfElements,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Jacobian matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
            !Use RHS variable for residual vector, otherwise first nonlinear variable if no RHS
            RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
            IF(ASSOCIATED(RHS_MAPPING)) THEN
              FIELD_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
            ELSE
              FIELD_VARIABLE=>NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
            ENDIF
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,FIELD_VARIABLE,ERR,ERROR,*999)
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL_CALCULATED=0
          ELSE
            CALL FlagError("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
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
            CALL FlagError("RHS mapping is not associated.",ERR,ERROR,*999)
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
              CALL FlagError("RHS mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_ELEMENT_INITIALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_ELEMENT_INITIALISE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
      IF(ASSOCIATED(EQUATIONS_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      CALL EquationsMatrices_NodalMatrixFinalise(EQUATIONS_MATRIX%NodalMatrix,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_MATRIX%TEMP_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(EQUATIONS_MATRIX%TEMP_VECTOR,ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRIX_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRIX_FINALISE",ERR,ERROR)
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

    ENTERS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES) THEN
        EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
          EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
              IF(ASSOCIATED(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
                LOCAL_ERROR="Equations matrix for dynamic matrix number "//TRIM(NumberToVString(MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is already associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                ALLOCATE(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate equations matrix.",ERR,ERROR,*999)
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
                CALL EquationsMatrices_ElementMatrixInitialise(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
                CALL EquationsMatrices_NodalMatrixInitialise(EQUATIONS_MATRIX%NodalMatrix,ERR,ERROR,*999)
                NULLIFY(EQUATIONS_MATRIX%TEMP_VECTOR)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FlagError("Dynamic matrices equations matrices is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified dynamic matrix number of "//TRIM(NumberToVString(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NumberToVString(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Dynamic matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRIX_FINALISE(DYNAMIC_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRIX_DYNAMIC_INITIALISE",ERR,ERROR)
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

    ENTERS("EQUATIONS_MATRIX_LINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(LINEAR_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES) THEN
        EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
          EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(LINEAR_MAPPING)) THEN
              IF(ASSOCIATED(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
                LOCAL_ERROR="Equations matrix for linear matrix number "//TRIM(NumberToVString(MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is already associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                ALLOCATE(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate equations matrix.",ERR,ERROR,*999)
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
                CALL EquationsMatrices_ElementMatrixInitialise(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
                CALL EquationsMatrices_NodalMatrixInitialise(EQUATIONS_MATRIX%NodalMatrix,ERR,ERROR,*999)
                NULLIFY(EQUATIONS_MATRIX%TEMP_VECTOR)
              ENDIF
            ELSE
              CALL FlagError("Equations mapping linear mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FlagError("Linear matrices equations matrices is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified linear matrix number of "//TRIM(NumberToVString(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NumberToVString(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Linear matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRIX_LINEAR_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRIX_FINALISE(LINEAR_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRIX_LINEAR_INITIALISE",ERR,ERROR)
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
     
    ENTERS("EQUATIONS_MATRICES_DYNAMIC_FINALISE",ERR,ERROR,*999)

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
    
    EXITS("EQUATIONS_MATRICES_DYNAMIC_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_DYNAMIC_FINALISE",ERR,ERROR)
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
     
    ENTERS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE",ERR,ERROR,*998)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
        IF(ASSOCIATED(DYNAMIC_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%DYNAMIC_MATRICES)) THEN
            CALL FlagError("Equations matrices dynamic matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices dynamic matrices.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%DYNAMIC_MATRICES%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
            EQUATIONS_MATRICES%DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES=DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
            ALLOCATE(EQUATIONS_MATRICES%DYNAMIC_MATRICES%MATRICES(DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices dynamic matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,DYNAMIC_MAPPING%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
              NULLIFY(EQUATIONS_MATRICES%DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL EQUATIONS_MATRIX_DYNAMIC_INITIALISE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
            NULLIFY(EQUATIONS_MATRICES%DYNAMIC_MATRICES%TEMP_VECTOR)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices equations mapping is not associated.",ERR,ERROR,*998)        
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_DYNAMIC_FINALISE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_DYNAMIC_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_DYNAMIC_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Adds the Jacobain matrices into the equations Jacobian.
  SUBROUTINE EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobian_matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD()")
#endif

    ENTERS("EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
      IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
        DO jacobian_matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(jacobian_matrix_idx)%PTR
          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
            IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
              !Add in Jacobian element matrices
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(JACOBIAN_MATRIX%JACOBIAN,JACOBIAN_MATRIX%ELEMENT_JACOBIAN%ROW_DOFS(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%COLUMN_DOFS(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS),JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1: &
                & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_ROWS,1:JACOBIAN_MATRIX%ELEMENT_JACOBIAN%NUMBER_OF_COLUMNS), &
                & ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Jacobian matrix for Jacobian matrix index "// &
              & TRIM(NumberToVString(jacobian_matrix_idx,"*",ERR,ERROR))//" is not associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !jacobian_matrix_idx
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD()")
#endif
    
    EXITS("EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_JACOBIAN_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Outputs the equations Jacobian matrices
  SUBROUTINE EQUATIONS_MATRICES_JACOBIAN_OUTPUT(ID,EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations Jacobian matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobian_matrix_idx
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("EQUATIONS_MATRICES_JACOBIAN_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"",ERR,ERROR,*999)
          CALL WRITE_STRING(ID,"Jacobian matrices:",ERR,ERROR,*999)
          DO jacobian_matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(jacobian_matrix_idx)%PTR
            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN        
              CALL WRITE_STRING(ID,"Jacobian matrix:",ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,JACOBIAN_MATRIX%JACOBIAN,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="Jacobian matrix for Jacobian matrix index "// &
                & TRIM(NumberToVString(jacobian_matrix_idx,"*",ERR,ERROR))//" is not associated."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !jacobian_matrix_idx
        ENDIF
      ELSE
        CALL FlagError("Equations matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_JACOBIAN_OUTPUT")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_JACOBIAN_OUTPUT",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_JACOBIAN_OUTPUT
  
  !
  !================================================================================================================================
  !

  !>Sets the Jacobian calculation types of the residual variables
  SUBROUTINE EquationsMatrices_JacobianTypesSet(equationsMatrices,jacobianTypes,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices !<A pointer to the equations matrices to set the Jacobian type for.
    INTEGER(INTG), INTENT(IN) :: jacobianTypes(:) !<jacobianTypes(matrix_idx). The Jacobian calculation type for the matrix_idx'th Jacobian matrix.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    INTEGER(INTG) :: matrixIdx,numberOfjacobians,jacobianType
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsMatrices_JacobianTypesSet",err,error,*999)

    IF(ASSOCIATED(equationsMatrices)) THEN
      IF(equationsMatrices%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",err,error,*999)
      ELSE
        nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
        IF(ASSOCIATED(nonlinearMatrices)) THEN
          numberOfJacobians=SIZE(jacobianTypes,1)
          IF(numberOfJacobians==nonlinearMatrices%NUMBER_OF_JACOBIANS) THEN
            DO matrixIdx=1,numberOfJacobians
              jacobianType=jacobianTypes(matrixIdx)
              SELECT CASE(jacobianType)
              CASE(EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED, &
                  & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
                nonlinearMatrices%JACOBIANS(matrixIdx)%PTR%JACOBIAN_CALCULATION_TYPE=jacobianType
              CASE DEFAULT
                localError="Invalid Jacobian calculation type of " &
                  & //TRIM(NumberToVString(jacobianType,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            END DO
          ELSE
            localError="Invalid number of Jacobian calculation types. The number of types " &
              & //TRIM(NumberToVString(numberOfJacobians,"*",err,error)) &
              & //" should be "//TRIM(NumberToVString(nonlinearMatrices%NUMBER_OF_JACOBIANS,"*",err,error))
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices nonlinear matrices are not associated",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices are not associated",err,error,*999)
    ENDIF

    EXITS("EquationsMatrices_JacobianTypesSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_JacobianTypesSet",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrices_JacobianTypesSet

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
     
    ENTERS("EQUATIONS_MATRICES_LINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_MATRICES)) THEN
      IF(ALLOCATED(LINEAR_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(LINEAR_MATRICES%MATRICES,1)
          CALL EQUATIONS_MATRIX_FINALISE(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(LINEAR_MATRICES%MATRICES)
      ENDIF
      DEALLOCATE(LINEAR_MATRICES)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_LINEAR_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_LINEAR_FINALISE",ERR,ERROR)
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
     
    ENTERS("EQUATIONS_MATRICES_LINEAR_INITIALISE",ERR,ERROR,*998)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
        IF(ASSOCIATED(LINEAR_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%LINEAR_MATRICES)) THEN
            CALL FlagError("Equations matrices linear matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%LINEAR_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices linear matrices.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%LINEAR_MATRICES%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
            EQUATIONS_MATRICES%LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES=LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
            ALLOCATE(EQUATIONS_MATRICES%LINEAR_MATRICES%MATRICES(LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices linear matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
              NULLIFY(EQUATIONS_MATRICES%LINEAR_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL EQUATIONS_MATRIX_LINEAR_INITIALISE(EQUATIONS_MATRICES%LINEAR_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices equations mapping is not associated.",ERR,ERROR,*998)        
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_LINEAR_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_LINEAR_FINALISE(EQUATIONS_MATRICES%LINEAR_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_LINEAR_INITIALISE",ERR,ERROR)
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
    INTEGER(INTG) :: matrix_idx
     
    ENTERS("EQUATIONS_MATRICES_NONLINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
      IF(ALLOCATED(NONLINEAR_MATRICES%JACOBIANS)) THEN
        DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          CALL EQUATIONS_JACOBIAN_FINALISE(NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO
        DEALLOCATE(NONLINEAR_MATRICES%JACOBIANS)
      ENDIF
      IF(ASSOCIATED(NONLINEAR_MATRICES%RESIDUAL)) CALL DISTRIBUTED_VECTOR_DESTROY(NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ERR,ERROR,*999)
      CALL EquationsMatrices_NodalVectorFinalise(NONLINEAR_MATRICES%NodalResidual,ERR,ERROR,*999)
      DEALLOCATE(NONLINEAR_MATRICES)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_NONLINEAR_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_NONLINEAR_FINALISE",ERR,ERROR)
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
    INTEGER(INTG) :: matrix_idx,DUMMY_ERR
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    ENTERS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
        IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%NONLINEAR_MATRICES)) THEN
            CALL FlagError("Equations matrices nonlinear matrices is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices nonlinear matrices.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%UPDATE_RESIDUAL=.TRUE.
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(EQUATIONS_MATRICES%NONLINEAR_MATRICES%RESIDUAL)
            CALL EquationsMatrices_ElementVectorInitialise(EQUATIONS_MATRICES%NONLINEAR_MATRICES%ELEMENT_RESIDUAL,ERR,ERROR,*999)
            CALL EquationsMatrices_NodalVectorInitialise(EQUATIONS_MATRICES%NONLINEAR_MATRICES%NodalResidual,ERR,ERROR,*999)
            EQUATIONS_MATRICES%NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS=NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
            ALLOCATE(EQUATIONS_MATRICES%NONLINEAR_MATRICES%JACOBIANS(NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices Jacobian matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,NONLINEAR_MAPPING%NUMBER_OF_RESIDUAL_VARIABLES
              NULLIFY(EQUATIONS_MATRICES%NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR)
              CALL EQUATIONS_JACOBIAN_INITIALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_NONLINEAR_FINALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_NONLINEAR_INITIALISE",ERR,ERROR)
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
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    
    ENTERS("EQUATIONS_MATRICES_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"",ERR,ERROR,*999)
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
              CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
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
              CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          CALL WRITE_STRING(ID,"Nonlinear vectors:",ERR,ERROR,*999)
          IF(ASSOCIATED(NONLINEAR_MATRICES%RESIDUAL)) THEN
            CALL WRITE_STRING(ID,"Residual vector:",ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,NONLINEAR_MATRICES%RESIDUAL,ERR,ERROR,*999)
          ELSE
            CALL FlagError("Nonlinear matrices residual is not associated.",ERR,ERROR,*999)
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
        CALL FlagError("Equations matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_OUTPUT")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_OUTPUT",ERR,ERROR)
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
     
    ENTERS("EQUATIONS_MATRICES_RHS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_VECTOR)) THEN
      IF(ASSOCIATED(RHS_VECTOR%VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(RHS_VECTOR%VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
      CALL EquationsMatrices_NodalVectorFinalise(RHS_VECTOR%NodalVector,ERR,ERROR,*999)
      DEALLOCATE(RHS_VECTOR)
    ENDIF      
     
    EXITS("EQUATIONS_MATRICES_RHS_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_RHS_FINALISE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_RHS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
        IF(ASSOCIATED(RHS_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%RHS_VECTOR)) THEN
            CALL FlagError("Equations matrices RHS vector is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%RHS_VECTOR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices RHS vector.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%RHS_VECTOR%UPDATE_VECTOR=.TRUE.
            EQUATIONS_MATRICES%RHS_VECTOR%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(EQUATIONS_MATRICES%RHS_VECTOR%VECTOR)
            CALL EquationsMatrices_ElementVectorInitialise(EQUATIONS_MATRICES%RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
            CALL EquationsMatrices_NodalVectorInitialise(EQUATIONS_MATRICES%RHS_VECTOR%NodalVector,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices equation mapping is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_RHS_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_RHS_FINALISE(EQUATIONS_MATRICES%RHS_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_RHS_INITIALISE",ERR,ERROR)
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
     
    ENTERS("EQUATIONS_MATRICES_SOURCE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOURCE_VECTOR)) THEN
      IF(ASSOCIATED(SOURCE_VECTOR%VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(SOURCE_VECTOR%VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(SOURCE_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
      CALL EquationsMatrices_NodalVectorFinalise(SOURCE_VECTOR%NodalVector,ERR,ERROR,*999)
      DEALLOCATE(SOURCE_VECTOR)
    ENDIF      
     
    EXITS("EQUATIONS_MATRICES_SOURCE_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_SOURCE_FINALISE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_SOURCE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES%SOURCE_VECTOR)) THEN
            CALL FlagError("Equations matrices source vector is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%SOURCE_VECTOR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations matrices source vector.",ERR,ERROR,*999)
            EQUATIONS_MATRICES%SOURCE_VECTOR%UPDATE_VECTOR=.TRUE.
            EQUATIONS_MATRICES%SOURCE_VECTOR%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(EQUATIONS_MATRICES%SOURCE_VECTOR%VECTOR)
            CALL EquationsMatrices_ElementVectorInitialise(EQUATIONS_MATRICES%SOURCE_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
            CALL EquationsMatrices_NodalVectorInitialise(EQUATIONS_MATRICES%SOURCE_VECTOR%NodalVector,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Equations matrices equation mapping is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_SOURCE_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_SOURCE_FINALISE(EQUATIONS_MATRICES%SOURCE_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_SOURCE_INITIALISE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have already been finished.",ERR,ERROR,*999)
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
                  LOCAL_ERROR="The specified lumping type of "//TRIM(NumberToVString(LUMPING_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the lumping type array ("//TRIM(NumberToVString(SIZE(LUMPING_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NumberToVString(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have already been finished.",ERR,ERROR,*999)
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
                  LOCAL_ERROR="The specified storage type of "//TRIM(NumberToVString(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the dynamic matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NumberToVString(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NumberToVString(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
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
                  LOCAL_ERROR="The specified storage type of "//TRIM(NumberToVString(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the linear matrix number "//TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NumberToVString(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of linear matrices ("// &
              & TRIM(NumberToVString(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet0(EQUATIONS_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("EquationsMatrices_NonlinearStorageTypeSet0",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          IF(SIZE(STORAGE_TYPE,1)==NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                SELECT CASE(STORAGE_TYPE(matrix_idx))
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
                  LOCAL_ERROR="The specified storage type of "//TRIM(NumberToVString(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                    & " for the Jacobian matrix is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Jacobian matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NumberToVString(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of Jacobian matrices ("// &
              & TRIM(NumberToVString(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NonlinearStorageTypeSet0")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NonlinearStorageTypeSet0",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of all nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet1(EQUATIONS_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE !<STORAGE_TYPE. The storage type for all Jacobian equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: STORAGE_TYPES(:)
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES

    ENTERS("EquationsMatrices_NonlinearStorageTypeSet1",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          ALLOCATE(STORAGE_TYPES(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate storage types.",ERR,ERROR,*999)
          STORAGE_TYPES=STORAGE_TYPE
          CALL EquationsMatrices_NonlinearStorageTypeSet0(EQUATIONS_MATRICES,STORAGE_TYPES,ERR,ERROR,*999)
          DEALLOCATE(STORAGE_TYPES)
        ELSE
          CALL FlagError("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("EquationsMatrices_NonlinearStorageTypeSet1")
    RETURN
999 ERRORSEXITS("EquationsMatrices_NonlinearStorageTypeSet1",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStorageTypeSet1

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the dynamic equations matrices
  SUBROUTINE EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
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

    ENTERS("EquationsMatrices_DynamicStructureTypeSet",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
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
                CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NODAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NumberToVString(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for dynamic matrix number "// &
                    & TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NumberToVString(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of dynamic matrices ("// &
              & TRIM(NumberToVString(DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EquationsMatrices_DynamicStructureTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_DynamicStructureTypeSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_DynamicStructureTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the linear equations matrices
  SUBROUTINE EquationsMatrices_LinearStructureTypeSet(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
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

    ENTERS("EquationsMatrices_LinearStructureTypeSet",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
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
                CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
                  EQUATIONS_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NODAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NumberToVString(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for linear matrix number "// &
                    & TRIM(NumberToVString(matrix_idx,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NumberToVString(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of linear matrices ("// &
              & TRIM(NumberToVString(LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EquationsMatrices_LinearStructureTypeSet")
    RETURN
999 ERRORSEXITS("EquationsMatrices_LinearStructureTypeSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_LinearStructureTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet0(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The structure type for the matrix_idx'th Jacobian equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("EquationsMatrices_NonlinearStructureTypeSet0",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS) THEN
            DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
              JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                SELECT CASE(STRUCTURE_TYPE(matrix_idx))
                CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
                CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_FEM_STRUCTURE
                CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_DIAGONAL_STRUCTURE
                CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
                  JACOBIAN_MATRIX%STRUCTURE_TYPE=EQUATIONS_MATRIX_NODAL_STRUCTURE
                CASE DEFAULT
                  LOCAL_ERROR="The specified strucutre type of "// &
                    & TRIM(NumberToVString(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for the Jacobian matrix is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NumberToVString(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of Jacobian matrices ("// &
              & TRIM(NumberToVString(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS,"*",ERR,ERROR))//")."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("EquationsMatrices_NonlinearStructureTypeSet0")
    RETURN
999 ERRORS("EquationsMatrices_NonlinearStructureTypeSet0",ERR,ERROR)
    EXITS("EquationsMatrices_NonlinearStructureTypeSet0")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet0

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of all nonlinear (Jacobian) equations matrices
  SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet1(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE !<The structure type for all Jacobian equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: STRUCTURE_TYPES(:)
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES

    ENTERS("EquationsMatrices_NonlinearStructureTypeSet1",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FlagError("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          ALLOCATE(STRUCTURE_TYPES(NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate storage types.",ERR,ERROR,*999)
          STRUCTURE_TYPES=STRUCTURE_TYPE
          CALL EquationsMatrices_NonlinearStructureTypeSet0(EQUATIONS_MATRICES,STRUCTURE_TYPES,ERR,ERROR,*999)
          DEALLOCATE(STRUCTURE_TYPES)
        ELSE
          CALL FlagError("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("EquationsMatrices_NonlinearStructureTypeSet1")
    RETURN
999 ERRORS("EquationsMatrices_NonlinearStructureTypeSet1",ERR,ERROR)
    EXITS("EquationsMatrices_NonlinearStructureTypeSet1")
    RETURN 1
    
  END SUBROUTINE EquationsMatrices_NonlinearStructureTypeSet1

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
   
    ENTERS("EQUATIONS_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      CALL EQUATIONS_MATRICES_DYNAMIC_FINALISE(EQUATIONS_MATRICES%DYNAMIC_MATRICES,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_LINEAR_FINALISE(EQUATIONS_MATRICES%LINEAR_MATRICES,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_NONLINEAR_FINALISE(EQUATIONS_MATRICES%NONLINEAR_MATRICES,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_RHS_FINALISE(EQUATIONS_MATRICES%RHS_VECTOR,ERR,ERROR,*999)      
      CALL EQUATIONS_MATRICES_SOURCE_FINALISE(EQUATIONS_MATRICES%SOURCE_VECTOR,ERR,ERROR,*999)      
      DEALLOCATE(EQUATIONS_MATRICES)
    ENDIF
       
    EXITS("EQUATIONS_MATRICES_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_FINALISE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MATRICES)) THEN
        CALL FlagError("Equations matrices is already associated for this equations.",ERR,ERROR,*998)
      ELSE
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
            ALLOCATE(EQUATIONS%EQUATIONS_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate equations equations matrices.",ERR,ERROR,*999)
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
            CALL FlagError("Equations mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("EQUATIONS_MATRICES_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_FINALISE(EQUATIONS%EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_MATRICES_INITIALISE",ERR,ERROR)
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
    
    ENTERS("EQUATIONS_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)
    
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
              CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
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
              CALL FlagError("Equations matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
      ENDIF
      IF(SELECTION_TYPE==EQUATIONS_MATRICES_ALL.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_NONLINEAR_ONLY.OR. &
        & SELECTION_TYPE==EQUATIONS_MATRICES_JACOBIAN_ONLY) THEN
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
          DO matrix_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(matrix_idx)%PTR
            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
              IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
                CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(JACOBIAN_MATRIX%JACOBIAN,VALUE,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Jacobian matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO
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
      CALL FlagError("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("EQUATIONS_MATRICES_VALUES_INITIALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_MATRICES_VALUES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_VALUES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a equations matrix.
  SUBROUTINE EquationsMatrix_StructureCalculate(equationsMatrix,numberOfNonZeros,rowIndices,columnIndices,list,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix !<A pointer to the equations matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    type(LinkedList),pointer :: list(:) 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) ::  columnIdx,dummyErr,elementIdx,globalColumn,localColumn,local_ny,matrixNumber,mk,mp,ne,nh,nh2,nn,nnk,np
    INTEGER(INTG) ::  numberOfColumns,nyy,nyyg,npg,nhg,local_cols,local_dof,mv
    INTEGER(INTG) :: dofIdx,nodeIdx,componentIdx,localDofIdx
    INTEGER(INTG) :: versionIdx,derivativeIdx,numberOfDerivatives,numberOfVersions
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: dependentDofsDomainMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: dependentDofsParamMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError
    ENTERS("EquationsMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    IF(ASSOCIATED(equationsMatrix)) THEN
      IF(.NOT.ASSOCIATED(rowIndices)) THEN
        IF(.NOT.ASSOCIATED(columnIndices)) THEN
          matrixNumber=equationsMatrix%MATRIX_NUMBER
          SELECT CASE(equationsMatrix%STRUCTURE_TYPE)
          CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
            CALL FlagError("There is no structure to calculate for a matrix with no structure.",err,error,*998)
          CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
            SELECT CASE(equationsMatrix%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              linearMatrices=>equationsMatrix%LINEAR_MATRICES
              dynamicMatrices=>equationsMatrix%DYNAMIC_MATRICES
              IF(ASSOCIATED(dynamicMatrices).OR.ASSOCIATED(linearMatrices)) THEN
                IF(ASSOCIATED(dynamicMatrices)) THEN
                  equationsMatrices=>dynamicMatrices%EQUATIONS_MATRICES
                ELSE
                  equationsMatrices=>linearMatrices%EQUATIONS_MATRICES
                ENDIF
                IF(ASSOCIATED(equationsMatrices)) THEN
                  equations=>equationsMatrices%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    equationsMapping=>equationsMatrices%EQUATIONS_MAPPING
                    IF(ASSOCIATED(equationsMapping)) THEN
                      dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
                      linearMapping=>equationsMapping%LINEAR_MAPPING
                      IF(ASSOCIATED(dynamicMapping).OR.ASSOCIATED(linearMapping)) THEN
                        equationsSet=>equations%EQUATIONS_SET
                        IF(ASSOCIATED(equationsSet)) THEN
                          dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            IF(ASSOCIATED(dynamicMatrices)) THEN
                              fieldVariable=>dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%VARIABLE
                            ELSE
                              fieldVariable=>linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%VARIABLE
                            ENDIF
                            IF(ASSOCIATED(fieldVariable)) THEN
                              dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
                              IF(ASSOCIATED(dependentDofsDomainMapping)) THEN
                                dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(dependentDofsParamMapping)) THEN
                                  !Allocate lists
                                  ALLOCATE(columnIndicesLists(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                  IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
                                  !Allocate row indices
                                  ALLOCATE(rowIndices(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                  IF(ERR/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
                                  rowIndices(1)=1
                                  
                                  !First, loop over the rows and calculate the number of non-zeros
                                  numberOfNonZeros=0
                                  DO local_ny=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                    IF(dependentDofsParamMapping%DOF_TYPE(1,local_ny)==FIELD_NODE_DOF_TYPE) THEN
                                      nyy=dependentDofsParamMapping%DOF_TYPE(2,local_ny)!value for a particular field dof (local_ny)
                                      np=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,nyy)!node number (np) of the field parameter
                                      nh=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,nyy)!component number (nh) of the field parameter
                                      domainNodes=>fieldVariable%COMPONENTS(nh)%DOMAIN%TOPOLOGY%NODES
                                      
                                      !Set up list
                                      NULLIFY(columnIndicesLists(local_ny)%PTR)
                                      CALL LIST_CREATE_START(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                      CALL LIST_DATA_TYPE_SET(columnIndicesLists(local_ny)%PTR,LIST_INTG_TYPE,err,error,*999)
                                      CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(local_ny)%PTR,domainNodes%NODES(np)% &
                                        & NUMBER_OF_SURROUNDING_ELEMENTS*fieldVariable%COMPONENTS(nh)% &
                                        & maxNumberElementInterpolationParameters,err,error,*999)
                                      CALL LIST_CREATE_FINISH(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                      !Loop over all elements containing the dof
                                      DO elementIdx=1,domainNodes%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                                        ne=domainNodes%NODES(np)%SURROUNDING_ELEMENTS(elementIdx)
                                        DO nh2=1,fieldVariable%NUMBER_OF_COMPONENTS
                                          domainElements=>fieldVariable%COMPONENTS(nh2)%DOMAIN%TOPOLOGY%ELEMENTS
                                          basis=>domainElements%ELEMENTS(ne)%BASIS
                                          DO nn=1,basis%NUMBER_OF_NODES
                                            mp=domainElements%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                            DO nnk=1,basis%NUMBER_OF_DERIVATIVES(nn)
                                              mk=domainElements%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                              mv=domainElements%ELEMENTS(ne)%elementVersions(nnk,nn)
                                              !Find the local and global column and add the global column to the indices list
                                              localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                                & NODES(mp)%DERIVATIVES(mk)%VERSIONS(mv)
                                              globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                          
                                              CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn,err,error,*999)
                                                
                                            ENDDO !mk
                                          ENDDO !nn
                                        ENDDO !nh2
                                      ENDDO !elementIdx
                                      CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                      CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_ny)%PTR,numberOfColumns, &
                                        & err,error,*999)
                                      numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                      rowIndices(local_ny+1)=numberOfNonZeros+1
                                    ELSE
                                      localError="Local dof number "//TRIM(NumberToVString(local_ny,"*",err,error))// &
                                        & " is not a node based dof."
                                      CALL FlagError(localError,err,error,*999)
                                    ENDIF
                                  ENDDO !local_ny
                                  
                                  
                                  !Allocate and setup the column locations
                                  ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)

                                  ALLOCATE(list(dependentDofsDomainMapping%NUMBER_OF_GLOBAL))

                                  IF(ERR/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
                                  DO local_ny=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                   
                                    CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(local_ny)%PTR,numberOfColumns,columns, &
                                      & err,error,*999)        
                                    DO columnIdx=1,numberOfColumns
                                      !COLUMNS store the list of nonzero column indices for each local row (local_ny)
                                      columnIndices(rowIndices(local_ny)+columnIdx-1)=columns(columnIdx) 

                                      ! global to local columns
                                       IF(ASSOCIATED(linearMapping).OR.ASSOCIATED(dynamicMapping)) THEN
                                         IF(ASSOCIATED(dynamicMatrices)) THEN
                                           local_cols=equationsMatrices%equations_mapping%dynamic_mapping &
                                             & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                             & (columns(columnIdx))%LOCAL_NUMBER(1)
                                           local_dof = local_cols
                                           ! Column to dof mapping?
                                           !local_dof=equationsMatrices%equations_mapping%dynamic_mapping% &
                                            ! & equations_matrix_to_var_maps(1)%column_to_dof_map(local_cols)
                                         ELSE
                                           local_cols=equationsMatrices%equations_mapping%linear_mapping &
                                             & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                             & (columns(columnIdx))%LOCAL_NUMBER(1)
                                           local_dof = local_cols
                                         ENDIF
                                       ENDIF
                                       nyyg=dependentDofsParamMapping%DOF_TYPE(2,local_dof)
                                       npg=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,nyyg)
                                       nhg=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,nyyg)
                                       domainNodes=>fieldVariable%COMPONENTS(nhg)%DOMAIN%TOPOLOGY%NODES
                            
                                      ! Check whether boundary node    
                                      IF(domainNodes%NODES(npg)%BOUNDARY_NODE)THEN
                                        CALL LinkedList_Add(list(columns(columnIdx)),local_ny,ERR,ERROR,*999)
                                      ENDIF
                                    
                                    ENDDO !columnIdx
                                    DEALLOCATE(columns)                                    
                                  ENDDO !local_ny

                                 
                                  IF(DIAGNOSTICS1) THEN
                                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix structure:",err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix number : ",matrixNumber, &
                                      & err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                      & dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                      & dependentDofsDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                      & numberOfNonZeros,err,error,*999)
                                    IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                      & dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                      sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
                                        & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
                                      CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ", &
                                        & sparsity,"F6.2",err,error,*999)
                                    ENDIF
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,dependentDofsDomainMapping% &
                                      & TOTAL_NUMBER_OF_LOCAL+1,8,8,rowIndices,'("  Row indices    :",8(X,I13))', &
                                      & '(18X,8(X,I13))',err,error,*999)
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
                                      & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Dependent dofs parameter mapping is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Dependent dofs domain mapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Either equations mapping dynamic mapping or linear mapping is not associated.", &
                          & err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Dynamic or linear matrices equations matrices is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Either equations matrix dynamic or linear matrices is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NumberToVString(equationsMatrix%STORAGE_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
            SELECT CASE(equationsMatrix%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              linearMatrices=>equationsMatrix%LINEAR_MATRICES
              dynamicMatrices=>equationsMatrix%DYNAMIC_MATRICES
              IF(ASSOCIATED(dynamicMatrices).OR.ASSOCIATED(linearMatrices)) THEN
                IF(ASSOCIATED(dynamicMatrices)) THEN
                  equationsMatrices=>dynamicMatrices%EQUATIONS_MATRICES
                ELSE
                  equationsMatrices=>linearMatrices%EQUATIONS_MATRICES
                ENDIF
                IF(ASSOCIATED(equationsMatrices)) THEN
                  equations=>equationsMatrices%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    equationsMapping=>equationsMatrices%EQUATIONS_MAPPING
                    IF(ASSOCIATED(equationsMapping)) THEN
                      dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
                      linearMapping=>equationsMapping%LINEAR_MAPPING
                      IF(ASSOCIATED(dynamicMapping).OR.ASSOCIATED(linearMapping)) THEN
                        equationsSet=>equations%EQUATIONS_SET
                        IF(ASSOCIATED(equationsSet)) THEN
                          dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            IF(ASSOCIATED(dynamicMatrices)) THEN
                              fieldVariable=>dynamicMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%VARIABLE
                            ELSE
                              fieldVariable=>linearMapping%EQUATIONS_MATRIX_TO_VAR_MAPS(matrixNumber)%VARIABLE
                            ENDIF
                            IF(ASSOCIATED(fieldVariable)) THEN
                              dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
                              IF(ASSOCIATED(dependentDofsDomainMapping)) THEN
                                dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(dependentDofsParamMapping)) THEN
                                  !Allocate lists
                                  ALLOCATE(columnIndicesLists(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                  IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
                                  !Allocate row indices
                                  ALLOCATE(rowIndices(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                  IF(ERR/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
                                  rowIndices(1)=1
                                  
                                  !First, loop over the rows and calculate the number of non-zeros
                                  numberOfNonZeros=0
                                  DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                    IF(dependentDofsParamMapping%DOF_TYPE(1,localDofIdx)==FIELD_NODE_DOF_TYPE) THEN
                                      dofIdx=dependentDofsParamMapping%DOF_TYPE(2,localDofIdx)!value for a particular field dof (localDofIdx)
                                      nodeIdx=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx)!node number (np) of the field parameter
                                      componentIdx=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx)!component number (nh) of the field parameter
                                      domainNodes=>fieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES
                                      
                                      !Set up list
                                      NULLIFY(columnIndicesLists(localDofIdx)%PTR)
                                      CALL LIST_CREATE_START(columnIndicesLists(localDofIdx)%PTR,err,error,*999)
                                      CALL LIST_DATA_TYPE_SET(columnIndicesLists(localDofIdx)%PTR,LIST_INTG_TYPE,err,error,*999)

                                      CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(localDofIdx)%PTR, &
                                        & fieldVariable%NUMBER_OF_COMPONENTS* &
                                        & fieldVariable%maxNumberElementInterpolationParameters,err,error,*999)

                                      CALL LIST_CREATE_FINISH(columnIndicesLists(localDofIdx)%PTR,err,error,*999)
                                      !Loop over all components,nodes,derivatives, and versions
                                      DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
                                        numberOfDerivatives=fieldVariable%components(componentIdx)%domain%topology% &
                                         & nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                                        DO derivativeIdx=1,numberOfDerivatives
                                          numberOfVersions=fieldVariable%components(componentIdx)%domain%topology% &
                                           & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                                          DO versionIdx=1,numberOfVersions
                                            localColumn=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                             & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                                            globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)

                                            CALL LIST_ITEM_ADD(columnIndicesLists(localDofIdx)%PTR,globalColumn,err,error,*999)
                                            
                                          ENDDO !versionIdx
                                        ENDDO !derivativeIdx
                                      ENDDO !componentIdx

                                      CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(localDofIdx)%PTR,err,error,*999)
                                      CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(localDofIdx)%PTR,numberOfColumns, &
                                        & err,error,*999)
                                      numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                      rowIndices(localDofIdx+1)=numberOfNonZeros+1
                                    ELSE
                                      localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))// &
                                        & " is not a node based dof."
                                      CALL FlagError(localError,err,error,*999)
                                    ENDIF
                                  ENDDO !localDofIdx
                                  
                                  !Allocate and setup the column locations
                                  ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)
                                  ALLOCATE(list(dependentDofsDomainMapping%NUMBER_OF_GLOBAL))

                                  IF(ERR/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
                                  DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                   
                                    CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(localDofIdx)%PTR,numberOfColumns,columns, &
                                      & err,error,*999)        
                                    DO columnIdx=1,numberOfColumns
                                      !columns store the list of nonzero column indices for each local row (localDofIdx)
                                      columnIndices(rowIndices(localDofIdx)+columnIdx-1)=columns(columnIdx) 

                                      ! global to local columns
                                       IF(ASSOCIATED(linearMapping).OR.ASSOCIATED(dynamicMapping)) THEN
                                         IF(ASSOCIATED(dynamicMatrices)) THEN
                                           local_cols=equationsMatrices%equations_mapping%dynamic_mapping &
                                             & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                             & (columns(columnIdx))%LOCAL_NUMBER(1)
                                           local_dof = local_cols
                                           ! Column to dof mapping?
                                           !local_dof=equationsMatrices%equations_mapping%dynamic_mapping% &
                                            ! & equations_matrix_to_var_maps(1)%column_to_dof_map(local_cols)
                                         ELSE
                                           local_cols=equationsMatrices%equations_mapping%linear_mapping &
                                             & %equations_matrix_to_var_maps(1)%column_dofs_mapping%global_to_local_map &
                                             & (columns(columnIdx))%LOCAL_NUMBER(1)
                                           local_dof = local_cols
                                         ENDIF
                                       ENDIF
                                       nyyg=dependentDofsParamMapping%DOF_TYPE(2,local_dof)
                                       npg=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,nyyg)
                                       nhg=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,nyyg)
                                       domainNodes=>fieldVariable%COMPONENTS(nhg)%DOMAIN%TOPOLOGY%NODES
                            
                                      ! Check whether boundary node    
                                      IF(domainNodes%NODES(npg)%BOUNDARY_NODE)THEN
                                        CALL LinkedList_Add(list(columns(columnIdx)),localDofIdx,err,error,*999)
                                      ENDIF
                                    
                                    ENDDO !columnIdx
                                    DEALLOCATE(columns)                                    
                                  ENDDO !localDofIdx
                                 
                                  IF(DIAGNOSTICS1) THEN
                                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix structure:",err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix number : ",matrixNumber, &
                                      & err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                      & dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                      & dependentDofsDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
                                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                      & numberOfNonZeros,err,error,*999)
                                    IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                      & dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                      sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
                                        & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
                                      CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ", &
                                        & sparsity,"F6.2",err,error,*999)
                                    ENDIF
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,dependentDofsDomainMapping% &
                                      & TOTAL_NUMBER_OF_LOCAL+1,8,8,rowIndices,'("  Row indices    :",8(X,I13))', &
                                      & '(18X,8(X,I13))',err,error,*999)
                                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices, &
                                      & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Dependent dofs parameter mapping is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Dependent dofs domain mapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Either equations mapping dynamic mapping or linear mapping is not associated.", &
                          & err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Dynamic or linear matrices equations matrices is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Either equations matrix dynamic or linear matrices is not associated.",err,error,*999)
              ENDIF

            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NumberToVString(equationsMatrix%STORAGE_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(EQUATIONS_MATRIX_DIAGONAL_STRUCTURE)
            CALL FlagError("There is not structure to calculate for a diagonal matrix.",err,error,*998)
          CASE DEFAULT
            localError="The matrix structure type of "// &
              & TRIM(NumberToVString(equationsMatrix%STRUCTURE_TYPE,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*998)
          END SELECT
        ELSE
          CALL FlagError("Column indices is already associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Row indieces is already associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Equations matrix is not associated.",err,error,*999)
    ENDIF
      
    EXITS("EquationsMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%PTR)) &
          & CALL LIST_DESTROY(columnIndicesLists(localDofIdx)%PTR,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("EquationsMatrix_StructureCalculate",err,error)
    RETURN 1
  END SUBROUTINE EquationsMatrix_StructureCalculate

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a Jacobian matrix.
  SUBROUTINE JacobianMatrix_StructureCalculate(jacobianMatrix,numberOfNonZeros,rowIndices,columnIndices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix !<A pointer to the Jacobian matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: numberOfNonZeros !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: rowIndices(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: columnIndices(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) ::  columnIdx,dummyErr,elementIdx,globalColumn,localColumn,local_ny,mk,mp,ne,nh,nh2,nn,nnk,np,mv, &
      & numberOfColumns,nyy,matrixNumber
    INTEGER(INTG) :: dofIdx,nodeIdx,componentIdx,versionIdx,derivativeIdx,numberOfVersions,numberOfDerivatives
    INTEGER(INTG) :: localDofIdx
    INTEGER(INTG), ALLOCATABLE :: columns(:)
    REAL(DP) :: sparsity
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: dependentDofsDomainMapping,rowDofsDomainMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONlinearMatrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: dependentDofsParamMapping,rowDofsParamMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,rowVariable
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: columnIndicesLists(:)
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("JacobianMatrix_StructureCalculate",err,error,*998)

    numberOfNonZeros=0
    IF(ASSOCIATED(jacobianMatrix)) THEN
      matrixNumber=jacobianMatrix%JACOBIAN_NUMBER
      IF(.NOT.ASSOCIATED(rowIndices)) THEN
        IF(.NOT.ASSOCIATED(columnIndices)) THEN
          SELECT CASE(jacobianMatrix%STRUCTURE_TYPE)
          CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
            CALL FlagError("Not implemented.",err,error,*998)
          CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
            SELECT CASE(jacobianMatrix%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              nonlinearMatrices=>jacobianMatrix%NONLINEAR_MATRICES
              IF(ASSOCIATED(nonlinearMatrices)) THEN
                equationsMatrices=>nonlinearMatrices%EQUATIONS_MATRICES
                IF(ASSOCIATED(equationsMatrices)) THEN
                  equations=>equationsMatrices%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    equationsMapping=>equationsMatrices%EQUATIONS_MAPPING
                    IF(ASSOCIATED(equationsMapping)) THEN
                      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
                      IF(ASSOCIATED(nonlinearMapping)) THEN
                        equationsSet=>equations%EQUATIONS_SET
                        IF(ASSOCIATED(equationsSet)) THEN
                          dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixNumber)%VARIABLE
                            IF(ASSOCIATED(fieldVariable)) THEN
                              dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
                              IF(ASSOCIATED(dependentDofsDomainMapping)) THEN
                                dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(dependentDofsParamMapping)) THEN
                                  !If RHS variable exists, use this for row DOFs, else use the first nonlinear variable
                                  IF(ASSOCIATED(equationsMapping%RHS_MAPPING)) THEN
                                    rowVariable=>equationsMapping%RHS_MAPPING%RHS_VARIABLE
                                  ELSE
                                    rowVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
                                  ENDIF
                                  IF(ASSOCIATED(rowVariable)) THEN
                                    rowDofsDomainMapping=>rowVariable%DOMAIN_MAPPING
                                    rowDofsParamMapping=>rowVariable%DOF_TO_PARAM_MAP
                                  ELSE
                                    CALL FlagError("RHS or first nonlinear variable is not associated",err,error,*999)
                                  ENDIF
                                  IF(ASSOCIATED(rowDofsDomainMapping)) THEN
                                    IF(ASSOCIATED(rowDofsParamMapping)) THEN
                                      !Allocate lists
                                      ALLOCATE(columnIndicesLists(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
                                      !Allocate row indices
                                      ALLOCATE(rowIndices(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
                                      rowIndices(1)=1
                                      !First, loop over the rows and calculate the number of non-zeros
                                      numberOfNonZeros=0
                                      DO local_ny=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        SELECT CASE(rowDofsParamMapping%DOF_TYPE(1,local_ny))
                                        CASE(FIELD_CONSTANT_INTERPOLATION)
                                          CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_NODE_DOF_TYPE)
                                          nyy=rowDofsParamMapping%DOF_TYPE(2,local_ny)
                                          np=rowDofsParamMapping%NODE_DOF2PARAM_MAP(3,nyy) !node number
                                          nh=rowDofsParamMapping%NODE_DOF2PARAM_MAP(4,nyy) !component number
                                          domainNodes=>rowVariable%COMPONENTS(nh)%DOMAIN%TOPOLOGY%NODES
                                          !Set up list
                                          NULLIFY(columnIndicesLists(local_ny)%PTR)
                                          CALL LIST_CREATE_START(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_DATA_TYPE_SET(columnIndicesLists(local_ny)%PTR,LIST_INTG_TYPE,err,error,*999)
                                          CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(local_ny)%PTR,domainNodes%NODES(np)% &
                                            & NUMBER_OF_SURROUNDING_ELEMENTS*rowVariable%COMPONENTS(nh)% &
                                            & maxNumberElementInterpolationParameters,err,error,*999)
                                          CALL LIST_CREATE_FINISH(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          !Loop over all elements containing the dof
                                          DO elementIdx=1,domainNodes%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                                            ne=domainNodes%NODES(np)%SURROUNDING_ELEMENTS(elementIdx)
                                            DO nh2=1,fieldVariable%NUMBER_OF_COMPONENTS
                                              SELECT CASE(fieldVariable%COMPONENTS(nh2)%INTERPOLATION_TYPE)
                                              CASE(FIELD_CONSTANT_INTERPOLATION)
                                                ! do nothing? this will probably never be encountered...?
                                              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                                localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                  & ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                                                globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                                CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn,err,error,*999)
                                              CASE(FIELD_NODE_BASED_INTERPOLATION)
                                                domainElements=>fieldVariable%COMPONENTS(nh2)%DOMAIN%TOPOLOGY%ELEMENTS
                                                basis=>domainElements%ELEMENTS(ne)%BASIS
                                                DO nn=1,basis%NUMBER_OF_NODES
                                                  mp=domainElements%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                                  DO nnk=1,basis%NUMBER_OF_DERIVATIVES(nn)
                                                    mk=domainElements%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                                    mv=domainElements%ELEMENTS(ne)%elementVersions(nnk,nn)
                                                    !Find the local and global column and add the global column to the indices list
                                                    localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                      & NODE_PARAM2DOF_MAP%NODES(mp)%DERIVATIVES(mk)%VERSIONS(mv)
                                                    globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                                    CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn, &
                                                      & err,error,*999)
                                                  ENDDO !mk
                                                ENDDO !nn
                                              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                                CALL FlagError("Grid point based interpolation is not implemented yet.",& 
                                                  & err,error,*999)
                                              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                                CALL FlagError("Gauss point based interpolation is not implemented yet.",&
                                                  & err,error,*999)
                                              CASE DEFAULT
                                                localError="Local dof number "//TRIM(NumberToVString(local_ny,"*",err,error))// &
                                                  & " has invalid interpolation type."
                                                CALL FlagError(localError,err,error,*999)
                                              END SELECT
                                            ENDDO !nh2
                                          ENDDO !elementIdx
                                          CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_ny)%PTR,numberOfColumns, &
                                            & err,error,*999)
                                          numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                          rowIndices(local_ny+1)=numberOfNonZeros+1
                                        CASE(FIELD_ELEMENT_DOF_TYPE)
                                          ! row corresponds to a variable that's element-wisely interpolated
                                          nyy=rowDofsParamMapping%DOF_TYPE(2,local_ny)          ! nyy = index in ELEMENT_DOF2PARAM_MAP
                                          ne=rowDofsParamMapping%ELEMENT_DOF2PARAM_MAP(1,nyy)   ! current element (i.e. corresponds to current dof)
                                          nh=rowDofsParamMapping%ELEMENT_DOF2PARAM_MAP(2,nyy)   ! current variable component
                                          domainElements=>rowVariable%COMPONENTS(nh)%DOMAIN%TOPOLOGY%ELEMENTS
                                          basis=>domainElements%ELEMENTS(ne)%BASIS
                                          !Set up list
                                          NULLIFY(columnIndicesLists(local_ny)%PTR)
                                          CALL LIST_CREATE_START(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_DATA_TYPE_SET(columnIndicesLists(local_ny)%PTR,LIST_INTG_TYPE,err,error,*999)
                                          CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(local_ny)%PTR, &
                                            & rowVariable%COMPONENTS(nh)%maxNumberElementInterpolationParameters+1, &
                                            & err,error,*999) ! size = all nodal dofs + itself
                                          CALL LIST_CREATE_FINISH(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          DO nh2=1,fieldVariable%NUMBER_OF_COMPONENTS
                                            SELECT CASE(fieldVariable%COMPONENTS(nh2)%INTERPOLATION_TYPE)
                                            CASE(FIELD_CONSTANT_INTERPOLATION)
                                              CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
                                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                              ! it's assumed that element-based variables arne't directly coupled
                                              ! put a diagonal entry
                                              localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP% &
                                                & ELEMENT_PARAM2DOF_MAP%ELEMENTS(ne)
                                              globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                              CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn,err,error,*999)
                                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                                              ! loop over all nodes in the element (and dofs belonging to them)
                                              DO nn=1,basis%NUMBER_OF_NODES
                                                mp=domainElements%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                                DO nnk=1,basis%NUMBER_OF_DERIVATIVES(nn)
                                                  mk=domainElements%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                                  mv=domainElements%ELEMENTS(ne)%elementVersions(nnk,nn)
                                                  !Find the local and global column and add the global column to the indices list
                                                  localColumn=fieldVariable%COMPONENTS(nh2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                                    & NODES(mp)%DERIVATIVES(mk)%VERSIONS(mv)
                                                  globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)
                                                  CALL LIST_ITEM_ADD(columnIndicesLists(local_ny)%PTR,globalColumn, &
                                                    & err,error,*999)
                                                ENDDO !mk
                                              ENDDO !nn
                                            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                              CALL FlagError("Grid point based interpolation is not implemented yet.", &
                                                & err,error,*999)
                                            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                              CALL FlagError("Gauss point based interpolation is not implemented yet.", &
                                                & err,error,*999)
                                            CASE DEFAULT
                                              localError="Local dof number "//TRIM(NumberToVString(local_ny,"*",err,error))// &
                                                & " has invalid interpolation type."
                                              CALL FlagError(localError,err,error,*999)
                                            END SELECT
                                          ENDDO !nh2
                                          ! clean up the list
                                          CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(local_ny)%PTR,err,error,*999)
                                          CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(local_ny)%PTR,numberOfColumns, &
                                            & err,error,*999)
                                          numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                          rowIndices(local_ny+1)=numberOfNonZeros+1
                                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                          CALL FlagError("Grid point based interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                          CALL FlagError("Gauss point based interpolation is not implemented yet.",err,error,*999)
                                        CASE DEFAULT
                                          localError="Local dof number "//TRIM(NumberToVString(local_ny,"*",err,error))// &
                                            & " has an invalid type."
                                          CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      ENDDO !local_ny
                                      !Allocate and setup the column locations
                                      ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
                                      DO local_ny=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(local_ny)%PTR,numberOfColumns,columns, &
                                          & err,error,*999)
                                        DO columnIdx=1,numberOfColumns
                                          columnIndices(rowIndices(local_ny)+columnIdx-1)=columns(columnIdx)
                                        ENDDO !columnIdx
                                        DEALLOCATE(columns)
                                      ENDDO !local_ny
                                      IF(DIAGNOSTICS1) THEN
                                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Jacobian matrix structure:",err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                          & dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                          & dependentDofsDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                          & numberOfNonZeros,err,error,*999)
                                        IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                          & dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                          sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
                                            & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
                                          CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ", &
                                            & sparsity,"F6.2",err,error,*999)
                                        ENDIF
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,rowDofsDomainMapping% &
                                          & TOTAL_NUMBER_OF_LOCAL+1,8,8,rowIndices,'("  Row indices    :",8(X,I13))', &
                                          & '(18X,8(X,I13))',err,error,*999)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices,&
                                          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Row dofs parameter mapping is not associated.",err,error,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Row dofs domain mapping is not associated.",err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Dependent dofs parameter mapping is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Dependent dofs domain mapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Nonlinear matrices equations matrices is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations matrix nonlinear matrices is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NumberToVString(jacobianMatrix%STORAGE_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(EQUATIONS_MATRIX_NODAL_STRUCTURE)
            SELECT CASE(jacobianMatrix%STORAGE_TYPE)
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              nonlinearMatrices=>jacobianMatrix%NONLINEAR_MATRICES
              IF(ASSOCIATED(nonlinearMatrices)) THEN
                equationsMatrices=>nonlinearMatrices%EQUATIONS_MATRICES
                IF(ASSOCIATED(equationsMatrices)) THEN
                  equations=>equationsMatrices%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    equationsMapping=>equationsMatrices%EQUATIONS_MAPPING
                    IF(ASSOCIATED(equationsMapping)) THEN
                      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
                      IF(ASSOCIATED(nonlinearMapping)) THEN
                        equationsSet=>equations%EQUATIONS_SET
                        IF(ASSOCIATED(equationsSet)) THEN
                          dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            fieldVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(matrixNumber)%VARIABLE
                            IF(ASSOCIATED(fieldVariable)) THEN
                              dependentDofsDomainMapping=>fieldVariable%DOMAIN_MAPPING
                              IF(ASSOCIATED(dependentDofsDomainMapping)) THEN
                                dependentDofsParamMapping=>fieldVariable%DOF_TO_PARAM_MAP
                                IF(ASSOCIATED(dependentDofsParamMapping)) THEN
                                  !If RHS variable exists, use this for row DOFs, else use the first nonlinear variable
                                  IF(ASSOCIATED(equationsMapping%RHS_MAPPING)) THEN
                                    rowVariable=>equationsMapping%RHS_MAPPING%RHS_VARIABLE
                                  ELSE
                                    rowVariable=>nonlinearMapping%JACOBIAN_TO_VAR_MAP(1)%VARIABLE
                                  ENDIF
                                  IF(ASSOCIATED(rowVariable)) THEN
                                    rowDofsDomainMapping=>rowVariable%DOMAIN_MAPPING
                                    rowDofsParamMapping=>rowVariable%DOF_TO_PARAM_MAP
                                  ELSE
                                    CALL FlagError("RHS or first nonlinear variable is not associated",err,error,*999)
                                  ENDIF
                                  IF(ASSOCIATED(rowDofsDomainMapping)) THEN
                                    IF(ASSOCIATED(rowDofsParamMapping)) THEN
                                      !Allocate lists
                                      ALLOCATE(columnIndicesLists(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",err,error,*999)
                                      !Allocate row indices
                                      ALLOCATE(rowIndices(rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate row indices.",err,error,*999)
                                      rowIndices(1)=1
                                      !First, loop over the rows and calculate the number of non-zeros
                                      numberOfNonZeros=0
                                      DO localDofIdx=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        SELECT CASE(rowDofsParamMapping%DOF_TYPE(1,localDofIdx))
                                        CASE(FIELD_CONSTANT_INTERPOLATION)
                                          CALL FlagError("Constant interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_NODE_DOF_TYPE)
                                          dofIdx=dependentDofsParamMapping%DOF_TYPE(2,localDofIdx)!value for a particular field dof (localDofIdx)
                                          nodeIdx=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(3,dofIdx)!node number (np) of the field parameter
                                          componentIdx=dependentDofsParamMapping%NODE_DOF2PARAM_MAP(4,dofIdx)!component number (nh) of the field parameter
                                          domainNodes=>fieldVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%NODES

                                          !Set up list
                                          NULLIFY(columnIndicesLists(localDofIdx)%PTR)
                                          CALL LIST_CREATE_START(columnIndicesLists(localDofIdx)%PTR,err,error,*999)
                                          CALL LIST_DATA_TYPE_SET(columnIndicesLists(localDofIdx)%PTR,LIST_INTG_TYPE,err,error,*999)

                                          CALL LIST_INITIAL_SIZE_SET(columnIndicesLists(localDofIdx)%PTR, &
                                            & fieldVariable%NUMBER_OF_COMPONENTS* &
                                            & fieldVariable%maxNumberElementInterpolationParameters,err,error,*999)

                                          CALL LIST_CREATE_FINISH(columnIndicesLists(localDofIdx)%PTR,err,error,*999)
                                          !Loop over all components,nodes,derivatives, and versions
                                          DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
                                            SELECT CASE(fieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
                                            CASE(FIELD_NODE_BASED_INTERPOLATION)
                                              numberOfDerivatives=fieldVariable%components(componentIdx)%domain%topology% &
                                               & nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                                              DO derivativeIdx=1,numberOfDerivatives
                                                numberOfVersions=fieldVariable%components(componentIdx)%domain%topology% &
                                                 & nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                                                DO versionIdx=1,numberOfVersions
                                                  localColumn=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                                   & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                   & VERSIONS(versionIdx)
                                                  globalColumn=fieldVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localColumn)

                                                  CALL LIST_ITEM_ADD(columnIndicesLists(localDofIdx)%PTR,globalColumn, &
                                                    & err,error,*999)

                                                ENDDO !versionIdx
                                              ENDDO !derivativeIdx
                                            CASE DEFAULT
                                              localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))// &
                                                & " has invalid interpolation type."
                                              CALL FlagError(localError,err,error,*999)
                                            END SELECT
                                          ENDDO !componentIdx

                                          CALL LIST_REMOVE_DUPLICATES(columnIndicesLists(localDofIdx)%PTR,err,error,*999)
                                          CALL LIST_NUMBER_OF_ITEMS_GET(columnIndicesLists(localDofIdx)%PTR,numberOfColumns, &
                                            & err,error,*999)
                                          numberOfNonZeros=numberOfNonZeros+numberOfColumns
                                          rowIndices(localDofIdx+1)=numberOfNonZeros+1
                                        CASE(FIELD_ELEMENT_DOF_TYPE)
                                          CALL FlagError("Element based interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                          CALL FlagError("Grid point based interpolation is not implemented yet.",err,error,*999)
                                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                          CALL FlagError("Gauss point based interpolation is not implemented yet.",err,error,*999)
                                        CASE DEFAULT
                                          localError="Local dof number "//TRIM(NumberToVString(localDofIdx,"*",err,error))// &
                                            & " has an invalid type."
                                          CALL FlagError(localError,err,error,*999)
                                        END SELECT
                                      ENDDO !localDofIdx
                                      !Allocate and setup the column locations
                                      ALLOCATE(columnIndices(numberOfNonZeros),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate column indices.",err,error,*999)
                                      DO localDofIdx=1,rowDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
                                        CALL LIST_DETACH_AND_DESTROY(columnIndicesLists(localDofIdx)%PTR,numberOfColumns,columns, &
                                          & err,error,*999)
                                        DO columnIdx=1,numberOfColumns
                                          columnIndices(rowIndices(localDofIdx)+columnIdx-1)=columns(columnIdx)
                                        ENDDO !columnIdx
                                        DEALLOCATE(columns)
                                      ENDDO !localDofIdx
                                      IF(DIAGNOSTICS1) THEN
                                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Jacobian matrix structure :",err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                          & dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL,err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                          & dependentDofsDomainMapping%NUMBER_OF_GLOBAL,err,error,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                          & numberOfNonZeros,err,error,*999)
                                        IF(dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL* &
                                          & dependentDofsDomainMapping%NUMBER_OF_GLOBAL/=0) THEN
                                          sparsity=(1.0_DP-REAL(numberOfNonZeros,DP)/REAL(dependentDofsDomainMapping% &
                                            & TOTAL_NUMBER_OF_LOCAL*dependentDofsDomainMapping%NUMBER_OF_GLOBAL,DP))*100.0_DP
                                          CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (% of zeros) = ", &
                                            & sparsity,"F6.2",err,error,*999)
                                        ENDIF
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,rowDofsDomainMapping% &
                                          & TOTAL_NUMBER_OF_LOCAL+1,8,8,rowIndices,'("  Row indices    :",8(X,I13))', &
                                          & '(18X,8(X,I13))',err,error,*999)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfNonZeros,8,8,columnIndices,&
                                          & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', err,error,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Row dofs parameter mapping is not associated.",err,error,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Row dofs domain mapping is not associated.",err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Dependent dofs parameter mapping is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Dependent dofs domain mapping is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations mapping nonlinear mapping is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Nonlinear matrices equations matrices is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations matrix nonlinear matrices is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NumberToVString(jacobianMatrix%STORAGE_TYPE,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The matrix structure type of "// &
              & TRIM(NumberToVString(jacobianMatrix%STRUCTURE_TYPE,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*998)
          END SELECT
        ELSE
          CALL FlagError("Column indices is already associated.",err,error,*998)
        ENDIF
      ELSE
        CALL FlagError("Row indices is already associated.",err,error,*998)
      ENDIF
    ELSE
      CALL FlagError("Jacobian matrix is not associated.",err,error,*999)
    ENDIF
      
    EXITS("JacobianMatrix_StructureCalculate")
    RETURN
999 IF(ASSOCIATED(rowIndices)) DEALLOCATE(rowIndices)
    IF(ASSOCIATED(columnIndices)) DEALLOCATE(columnIndices)
    IF(ALLOCATED(columns)) DEALLOCATE(columns)
    IF(ALLOCATED(columnIndicesLists)) THEN
      DO localDofIdx=1,dependentDofsDomainMapping%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(columnIndicesLists(localDofIdx)%PTR)) &
          & CALL LIST_DESTROY(columnIndicesLists(localDofIdx)%PTR,dummyErr,dummyError,*998)
      ENDDO !localDofIdx
      DEALLOCATE(columnIndicesLists)
    ENDIF
998 ERRORSEXITS("JacobianMatrix_StructureCalculate",err,error)
    RETURN 1
  END SUBROUTINE JacobianMatrix_StructureCalculate

  !
  !================================================================================================================================
  !
 
END MODULE EQUATIONS_MATRICES_ROUTINES
