!> \file
!> \author Chris Bradley
!> \brief This module contains all interface matrices routines.
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

!>This module contains all interface matrices routines.
MODULE INTERFACE_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_MATRICES_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup INTERFACE_MATRICES_ROUTINES_InterfaceMatrixStructureTypes INTERFACE_MATRICES_ROUTINES::InterfaceMatrixStructureTypes
  !> \brief Interface matrices structure (sparsity) types
  !> \see INTERFACE_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see INTERFACE_MATRICES_ROUTINES_InterfaceMatrixStructureTypes,INTERFACE_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see INTERFACE_MATRICES_ROUTINES_InterfaceMatrixStructureTypes,INTERFACE_MATRICES_ROUTINES 
  !>@}

  !> \addtogroup INTERFACE_MATRICES_ROUTINES_InterfaceMatricesSparsityTypes INTERFACE_MATRICES_ROUTINES::InterfaceMatricesSparsityTypes
  !> \brief Interface matrices sparsity types
  !> \see INTERFACE_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRICES_SPARSE_MATRICES=1 !<Use sparse interface matrices \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesSparsityTypes,INTERFACE_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRICES_FULL_MATRICES=2 !<Use fully populated interface matrices \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesSparsityTypes,INTERFACE_MATRICES_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC INTERFACE_MATRIX_NO_STRUCTURE,INTERFACE_MATRIX_FEM_STRUCTURE

  PUBLIC INTERFACE_MATRICES_SPARSE_MATRICES,INTERFACE_MATRICES_FULL_MATRICES

  PUBLIC INTERFACE_MATRICES_CREATE_FINISH,INTERFACE_MATRICES_CREATE_START

  PUBLIC INTERFACE_MATRICES_DESTROY

  PUBLIC INTERFACE_MATRICES_ELEMENT_ADD

  PUBLIC InterfaceMatrices_ElementCalculate

  PUBLIC INTERFACE_MATRICES_ELEMENT_FINALISE,InterfaceMatrices_ElementInitialise

  PUBLIC INTERFACE_MATRICES_OUTPUT

  PUBLIC INTERFACE_MATRICES_STORAGE_TYPE_SET

  PUBLIC INTERFACE_MATRICES_STRUCTURE_TYPE_SET

  PUBLIC INTERFACE_MATRICES_VALUES_INITIALISE
  
  PUBLIC InterfaceMatrix_TimeDependenceTypeSet,InterfaceMatrix_TimeDependenceTypeGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalise a interface matrix and deallocate all memory
  SUBROUTINE INTERFACE_MATRIX_FINALISE(INTERFACE_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("INTERFACE_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
      IF(ASSOCIATED(INTERFACE_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE_MATRIX%MATRIX_TRANSPOSE)) CALL DISTRIBUTED_MATRIX_DESTROY(INTERFACE_MATRIX%MATRIX_TRANSPOSE, &
        & ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_MATRIX)
    ENDIF
    
    EXITS("INTERFACE_MATRIX_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRIX_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an interface matrix.
  SUBROUTINE INTERFACE_MATRIX_INITIALISE(INTERFACE_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices to initialise the interface matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The matrix number in the interface matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("INTERFACE_MATRIX_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
        INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          IF(ASSOCIATED(INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
            LOCAL_ERROR="Interface matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
              & " is already associated."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
          ELSE
            ALLOCATE(INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface matrix.",ERR,ERROR,*999)
            INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
            INTERFACE_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
            INTERFACE_MATRIX%INTERFACE_MATRICES=>INTERFACE_MATRICES
            INTERFACE_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
            INTERFACE_MATRIX%STRUCTURE_TYPE=INTERFACE_MATRIX_NO_STRUCTURE
            INTERFACE_MATRIX%UPDATE_MATRIX=.TRUE.
            INTERFACE_MATRIX%FIRST_ASSEMBLY=.TRUE.
            INTERFACE_MATRIX%HAS_TRANSPOSE=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%HAS_TRANSPOSE
            INTERFACE_MATRIX%NUMBER_OF_ROWS=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_ROWS
            INTERFACE_MATRIX%TOTAL_NUMBER_OF_ROWS=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)% &
              & TOTAL_NUMBER_OF_ROWS
            INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%INTERFACE_MATRIX=>INTERFACE_MATRIX
            NULLIFY(INTERFACE_MATRIX%MATRIX)
            NULLIFY(INTERFACE_MATRIX%MATRIX_TRANSPOSE)
            NULLIFY(INTERFACE_MATRIX%TEMP_VECTOR)
            NULLIFY(INTERFACE_MATRIX%TEMP_TRANSPOSE_VECTOR)
            CALL EquationsMatrices_ElementMatrixInitialise(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified interface matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//"."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("INTERFACE_MATRIX_INITIALISE")
    RETURN
999 CALL INTERFACE_MATRIX_FINALISE(INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_MATRIX_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds the element matrices into the interface matrices.
  SUBROUTINE INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_ADD()")
#endif

    ENTERS("INTERFACE_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      !Add the element matrices
      DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
        INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
          IF(INTERFACE_MATRIX%UPDATE_MATRIX) THEN
            !Add the element matrix into the distributed interface equations matrix
            CALL DISTRIBUTED_MATRIX_VALUES_ADD(INTERFACE_MATRIX%MATRIX,INTERFACE_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
              & INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),INTERFACE_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
              & INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
              & INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
              & ERR,ERROR,*999)
            !If the interface matrix has a transpose add it
            IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_VALUES_ADD(INTERFACE_MATRIX%MATRIX_TRANSPOSE,INTERFACE_MATRIX%ELEMENT_MATRIX%COLUMN_DOFS(1: &
                & INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS),INTERFACE_MATRIX%ELEMENT_MATRIX%ROW_DOFS(1: &
                & INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS),TRANSPOSE(INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(1: &
                & INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:INTERFACE_MATRIX%ELEMENT_MATRIX%NUMBER_OF_COLUMNS)), &
                & ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="Interface matrix for interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      RHS_VECTOR=>INTERFACE_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        IF(RHS_VECTOR%UPDATE_VECTOR) THEN
          !Add the rhs element vector
          CALL DISTRIBUTED_VECTOR_VALUES_ADD(RHS_VECTOR%RHS_VECTOR,RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS(1: &
            & RHS_VECTOR%ELEMENT_VECTOR%NUMBER_OF_ROWS),RHS_VECTOR%ELEMENT_VECTOR%VECTOR(1:RHS_VECTOR% &
            & ELEMENT_VECTOR%NUMBER_OF_ROWS),ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_ADD()")
#endif
    
    EXITS("INTERFACE_MATRICES_ELEMENT_ADD")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_ELEMENT_ADD",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Calculate the positions of the element matrices in the interface matrices. 
  SUBROUTINE InterfaceMatrices_ElementCalculate(interfaceMatrices,interfaceElementNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: interfaceMatrices !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,rowsElementNumber,rowsMeshIdx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable,rowsFieldVariable
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: interfaceMapping
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: interfaceMatrix
    TYPE(INTERFACE_RHS_TYPE), POINTER :: rhsVector
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: meshConnectivity
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(VARYING_STRING) :: localError

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("InterfaceMatrices_ElementCalculate()")
#endif

    ENTERS("InterfaceMatrices_ElementCalculate",err,error,*999)

    IF(ASSOCIATED(interfaceMatrices)) THEN
      interfaceMapping=>interfaceMatrices%INTERFACE_MAPPING
      IF(ASSOCIATED(interfaceMapping)) THEN
        interfaceEquations=>interfaceMapping%INTERFACE_EQUATIONS
        IF(ASSOCIATED(interfaceEquations)) THEN
          interfaceCondition=>interfaceEquations%INTERFACE_CONDITION
          IF(ASSOCIATED(interfaceCondition)) THEN
            interface=>interfaceCondition%INTERFACE
            IF(ASSOCIATED(interface)) THEN
              SELECT CASE(interfaceCondition%integrationType)
              CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
                meshConnectivity=>interface%MESH_CONNECTIVITY
                IF(ASSOCIATED(meshConnectivity)) THEN
                  IF(ALLOCATED(meshConnectivity%ELEMENT_CONNECTIVITY)) THEN
                    !Calculate the row and columns for the interface equations matrices
                    DO matrixIdx=1,interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES
                      interfaceMatrix=>interfaceMatrices%MATRICES(matrixIdx)%PTR
                      IF(ASSOCIATED(interfaceMatrix)) THEN
                        rowsFieldVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                        colsFieldVariable=>interfaceMapping%LAGRANGE_VARIABLE !\todo: TEMPORARY: Needs generalising
                        rowsMeshIdx=interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%MESH_INDEX
                        IF(ASSOCIATED(rowsFieldVariable,colsFieldVariable)) THEN
                          ! If the rows and column variables are both the Lagrange variable (this is the diagonal matrix)
                          rowsElementNumber=InterfaceElementNumber
                        ELSE
                          rowsElementNumber=meshConnectivity%ELEMENT_CONNECTIVITY(InterfaceElementNumber,rowsMeshIdx)% &
                            & COUPLED_MESH_ELEMENT_NUMBER
                        ENDIF
                        CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(interfaceMatrix%ELEMENT_MATRIX, &
                          & interfaceMatrix%UPDATE_MATRIX,[rowsElementNumber],[interfaceElementNumber],rowsFieldVariable, &
                          & colsFieldVariable,err,error,*999)
                      ELSE
                        localError="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                          & " is not associated."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO !matrixIdx
                  ELSE
                    CALL FlagError("Interface element connectivity is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Interface mesh connectivity is not associated.",err,error,*999)              
                ENDIF
              CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
                pointsConnectivity=>interface%pointsConnectivity
                IF(ASSOCIATED(pointsConnectivity)) THEN
                  IF(ALLOCATED(pointsConnectivity%coupledElements)) THEN
                    DO matrixIdx=1,interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES
                      interfaceMatrix=>interfaceMatrices%MATRICES(matrixIdx)%PTR
                      IF(ASSOCIATED(interfaceMatrix)) THEN 
                        IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                            matrixIdx==interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES) THEN
                          rowsFieldVariable=>interfaceMapping%LAGRANGE_VARIABLE
                          colsFieldVariable=>interfaceMapping%LAGRANGE_VARIABLE
                          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(interfaceMatrix%ELEMENT_MATRIX, &
                            & interfaceMatrix%UPDATE_MATRIX,[InterfaceElementNumber],[InterfaceElementNumber], &
                            & rowsFieldVariable,colsFieldVariable,err,error,*999)
                        ELSE
                          rowsFieldVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                          colsFieldVariable=>interfaceMapping%LAGRANGE_VARIABLE !\todo: TEMPORARY: Needs generalising
                          rowsMeshIdx=interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%MESH_INDEX
                          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(interfaceMatrix%ELEMENT_MATRIX, &
                            & interfaceMatrix%UPDATE_MATRIX,pointsConnectivity%coupledElements(InterfaceElementNumber, &
                            & rowsMeshIdx)%elementNumbers,[InterfaceElementNumber],rowsFieldVariable,colsFieldVariable, &
                            & err,error,*999)
                        ENDIF
                      ELSE
                        localError="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                          & " is not associated."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO !matrixIdx
                  ELSE
                    CALL FlagError("Interface points connectivity coupled elements is not allocated.",err,error,*999)             
                  ENDIF
                ELSE
                  CALL FlagError("Interface points connectivity is not associated.",err,error,*999)              
                ENDIF
              CASE DEFAULT
                localError="The interface condition integration type of "// &
                  & TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(localError,ERR,ERROR,*999)
              END SELECT
              !RHS element matrix dofs are the same for both mesh and points connectivity, right now
              rhsVector=>interfaceMatrices%RHS_VECTOR
              IF(ASSOCIATED(rhsVector)) THEN
                rhsMapping=>interfaceMapping%RHS_MAPPING
                IF(ASSOCIATED(rhsMapping)) THEN
                  !Calculate the rows  for the equations RHS
                  rowsFieldVariable=>rhsMapping%RHS_VARIABLE
                  CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_CALCULATE(rhsVector%ELEMENT_VECTOR,rhsVector%UPDATE_VECTOR, &
                    & interfaceElementNumber,rowsFieldVariable,err,error,*999)
                ELSE
                  CALL FlagError("Interface mapping rhs mapping is not associated.",err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Interface condition interface is not associated.",err,error,*999)            
            ENDIF
          ELSE
            CALL FlagError("Interface equations interface condition is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mapping interface equations is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface mapping is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not allocated",err,error,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("InterfaceMatrices_ElementCalculate()")
#endif
    
    EXITS("InterfaceMatrices_ElementCalculate")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ElementCalculate",err,error)
    RETURN 1
  END SUBROUTINE InterfaceMatrices_ElementCalculate

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information for interface matrices and deallocate all memory
  SUBROUTINE INTERFACE_MATRICES_ELEMENT_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<The interface matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("INTERFACE_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
        INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is not associated."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      RHS_VECTOR=>INTERFACE_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        !Finalise the interface element vector
        RHS_VECTOR%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%ROW_DOFS)
        IF(ALLOCATED(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(RHS_VECTOR%ELEMENT_VECTOR%VECTOR)
      ENDIF

      ENDDO !matrix_idx
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MATRICES_ELEMENT_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_ELEMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the interface matrices
  SUBROUTINE InterfaceMatrices_ElementInitialise(interfaceMatrices,err,error,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: interfaceMatrices !The interface matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,rowsMeshIdx
    INTEGER(INTG) :: rowsNumberOfElements,colsNumberOfElements !Number of elements in the row and col variables whose dofs are present in interface element matrix
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: interfaceMapping
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations
    TYPE(INTERFACE_TYPE), POINTER :: interface
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: interfaceMatrix
    TYPE(INTERFACE_RHS_TYPE), POINTER :: rhsVector
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: rhsMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: colsFieldVariable,rowsFieldVariable
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("InterfaceMatrices_ElementInitialise",err,error,*999)

    IF(ASSOCIATED(interfaceMatrices)) THEN
      interfaceMapping=>interfaceMatrices%INTERFACE_MAPPING
      IF(ASSOCIATED(interfaceMapping)) THEN
        interfaceEquations=>interfaceMapping%INTERFACE_EQUATIONS
        IF(ASSOCIATED(interfaceEquations)) THEN
          interfaceCondition=>interfaceEquations%INTERFACE_CONDITION
          IF(ASSOCIATED(interfaceCondition)) THEN
            SELECT CASE(interfaceCondition%integrationType)
              CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
              DO matrixIdx=1,interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES
                interfaceMatrix=>interfaceMatrices%MATRICES(matrixIdx)%PTR
                IF(ASSOCIATED(interfaceMatrix)) THEN
                  rowsNumberOfElements=1
                  colsNumberOfElements=1
                  rowsFieldVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                  colsFieldVariable=>interfaceMapping%LAGRANGE_VARIABLE !TEMPORARY: Needs generalising
                  CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(interfaceMatrix%ELEMENT_MATRIX,rowsFieldVariable, &
                    & colsFieldVariable,rowsNumberOfElements,colsNumberOfElements,err,error,*999)
                ELSE
                  localError="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                    & " is not associated."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !matrixIdx
            CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION) 
              interface=>interfaceCondition%INTERFACE
              IF(ASSOCIATED(interface))THEN
                pointsConnectivity=>interface%pointsConnectivity
                IF(ASSOCIATED(pointsConnectivity)) THEN
                  IF(ALLOCATED(pointsConnectivity%coupledElements)) THEN
                    DO matrixIdx=1,interfaceMatrices%NUMBER_OF_INTERFACE_MATRICES !\todo: Need to separate the case for penalty matrix
                      interfaceMatrix=>interfaceMatrices%MATRICES(matrixIdx)%PTR
                      IF(ASSOCIATED(interfaceMatrix)) THEN
                        colsNumberOfElements=1
                        rowsFieldVariable=>interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%VARIABLE
                        colsFieldVariable=>interfaceMapping%LAGRANGE_VARIABLE !TEMPORARY: Needs generalising
                        rowsMeshIdx=interfaceMapping%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrixIdx)%MESH_INDEX
                        CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(interfaceMatrix%ELEMENT_MATRIX,rowsFieldVariable, &
                          & colsFieldVariable,pointsConnectivity%maxNumberOfCoupledElements(rowsMeshIdx), &
                          & colsNumberOfElements,err,error,*999)
                      ELSE
                        localError="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrixIdx,"*",err,error))// &
                          & " is not associated."
                        CALL FlagError(localError,err,error,*999)
                      ENDIF
                    ENDDO !matrixIdx
                  ELSE
                    CALL FlagError("Interface points connectivity coupled elements is not allocated.",err,error,*999)             
                  ENDIF
                ELSE
                  CALL FlagError("Interface points connectivity is not associated.",err,error,*999)              
                ENDIF
              ELSE
                CALL FlagError("Interface is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The interface condition integration type of "// &
                & TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType,"*",ERR,ERROR))//" is invalid."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Interface condition is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface equations is not associated.",err,error,*999)
        ENDIF
        rhsVector=>interfaceMatrices%RHS_VECTOR
        IF(ASSOCIATED(rhsVector)) THEN
          !Initialise the RHS element vector
          rhsMapping=>interfaceMapping%RHS_MAPPING
          IF(ASSOCIATED(rhsMapping)) THEN
            rowsFieldVariable=>rhsMapping%RHS_VARIABLE
            CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_SETUP(rhsVector%ELEMENT_VECTOR,rowsFieldVariable,err,error,*999)
          ELSE
            CALL FlagError("RHS mapping is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Interface matrices mapping is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",err,error,*999)
    ENDIF
    
    EXITS("InterfaceMatrices_ElementInitialise")
    RETURN
999 ERRORSEXITS("InterfaceMatrices_ElementInitialise",err,error)
    RETURN 1
  END SUBROUTINE InterfaceMatrices_ElementInitialise

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for an interface matrix.
  SUBROUTINE INTERFACE_MATRIX_STRUCTURE_CALCULATE(INTERFACE_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
    & TRANSPOSE_ROW_INDICES,TRANSPOSE_COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return, the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return, a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return, a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: TRANSPOSE_ROW_INDICES(:) !<On return, if the interface matrix has a transpose a pointer to transpose row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: TRANSPOSE_COLUMN_INDICES(:) !<On return, if the interface matrix has a transpose a pointer to the transpose column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   INTEGER(INTG) :: column_version,column_derivative,column_idx,column_component_idx,column_local_derivative_idx, &
     & column_local_node_idx, column_node,DUMMY_ERR,domain_element,global_column,global_row,interface_element_idx, &
     & INTERFACE_MESH_INDEX,local_column,local_row,MATRIX_NUMBER,NUMBER_OF_COLUMNS,NUMBER_OF_ROWS,row_component_idx, &
     & row_version,row_derivative,row_local_derivative_idx,row_idx,row_local_node_idx,row_node,TRANSPOSE_NUMBER_OF_NON_ZEROS
    INTEGER(INTG), ALLOCATABLE :: COLUMNS(:),TRANSPOSE_COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: COLUMN_BASIS,ROW_BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOFS_DOMAIN_MAPPING,ROW_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: COLUMN_DOMAIN_ELEMENTS,ROW_DOMAIN_ELEMENTS
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: MESH_CONNECTIVITY
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: COLUMN_DOFS_PARAM_MAPPING,ROW_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLUMN_VARIABLE,ROW_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: TRANSPOSE_COLUMN_INDICES_LISTS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    ENTERS("INTERFACE_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*999)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
      IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
        IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
          IF(.NOT.ASSOCIATED(TRANSPOSE_ROW_INDICES)) THEN
            IF(.NOT.ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) THEN
              MATRIX_NUMBER=INTERFACE_MATRIX%MATRIX_NUMBER
              SELECT CASE(INTERFACE_MATRIX%STRUCTURE_TYPE)
              CASE(INTERFACE_MATRIX_NO_STRUCTURE)
                CALL FlagError("There is no structure to calculate for a matrix with no structure.",ERR,ERROR,*998)
              CASE(INTERFACE_MATRIX_FEM_STRUCTURE)
                SELECT CASE(INTERFACE_MATRIX%STORAGE_TYPE)
                CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                  INTERFACE_MATRICES=>INTERFACE_MATRIX%INTERFACE_MATRICES
                  IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
                    INTERFACE_EQUATIONS=>INTERFACE_MATRICES%INTERFACE_EQUATIONS
                    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                      INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
                      IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
                        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                          INTERFACE=>INTERFACE_CONDITION%INTERFACE
                          IF(ASSOCIATED(INTERFACE)) THEN
                            MESH_CONNECTIVITY=>INTERFACE%MESH_CONNECTIVITY
                            IF(ASSOCIATED(MESH_CONNECTIVITY)) THEN
                              ROW_VARIABLE=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%VARIABLE
                              INTERFACE_MESH_INDEX=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%MESH_INDEX
                              IF(ASSOCIATED(ROW_VARIABLE)) THEN
                                COLUMN_VARIABLE=>INTERFACE_MAPPING%LAGRANGE_VARIABLE
                                IF(ASSOCIATED(COLUMN_VARIABLE)) THEN
                                  ROW_DOFS_DOMAIN_MAPPING=>ROW_VARIABLE%DOMAIN_MAPPING
                                  IF(ASSOCIATED(ROW_DOFS_DOMAIN_MAPPING)) THEN
                                    COLUMN_DOFS_DOMAIN_MAPPING=>COLUMN_VARIABLE%DOMAIN_MAPPING
                                    IF(ASSOCIATED(COLUMN_DOFS_DOMAIN_MAPPING)) THEN
                                      ROW_DOFS_PARAM_MAPPING=>ROW_VARIABLE%DOF_TO_PARAM_MAP
                                      IF(ASSOCIATED(ROW_DOFS_PARAM_MAPPING)) THEN
                                        COLUMN_DOFS_PARAM_MAPPING=>COLUMN_VARIABLE%DOF_TO_PARAM_MAP
                                        IF(ASSOCIATED(COLUMN_DOFS_PARAM_MAPPING)) THEN
                                          !Allocate lists
                                          ALLOCATE(COLUMN_INDICES_LISTS(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                          IF(ERR/=0) CALL FlagError("Could not allocate column indices lists.",ERR,ERROR,*999)
                                          DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            !Set up list
                                            NULLIFY(COLUMN_INDICES_LISTS(local_row)%PTR)
                                            CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                            CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_row)%PTR,LIST_INTG_TYPE, &
                                              & ERR,ERROR,*999)
                                            CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_row)%PTR,50,ERR,ERROR,*999)
                                            CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                          ENDDO !local_row
                                          !Allocate row indices
                                          ALLOCATE(ROW_INDICES(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                          IF(ERR/=0) CALL FlagError("Could not allocate row indices.",ERR,ERROR,*999)
                                          IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                            !Allocate transpose lists
                                            ALLOCATE(TRANSPOSE_COLUMN_INDICES_LISTS(COLUMN_DOFS_DOMAIN_MAPPING% &
                                              & TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                                            IF(ERR/=0) CALL FlagError("Could not allocate transpose column indices lists.", &
                                              & ERR,ERROR,*999)
                                            DO local_column=1,COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                              !Set up list
                                              NULLIFY(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR)
                                              CALL LIST_CREATE_START(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & ERR,ERROR,*999)
                                              CALL LIST_DATA_TYPE_SET(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & LIST_INTG_TYPE,ERR,ERROR,*999)
                                              CALL LIST_INITIAL_SIZE_SET(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR,50, &
                                                & ERR,ERROR,*999)
                                              CALL LIST_CREATE_FINISH(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & ERR,ERROR,*999)
                                            ENDDO !local_column
                                            !Allocate transpose row indices
                                            ALLOCATE(TRANSPOSE_ROW_INDICES(COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1), &
                                              & STAT=ERR)
                                            IF(ERR/=0) CALL FlagError("Could not allocate transpose row indices.",ERR,ERROR,*999)
                                          ENDIF
                                          !Loop over the number of components in the Lagrange multipler variable
                                          DO column_component_idx=1,COLUMN_VARIABLE%NUMBER_OF_COMPONENTS
                                            IF(COLUMN_VARIABLE%COMPONENTS(column_component_idx)%INTERPOLATION_TYPE== &
                                              & FIELD_NODE_BASED_INTERPOLATION) THEN
                                              !Loop over the elements in the interface mesh
                                              COLUMN_DOMAIN_ELEMENTS=>COLUMN_VARIABLE%COMPONENTS(column_component_idx)%DOMAIN% &
                                                & TOPOLOGY%ELEMENTS
                                              DO interface_element_idx=1,COLUMN_DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                                                COLUMN_BASIS=>COLUMN_DOMAIN_ELEMENTS%ELEMENTS(interface_element_idx)%BASIS
                                                !Loop over the column DOFs in the element
                                                DO column_local_node_idx=1,COLUMN_BASIS%NUMBER_OF_NODES
                                                  column_node=COLUMN_DOMAIN_ELEMENTS%ELEMENTS(interface_element_idx)% &
                                                    & ELEMENT_NODES(column_local_node_idx)
                                                  DO column_local_derivative_idx=1,COLUMN_BASIS% &
                                                    & NUMBER_OF_DERIVATIVES(column_local_node_idx)
                                                    column_derivative=COLUMN_DOMAIN_ELEMENTS%ELEMENTS(interface_element_idx)% &
                                                      & ELEMENT_DERIVATIVES(column_local_derivative_idx,column_local_node_idx)
                                                    column_version=COLUMN_DOMAIN_ELEMENTS%ELEMENTS(interface_element_idx)% &
                                                      & elementVersions(column_local_derivative_idx,column_local_node_idx)
                                                    local_column=COLUMN_VARIABLE%COMPONENTS(column_component_idx)% &
                                                      & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(column_node)% &
                                                      & DERIVATIVES(column_derivative)%VERSIONS(column_version)
                                                    global_column=COLUMN_DOFS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                                    !Loop over the components in the dependent variable
                                                    DO row_component_idx=1,ROW_VARIABLE%NUMBER_OF_COMPONENTS
                                                      SELECT CASE(ROW_VARIABLE%COMPONENTS(row_component_idx)%INTERPOLATION_TYPE)
                                                      CASE(FIELD_CONSTANT_INTERPOLATION)
                                                        local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                          & CONSTANT_PARAM2DOF_MAP
                                                        CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                          & ERR,ERROR,*999)
                                                        IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                                          global_row=ROW_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                          CALL LIST_ITEM_ADD(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                            & global_row,ERR,ERROR,*999)
                                                        ENDIF
                                                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                                        domain_element=MESH_CONNECTIVITY% &
                                                          & ELEMENT_CONNECTIVITY(interface_element_idx,INTERFACE_MESH_INDEX)% &
                                                          & COUPLED_MESH_ELEMENT_NUMBER
                                                        local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                          & ELEMENT_PARAM2DOF_MAP%ELEMENTS(domain_element)
                                                        CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                          & ERR,ERROR,*999)
                                                        IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                                          global_row=ROW_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                          CALL LIST_ITEM_ADD(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                            & global_row,ERR,ERROR,*999)
                                                        ENDIF
                                                      CASE(FIELD_NODE_BASED_INTERPOLATION)
                                                        ROW_DOMAIN_ELEMENTS=>ROW_VARIABLE%COMPONENTS(row_component_idx)%DOMAIN% &
                                                          & TOPOLOGY%ELEMENTS
                                                        domain_element=MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY( &
                                                          & interface_element_idx,INTERFACE_MESH_INDEX)%COUPLED_MESH_ELEMENT_NUMBER
                                                        ROW_BASIS=>ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)%BASIS
                                                        !Loop over the row DOFs in the domain mesh element
                                                        DO row_local_node_idx=1,ROW_BASIS%NUMBER_OF_NODES
                                                          row_node=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                            & ELEMENT_NODES(row_local_node_idx)
                                                          DO row_local_derivative_idx=1,ROW_BASIS% &
                                                            & NUMBER_OF_DERIVATIVES(row_local_node_idx)
                                                            row_derivative=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                              & ELEMENT_DERIVATIVES(row_local_derivative_idx,row_local_node_idx)
                                                            row_version=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                              & elementVersions(row_local_derivative_idx,row_local_node_idx)
                                                            local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)% &
                                                              & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(row_node)% &
                                                              & DERIVATIVES(row_derivative)%VERSIONS(row_version)
                                                            CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                              & ERR,ERROR,*999)
                                                            IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                                              global_row=ROW_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_row)
                                                              CALL LIST_ITEM_ADD(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                                & global_row,ERR,ERROR,*999)
                                                            ENDIF
                                                          ENDDO !row_local_derivative_idx
                                                        ENDDO !row_local_node_idx
                                                      CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                                        CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                      CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                                        CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                      CASE DEFAULT
                                                        LOCAL_ERROR="The row variable interpolation type of "// &
                                                          & TRIM(NUMBER_TO_VSTRING(ROW_VARIABLE%COMPONENTS(row_component_idx)% &
                                                          INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid."
                                                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                                      END SELECT
                                                    ENDDO !row_component_idx
                                                  ENDDO !column_local_derivative_idx
                                                ENDDO !column_local_node_idx
                                              ENDDO !interface_element_idx
                                            ELSE
                                              CALL FlagError("Only node based fields implemented.",ERR,ERROR,*999)
                                            ENDIF
                                          ENDDO !column_component_idx
                                          ROW_INDICES(1)=1
                                          DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                            CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                              & ERR,ERROR,*999)
                                            NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                            ROW_INDICES(local_row+1)=NUMBER_OF_NON_ZEROS+1
                                          ENDDO !local_row
                                          IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                            TRANSPOSE_NUMBER_OF_NON_ZEROS=0
                                            TRANSPOSE_ROW_INDICES(1)=1
                                            DO local_column=1,COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                              CALL LIST_REMOVE_DUPLICATES(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & ERR,ERROR,*999)
                                              CALL LIST_NUMBER_OF_ITEMS_GET(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                                              TRANSPOSE_NUMBER_OF_NON_ZEROS=TRANSPOSE_NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                              TRANSPOSE_ROW_INDICES(local_column+1)=TRANSPOSE_NUMBER_OF_NON_ZEROS+1
                                            ENDDO !local_column
                                            !Sanity check - the number of non-zeros should be the same
                                            IF(TRANSPOSE_NUMBER_OF_NON_ZEROS/=NUMBER_OF_NON_ZEROS) THEN
                                              LOCAL_ERROR="Invalid number of non-zeros. The number of non-zeros in the "// &
                                                & "transposed matrix ("//TRIM(NUMBER_TO_VSTRING(TRANSPOSE_NUMBER_OF_NON_ZEROS, &
                                                & "*",ERR,ERROR))//") does not match the number of non-zeros in the interface "// &
                                                & "matrix ("//TRIM(NUMBER_TO_VSTRING(TRANSPOSE_NUMBER_OF_NON_ZEROS,"*",ERR, &
                                                & ERROR))//")."
                                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                            ENDIF
                                          ENDIF
                                          !Allocate and setup the column locations
                                          ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                          IF(ERR/=0) CALL FlagError("Could not allocate column indices.",ERR,ERROR,*999)
                                          DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                              & COLUMNS,ERR,ERROR,*999)
                                            DO column_idx=1,NUMBER_OF_COLUMNS
                                              COLUMN_INDICES(ROW_INDICES(local_row)+column_idx-1)=COLUMNS(column_idx)
                                            ENDDO !column_idx
                                            DEALLOCATE(COLUMNS)
                                          ENDDO !local_row
                                          IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                            !Allocate and setup the column locations
                                            ALLOCATE(TRANSPOSE_COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                            IF(ERR/=0) &
                                              & CALL FlagError("Could not allocate transpose column indices.",ERR,ERROR,*999)
                                            DO local_column=1,COLUMN_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                              CALL LIST_DETACH_AND_DESTROY(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR, &
                                                & NUMBER_OF_ROWS,TRANSPOSE_COLUMNS,ERR,ERROR,*999)
                                              DO row_idx=1,NUMBER_OF_ROWS
                                                TRANSPOSE_COLUMN_INDICES(TRANSPOSE_ROW_INDICES(local_column)+row_idx-1)= &
                                                  & TRANSPOSE_COLUMNS(row_idx)
                                              ENDDO !row_idx
                                              DEALLOCATE(TRANSPOSE_COLUMNS)
                                            ENDDO !local_column
                                          ENDIF
                                          IF(DIAGNOSTICS1) THEN
                                            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interface matrix structure:",ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Interface matrix number : ", &
                                              & MATRIX_NUMBER,ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                              & ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                              & COLUMN_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                                            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ", &
                                              & NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                                            IF(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                                              & COLUMN_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL/=0) THEN
                                              SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(ROW_DOFS_DOMAIN_MAPPING% &
                                                & TOTAL_NUMBER_OF_LOCAL*COLUMN_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,DP)*100.0_DP
                                              CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY, &
                                                & "F6.2",ERR,ERROR,*999)
                                            ENDIF
                                            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ROW_DOFS_DOMAIN_MAPPING% &
                                              & TOTAL_NUMBER_OF_LOCAL+1,5,5,ROW_INDICES, &
                                              & '("  Row indices              :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
                                            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                              & COLUMN_INDICES,'("  Column indices           :",5(X,I13))','(28X,5(X,I13))', &
                                              & ERR,ERROR,*999)
                                            IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN 
                                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,COLUMN_DOFS_DOMAIN_MAPPING% &
                                                & TOTAL_NUMBER_OF_LOCAL+1,5,5,TRANSPOSE_ROW_INDICES, &
                                                & '("  Transpose row indices    :",5(X,I13))','(28X,5(X,I13))',ERR,ERROR,*999)
                                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                                & TRANSPOSE_COLUMN_INDICES,'("  Transpose column indices :",5(X,I13))', &
                                                & '(28X,5(X,I13))',ERR,ERROR,*999)
                                            ENDIF
                                          ENDIF
                                        ELSE
                                          CALL FlagError("Column dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FlagError("Row dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Column dofs domain mapping is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Row dofs domain mapping is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Column field variable is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Row field variable is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Interface condition interface is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE DEFAULT
                  LOCAL_ERROR="The matrix storage type of "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE DEFAULT
                LOCAL_ERROR="The matrix structure type of "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRIX%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*998)
              END SELECT
            ELSE
              CALL FlagError("Transpose column indices is already associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FlagError("Transpose row indieces is already associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FlagError("Column indices is already associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FlagError("Row indieces is already associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*998)
    ENDIF
     
    EXITS("INTERFACE_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(TRANSPOSE_ROW_INDICES)) DEALLOCATE(TRANSPOSE_ROW_INDICES)
    IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) DEALLOCATE(TRANSPOSE_COLUMN_INDICES)
    IF(ALLOCATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(TRANSPOSE_COLUMNS)) DEALLOCATE(TRANSPOSE_COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO local_row=1,SIZE(COLUMN_INDICES_LISTS,1)
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(local_row)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(local_row)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_row
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
    IF(ALLOCATED(TRANSPOSE_COLUMN_INDICES_LISTS)) THEN
      DO local_column=1,SIZE(TRANSPOSE_COLUMN_INDICES_LISTS,1)
        IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR)) &
          & CALL LIST_DESTROY(TRANSPOSE_COLUMN_INDICES_LISTS(local_column)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_row
      DEALLOCATE(TRANSPOSE_COLUMN_INDICES_LISTS)
    ENDIF
998 ERRORSEXITS("INTERFACE_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRIX_STRUCTURE_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the interface matrices for the interface equations
  SUBROUTINE INTERFACE_MATRICES_CREATE_FINISH(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<The pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:),TRANSPOSE_ROW_INDICES(:),TRANSPOSE_COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)
    NULLIFY(TRANSPOSE_ROW_INDICES)
    NULLIFY(TRANSPOSE_COLUMN_INDICES)

    NULLIFY(ROW_DOMAIN_MAP)
    NULLIFY(COLUMN_DOMAIN_MAP)

    ENTERS("INTERFACE_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FlagError("Interface matrices have already been finished.",ERR,ERROR,*998)
      ELSE
        INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          COLUMN_DOMAIN_MAP=>INTERFACE_MAPPING%COLUMN_DOFS_MAPPING
          IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
            !Now create the individual interface matrices
            DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
              INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
              IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
                ROW_DOMAIN_MAP=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%ROW_DOFS_MAPPING
                IF(ASSOCIATED(ROW_DOMAIN_MAP)) THEN
                  !Create the distributed equations matrix
                  CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,INTERFACE_MATRICES% &
                    & MATRICES(matrix_idx)%PTR%MATRIX,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(INTERFACE_MATRIX%MATRIX,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(INTERFACE_MATRIX%MATRIX,INTERFACE_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                  IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                    CALL DISTRIBUTED_MATRIX_CREATE_START(COLUMN_DOMAIN_MAP,ROW_DOMAIN_MAP,INTERFACE_MATRICES% &
                      & MATRICES(matrix_idx)%PTR%MATRIX_TRANSPOSE,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(INTERFACE_MATRIX%MATRIX_TRANSPOSE,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(INTERFACE_MATRIX%MATRIX_TRANSPOSE,INTERFACE_MATRIX%STORAGE_TYPE, &
                      & ERR,ERROR,*999)
                  ENDIF
                  !Calculate and set the matrix structure/sparsity pattern
                  IF(INTERFACE_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                    & INTERFACE_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                    CALL INTERFACE_MATRIX_STRUCTURE_CALCULATE(INTERFACE_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                      & TRANSPOSE_ROW_INDICES,TRANSPOSE_COLUMN_INDICES,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(INTERFACE_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(INTERFACE_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                      & ERR,ERROR,*999)
                    IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                      CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(INTERFACE_MATRIX%MATRIX_TRANSPOSE,NUMBER_OF_NON_ZEROS, &
                        & ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(INTERFACE_MATRIX%MATRIX_TRANSPOSE,TRANSPOSE_ROW_INDICES, &
                        & TRANSPOSE_COLUMN_INDICES,ERR,ERROR,*999)
                     ENDIF
                    IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                    IF(ASSOCIATED(TRANSPOSE_ROW_INDICES)) DEALLOCATE(TRANSPOSE_ROW_INDICES)
                    IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) DEALLOCATE(TRANSPOSE_COLUMN_INDICES)
                  ENDIF
                  CALL DISTRIBUTED_MATRIX_CREATE_FINISH(INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
                  IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                    CALL DISTRIBUTED_MATRIX_CREATE_FINISH(INTERFACE_MATRIX%MATRIX_TRANSPOSE,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Row domain map for interface matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Interface matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
            RHS_VECTOR=>INTERFACE_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              !Set up the interface RHS vector
              CALL DISTRIBUTED_VECTOR_CREATE_START(COLUMN_DOMAIN_MAP,INTERFACE_MATRICES%RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(RHS_VECTOR%RHS_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
            ENDIF
            !Finish up
            INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED=.TRUE.
          ELSE
            CALL FlagError("Column domain map is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("INTERFACE_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(TRANSPOSE_ROW_INDICES)) DEALLOCATE(TRANSPOSE_ROW_INDICES)
    IF(ASSOCIATED(TRANSPOSE_COLUMN_INDICES)) DEALLOCATE(TRANSPOSE_COLUMN_INDICES)
    CALL INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_MATRICES_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of the interface matrices and rhs for the interface equations
  SUBROUTINE INTERFACE_MATRICES_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<The pointer to the interface equations to create the interface equations matrices for
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<On return, a pointer to the interface matrices being created. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables

    ENTERS("INTERFACE_MATRICES_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN      
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
          CALL FlagError("Interface matrices is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(INTERFACE_MATRICES)
          !Initialise the interface matrices
          CALL INTERFACE_MATRICES_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
        ENDIF
      ELSE
        CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MATRICES_CREATE_START")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_CREATE_START",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRICES_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the interface matrices
  SUBROUTINE INTERFACE_MATRICES_DESTROY(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer the interface matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      CALL INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
        
    EXITS("INTERFACE_MATRICES_DESTROY")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_DESTROY",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE INTERFACE_MATRICES_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the interface matrices and deallocate all memory.
  SUBROUTINE INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
   
    ENTERS("INTERFACE_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(ALLOCATED(INTERFACE_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(INTERFACE_MATRICES%MATRICES,1)
          CALL INTERFACE_MATRIX_FINALISE(INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(INTERFACE_MATRICES%MATRICES)
      ENDIF
      CALL INTERFACE_MATRICES_RHS_FINALISE(INTERFACE_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_MATRICES)
    ENDIF
       
    EXITS("INTERFACE_MATRICES_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the interface matrices for the interface equations.
  SUBROUTINE INTERFACE_MATRICES_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)
    
     !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to initialise the interface matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    ENTERS("INTERFACE_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MATRICES)) THEN
        CALL FlagError("Interface matrices is already associated for this interface equations.",ERR,ERROR,*998)
      ELSE
        INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
            ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface equations interface matrices.",ERR,ERROR,*999)
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED=.FALSE.
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_MAPPING=>INTERFACE_MAPPING
            NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MATRICES%SOLVER_MAPPING)
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_COLUMNS=INTERFACE_MAPPING%NUMBER_OF_COLUMNS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%TOTAL_NUMBER_OF_COLUMNS=INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_GLOBAL_COLUMNS=INTERFACE_MAPPING%NUMBER_OF_GLOBAL_COLUMNS
            NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MATRICES%RHS_VECTOR)
            !Allocate and initialise the matrices
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
            ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(INTERFACE_EQUATIONS%INTERFACE_MATRICES% &
              & NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
              NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL INTERFACE_MATRIX_INITIALISE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
            CALL INTERFACE_MATRICES_RHS_INITIALISE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,ERR,ERROR,*999)
          ELSE
            CALL FlagError("Interface mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface equations interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("INTERFACE_MATRICES_INITIALISE")
    RETURN
999 CALL INTERFACE_MATRICES_FINALISE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_MATRICES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Outputs the interface matrices
  SUBROUTINE INTERFACE_MATRICES_OUTPUT(ID,INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR
    
    ENTERS("INTERFACE_MATRICES_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"Interface matrices:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Number of interface matrices = ",INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES, &
          & ERR,ERROR,*999)
        DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
            CALL WRITE_STRING_VALUE(ID,"Interface matrix : ",matrix_idx,ERR,ERROR,*999)
            CALL WRITE_STRING(ID,"Standard matrix:",ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_OUTPUT(ID,INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
            IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
              CALL WRITE_STRING(ID,"Transposed matrix:",ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,INTERFACE_MATRIX%MATRIX_TRANSPOSE,ERR,ERROR,*999)
            ENDIF
         ELSE
            CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
        RHS_VECTOR=>INTERFACE_MATRICES%RHS_VECTOR
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          CALL WRITE_STRING(ID,"Interface RHS vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MATRICES_OUTPUT")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_OUTPUT",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRICES_OUTPUT
  
  !
  !================================================================================================================================
  !

  !>Finalises the interface matrices RHS vector and deallocates all memory
  SUBROUTINE INTERFACE_MATRICES_RHS_FINALISE(RHS_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR !<A pointer to the equation matrices RHS vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    ENTERS("INTERFACE_MATRICES_RHS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(RHS_VECTOR)) THEN
      IF(ASSOCIATED(RHS_VECTOR%RHS_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(RHS_VECTOR%RHS_VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(RHS_VECTOR)
    ENDIF      
     
    EXITS("INTERFACE_MATRICES_RHS_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_RHS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_RHS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the interface matrices RHS vector
  SUBROUTINE INTERFACE_MATRICES_RHS_INITIALISE(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the equation matrices to initialise the rhs vector for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    ENTERS("INTERFACE_MATRICES_RHS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
      IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
        RHS_MAPPING=>INTERFACE_MAPPING%RHS_MAPPING
        IF(ASSOCIATED(RHS_MAPPING)) THEN
          IF(ASSOCIATED(INTERFACE_MATRICES%RHS_VECTOR)) THEN
            CALL FlagError("Interface matrices RHS vector is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(INTERFACE_MATRICES%RHS_VECTOR,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface matrices RHS vector.",ERR,ERROR,*999)
            INTERFACE_MATRICES%RHS_VECTOR%UPDATE_VECTOR=.TRUE.
            INTERFACE_MATRICES%RHS_VECTOR%FIRST_ASSEMBLY=.TRUE.
            NULLIFY(INTERFACE_MATRICES%RHS_VECTOR%RHS_VECTOR)
            CALL EquationsMatrices_ElementVectorInitialise(INTERFACE_MATRICES%RHS_VECTOR%ELEMENT_VECTOR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Interface matrices equation mapping is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("INTERFACE_MATRICES_RHS_INITIALISE")
    RETURN
999 CALL INTERFACE_MATRICES_RHS_FINALISE(INTERFACE_MATRICES%RHS_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_MATRICES_RHS_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_RHS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the interface matrices
  SUBROUTINE INTERFACE_MATRICES_STORAGE_TYPE_SET(INTERFACE_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th inteface matrices. \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesSparsityTypes,INTERFACE_MATRICES_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("INTERFACE_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FlagError("Interface matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STORAGE_TYPE,1)==INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
          DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
            INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
              SELECT CASE(STORAGE_TYPE(matrix_idx))
              CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE        
              CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
              CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                INTERFACE_MATRIX%STORAGE_TYPE=DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE
              CASE DEFAULT
                LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                  & " for interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of interface matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRICES_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the interface matrices.
  SUBROUTINE INTERFACE_MATRICES_STRUCTURE_TYPE_SET(INTERFACE_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The structure type for the  matrix_idx'th interface matrix \see INTERFACE_MATRICES_ROUTINES_InterfaceMatrixStructureTypes,INTERFACE_MATRICES_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FlagError("Interface matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STRUCTURE_TYPE,1)==INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
          DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
            INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
              SELECT CASE(STRUCTURE_TYPE(matrix_idx))
              CASE(INTERFACE_MATRIX_NO_STRUCTURE)
                INTERFACE_MATRIX%STRUCTURE_TYPE=INTERFACE_MATRIX_NO_STRUCTURE
              CASE(INTERFACE_MATRIX_FEM_STRUCTURE)
                INTERFACE_MATRIX%STRUCTURE_TYPE=INTERFACE_MATRIX_FEM_STRUCTURE
              CASE DEFAULT
                LOCAL_ERROR="The specified strucutre type of "// &
                  & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for interface matrix number "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of interface matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_STRUCTURE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Initialise the values of the interface matrices to the given value e.g., 0.0_DP
  SUBROUTINE INTERFACE_MATRICES_VALUES_INITIALISE(INTERFACE_MATRICES,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices to initialise the values for
    REAL(DP), INTENT(IN) :: VALUE !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: RHS_VECTOR
    
    ENTERS("INTERFACE_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
        INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
          IF(INTERFACE_MATRIX%UPDATE_MATRIX) THEN
            CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(INTERFACE_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
            IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
              CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(INTERFACE_MATRIX%MATRIX_TRANSPOSE,VALUE,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("Interface matrix is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      RHS_VECTOR=>INTERFACE_MATRICES%RHS_VECTOR
      IF(ASSOCIATED(RHS_VECTOR)) THEN
        IF(RHS_VECTOR%UPDATE_VECTOR) THEN
          CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(RHS_VECTOR%RHS_VECTOR,VALUE,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("INTERFACE_MATRICES_VALUES_INITIALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_MATRICES_VALUES_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_VALUES_INITIALISE

  !
  !================================================================================================================================
  !
  
  SUBROUTINE InterfaceMatrix_TimeDependenceTypeSet(InterfaceCondition, &
    & interfaceMatrixIndex,IsTranspose,TimeDependenceType,Err,Error,*)
    
    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: InterfaceCondition
    INTEGER(INTG), INTENT(IN) :: InterfaceMatrixIndex
    LOGICAL, INTENT(IN) :: IsTranspose
    INTEGER(INTG), INTENT(IN) :: TimeDependenceType
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: InterfaceEquations
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: InterfaceMatrices
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: InterfaceMatrix
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("InterfaceMatrix_TimeDependenceTypeSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(InterfaceCondition)) THEN
      InterfaceEquations=>InterfaceCondition%INTERFACE_EQUATIONS
      IF(ASSOCIATED(InterfaceEquations)) THEN
        InterfaceMatrices=>InterfaceEquations%INTERFACE_MATRICES
        IF(ASSOCIATED(InterfaceMatrices)) THEN
          InterfaceMatrix=>InterfaceMatrices%MATRICES(InterfaceMatrixIndex)%PTR
          IF(ASSOCIATED(InterfaceMatrix)) THEN
            IF(.NOT.IsTranspose) THEN
              InterfaceMatrix%INTERFACE_MATRIX_TIME_DEPENDENCE_TYPE=TimeDependenceType
            ELSE
              IF(InterfaceMatrix%HAS_TRANSPOSE) THEN
                InterfaceMatrix%INTERFACE_MATRIX_TRANSPOSE_TIME_DEPENDENCE_TYPE=TimeDependenceType
              ELSE
                LOCAL_ERROR="Interface matrices has_transpose flag is .false. but interface matrix type is transpose."
                CALL FlagError(LOCAL_ERROR,Err,Error,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FlagError("Interface matrix is not associated",Err,Error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface matrices not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface equations not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface condition is not associated.",Err,Error,*999)
    ENDIF
    
    EXITS("InterfaceMatrix_TimeDependenceTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_TimeDependenceTypeSet",Err,Error)
    RETURN 1
  END SUBROUTINE InterfaceMatrix_TimeDependenceTypeSet

  !
  !================================================================================================================================
  !
  
  SUBROUTINE InterfaceMatrix_TimeDependenceTypeGet(InterfaceCondition, &
    & interfaceMatrixIndex,IsTranspose,TimeDependenceType,Err,Error,*)
    
    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: InterfaceCondition
    INTEGER(INTG), INTENT(IN) :: InterfaceMatrixIndex
    LOGICAL, INTENT(IN) :: IsTranspose
    INTEGER(INTG), INTENT(OUT) :: TimeDependenceType
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Local variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: InterfaceEquations
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: InterfaceMatrices
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: InterfaceMatrix
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("InterfaceMatrix_TimeDependenceTypeGet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(InterfaceCondition)) THEN
      InterfaceEquations=>InterfaceCondition%INTERFACE_EQUATIONS
      IF(ASSOCIATED(InterfaceEquations)) THEN
        InterfaceMatrices=>InterfaceEquations%INTERFACE_MATRICES
        IF(ASSOCIATED(InterfaceMatrices)) THEN
          InterfaceMatrix=>InterfaceMatrices%MATRICES(InterfaceMatrixIndex)%PTR
          IF(ASSOCIATED(InterfaceMatrix)) THEN
            IF(.NOT.IsTranspose) THEN
              TimeDependenceType=InterfaceMatrix%INTERFACE_MATRIX_TIME_DEPENDENCE_TYPE
            ELSE
              IF(InterfaceMatrix%HAS_TRANSPOSE) THEN
                TimeDependenceType=InterfaceMatrix%INTERFACE_MATRIX_TRANSPOSE_TIME_DEPENDENCE_TYPE
              ELSE
                LOCAL_ERROR="Interface matrices has_transpose flag is .false. but interface matrix type is transpose."
                CALL FlagError(LOCAL_ERROR,Err,Error,*999)
              ENDIF
            ENDIF
            !Sanity check
            IF(.NOT.(TimeDependenceType>0.AND.TimeDependenceType<=NUMBER_OF_INTERFACE_MATRIX_TYPES)) THEN
              LOCAL_ERROR="Invalid time dependence type of "//TRIM(NUMBER_TO_VSTRING(TimeDependenceType,"*",ERR,ERROR))// &
                & ". Must be > 0 and <= "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_INTERFACE_MATRIX_TYPES,"*",ERR,ERROR))// &
                & "."
              CALL FlagError(LOCAL_ERROR,Err,Error,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface matrix is not associated",Err,Error,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface matrices not associated.",Err,Error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface equations not associated.",Err,Error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface condition is not associated.",Err,Error,*999)
    ENDIF
    
    EXITS("InterfaceMatrix_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceMatrix_TimeDependenceTypeGet",Err,Error)
    RETURN 1
  END SUBROUTINE InterfaceMatrix_TimeDependenceTypeGet

  !
  !================================================================================================================================
  !


END MODULE INTERFACE_MATRICES_ROUTINES
