!> \file
!> $Id: equations_matrices_routines.f90 28 2007-07-27 08:35:14Z cpb $
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

!> This module handles all equations matrix and rhs routines.
MODULE EQUATIONS_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes EQUATIONS_MATRICES_ROUTINES::EquationsMatrixStructureTypes
  !> \brief Equations matrices structure (sparsity) types
  !> \see EQUATIONS_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see EQUATIONS_MATRICES_ROUTINES_EquationsMatrixStructureTypes,EQUATIONS_MATRICES_ROUTINES
  !>@}

  !> \addtogroup EQUATIONS_MATRICES_ROUTINES_EquationsMatricesSparsityTypes EQUATIONS_MATRICES_ROUTINES::EquationsMatricesSparsityTypes
  !> \brief Equations matrices sparsity types
  !> \see EQUATIONS_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_SPARSE_MATRICES=1 !<Use sparse equations matrices \see EQUATIONS_MATRICES_ROUTINES_EquationsMatricesSparsityTypes,EQUATIONS_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRICES_FULL_MATRICES=2 !<Use fully populated equation matrices \see EQUATIONS_MATRICES_ROUTINES_EquationsMatricesSparsityTypes,EQUATIONS_MATRICES_ROUTINES
 !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATIONS_MATRIX_NO_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE

  PUBLIC EQUATIONS_MATRICES_SPARSE_MATRICES,EQUATIONS_MATRICES_FULL_MATRICES

  PUBLIC EQUATIONS_MATRICES_CREATE_FINISH,EQUATIONS_MATRICES_CREATE_START,EQUATIONS_MATRICES_DESTROY

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC EQUATIONS_MATRICES_ELEMENT_ADD,EQUATIONS_MATRICES_ELEMENT_CALCULATE,EQUATIONS_MATRICES_ELEMENT_INITIALISE, &
    & EQUATIONS_MATRICES_ELEMENT_FINALISE,EQUATIONS_MATRICES_VALUES_INITIALISE

  PUBLIC EQUATIONS_MATRICES_OUTPUT

  PUBLIC EQUATIONS_MATRICES_STORAGE_TYPE_SET,EQUATIONS_MATRICES_STRUCTURE_TYPE_SET

CONTAINS

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
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: VARIABLE_DOMAIN_MAP
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: VARIABLE
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)

    CALL ENTERS("EQUATIONS_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have already been finished.",ERR,ERROR,*998)
      ELSE
        EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          !Now create the individual equations matrices
          DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
            VARIABLE=>EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE
            VARIABLE_DOMAIN_MAP=>VARIABLE%DOMAIN_MAPPING
            !Create the distributed equations matrix
            CALL DISTRIBUTED_MATRIX_CREATE_START(VARIABLE_DOMAIN_MAP,EQUATIONS_MATRICES%MATRICES(matrix_idx)% &
              & MATRIX,ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX, &
              & MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX, &
              & EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE,ERR,ERROR,*999)
            !Calculate and set the matrix structure/sparsity pattern
            IF(EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE/=MATRIX_BLOCK_STORAGE_TYPE) THEN
              CALL EQUATIONS_MATRIX_STRUCTURE_CALCULATE(EQUATIONS_MATRICES%MATRICES(matrix_idx),NUMBER_OF_NON_ZEROS, &
                & ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                & NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX,ROW_INDICES, &
                & COLUMN_INDICES,ERR,ERROR,*999)
              IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
              IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
            ENDIF
            CALL DISTRIBUTED_MATRIX_CREATE_FINISH(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
          ENDDO !matrix_idx
          !Finish setting up the equations RHS vector
          VARIABLE=>EQUATIONS_MAPPING%RHS_VARIABLE
          IF(ASSOCIATED(VARIABLE)) THEN
            VARIABLE_DOMAIN_MAP=>VARIABLE%DOMAIN_MAPPING
            CALL DISTRIBUTED_VECTOR_CREATE_START(VARIABLE_DOMAIN_MAP,EQUATIONS_MATRICES%VECTOR,ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(EQUATIONS_MATRICES%VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_CREATE_FINISH(EQUATIONS_MATRICES%VECTOR,ERR,ERROR,*999)
          ENDIF
          EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED=.TRUE.
        ELSE
          CALL FLAG_ERROR("Equations mapping is not associated",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*998)
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
  SUBROUTINE EQUATIONS_MATRICES_CREATE_START(EQUATIONS_MAPPING,EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING !<The pointer to the equations mapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<On return, a pointer to the equations matrices being created.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: DUMMY_ERROR    

    CALL ENTERS("EQUATIONS_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
      IF(EQUATIONS_MAPPING%EQUATIONS_MAPPING_FINISHED) THEN
        EQUATIONS=>EQUATIONS_MAPPING%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            CALL FLAG_ERROR("Equations matrices is already associated",ERR,ERROR,*998)
          ELSE
            NULLIFY(EQUATIONS_MATRICES)
            !Initialise the equations matrices
            CALL EQUATIONS_MATRICES_INITIALISE(EQUATIONS,ERR,ERROR,*999)
            EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations mapping has not been finished",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*998)
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

  !>Adds the element matrices and rhs vector into the equations matrices and rhs vector.
  SUBROUTINE EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      !Add the element matrices
      DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
        IF(EQUATIONS_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX) THEN
          CALL DISTRIBUTED_MATRIX_VALUES_ADD(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX,EQUATIONS_MATRICES% &
            & MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS(1:EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX% &
            & NUMBER_OF_ROWS),EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS(1:EQUATIONS_MATRICES%MATRICES( &
            & matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS),EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX(1: &
            & EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS,1:EQUATIONS_MATRICES%MATRICES(matrix_idx)% &
            & ELEMENT_MATRIX%NUMBER_OF_COLUMNS),ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      !Add the rhs
      IF(EQUATIONS_MATRICES%UPDATE_VECTOR) THEN
        CALL DISTRIBUTED_VECTOR_VALUES_ADD(EQUATIONS_MATRICES%VECTOR,EQUATIONS_MATRICES%ELEMENT_VECTOR%ROW_DOFS(1: &
          & EQUATIONS_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS),EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR(1:EQUATIONS_MATRICES% &
          & ELEMENT_VECTOR%NUMBER_OF_ROWS),ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
    
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
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,global_ny,local_ny,matrix_idx,node,node_idx
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        !Calculate the row and columns for the equations matrices
        DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
          EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS=0
          EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
          IF(EQUATIONS_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX) THEN
            FIELD_VARIABLE=>EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE
            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
                IF(ELEMENT_NUMBER>=1.AND.ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                  BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
!!TODO :: CHANGE TO LOOP OVER ELEMENT PARAMETERS (ns)
                  DO node_idx=1,BASIS%NUMBER_OF_NODES
                    node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                    DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                      derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node,1)
                      global_ny=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                      EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS= &
                        & EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                      EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS= &
                        & EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                      EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS(EQUATIONS_MATRICES%MATRICES(matrix_idx)% &
                        & ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                      EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS(EQUATIONS_MATRICES%MATRICES(matrix_idx)% &
                        & ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                    ENDDO !derivative_idx
                  ENDDO !node_idx
                ELSE
                  LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                    & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                    & ". The element number must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The interpolation type for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                  & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                  & " is not nodally based."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
              ENDIF
            ENDDO !component_idx
            EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX=0.0_DP
          ENDIF
        ENDDO !matrix_idx
        !Calculate the row and columns for the equations RHS
        EQUATIONS_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS=0
        IF(EQUATIONS_MATRICES%UPDATE_VECTOR) THEN
          FIELD_VARIABLE=>EQUATIONS_MAPPING%RHS_VARIABLE
          DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
              IF(ELEMENT_NUMBER>=1.AND.ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
                BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
                DO node_idx=1,BASIS%NUMBER_OF_NODES
                  node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                  DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                    derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node,1)
                    EQUATIONS_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS=EQUATIONS_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                    EQUATIONS_MATRICES%ELEMENT_VECTOR%ROW_DOFS(EQUATIONS_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
                  ENDDO !derivative_idx
                ENDDO !node_idx
              ELSE
                LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                  & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                  & ". The element number must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The interpolation type for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                & " is not nodally based."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
            ENDIF
          ENDDO !component_idx
          EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR=0.0_DP
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
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
    
    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      !Finalise the element matrices
      DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
        EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
        EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
        IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS)) &
          & DEALLOCATE(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS)
        IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS)) &
          & DEALLOCATE(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS)
        IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX)) &
          & DEALLOCATE(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX)        
      ENDDO !matrix_idx
      !Finalise the element vector
      EQUATIONS_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
      IF(ALLOCATED(EQUATIONS_MATRICES%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(EQUATIONS_MATRICES%ELEMENT_VECTOR%ROW_DOFS)
      IF(ALLOCATED(EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR)
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
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_ELEMENT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
        !Finalise the element matrices
        DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
          FIELD_VARIABLE=>EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(matrix_idx)%VARIABLE
          EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS= &
            & FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
          EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS= &
            & FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
          IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS)) THEN
            CALL FLAG_ERROR("Element matrix row dofs already allocated.",ERR,ERROR,*999)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS(EQUATIONS_MATRICES%MATRICES(matrix_idx)% &
              & ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix row dofs.",ERR,ERROR,*999)
          ENDIF
          IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS)) THEN
            CALL FLAG_ERROR("Element matrix column dofs already allocated.",ERR,ERROR,*999)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS(EQUATIONS_MATRICES%MATRICES(matrix_idx)% &
              & ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix column dofs.",ERR,ERROR,*999)
          ENDIF
          IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX)) THEN
            CALL FLAG_ERROR("Element matrix already allocated.",ERR,ERROR,*999)
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX( &
              & EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
              & EQUATIONS_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix.",ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
        !Initialise the element vector
        FIELD_VARIABLE=>EQUATIONS_MAPPING%RHS_VARIABLE
        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          EQUATIONS_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
          IF(ALLOCATED(EQUATIONS_MATRICES%ELEMENT_VECTOR%ROW_DOFS)) THEN
            CALL FLAG_ERROR("Element vector row dofs already allocated.",ERR,ERROR,*999)        
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%ELEMENT_VECTOR%ROW_DOFS(EQUATIONS_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector row dofs.",ERR,ERROR,*999)
          ENDIF
          IF(ALLOCATED(EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR)) THEN
            CALL FLAG_ERROR("Element vector already allocated.",ERR,ERROR,*999)        
          ELSE
            ALLOCATE(EQUATIONS_MATRICES%ELEMENT_VECTOR%VECTOR(EQUATIONS_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector.",ERR,ERROR,*999)
          ENDIF
        ELSE
          EQUATIONS_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not allocated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_INITIALISE")
    RETURN
999 CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_MATRICES_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise a equations matrix and deallocate all memory
  SUBROUTINE EQUATIONS_MATRIX_FINALISE(EQUATIONS_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TYPE):: EQUATIONS_MATRIX !<The equations matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("EQUATIONS_MATRIX_FINALISE",ERR,ERROR,*999)

    CALL DISTRIBUTED_MATRIX_DESTROY(EQUATIONS_MATRIX%MATRIX,ERR,ERROR,*999)
    CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(EQUATIONS_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
    
    CALL EXITS("EQUATIONS_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the equations matrix.
  SUBROUTINE EQUATIONS_MATRIX_INITIALISE(EQUATIONS_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices to initialise the equations matrix for
    INTEGER(INTG) :: MATRIX_NUMBER !<The matrix number in the equations matrices to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRIX_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=EQUATIONS_MATRICES%NUMBER_OF_MATRICES) THEN
        EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%MATRIX_NUMBER=MATRIX_NUMBER
          EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%EQUATIONS_MATRICES=>EQUATIONS_MATRICES
          EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
          EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
          EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%UPDATE_MATRIX=.TRUE.
          EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%NUMBER_OF_COLUMNS=EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS( &
            & MATRIX_NUMBER)%NUMBER_OF_COLUMNS
          NULLIFY(EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%MATRIX)
          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(EQUATIONS_MATRICES%MATRICES(MATRIX_NUMBER)%ELEMENT_MATRIX, &
            & ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRIX_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRIX_INITIALISE

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
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    
    CALL ENTERS("EQUATIONS_MATRICES_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          CALL WRITE_STRING(ID,"Equations matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of matrices = ",EQUATIONS_MATRICES%NUMBER_OF_MATRICES,ERR,ERROR,*999)
          DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
            CALL WRITE_STRING_VALUE(ID,"Equations matrix : ",matrix_idx,ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_OUTPUT(ID,EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
          ENDDO !matrix_idx
          IF(EQUATIONS_MAPPING%RHS_VARIABLE_TYPE/=0) THEN
            CALL WRITE_STRING(ID,"Equations RHS vector:",ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,EQUATIONS_MATRICES%VECTOR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
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

  !>Sets the storage type (sparsity) of the equations matrices
  SUBROUTINE EQUATIONS_MATRICES_STORAGE_TYPE_SET(EQUATIONS_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the eqautions matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STORAGE_TYPE,1)==EQUATIONS_MATRICES%NUMBER_OF_MATRICES) THEN
          DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
            SELECT CASE(STORAGE_TYPE(matrix_idx))
            CASE(MATRIX_BLOCK_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
            CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_DIAGONAL_STORAGE_TYPE        
            CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
            CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_ROW_MAJOR_STORAGE_TYPE
            CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
            CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
            CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=MATRIX_ROW_COLUMN_STORAGE_TYPE
            CASE DEFAULT
              LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                & " for the matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the equations matrices
  SUBROUTINE EQUATIONS_MATRICES_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The storage type for the matrix_idx'th equations matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Equations matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STRUCTURE_TYPE,1)==EQUATIONS_MATRICES%NUMBER_OF_MATRICES) THEN
          DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
            SELECT CASE(STRUCTURE_TYPE(matrix_idx))
            CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE=EQUATIONS_MATRIX_NO_STRUCTURE
            CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
              EQUATIONS_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE=EQUATIONS_MATRIX_FEM_STRUCTURE
            CASE DEFAULT
              LOCAL_ERROR="The specified strucutre type of "// &
                & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for the matrix number "// &
                & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_MATRICES_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_MATRICES_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_MATRICES_STRUCTURE_TYPE_SET
  
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
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("EQUATIONS_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      IF(ALLOCATED(EQUATIONS_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(EQUATIONS_MATRICES%MATRICES,1)
          CALL EQUATIONS_MATRIX_FINALISE(EQUATIONS_MATRICES%MATRICES(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(EQUATIONS_MATRICES%MATRICES)
      ENDIF
      CALL DISTRIBUTED_VECTOR_DESTROY(EQUATIONS_MATRICES%VECTOR,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_FINALISE(EQUATIONS_MATRICES%ELEMENT_VECTOR,ERR,ERROR,*999)
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
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MATRICES)) THEN
        CALL FLAG_ERROR("Equations matrices is already associated for this equations.",ERR,ERROR,*998)
      ELSE
        EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
        IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
          ALLOCATE(EQUATIONS%EQUATIONS_MATRICES,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations equations matrices.",ERR,ERROR,*999)
          EQUATIONS%EQUATIONS_MATRICES%EQUATIONS=>EQUATIONS
          EQUATIONS%EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED=.FALSE.
          EQUATIONS%EQUATIONS_MATRICES%EQUATIONS_MAPPING=>EQUATIONS_MAPPING
          NULLIFY(EQUATIONS%EQUATIONS_MATRICES%SOLUTION_MAPPING)
          EQUATIONS%EQUATIONS_MATRICES%NUMBER_OF_ROWS=EQUATIONS_MAPPING%NUMBER_OF_ROWS
          EQUATIONS%EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS=EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
          EQUATIONS%EQUATIONS_MATRICES%NUMBER_OF_MATRICES=EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
          ALLOCATE(EQUATIONS%EQUATIONS_MATRICES%MATRICES(EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations matrices matrices.",ERR,ERROR,*999)
          DO matrix_idx=1,EQUATIONS_MAPPING%NUMBER_OF_EQUATIONS_MATRICES
            CALL EQUATIONS_MATRIX_INITIALISE(EQUATIONS%EQUATIONS_MATRICES,matrix_idx,ERR,ERROR,*999)
          ENDDO !matrix_idx
          EQUATIONS%EQUATIONS_MATRICES%UPDATE_VECTOR=.TRUE.
          NULLIFY(EQUATIONS%EQUATIONS_MATRICES%VECTOR)
          CALL EQUATIONS_MATRICES_ELEMENT_VECTOR_INITIALISE(EQUATIONS%EQUATIONS_MATRICES%ELEMENT_VECTOR,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*998)
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

  !>Initialise the values of the equations matrices and rhs to the given value e.g., 0.0_DP
  SUBROUTINE EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES !<A pointer to the equations matrices to initialise the values for
    REAL(DP), INTENT(IN) :: VALUE !<The value to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("EQUATIONS_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
      DO matrix_idx=1,EQUATIONS_MATRICES%NUMBER_OF_MATRICES
        IF(EQUATIONS_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX) THEN
          CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(EQUATIONS_MATRICES%MATRICES(matrix_idx)%MATRIX,VALUE,ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      IF(EQUATIONS_MATRICES%UPDATE_VECTOR) THEN
        CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(EQUATIONS_MATRICES%VECTOR,VALUE,ERR,ERROR,*999)
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
  SUBROUTINE EQUATIONS_MATRIX_STRUCTURE_CALCULATE(EQUATIONS_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_MATRIX_TYPE) :: EQUATIONS_MATRIX !<The equations matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::  column_idx,DUMMY_ERR,elem_idx,global_column,local_column,local_ny,MATRIX_NUMBER,mk,mp,ne,nh,nn,nnk,np, &
      & NUMBER_OF_COLUMNS,nyy
    INTEGER(INTG), POINTER :: COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: DEPENDENT_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(COLUMNS)
    
    CALL ENTERS("EQUATIONS_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*998)

    NUMBER_OF_NON_ZEROS=0
    IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
      IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
        MATRIX_NUMBER=EQUATIONS_MATRIX%MATRIX_NUMBER
        SELECT CASE(EQUATIONS_MATRIX%STRUCTURE_TYPE)
        CASE(EQUATIONS_MATRIX_NO_STRUCTURE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*998)
        CASE(EQUATIONS_MATRIX_FEM_STRUCTURE)
          SELECT CASE(EQUATIONS_MATRIX%STORAGE_TYPE)
          CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
            EQUATIONS_MATRICES=>EQUATIONS_MATRIX%EQUATIONS_MATRICES
            IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
              EQUATIONS=>EQUATIONS_MATRICES%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_MAPPING=>EQUATIONS_MATRICES%EQUATIONS_MAPPING
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                  EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                      FIELD_VARIABLE=>EQUATIONS_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(MATRIX_NUMBER)%VARIABLE
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        DEPENDENT_DOFS_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                        IF(ASSOCIATED(DEPENDENT_DOFS_DOMAIN_MAPPING)) THEN
                          DEPENDENT_DOFS_PARAM_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOF_TO_PARAM_MAP
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
                                  BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                                  DO nn=1,BASIS%NUMBER_OF_NODES
                                    mp=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                    DO nnk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                                      mk=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(nnk,nn)
                                      !Find the local and global column and add the global column to the indices list
                                      local_column=FIELD_VARIABLE%COMPONENTS(nh)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(mk,mp,1)
                                      global_column=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                      CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_ny)%PTR,global_column,ERR,ERROR,*999)
                                    ENDDO !mk
                                  ENDDO !nn
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
                              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix structure:",ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Equations matrix number : ",MATRIX_NUMBER, &
                                & ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ", &
                                & DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ", &
                                & DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",NUMBER_OF_NON_ZEROS, &
                                & ERR,ERROR,*999)
                              IF(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                                & DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL/=0) THEN
                                SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL* &
                                  DEPENDENT_DOFS_DOMAIN_MAPPING%NUMBER_OF_GLOBAL,DP)*100.0_DP
                                CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F5.2", &
                                  & ERR,ERROR,*999)
                              ENDIF
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DEPENDENT_DOFS_DOMAIN_MAPPING% &
                                & TOTAL_NUMBER_OF_LOCAL+1,8,8,ROW_INDICES,'("  Row indices    :",8(X,I13))','(18X,8(X,I13))', &
                                & ERR,ERROR,*999)
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
                  CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The matrix storage type of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
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
    
    CALL EXITS("EQUATIONS_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(COLUMNS)) DEALLOCATE(COLUMNS)
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
  
END MODULE EQUATIONS_MATRICES_ROUTINES
