!> \file
!> $Id: interface_matrices_routines.f90 690 2009-09-30 23:27:16Z chrispbradley $
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

!>This module contains all interface matrices routines.
MODULE INTERFACE_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES

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

  PUBLIC INTERFACE_MATRICES_ELEMENT_CALCULATE

  PUBLIC INTERFACE_MATRICES_ELEMENT_FINALISE,INTERFACE_MATRICES_ELEMENT_INITIALISE

  PUBLIC INTERFACE_MATRICES_OUTPUT

  PUBLIC INTERFACE_MATRICES_STORAGE_TYPE_SET

  PUBLIC INTERFACE_MATRICES_STRUCTURE_TYPE_SET

  PUBLIC INTERFACE_MATRICES_VALUES_INITIALISE

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
    
    CALL ENTERS("INTERFACE_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
      IF(ASSOCIATED(INTERFACE_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
      CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_MATRIX)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRIX_FINALISE")
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

    CALL ENTERS("INTERFACE_MATRIX_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
        INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          IF(ASSOCIATED(INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
            LOCAL_ERROR="Interface matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
              & " is already associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          ELSE
            ALLOCATE(INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface matrix.",ERR,ERROR,*999)
            INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
            INTERFACE_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
            INTERFACE_MATRIX%INTERFACE_MATRICES=>INTERFACE_MATRICES
            INTERFACE_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
            INTERFACE_MATRIX%STRUCTURE_TYPE=INTERFACE_MATRIX_NO_STRUCTURE
            INTERFACE_MATRIX%UPDATE_MATRIX=.TRUE.
            INTERFACE_MATRIX%FIRST_ASSEMBLY=.TRUE.
            INTERFACE_MATRIX%NUMBER_OF_ROWS=INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%NUMBER_OF_ROWS
            INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MATRIX_NUMBER)%INTERFACE_MATRIX=>INTERFACE_MATRIX
            NULLIFY(INTERFACE_MATRIX%MATRIX)
            CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_INITIALISE(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified interface matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The matrix number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRIX_INITIALISE")
    RETURN
999 CALL INTERFACE_MATRIX_FINALISE(INTERFACE_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRIX_INITIALISE")
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
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_ADD()")
#endif

    CALL ENTERS("INTERFACE_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

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
          ENDIF
        ELSE
          LOCAL_ERROR="Interface matrix for interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
    ELSE
      CALL FLAG_ERROR("Interface matrices is not allocated.",ERR,ERROR,*999)
    ENDIF
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_ADD()")
#endif
    
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_ADD")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_ELEMENT_ADD",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_ADD")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the interface matrices of the element matrices.
  SUBROUTINE INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !<A pointer to the interface matrices
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx 
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLS_FIELD_VARIABLE,ROWS_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("INTERFACE_MATRICES_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
      IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
        !Calculate the row and columns for the interface equations matrices
        DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
            ROWS_FIELD_VARIABLE=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE
            COLS_FIELD_VARIABLE=>INTERFACE_MAPPING%LAGRANGE_VARIABLE !TEMPORARY: Needs generalising
            CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_CALCULATE(INTERFACE_MATRIX%ELEMENT_MATRIX, &
              & INTERFACE_MATRIX%UPDATE_MATRIX,ELEMENT_NUMBER,ROWS_FIELD_VARIABLE,COLS_FIELD_VARIABLE, &
              & ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ELSE
        CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_CALCULATE()")
#endif
    
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_ELEMENT_CALCULATE

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
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
        INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
          CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_FINALISE(INTERFACE_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
            & " is not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the interface matrices
  SUBROUTINE INTERFACE_MATRICES_ELEMENT_INITIALISE(INTERFACE_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES !The interface matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: COLS_FIELD_VARIABLE,ROWS_FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("INTERFACE_MATRICES_ELEMENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      INTERFACE_MAPPING=>INTERFACE_MATRICES%INTERFACE_MAPPING
      IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
        DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
            ROWS_FIELD_VARIABLE=>INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(matrix_idx)%VARIABLE
            COLS_FIELD_VARIABLE=>INTERFACE_MAPPING%LAGRANGE_VARIABLE !TEMPORARY: Needs generalising
            CALL EQUATIONS_MATRICES_ELEMENT_MATRIX_SETUP(INTERFACE_MATRIX%ELEMENT_MATRIX,ROWS_FIELD_VARIABLE, &
              & COLS_FIELD_VARIABLE,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Interface matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
              & " is not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ELSE
        CALL FLAG_ERROR("Interface matrices mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for an interface matrix.
  SUBROUTINE INTERFACE_MATRIX_STRUCTURE_CALCULATE(INTERFACE_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX !<A pointer to the interface matrix to calculate the strucute for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   INTEGER(INTG) ::  column_derivative,column_idx,column_component_idx,column_local_derivative_idx,column_local_node_idx, &
     & column_node,DUMMY_ERR,domain_element,domain_element_idx,global_column,interface_element_idx,INTERFACE_MESH_INDEX, &
     & local_column,local_row,MATRIX_NUMBER,NUMBER_OF_COLUMNS,row_component_idx,row_derivative,row_local_derivative_idx, &
     & row_local_node_idx,row_node
    INTEGER(INTG), POINTER :: COLUMNS(:)
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
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(COLUMNS)
       
    CALL ENTERS("INTERFACE_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*999)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
      IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
        IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
          MATRIX_NUMBER=INTERFACE_MATRIX%MATRIX_NUMBER
          SELECT CASE(INTERFACE_MATRIX%STRUCTURE_TYPE)
          CASE(INTERFACE_MATRIX_NO_STRUCTURE)
            CALL FLAG_ERROR("There is no structure to calculate for a matrix with no structure.",ERR,ERROR,*998)
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
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",ERR,ERROR,*999)
                                      DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                        !Set up list
                                        NULLIFY(COLUMN_INDICES_LISTS(local_row)%PTR)
                                        CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                        CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(local_row)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                                        CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(local_row)%PTR,50,ERR,ERROR,*999)
                                        CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(local_row)%PTR,ERR,ERROR,*999)
                                      ENDDO !local_row
                                      !Allocate row indices
                                      ALLOCATE(ROW_INDICES(ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",ERR,ERROR,*999)
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
                                                local_column=COLUMN_VARIABLE%COMPONENTS(column_component_idx)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP(column_derivative,column_node)
                                                global_column=COLUMN_DOFS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_column)
                                                !Loop over the components in the dependent variable
                                                DO row_component_idx=1,ROW_VARIABLE%NUMBER_OF_COMPONENTS
                                                  SELECT CASE(ROW_VARIABLE%COMPONENTS(row_component_idx)%INTERPOLATION_TYPE)
                                                  CASE(FIELD_CONSTANT_INTERPOLATION)
                                                    local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                      & CONSTANT_PARAM2DOF_MAP
                                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                      & ERR,ERROR,*999)                                               
                                                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                                                    domain_element=MESH_CONNECTIVITY% &
                                                      & ELEMENTS_CONNECTIVITY(interface_element_idx,INTERFACE_MESH_INDEX)% &
                                                      & COUPLED_MESH_ELEMENT_NUMBER
                                                    local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                      & ELEMENT_PARAM2DOF_MAP(domain_element)
                                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                      & ERR,ERROR,*999)
                                                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                                                    ROW_DOMAIN_ELEMENTS=>ROW_VARIABLE%COMPONENTS(row_component_idx)%DOMAIN% &
                                                      & TOPOLOGY%ELEMENTS
                                                    domain_element=MESH_CONNECTIVITY%ELEMENTS_CONNECTIVITY(interface_element_idx, &
                                                      & INTERFACE_MESH_INDEX)%COUPLED_MESH_ELEMENT_NUMBER
                                                    ROW_BASIS=>ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)%BASIS
                                                    !Loop over the row DOFs in the domain mesh element
                                                    DO row_local_node_idx=1,ROW_BASIS%NUMBER_OF_NODES
                                                      row_node=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                        & ELEMENT_NODES(row_local_node_idx)
                                                      DO row_local_derivative_idx=1,ROW_BASIS% &
                                                        & NUMBER_OF_DERIVATIVES(row_local_node_idx)
                                                        row_derivative=ROW_DOMAIN_ELEMENTS%ELEMENTS(domain_element)% &
                                                          & ELEMENT_DERIVATIVES(row_local_derivative_idx,row_local_node_idx)
                                                        local_row=ROW_VARIABLE%COMPONENTS(row_component_idx)%PARAM_TO_DOF_MAP% &
                                                          & NODE_PARAM2DOF_MAP(row_derivative,row_node)
                                                        CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(local_row)%PTR,global_column, &
                                                          & ERR,ERROR,*999)                                                      
                                                      ENDDO !row_local_derivative_idx
                                                    ENDDO !row_local_node_idx
                                                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                                                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                                                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                                  CASE DEFAULT
                                                    LOCAL_ERROR="The row variable interpolation type of "// &
                                                      & TRIM(NUMBER_TO_VSTRING(ROW_VARIABLE%COMPONENTS(row_component_idx)% &
                                                      INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid."
                                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                                  END SELECT
                                                ENDDO !row_component_idx
                                              ENDDO !column_local_derivative_idx
                                            ENDDO !column_local_node_idx
                                          ENDDO !interface_element_idx
                                        ELSE
                                          CALL FLAG_ERROR("Only node based fields implemented.",ERR,ERROR,*999)
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
                                      !Allocate and setup the column locations
                                      ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices.",ERR,ERROR,*999)
                                      DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                        CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_row)%PTR,NUMBER_OF_COLUMNS, &
                                          & COLUMNS,ERR,ERROR,*999)
                                        DO column_idx=1,NUMBER_OF_COLUMNS
                                          COLUMN_INDICES(ROW_INDICES(local_row)+column_idx-1)=COLUMNS(column_idx)
                                        ENDDO !column_idx
                                        DEALLOCATE(COLUMNS)
                                      ENDDO !local_ny
                                      IF(DIAGNOSTICS1) THEN
                                        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interface matrix structure:",ERR,ERROR,*999)
                                        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Interface matrix number : ",MATRIX_NUMBER, &
                                          & ERR,ERROR,*999)
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
                                          CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F6.2", &
                                            & ERR,ERROR,*999)
                                        ENDIF
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ROW_DOFS_DOMAIN_MAPPING% &
                                          & TOTAL_NUMBER_OF_LOCAL+1,8,8,ROW_INDICES,'("  Row indices    :",8(X,I13))', &
                                          & '(18X,8(X,I13))',ERR,ERROR,*999)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8, &
                                          & COLUMN_INDICES,'("  Column indices :",8(X,I13))','(18X,8(X,I13))', ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Column dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Row dofs parameter mapping is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Column dofs domain mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Row dofs domain mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Column field variable is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Row field variable is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Interface mesh connectivity is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Interface condition interface is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The matrix storage type of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The matrix structure type of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRIX%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Column indices is already associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indieces is already associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrix is not associated.",ERR,ERROR,*998)
    ENDIF
     
    CALL EXITS("INTERFACE_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO local_row=1,ROW_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(local_row)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(local_row)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !local_row
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
998 CALL ERRORS("INTERFACE_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRIX_STRUCTURE_CALCULATE")
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
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)

    NULLIFY(ROW_DOMAIN_MAP)
    NULLIFY(COLUMN_DOMAIN_MAP)

    CALL ENTERS("INTERFACE_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Interface matrices have already been finished.",ERR,ERROR,*998)
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
                  !Calculate and set the matrix structure/sparsity pattern
                  IF(INTERFACE_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                    & INTERFACE_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                    CALL INTERFACE_MATRIX_STRUCTURE_CALCULATE(INTERFACE_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(INTERFACE_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(INTERFACE_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                      & ERR,ERROR,*999)
                    IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                  ENDIF
                  CALL DISTRIBUTED_MATRIX_CREATE_FINISH(INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="Row domain map for interface matrix number "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is not associated."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Interface matrix for matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx                
            !Finish up
            INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Column domain map is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    CALL INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MATRICES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_CREATE_FINISH")
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

    CALL ENTERS("INTERFACE_MATRICES_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN      
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
          CALL FLAG_ERROR("Interface matrices is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(INTERFACE_MATRICES)
          !Initialise the interface matrices
          CALL INTERFACE_MATRICES_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_CREATE_START")
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

    CALL ENTERS("INTERFACE_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      CALL INTERFACE_MATRICES_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("INTERFACE_MATRICES_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_DESTROY",ERR,ERROR)    
    CALL EXITS("INTERFACE_MATRICES_DESTROY")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MATRICES_DESTROY
  !
  !================================================================================================================================
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
   
    CALL ENTERS("INTERFACE_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(ALLOCATED(INTERFACE_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(INTERFACE_MATRICES%MATRICES,1)
          CALL INTERFACE_MATRIX_FINALISE(INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(INTERFACE_MATRICES%MATRICES)
      ENDIF
      DEALLOCATE(INTERFACE_MATRICES)
    ENDIF
       
    CALL EXITS("INTERFACE_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_FINALISE")
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
    
    CALL ENTERS("INTERFACE_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MATRICES)) THEN
        CALL FLAG_ERROR("Interface matrices is already associated for this interface equations.",ERR,ERROR,*998)
      ELSE
        INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          IF(INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED) THEN
            ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface equations interface matrices.",ERR,ERROR,*999)
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED=.FALSE.
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%INTERFACE_MAPPING=>INTERFACE_MAPPING
            NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MATRICES%SOLVER_MAPPING)
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_COLUMNS=INTERFACE_MAPPING%NUMBER_OF_COLUMNS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%TOTAL_NUMBER_OF_COLUMNS=INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_GLOBAL_COLUMNS=INTERFACE_MAPPING%NUMBER_OF_GLOBAL_COLUMNS
            !Allocate and initialise the matrices
            INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
            ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(INTERFACE_EQUATIONS%INTERFACE_MATRICES% &
              & NUMBER_OF_INTERFACE_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface matrices matrices.",ERR,ERROR,*999)
            DO matrix_idx=1,INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
              NULLIFY(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR)
              CALL INTERFACE_MATRIX_INITIALISE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,matrix_idx,ERR,ERROR,*999)
            ENDDO !matrix_idx
          ELSE
            CALL FLAG_ERROR("Interface mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface equations interface mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MATRICES_INITIALISE")
    RETURN
999 CALL INTERFACE_MATRICES_FINALISE(INTERFACE_EQUATIONS%INTERFACE_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MATRICES_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_INITIALISE")
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
    
    CALL ENTERS("INTERFACE_MATRICES_OUTPUT",ERR,ERROR,*999)
    
    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"Interface matrices:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Number of interface matrices = ",INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,ERR,ERROR,*999)
        DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
          IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
            CALL WRITE_STRING_VALUE(ID,"Interface matrix : ",matrix_idx,ERR,ERROR,*999)
            CALL DISTRIBUTED_MATRIX_OUTPUT(ID,INTERFACE_MATRIX%MATRIX,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Interface matrix is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !matrix_idx
      ELSE
        CALL FLAG_ERROR("Interface matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_OUTPUT")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_OUTPUT",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_OUTPUT")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MATRICES_OUTPUT
  
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
    
    CALL ENTERS("INTERFACE_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Interface matrices have been finished.",ERR,ERROR,*999)
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
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of interface matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_STORAGE_TYPE_SET")
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

    CALL ENTERS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      IF(INTERFACE_MATRICES%INTERFACE_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Interface matrices have been finished.",ERR,ERROR,*999)
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
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Interface matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of interface matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_STRUCTURE_TYPE_SET")
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
    
    CALL ENTERS("INTERFACE_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
      DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
        INTERFACE_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR
        IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
          IF(INTERFACE_MATRIX%UPDATE_MATRIX) THEN
            CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(INTERFACE_MATRIX%MATRIX,VALUE,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface matrix is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
    ELSE
      CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("INTERFACE_MATRICES_VALUES_INITIALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MATRICES_VALUES_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MATRICES_VALUES_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MATRICES_VALUES_INITIALISE

  !
  !================================================================================================================================
  !


END MODULE INTERFACE_MATRICES_ROUTINES
