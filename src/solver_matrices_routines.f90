!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all solver matrix and rhs routines.
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

!> This module handles all solver matrix and rhs routines.
MODULE SOLVER_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters


  !> \addtogroup SOLVER_MATRICES_ROUTINES_SelectMatricesTypes SOLVER_MATRICES_ROUTINES::SelectMatricesTypes
  !> \brief The types of selection available for the solver matrices
  !> \see SOLVER_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_ALL=1 !<Select all the solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
!  redundant when introducing dynamic nonlinear equations
!  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_LINEAR_ONLY=3 !<Select only the linear solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear solver matrices and vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian solver matrix \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual solver vector \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_ONLY=7 !<Select only the RHS solver vector \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_RESIDUAL_ONLY=8 !<Select only the residual and RHS solver vectors \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_MATRICES_ALL,SOLVER_MATRICES_LINEAR_ONLY,SOLVER_MATRICES_NONLINEAR_ONLY, &
    & SOLVER_MATRICES_JACOBIAN_ONLY,SOLVER_MATRICES_RESIDUAL_ONLY,SOLVER_MATRICES_RHS_ONLY, & 
    & SOLVER_MATRICES_RHS_RESIDUAL_ONLY !,SOLVER_MATRICES_DYNAMIC_ONLY

  PUBLIC SOLVER_MATRIX_EQUATIONS_MATRIX_ADD,SOLVER_MATRIX_JACOBIAN_MATRIX_ADD
  
  PUBLIC SOLVER_MATRICES_CREATE_FINISH,SOLVER_MATRICES_CREATE_START,SOLVER_MATRICES_DESTROY,SOLVER_MATRICES_LIBRARY_TYPE_SET, &
    & SOLVER_MATRICES_OUTPUT,SOLVER_MATRICES_STORAGE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating the solver matrices
  SUBROUTINE SOLVER_MATRICES_CREATE_FINISH(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(COLUMN_INDICES)
    NULLIFY(ROW_INDICES)
    
    CALL ENTERS("SOLVER_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Solver matrices have already been finished",ERR,ERROR,*998)
      ELSE
        SOLVER_EQUATIONS=>SOLVER_MATRICES%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            !Now create the individual solver matrices
            ROW_DOMAIN_MAP=>SOLVER_MAPPING%ROW_DOFS_MAPPING
            IF(ASSOCIATED(ROW_DOMAIN_MAP)) THEN
              DO matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  COLUMN_DOMAIN_MAP=>SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(matrix_idx)%COLUMN_DOFS_MAPPING
                  IF(ASSOCIATED(COLUMN_DOMAIN_MAP)) THEN
                    !Create the distributed solver matrix
                    CALL DISTRIBUTED_MATRIX_CREATE_START(ROW_DOMAIN_MAP,COLUMN_DOMAIN_MAP,SOLVER_MATRICES%MATRICES(matrix_idx)% &
                         & PTR%MATRIX,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_LIBRARY_TYPE_SET(SOLVER_MATRIX%MATRIX,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(SOLVER_MATRIX%MATRIX,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(SOLVER_MATRIX%MATRIX,SOLVER_MATRIX%STORAGE_TYPE,ERR,ERROR,*999)
                    !Calculate and set the matrix structure/sparsity pattern
                    IF(SOLVER_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE.AND. &
                      & SOLVER_MATRIX%STORAGE_TYPE/=DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
                      CALL SOLVER_MATRIX_STRUCTURE_CALCULATE(SOLVER_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES, &
                        & COLUMN_INDICES,ERR,ERROR,*999)                  
                      CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(SOLVER_MATRIX%MATRIX,NUMBER_OF_NON_ZEROS, &
                        & ERR,ERROR,*999)
                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(SOLVER_MATRIX%MATRIX,ROW_INDICES,COLUMN_INDICES, &
                        & ERR,ERROR,*999)
                      IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                      IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                    ENDIF
                    CALL DISTRIBUTED_MATRIX_CREATE_FINISH(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
                    !Allocate the distributed solver vector
                    CALL DISTRIBUTED_VECTOR_CREATE_START(COLUMN_DOMAIN_MAP,SOLVER_MATRICES%MATRICES(matrix_idx)% &
                         & PTR%SOLVER_VECTOR,ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRIX%SOLVER_VECTOR,SOLVER_MATRICES%LIBRARY_TYPE, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRIX%SOLVER_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                    CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Column domain mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !matrix_idx
              IF(SOLVER_EQUATIONS%LINEARITY==PROBLEM_SOLVER_NONLINEAR) THEN
                !Allocate the nonlinear matrices and vectors                  
                !Allocate the distributed residual vector
                CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRICES%RESIDUAL,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRICES%RESIDUAL,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)                  
              ENDIF
!!TODO: what to do if there is no RHS
              !Allocate the distributed rhs vector
              CALL DISTRIBUTED_VECTOR_CREATE_START(ROW_DOMAIN_MAP,SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
              !Finish up
              SOLVER_MATRICES%SOLVER_MATRICES_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("Row domain mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver matrices solver equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    CALL SOLVER_MATRICES_FINALISE(SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRICES_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE SOLVER_MATRICES_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the solver matrices
  SUBROUTINE SOLVER_MATRICES_CREATE_START(SOLVER_EQUATIONS,SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to create the solver matrices for
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<On return, a pointer to the solver matrices. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_MATRICES)) THEN
          CALL FLAG_ERROR("Solver matrices is already associated",ERR,ERROR,*998)
        ELSE
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES)
          CALL SOLVER_MATRICES_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
        ENDIF
     ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_CREATE_START")
    RETURN
999 CALL SOLVER_MATRICES_FINALISE(SOLVER_EQUATIONS%SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRICES_CREATE_START",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_CREATE_START")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_CREATE_START
        
  !
  !================================================================================================================================
  !

  !>Destroy the solver matrices
  SUBROUTINE SOLVER_MATRICES_DESTROY(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer the solver matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      CALL SOLVER_MATRICES_FINALISE(SOLVER_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_DESTROY")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_DESTROY",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_DESTROY")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the solver matrices and deallocates all memory
  SUBROUTINE SOLVER_MATRICES_FINALISE(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("SOLVER_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(ALLOCATED(SOLVER_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(SOLVER_MATRICES%MATRICES,1)
          CALL SOLVER_MATRIX_FINALISE(SOLVER_MATRICES%MATRICES(matrix_idx)%PTR,ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(SOLVER_MATRICES%MATRICES)
      ENDIF
      IF(ASSOCIATED(SOLVER_MATRICES%RESIDUAL)) CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER_MATRICES%RHS_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(SOLVER_MATRICES)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Initialises the solver matrices for solver equations
  SUBROUTINE SOLVER_MATRICES_INITIALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to initialise the solver matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,equations_matrix_idx,equations_set_idx,matrix_idx
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MATRICES)) THEN
        CALL FLAG_ERROR("Solver matrices is already associated for this solver equations.",ERR,ERROR,*998)
      ELSE
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          ALLOCATE(SOLVER_EQUATIONS%SOLVER_MATRICES,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrices.",ERR,ERROR,*999)
          SOLVER_EQUATIONS%SOLVER_MATRICES%SOLVER_EQUATIONS=>SOLVER_EQUATIONS
          SOLVER_EQUATIONS%SOLVER_MATRICES%SOLVER_MATRICES_FINISHED=.FALSE.
          SOLVER_EQUATIONS%SOLVER_MATRICES%SOLVER_MAPPING=>SOLVER_MAPPING
          SOLVER_EQUATIONS%SOLVER_MATRICES%NUMBER_OF_ROWS=SOLVER_MAPPING%NUMBER_OF_ROWS
          SOLVER_EQUATIONS%SOLVER_MATRICES%NUMBER_OF_GLOBAL_ROWS=SOLVER_MAPPING%NUMBER_OF_GLOBAL_ROWS
          SOLVER_EQUATIONS%SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
          SOLVER_EQUATIONS%SOLVER_MATRICES%NUMBER_OF_MATRICES=SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
          ALLOCATE(SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrices matrices.",ERR,ERROR,*999)
          DO matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
            NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(matrix_idx)%PTR)
            CALL SOLVER_MATRIX_INITIALISE(SOLVER_EQUATIONS%SOLVER_MATRICES,matrix_idx,ERR,ERROR,*999)
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              IF(ALLOCATED(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                & matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS)) THEN
                DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                  !Add the solver matrix to the solvers mapping
                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                    & matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR%SOLVER_MATRIX=> &
                    & SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
                ENDDO !equations_matrix_idx
              ELSE
                IF(ASSOCIATED(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                  & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP)) THEN
                  SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP%SOLVER_MATRIX=> &
                    & SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
                ELSE
                  DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                    !Add the solver matrix to the solvers mapping
                    SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM( &
                      & matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS(equations_matrix_idx)%PTR%SOLVER_MATRIX=> &
                      & SOLVER_EQUATIONS%SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
                  ENDDO !equations_matrix_idx
                ENDIF
              ENDIF
            ENDDO !equations_set_idx
          ENDDO !matrix_idx
          IF(SOLVER_EQUATIONS%LINEARITY==PROBLEM_SOLVER_NONLINEAR) THEN
            SOLVER_EQUATIONS%SOLVER_MATRICES%UPDATE_RESIDUAL=.TRUE.
          ELSE
            SOLVER_EQUATIONS%SOLVER_MATRICES%UPDATE_RESIDUAL=.FALSE.
          ENDIF
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES%RESIDUAL)
          SOLVER_EQUATIONS%SOLVER_MATRICES%UPDATE_RHS_VECTOR=.TRUE.
          NULLIFY(SOLVER_EQUATIONS%SOLVER_MATRICES%RHS_VECTOR)
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_INITIALISE")
    RETURN
999 CALL SOLVER_MATRICES_FINALISE(SOLVER_EQUATIONS%SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRICES_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_INITIALISE
  
  !
  !================================================================================================================================
  !
  
  !>Gets the library type for the solver matrices (and vectors)
  SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_GET(SOLVER_MATRICES,LIBRARY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(OUT) :: LIBRARY_TYPE !<On return, the library type of the specified solver matrices \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MATRICES_LIBRARY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        LIBRARY_TYPE=SOLVER_MATRICES%LIBRARY_TYPE
      ELSE
        CALL FLAG_ERROR("Solver matrices has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_LIBRARY_TYPE_GET")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_LIBRARY_TYPE_GET",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_LIBRARY_TYPE_GET")
    RETURN
  END SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_GET
  
        
  !
  !================================================================================================================================
  !

  !>Sets the library type for the solver matrices (and vectors)
  SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,LIBRARY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices.
    INTEGER(INTG), INTENT(IN) :: LIBRARY_TYPE !<The library type to set \see DISTRIBUTED_MATRIX_VECTOR_LibraryTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_MATRICES_LIBRARY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Solver matrices has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(LIBRARY_TYPE)
        CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
          SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE
        CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
          SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The library type of "// TRIM(NUMBER_TO_VSTRING(LIBRARY_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_LIBRARY_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_LIBRARY_TYPE_SET",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_LIBRARY_TYPE_SET")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_LIBRARY_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Outputs the solver matrices
  SUBROUTINE SOLVER_MATRICES_OUTPUT(ID,SELECTION_TYPE,SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The type of matrix selection \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    
    CALL ENTERS("SOLVER_MATRICES_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
!           & SELECTION_TYPE==SOLVER_MATRICES_DYNAMIC_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
          CALL WRITE_STRING(ID,"Solver matrices:",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(ID,"Number of matrices = ",SOLVER_MATRICES%NUMBER_OF_MATRICES,ERR,ERROR,*999)
          DO matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              CALL WRITE_STRING_VALUE(ID,"Solver matrix : ",matrix_idx,ERR,ERROR,*999) 
              CALL DISTRIBUTED_MATRIX_OUTPUT(ID,SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)    
            ELSE
              CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ENDIF
        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
          IF(ASSOCIATED(SOLVER_MATRICES%RESIDUAL)) THEN
            CALL WRITE_STRING(ID,"Solver residual vector:",ERR,ERROR,*999)     
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,SOLVER_MATRICES%RESIDUAL,ERR,ERROR,*999)  
          ENDIF
        ENDIF
        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
!          & SELECTION_TYPE==SOLVER_MATRICES_DYNAMIC_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RHS_ONLY.OR. &
          & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
          IF(ASSOCIATED(SOLVER_MATRICES%RHS_VECTOR)) THEN
            CALL WRITE_STRING(ID,"Solver RHS vector:",ERR,ERROR,*999)     
            CALL DISTRIBUTED_VECTOR_OUTPUT(ID,SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver matrices have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_OUTPUT")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_OUTPUT",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_OUTPUT")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_OUTPUT
  
  !
  !================================================================================================================================
  !
  
  !>Gets the storage type (sparsity) of the solver matrices
  SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_GET(SOLVER_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). On return, the storage type for the matrix_idx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_MATRICES_STORAGE_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        IF(SIZE(STORAGE_TYPE,1)>=SOLVER_MATRICES%NUMBER_OF_MATRICES) THEN
          DO matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              STORAGE_TYPE(matrix_idx)=SOLVER_MATRIX%STORAGE_TYPE
            ELSE
              CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of STORAGE_TYPE is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver matrices have not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_STORAGE_TYPE_GET")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_STORAGE_TYPE_GET",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_STORAGE_TYPE_GET")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the solver matrices
  SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th solver matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Solver matrices have been finished.",ERR,ERROR,*999)
      ELSE
        IF(SIZE(STORAGE_TYPE,1)==SOLVER_MATRICES%NUMBER_OF_MATRICES) THEN
          DO matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrix_idx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              SELECT CASE(STORAGE_TYPE(matrix_idx))
              CASE(MATRIX_BLOCK_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
              CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_DIAGONAL_STORAGE_TYPE        
              CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
              CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_ROW_MAJOR_STORAGE_TYPE
              CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
              CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
              CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
                SOLVER_MATRIX%STORAGE_TYPE=MATRIX_ROW_COLUMN_STORAGE_TYPE
              CASE DEFAULT
                LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                  & " for the matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !matrix_idx
        ELSE
          LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
            & ") is not equal to the number of matrices ("// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Adds alpha times the equations matrix into the solver matrix
  SUBROUTINE SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx,ALPHA,EQUATIONS_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: equations_set_idx !<The equations set index in the solver mapping that contains the equations matrix to add
    REAL(DP), INTENT(IN) :: ALPHA !<The multiplicative factor for the equations matrix
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX !<A pointer to the equations matrix to add    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,equations_row_number,EQUATIONS_STORAGE_TYPE, &
      & solver_column_idx,solver_column_number,solver_row_idx,solver_row_number
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(DP) :: column_coupling_coefficient,row_coupling_coefficient,VALUE
    REAL(DP), POINTER :: EQUATIONS_MATRIX_DATA(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
        IF(ABS(ALPHA)>ZERO_TOLERANCE) THEN
          SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
              SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
                DYNAMIC_MATRICES=>EQUATIONS_MATRIX%DYNAMIC_MATRICES
                IF(ASSOCIATED(DYNAMIC_MATRICES).OR.ASSOCIATED(LINEAR_MATRICES)) THEN
                  IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                    EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
                  ELSE
                    EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
                  ENDIF
                  IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                    IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
                      IF(equations_set_idx>0.AND.equations_set_idx<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS) THEN
                        EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_EM(EQUATIONS_MATRIX%MATRIX_NUMBER)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS(SOLVER_MATRIX%MATRIX_NUMBER)%PTR
                        IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                          SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                          IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN
                            EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                            IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                              CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE, &
                                & ERR,ERROR,*999)
                              CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA,ERR,ERROR,*999)
                              SELECT CASE(EQUATIONS_STORAGE_TYPE)
                              CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                !Loop over the rows of the equations matrix
                                DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this equations row is mapped to
                                  DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    !Loop over the columns of the equations matrix
                                    DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                      !Loop over the solution columns this equations column is mapped to
                                      DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                        & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                        solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                          & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                        column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                          & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                          & COUPLING_COEFFICIENTS(solver_column_idx)
                                        !Add in the solver matrix value
                                        VALUE=ALPHA*EQUATIONS_MATRIX_DATA(equations_row_number+ &
                                          & (equations_column_number-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                                          & row_coupling_coefficient*column_coupling_coefficient
                                        CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX, &
                                          & solver_row_number,solver_column_number,VALUE,ERR,ERROR,*999)
                                      ENDDO !solver_column_idx
                                    ENDDO !equations_column_number
                                  ENDDO !solver_row_idx
                                ENDDO !equations_row_number
                              CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                !Loop over the rows of the equations matrix
                                DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this equations row is mapped to
                                  DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    equations_column_number=equations_row_number
                                    !Loop over the solution columns this equations column is mapped to
                                    DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                      & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                      solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                        & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                      column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                        & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                        & COUPLING_COEFFICIENTS(solver_column_idx)
                                      !Add in the solver matrix value
                                      VALUE=ALPHA*EQUATIONS_MATRIX_DATA(equations_row_number)* &
                                        & row_coupling_coefficient*column_coupling_coefficient
                                      CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX, &
                                        & solver_row_number,solver_column_number,VALUE,ERR,ERROR,*999)
                                    ENDDO !solver_column_idx
                                  ENDDO !solver_row_idx
                                ENDDO !equations_row_number
                              CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                              CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                              CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                  & ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*999)
                                !Loop over the rows of the equations matrix
                                DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this equations row is mapped to
                                  DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    !Loop over the columns of the equations matrix
                                    DO equations_column_idx=ROW_INDICES(equations_row_number),ROW_INDICES(equations_row_number+1)-1
                                      equations_column_number=COLUMN_INDICES(equations_column_idx)
                                      !Loop over the solution columns this equations column is mapped to
                                      DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                        & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                        solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                          & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                        column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                          & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                          & COUPLING_COEFFICIENTS(solver_column_idx)
                                        !Add in the solver matrix value
                                        VALUE=ALPHA*EQUATIONS_MATRIX_DATA(equations_column_idx)*row_coupling_coefficient* &
                                          & column_coupling_coefficient
                                        CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX, &
                                          & solver_row_number,solver_column_number,VALUE,ERR,ERROR,*999)
                                      ENDDO !solution_column_idx
                                    ENDDO !equations_column_idx
                                  ENDDO !solution_row_idx
                                ENDDO !equations_row_number
                              CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                              CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                              CASE DEFAULT
                                LOCAL_ERROR="The equations matrix storage type of "// &
                                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                              CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                & ERR,ERROR,*999)
                            ELSE
                              CALL FLAG_ERROR("The equations matrix distributed matrix is not associated",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations to solver map is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The specified equations set index of "// &
                          & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                          & " is invalid. The equations set index needs to be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations matrices have not been finished.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Dynamic or linear matrices equations matrices is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations matrix dynamic or linear matrices is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrices solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver matrices have not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver matrix solver matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD")
    RETURN
999 CALL ERRORS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD",ERR,ERROR)
    CALL EXITS("SOLVER_MATRIX_EQUATIONS_MATRIX_ADD")
    RETURN 1
  END SUBROUTINE SOLVER_MATRIX_EQUATIONS_MATRIX_ADD

  !
  !================================================================================================================================
  !

  !>Adds alpha times the Jacobian matrix into the solver matrix
  SUBROUTINE SOLVER_MATRIX_JACOBIAN_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx,ALPHA,JACOBIAN_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix
    INTEGER(INTG), INTENT(IN) :: equations_set_idx !<The equations set index in the solver mapping that contains the Jacobian matrix to add
    REAL(DP), INTENT(IN) :: ALPHA !<The multiplicative factor for the Jacobian matrix
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX !<A pointer to the Jacobian matrix to add    
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobian_column_idx,jacobian_column_number,jacobian_row_number,JACOBIAN_STORAGE_TYPE, &
      & solver_column_idx,solver_column_number,solver_row_idx,solver_row_number
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(DP) :: column_coupling_coefficient,row_coupling_coefficient,VALUE
    REAL(DP), POINTER :: JACOBIAN_MATRIX_DATA(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: JACOBIAN_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      NULLIFY(SOLVER_MATRICES)
      NULLIFY(SOLVER_MAPPING)
      NULLIFY(NONLINEAR_MATRICES)
      NULLIFY(EQUATIONS_MATRICES)
      NULLIFY(JACOBIAN_TO_SOLVER_MAP)
      NULLIFY(SOLVER_DISTRIBUTED_MATRIX)
      NULLIFY(JACOBIAN_DISTRIBUTED_MATRIX)

      IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
        IF(ABS(ALPHA)>ZERO_TOLERANCE) THEN
          SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
              SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
                IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                  EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
                  IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                    IF(EQUATIONS_MATRICES%EQUATIONS_MATRICES_FINISHED) THEN
                      IF(equations_set_idx>0.AND.equations_set_idx<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS) THEN
                        JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(SOLVER_MATRIX%MATRIX_NUMBER)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                        IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                          SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                          IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN
                            JACOBIAN_DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
                            IF(ASSOCIATED(JACOBIAN_DISTRIBUTED_MATRIX)) THEN
                              CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_STORAGE_TYPE, &
                                & ERR,ERROR,*999)
                              CALL DISTRIBUTED_MATRIX_DATA_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA,ERR,ERROR,*999)
                              SELECT CASE(JACOBIAN_STORAGE_TYPE)
                              CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                !Loop over the rows of the Jacobian matrix
                                DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this Jacobian row is mapped to
                                  DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    !Loop over the columns of the Jacobian matrix
                                    DO jacobian_column_number=1,JACOBIAN_MATRIX%NUMBER_OF_COLUMNS
                                      !Loop over the solution columns this Jacobian column is mapped to
                                      DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                        & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                        solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                          & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                        column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                          & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                          & COUPLING_COEFFICIENTS(solver_column_idx)
                                        !Add in the solver matrix value
                                        VALUE=ALPHA*JACOBIAN_MATRIX_DATA(jacobian_row_number+(jacobian_column_number-1)* &
                                          & EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)*row_coupling_coefficient* &
                                          & column_coupling_coefficient
                                        CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                          & solver_column_number,VALUE,ERR,ERROR,*999)
                                      ENDDO !solver_column_idx
                                    ENDDO !jacobian_column_number
                                  ENDDO !solver_row_idx
                                ENDDO !jacobian_row_number
                              CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                !Loop over the rows of the Jacobian matrix
                                DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this Jacobian row is mapped to
                                  DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    jacobian_column_number=jacobian_row_number
                                    !Loop over the solution columns this Jacobian column is mapped to
                                    DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                      & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                      solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                        & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                      column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                        & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                        & COUPLING_COEFFICIENTS(solver_column_idx)
                                      !Add in the solver matrix value
                                      VALUE=ALPHA*JACOBIAN_MATRIX_DATA(jacobian_row_number)* &
                                        & row_coupling_coefficient*column_coupling_coefficient
                                      CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX, &
                                        & solver_row_number,solver_column_number,VALUE,ERR,ERROR,*999)
                                    ENDDO !solver_column_idx
                                  ENDDO !solver_row_idx
                                ENDDO !jacobian_row_number
                              CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                              CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                              CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(JACOBIAN_DISTRIBUTED_MATRIX,ROW_INDICES, &
                                  & COLUMN_INDICES,ERR,ERROR,*999)
                                !Loop over the rows of the Jacobian matrix
                                DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this Jacobian row is mapped to
                                  DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    !Loop over the columns of the Jacobian matrix
                                    DO jacobian_column_idx=ROW_INDICES(jacobian_row_number), &
                                      & ROW_INDICES(jacobian_row_number+1)-1
                                      jacobian_column_number=COLUMN_INDICES(jacobian_column_idx)
                                      !Loop over the solution columns this equations column is mapped to
                                      DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                        & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                        solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                          & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                        column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                          & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                          & COUPLING_COEFFICIENTS(solver_column_idx)
                                        !Add in the solver matrix value
                                        VALUE=ALPHA*JACOBIAN_MATRIX_DATA(jacobian_column_idx)*row_coupling_coefficient* &
                                          & column_coupling_coefficient
                                        CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                          & solver_column_number,VALUE,ERR,ERROR,*999)                                    
                                      ENDDO !solution_column_idx
                                    ENDDO !jacobian_column_idx
                                  ENDDO !solution_row_idx
                                ENDDO !jacobian_row_number
                              CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                              CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                              CASE DEFAULT
                                LOCAL_ERROR="The Jacobian matrix storage type of "// &
                                  & TRIM(NUMBER_TO_VSTRING(JACOBIAN_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                              CALL DISTRIBUTED_MATRIX_DATA_RESTORE(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA, &
                                & ERR,ERROR,*999)
                            ELSE
                              CALL FLAG_ERROR("The Jacobian matrix distributed matrix is not associated",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Jacobian to solver map is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="The specified equations set index of "// &
                          & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                          & " is invalid. The equations set index needs to be between 1 and "// &
                          & TRIM(NUMBER_TO_VSTRING(SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations matrices have not been finished.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Jacobian matrix nonlinear matrices is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrices solver mapping is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver matrices have not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver matrix solver matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD")
    RETURN
999 CALL ERRORS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD",ERR,ERROR)
    CALL EXITS("SOLVER_MATRIX_JACOBIAN_MATRIX_ADD")
    RETURN 1
  END SUBROUTINE SOLVER_MATRIX_JACOBIAN_MATRIX_ADD

  !
  !================================================================================================================================
  !

  !>Calculates the structure (sparsity) of the solver matrix from the soluton mapping.
  SUBROUTINE SOLVER_MATRIX_STRUCTURE_CALCULATE(SOLVER_MATRIX,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix to calculate the structure for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the solver matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices in compressed row format. The pointers must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices in compressed row format. The pointers must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,DUMMY_ERR,equations_matrix_idx,equations_row_number, &
      & equations_set_idx,EQUATIONS_STORAGE_TYPE,jacobian_column_idx,jacobian_column_number,jacobian_row_number, &
      & MAX_COLUMN_INDICES,MAX_COLUMNS_PER_ROW,NUMBER_OF_COLUMNS,solver_column_idx,solver_column_number,solver_matrix_idx, &
      & solver_row_idx,solver_row_number
    INTEGER(INTG), POINTER :: COLUMNS(:),EQUATIONS_ROW_INDICES(:),EQUATIONS_COLUMN_INDICES(:)
    REAL(DP) :: SPARSITY
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    NULLIFY(COLUMNS)
 
    CALL ENTERS("SOLVER_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*999)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
        IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
          SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
          IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN
            IF(SOLVER_DISTRIBUTED_MATRIX%MATRIX_FINISHED) THEN
              CALL FLAG_ERROR("The solver distributed matrix has already been finished.",ERR,ERROR,*998)        
            ELSE
              SOLVER_MATRICES=>SOLVER_MATRIX%SOLVER_MATRICES             
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  SELECT CASE(SOLVER_MATRIX%STORAGE_TYPE)
                  CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                    CALL FLAG_ERROR("Can not calcualte the structure for a block storage matrix.",ERR,ERROR,*999)
                  CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                    CALL FLAG_ERROR("Can not calcualte the structure for a diagonal matrix.",ERR,ERROR,*999)
                  CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                  CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                    solver_matrix_idx=SOLVER_MATRIX%MATRIX_NUMBER
                    !Allocate lists
                    ALLOCATE(COLUMN_INDICES_LISTS(SOLVER_MAPPING%NUMBER_OF_ROWS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists.",ERR,ERROR,*999)
                    !Allocate row indices
                    ALLOCATE(ROW_INDICES(SOLVER_MAPPING%NUMBER_OF_ROWS+1),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices.",ERR,ERROR,*999)
                    ROW_INDICES(1)=1
                    !Find the maximum number of column indices
                    MAX_COLUMN_INDICES=0
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%PTR
                          IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                            EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                              DYNAMIC_MATRICES=>EQUATIONS_MATRIX%DYNAMIC_MATRICES
                              IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                                EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                  IF(ASSOCIATED(DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,MAX_COLUMNS_PER_ROW, &
                                      & ERR,ERROR,*999)
                                    MAX_COLUMN_INDICES=MAX_COLUMN_INDICES+MAX_COLUMNS_PER_ROW
                                  ELSE
                                    CALL FLAG_ERROR("Equations matrix distributed matrix is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Dynamic matrices equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations matrix dynamic matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Equations matrix is not assocaited.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations to solver matrix map is not assocaited.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !equations_matrix_idx
                      ELSE
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%PTR
                          IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                            EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                              LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
                              IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                                EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                  IF(ASSOCIATED(DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,MAX_COLUMNS_PER_ROW, &
                                      & ERR,ERROR,*999)
                                    MAX_COLUMN_INDICES=MAX_COLUMN_INDICES+MAX_COLUMNS_PER_ROW
                                  ELSE
                                    CALL FLAG_ERROR("Equations matrix distributed matrix is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Linear matrices equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations matrix linear matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Equations matrix is not assocaited.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations to solver matrix map is not assocaited.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !equations_matrix_idx
                        JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                        IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                          JACOBIAN_MATRIX=>JACOBIAN_TO_SOLVER_MAP%JACOBIAN_MATRIX
                          IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                            NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
                            IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                              EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
                              IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
                                IF(ASSOCIATED(DISTRIBUTED_MATRIX)) THEN
                                  CALL DISTRIBUTED_MATRIX_MAX_COLUMNS_PER_ROW_GET(DISTRIBUTED_MATRIX,MAX_COLUMNS_PER_ROW, &
                                    & ERR,ERROR,*999)
                                  MAX_COLUMN_INDICES=MAX_COLUMN_INDICES+MAX_COLUMNS_PER_ROW
                                ELSE
                                  CALL FLAG_ERROR("Jacobian distributed matrix is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Jacobian matrix nonlinear matrices is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Jacobian matrix is not assocaited.",ERR,ERROR,*999)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDDO !equations_set_idx
                    !Set up the column indicies lists
                    DO solver_row_number=1,SOLVER_MAPPING%NUMBER_OF_ROWS
                      NULLIFY(COLUMN_INDICES_LISTS(solver_row_number)%PTR)
                      CALL LIST_CREATE_START(COLUMN_INDICES_LISTS(solver_row_number)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(COLUMN_INDICES_LISTS(solver_row_number)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(COLUMN_INDICES_LISTS(solver_row_number)%PTR,MAX_COLUMN_INDICES,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(COLUMN_INDICES_LISTS(solver_row_number)%PTR,ERR,ERROR,*999)
                    ENDDO !solver_row_number
                    !Loop over the equations sets
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES>0) THEN
                        !Loop over the dynamic equations matrices mapped to the solver matrix and calculate the col indices by row.
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_DYNAMIC_EQUATIONS_MATRICES
                          !Note: pointers have been checked above
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%DYNAMIC_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%PTR
                          EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                          DYNAMIC_MATRICES=>EQUATIONS_MATRIX%DYNAMIC_MATRICES
                          EQUATIONS_MATRICES=>DYNAMIC_MATRICES%EQUATIONS_MATRICES
                          DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE,ERR,ERROR,*999)
                          SELECT CASE(EQUATIONS_STORAGE_TYPE)
                          CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                            !Loop over the rows of the equations matrix
                            DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)                                              
                                !Loop over the columns of the equations matrix
                                DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                  !Loop over the solver columns this equations column is mapped to
                                  DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                    & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                    solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                      & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                      & ERR,ERROR,*999)
                                  ENDDO !solver_column_idx
                                ENDDO !equations_column_number
                              ENDDO !solver_row_idx
                            ENDDO !equations_row_number                          
                          CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                            !Loop over the rows of the equations matrix
                            DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)                                              
                                equations_column_number=equations_row_number
                                !Loop over the solver columns this equations column is mapped to
                                DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                  & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                  solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                    & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                  CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                    & ERR,ERROR,*999)
                                ENDDO !solver_column_idx
                              ENDDO !solver_row_idx
                            ENDDO !equations_row_number                          
                          CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                          CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX,EQUATIONS_ROW_INDICES, &
                              & EQUATIONS_COLUMN_INDICES,ERR,ERROR,*999)
                            !Loop over the rows of the equations matrix
                            DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)
                                !Loop over the columns of the equations matrix
                                DO equations_column_idx=EQUATIONS_ROW_INDICES(equations_row_number), &
                                  & EQUATIONS_ROW_INDICES(equations_row_number+1)-1
                                  equations_column_number=EQUATIONS_COLUMN_INDICES(equations_column_idx)
                                  !Loop over the solver columns this equations column is mapped to
                                  DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                    & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                    solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                      & equations_column_number)%SOLVER_COLS(solver_column_idx)                                   
                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                      & ERR,ERROR,*999)
                                  ENDDO !solver_column_idx
                                ENDDO !equations_column_idx
                              ENDDO !equations_row_idx
                            ENDDO !equations_row_number                          
                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                          CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            LOCAL_ERROR="The matrix storage type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ENDDO !equations_matrix_idx
                      ELSE
                        !Loop over the linear equations matrices mapped to the solver matrix and calculate the col indices by row.
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          !Note: pointers have been checked above
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%PTR
                          EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                          LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
                          EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
                          DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE,ERR,ERROR,*999)
                          SELECT CASE(EQUATIONS_STORAGE_TYPE)
                          CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                            !Loop over the rows of the equations matrix
                            DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)                                              
                                !Loop over the columns of the equations matrix
                                DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                  !Loop over the solver columns this equations column is mapped to
                                  DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                    & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                    solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                      & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                      & ERR,ERROR,*999)
                                  ENDDO !solver_column_idx
                                ENDDO !equations_column_number
                              ENDDO !solver_row_idx
                            ENDDO !equations_row_number                          
                          CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                            !Loop over the rows of the equations matrix
                            DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)                                              
                                equations_column_number=equations_row_number
                                !Loop over the solver columns this equations column is mapped to
                                DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                  & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                  solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                    & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                  CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                    & ERR,ERROR,*999)
                                ENDDO !solver_column_idx
                              ENDDO !solver_row_idx
                            ENDDO !equations_row_number 
                          CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                          CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX,EQUATIONS_ROW_INDICES, &
                              & EQUATIONS_COLUMN_INDICES,ERR,ERROR,*999)
                            !Loop over the rows of the equations matrix
                            DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)
                                !Loop over the columns of the equations matrix
                                DO equations_column_idx=EQUATIONS_ROW_INDICES(equations_row_number), &
                                  & EQUATIONS_ROW_INDICES(equations_row_number+1)-1
                                  equations_column_number=EQUATIONS_COLUMN_INDICES(equations_column_idx)
                                  !Loop over the solver columns this equations column is mapped to
                                  DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                    & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                    solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                      & equations_column_number)%SOLVER_COLS(solver_column_idx)                                   
                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                      & ERR,ERROR,*999)
                                  ENDDO !solver_column_idx
                                ENDDO !equations_column_idx
                              ENDDO !equations_row_idx
                            ENDDO !equations_row_number                          
                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                          CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            LOCAL_ERROR="The matrix storage type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ENDDO !equations_matrix_idx
                        !Now add any columns from the Jacobian
                        JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                        IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                          !Note: pointers have been checked above
                          JACOBIAN_MATRIX=>JACOBIAN_TO_SOLVER_MAP%JACOBIAN_MATRIX
                          NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
                          EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
                          DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
                          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE,ERR,ERROR,*999)
                          SELECT CASE(EQUATIONS_STORAGE_TYPE)
                          CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                            !Loop over the rows of the Jacobian matrix
                            DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)                                              
                                !Loop over the columns of the Jacobian
                                DO jacobian_column_number=1,JACOBIAN_MATRIX%NUMBER_OF_COLUMNS
                                  !Loop over the solver columns this equations column is mapped to
                                  DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                    & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                    solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                      & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                      & ERR,ERROR,*999)
                                  ENDDO !solver_column_idx
                                ENDDO !jacobian_column_number
                              ENDDO !solver_row_idx
                            ENDDO !jacobian_row_number                          
                          CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                            !Loop over the rows of the Jacobian matrix
                            DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)                                              
                                jacobian_column_number=jacobian_row_number
                                !Loop over the solver columns this equations column is mapped to
                                DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                  & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                  solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                    & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                  CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                    & ERR,ERROR,*999)
                                ENDDO !solver_column_idx
                              ENDDO !solver_row_idx
                            ENDDO !jacobian_row_number                          
                          CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                          CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(DISTRIBUTED_MATRIX,EQUATIONS_ROW_INDICES, &
                              & EQUATIONS_COLUMN_INDICES,ERR,ERROR,*999)
                            !Loop over the rows of the Jacobian matrix
                            DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                              !Loop over the solver rows this equations row is mapped to
                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                  & SOLVER_ROWS(solver_row_idx)
                                !Loop over the columns of the Jacobian matrix
                                DO jacobian_column_idx=EQUATIONS_ROW_INDICES(jacobian_row_number), &
                                  & EQUATIONS_ROW_INDICES(jacobian_row_number+1)-1
                                  jacobian_column_number=EQUATIONS_COLUMN_INDICES(jacobian_column_idx)
                                  !Loop over the solver columns this Jacobian column is mapped to
                                  DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                    & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                    solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                      & jacobian_column_number)%SOLVER_COLS(solver_column_idx)                                   
                                    CALL LIST_ITEM_ADD(COLUMN_INDICES_LISTS(solver_row_number)%PTR,solver_column_number, &
                                      & ERR,ERROR,*999)
                                  ENDDO !solver_column_idx
                                ENDDO !jacobian_column_idx
                              ENDDO !solver_row_idx
                            ENDDO !jacobian_row_number                          
                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                          CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            LOCAL_ERROR="The Jacobian storage type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ENDIF
                      ENDIF
                    ENDDO !equations_set_idx
                    !Loop over the rows to calculate the number of non-zeros and setup the row indicces
                    DO solver_row_number=1,SOLVER_MAPPING%NUMBER_OF_ROWS
                      CALL LIST_REMOVE_DUPLICATES(COLUMN_INDICES_LISTS(solver_row_number)%PTR,ERR,ERROR,*999)
                      CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(solver_row_number)%PTR,NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                      NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                      ROW_INDICES(solver_row_number+1)=NUMBER_OF_NON_ZEROS+1
                    ENDDO !solver_row_number
                    !Allocate and setup the column locations
                    ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices",ERR,ERROR,*999)
                    DO solver_row_number=1,SOLVER_MAPPING%NUMBER_OF_ROWS
                      CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(solver_row_number)%PTR,NUMBER_OF_COLUMNS,COLUMNS, &
                        & ERR,ERROR,*999)
                      DO solver_column_idx=1,NUMBER_OF_COLUMNS
                        COLUMN_INDICES(ROW_INDICES(solver_row_number)+solver_column_idx-1)=COLUMNS(solver_column_idx)
                      ENDDO !solver_column_idx
                      DEALLOCATE(COLUMNS)
                    ENDDO !solver_row_idx
                  CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                  CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                  CASE DEFAULT
                    LOCAL_ERROR="The matrix storage type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRIX%STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT

                  IF(DIAGNOSTICS1) THEN
                    CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix structure:",ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Solver matrix number : ",SOLVER_MATRIX%MATRIX_NUMBER, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",SOLVER_MATRICES%NUMBER_OF_ROWS, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of columns = ",SOLVER_MATRIX%NUMBER_OF_COLUMNS, &
                      & ERR,ERROR,*999)
                    CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of non zeros = ",NUMBER_OF_NON_ZEROS,ERR,ERROR,*999)
                    IF(SOLVER_MATRICES%NUMBER_OF_ROWS*SOLVER_MATRIX%NUMBER_OF_COLUMNS/=0) THEN
                      SPARSITY=REAL(NUMBER_OF_NON_ZEROS,DP)/REAL(SOLVER_MATRICES%NUMBER_OF_ROWS* &
                        & SOLVER_MATRIX%NUMBER_OF_COLUMNS,DP)*100.0_DP
                      CALL WRITE_STRING_FMT_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sparsity (%) = ",SPARSITY,"F6.2", ERR,ERROR,*999)
                    ENDIF
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLVER_MATRICES%NUMBER_OF_ROWS+1,8,8,ROW_INDICES, &
                      & '("  Row indices    :",8(X,I13))','(18X,8(X,I13))',ERR,ERROR,*999)
                    CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_NON_ZEROS,8,8,COLUMN_INDICES, &
                      & '("  Column indices :",8(X,I13))','(18X,8(X,I13))', ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrices solver mapping is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrix solver matrices is not associated",ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver matrix distributed matrix is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Column indices is already associated",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Row indices is already associated",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("SOLVER_MATRIX_STRUCTURE_CALCULATE")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    IF(ASSOCIATED(COLUMNS)) DEALLOCATE(COLUMNS)
    IF(ALLOCATED(COLUMN_INDICES_LISTS)) THEN
      DO solver_row_number=1,SOLVER_MAPPING%NUMBER_OF_ROWS
        IF(ASSOCIATED(COLUMN_INDICES_LISTS(solver_row_number)%PTR)) &
          & CALL LIST_DESTROY(COLUMN_INDICES_LISTS(solver_row_number)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO !solver_row_number
      DEALLOCATE(COLUMN_INDICES_LISTS)
    ENDIF
998 CALL ERRORS("SOLVER_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRIX_STRUCTURE_CALCULATE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_STRUCTURE_CALCULATE
        
  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix and deallocates all memory
  SUBROUTINE SOLVER_MATRIX_FINALISE(SOLVER_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      IF(ASSOCIATED(SOLVER_MATRIX%MATRIX)) CALL DISTRIBUTED_MATRIX_DESTROY(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER_MATRIX%SOLVER_VECTOR)) CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
      DEALLOCATE(SOLVER_MATRIX)
    ENDIF
    
    CALL EXITS("SOLVER_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MATRIX_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRIX_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Forms a solver matrix by initialising the structure of the matrix to zero.
  SUBROUTINE SOLVER_MATRIX_FORM(SOLVER_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX !<A pointer to the solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MATRIX_FORM",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
      CALL DISTRIBUTED_MATRIX_FORM(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRIX_FORM")
    RETURN
999 CALL ERRORS("SOLVER_MATRIX_FORM",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRIX_FORM")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FORM
        
    !
  !================================================================================================================================
  !

  !>Initialises a solver matrix
  SUBROUTINE SOLVER_MATRIX_INITIALISE(SOLVER_MATRICES,MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices to initialise
    INTEGER(INTG), INTENT(IN) :: MATRIX_NUMBER !<The matrix number in the solver matrices to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_MATRIX_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(MATRIX_NUMBER>0.AND.MATRIX_NUMBER<=SOLVER_MATRICES%NUMBER_OF_MATRICES) THEN
        SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          IF(ASSOCIATED(SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR)) THEN
            CALL FLAG_ERROR("Solver matrix is already associated.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrix.",ERR,ERROR,*999)
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR
            SOLVER_MATRIX%MATRIX_NUMBER=MATRIX_NUMBER
            SOLVER_MATRIX%SOLVER_MATRICES=>SOLVER_MATRICES
            SOLVER_MATRIX%STORAGE_TYPE=MATRIX_BLOCK_STORAGE_TYPE
            SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
            SOLVER_MATRIX%NUMBER_OF_COLUMNS=SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(MATRIX_NUMBER)%NUMBER_OF_COLUMNS
            SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(MATRIX_NUMBER)%SOLVER_MATRIX=>SOLVER_MATRIX
            NULLIFY(SOLVER_MATRIX%SOLVER_VECTOR)
            NULLIFY(SOLVER_MATRIX%MATRIX)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified matrix number of "//TRIM(NUMBER_TO_VSTRING(MATRIX_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The number must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_MATRIX_INITIALISE")
    RETURN
999 CALL SOLVER_MATRIX_FINALISE(SOLVER_MATRICES%MATRICES(MATRIX_NUMBER)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRIX_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRIX_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_INITIALISE
        
  !
  !================================================================================================================================
  !

END MODULE SOLVER_MATRICES_ROUTINES
