!> \file
!> $Id: solver_matrices_routines.f90 28 2007-07-27 08:35:14Z cpb $
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

!> This module handles all solver matrix and rhs routines.
MODULE SOLVER_MATRICES_ROUTINES

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_MATRICES_CREATE_FINISH,SOLVER_MATRICES_CREATE_START,SOLVER_MATRICES_DESTROY,SOLVER_MATRICES_LIBRARY_TYPE_SET, &
    & SOLVER_MATRICES_OUTPUT

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
    INTEGER(INTG) :: matrix_idx,NUMBER_OF_NON_ZEROS
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER

    NULLIFY(COLUMN_INDICES)
    NULLIFY(ROW_INDICES)
    
    CALL ENTERS("SOLVER_MATRICES_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Solver matrices have already been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLVER_MATRICES%CREATE_VALUES_CACHE)) THEN
          SOLVER=>SOLVER_MATRICES%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            SOLUTION_MAPPING=>SOLVER%SOLUTION_MAPPING
            IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
              !Allocate the solver matrices
              ALLOCATE(SOLVER_MATRICES%MATRICES(SOLVER_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrices matrices",ERR,ERROR,*999)
              DO matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
                !Create the global matrix
                CALL SOLVER_MATRIX_INITIALISE(SOLVER_MATRICES%MATRICES(matrix_idx),ERR,ERROR,*999)
                SOLVER_MATRICES%MATRICES(matrix_idx)%SOLVER_MATRICES=>SOLVER_MATRICES
                SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx)%SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(matrix_idx)
                SOLVER_MATRICES%MATRICES(matrix_idx)%SOLVER_MATRIX_MAP=>SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx)
                !Set up the matrix parameters
                SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX_NUMBER=matrix_idx
                SOLVER_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=SOLVER_MATRICES%CREATE_VALUES_CACHE% &
                    & MATRIX_STORAGE_TYPE(matrix_idx)
                !Calculate the matrix mappings
                SOLVER_MATRICES%MATRICES(matrix_idx)%NUMBER_OF_ROWS=SOLUTION_MAPPING%ROW_DOFS_MAPPING%NUMBER_OF_LOCAL
                SOLVER_MATRICES%MATRICES(matrix_idx)%NUMBER_OF_COLUMNS=SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx)% &
                  & COLUMN_DOFS_MAPPING%NUMBER_OF_GLOBAL
!!TODO: Fix this
                SOLVER_MATRICES%NUMBER_OF_ROWS=SOLUTION_MAPPING%ROW_DOFS_MAPPING%NUMBER_OF_LOCAL
                SOLVER_MATRICES%TOTAL_NUMBER_OF_ROWS=SOLUTION_MAPPING%ROW_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
                !Allocate the distributed solver matrix
                CALL DISTRIBUTED_MATRIX_CREATE_START(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING, &
                  & SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
                CALL DISTRIBUTED_MATRIX_LIBRARY_TYPE_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                  & SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                  & MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                  & SOLVER_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE,ERR,ERROR,*999)                
                !Calculate and set the matrix structure/sparsity pattern
                IF(SOLVER_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE/=MATRIX_BLOCK_STORAGE_TYPE) THEN
                  CALL DISTRIBUTED_MATRIX_SOLVER_STRUCTURE_CALC_FROM_MAP(SOLVER_MATRICES%MATRICES(matrix_idx), &
                    & NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*999)                  
                  CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX,NUMBER_OF_NON_ZEROS, &
                    & ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX,ROW_INDICES, &
                    & COLUMN_INDICES,ERR,ERROR,*999)
                  IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                  IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                ENDIF
                CALL DISTRIBUTED_MATRIX_CREATE_FINISH(SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
                !Allocate the distributed solution vector
                CALL DISTRIBUTED_VECTOR_CREATE_START(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx)%COLUMN_DOFS_MAPPING, &
                  & SOLVER_MATRICES%MATRICES(matrix_idx)%SOLVER_VECTOR,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%SOLVER_VECTOR, &
                  & SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRICES%MATRICES(matrix_idx)%SOLVER_VECTOR, &
                  & MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRICES%MATRICES(matrix_idx)%SOLVER_VECTOR,ERR,ERROR,*999)
              ENDDO !matrix_idx
              !Allocate the distributed rhs vector
              CALL DISTRIBUTED_VECTOR_CREATE_START(SOLUTION_MAPPING%ROW_DOFS_MAPPING,SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_LIBRARY_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR,SOLVER_MATRICES%LIBRARY_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(SOLVER_MATRICES%RHS_VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
              !Finish up
              CALL SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE(SOLVER_MATRICES%CREATE_VALUES_CACHE,ERR,ERROR,*999)
              SOLVER_MATRICES%SOLVER_MATRICES_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("Solver solution mapping is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver matrices solver is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver matrices create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the solver matrices
  SUBROUTINE SOLVER_MATRICES_CREATE_START(SOLVER,SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to create the solver matrices for
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLVER_MATRICES)) THEN
          CALL FLAG_ERROR("Solver matrices is already associated",ERR,ERROR,*998)
        ELSE
          NULLIFY(SOLVER_MATRICES)
          CALL SOLVER_MATRICES_INITIALISE(SOLVER,ERR,ERROR,*999)
          SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_CREATE_START")
    RETURN
999 CALL SOLVER_MATRICES_FINALISE(SOLVER%SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRICES_CREATE_START",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_CREATE_START")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_CREATE_START
        
  !
  !================================================================================================================================
  !

  !>Finalises a solver matrices create value cache and deallocates all memory
  SUBROUTINE SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the the solver matrices create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Intialises a solver matrices create value cache and deallocates all memory
  SUBROUTINE SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE(SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the the solver matrices 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(ASSOCIATED(SOLVER_MATRICES%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Solver matrices create values cache is already associated",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER_MATRICES%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrices create values cache",ERR,ERROR,*999)
        ALLOCATE(SOLVER_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(SOLVER_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrices create value cache matrix storage type",ERR,ERROR,*999)
        SOLVER_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE(SOLVER_MATRICES%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE
        
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
          CALL SOLVER_MATRIX_FINALISE(SOLVER_MATRICES%MATRICES(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(SOLVER_MATRICES%MATRICES)
      ENDIF
      CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
      CALL SOLVER_MATRICES_CREATE_VALUES_CACHE_FINALISE(SOLVER_MATRICES%CREATE_VALUES_CACHE,ERR,ERROR,*999)
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

  !>Initialises the solver matrices for a solver
  SUBROUTINE SOLVER_MATRICES_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the problem solver to initialise the solver matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%SOLVER_MATRICES)) THEN
        CALL FLAG_ERROR("Solver matrices is already associated for this solver",ERR,ERROR,*998)
      ELSE
        SOLUTION_MAPPING=>SOLVER%SOLUTION_MAPPING
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
          ALLOCATE(SOLVER%SOLVER_MATRICES,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrices",ERR,ERROR,*999)
          SOLVER%SOLVER_MATRICES%SOLVER=>SOLVER
          SOLVER%SOLVER_MATRICES%SOLVER_MATRICES_FINISHED=.FALSE.
          SOLVER%SOLVER_MATRICES%SOLUTION_MAPPING=>SOLVER%SOLUTION_MAPPING
          SOLVER%SOLVER_MATRICES%NUMBER_OF_ROWS=SOLUTION_MAPPING%NUMBER_OF_ROWS
          SOLVER%SOLVER_MATRICES%TOTAL_NUMBER_OF_ROWS=SOLUTION_MAPPING%TOTAL_NUMBER_OF_ROWS
          SOLVER%SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
          SOLVER%SOLVER_MATRICES%NUMBER_OF_MATRICES=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
          SOLVER%SOLVER_MATRICES%UPDATE_RHS_VECTOR=.FALSE.
          NULLIFY(SOLVER%SOLVER_MATRICES%RHS_VECTOR)
          NULLIFY(SOLVER%SOLVER_MATRICES%CREATE_VALUES_CACHE)
          CALL SOLVER_MATRICES_CREATE_VALUES_CACHE_INITIALISE(SOLVER%SOLVER_MATRICES,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Solver solution mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_MATRICES_INITIALISE")
    RETURN
999 CALL SOLVER_MATRICES_FINALISE(SOLVER%SOLVER_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_MATRICES_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRICES_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRICES_INITIALISE
        
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
        CALL FLAG_ERROR("Solver matrices is finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(LIBRARY_TYPE)
        CASE(DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE)
          SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_CMISS_TYPE
        CASE(DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE)
          SOLVER_MATRICES%LIBRARY_TYPE=DISTRIBUTED_MATRIX_VECTOR_PETSC_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The library type of "// TRIM(NUMBER_TO_VSTRING(LIBRARY_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_MATRICES_OUTPUT(ID,SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES !<A pointer to the solver matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("SOLVER_MATRICES_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_MATRICES)) THEN
      IF(SOLVER_MATRICES%SOLVER_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"Solver matrices:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Number of matrices = ",SOLVER_MATRICES%NUMBER_OF_MATRICES,ERR,ERROR,*999)
        DO matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
          CALL WRITE_STRING_VALUE(ID,"Solver matrix : ",matrix_idx,ERR,ERROR,*999)
          CALL DISTRIBUTED_MATRIX_OUTPUT(ID,SOLVER_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
        ENDDO !matrix_idx
        IF(ASSOCIATED(SOLVER_MATRICES%RHS_VECTOR)) THEN
          CALL WRITE_STRING(ID,"Solver RHS vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,SOLVER_MATRICES%RHS_VECTOR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver matrices have not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrices is not associated",ERR,ERROR,*999)
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

  !>Finalises a solver matrix and deallocates all memory
  SUBROUTINE SOLVER_MATRIX_FINALISE(SOLVER_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE) :: SOLVER_MATRIX !<The solver matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_MATRIX_FINALISE",ERR,ERROR,*999)

    CALL DISTRIBUTED_MATRIX_DESTROY(SOLVER_MATRIX%MATRIX,ERR,ERROR,*999)
    CALL DISTRIBUTED_VECTOR_DESTROY(SOLVER_MATRIX%SOLVER_VECTOR,ERR,ERROR,*999)
        
    CALL EXITS("SOLVER_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_MATRIX_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRIX_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Initialises a solver matrix
  SUBROUTINE SOLVER_MATRIX_INITIALISE(SOLVER_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_TYPE) :: SOLVER_MATRIX !<The solver matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_MATRIX_INITIALISE",ERR,ERROR,*999)

    SOLVER_MATRIX%MATRIX_NUMBER=0
    NULLIFY(SOLVER_MATRIX%SOLVER_MATRICES)
    NULLIFY(SOLVER_MATRIX%SOLVER_MATRIX_MAP)
    SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
    SOLVER_MATRIX%NUMBER_OF_ROWS=0
    SOLVER_MATRIX%NUMBER_OF_COLUMNS=0
    NULLIFY(SOLVER_MATRIX%MATRIX)
    NULLIFY(SOLVER_MATRIX%SOLVER_VECTOR)
    
    CALL EXITS("SOLVER_MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_MATRIX_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_MATRIX_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_MATRIX_INITIALISE
        
  !
  !================================================================================================================================
  !

END MODULE SOLVER_MATRICES_ROUTINES
