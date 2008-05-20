!> \file
!> $Id: global_matrices_routines.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief This module handles all global matrix and rhs routines.
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

!> This module handles all global matrix and rhs routines.
MODULE GLOBAL_MATRICES_ROUTINES

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

  !> \addtogroup GLOBAL_SOLVER_MATRICES_ROUTINES_GlobalMatrixStructureTypes GLOBAL_SOLVER_MATRICES_ROUTINES::GlobalMatrixStructureTypes
  !> \brief Global matrices structure (sparsity) types
  !> \see GLOBAL_SOLVER_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: GLOBAL_MATRIX_NO_STRUCTURE=1 !<No matrix structure - all elements can contain a value. \see GLOBAL_SOLVER_MATRICES_ROUTINES_GlobalMatrixStructureTypes,GLOBAL_SOLVER_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: GLOBAL_MATRIX_FEM_STRUCTURE=2 !<Finite element matrix structure. \see GLOBAL_SOLVER_MATRICES_ROUTINES_GlobalMatrixStructureTypes,GLOBAL_SOLVER_MATRICES_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC GLOBAL_MATRICES_CREATE_FINISH,GLOBAL_MATRICES_CREATE_START,GLOBAL_MATRICES_DESTROY

  !!TODO check if the elements should be create/destroy rather than initialise/finalise
  PUBLIC GLOBAL_MATRICES_ELEMENT_ADD,GLOBAL_MATRICES_ELEMENT_CALCULATE,GLOBAL_MATRICES_ELEMENT_INITIALISE, &
    & GLOBAL_MATRICES_ELEMENT_FINALISE,GLOBAL_MATRICES_VALUES_INITIALISE

  PUBLIC GLOBAL_MATRICES_OUTPUT

  PUBLIC GLOBAL_MATRICES_COEFFICIENTS_SET,GLOBAL_MATRICES_NUMBER_SET,GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET, &
    & GLOBAL_MATRICES_STORAGE_TYPE_SET,GLOBAL_MATRICES_STRUCTURE_TYPE_SET,GLOBAL_MATRICES_VARIABLE_TYPES_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the creation of the global matrices and rhs for the the problem solution
  SUBROUTINE GLOBAL_MATRICES_CREATE_FINISH(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,NUMBER_OF_NON_ZEROS,variable_type
    INTEGER(INTG), POINTER :: ROW_INDICES(:),COLUMN_INDICES(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_VARIABLE_DOMAIN_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(ROW_INDICES)
    NULLIFY(COLUMN_INDICES)

    CALL ENTERS("GLOBAL_MATRICES_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices have already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
          PROBLEM_SOLUTION=>GLOBAL_MATRICES%PROBLEM_SOLUTION
          IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
            PROBLEM=>PROBLEM_SOLUTION%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                !Calculate the number of variable type maps and initialise
                DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
                  variable_type=GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)
                  GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES= &
                    & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES+1
                ENDDO !matrix_idx
                DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                  ALLOCATE(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(GLOBAL_MATRICES% &
                    & VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable type maps matrix map",ERR,ERROR,*999)
                  DO matrix_idx=1,GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES
                    NULLIFY(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(matrix_idx)%PTR)
                  ENDDO !matrix_idx
                  GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES=0
                ENDDO !variable_type
                !Allocate the global matrices
                ALLOCATE(GLOBAL_MATRICES%MATRICES(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global matrices",ERR,ERROR,*999)
                !Now create the individual global matrices
                DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
                  !Create the global matrix
                  CALL GLOBAL_MATRIX_INITIALISE(GLOBAL_MATRICES%MATRICES(matrix_idx),ERR,ERROR,*999)
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%GLOBAL_MATRICES=>GLOBAL_MATRICES
                  !Set up the matrix parameters
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX_NUMBER=matrix_idx
                  variable_type=GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE_TYPE=variable_type
                  GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES= &
                    & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES+1
                  GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)% &
                    & NUMBER_OF_GLOBAL_MATRICES)%PTR=>GLOBAL_MATRICES%MATRICES(matrix_idx)
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(GLOBAL_MATRICES% &
                    & MATRICES(matrix_idx)%VARIABLE_TYPE)%PTR
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX_COEFFICIENT=GLOBAL_MATRICES%CREATE_VALUES_CACHE% &
                    & MATRIX_COEFFICIENTS(matrix_idx)
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE=GLOBAL_MATRICES%CREATE_VALUES_CACHE% &
                    & MATRIX_STORAGE_TYPE(matrix_idx)
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE=GLOBAL_MATRICES%CREATE_VALUES_CACHE% &
                    & MATRIX_STRUCTURE_TYPE(matrix_idx)
                  !Calculate the matrix mappings
                  DEPENDENT_VARIABLE_DOMAIN_MAP=>GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE%DOMAIN_MAPPING
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%NUMBER_OF_ROWS=DEPENDENT_VARIABLE_DOMAIN_MAP%TOTAL_NUMBER_OF_LOCAL
                  GLOBAL_MATRICES%MATRICES(matrix_idx)%NUMBER_OF_COLUMNS=DEPENDENT_VARIABLE_DOMAIN_MAP%NUMBER_OF_GLOBAL
!!TODO: Fix this
                  GLOBAL_MATRICES%NUMBER_OF_ROWS=DEPENDENT_VARIABLE_DOMAIN_MAP%TOTAL_NUMBER_OF_LOCAL
                  GLOBAL_MATRICES%TOTAL_NUMBER_OF_ROWS=DEPENDENT_VARIABLE_DOMAIN_MAP%NUMBER_OF_GLOBAL
                  !Create the distributed global matrix
                  CALL DISTRIBUTED_MATRIX_CREATE_START(DEPENDENT_VARIABLE_DOMAIN_MAP,GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                    & ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_DATA_TYPE_SET(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                    & MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_SET(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX, &
                    & GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE,ERR,ERROR,*999)
                  !Calculate and set the matrix structure/sparsity pattern
                  IF(GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE/=MATRIX_BLOCK_STORAGE_TYPE) THEN
                    CALL GLOBAL_MATRIX_STRUCTURE_CALCULATE(GLOBAL_MATRICES,matrix_idx,NUMBER_OF_NON_ZEROS,ROW_INDICES, &
                      & COLUMN_INDICES,ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_NUMBER_NON_ZEROS_SET(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX,NUMBER_OF_NON_ZEROS, &
                      & ERR,ERROR,*999)
                    CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_SET(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX,ROW_INDICES, &
                      & COLUMN_INDICES,ERR,ERROR,*999)
                    IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
                    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
                  ENDIF
                  CALL DISTRIBUTED_MATRIX_CREATE_FINISH(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
                ENDDO !matrix_idx
                !Finish setting up the global RHS vector
                IF(GLOBAL_MATRICES%RHS_VARIABLE_TYPE==0) THEN
                  NULLIFY(GLOBAL_MATRICES%RHS_VARIABLE)
                ELSE
                  GLOBAL_MATRICES%RHS_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(GLOBAL_MATRICES%RHS_VARIABLE_TYPE)%PTR
                  DEPENDENT_VARIABLE_DOMAIN_MAP=>GLOBAL_MATRICES%RHS_VARIABLE%DOMAIN_MAPPING
                  CALL DISTRIBUTED_VECTOR_CREATE_START(DEPENDENT_VARIABLE_DOMAIN_MAP,GLOBAL_MATRICES%VECTOR,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(GLOBAL_MATRICES%VECTOR,MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_CREATE_FINISH(GLOBAL_MATRICES%VECTOR,ERR,ERROR,*999)
                ENDIF
                !Finish the create values cache
                CALL GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE(GLOBAL_MATRICES%CREATE_VALUES_CACHE,ERR,ERROR,*999)
                GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED=.TRUE.
              ELSE
                CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Problem solution problem is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Global matrices problem solution is not associated",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices create values cache is not associated",ERR,ERROR,*998)
        ENDIF        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(ROW_INDICES)) DEALLOCATE(ROW_INDICES)
    IF(ASSOCIATED(COLUMN_INDICES)) DEALLOCATE(COLUMN_INDICES)
    CALL GLOBAL_MATRICES_FINALISE(GLOBAL_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GLOBAL_MATRICES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of the global matrices and rhs for the the problem solution
  SUBROUTINE GLOBAL_MATRICES_CREATE_START(PROBLEM_SOLUTION,GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<The pointer to the problem solution
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<On return, a pointer to the global matrices being created.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string  
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR    

    CALL ENTERS("GLOBAL_MATRICES_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(PROBLEM_SOLUTION%SOLUTION_FINISHED) THEN
        CALL FLAG_ERROR("Problem solution has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
          CALL FLAG_ERROR("Global matrices is already associated",ERR,ERROR,*998)
        ELSE
          NULLIFY(GLOBAL_MATRICES)
          IF(ASSOCIATED(PROBLEM_SOLUTION%PROBLEM)) THEN
            !Initialise the global matrices
            CALL GLOBAL_MATRICES_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*999)
            GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
          ELSE
            CALL FLAG_ERROR("Problem solution problem is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_CREATE_START")
    RETURN
999 CALL GLOBAL_MATRICES_FINALISE(PROBLEM_SOLUTION%GLOBAL_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GLOBAL_MATRICES_CREATE_START",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_CREATE_START")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the global matrices
  SUBROUTINE GLOBAL_MATRICES_DESTROY(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer the global matrices to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GLOBAL_MATRICES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      CALL GLOBAL_MATRICES_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("GLOBAL_MATRICES_DESTROY")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_DESTROY",ERR,ERROR)    
    CALL EXITS("GLOBAL_MATRICES_DESTROY")
    RETURN 1
   
  END SUBROUTINE GLOBAL_MATRICES_DESTROY

  !
  !================================================================================================================================
  !
  
  !>Finalises the DOF to global matrices map for a problem and deallocates all memory
  SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE(DOF_TO_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(DOF_TO_GLOBAL_MAP_TYPE) :: DOF_TO_GLOBAL_MAP !<The DOF to global map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(DOF_TO_GLOBAL_MAP%MATRIX_MAPPINGS)) DEALLOCATE(DOF_TO_GLOBAL_MAP%MATRIX_MAPPINGS)
    IF(ALLOCATED(DOF_TO_GLOBAL_MAP%COLUMN_MAPPINGS)) DEALLOCATE(DOF_TO_GLOBAL_MAP%COLUMN_MAPPINGS)
      
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE

  !
  !================================================================================================================================
  !
  
  !>Initialises the DOF to global matrices map.
  SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE(DOF_TO_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(DOF_TO_GLOBAL_MAP_TYPE) :: DOF_TO_GLOBAL_MAP !<The DOF to global map to initalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE",ERR,ERROR,*999)

    DOF_TO_GLOBAL_MAP%ROW_MAPPING=0
    DOF_TO_GLOBAL_MAP%NUMBER_OF_MATRICES=0
      
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Finalises the DOF to global matrices map for a problem and deallocates all memory
  SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD

    CALL ENTERS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(ALLOCATED(GLOBAL_MATRICES%DOF_TO_GLOBAL_MAPS)) THEN
        DEPENDENT_FIELD=>GLOBAL_MATRICES%PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          DO ny=1,DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING%NUMBER_OF_LOCAL
            CALL GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_FINALISE(GLOBAL_MATRICES%DOF_TO_GLOBAL_MAPS(ny),ERR,ERROR,*999)
          ENDDO !ny
        ENDIF
        DEALLOCATE(GLOBAL_MATRICES%DOF_TO_GLOBAL_MAPS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the DOF to global matrices map information for the global matrices
  SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,ny
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(ALLOCATED(GLOBAL_MATRICES%DOF_TO_GLOBAL_MAPS)) THEN
        CALL FLAG_ERROR("DOF to global mapping is already associated for the global matrices",ERR,ERROR,*998)
      ELSE
        PROBLEM_SOLUTION=>GLOBAL_MATRICES%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          PROBLEM=>PROBLEM_SOLUTION%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            DEPENDENT_FIELD=>GLOBAL_MATRICES%PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              ALLOCATE(GLOBAL_MATRICES%DOF_TO_GLOBAL_MAPS(DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING%NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global matrices dof to global maps",ERR,ERROR,*999)
              DO ny=1,DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING%NUMBER_OF_LOCAL
                CALL GLOBAL_MATRICES_DOF_TO_GLOBAL_MAP_INITIALISE(GLOBAL_MATRICES%DOF_TO_GLOBAL_MAPS(ny),ERR,ERROR,*999)
              ENDDO !ny
            ELSE
              CALL FLAG_ERROR("Problem dependent field is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Problem solution problem is not associated",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices problem solution is not associated",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE")
    RETURN
999 CALL GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE(GLOBAL_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)    
998 CALL ERRORS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise an element matrix and deallocate all memory
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE):: ELEMENT_MATRIX !<The element matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT_MATRIX%ROW_DOFS)) DEALLOCATE(ELEMENT_MATRIX%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%COLUMN_DOFS)) DEALLOCATE(ELEMENT_MATRIX%COLUMN_DOFS)
    IF(ALLOCATED(ELEMENT_MATRIX%MATRIX)) DEALLOCATE(ELEMENT_MATRIX%MATRIX)
    
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element matrix.
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE(ELEMENT_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_MATRIX_TYPE) :: ELEMENT_MATRIX !The element matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE",ERR,ERROR,*999)

    ELEMENT_MATRIX%GLOBAL_MATRIX_NUMBER=0
    ELEMENT_MATRIX%NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
    ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
       
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise an element vector and deallocate all memory
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE):: ELEMENT_VECTOR !<The element vector to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(ELEMENT_VECTOR%ROW_DOFS)
    IF(ALLOCATED(ELEMENT_VECTOR%VECTOR)) DEALLOCATE(ELEMENT_VECTOR%VECTOR)
    
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element vector
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE(ELEMENT_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR !The element vector to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE",ERR,ERROR,*999)

    ELEMENT_VECTOR%NUMBER_OF_ROWS=0
    ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
       
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Adds the element matrices and rhs vector into the global matrices and rhs vector.
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_ADD(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      !Add the element matrices
      DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
        IF(GLOBAL_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX) THEN
          CALL DISTRIBUTED_MATRIX_VALUES_ADD(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX,GLOBAL_MATRICES%MATRICES(matrix_idx)% &
            & ELEMENT_MATRIX%ROW_DOFS(1:GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS),GLOBAL_MATRICES% &
            & MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS(1:GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX% &
            & NUMBER_OF_ROWS),GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX(1:GLOBAL_MATRICES%MATRICES(matrix_idx)% &
            & ELEMENT_MATRIX%NUMBER_OF_ROWS,1:GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS), &
            & ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      !Add the rhs
      IF(GLOBAL_MATRICES%UPDATE_VECTOR) THEN
        CALL DISTRIBUTED_VECTOR_VALUES_ADD(GLOBAL_MATRICES%VECTOR,GLOBAL_MATRICES%ELEMENT_VECTOR%ROW_DOFS(1:GLOBAL_MATRICES% &
          & ELEMENT_VECTOR%NUMBER_OF_ROWS),GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR(1:GLOBAL_MATRICES%ELEMENT_VECTOR% &
          & NUMBER_OF_ROWS),ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_ADD")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_ADD",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_ADD")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_ADD

  !
  !================================================================================================================================
  !

  !>Calculate the positions in the global matrices and rhs of the element matrices and rhs vector. Old CMISS name MELGE.
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_CALCULATE(GLOBAL_MATRICES,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative,derivative_idx,global_ny,local_ny,matrix_idx,node,node_idx
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      !Calculate the row and columns for the global matrices
      DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
        GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS=0
        GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS=0
        IF(GLOBAL_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX) THEN
          FIELD_VARIABLE=>GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE
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
                    GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS= &
                      & GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_ROWS+1
                    GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS= &
                      & GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%NUMBER_OF_COLUMNS+1
                    GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS(GLOBAL_MATRICES%MATRICES(matrix_idx)% &
                      & ELEMENT_MATRIX%NUMBER_OF_ROWS)=local_ny
                    GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS(GLOBAL_MATRICES%MATRICES(matrix_idx)% &
                      & ELEMENT_MATRIX%NUMBER_OF_COLUMNS)=global_ny
                  ENDDO !derivative_idx
                ENDDO !node_idx
              ELSE
                LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                  & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                  & ". The element number must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The interpolation type for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                & " is not nodally based"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
            ENDIF
          ENDDO !component_idx
          GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX=0.0_DP
        ENDIF
      ENDDO !matrix_idx
      !Calculate the row and columns for the global RHS
      GLOBAL_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS=0
      IF(GLOBAL_MATRICES%UPDATE_VECTOR) THEN
        FIELD_VARIABLE=>GLOBAL_MATRICES%RHS_VARIABLE
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
                  GLOBAL_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS=GLOBAL_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS+1
                  GLOBAL_MATRICES%ELEMENT_VECTOR%ROW_DOFS(GLOBAL_MATRICES%ELEMENT_VECTOR%NUMBER_OF_ROWS)=local_ny
                ENDDO !derivative_idx
              ENDDO !node_idx
            ELSE
              LOCAL_ERROR="Element number "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
                & ". The element number must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The interpolation type for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " of dependent variable number "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%VARIABLE_NUMBER,"*",ERR,ERROR))// &
              & " is not nodally based"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
          ENDIF
        ENDDO !component_idx
        GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR=0.0_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalise the element calculation information and deallocate all memory
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The global matrices for which to finalise the elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      !Finalise the element matrices
      DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
        GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS=0
        GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS=0
        IF(ALLOCATED(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS)) &
          & DEALLOCATE(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS)
        IF(ALLOCATED(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS)) &
          & DEALLOCATE(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS)
        IF(ALLOCATED(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX)) &
          & DEALLOCATE(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX)        
      ENDDO !matrix_idx
      !Finalise the element vector
      GLOBAL_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
      IF(ALLOCATED(GLOBAL_MATRICES%ELEMENT_VECTOR%ROW_DOFS)) DEALLOCATE(GLOBAL_MATRICES%ELEMENT_VECTOR%ROW_DOFS)
      IF(ALLOCATED(GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR)) DEALLOCATE(GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR)
    ELSE
      CALL FLAG_ERROR("Global matrices is not allocated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the element calculation information for the global matrices
  SUBROUTINE GLOBAL_MATRICES_ELEMENT_INITIALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !The global matrices to initialise the element information for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("GLOBAL_MATRICES_ELEMENT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      !Finalise the element matrices
      DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
        FIELD_VARIABLE=>GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE
        GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS= &
          & FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
        GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS= &
          & FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
        IF(ALLOCATED(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS)) THEN
          CALL FLAG_ERROR("Element matrix row dofs already allocated",ERR,ERROR,*999)
        ELSE
          ALLOCATE(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%ROW_DOFS(GLOBAL_MATRICES%MATRICES(matrix_idx)% &
            & ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix row dofs",ERR,ERROR,*999)
        ENDIF
        IF(ALLOCATED(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS)) THEN
          CALL FLAG_ERROR("Element matrix column dofs already allocated",ERR,ERROR,*999)
        ELSE
          ALLOCATE(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%COLUMN_DOFS(GLOBAL_MATRICES%MATRICES(matrix_idx)% &
            & ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix column dofs",ERR,ERROR,*999)
       ENDIF
       IF(ALLOCATED(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX)) THEN
         CALL FLAG_ERROR("Element matrix already allocated",ERR,ERROR,*999)
       ELSE
         ALLOCATE(GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MATRIX( &
           & GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
           & GLOBAL_MATRICES%MATRICES(matrix_idx)%ELEMENT_MATRIX%MAX_NUMBER_OF_COLUMNS),STAT=ERR)
         IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element matrix",ERR,ERROR,*999)
       ENDIF
      ENDDO !matrix_idx
      !Initialise the element vector
      FIELD_VARIABLE=>GLOBAL_MATRICES%RHS_VARIABLE
      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
        GLOBAL_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
        IF(ALLOCATED(GLOBAL_MATRICES%ELEMENT_VECTOR%ROW_DOFS)) THEN
          CALL FLAG_ERROR("Element vector row dofs already allocated",ERR,ERROR,*999)        
        ELSE
          ALLOCATE(GLOBAL_MATRICES%ELEMENT_VECTOR%ROW_DOFS(GLOBAL_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector row dofs",ERR,ERROR,*999)
        ENDIF
        IF(ALLOCATED(GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR)) THEN
          CALL FLAG_ERROR("Element vector already allocated",ERR,ERROR,*999)        
        ELSE
          ALLOCATE(GLOBAL_MATRICES%ELEMENT_VECTOR%VECTOR(GLOBAL_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element vector",ERR,ERROR,*999)
        ENDIF
      ELSE
        GLOBAL_MATRICES%ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not allocated",ERR,ERROR,*998)
    ENDIF    
    
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_INITIALISE")
    RETURN
999 CALL GLOBAL_MATRICES_ELEMENT_FINALISE(GLOBAL_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GLOBAL_MATRICES_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise a global matrix and deallocate all memory
  SUBROUTINE GLOBAL_MATRIX_FINALISE(GLOBAL_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRIX_TYPE):: GLOBAL_MATRIX !<The global matrix to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("GLOBAL_MATRIX_FINALISE",ERR,ERROR,*999)

    CALL DISTRIBUTED_MATRIX_DESTROY(GLOBAL_MATRIX%MATRIX,ERR,ERROR,*999)
    CALL GLOBAL_MATRICES_ELEMENT_MATRIX_FINALISE(GLOBAL_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
    
    CALL EXITS("GLOBAL_MATRIX_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRIX_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRIX_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRIX_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the global matrix.
  SUBROUTINE GLOBAL_MATRIX_INITIALISE(GLOBAL_MATRIX,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRIX_TYPE) :: GLOBAL_MATRIX !<The global matrix to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("GLOBAL_MATRIX_INITIALISE",ERR,ERROR,*999)

    GLOBAL_MATRIX%MATRIX_NUMBER=0
    NULLIFY(GLOBAL_MATRIX%GLOBAL_MATRICES)
    GLOBAL_MATRIX%VARIABLE_TYPE=0
    NULLIFY(GLOBAL_MATRIX%VARIABLE)
    GLOBAL_MATRIX%MATRIX_COEFFICIENT=1.0_DP !Matrices in an equation set are added by default
    GLOBAL_MATRIX%STORAGE_TYPE=0
    GLOBAL_MATRIX%STRUCTURE_TYPE=0
    GLOBAL_MATRIX%UPDATE_MATRIX=.TRUE.
    GLOBAL_MATRIX%NUMBER_OF_ROWS=0
    GLOBAL_MATRIX%NUMBER_OF_COLUMNS=0
    NULLIFY(GLOBAL_MATRIX%MATRIX)
    CALL GLOBAL_MATRICES_ELEMENT_MATRIX_INITIALISE(GLOBAL_MATRIX%ELEMENT_MATRIX,ERR,ERROR,*999)
       
    CALL EXITS("GLOBAL_MATRIX_INITIALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRIX_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRIX_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRIX_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the coefficients for the global matrices in an equation set
  SUBROUTINE GLOBAL_MATRICES_COEFFICIENTS_SET(GLOBAL_MATRICES,MATRIX_COEFFICIENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices.
    REAL(DP), INTENT(IN) :: MATRIX_COEFFICIENTS(:) !<The matrix coefficients
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GLOBAL_MATRICES_COEFFICIENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices is finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
          IF(SIZE(MATRIX_COEFFICIENTS,1)==GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN            
            GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS=MATRIX_COEFFICIENTS
          ELSE
            LOCAL_ERROR="Invalid size of matrix coefficeints. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_COEFFICIENTS,1),"*",ERR,ERROR))// &
              & ") must match the number of global matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_COEFFICIENTS_SET")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_COEFFICIENTS_SET",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_COEFFICIENTS_SET")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_COEFFICIENTS_SET

  !
  !================================================================================================================================
  !

  !>Finalises the create values cache for the global matrices
  SUBROUTINE GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string   
    !Local Variables

    CALL ENTERS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the create values cache for the global matrices
  SUBROUTINE GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The pointer to the global matrices.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string 
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Global matrices create values cache is already associated",ERR,ERROR,*998)
      ELSE
        DEPENDENT_FIELD=>GLOBAL_MATRICES%PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
        !Create as many global matrices and rhs vector as there are dependent field variables. By default the second dependent
        !field variable will be the global rhs variable. 
        IF(DEPENDENT_FIELD%NUMBER_OF_VARIABLES==1) THEN
          CALL FLAG_ERROR("Dependent field only has one variable which cannot be mapped to both a global matrix and rhs vector", &
            & ERR,ERROR,*999)
        ELSE IF(DEPENDENT_FIELD%NUMBER_OF_VARIABLES==2) THEN
          GLOBAL_MATRICES%NUMBER_OF_MATRICES=1
          GLOBAL_MATRICES%RHS_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(2)%PTR%VARIABLE_TYPE
        ELSE
          GLOBAL_MATRICES%NUMBER_OF_MATRICES=DEPENDENT_FIELD%NUMBER_OF_VARIABLES-1
          GLOBAL_MATRICES%RHS_VARIABLE_TYPE=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(2)%PTR%VARIABLE_TYPE
        ENDIF
        ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global matrices create values cache",ERR,ERROR,*999)
        ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix variable types",ERR,ERROR,*999)
        ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix storage type",ERR,ERROR,*999)
        ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix structure type",ERR,ERROR,*999)                
        ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache matrix coefficients",ERR,ERROR,*999)                
        GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1)=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(1)%PTR%VARIABLE_TYPE
        GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(1)=MATRIX_BLOCK_STORAGE_TYPE
        GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(1)=GLOBAL_MATRIX_NO_STRUCTURE
        GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1)=1.0_DP
        DO matrix_idx=2,GLOBAL_MATRICES%NUMBER_OF_MATRICES
          GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)= &
            & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(matrix_idx+1)%PTR%VARIABLE_TYPE
          GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_BLOCK_STORAGE_TYPE
          GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(matrix_idx)=GLOBAL_MATRIX_NO_STRUCTURE
          GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(matrix_idx)=1.0_DP
        ENDDO !matrix_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE(GLOBAL_MATRICES%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variables and the global matrices
  SUBROUTINE GLOBAL_MATRICES_NUMBER_SET(GLOBAL_MATRICES,NUMBER_OF_GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_GLOBAL_MATRICES !<The number of global matrices for the problem.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_VARIABLE_TYPES(:),OLD_MATRIX_STORAGE_TYPE(:),OLD_MATRIX_STRUCTURE_TYPE(:)
    REAL(DP), ALLOCATABLE :: OLD_MATRIX_COEFFICIENTS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GLOBAL_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices have been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
          !Check number of matrices to create is valid
          IF(GLOBAL_MATRICES%RHS_VARIABLE_TYPE==0) THEN
            IF(NUMBER_OF_GLOBAL_MATRICES<1.OR.NUMBER_OF_GLOBAL_MATRICES>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
              LOCAL_ERROR="The requested number of matrices to create ("// &
                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GLOBAL_MATRICES,"*",ERR,ERROR))// &
                & ") is invalid. For problems without a global RHS the number must be between >= 1 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE            
            IF(NUMBER_OF_GLOBAL_MATRICES<1.OR.NUMBER_OF_GLOBAL_MATRICES>=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
              LOCAL_ERROR="The requested number of matrices to create ("// &
                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GLOBAL_MATRICES,"*",ERR,ERROR))// &
                & ") is invalid. For problems with a global RHS the number must be between >= 1 and < "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
          !If we need to reallocate and reset all the create_values cache arrays and change the number of matrices
          IF(NUMBER_OF_GLOBAL_MATRICES/=GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
            ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types",ERR,ERROR,*999)
            ALLOCATE(OLD_MATRIX_STORAGE_TYPE(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix storage type",ERR,ERROR,*999)
            ALLOCATE(OLD_MATRIX_STRUCTURE_TYPE(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix structure type",ERR,ERROR,*999)
            ALLOCATE(OLD_MATRIX_COEFFICIENTS(GLOBAL_MATRICES%NUMBER_OF_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix coefficients",ERR,ERROR,*999)
            OLD_MATRIX_VARIABLE_TYPES=GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
            OLD_MATRIX_STORAGE_TYPE=GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE
            OLD_MATRIX_STRUCTURE_TYPE=GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE
            OLD_MATRIX_COEFFICIENTS=GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS
            DEALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
            DEALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE)
            DEALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE)
            DEALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
            ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types",ERR,ERROR,*999)
            ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix storage type",ERR,ERROR,*999)
            ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix structure type",ERR,ERROR,*999)
            ALLOCATE(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix coefficients",ERR,ERROR,*999)
            IF(NUMBER_OF_GLOBAL_MATRICES>GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:GLOBAL_MATRICES%NUMBER_OF_MATRICES)= &
                & OLD_MATRIX_VARIABLE_TYPES
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(GLOBAL_MATRICES%NUMBER_OF_MATRICES+1: &
                & NUMBER_OF_GLOBAL_MATRICES)=OLD_MATRIX_VARIABLE_TYPES(1)
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(1:GLOBAL_MATRICES%NUMBER_OF_MATRICES)= &
                & OLD_MATRIX_STORAGE_TYPE
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(GLOBAL_MATRICES%NUMBER_OF_MATRICES+1: &
                & NUMBER_OF_GLOBAL_MATRICES)=OLD_MATRIX_STORAGE_TYPE(1)
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(1:GLOBAL_MATRICES%NUMBER_OF_MATRICES)= &
                & OLD_MATRIX_STRUCTURE_TYPE
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(GLOBAL_MATRICES%NUMBER_OF_MATRICES+1: &
                & NUMBER_OF_GLOBAL_MATRICES)=OLD_MATRIX_STRUCTURE_TYPE(1)
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:GLOBAL_MATRICES%NUMBER_OF_MATRICES)= &
                & OLD_MATRIX_COEFFICIENTS
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(GLOBAL_MATRICES%NUMBER_OF_MATRICES+1: &
                & NUMBER_OF_GLOBAL_MATRICES)=OLD_MATRIX_COEFFICIENTS(1)
            ELSE
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:NUMBER_OF_GLOBAL_MATRICES)= &
                & OLD_MATRIX_VARIABLE_TYPES(1:NUMBER_OF_GLOBAL_MATRICES)
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(1:NUMBER_OF_GLOBAL_MATRICES)= &
                & OLD_MATRIX_STORAGE_TYPE(1:NUMBER_OF_GLOBAL_MATRICES)
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(1:NUMBER_OF_GLOBAL_MATRICES)= &
                & OLD_MATRIX_STRUCTURE_TYPE(1:NUMBER_OF_GLOBAL_MATRICES)
              GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS(1:NUMBER_OF_GLOBAL_MATRICES)= &
                & OLD_MATRIX_COEFFICIENTS(1:NUMBER_OF_GLOBAL_MATRICES)
            ENDIF
            GLOBAL_MATRICES%NUMBER_OF_MATRICES=NUMBER_OF_GLOBAL_MATRICES
            IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
            IF(ALLOCATED(OLD_MATRIX_STORAGE_TYPE)) DEALLOCATE(OLD_MATRIX_STORAGE_TYPE)
            IF(ALLOCATED(OLD_MATRIX_STRUCTURE_TYPE)) DEALLOCATE(OLD_MATRIX_STRUCTURE_TYPE)
            IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
    IF(ALLOCATED(OLD_MATRIX_STORAGE_TYPE)) DEALLOCATE(OLD_MATRIX_STORAGE_TYPE)
    IF(ALLOCATED(OLD_MATRIX_STRUCTURE_TYPE)) DEALLOCATE(OLD_MATRIX_STRUCTURE_TYPE)
    IF(ALLOCATED(OLD_MATRIX_COEFFICIENTS)) DEALLOCATE(OLD_MATRIX_COEFFICIENTS)
    CALL ERRORS("GLOBAL_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Outputs the global matrices
  SUBROUTINE GLOBAL_MATRICES_OUTPUT(ID,GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the ouptut stream
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("GLOBAL_MATRICES_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL WRITE_STRING(ID,"Global matrices:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Number of matrices = ",GLOBAL_MATRICES%NUMBER_OF_MATRICES,ERR,ERROR,*999)
        DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
          CALL WRITE_STRING_VALUE(ID,"Global matrix : ",matrix_idx,ERR,ERROR,*999)
          CALL DISTRIBUTED_MATRIX_OUTPUT(ID,GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX,ERR,ERROR,*999)
        ENDDO !matrix_idx
        IF(GLOBAL_MATRICES%RHS_VARIABLE_TYPE/=0) THEN
          CALL WRITE_STRING(ID,"Global RHS vector:",ERR,ERROR,*999)
          CALL DISTRIBUTED_VECTOR_OUTPUT(ID,GLOBAL_MATRICES%VECTOR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Global matrices have not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_OUTPUT")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_OUTPUT",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_OUTPUT")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_OUTPUT

  !
  !================================================================================================================================
  !

  !>Sets the storage type (sparsity) of the global matrices
  SUBROUTINE GLOBAL_MATRICES_STORAGE_TYPE_SET(GLOBAL_MATRICES,STORAGE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(IN) :: STORAGE_TYPE(:) !<STORAGE_TYPE(matrix_idx). The storage type for the matrix_idx'th global matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GLOBAL_MATRICES_STORAGE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices have been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
          IF(SIZE(STORAGE_TYPE,1)==GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
            DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
              SELECT CASE(STORAGE_TYPE(matrix_idx))
              CASE(MATRIX_BLOCK_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_BLOCK_STORAGE_TYPE
              CASE(MATRIX_DIAGONAL_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_DIAGONAL_STORAGE_TYPE        
              CASE(MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_COLUMN_MAJOR_STORAGE_TYPE
              CASE(MATRIX_ROW_MAJOR_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_ROW_MAJOR_STORAGE_TYPE
              CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_COMPRESSED_ROW_STORAGE_TYPE
              CASE(MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE
              CASE(MATRIX_ROW_COLUMN_STORAGE_TYPE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STORAGE_TYPE(matrix_idx)=MATRIX_ROW_COLUMN_STORAGE_TYPE
              CASE DEFAULT
                LOCAL_ERROR="The specified storage type of "//TRIM(NUMBER_TO_VSTRING(STORAGE_TYPE(matrix_idx),"*",ERR,ERROR))// &
                  & " for the matrix number "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the storage type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STORAGE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRICES_STORAGE_TYPE_SET")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_STORAGE_TYPE_SET",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_STORAGE_TYPE_SET")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_STORAGE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets the structure (sparsity) of the global matrices
  SUBROUTINE GLOBAL_MATRICES_STRUCTURE_TYPE_SET(GLOBAL_MATRICES,STRUCTURE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(IN) :: STRUCTURE_TYPE(:) !<STRUCTURE_TYPE(matrix_idx). The storage type for the matrix_idx'th global matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GLOBAL_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices have been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
          IF(SIZE(STRUCTURE_TYPE,1)==GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
            DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
              SELECT CASE(STRUCTURE_TYPE(matrix_idx))
              CASE(GLOBAL_MATRIX_NO_STRUCTURE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(matrix_idx)=GLOBAL_MATRIX_NO_STRUCTURE
              CASE(GLOBAL_MATRIX_FEM_STRUCTURE)
                GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_STRUCTURE_TYPE(matrix_idx)=GLOBAL_MATRIX_FEM_STRUCTURE
              CASE DEFAULT
                LOCAL_ERROR="The specified strucutre type of "// &
                  & TRIM(NUMBER_TO_VSTRING(STRUCTURE_TYPE(matrix_idx),"*",ERR,ERROR))//" for the matrix number "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDDO !matrix_idx
          ELSE
            LOCAL_ERROR="The size of the structure type array ("//TRIM(NUMBER_TO_VSTRING(SIZE(STRUCTURE_TYPE,1),"*",ERR,ERROR))// &
              & ") is not equal to the number of matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices create value cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_STRUCTURE_TYPE_SET")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_STRUCTURE_TYPE_SET",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_STRUCTURE_TYPE_SET")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_STRUCTURE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Finalise the global matrices for the problem solution and deallocate all memory.
  SUBROUTINE GLOBAL_MATRICES_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("GLOBAL_MATRICES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      CALL GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*999)
      CALL GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*999)
      IF(ALLOCATED(GLOBAL_MATRICES%GLOBAL_ROW_TO_DOF_MAP)) DEALLOCATE(GLOBAL_MATRICES%GLOBAL_ROW_TO_DOF_MAP)
      IF(ALLOCATED(GLOBAL_MATRICES%MATRICES)) THEN
        DO matrix_idx=1,SIZE(GLOBAL_MATRICES%MATRICES,1)
          CALL GLOBAL_MATRIX_FINALISE(GLOBAL_MATRICES%MATRICES(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(GLOBAL_MATRICES%MATRICES)
      ENDIF
      CALL DISTRIBUTED_VECTOR_DESTROY(GLOBAL_MATRICES%VECTOR,ERR,ERROR,*999)
      CALL GLOBAL_MATRICES_ELEMENT_VECTOR_FINALISE(GLOBAL_MATRICES%ELEMENT_VECTOR,ERR,ERROR,*999)
      CALL GLOBAL_MATRICES_CREATE_VALUES_CACHE_FINALISE(GLOBAL_MATRICES%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(GLOBAL_MATRICES)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the global matrices for the problem solution.
  SUBROUTINE GLOBAL_MATRICES_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

     !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution to initialise the global matrices for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("GLOBAL_MATRICES_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%GLOBAL_MATRICES)) THEN
        CALL FLAG_ERROR("Global matrices is already associated for this problem solution",ERR,ERROR,*998)
      ELSE
        ALLOCATE(PROBLEM_SOLUTION%GLOBAL_MATRICES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution global matrices",ERR,ERROR,*999)
        PROBLEM_SOLUTION%GLOBAL_MATRICES%PROBLEM_SOLUTION=>PROBLEM_SOLUTION
        PROBLEM_SOLUTION%GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED=.FALSE.
        NULLIFY(PROBLEM_SOLUTION%SOLUTION_MAPPING)
        PROBLEM_SOLUTION%GLOBAL_MATRICES%NUMBER_OF_ROWS=0
        PROBLEM_SOLUTION%GLOBAL_MATRICES%TOTAL_NUMBER_OF_ROWS=0
        PROBLEM_SOLUTION%GLOBAL_MATRICES%NUMBER_OF_MATRICES=0
        PROBLEM_SOLUTION%GLOBAL_MATRICES%RHS_VARIABLE_TYPE=0
        NULLIFY(PROBLEM_SOLUTION%GLOBAL_MATRICES%RHS_VARIABLE)
        CALL GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE(PROBLEM_SOLUTION%GLOBAL_MATRICES,ERR,ERROR,*999)
        CALL GLOBAL_MATRICES_DOF_TO_GLOBAL_MAPS_INITIALISE(PROBLEM_SOLUTION%GLOBAL_MATRICES,ERR,ERROR,*999)
        PROBLEM_SOLUTION%GLOBAL_MATRICES%UPDATE_VECTOR=.TRUE.
        NULLIFY(PROBLEM_SOLUTION%GLOBAL_MATRICES%VECTOR)
        CALL GLOBAL_MATRICES_ELEMENT_VECTOR_INITIALISE(PROBLEM_SOLUTION%GLOBAL_MATRICES%ELEMENT_VECTOR,ERR,ERROR,*999)
        NULLIFY(PROBLEM_SOLUTION%GLOBAL_MATRICES%CREATE_VALUES_CACHE)
        CALL GLOBAL_MATRICES_CREATE_VALUES_CACHE_INITIALISE(PROBLEM_SOLUTION%GLOBAL_MATRICES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_INITIALISE")
    RETURN
999 CALL GLOBAL_MATRICES_FINALISE(PROBLEM_SOLUTION%GLOBAL_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("GLOBAL_MATRICES_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between a dependent field variable and the global rhs vector.
  SUBROUTINE GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET(GLOBAL_MATRICES,RHS_VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices to set
    INTEGER(INTG), INTENT(IN) :: RHS_VARIABLE_TYPE !<The variable type associated with the global rhs vector. If the problem does not have a rhs vector then the variable type on input should be zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices have been finished",ERR,ERROR,*999)
      ELSE
        IF(RHS_VARIABLE_TYPE==0) THEN
          GLOBAL_MATRICES%RHS_VARIABLE_TYPE=0
        ELSE
          IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
            !Check that the given RHS variable number is not already being used for a matrix
            DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
              IF(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)==RHS_VARIABLE_TYPE) THEN
                LOCAL_ERROR="The specified RHS variable type of "//TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " is the same as the matrix variable type for matrix number "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))
              ENDIF
            ENDDO !matrix_idx          
            DEPENDENT_FIELD=>GLOBAL_MATRICES%PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              !Check the RHS variable number is defined on the dependent field
              IF(RHS_VARIABLE_TYPE>=1.AND.RHS_VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                IF(ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(RHS_VARIABLE_TYPE)%PTR)) THEN
                  GLOBAL_MATRICES%RHS_VARIABLE_TYPE=RHS_VARIABLE_TYPE
                ELSE
                  LOCAL_ERROR="The specified RHS variable type of "//TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " is not defined on the dependent field"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The specified RHS variable type of "//TRIM(NUMBER_TO_VSTRING(RHS_VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The number must either be zero or >= 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Global matrices create values cache is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_RHS_VARIABLE_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Initialise the values of the global matrices and rhs to the given value e.g., 0.0_DP
  SUBROUTINE GLOBAL_MATRICES_VALUES_INITIALISE(GLOBAL_MATRICES,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    REAL(DP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    
    CALL ENTERS("GLOBAL_MATRICES_VALUES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
        IF(GLOBAL_MATRICES%MATRICES(matrix_idx)%UPDATE_MATRIX) THEN
          CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(GLOBAL_MATRICES%MATRICES(matrix_idx)%MATRIX,VALUE,ERR,ERROR,*999)
        ENDIF
      ENDDO !matrix_idx
      IF(GLOBAL_MATRICES%UPDATE_VECTOR) THEN
        CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(GLOBAL_MATRICES%VECTOR,VALUE,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_VALUES_INITIALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_VALUES_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_VALUES_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_VALUES_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Finalises the variable types map for global matrices deallocates all memory
  SUBROUTINE GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx

    CALL ENTERS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(ALLOCATED(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS)) THEN
        DO variable_idx=1,SIZE(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS,1)
          !Just handle the finalisation of GLOBAL_MATRICES_VARIABLE_MAP_TYPE here as there is only one item.
          IF(ALLOCATED(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_idx)%MATRIX_MAP)) &
            & DEALLOCATE(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_idx)%MATRIX_MAP)
        ENDDO !variable_idx
        DEALLOCATE(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the variable type map information for the global matrices
  SUBROUTINE GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE(GLOBAL_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<The pointer to the global matrices
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(ALLOCATED(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS)) THEN
        CALL FLAG_ERROR("Variable type maps is already associated for the global matrices",ERR,ERROR,*998)
      ELSE
        ALLOCATE(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global matrices variable type maps",ERR,ERROR,*999)
        DO variable_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
          !Just handle the initialisation of GLOBAL_MATRICES_VARIABLE_MAP_TYPE here as there is only one item.
          GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_idx)%NUMBER_OF_GLOBAL_MATRICES=0
        ENDDO !variable_idx 
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE")
    RETURN
999 CALL GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_FINALISE(GLOBAL_MATRICES,DUMMY_ERR,DUMMY_ERROR,*998)    
998 CALL ERRORS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_VARIABLE_TYPE_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the mapping between the dependent field variable types and the global matrices
  SUBROUTINE GLOBAL_MATRICES_VARIABLE_TYPES_SET(GLOBAL_MATRICES,MATRIX_VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(IN) :: MATRIX_VARIABLE_TYPES(:) !<The matrix variable types to map to each global matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("GLOBAL_MATRICES_VARIABLE_TYPES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        CALL FLAG_ERROR("Global matrices is finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(GLOBAL_MATRICES%CREATE_VALUES_CACHE)) THEN
          IF(SIZE(MATRIX_VARIABLE_TYPES,1)==GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
            DEPENDENT_FIELD=>GLOBAL_MATRICES%PROBLEM_SOLUTION%PROBLEM%DEPENDENT%DEPENDENT_FIELD
            !Check input values
            DO matrix_idx=1,GLOBAL_MATRICES%NUMBER_OF_MATRICES
              !Check to see if the variable numbers are defined on the dependent field
              IF(MATRIX_VARIABLE_TYPES(matrix_idx)>=1.OR.MATRIX_VARIABLE_TYPES(matrix_idx)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                IF(.NOT.ASSOCIATED(DEPENDENT_FIELD%VARIABLE_TYPE_MAP(MATRIX_VARIABLE_TYPES(matrix_idx))%PTR)) THEN
                  LOCAL_ERROR="The matrix variable type of "// &
                    & TRIM(NUMBER_TO_VSTRING(MATRIX_VARIABLE_TYPES(matrix_idx),"*",ERR,ERROR))//" for matrix "// &
                    & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                    & " is not defined on the dependent field"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The matrix variable type of "// &
                  & TRIM(NUMBER_TO_VSTRING(MATRIX_VARIABLE_TYPES(matrix_idx),"*",ERR,ERROR))//" for matrix "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
                  & " is invalid. The variable types must be >= 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              !Check to see if the matrix variable number is not already being used for the RHS vector
              IF(GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(matrix_idx)==GLOBAL_MATRICES%RHS_VARIABLE_TYPE) THEN
                LOCAL_ERROR="The matrix variable type of "// &
                  & TRIM(NUMBER_TO_VSTRING(MATRIX_VARIABLE_TYPES(matrix_idx),"*",ERR,ERROR))//" for matrix  "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))//" is the same as the current RHS global variable type"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
              ENDIF
            ENDDO !matrix_idx
            GLOBAL_MATRICES%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=MATRIX_VARIABLE_TYPES
          ELSE
            LOCAL_ERROR="Invalid size of matrix variable types. The size of the supplied array ("// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(MATRIX_VARIABLE_TYPES,1),"*",ERR,ERROR))// &
              & ") must match the number of global matrices ("// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//")"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("GLOBAL_MATRICES_VARIABLE_TYPES_SET")
    RETURN
999 CALL ERRORS("GLOBAL_MATRICES_VARIABLE_TYPES_SET",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRICES_VARIABLE_TYPES_SET")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRICES_VARIABLE_TYPES_SET

  !
  !================================================================================================================================
  !

  !>Caclulates the matrix structure (sparsity) for a global matrix.
  SUBROUTINE GLOBAL_MATRIX_STRUCTURE_CALCULATE(GLOBAL_MATRICES,matrix_idx,NUMBER_OF_NON_ZEROS,ROW_INDICES,COLUMN_INDICES, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices
    INTEGER(INTG), INTENT(IN) :: matrix_idx !<The matrix number to calculate the structure of
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NON_ZEROS !<On return the number of non-zeros in the matrix
    INTEGER(INTG), POINTER :: ROW_INDICES(:) !<On return a pointer to row location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:) !<On return a pointer to the column location indices in compressed row format. The pointer must be NULL on entry and the calling routine is responsible for deallocation.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::  column_idx,DUMMY_ERR,elem_idx,global_column,local_column,local_ny,mk,mp,ne,nh,nn,nnk,np, &
      & NUMBER_OF_COLUMNS,nyy
    INTEGER(INTG), POINTER :: COLUMNS(:)
    REAL(DP) :: SPARSITY
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_DOMAIN_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES    
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE), POINTER :: DEPENDENT_DOFS_PARAM_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: COLUMN_INDICES_LISTS(:)
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(COLUMNS)
    
    CALL ENTERS("GLOBAL_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR,*998)

    NUMBER_OF_NON_ZEROS=0
    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(matrix_idx>=1.AND.matrix_idx<=GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN            
        IF(.NOT.ASSOCIATED(ROW_INDICES)) THEN
          IF(.NOT.ASSOCIATED(COLUMN_INDICES)) THEN
            SELECT CASE(GLOBAL_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE)
            CASE(GLOBAL_MATRIX_NO_STRUCTURE)
              CALL FLAG_ERROR("Not implemented",ERR,ERROR,*998)
            CASE(GLOBAL_MATRIX_FEM_STRUCTURE)
              SELECT CASE(GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE)
              CASE(MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                PROBLEM_SOLUTION=>GLOBAL_MATRICES%PROBLEM_SOLUTION
                IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
                  PROBLEM=>PROBLEM_SOLUTION%PROBLEM
                  IF(ASSOCIATED(PROBLEM)) THEN
                    DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                      FIELD_VARIABLE=>GLOBAL_MATRICES%MATRICES(matrix_idx)%VARIABLE
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        DEPENDENT_DOFS_DOMAIN_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
                        IF(ASSOCIATED(DEPENDENT_DOFS_DOMAIN_MAPPING)) THEN
                          DEPENDENT_DOFS_PARAM_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOF_TO_PARAM_MAP
                          IF(ASSOCIATED(DEPENDENT_DOFS_PARAM_MAPPING)) THEN
                            !Allocate lists
                            ALLOCATE(COLUMN_INDICES_LISTS(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices lists",ERR,ERROR,*999)
                            !Allocate row indices
                            ALLOCATE(ROW_INDICES(DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL+1),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate row indices",ERR,ERROR,*999)
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
                                CALL LIST_NUMBER_OF_ITEMS_GET(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                                NUMBER_OF_NON_ZEROS=NUMBER_OF_NON_ZEROS+NUMBER_OF_COLUMNS
                                ROW_INDICES(local_ny+1)=NUMBER_OF_NON_ZEROS+1
                              ELSE
                                LOCAL_ERROR="Local dof number "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                                  & " is not a node based dof"
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !local_ny
                            !Allocate and setup the column locations
                            ALLOCATE(COLUMN_INDICES(NUMBER_OF_NON_ZEROS),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate column indices",ERR,ERROR,*999)
                            DO local_ny=1,DEPENDENT_DOFS_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                              CALL LIST_DETACH_AND_DESTROY(COLUMN_INDICES_LISTS(local_ny)%PTR,NUMBER_OF_COLUMNS,COLUMNS, &
                                & ERR,ERROR,*999)
                              DO column_idx=1,NUMBER_OF_COLUMNS
                                COLUMN_INDICES(ROW_INDICES(local_ny)+column_idx-1)=COLUMNS(column_idx)
                              ENDDO !column_idx
                              DEALLOCATE(COLUMNS)
                            ENDDO !local_ny
                            IF(DIAGNOSTICS1) THEN
                              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Global matrix structure:",ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Global matrix number : ",matrix_idx,ERR,ERROR,*999)
                              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",DEPENDENT_DOFS_DOMAIN_MAPPING% &
                                & TOTAL_NUMBER_OF_LOCAL,ERR,ERROR,*999)
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
                            CALL FLAG_ERROR("Dependent dofs parameter mapping is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Dependent dofs domain mapping is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Dependent field variable is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Problem dependent field is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Problem solution problem is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Global matrices problem solution is not associated",ERR,ERROR,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="The matrix storage type of "// &
                  & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%MATRICES(matrix_idx)%STORAGE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The matrix structure type of "// &
                & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%MATRICES(matrix_idx)%STRUCTURE_TYPE,"*",ERR,ERROR))//" is invalid"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Column indices is already associated",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Row indieces is already associated",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The matrix index of "//TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))// &
          & " is invalid. The index must be >= 1 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("GLOBAL_MATRIX_STRUCTURE_CALCULATE")
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
998 CALL ERRORS("GLOBAL_MATRIX_STRUCTURE_CALCULATE",ERR,ERROR)
    CALL EXITS("GLOBAL_MATRIX_STRUCTURE_CALCULATE")
    RETURN 1
  END SUBROUTINE GLOBAL_MATRIX_STRUCTURE_CALCULATE

  !
  !================================================================================================================================
  !
  
END MODULE GLOBAL_MATRICES_ROUTINES
