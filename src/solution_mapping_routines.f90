!> \file
!> $Id: solution_mapping_routines.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief This module handles all solution mapping routines.
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

!> This module handles all solution mapping routines.
MODULE SOLUTION_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE COMP_ENVIRONMENT
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
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

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLUTION_MAPPING_CREATE_FINISH,SOLUTION_MAPPING_CREATE_START,SOLUTION_MAPPING_DESTROY, &
    & SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET,SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solution mapping
  SUBROUTINE SOLUTION_MAPPING_CREATE_FINISH(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mapping has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
          CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE(SOLUTION_MAPPING,ERR,ERROR,*999)
          CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE(SOLUTION_MAPPING,ERR,ERROR,*999)
!!TODO set and check row numbers
          CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
          SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED=.TRUE.            
        ELSE
          CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_FINISH")
    RETURN
999 CALL SOLUTION_MAPPING_FINALISE(SOLUTION_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solution mapping for a problem solution
  SUBROUTINE SOLUTION_MAPPING_CREATE_START(GLOBAL_MATRICES,SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES !<A pointer to the global matrices to create the solution mapping from.
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    
    CALL ENTERS("SOLUTION_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
      IF(GLOBAL_MATRICES%GLOBAL_MATRICES_FINISHED) THEN
        PROBLEM_SOLUTION=>GLOBAL_MATRICES%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          IF(PROBLEM_SOLUTION%SOLUTION_FINISHED) THEN
            CALL FLAG_ERROR("Problem solution has already been finished",ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
              CALL FLAG_ERROR("Solution mapping is already assocaited",ERR,ERROR,*999)
            ELSE
              NULLIFY(SOLUTION_MAPPING)
              CALL SOLUTION_MAPPING_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*999)
              SOLUTION_MAPPING=>PROBLEM_SOLUTION%SOLUTION_MAPPING
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("Global matrices problem solution is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Global matrices have not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Global matrices are not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_START")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises a solution mapping create values cache and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solution mapping create values cache 
  SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the create values cache
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,variable_idx,variable_type
    LOGICAL :: MATRIX_DONE
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Solution mapping create values cache is already associated",ERR,ERROR,*998)
      ELSE
        PROBLEM_SOLUTION=>SOLUTION_MAPPING%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
          IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
            ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache",ERR,ERROR,*999)
            ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
              & SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping create values cache matrix variable types", &
              & ERR,ERROR,*999)
            SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES=0
            !Map the first variable found to the first solver matrix, the second variable found to the second, etc.
            variable_type=1
            DO matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
              MATRIX_DONE=.FALSE.
              DO WHILE(variable_type<=FIELD_NUMBER_OF_VARIABLE_TYPES.AND..NOT.MATRIX_DONE)
                IF(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES>0) THEN                  
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,matrix_idx)=1
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1,matrix_idx)=variable_type
                  MATRIX_DONE=.TRUE.
                ELSE
                  variable_type=variable_type+1
                ENDIF
              ENDDO
              IF(.NOT.MATRIX_DONE) THEN
                !Error - could not find any more variables to map to this solver matrix
                LOCAL_ERROR="Could not find any unmapped variables for solver matrix "// &
                  & TRIM(NUMBER_TO_VSTRING(matrix_idx,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !matrix_idx
            !Check if there are still unmapped variables.
            DO variable_idx=variable_type+1,FIELD_NUMBER_OF_VARIABLE_TYPES
              IF(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_idx)%NUMBER_OF_GLOBAL_MATRICES>0) THEN
                LOCAL_ERROR="Variable type "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                  & " has not been mapped to any solver matrices"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !variable_idx            
          ELSE
            CALL FLAG_ERROR("The problem solution global matrices is not associated",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solution mapping problem solution is not associated",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroy a solution mapping.
  SUBROUTINE SOLUTION_MAPPING_DESTROY(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer the solution mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      CALL SOLUTION_MAPPING_FINALISE(SOLUTION_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLUTION_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_DESTROY",ERR,ERROR)    
    CALL EXITS("SOLUTION_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE SOLUTION_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the solution mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_FINALISE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("SOLUTION_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ALLOCATED(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS)) THEN
        DO matrix_idx=1,SIZE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS,1)
          CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS)        
      ENDIF
      CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_FINALISE(SOLUTION_MAPPING%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      DEALLOCATE(SOLUTION_MAPPING)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the mapping of global variables to a solver matrix for the solution mapping
  SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET(SOLUTION_MAPPING,SOLVER_MATRIX,VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(IN) :: SOLVER_MATRIX !<The solver matrix number to set the global matrices for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPES(:) !<The variable types to map to the solver matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mappings has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
          PROBLEM_SOLUTION=>SOLUTION_MAPPING%PROBLEM_SOLUTION
          IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
            GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
            IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
              IF(SOLVER_MATRIX>=1.AND.SOLVER_MATRIX<=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
                IF(SIZE(VARIABLE_TYPES,1)>=1.AND.SIZE(VARIABLE_TYPES,1)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                  DO variable_idx=1,SIZE(VARIABLE_TYPES,1)
                    !!TODO: CHECK THAT THE VARIABLE TYPE IS NOT REPEATED
                    IF(VARIABLE_TYPES(variable_idx)<1.OR. &
                      & VARIABLE_TYPES(variable_idx)>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                      LOCAL_ERROR="The variable type number of "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                        & " at position "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                        & " in the array is invalid. The number must be >=1 and <= "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                    IF(GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(VARIABLE_TYPES(variable_idx))%NUMBER_OF_GLOBAL_MATRICES==0) THEN
                      LOCAL_ERROR="The variable type number of "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                        & " at position "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                        & " in the array is invalid. That variable type is not mapped to any global matrices"
                    ENDIF
                  ENDDO !matrix_idx
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,SOLVER_MATRIX)=SIZE(VARIABLE_TYPES,1)
                  SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1:SIZE(VARIABLE_TYPES,1),SOLVER_MATRIX)= &
                    & VARIABLE_TYPES
                ELSE
                  LOCAL_ERROR="The supplied size of variable types array of "// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(VARIABLE_TYPES,1),"*",ERR,ERROR))// &
                    & " is invalid. The size must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The solver matrix number of "//TRIM(NUMBER_TO_VSTRING(SOLVER_MATRIX,"*",ERR,ERROR))// &
                  & " is invalid. The number must be >= 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Problem solution global matrices are not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solution mapping problem solution is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_VARIABLES_SET

  !
  !================================================================================================================================
  !

  !>Initialises the solution mapping and deallocates all memory.
  SUBROUTINE SOLUTION_MAPPING_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem_solution
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%SOLUTION_MAPPING)) THEN
        CALL FLAG_ERROR("Problem solution solution mapping is already associated",ERR,ERROR,*998)
      ELSE
        PROBLEM=>PROBLEM_SOLUTION%PROBLEM
        IF(ASSOCIATED(PROBLEM)) THEN
          DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
            ALLOCATE(PROBLEM_SOLUTION%SOLUTION_MAPPING,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution solution mapping",ERR,ERROR,*999)
            PROBLEM_SOLUTION%SOLUTION_MAPPING%PROBLEM_SOLUTION=>PROBLEM_SOLUTION
            PROBLEM_SOLUTION%SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED=.FALSE.
            PROBLEM_SOLUTION%SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES=1
            CALL SOLUTION_MAPPING_CREATE_VALUES_CACHE_INITIALISE(PROBLEM_SOLUTION%SOLUTION_MAPPING,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("The problem dependent field is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The problem solution problem is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_FINALISE(PROBLEM_SOLUTION%SOLUTION_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solver matrices for the solution mapping
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLUTION_MAPPING,NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SOLVER_MATRICES !<The number of solver matrices for the problem.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    INTEGER(INTG), ALLOCATABLE :: OLD_MATRIX_VARIABLE_TYPES(:,:)
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
        CALL FLAG_ERROR("Solution mappings has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
          PROBLEM_SOLUTION=>SOLUTION_MAPPING%PROBLEM_SOLUTION
          IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
            GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
            IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
              !Check number of matrices to set is valid
              IF(NUMBER_OF_SOLVER_MATRICES>0.AND.NUMBER_OF_SOLVER_MATRICES<=GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
                !If we need to reallocate and reset all the create values cache arrays and change the number of matrices
                IF(NUMBER_OF_SOLVER_MATRICES/=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
                  ALLOCATE(OLD_MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                    & SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old matrix variable types",ERR,ERROR,*999)
                  OLD_MATRIX_VARIABLE_TYPES=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES
                  DEALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES)
                  ALLOCATE(SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0:FIELD_NUMBER_OF_VARIABLE_TYPES, &
                    & NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate matrix variable types",ERR,ERROR,*999)
                  IF(NUMBER_OF_SOLVER_MATRICES>SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
                    SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,1:SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES)= &
                      & OLD_MATRIX_VARIABLE_TYPES(:,1:SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES)
                    DO matrix_idx=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES+1,NUMBER_OF_SOLVER_MATRICES
                      SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,matrix_idx)=1
                      SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(1,matrix_idx)=FIELD_STANDARD_VARIABLE_TYPE
                      SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(2:FIELD_NUMBER_OF_VARIABLE_TYPES,matrix_idx)=0
                    ENDDO !matrix_idx
                  ELSE
                    SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(:,1:NUMBER_OF_SOLVER_MATRICES)= &
                      & OLD_MATRIX_VARIABLE_TYPES(:,1:NUMBER_OF_SOLVER_MATRICES)
                  ENDIF
                  SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES=NUMBER_OF_SOLVER_MATRICES
                  IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
                ENDIF
              ELSE
                LOCAL_ERROR="The specified number of solver matrices of "// &
                  & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
                  & " is invalid. The number must be >= 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                  & " (the number of global matrices)"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Problem solution global matrices is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solution mapping problem solution is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_MATRIX_VARIABLE_TYPES)) DEALLOCATE(OLD_MATRIX_VARIABLE_TYPES)
    CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Finalises a solver to global map and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE(SOLVER_TO_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TO_GLOBAL_MAP_TYPE) :: SOLVER_TO_GLOBAL_MAP !<The solver matrix map to finalise the solver to global maps for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_TO_GLOBAL_MAP%GLOBAL_MATRIX_NUMBERS)) DEALLOCATE(SOLVER_TO_GLOBAL_MAP%GLOBAL_MATRIX_NUMBERS)
    IF(ALLOCATED(SOLVER_TO_GLOBAL_MAP%GLOBAL_DOFS)) DEALLOCATE(SOLVER_TO_GLOBAL_MAP%GLOBAL_DOFS)
    IF(ALLOCATED(SOLVER_TO_GLOBAL_MAP%COUPLING_COEFFICIENTS)) DEALLOCATE(SOLVER_TO_GLOBAL_MAP%COUPLING_COEFFICIENTS)
               
    CALL EXITS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solver to global map
  SUBROUTINE SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE(SOLVER_TO_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TO_GLOBAL_MAP_TYPE) :: SOLVER_TO_GLOBAL_MAP !<The solver matrix map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE",ERR,ERROR,*999)

    SOLVER_TO_GLOBAL_MAP%NUMBER_OF_GLOBAL_DOFS=0
               
    CALL EXITS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix map solver to globals map and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE(SOLVER_TO_GLOBALS_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TO_GLOBALS_MAP_TYPE) :: SOLVER_TO_GLOBALS_MAP !<The solver to globals map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,row_idx

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP)) THEN
      DO row_idx=1,SIZE(SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP,1)
        CALL SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE(SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx),ERR,ERROR,*999)
      ENDDO !row_idx
      DEALLOCATE(SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP)
    ENDIF
    IF(ALLOCATED(SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP)) THEN
      DO column_idx=1,SIZE(SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP,1)
        CALL SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_FINALISE(SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx),ERR,ERROR,*999)
      ENDDO !column_idx
      DEALLOCATE(SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP)
    ENDIF
              
    CALL EXITS("SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises a global to solver map and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE(GLOBAL_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_TO_SOLVER_MAP_TYPE), INTENT(OUT) :: GLOBAL_TO_SOLVER_MAP !<The global to solver map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(GLOBAL_TO_SOLVER_MAP%SOLUTION_DOFS)) DEALLOCATE(GLOBAL_TO_SOLVER_MAP%SOLUTION_DOFS)
    IF(ALLOCATED(GLOBAL_TO_SOLVER_MAP%COUPLING_COEFFICIENTS)) DEALLOCATE(GLOBAL_TO_SOLVER_MAP%COUPLING_COEFFICIENTS)
        
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a global to solver map 
  SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE(GLOBAL_TO_SOLVER_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_TO_SOLVER_MAP_TYPE) :: GLOBAL_TO_SOLVER_MAP !<The global to solver map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE",ERR,ERROR,*999)

    GLOBAL_TO_SOLVER_MAP%NUMBER_OF_SOLUTION_DOFS=0
        
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a global to solver maps and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE(GLOBAL_TO_SOLVER_MAPS,ERR,ERROR,*)

    !Argument variables
    TYPE(GLOBAL_TO_SOLVER_MAPS_TYPE) :: GLOBAL_TO_SOLVER_MAPS !<The global to solver maps to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,row_idx
    
    CALL ENTERS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(GLOBAL_TO_SOLVER_MAPS%GLOBAL_ROW_MAP)) THEN
      DO row_idx=1,SIZE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_ROW_MAP,1)
        CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_ROW_MAP(row_idx),ERR,ERROR,*999)
      ENDDO !row_idx
      DEALLOCATE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_ROW_MAP)
    ENDIF
    IF(ALLOCATED(GLOBAL_TO_SOLVER_MAPS%GLOBAL_COLUMN_MAP)) THEN
      DO column_idx=1,SIZE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_COLUMN_MAP,1)
        CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_FINALISE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_COLUMN_MAP(column_idx),ERR,ERROR,*999)
      ENDDO !column_idx
      DEALLOCATE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_COLUMN_MAP)
    ENDIF
    IF(ALLOCATED(GLOBAL_TO_SOLVER_MAPS%GLOBAL_COLUMN_CONSTANTS)) DEALLOCATE(GLOBAL_TO_SOLVER_MAPS%GLOBAL_COLUMN_CONSTANTS)
       
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a global to solver map
  SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE(SOLVER_MATRIX_MAPPING,GLOBAL_MATRIX_NUMBER,MATRIX_NUMBER, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_MAPPING_TYPE), TARGET :: SOLVER_MATRIX_MAPPING !<The solver matrix mapping to initialise the global to solver map for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_MATRIX_NUMBER !<The global matrix number corresponding to the matrix number to initialise
    INTEGER(INTG), INTENT(IN) :: MATRIX_NUMBER !<The matrix number to initialise the global to solver maps for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,DUMMY_ERR,NUMBER_OF_ROWS,NUMBER_OF_COLUMNS,row_idx
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR,*998)

    SOLUTION_MAPPING=>SOLVER_MATRIX_MAPPING%SOLUTION_MAPPING
    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
        PROBLEM_SOLUTION=>SOLUTION_MAPPING%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
          IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
            IF(GLOBAL_MATRIX_NUMBER>=1.AND.GLOBAL_MATRIX_NUMBER<=GLOBAL_MATRICES%NUMBER_OF_MATRICES) THEN
              SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%MATRIX_NUMBER=MATRIX_NUMBER
              SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%GLOBAL_MATRIX_NUMBER=GLOBAL_MATRIX_NUMBER
              SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%GLOBAL_MATRIX=> &
                & GLOBAL_MATRICES%MATRICES(GLOBAL_MATRIX_NUMBER)
              NUMBER_OF_ROWS=GLOBAL_MATRICES%MATRICES(GLOBAL_MATRIX_NUMBER)%NUMBER_OF_ROWS
              NUMBER_OF_COLUMNS=GLOBAL_MATRICES%MATRICES(GLOBAL_MATRIX_NUMBER)%NUMBER_OF_COLUMNS
              ALLOCATE(SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%GLOBAL_ROW_MAP(NUMBER_OF_ROWS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global to solver map global row map",ERR,ERROR,*999)
              ALLOCATE(SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%GLOBAL_COLUMN_MAP(NUMBER_OF_COLUMNS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global to solver map global column map",ERR,ERROR,*999)
              ALLOCATE(SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%GLOBAL_COLUMN_CONSTANTS(NUMBER_OF_COLUMNS), &
                & STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global to solver map global column constants",ERR,ERROR,*999)
              DO row_idx=1,NUMBER_OF_ROWS
                CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE(SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)% &
                  & GLOBAL_ROW_MAP(row_idx),ERR,ERROR,*999)
              ENDDO !row_idx
              DO column_idx=1,NUMBER_OF_COLUMNS
                CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAP_INITIALISE(SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)% &
                  & GLOBAL_COLUMN_MAP(column_idx),ERR,ERROR,*999)
              ENDDO !column_idx
              SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER)%GLOBAL_COLUMN_CONSTANTS=0.0_DP
            ELSE
              LOCAL_ERROR="The specified global matrix number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRIX_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The number must be >= 1 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(GLOBAL_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Problem solution global matrices is not associated",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution matpping problem solution is not associated",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver matrix mapping solution mapping is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE(SOLVER_MATRIX_MAPPING%GLOBAL_TO_SOLVER_MAPS(MATRIX_NUMBER),DUMMY_ERR, &
      & DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix map and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE(SOLVER_MATRIX_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_MATRIX_MAPPING_TYPE) :: SOLVER_MATRIX_MAP !<The solver matrix map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(SOLVER_MATRIX_MAP%GLOBAL_TO_SOLVER_MAPS)) THEN
      DO matrix_idx=1,SIZE(SOLVER_MATRIX_MAP%GLOBAL_TO_SOLVER_MAPS,1)
        CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_FINALISE(SOLVER_MATRIX_MAP%GLOBAL_TO_SOLVER_MAPS(matrix_idx),ERR,ERROR,*999)
      ENDDO !matrix_idx
      DEALLOCATE(SOLVER_MATRIX_MAP%GLOBAL_TO_SOLVER_MAPS)
    ENDIF
    CALL SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_FINALISE(SOLVER_MATRIX_MAP%SOLVER_TO_GLOBALS_MAP,ERR,ERROR,*999)
       
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solution mapping solver matrix maps 
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE(SOLUTION_MAPPING,SOLVER_MATRIX_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping containing the solver matrix maps
    INTEGER(INTG), INTENT(IN) :: SOLVER_MATRIX_NUMBER !<The solver matrix number of the solution mapping to intialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,GLOBAL_MATRIX_NUMBER,matrix_idx,variable_idx,variable_type
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
        PROBLEM_SOLUTION=>SOLUTION_MAPPING%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          PROBLEM=>PROBLEM_SOLUTION%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
              IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
                IF(SOLVER_MATRIX_NUMBER>=1.AND.SOLVER_MATRIX_NUMBER<=SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES) THEN
                  !Initialise preliminaries
                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%SOLVER_MATRIX_NUMBER=SOLVER_MATRIX_NUMBER
                  NULLIFY(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%SOLVER_MATRIX)
                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%SOLUTION_MAPPING=>SOLUTION_MAPPING
                  !Calculate the number of variables and global matrices mapped to this solver matrix
                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%NUMBER_OF_VARIABLES=SOLUTION_MAPPING% &
                    & CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,SOLVER_MATRIX_NUMBER)
                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%NUMBER_OF_GLOBAL_MATRICES=0
                  DO variable_idx=1,SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,SOLVER_MATRIX_NUMBER)
                    variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,SOLVER_MATRIX_NUMBER)
                    SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%NUMBER_OF_GLOBAL_MATRICES= &
                      & SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%NUMBER_OF_GLOBAL_MATRICES + &
                      & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES                    
                  ENDDO !variable_idx
                  !Allocate and initialise variables and global to solver maps
                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%VARIABLES(SOLUTION_MAPPING% &
                    & SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%NUMBER_OF_VARIABLES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrix maps variables",ERR,ERROR,*999)
                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%GLOBAL_TO_SOLVER_MAPS(SOLUTION_MAPPING% &
                    & SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrix maps global to solver maps",ERR,ERROR,*999)
                  GLOBAL_MATRIX_NUMBER=0
                  DO variable_idx=1,SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,SOLVER_MATRIX_NUMBER)
                    variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,SOLVER_MATRIX_NUMBER)
                    SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%VARIABLES(variable_idx)%PTR=> &
                      & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                    DO matrix_idx=1,GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES
                      GLOBAL_MATRIX_NUMBER=GLOBAL_MATRIX_NUMBER+1
                      CALL SOLUTION_MAPPING_GLOBAL_TO_SOLVER_MAPS_INITIALISE(SOLUTION_MAPPING% &
                        & SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER),GLOBAL_MATRIX_NUMBER,matrix_idx,ERR,ERROR,*999)
                    ENDDO !matrix_idx
                  ENDDO !variable_type
                  !Initialise the solver to globals map
                  !It does nothing at the moment
                  !CALL SOLUTION_MAPPING_SOLVER_TO_GLOBALS_MAP_INITIALISE(SOLUTION_MAPPING% &
                  !  & SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%SOLVER_TO_GLOBALS_MAP,ERR,ERROR,*999)
                  !Allocate and initialise the solver dof domain mapping
                  !!TODO: will eventually need row and column domain mappings
                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)%COLUMN_DOFS_MAPPING,STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocated solver matrix maps column domain mapping",ERR,ERROR,*999)
                  CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER)% &
                    & COLUMN_DOFS_MAPPING,DEPENDENT_FIELD%DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="Solver matrix number "//TRIM(NUMBER_TO_VSTRING(SOLVER_MATRIX_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for this solution mapping which has "// &
                    & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))//" solver matrices"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Problem solution global matrices is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Problem dependent field is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Problem solution problem is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping problem solution is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE                
        CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLVER_MATRIX_NUMBER), &
      & DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the solver matrix maps for a solution mapping
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping containing the solver matrix maps to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: column_idx,field_ny,GLOBAL_MATRIX_NUMBER,GLOBAL_MATRIX_OFFSET,global_ny,local_ny,matrix_idx, &
      & NUMBER_OF_GLOBAL_MATRICES,NUMBER_OF_GLOBAL_SOLVER_DOFS,NUMBER_OF_LOCAL_SOLVER_DOFS,rank,rank_idx,row_idx, &
      & solver_matrix_idx,variable_idx,variable_type
    INTEGER(INTG) :: VARIABLE_TYPE_LIST(FIELD_NUMBER_OF_VARIABLE_TYPES)
    LOGICAL :: RANK_DOF
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(GLOBAL_MATRICES_TYPE), POINTER :: GLOBAL_MATRICES
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION    

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ASSOCIATED(SOLUTION_MAPPING%CREATE_VALUES_CACHE)) THEN
        PROBLEM_SOLUTION=>SOLUTION_MAPPING%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          PROBLEM=>PROBLEM_SOLUTION%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            GLOBAL_MATRICES=>PROBLEM_SOLUTION%GLOBAL_MATRICES
            IF(ASSOCIATED(GLOBAL_MATRICES)) THEN
              FIXED_CONDITIONS=>PROBLEM%FIXED_CONDITIONS
              IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
                DEPENDENT_FIELD=>PROBLEM%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  DO solver_matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
                    !Calculate the number of solution variables
                    NUMBER_OF_GLOBAL_SOLVER_DOFS=0
                    NUMBER_OF_LOCAL_SOLVER_DOFS=0
!!TODO: see how slow this is. At the moment we go through number of ranks*number of globals. We could presort the global nys into a list for each rank. This would take additional memory.
                    DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
                      DO variable_idx=1,SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,solver_matrix_idx)
                        variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,solver_matrix_idx)
                        DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                        DO global_ny=1,DEPENDENT_VARIABLE%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                          RANK_DOF=.FALSE.
                          DO rank_idx=1,DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%NUMBER_OF_DOMAINS
                            IF(DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%DOMAIN_NUMBER(rank_idx)==rank &
                              & .AND.DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_TYPE(rank_idx)/= &
                              & DOMAIN_LOCAL_GHOST) THEN
                              RANK_DOF=.TRUE.
                              EXIT
                            ENDIF
                          ENDDO !rank_idx
                          IF(RANK_DOF) THEN
                            field_ny=DEPENDENT_VARIABLE%GLOBAL_DOF_LIST(global_ny)
!!TODO replace the 0 here with a named constant. Needs equations moved into its own module first.
                            IF(FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(field_ny)==0) THEN
                              NUMBER_OF_GLOBAL_SOLVER_DOFS=NUMBER_OF_GLOBAL_SOLVER_DOFS+1
                              IF(rank==COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER) &
                                & NUMBER_OF_LOCAL_SOLVER_DOFS=NUMBER_OF_LOCAL_SOLVER_DOFS+1
                            ENDIF
                          ENDIF
                        ENDDO !global_ny
                      ENDDO !variable_idx
                    ENDDO !rank
                    !Allocate memory
                    DOMAIN_MAPPING=>SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%COLUMN_DOFS_MAPPING
                    ALLOCATE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping dofs mapping global to local map", &
                      & ERR,ERROR,*999)
                    DOMAIN_MAPPING%NUMBER_OF_GLOBAL=NUMBER_OF_GLOBAL_SOLVER_DOFS
                    ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                      & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrix solver to globals solver row map",ERR,ERROR,*999)
                    ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                      & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver matrix solver to globals solver column map", &
                      & ERR,ERROR,*999)
                    !Set matrix sizes
                    SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%NUMBER_OF_ROWS=NUMBER_OF_LOCAL_SOLVER_DOFS
                    SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%TOTAL_NUMBER_OF_ROWS=NUMBER_OF_GLOBAL_SOLVER_DOFS
                    SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%NUMBER_OF_COLUMNS=NUMBER_OF_GLOBAL_SOLVER_DOFS
!!TODO: fix this
                    SOLUTION_MAPPING%NUMBER_OF_ROWS=NUMBER_OF_LOCAL_SOLVER_DOFS
                    SOLUTION_MAPPING%TOTAL_NUMBER_OF_ROWS=NUMBER_OF_GLOBAL_SOLVER_DOFS
                    !Calculate mappings
                    NUMBER_OF_GLOBAL_SOLVER_DOFS=0
                    !Loop over the ranks to ensure that the lowest ranks have the lowest numbered solver variables
!!TODO: see how slow this is. At the moment we go through number of ranks*number of globals. We could presort the global nys into a list for each rank. This would take additional memory.
                    DO rank=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
                      NUMBER_OF_LOCAL_SOLVER_DOFS=0
                      !Loop over the variables
                      GLOBAL_MATRIX_OFFSET=0
                      DO variable_idx=1,SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(0,solver_matrix_idx)
                        variable_type=SOLUTION_MAPPING%CREATE_VALUES_CACHE%MATRIX_VARIABLE_TYPES(variable_idx,solver_matrix_idx)
                        DEPENDENT_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                        DO global_ny=1,DEPENDENT_VARIABLE%DOMAIN_MAPPING%NUMBER_OF_GLOBAL
                          field_ny=DEPENDENT_VARIABLE%GLOBAL_DOF_LIST(global_ny)
                          RANK_DOF=.FALSE.
                          DO rank_idx=1,DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%NUMBER_OF_DOMAINS
                            IF(DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%DOMAIN_NUMBER(rank_idx)==rank &
                              & .AND.DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_TYPE(rank_idx)/= &
                              & DOMAIN_LOCAL_GHOST) THEN
                              local_ny=DEPENDENT_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_NUMBER(rank_idx)
                              RANK_DOF=.TRUE.
                              EXIT
                            ENDIF
                          ENDDO !rank_idx
                          IF(RANK_DOF) THEN
                            NUMBER_OF_GLOBAL_MATRICES=GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%NUMBER_OF_GLOBAL_MATRICES
                            IF(FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(field_ny)==0) THEN
                              !Map the global dof to a new solution dof
                              NUMBER_OF_GLOBAL_SOLVER_DOFS=NUMBER_OF_GLOBAL_SOLVER_DOFS+1
                              NUMBER_OF_LOCAL_SOLVER_DOFS=NUMBER_OF_LOCAL_SOLVER_DOFS+1
                              CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                                & NUMBER_OF_GLOBAL_SOLVER_DOFS),ERR,ERROR,*999)
                              !Solution dofs are not ghosted
                              ALLOCATE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%LOCAL_NUMBER(1),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dofs global to local map local number", &
                                & ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%DOMAIN_NUMBER(1),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dofs global to local map domain number", &
                                & ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%LOCAL_TYPE(1),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dofs global to local map domain number", &
                                & ERR,ERROR,*999)
                              DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%NUMBER_OF_DOMAINS=1
                              DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%LOCAL_NUMBER(1)= &
                                & NUMBER_OF_LOCAL_SOLVER_DOFS
                              DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%DOMAIN_NUMBER(1)=rank
                              DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                              !If this is my rank then set up the solver->global and global->solver mappings
                              IF(rank==COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER) THEN
                                !First set up the solver -> global matrix mappings for the new solver dof
                                !Initialise the solver to global row map
                                CALL SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
                                  & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS), &
                                  & ERR,ERROR,*999)
                                !Allocate the solver to global row map items
                                ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%GLOBAL_MATRIX_NUMBERS(NUMBER_OF_GLOBAL_MATRICES), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver row map global matrix numbers", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%GLOBAL_DOFS(NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver row map global dofs",ERR,ERROR,*999)
                                ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%COUPLING_COEFFICIENTS(NUMBER_OF_GLOBAL_MATRICES), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver row map coupling coeffficients", &
                                  & ERR,ERROR,*999)
                                SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%NUMBER_OF_GLOBAL_DOFS=NUMBER_OF_GLOBAL_MATRICES
                                !Initialise the solver to global column map
                                CALL SOLUTION_MAPPING_SOLVER_TO_GLOBAL_MAP_INITIALISE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
                                  & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS), &
                                  & ERR,ERROR,*999)
                                !Allocate the solver to global column map items
                                ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%GLOBAL_MATRIX_NUMBERS( &
                                  & NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global column map global matrix numbers", &
                                  & ERR,ERROR,*999)
                                ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%GLOBAL_DOFS(NUMBER_OF_GLOBAL_MATRICES), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global column map global dofs",ERR,ERROR,*999)
                                ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%COUPLING_COEFFICIENTS( &
                                  & NUMBER_OF_GLOBAL_MATRICES),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global column map coupling coeffficients", &
                                  & ERR,ERROR,*999)
                                SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                  & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%NUMBER_OF_GLOBAL_DOFS= &
                                  & NUMBER_OF_GLOBAL_MATRICES
                                !Loop over the global matrices associated with the variable and set the row and column maps
                                DO matrix_idx=1,NUMBER_OF_GLOBAL_MATRICES
                                  !Set the row map
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                    & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%GLOBAL_MATRIX_NUMBERS(matrix_idx)= &
                                    & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(matrix_idx)%PTR% &
                                    & MATRIX_NUMBER
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                    & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%GLOBAL_DOFS(matrix_idx)=local_ny
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                    & SOLVER_ROW_MAP(NUMBER_OF_LOCAL_SOLVER_DOFS)%COUPLING_COEFFICIENTS(matrix_idx)=1.0_DP
                                  !Set the column map
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                    & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%GLOBAL_MATRIX_NUMBERS(matrix_idx)= &
                                    & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(matrix_idx)%PTR% &
                                    & MATRIX_NUMBER
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                    & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%GLOBAL_DOFS(matrix_idx)=global_ny
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP% &
                                    & SOLVER_COLUMN_MAP(NUMBER_OF_GLOBAL_SOLVER_DOFS)%COUPLING_COEFFICIENTS(matrix_idx)= &
                                    & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(matrix_idx)%PTR% &
                                    & MATRIX_COEFFICIENT
                                ENDDO !matrix idx
                                !Secondly, set up the global -> solver mapping                                
                                DO matrix_idx=1,NUMBER_OF_GLOBAL_MATRICES
                                  !Allocate the global to solver map row items
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%NUMBER_OF_SOLUTION_DOFS=1
                                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%SOLUTION_DOFS(SOLUTION_MAPPING% &
                                    & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%NUMBER_OF_SOLUTION_DOFS),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global row map solution dofs.",ERR,ERROR,*999)
                                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%COUPLING_COEFFICIENTS( &
                                    & SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%NUMBER_OF_SOLUTION_DOFS),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global row map coupling coefficients.", &
                                    & ERR,ERROR,*999)
                                  !Allocate the global to solver map column items
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)%NUMBER_OF_SOLUTION_DOFS=1
                                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)%SOLUTION_DOFS( &
                                    & SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)% &
                                    & NUMBER_OF_SOLUTION_DOFS),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global column map solution dofs.",ERR,ERROR,*999)
                                  ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)%COUPLING_COEFFICIENTS( &
                                    & SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)% &
                                    & NUMBER_OF_SOLUTION_DOFS),STAT=ERR)
                                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global column map coupling coefficients.", &
                                    & ERR,ERROR,*999)
                                  !Do the row mapping
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%SOLUTION_DOFS(1)= &
                                    & NUMBER_OF_LOCAL_SOLVER_DOFS
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%COUPLING_COEFFICIENTS(1)= &
                                    & 1.0_DP
                                  !Do the column mapping
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)%SOLUTION_DOFS(1)= &
                                    & NUMBER_OF_GLOBAL_SOLVER_DOFS
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)%COUPLING_COEFFICIENTS(1)= &
                                    & GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)%MATRIX_MAP(matrix_idx)%PTR% &
                                    & MATRIX_COEFFICIENT
                                ENDDO !matrix_idx
                              ENDIF
                            ELSE
                              IF(rank==COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER) THEN
                                !Set up the global -> solver mapping
                                DO matrix_idx=1,NUMBER_OF_GLOBAL_MATRICES
                                  !Do the row mapping
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_ROW_MAP(local_ny)%NUMBER_OF_SOLUTION_DOFS=0
                                  !Do the column mapping
                                  SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS( &
                                    & GLOBAL_MATRIX_OFFSET+matrix_idx)%GLOBAL_COLUMN_MAP(global_ny)%NUMBER_OF_SOLUTION_DOFS=0
                                ENDDO !matrix_idx
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO !global_ny
                        GLOBAL_MATRIX_OFFSET=GLOBAL_MATRIX_OFFSET+GLOBAL_MATRICES%VARIABLE_TYPE_MAPS(variable_type)% &
                          & NUMBER_OF_GLOBAL_MATRICES
                      ENDDO !variable_idx
                    ENDDO !rank
                    CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOMAIN_MAPPING,ERR,ERROR,*999)
!!TODO: do this properly!!!
                    SOLUTION_MAPPING%ROW_DOFS_MAPPING=>SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%COLUMN_DOFS_MAPPING
                  ENDDO !solver_matrix_idx
                ELSE
                  CALL FLAG_ERROR("The problem solution dependent field is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solution fixed conditions is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solution global matrices is not associated",ERR,ERROR,*999)              
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem solution problem is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solution mapping problem solution is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution mapping create values cache is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Solution mappings:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of rows = ",SOLUTION_MAPPING%NUMBER_OF_ROWS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of rows = ",SOLUTION_MAPPING%TOTAL_NUMBER_OF_ROWS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of solver matrices = ",SOLUTION_MAPPING% &
        & NUMBER_OF_SOLVER_MATRICES,ERR,ERROR,*999)
      DO solver_matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Solver matrix : ",solver_matrix_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of rows = ",SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
          & solver_matrix_idx)%NUMBER_OF_ROWS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Total number of rows = ",SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
          & solver_matrix_idx)%TOTAL_NUMBER_OF_ROWS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of columns = ",SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
          & solver_matrix_idx)%NUMBER_OF_COLUMNS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of variables = ",SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
          & solver_matrix_idx)%NUMBER_OF_VARIABLES,ERR,ERROR,*999)
        VARIABLE_TYPE_LIST=0
        DO variable_idx=1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%NUMBER_OF_VARIABLES
          VARIABLE_TYPE_LIST(variable_idx)=SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%VARIABLES(variable_idx)%PTR% &
            & VARIABLE_TYPE
        ENDDO !variable_idx
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
          & NUMBER_OF_VARIABLES,4,4,VARIABLE_TYPE_LIST,'("    Variable types :",4(X,I8))','(20X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Global to solver mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of global matrices = ",SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
          & solver_matrix_idx)%NUMBER_OF_GLOBAL_MATRICES,ERR,ERROR,*999)
        DO matrix_idx=1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%NUMBER_OF_GLOBAL_MATRICES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global matrix : ",matrix_idx,ERR,ERROR,*999)
          GLOBAL_MATRIX_NUMBER=SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)% &
            & GLOBAL_MATRIX_NUMBER
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Global matrix number = ",GLOBAL_MATRIX_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of rows = ",GLOBAL_MATRICES%NUMBER_OF_ROWS,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Total number of rows = ",GLOBAL_MATRICES%TOTAL_NUMBER_OF_ROWS, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of columns = ",GLOBAL_MATRICES%MATRICES(matrix_idx)% &
            & NUMBER_OF_COLUMNS  ,ERR,ERROR,*999)
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"        Row maps:",ERR,ERROR,*999)
          DO row_idx=1,GLOBAL_MATRICES%NUMBER_OF_ROWS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Row number : ",row_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of solution dofs = ",SOLUTION_MAPPING% &
              & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_ROW_MAP(row_idx)% &
              & NUMBER_OF_SOLUTION_DOFS,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
              & GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_ROW_MAP(row_idx)%NUMBER_OF_SOLUTION_DOFS,5,5,SOLUTION_MAPPING% &
              & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_ROW_MAP(row_idx)%SOLUTION_DOFS, &
              & '("          Solution dof numbers  :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999) 
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
              & GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_ROW_MAP(row_idx)%NUMBER_OF_SOLUTION_DOFS,5,5,SOLUTION_MAPPING% &
              & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_ROW_MAP(row_idx)% &
              & COUPLING_COEFFICIENTS,'("          Coupling coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999) 
          ENDDO !row_idx
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"        Column maps:",ERR,ERROR,*999)
          DO column_idx=1,GLOBAL_MATRICES%MATRICES(matrix_idx)%NUMBER_OF_COLUMNS           
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Column number : ",column_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"            Number of solution dofs = ",SOLUTION_MAPPING% &
              & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_COLUMN_MAP(column_idx)% &
              & NUMBER_OF_SOLUTION_DOFS,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
              & GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_COLUMN_MAP(column_idx)%NUMBER_OF_SOLUTION_DOFS,5,5,SOLUTION_MAPPING% &
              & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_COLUMN_MAP(column_idx)% &
              & SOLUTION_DOFS,'("          Solution dof numbers  :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999) 
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
              & GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_COLUMN_MAP(column_idx)%NUMBER_OF_SOLUTION_DOFS,5,5,SOLUTION_MAPPING% &
              & SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)%GLOBAL_COLUMN_MAP(column_idx)% &
              & COUPLING_COEFFICIENTS,'("          Coupling coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999) 
          ENDDO !column_idx
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & NUMBER_OF_COLUMNS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%GLOBAL_TO_SOLVER_MAPS(matrix_idx)% &
            & GLOBAL_COLUMN_CONSTANTS,'("        Global column constants :",5(X,E13.6))','(33X,5(X,E13.6))',ERR,ERROR,*999)
        ENDDO !matrix_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Solver to global mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Row maps:",ERR,ERROR,*999)
        DO row_idx=1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%NUMBER_OF_ROWS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Row number : ",row_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of global dofs = ",SOLUTION_MAPPING% &
            & SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)% &
            & NUMBER_OF_GLOBAL_DOFS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)%NUMBER_OF_GLOBAL_DOFS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
            & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)%GLOBAL_MATRIX_NUMBERS, &
            & '("          Global matrix numbers :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)%NUMBER_OF_GLOBAL_DOFS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
            & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)%GLOBAL_DOFS, &
            & '("          Global dof numbers    :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)%NUMBER_OF_GLOBAL_DOFS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
            & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_ROW_MAP(row_idx)%COUPLING_COEFFICIENTS, &
            & '("          Coupling coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999)
        ENDDO !row_idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"      Column maps:",ERR,ERROR,*999)
        DO column_idx=1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)%NUMBER_OF_COLUMNS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Column number : ",column_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Number of global dofs = ",SOLUTION_MAPPING% &
            & SOLVER_MATRIX_MAPS(solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)% &
            & NUMBER_OF_GLOBAL_DOFS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)%NUMBER_OF_GLOBAL_DOFS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
            & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)%GLOBAL_MATRIX_NUMBERS, &
            & '("          Global matrix numbers :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)%NUMBER_OF_GLOBAL_DOFS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
            & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)%GLOBAL_DOFS, &
            & '("          Global dof numbers    :",5(X,I13))','(31X,5(X,I13))',ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(solver_matrix_idx)% &
            & SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)%NUMBER_OF_GLOBAL_DOFS,5,5,SOLUTION_MAPPING%SOLVER_MATRIX_MAPS( &
            & solver_matrix_idx)%SOLVER_TO_GLOBALS_MAP%SOLVER_COLUMN_MAP(column_idx)%COUPLING_COEFFICIENTS, &
            & '("          Coupling coefficients :",5(X,E13.6))','(31X,5(X,E13.6))',ERR,ERROR,*999)
        ENDDO !column_idx
      ENDDO !solver_matrix_idx
    ENDIF
    
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises a solver matrix mapping and deallocates all memory
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping containing the solver matrix mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      IF(ALLOCATED(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS)) THEN
        DO matrix_idx=1,SIZE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS,1)
          CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAP_FINALISE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(matrix_idx),ERR,ERROR,*999)
        ENDDO !matrix_idx
        DEALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS)
      ENDIF
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a solution mapping solver matrix maps 
  SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE(SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<A pointer to the solution mapping containing the solver matrix maps
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
      ALLOCATE(SOLUTION_MAPPING%SOLVER_MATRIX_MAPS(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution mapping solver matrix mappings",ERR,ERROR,*999)
      DO matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
        CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAP_INITIALISE(SOLUTION_MAPPING,matrix_idx,ERR,ERROR,*999)
      ENDDO !matrix_idx      
    ELSE
      CALL FLAG_ERROR("Solution mapping is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE")
    RETURN
999 CALL SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_FINALISE(SOLUTION_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_MAPPING_SOLVER_MATRIX_MAPS_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE SOLUTION_MAPPING_ROUTINES
