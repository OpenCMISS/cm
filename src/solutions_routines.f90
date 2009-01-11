!> \file
!> $Id: solutions_routines.f90 248 2008-11-28 11:14:17Z chrispbradley $
!> \author Chris Bradley
!> \brief This module handles all solutions routines.
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

!> This module handles all solutions routines.
MODULE SOLUTIONS_ROUTINES

  USE BASE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLUTION_SOLUTION_MAPPING_GET,SOLUTION_SOLVER_GET,SOLUTIONS_CREATE_FINISH,SOLUTIONS_CREATE_START,SOLUTIONS_DESTROY, &
    & SOLUTIONS_LINEARITY_SET,SOLUTIONS_NUMBER_SET,SOLUTIONS_SOLUTION_GET,SOLUTIONS_TIME_DEPENDENCE_SET


CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalises a solution and deallocates all memory.
  SUBROUTINE SOLUTION_FINALISE(SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLUTION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(ASSOCIATED(SOLUTION%SOLUTION_MAPPING)) CALL SOLUTION_MAPPING_DESTROY(SOLUTION%SOLUTION_MAPPING,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLUTION%SOLVER)) CALL SOLVER_DESTROY(SOLUTION%SOLVER,ERR,ERROR,*999)
      DEALLOCATE(SOLUTION)
    ENDIF
       
    CALL EXITS("SOLUTION_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTION_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solution for a problem.
  SUBROUTINE SOLUTION_INITIALISE(SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLUTION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      SOLUTION%SOLUTION_NUMBER=0
      NULLIFY(SOLUTION%CONTROL_LOOP)
      SOLUTION%SOLUTION_FINISHED=.FALSE.
      SOLUTION%LINEARITY=PROBLEM_SOLUTION_LINEAR
      SOLUTION%TIME_DEPENDENCE=PROBLEM_SOLUTION_STATIC
      NULLIFY(SOLUTION%SOLUTION_MAPPING)
      NULLIFY(SOLUTION%SOLVER)
    ELSE
      CALL FLAG_ERROR("Solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLUTION_INITIALISE")
    RETURN
999 CALL ERRORS("SOLUTION_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTION_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTION_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solution mapping for a solution.
  SUBROUTINE SOLUTION_SOLUTION_MAPPING_GET(SOLUTION,SOLUTION_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to get the solution mapping for
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING !<On exit, a pointer to the solution mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLUTION_SOLUTION_MAPPING_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(SOLUTION%SOLUTION_FINISHED) THEN
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
          CALL FLAG_ERROR("Solution mapping is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(SOLUTION_MAPPING)
          SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
          IF(.NOT.ASSOCIATED(SOLUTION_MAPPING)) CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTION_SOLUTION_MAPPING_GET")
    RETURN
999 NULLIFY(SOLUTION_MAPPING)
998 CALL ERRORS("SOLUTION_SOLUTION_MAPPING_GET",ERR,ERROR)
    CALL EXITS("SOLUTION_SOLUTION_MAPPING_GET")
    RETURN 1
  END SUBROUTINE SOLUTION_SOLUTION_MAPPING_GET
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a solution.
  SUBROUTINE SOLUTION_SOLVER_GET(SOLUTION,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to get the solver for
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On exit, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLUTION_SOLVER_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(SOLUTION%SOLUTION_FINISHED) THEN
        IF(ASSOCIATED(SOLVER)) THEN
          CALL FLAG_ERROR("Solver is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(SOLVER)
          SOLVER=>SOLUTION%SOLVER
          IF(.NOT.ASSOCIATED(SOLVER)) CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLUTION_SOLVER_GET")
    RETURN
999 NULLIFY(SOLVER)
998 CALL ERRORS("SOLUTION_SOLVER_GET",ERR,ERROR)
    CALL EXITS("SOLUTION_SOLVER_GET")
    RETURN 1
  END SUBROUTINE SOLUTION_SOLVER_GET
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of solutions.
  SUBROUTINE SOLUTIONS_CREATE_FINISH(SOLUTIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
   
    CALL ENTERS("SOLUTIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      IF(SOLUTIONS%SOLUTIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solutions has already been finished.",ERR,ERROR,*999)
      ELSE        
        CONTROL_LOOP=>SOLUTIONS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN          
          !Finish the problem solution creation
          IF(ALLOCATED(SOLUTIONS%SOLUTIONS)) THEN
            DO solution_idx=1,SOLUTIONS%NUMBER_OF_SOLUTIONS
              SOLUTION=>SOLUTIONS%SOLUTIONS(solution_idx)%PTR
              IF(ASSOCIATED(SOLUTION)) THEN
                SOLUTION%SOLUTION_FINISHED=.TRUE.
              ELSE
                CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !solution_idx            
            SOLUTIONS%SOLUTIONS_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Solutions solutions is not allocated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solutions control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solutions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLUTIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLUTIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("SOLUTIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE SOLUTIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a solution for the control loop. 
  SUBROUTINE SOLUTIONS_CREATE_START(CONTROL_LOOP,SOLUTIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to create the solutions for
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<On exit, a pointer to the solutions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLUTIONS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
          IF(ASSOCIATED(SOLUTIONS)) THEN
            CALL FLAG_ERROR("Solutions is already associated.",ERR,ERROR,*999)
          ELSE
            NULLIFY(SOLUTIONS)
            !Initialise the solutions
            CALL SOLUTIONS_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
            !Return the pointer
            SOLUTIONS=>CONTROL_LOOP%SOLUTIONS
          ENDIF
        ELSE
          LOCAL_ERROR="Invalid control loop setup. The specified control loop has "// &
            & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))// &
            & " sub loops. To create solutions the control loop must have 0 sub loops."          
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Control loop has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLUTIONS_CREATE_START")
    RETURN
999 CALL ERRORS("SOLUTIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLUTIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE SOLUTIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys the solutions
  SUBROUTINE SOLUTIONS_DESTROY(SOLUTIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("SOLUTIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      CALL SOLUTIONS_FINALISE(SOLUTIONS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solutions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLUTIONS_DESTROY")
    RETURN
999 CALL ERRORS("SOLUTIONS_DESTROY",ERR,ERROR)
    CALL EXITS("SOLUTIONS_DESTROY")
    RETURN 1
  END SUBROUTINE SOLUTIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the solutions and deallocates all memory
  SUBROUTINE SOLUTIONS_FINALISE(SOLUTIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
 
    CALL ENTERS("SOLUTIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      IF(ALLOCATED(SOLUTIONS%SOLUTIONS)) THEN
        DO solution_idx=1,SIZE(SOLUTIONS%SOLUTIONS,1)
          CALL SOLUTION_FINALISE(SOLUTIONS%SOLUTIONS(solution_idx)%PTR,ERR,ERROR,*999)
        ENDDO ! solution_idx
        DEALLOCATE(SOLUTIONS%SOLUTIONS)
      ENDIF
      DEALLOCATE(SOLUTIONS)
    ENDIF
       
    CALL EXITS("SOLUTIONS_FINALISE")
    RETURN
999 CALL ERRORS("SOLUTIONS_FINALISE",ERR,ERROR)
    CALL EXITS("SOLUTIONS_FINALISE")
    RETURN 1
  END SUBROUTINE SOLUTIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the solutions for a problem.
  SUBROUTINE SOLUTIONS_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to initialise the solutions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,solution_idx
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLUTIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      SOLUTIONS=>CONTROL_LOOP%SOLUTIONS
      IF(ASSOCIATED(SOLUTIONS)) THEN
        CALL FLAG_ERROR("Solutions is already allocated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%SOLUTIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate control loop solutions.",ERR,ERROR,*999)
        CONTROL_LOOP%SOLUTIONS%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%SOLUTIONS%SOLUTIONS_FINISHED=.FALSE.
        CONTROL_LOOP%SOLUTIONS%NUMBER_OF_SOLUTIONS=1
        ALLOCATE(CONTROL_LOOP%SOLUTIONS%SOLUTIONS(CONTROL_LOOP%SOLUTIONS%NUMBER_OF_SOLUTIONS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solutions solutions.",ERR,ERROR,*999)
        DO solution_idx=1,CONTROL_LOOP%SOLUTIONS%NUMBER_OF_SOLUTIONS
          ALLOCATE(CONTROL_LOOP%SOLUTIONS%SOLUTIONS(solution_idx)%PTR,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution.",ERR,ERROR,*999)
          CALL SOLUTION_INITIALISE(CONTROL_LOOP%SOLUTIONS%SOLUTIONS(solution_idx)%PTR,ERR,ERROR,*999)
          CONTROL_LOOP%SOLUTIONS%SOLUTIONS(solution_idx)%PTR%CONTROL_LOOP=>CONTROL_LOOP
        ENDDO !solution_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTIONS_INITIALISE")
    RETURN
999 CALL SOLUTIONS_FINALISE(SOLUTIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLUTIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLUTIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE SOLUTIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity of a solution.
  SUBROUTINE SOLUTIONS_LINEARITY_SET(SOLUTIONS,SOLUTION_INDEX,LINEARITY,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to set the linearity for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_INDEX !<The index of the solution to set the linearity for
    INTEGER(INTG), INTENT(IN) :: LINEARITY !<The solution linearity to set \see PROBLEM_CONSTANTS_SolutionLinearityTypes,PROBLEM_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("SOLUTIONS_LINEARITY_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      IF(SOLUTIONS%SOLUTIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solutions have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLUTION_INDEX>0.AND.SOLUTION_INDEX<=SOLUTIONS%NUMBER_OF_SOLUTIONS) THEN
          SOLUTION=>SOLUTIONS%SOLUTIONS(SOLUTION_INDEX)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SELECT CASE(LINEARITY)
            CASE(PROBLEM_SOLUTION_LINEAR)
              SOLUTION%LINEARITY=PROBLEM_SOLUTION_LINEAR
            CASE(PROBLEM_SOLUTION_NONLINEAR)
              SOLUTION%LINEARITY=PROBLEM_SOLUTION_NONLINEAR
            CASE DEFAULT
              LOCAL_ERROR="The specified linearity of "//TRIM(NUMBER_TO_VSTRING(LINEARITY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
          ENDIF          
        ELSE
          LOCAL_ERROR="The specified solution index of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The number of solutions must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(SOLUTIONS%NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solutions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLUTIONS_LINEARITY_SET")
    RETURN
999 CALL ERRORS("SOLUTIONS_LINEARITY_SET",ERR,ERROR)
    CALL EXITS("SOLUTIONS_LINEARITY_SET")
    RETURN 1
  END SUBROUTINE SOLUTIONS_LINEARITY_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solutions.
  SUBROUTINE SOLUTIONS_NUMBER_SET(SOLUTIONS,NUMBER_OF_SOLUTIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to set the number for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SOLUTIONS !<The number of solutions to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    TYPE(SOLUTION_PTR_TYPE), ALLOCATABLE :: OLD_SOLUTIONS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("SOLUTIONS_NUMBER_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      IF(SOLUTIONS%SOLUTIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solutions have already been finished.",ERR,ERROR,*998)
      ELSE
        IF(NUMBER_OF_SOLUTIONS>0) THEN
          IF(NUMBER_OF_SOLUTIONS/=SOLUTIONS%NUMBER_OF_SOLUTIONS) THEN
            ALLOCATE(OLD_SOLUTIONS(SOLUTIONS%NUMBER_OF_SOLUTIONS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old solutions.",ERR,ERROR,*999)
            DO solution_idx=1,SOLUTIONS%NUMBER_OF_SOLUTIONS
              OLD_SOLUTIONS(solution_idx)%PTR=>SOLUTIONS%SOLUTIONS(solution_idx)%PTR
            ENDDO !solution_idx
            IF(ALLOCATED(SOLUTIONS%SOLUTIONS)) DEALLOCATE(SOLUTIONS%SOLUTIONS)
            ALLOCATE(SOLUTIONS%SOLUTIONS(NUMBER_OF_SOLUTIONS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solutions.",ERR,ERROR,*999)
            IF(NUMBER_OF_SOLUTIONS>SOLUTIONS%NUMBER_OF_SOLUTIONS) THEN
              DO solution_idx=1,SOLUTIONS%NUMBER_OF_SOLUTIONS
                SOLUTIONS%SOLUTIONS(solution_idx)%PTR=>OLD_SOLUTIONS(solution_idx)%PTR
              ENDDO !solution_idx
              DO solution_idx=SOLUTIONS%NUMBER_OF_SOLUTIONS+1,NUMBER_OF_SOLUTIONS
                ALLOCATE(SOLUTIONS%SOLUTIONS(solution_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution.",ERR,ERROR,*999)
                CALL SOLUTION_INITIALISE(SOLUTIONS%SOLUTIONS(solution_idx)%PTR,ERR,ERROR,*999)
              ENDDO !solution_idx
            ELSE
              DO solution_idx=1,NUMBER_OF_SOLUTIONS
                SOLUTIONS%SOLUTIONS(solution_idx)%PTR=>OLD_SOLUTIONS(solution_idx)%PTR
              ENDDO !solution_idx
              DO solution_idx=NUMBER_OF_SOLUTIONS+1,SOLUTIONS%NUMBER_OF_SOLUTIONS
                CALL SOLUTION_FINALISE(OLD_SOLUTIONS(solution_idx)%PTR,ERR,ERROR,*999)
              ENDDO !solution_idx
            ENDIF
            SOLUTIONS%NUMBER_OF_SOLUTIONS=NUMBER_OF_SOLUTIONS
          ENDIF
        ELSE
          LOCAL_ERROR="The specified number of solutions of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))// &
            & " is invalid. The number of solutions must be > 0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solutions is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTIONS_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_SOLUTIONS)) DEALLOCATE(OLD_SOLUTIONS)
998 CALL ERRORS("SOLUTIONS_NUMBER_SET",ERR,ERROR)
    CALL EXITS("SOLUTIONS_NUMBER_SET")
    RETURN 1
  END SUBROUTINE SOLUTIONS_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified solution in the list of solutions.
  SUBROUTINE SOLUTIONS_SOLUTION_GET(SOLUTIONS,SOLUTION_INDEX,SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to get the solution for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_INDEX !<The specified solution to get
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<On exit, a pointer to the specified solution. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("SOLUTIONS_SOLUTION_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      IF(SOLUTIONS%SOLUTIONS_FINISHED) THEN
        IF(ASSOCIATED(SOLUTION)) THEN
          CALL FLAG_ERROR("Solution is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(SOLUTION)
          IF(SOLUTION_INDEX>0.AND.SOLUTION_INDEX<=SOLUTIONS%NUMBER_OF_SOLUTIONS) THEN
            IF(ALLOCATED(SOLUTIONS%SOLUTIONS)) THEN
              SOLUTION=>SOLUTIONS%SOLUTIONS(SOLUTION_INDEX)%PTR
              IF(.NOT.ASSOCIATED(SOLUTION)) CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Solutions solutions is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified solution index of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The solution index must be >= 1 and <= "// &
              & TRIM(NUMBER_TO_VSTRING(SOLUTIONS%NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solutions has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solutions is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLUTIONS_SOLUTION_GET")
    RETURN
999 NULLIFY(SOLUTION)
998 CALL ERRORS("SOLUTIONS_SOLUTION_GET",ERR,ERROR)
    CALL EXITS("SOLUTIONS_SOLUTION_GET")
    RETURN 1
  END SUBROUTINE SOLUTIONS_SOLUTION_GET
  
  !
  !================================================================================================================================
  !
 
  !>Sets/changes the time dependence of a solution.
  SUBROUTINE SOLUTIONS_TIME_DEPENDENCE_SET(SOLUTIONS,SOLUTION_INDEX,TIME_DEPENDENCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<A pointer to the solutions to set the time dependence for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_INDEX !<The index of the solution to set the time dependence for
    INTEGER(INTG), INTENT(IN) :: TIME_DEPENDENCE !<The solution time dependence to set \see PROBLEM_CONSTANTS_SolutionTimeDependenceTypes,PROBLEM_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("SOLUTIONS_TIME_DEPENDENCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTIONS)) THEN
      IF(SOLUTIONS%SOLUTIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solutions have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLUTION_INDEX>0.AND.SOLUTION_INDEX<=SOLUTIONS%NUMBER_OF_SOLUTIONS) THEN
          SOLUTION=>SOLUTIONS%SOLUTIONS(SOLUTION_INDEX)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SELECT CASE(TIME_DEPENDENCE)
            CASE(PROBLEM_SOLUTION_STATIC)
              SOLUTION%TIME_DEPENDENCE=PROBLEM_SOLUTION_STATIC
            CASE(PROBLEM_SOLUTION_QUASISTATIC)
              SOLUTION%TIME_DEPENDENCE=PROBLEM_SOLUTION_QUASISTATIC
            CASE(PROBLEM_SOLUTION_FIRST_ORDER_DYNAMIC)
              SOLUTION%TIME_DEPENDENCE=PROBLEM_SOLUTION_FIRST_ORDER_DYNAMIC
            CASE(PROBLEM_SOLUTION_SECOND_ORDER_DYNAMIC)
              SOLUTION%TIME_DEPENDENCE=PROBLEM_SOLUTION_SECOND_ORDER_DYNAMIC
            CASE(PROBLEM_SOLUTION_THIRD_ORDER_DYNAMIC)
              SOLUTION%TIME_DEPENDENCE=PROBLEM_SOLUTION_THIRD_ORDER_DYNAMIC
            CASE DEFAULT
              LOCAL_ERROR="The specified time dependence of "//TRIM(NUMBER_TO_VSTRING(TIME_DEPENDENCE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
          ENDIF          
        ELSE
          LOCAL_ERROR="The specified solution index of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The number of solutions must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(SOLUTIONS%NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solutions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLUTIONS_TIME_DEPENDENCE_SET")
    RETURN
999 CALL ERRORS("SOLUTIONS_TIME_DEPENDENCE_SET",ERR,ERROR)
    CALL EXITS("SOLUTIONS_TIME_DEPENDENCE_SET")
    RETURN 1
  END SUBROUTINE SOLUTIONS_TIME_DEPENDENCE_SET
  
  !
  !================================================================================================================================
  !

END MODULE SOLUTIONS_ROUTINES

