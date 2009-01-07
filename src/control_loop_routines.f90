!> \file
!> $Id: control_loop_routines.f90 248 2008-11-28 11:14:17Z chrispbradley $
!> \author Chris Bradley
!> \brief This module handles all control loop routines.
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

!> This module handles all control loop routines.
MODULE CONTROL_LOOP_ROUTINES

  USE BASE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE SOLUTIONS_ROUTINES
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup CONTROL_LOOP_ROUTINES_ControlLoopIdentifiers CONTROL_LOOP_ROUTINES::ControlLoopIdentifiers
  !> \brief The linearity type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_NODE=0 !<The identifier for a each "leaf" node in a control loop. \see CONTROL_LOOP_ROUTINES_ControlLoopIdentifiers,CONTROL_LOOP_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  INTERFACE CONTROL_LOOP_GET
    MODULE PROCEDURE CONTROL_LOOP_GET_0
    MODULE PROCEDURE CONTROL_LOOP_GET_1
  END INTERFACE !CONTROL_LOOP_GET

  PUBLIC CONTROL_LOOP_NODE

  PUBLIC CONTROL_LOOP_CREATE_FINISH,CONTROL_LOOP_CREATE_START,CONTROL_LOOP_DESTROY,CONTROL_LOOP_GET,CONTROL_LOOP_ITERATIONS_SET, &
    & CONTROL_LOOP_MAXIMUM_ITERATIONS_SET,CONTROL_LOOP_NUMBER_SUB_LOOPS_GET,CONTROL_LOOP_NUMBER_SUB_LOOPS_SET, &
    & CONTROL_LOOP_SUB_LOOP_GET,CONTROL_LOOP_SOLUTIONS_DESTROY,CONTROL_LOOP_SOLUTIONS_GET,CONTROL_LOOP_SOLVER_DESTROY, &
    & CONTROL_LOOP_TIMES_SET,CONTROL_LOOP_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the process of creating a control loop
  RECURSIVE SUBROUTINE CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to finish.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
    
    CALL ENTERS("CONTROL_LOOP_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        !Finish the sub-loops first
        IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
          DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
            CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP2,ERR,ERROR,*999)
          ENDDO !loop_idx
        ENDIF
        !Finish this control loop
        CONTROL_LOOP%CONTROL_LOOP_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_CREATE_FINISH")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the process of creating a control loop for a problem.
  SUBROUTINE CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to initialise the control for.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On exit, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("CONTROL_LOOP_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("Control loop is already associated.",ERR,ERROR,*998)
      ELSE
        NULLIFY(CONTROL_LOOP)
        CALL CONTROL_LOOP_INITIALISE(PROBLEM,ERR,ERROR,*999)
        CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*998)
    ENDIF
              
    CALL EXITS("CONTROL_LOOP_CREATE_START")
    RETURN
999 CALL CONTROL_LOOP_FINALISE(PROBLEM%CONTROL_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONTROL_LOOP_CREATE_START",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_CREATE_START")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy a control loop
  SUBROUTINE CONTROL_LOOP_DESTROY(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   
    CALL ENTERS("CONTROL_LOOP_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_FINALISE(CONTROL_LOOP,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_DESTROY")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_DESTROY",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_DESTROY")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise a control loop and deallocate all memory
  RECURSIVE SUBROUTINE CONTROL_LOOP_FINALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
 
    CALL ENTERS("CONTROL_LOOP_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      !Finalise any sub control loops first
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
        DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
          CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
          CALL CONTROL_LOOP_FINALISE(CONTROL_LOOP2,ERR,ERROR,*999)
        ENDDO !loop_idx
        DEALLOCATE(CONTROL_LOOP%SUB_LOOPS)
      ENDIF
      !Finalise any solutions
      IF(ASSOCIATED(CONTROL_LOOP%SOLUTIONS)) CALL SOLUTIONS_DESTROY(CONTROL_LOOP%SOLUTIONS,ERR,ERROR,*999)
      !Now finalise this control loop
      CALL CONTROL_LOOP_SIMPLE_FINALISE(CONTROL_LOOP%SIMPLE_LOOP,ERR,ERROR,*999)
      CALL CONTROL_LOOP_FIXED_FINALISE(CONTROL_LOOP%FIXED_LOOP,ERR,ERROR,*999)
      CALL CONTROL_LOOP_TIME_FINALISE(CONTROL_LOOP%TIME_LOOP,ERR,ERROR,*999)
      CALL CONTROL_LOOP_WHILE_FINALISE(CONTROL_LOOP%WHILE_LOOP,ERR,ERROR,*999)
      DEALLOCATE(CONTROL_LOOP)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_FINALISE")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_FINALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_FINALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise the control for a problem.
  SUBROUTINE CONTROL_LOOP_INITIALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to initialise the control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("CONTROL_LOOP_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("Control loop is already associated for this problem.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(PROBLEM%CONTROL_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem control loop.",ERR,ERROR,*999)
        PROBLEM%CONTROL_LOOP%PROBLEM=>PROBLEM
        NULLIFY(PROBLEM%CONTROL_LOOP%PARENT_LOOP)
        PROBLEM%CONTROL_LOOP%CONTROL_LOOP_FINISHED=.FALSE.
        PROBLEM%CONTROL_LOOP%LOOP_TYPE=PROBLEM_CONTROL_SIMPLE_TYPE
        PROBLEM%CONTROL_LOOP%CONTROL_LOOP_LEVEL=1
        NULLIFY(PROBLEM%CONTROL_LOOP%SIMPLE_LOOP)
        NULLIFY(PROBLEM%CONTROL_LOOP%FIXED_LOOP)
        NULLIFY(PROBLEM%CONTROL_LOOP%TIME_LOOP)
        NULLIFY(PROBLEM%CONTROL_LOOP%WHILE_LOOP)
        PROBLEM%CONTROL_LOOP%NUMBER_OF_SUB_LOOPS=0
        NULLIFY(PROBLEM%CONTROL_LOOP%SOLUTIONS)
        CALL CONTROL_LOOP_SIMPLE_INITIALISE(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*998)
    ENDIF
              
    CALL EXITS("CONTROL_LOOP_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_FINALISE(PROBLEM%CONTROL_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONTROL_LOOP_INITIALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_INITIALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a fixed control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_FIXED_FINALISE(FIXED_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: FIXED_LOOP !<A pointer to the fixed control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_FIXED_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIXED_LOOP)) THEN
      DEALLOCATE(FIXED_LOOP)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_FIXED_FINALISE")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_FIXED_FINALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_FIXED_FINALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_FIXED_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  SUBROUTINE CONTROL_LOOP_GET_0(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT !<A pointer to the control loop to root
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On exit, the specified control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_GET_0",ERR,ERROR,*999)

    CALL CONTROL_LOOP_GET_1(CONTROL_LOOP_ROOT,(/CONTROL_LOOP_IDENTIFIER/),CONTROL_LOOP,ERR,ERROR,*999)
       
    CALL EXITS("CONTROL_LOOP_GET_0")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_GET_0",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_GET_0")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_GET_0

  !
  !================================================================================================================================
  !

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  SUBROUTINE CONTROL_LOOP_GET_1(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT !<A pointer to the control loop to root
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On exit, the specified control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: control_loop_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONTROL_LOOP_GET_1",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
      !IF(CONTROL_LOOP_ROOT%CONTROL_LOOP_FINISHED) THEN
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          CALL FLAG_ERROR("Control loop is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(CONTROL_LOOP)
          IF(COUNT(CONTROL_LOOP_IDENTIFIER==CONTROL_LOOP_NODE)==1) THEN
            IF(CONTROL_LOOP_IDENTIFIER(SIZE(CONTROL_LOOP_IDENTIFIER,1))==CONTROL_LOOP_NODE) THEN
              CONTROL_LOOP=>CONTROL_LOOP_ROOT
              DO control_loop_idx=1,SIZE(CONTROL_LOOP_IDENTIFIER,1)
                IF(CONTROL_LOOP_IDENTIFIER(control_loop_idx)==CONTROL_LOOP_NODE) THEN
                  EXIT
                ELSE
                  IF(CONTROL_LOOP_IDENTIFIER(control_loop_idx)>0.AND. &
                    & CONTROL_LOOP_IDENTIFIER(control_loop_idx)<=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
                    CONTROL_LOOP=>CONTROL_LOOP%SUB_LOOPS(CONTROL_LOOP_IDENTIFIER(control_loop_idx))%PTR
                    IF(.NOT.ASSOCIATED(CONTROL_LOOP)) THEN
                      LOCAL_ERROR="Control sub loop number "// &
                        & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP_IDENTIFIER(control_loop_idx),"*",ERR,ERROR))// &
                        & " at identifier index "//TRIM(NUMBER_TO_VSTRING(control_loop_idx,"*",ERR,ERROR))// &
                        & " is not associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Invalid control loop identifier. The identifier at index "// &
                      & TRIM(NUMBER_TO_VSTRING(control_loop_idx,"*",ERR,ERROR))//" is "// &
                      & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP_IDENTIFIER(control_loop_idx),"*",ERR,ERROR))// &
                      & ". The identifier must be between 1 and "// &
                      & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !control_loop_idx
            ELSE
              LOCAL_ERROR="Invalid control loop identifier. The last value in the identifier vector is "// &
                & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP_IDENTIFIER(SIZE(CONTROL_LOOP_IDENTIFIER,1)),"*",ERR,ERROR))// &
                & " and it should be "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP_NODE,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="Invalid control loop identifier. The control loop identifier has "// &
              & TRIM(NUMBER_TO_VSTRING(COUNT(CONTROL_LOOP_IDENTIFIER==CONTROL_LOOP_NODE),"*",ERR,ERROR))// &
              & " control loop node identifiers and it should only have 1."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      !ELSE
      !  CALL FLAG_ERROR("Control loop root has not been finished.",ERR,ERROR,*998)
      !ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop root is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_GET_1")
    RETURN
999 NULLIFY(CONTROL_LOOP)
998 CALL ERRORS("CONTROL_LOOP_GET_1",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_GET_1")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_GET_1

  !
  !================================================================================================================================
  !

  !>Initialises a fixed loop for a control loop.
  SUBROUTINE CONTROL_LOOP_FIXED_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to initialise the fixed loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONTROL_LOOP_FIXED_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%FIXED_LOOP)) THEN
        CALL FLAG_ERROR("The fixed loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%FIXED_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate fixed loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%FIXED_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%FIXED_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%FIXED_LOOP%START_ITERATION=1
        CONTROL_LOOP%FIXED_LOOP%STOP_ITERATION=100
        CONTROL_LOOP%FIXED_LOOP%ITERATION_INCREMENT=1
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_FIXED_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_FIXED_FINALISE(CONTROL_LOOP%FIXED_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONTROL_LOOP_FIXED_INITIALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_FIXED_INITIALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_FIXED_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the iteration parameters for a fixed control loop.
  SUBROUTINE CONTROL_LOOP_ITERATIONS_SET(CONTROL_LOOP,START_ITERATION,STOP_ITERATION,ITERATION_INCREMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to fixed control loop to set the iterations for
    INTEGER(INTG), INTENT(IN) :: START_ITERATION !<The start iteration for the fixed control loop.
    INTEGER(INTG), INTENT(IN) :: STOP_ITERATION !<The stop iteration for the fixed control loop.
    INTEGER(INTG), INTENT(IN) :: ITERATION_INCREMENT !<The iteration increment for the fixed control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: FIXED_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONTROL_LOOP_TIMES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_FIXED_LOOP_TYPE) THEN
          FIXED_LOOP=>CONTROL_LOOP%FIXED_LOOP
          IF(ASSOCIATED(FIXED_LOOP)) THEN
            IF(ITERATION_INCREMENT==0) THEN
              LOCAL_ERROR="The specified time increment of "//TRIM(NUMBER_TO_VSTRING(ITERATION_INCREMENT,"*",ERR,ERROR))// &
                & " is invalid. The iteration increment must not be zero."          
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              IF(ITERATION_INCREMENT>0) THEN
                IF(STOP_ITERATION<=START_ITERATION) THEN
                  LOCAL_ERROR="The specified stop iteration of "//TRIM(NUMBER_TO_VSTRING(STOP_ITERATION,"*",ERR,ERROR))// &
                    & " is incompatiable with a specified start increment of "// &
                    & TRIM(NUMBER_TO_VSTRING(START_ITERATION,"*",ERR,ERROR))// &
                    & ". For a positive iteration increment the stop iteration must be > than the start iteration."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                IF(START_ITERATION<=STOP_ITERATION) THEN
                  LOCAL_ERROR="The specified start iteration of "//TRIM(NUMBER_TO_VSTRING(START_ITERATION,"*",ERR,ERROR))// &
                    & " is incompatiable with a specified stop iteration of "// &
                    & TRIM(NUMBER_TO_VSTRING(STOP_ITERATION,"*",ERR,ERROR))// &
                    & ". For a negative iteration increment the stop iteration must be < than the start iteration."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDIF
            FIXED_LOOP%START_ITERATION=START_ITERATION
            FIXED_LOOP%STOP_ITERATION=STOP_ITERATION
            FIXED_LOOP%ITERATION_INCREMENT=ITERATION_INCREMENT
          ELSE
            CALL FLAG_ERROR("Control loop fixed loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified control loop is not a fixed control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_ITERATIONS_SET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_ITERATIONS_SET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_ITERATIONS_SET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_ITERATIONS_SET
  
  !
  !================================================================================================================================
  !

  !>Sets the maximum number of iterations for a while control loop.
  SUBROUTINE CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MAXIMUM_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to while control loop to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_ITERATIONS !<The maximum number of iterations for a while control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: WHILE_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONTROL_LOOP_TIMES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
          WHILE_LOOP=>CONTROL_LOOP%WHILE_LOOP
          IF(ASSOCIATED(WHILE_LOOP)) THEN
            IF(MAXIMUM_ITERATIONS<=0) THEN
              LOCAL_ERROR="The specified maximum number of iterations of "// &
                & TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                & " is invalid. The maximum number of iterations must be greater than zero."          
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
            ENDIF
            WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
          ELSE
            CALL FLAG_ERROR("Control loop while loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified control loop is not a while control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_MAXIMUM_ITERATIONS_SET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_MAXIMUM_ITERATIONS_SET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_MAXIMUM_ITERATIONS_SET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_MAXIMUM_ITERATIONS_SET

  !
  !================================================================================================================================
  !

  !>Gets the number of sub loops for a control loop.
  SUBROUTINE CONTROL_LOOP_NUMBER_SUB_LOOPS_GET(CONTROL_LOOP,NUMBER_OF_SUB_LOOPS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to get the number of sub loops for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_SUB_LOOPS !<On return, the number of sub loops for the specified control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_NUMBER_SUB_LOOPS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_SUB_LOOPS=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_NUMBER_SUB_LOOPS_GET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_NUMBER_SUB_LOOPS_GET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_NUMBER_SUB_LOOPS_GET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_NUMBER_SUB_LOOPS_GET

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the number of sub loops in a control loop.
  SUBROUTINE CONTROL_LOOP_NUMBER_SUB_LOOPS_SET(CONTROL_LOOP,NUMBER_OF_SUB_LOOPS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to set the number of sub loops for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SUB_LOOPS !<The number of sub loops to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_PTR_TYPE), ALLOCATABLE :: OLD_SUB_LOOPS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONTROL_LOOP_NUMBER_SUB_LOOPS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(NUMBER_OF_SUB_LOOPS>=0) THEN
          IF(NUMBER_OF_SUB_LOOPS/=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
            IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
              ALLOCATE(OLD_SUB_LOOPS(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old sub loops.",ERR,ERROR,*999)
              DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                OLD_SUB_LOOPS(loop_idx)%PTR=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
              ENDDO !loop_idx
              DEALLOCATE(CONTROL_LOOP%SUB_LOOPS)
            ENDIF
            ALLOCATE(CONTROL_LOOP%SUB_LOOPS(NUMBER_OF_SUB_LOOPS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate control loop sub loops.",ERR,ERROR,*999)
            IF(NUMBER_OF_SUB_LOOPS>CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
              DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR=>OLD_SUB_LOOPS(loop_idx)%PTR
              ENDDO !loop_idx
              DO loop_idx=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS+1,NUMBER_OF_SUB_LOOPS
                ALLOCATE(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sub loops control loop.",ERR,ERROR,*999)
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%PROBLEM=>CONTROL_LOOP%PROBLEM
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%PARENT_LOOP=>CONTROL_LOOP
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%CONTROL_LOOP_FINISHED=.FALSE.
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%LOOP_TYPE=PROBLEM_CONTROL_SIMPLE_TYPE
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%CONTROL_LOOP_LEVEL=CONTROL_LOOP%CONTROL_LOOP_LEVEL+1
                NULLIFY(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%SIMPLE_LOOP)
                NULLIFY(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%FIXED_LOOP)
                NULLIFY(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%TIME_LOOP)
                NULLIFY(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%WHILE_LOOP)
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%NUMBER_OF_SUB_LOOPS=0
                NULLIFY(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%SOLUTIONS)
                CALL CONTROL_LOOP_SIMPLE_INITIALISE(CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR,ERR,ERROR,*999)
              ENDDO !loop_idx
            ELSE
              DO loop_idx=1,NUMBER_OF_SUB_LOOPS
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR=>OLD_SUB_LOOPS(loop_idx)%PTR
              ENDDO !loop_idx
              DO loop_idx=NUMBER_OF_SUB_LOOPS+1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                CALL CONTROL_LOOP_FINALISE(OLD_SUB_LOOPS(loop_idx)%PTR,ERR,ERROR,*999)
              ENDDO !loop_idx
            ENDIF
            IF(ALLOCATED(OLD_SUB_LOOPS)) DEALLOCATE(OLD_SUB_LOOPS)
            CONTROL_LOOP%NUMBER_OF_SUB_LOOPS=NUMBER_OF_SUB_LOOPS
          ENDIF
        ELSE
          LOCAL_ERROR="The given number of sub loops of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))// &
            & " is invalid. The number of sub loops must be >= 0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_NUMBER_SUB_LOOPS_SET")
    RETURN
999 IF(ALLOCATED(OLD_SUB_LOOPS)) DEALLOCATE(OLD_SUB_LOOPS)
    CALL ERRORS("CONTROL_LOOP_NUMBER_SUB_LOOPS_SET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_NUMBER_SUB_LOOPS_SET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_NUMBER_SUB_LOOPS_SET

  !
  !================================================================================================================================
  !

  !>Finalises a simple control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_SIMPLE_FINALISE(SIMPLE_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_SIMPLE_TYPE), POINTER :: SIMPLE_LOOP !<A pointer to the simple control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_SIMPLE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SIMPLE_LOOP)) THEN
      DEALLOCATE(SIMPLE_LOOP)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_SIMPLE_FINALISE")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_SIMPLE_FINALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_SIMPLE_FINALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SIMPLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a simple loop for a control loop.
  SUBROUTINE CONTROL_LOOP_SIMPLE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to initialise the simple loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONTROL_LOOP_SIMPLE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%SIMPLE_LOOP)) THEN
        CALL FLAG_ERROR("The simple loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%SIMPLE_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate simple loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%SIMPLE_LOOP%CONTROL_LOOP=>CONTROL_LOOP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_SIMPLE_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_SIMPLE_FINALISE(CONTROL_LOOP%SIMPLE_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONTROL_LOOP_SIMPLE_INITIALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_SIMPLE_INITIALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SIMPLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Recursively destroys the solutions for a control loop and all sub control loops.
  RECURSIVE SUBROUTINE CONTROL_LOOP_SOLUTIONS_DESTROY(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to destroy the solutions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2

    CALL ENTERS("CONTROL_LOOP_SOLUTIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      !Destroy the solutions in any sub control loops first
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
        DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
          CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
          CALL CONTROL_LOOP_SOLUTIONS_DESTROY(CONTROL_LOOP2,ERR,ERROR,*999)
        ENDDO !loop_idx
      ENDIF
      !Destroy the solutions in this control loop
      IF(ASSOCIATED(CONTROL_LOOP%SOLUTIONS)) CALL SOLUTIONS_DESTROY(CONTROL_LOOP%SOLUTIONS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_SOLUTIONS_DESTROY")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_SOLUTIONS_DESTROY",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_SOLUTIONS_DESTROY")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SOLUTIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solutions for a control loop.
  RECURSIVE SUBROUTINE CONTROL_LOOP_SOLUTIONS_GET(CONTROL_LOOP,SOLUTIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to get the solutions for.
    TYPE(SOLUTIONS_TYPE), POINTER :: SOLUTIONS !<On exit, a pointer to the control loop solutions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_SOLUTIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLUTIONS)) THEN
        CALL FLAG_ERROR("Solutions is already associated.",ERR,ERROR,*999)
      ELSE
        SOLUTIONS=>CONTROL_LOOP%SOLUTIONS
        IF(.NOT.ASSOCIATED(SOLUTIONS)) CALL FLAG_ERROR("Solutions is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_SOLUTIONS_GET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_SOLUTIONS_GET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_SOLUTIONS_GET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SOLUTIONS_GET

  !
  !================================================================================================================================
  !

  !>Recursively destroys the solver for a control loop and all sub control loops. \todo Create solutions_solver_destory and call?
  RECURSIVE SUBROUTINE CONTROL_LOOP_SOLVER_DESTROY(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to destroy the solver for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx,solution_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
 
    CALL ENTERS("CONTROL_LOOP_SOLVER_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      !Destroy the solvers in any sub control loops first
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
        DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
          CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
          CALL CONTROL_LOOP_SOLVER_DESTROY(CONTROL_LOOP2,ERR,ERROR,*999)
        ENDDO !loop_idx
      ENDIF
      !Destroy the solutions in this control loop
      IF(ASSOCIATED(CONTROL_LOOP%SOLUTIONS)) THEN
        DO solution_idx=1,CONTROL_LOOP%SOLUTIONS%NUMBER_OF_SOLUTIONS
          SOLUTION=>CONTROL_LOOP%SOLUTIONS%SOLUTIONS(solution_idx)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SOLVER=>SOLUTION%SOLVER
            IF(ASSOCIATED(SOLVER)) CALL SOLVER_DESTROY(SOLVER,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !solution_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_SOLVER_DESTROY")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_SOLVER_DESTROY",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_SOLVER_DESTROY")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SOLVER_DESTROY

  !
  !================================================================================================================================
  !

  !>Gets/returns a pointer to the sub loops as specified by the sub loop index for a control loop.
  SUBROUTINE CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,SUB_LOOP_INDEX,SUB_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to get the sub loop for
    INTEGER(INTG), INTENT(IN) :: SUB_LOOP_INDEX !<The sub loop index to get
    TYPE(CONTROL_LOOP_TYPE), POINTER :: SUB_LOOP !<On return, a pointer to the specified sub loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("CONTROL_LOOP_SUB_LOOP_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SUB_LOOP)) THEN
        CALL FLAG_ERROR("Sub loop is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(SUB_LOOP)
        IF(SUB_LOOP_INDEX>0.AND.SUB_LOOP_INDEX<=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
          SUB_LOOP=>CONTROL_LOOP%SUB_LOOPS(SUB_LOOP_INDEX)%PTR
        ELSE
          LOCAL_ERROR="The specified sub loop index of "//TRIM(NUMBER_TO_VSTRING(SUB_LOOP_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The sub loop index must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_SUB_LOOP_GET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_SUB_LOOP_GET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_SUB_LOOP_GET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SUB_LOOP_GET

  !
  !================================================================================================================================
  !

  !>Finalises a time control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_TIME_FINALISE(TIME_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP !<A pointer to the time control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_TIME_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_LOOP)) THEN
      DEALLOCATE(TIME_LOOP)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_TIME_FINALISE")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_TIME_FINALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_TIME_FINALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIME_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a time loop for a control loop.
  SUBROUTINE CONTROL_LOOP_TIME_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to initialise the time loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONTROL_LOOP_TIME_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%TIME_LOOP)) THEN
        CALL FLAG_ERROR("The time loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%TIME_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate time loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%TIME_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%TIME_LOOP%CURRENT_TIME=0.0_DP
        CONTROL_LOOP%TIME_LOOP%START_TIME=0.0_DP
        CONTROL_LOOP%TIME_LOOP%STOP_TIME=1.0_DP
        CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT=0.01_DP        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_TIME_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_TIME_FINALISE(CONTROL_LOOP%TIME_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONTROL_LOOP_TIME_INITIALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_TIME_INITIALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIME_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the time parameters for a time control loop.
  SUBROUTINE CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,START_TIME,STOP_TIME,TIME_INCREMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to set the times for
    REAL(DP), INTENT(IN) :: START_TIME !<The start time for the time control loop.
    REAL(DP), INTENT(IN) :: STOP_TIME !<The stop time for the time control loop.
    REAL(DP), INTENT(IN) :: TIME_INCREMENT !<The time increment for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONTROL_LOOP_TIMES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            IF(ABS(TIME_INCREMENT)<=ZERO_TOLERANCE) THEN
              LOCAL_ERROR="The specified time increment of "//TRIM(NUMBER_TO_VSTRING(TIME_INCREMENT,"*",ERR,ERROR))// &
                & " is invalid. The time increment must not be zero."          
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              IF(TIME_INCREMENT>0.0_DP) THEN
                IF(STOP_TIME<=START_TIME) THEN
                  LOCAL_ERROR="The specified stop time of "//TRIM(NUMBER_TO_VSTRING(STOP_TIME,"*",ERR,ERROR))// &
                    & " is incompatiable with a specified start time of "//TRIM(NUMBER_TO_VSTRING(START_TIME,"*",ERR,ERROR))// &
                    & ". For a positive time increment the stop time must be > than the start time."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                IF(START_TIME<=STOP_TIME) THEN
                  LOCAL_ERROR="The specified start time of "//TRIM(NUMBER_TO_VSTRING(START_TIME,"*",ERR,ERROR))// &
                    & " is incompatiable with a specified stop time of "//TRIM(NUMBER_TO_VSTRING(STOP_TIME,"*",ERR,ERROR))// &
                    & ". For a negative time increment the stop time must be < than the start time."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDIF
            TIME_LOOP%START_TIME=START_TIME
            TIME_LOOP%STOP_TIME=STOP_TIME
            TIME_LOOP%TIME_INCREMENT=TIME_INCREMENT
          ELSE
            CALL FLAG_ERROR("Control loop time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified control loop is not a time control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_TIMES_SET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_TIMES_SET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_TIMES_SET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIMES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the control loop type.
  SUBROUTINE CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,LOOP_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to control loop to set the type of
    INTEGER(INTG), INTENT(IN) :: LOOP_TYPE !<The type of loop type to set \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("CONTROL_LOOP_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FLAG_ERROR("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(LOOP_TYPE/=CONTROL_LOOP%LOOP_TYPE) THEN
          !Initialise the new loop type
          SELECT CASE(LOOP_TYPE)
          CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
            CALL CONTROL_LOOP_SIMPLE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
            CALL CONTROL_LOOP_FIXED_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            CALL CONTROL_LOOP_TIME_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL CONTROL_LOOP_WHILE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The loop type of "//TRIM(NUMBER_TO_VSTRING(LOOP_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Finialise the old loop type
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
            CALL CONTROL_LOOP_SIMPLE_FINALISE(CONTROL_LOOP%SIMPLE_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
            CALL CONTROL_LOOP_FIXED_FINALISE(CONTROL_LOOP%FIXED_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            CALL CONTROL_LOOP_TIME_FINALISE(CONTROL_LOOP%TIME_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL CONTROL_LOOP_WHILE_FINALISE(CONTROL_LOOP%WHILE_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The control loop type of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CONTROL_LOOP%LOOP_TYPE=LOOP_TYPE
        ENDIF
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_TYPE_SET")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_TYPE_SET",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_TYPE_SET")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Finalises a while control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_WHILE_FINALISE(WHILE_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: WHILE_LOOP !<A pointer to the while control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("CONTROL_LOOP_WHILE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(WHILE_LOOP)) THEN
      DEALLOCATE(WHILE_LOOP)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_WHILE_FINALISE")
    RETURN
999 CALL ERRORS("CONTROL_LOOP_WHILE_FINALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_WHILE_FINALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_WHILE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a while loop for a control loop.
  SUBROUTINE CONTROL_LOOP_WHILE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to initialise the while loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("CONTROL_LOOP_WHILE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%WHILE_LOOP)) THEN
        CALL FLAG_ERROR("The while loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%WHILE_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate while loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%WHILE_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=100
        CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("CONTROL_LOOP_WHILE_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_WHILE_FINALISE(CONTROL_LOOP%WHILE_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("CONTROL_LOOP_WHILE_INITIALISE",ERR,ERROR)
    CALL EXITS("CONTROL_LOOP_WHILE_INITIALISE")
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_WHILE_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE CONTROL_LOOP_ROUTINES

