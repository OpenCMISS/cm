!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all problem routines.
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

!> This module handles all problem routines.
MODULE PROBLEM_ROUTINES

  USE BASE_ROUTINES
  USE CLASSICAL_FIELD_ROUTINES
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ELASTICITY_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE SOLVER_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  TYPE(PROBLEMS_TYPE), TARGET :: PROBLEMS
  
  !Interfaces

  INTERFACE PROBLEM_CONTROL_LOOP_GET
    MODULE PROCEDURE PROBLEM_CONTROL_LOOP_GET_0
    MODULE PROCEDURE PROBLEM_CONTROL_LOOP_GET_1
  END INTERFACE !PROBLEM_CONTROL_LOOP_GET

  INTERFACE PROBLEM_DESTROY
    MODULE PROCEDURE PROBLEM_DESTROY_NUMBER
    MODULE PROCEDURE PROBLEM_DESTROY_PTR
  END INTERFACE !PROBLEM_DESTROY
  
  INTERFACE PROBLEM_SPECIFICATION_GET
    MODULE PROCEDURE PROBLEM_SPECIFICATION_GET_NUMBER
    MODULE PROCEDURE PROBLEM_SPECIFICATION_GET_PTR
  END INTERFACE !PROBLEM_SPECIFICATION_GET

  INTERFACE PROBLEM_SPECIFICATION_SET
    MODULE PROCEDURE PROBLEM_SPECIFICATION_SET_NUMBER
    MODULE PROCEDURE PROBLEM_SPECIFICATION_SET_PTR
  END INTERFACE !PROBLEM_SPECIFICATION_SET

  INTERFACE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD
    MODULE PROCEDURE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0
    MODULE PROCEDURE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1
  END INTERFACE !PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD

  INTERFACE PROBLEM_SOLVER_EQUATIONS_GET
    MODULE PROCEDURE PROBLEM_SOLVER_EQUATIONS_GET_0
    MODULE PROCEDURE PROBLEM_SOLVER_EQUATIONS_GET_1
  END INTERFACE !PROBLEM_SOLVER_EQUATIONS_GET

  INTERFACE PROBLEM_SOLVER_GET
    MODULE PROCEDURE PROBLEM_SOLVER_GET_0
    MODULE PROCEDURE PROBLEM_SOLVER_GET_1
  END INTERFACE !PROBLEM_SOLVER_GET
  
  PUBLIC PROBLEMS_INITIALISE,PROBLEMS_FINALISE
  
  PUBLIC PROBLEM_CREATE_START,PROBLEM_CREATE_FINISH,PROBLEM_DESTROY

  PUBLIC PROBLEM_SPECIFICATION_SET

  PUBLIC PROBLEM_CONTROL_LOOP_CREATE_START,PROBLEM_CONTROL_LOOP_CREATE_FINISH,PROBLEM_CONTROL_LOOP_DESTROY,PROBLEM_CONTROL_LOOP_GET
  
  PUBLIC PROBLEM_SOLVER_EQUATIONS_CREATE_START,PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH,PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD

  PUBLIC PROBLEM_SOLVER_JACOBIAN_EVALUATE,PROBLEM_SOLVER_RESIDUAL_EVALUATE
  
  PUBLIC PROBLEM_SOLVERS_CREATE_START,PROBLEM_SOLVERS_CREATE_FINISH,PROBLEM_SOLVERS_DESTROY,PROBLEM_SOLVER_GET
   
  PUBLIC PROBLEM_SOLVE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Solves a problem.
  RECURSIVE SUBROUTINE PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: iteration_idx,loop_idx,solver_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: FIXED_LOOP
    TYPE(CONTROL_LOOP_SIMPLE_TYPE), POINTER :: SIMPLE_LOOP
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: WHILE_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_CONTROL_LOOP_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        !Solve this control loop
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          SIMPLE_LOOP=>CONTROL_LOOP%SIMPLE_LOOP
          IF(ASSOCIATED(SIMPLE_LOOP)) THEN
            IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
              !If there are no sub loops then solve.
              SOLVERS=>CONTROL_LOOP%SOLVERS
              IF(ASSOCIATED(SOLVERS)) THEN
                DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                  SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                  IF(ASSOCIATED(SOLVER)) THEN
                    CALL PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER%SOLVER_EQUATIONS,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !solver_idx
              ELSE
                CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              !If there are sub loops the recursively solve those control loops
              DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
              ENDDO !loop_idx
            ENDIF
          ELSE
            CALL FLAG_ERROR("Control loop simple loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
          FIXED_LOOP=>CONTROL_LOOP%FIXED_LOOP
          IF(ASSOCIATED(FIXED_LOOP)) THEN
            DO iteration_idx=FIXED_LOOP%START_ITERATION,FIXED_LOOP%STOP_ITERATION,FIXED_LOOP%ITERATION_INCREMENT
              FIXED_LOOP%ITERATION_NUMBER=iteration_idx
              IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                !If there are no sub loops then solve
                SOLVERS=>CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(SOLVERS)) THEN
                  DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                    SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                    IF(ASSOCIATED(SOLVER)) THEN
                      CALL PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER%SOLVER_EQUATIONS,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !solver_idx
                ELSE
                  CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !If there are sub loops the recursively solve those control loops
                DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                  CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                  CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                ENDDO !loop_idx
              ENDIF
            ENDDO !iteration_idx
          ELSE
            CALL FLAG_ERROR("Control loop fixed loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            TIME_LOOP%CURRENT_TIME=TIME_LOOP%START_TIME
            TIME_LOOP%ITERATION_NUMBER=0
            DO WHILE(TIME_LOOP%CURRENT_TIME<=TIME_LOOP%STOP_TIME)
              TIME_LOOP%ITERATION_NUMBER=TIME_LOOP%ITERATION_NUMBER+1
              IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                !If there are no sub loops then solve.
                SOLVERS=>CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(SOLVERS)) THEN
                  DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                    SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                    IF(ASSOCIATED(SOLVER)) THEN
                      CALL PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER%SOLVER_EQUATIONS,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !solver_idx
                ELSE
                  CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !If there are sub loops the recursively solve those control loops
                DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                  CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                  CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                ENDDO !loop_idx
              ENDIF
              TIME_LOOP%CURRENT_TIME=TIME_LOOP%CURRENT_TIME+TIME_LOOP%TIME_INCREMENT
            ENDDO !time loop
          ELSE
            CALL FLAG_ERROR("Control loop time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          WHILE_LOOP=>CONTROL_LOOP%WHILE_LOOP
          IF(ASSOCIATED(WHILE_LOOP)) THEN
            WHILE_LOOP%ITERATION_NUMBER=0
            DO WHILE(WHILE_LOOP%CONTINUE_LOOP.AND.WHILE_LOOP%ITERATION_NUMBER<=WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS)
              WHILE_LOOP%ITERATION_NUMBER=WHILE_LOOP%ITERATION_NUMBER+1
              IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
                !If there are no sub loops then solve
                SOLVERS=>CONTROL_LOOP%SOLVERS
                IF(ASSOCIATED(SOLVERS)) THEN
                  DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                    SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
                    IF(ASSOCIATED(SOLVER)) THEN
                      CALL PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER%SOLVER_EQUATIONS,ERR,ERROR,*999)
                    ELSE
                      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !solver_idx
                ELSE
                  CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                !If there are sub loops the recursively solve those control loops
                DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                  CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
                  CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP2,ERR,ERROR,*999)
                ENDDO !loop_idx
              ENDIF
            ENDDO !while loop
          ELSE
            CALL FLAG_ERROR("Control loop while loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The control loop loop type of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Control loop has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_SOLVE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a problem.
  SUBROUTINE PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx

    CALL ENTERS("PROBLEM_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Finish the problem specific setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_INITIAL_TYPE,PROBLEM_SETUP_FINISH_ACTION,ERR,ERROR,*999)
      !Finish the problem creation
      PROBLEM%PROBLEM_FINISHED=.TRUE.
    ELSE        
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of problems = ",PROBLEMS%NUMBER_OF_PROBLEMS,ERR,ERROR,*999)
      DO problem_idx=1,PROBLEMS%NUMBER_OF_PROBLEMS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Problem number  = ",problem_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  User number     = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global number   = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem class   = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%CLASS, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem type    = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%TYPE, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Problem subtype = ",PROBLEMS%PROBLEMS(problem_idx)%PTR%SUBTYPE, &
          & ERR,ERROR,*999)
      ENDDO !problem_idx    
    ENDIF
    
    CALL EXITS("PROBLEM_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("PROBLEM_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE PROBLEM_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a problem defined by USER_NUMBER.
  !>The default values of the PROBLEM attributes are:
  !>- CLASS: 4 (PROBLEM_CLASSICAL_FIELD_CLASS)
  !>- TYPE: 1 (PROBLEM_LAPLACE_EQUATION_TYPE)
  !>- SUBTYPE: 1 (PROBLEM_STANDARD_LAPLACE_SUBTYPE)
  SUBROUTINE PROBLEM_CREATE_START(USER_NUMBER,PROBLEM,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the problem to create
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<On return, a pointer to the created problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx
    TYPE(PROBLEM_TYPE), POINTER :: NEW_PROBLEM
    TYPE(PROBLEM_PTR_TYPE), POINTER :: NEW_PROBLEMS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    NULLIFY(NEW_PROBLEM)
    NULLIFY(NEW_PROBLEMS)

    CALL ENTERS("PROBLEM_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CALL FLAG_ERROR("Problem is already associated.",ERR,ERROR,*999)
    ELSE
      NULLIFY(PROBLEM)
      CALL PROBLEM_USER_NUMBER_FIND(USER_NUMBER,PROBLEM,ERR,ERROR,*999)
      IF(ASSOCIATED(PROBLEM)) THEN
        LOCAL_ERROR="Problem number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))//" has already been created."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        !Allocate the new problem
        ALLOCATE(NEW_PROBLEM,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problem.",ERR,ERROR,*999)
        !Initalise problem
        CALL PROBLEM_INITIALISE(NEW_PROBLEM,ERR,ERROR,*999)
        !Set default problem values
        NEW_PROBLEM%USER_NUMBER=USER_NUMBER
        NEW_PROBLEM%GLOBAL_NUMBER=PROBLEMS%NUMBER_OF_PROBLEMS+1
        NEW_PROBLEM%PROBLEMS=>PROBLEMS
        !Default to a standardised Laplace.
        NEW_PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        NEW_PROBLEM%TYPE=PROBLEM_LAPLACE_EQUATION_TYPE
        NEW_PROBLEM%SUBTYPE=PROBLEM_STANDARD_LAPLACE_SUBTYPE
        NEW_PROBLEM%PROBLEM_FINISHED=.FALSE.
        !Start problem specific setup
        CALL PROBLEM_SETUP(NEW_PROBLEM,PROBLEM_SETUP_INITIAL_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
        !Add new problem into list of problems
        ALLOCATE(NEW_PROBLEMS(PROBLEMS%NUMBER_OF_PROBLEMS+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problems.",ERR,ERROR,*999)
        DO problem_idx=1,PROBLEMS%NUMBER_OF_PROBLEMS
          NEW_PROBLEMS(problem_idx)%PTR=>PROBLEMS%PROBLEMS(problem_idx)%PTR
        ENDDO !problem_idx
        NEW_PROBLEMS(PROBLEMS%NUMBER_OF_PROBLEMS+1)%PTR=>NEW_PROBLEM
        IF(ASSOCIATED(PROBLEMS%PROBLEMS)) DEALLOCATE(PROBLEMS%PROBLEMS)
        PROBLEMS%PROBLEMS=>NEW_PROBLEMS
        PROBLEMS%NUMBER_OF_PROBLEMS=PROBLEMS%NUMBER_OF_PROBLEMS+1
        PROBLEM=>NEW_PROBLEM
      ENDIF
    ENDIF
    
    CALL EXITS("PROBLEM_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CREATE_START")
    RETURN 1   
  END SUBROUTINE PROBLEM_CREATE_START
  
  !
  !================================================================================================================================
  !


  !>Destroys a problem identified by a user number.
  SUBROUTINE PROBLEM_DESTROY_NUMBER(USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the problem to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_position
    LOGICAL :: FOUND
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEMS%PROBLEMS)) THEN
      
      !Find the problem identified by the user number
      FOUND=.FALSE.
      problem_position=0
      DO WHILE(problem_position<PROBLEMS%NUMBER_OF_PROBLEMS.AND..NOT.FOUND)
        problem_position=problem_position+1
        IF(PROBLEMS%PROBLEMS(problem_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
      ENDDO
      
      IF(FOUND) THEN
        
        PROBLEM=>PROBLEMS%PROBLEMS(problem_position)%PTR
        
        !Destroy all the problem components
        CALL PROBLEM_DESTROY(PROBLEM,ERR,ERROR,*999)
        
      ELSE
        LOCAL_ERROR="Problem number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))//" has not been created."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem problems is not associated.",ERR,ERROR,*999)
    ENDIF    

    CALL EXITS("PROBLEM_DESTROY_NUMBER")
    RETURN
999 CALL ERRORS("PROBLEM_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("PROBLEM_DESTROY_NUMBER")
    RETURN 1   
  END SUBROUTINE PROBLEM_DESTROY_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Destroys a problem identified by a pointer.
  SUBROUTINE PROBLEM_DESTROY_PTR(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx,problem_position
    TYPE(PROBLEM_PTR_TYPE), POINTER :: NEW_PROBLEMS(:)

    NULLIFY(NEW_PROBLEMS)

    CALL ENTERS("PROBLEM_DESTROY_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEMS%PROBLEMS)) THEN
        
        problem_position=PROBLEM%GLOBAL_NUMBER
      
        !Destroy all the problem components
        CALL PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*999)
        
        !Remove the problem from the list of problems
        IF(PROBLEMS%NUMBER_OF_PROBLEMS>1) THEN
          ALLOCATE(NEW_PROBLEMS(PROBLEMS%NUMBER_OF_PROBLEMS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new problems.",ERR,ERROR,*999)
          DO problem_idx=1,PROBLEMS%NUMBER_OF_PROBLEMS
            IF(problem_idx<problem_position) THEN
              NEW_PROBLEMS(problem_idx)%PTR=>PROBLEMS%PROBLEMS(problem_idx)%PTR
            ELSE IF(problem_idx>problem_position) THEN
              PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER=PROBLEMS%PROBLEMS(problem_idx)%PTR%GLOBAL_NUMBER-1
              NEW_PROBLEMS(problem_idx-1)%PTR=>PROBLEMS%PROBLEMS(problem_idx)%PTR
            ENDIF
          ENDDO !problem_idx
          DEALLOCATE(PROBLEMS%PROBLEMS)
          PROBLEMS%PROBLEMS=>NEW_PROBLEMS
          PROBLEMS%NUMBER_OF_PROBLEMS=PROBLEMS%NUMBER_OF_PROBLEMS-1
        ELSE
          DEALLOCATE(PROBLEMS%PROBLEMS)
          PROBLEMS%NUMBER_OF_PROBLEMS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Problem problems is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("PROBLEM_DESTROY_PTR")
    RETURN
999 IF(ASSOCIATED(NEW_PROBLEMS)) DEALLOCATE(NEW_PROBLEMS)
998 CALL ERRORS("PROBLEM_DESTROY_PTR",ERR,ERROR)
    CALL EXITS("PROBLEM_DESTROY_PTR")
    RETURN 1   
  END SUBROUTINE PROBLEM_DESTROY_PTR
  
  !
  !================================================================================================================================
  !

  !>Finalise the problem equations add and deallocate all memory.
  SUBROUTINE PROBLEM_EQUATIONS_ADD_FINALISE(EQUATIONS_TO_ADD,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_EQUATIONS_ADD_TYPE), POINTER :: EQUATIONS_TO_ADD !<A pointer to the problem equations add to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_EQUATIONS_ADD_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_TO_ADD)) THEN
      IF(ALLOCATED(EQUATIONS_TO_ADD%CONTROL_LOOP_IDENTIFIER)) DEALLOCATE(EQUATIONS_TO_ADD%CONTROL_LOOP_IDENTIFIER)
      NULLIFY(EQUATIONS_TO_ADD%EQUATIONS_SET_TO_ADD)
      DEALLOCATE(EQUATIONS_TO_ADD)
    ENDIF
       
    CALL EXITS("PROBLEM_EQUATIONS_ADD_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_EQUATIONS_ADD_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_EQUATIONS_ADD_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_EQUATIONS_ADD_FINALISE

 !
  !================================================================================================================================
  !

  !>Initialise the problem equations add and deallocate all memory.
  SUBROUTINE PROBLEM_EQUATIONS_ADD_INITIALISE(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to initialise the equations add for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier of the added equations set
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index of the added equations set
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_EQUATIONS_ADD_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%EQUATIONS_TO_ADD)) THEN
        CALL FLAG_ERROR("Problem equations to add is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM%EQUATIONS_TO_ADD,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations to add.",ERR,ERROR,*999)
        ALLOCATE(PROBLEM%EQUATIONS_TO_ADD%CONTROL_LOOP_IDENTIFIER(SIZE(CONTROL_LOOP_IDENTIFIER,1)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equatiions to add control loop identifier.",ERR,ERROR,*999)
        PROBLEM%EQUATIONS_TO_ADD%CONTROL_LOOP_IDENTIFIER=CONTROL_LOOP_IDENTIFIER
        PROBLEM%EQUATIONS_TO_ADD%SOLVER_INDEX=SOLVER_INDEX
        PROBLEM%EQUATIONS_TO_ADD%EQUATIONS_SET_TO_ADD=>EQUATIONS_SET
        PROBLEM%EQUATIONS_TO_ADD%EQUATIONS_SET_ADDED_INDEX=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_EQUATIONS_ADD_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_EQUATIONS_ADD_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_EQUATIONS_ADD_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_EQUATIONS_ADD_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise the problem and deallocate all memory.
  SUBROUTINE PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) CALL CONTROL_LOOP_DESTROY(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      DEALLOCATE(PROBLEM)
    ENDIF
       
    CALL EXITS("PROBLEM_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a problem.
  SUBROUTINE PROBLEM_INITIALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<The pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The errror code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      PROBLEM%USER_NUMBER=0
      PROBLEM%GLOBAL_NUMBER=0
      PROBLEM%PROBLEM_FINISHED=.FALSE.
      NULLIFY(PROBLEM%PROBLEMS)
      PROBLEM%CLASS=PROBLEM_NO_CLASS
      PROBLEM%TYPE=PROBLEM_NO_TYPE
      PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE
      NULLIFY(PROBLEM%CONTROL_LOOP)
      NULLIFY(PROBLEM%EQUATIONS_TO_ADD)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finish the creation of the control for the problem.
  SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the control for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_CONTROL_LOOP_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN
        IF(PROBLEM%CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
          CALL FLAG_ERROR("Problem control loop has already been finished.",ERR,ERROR,*999)
        ELSE
          !Finish problem specific startup
          CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_CONTROL_TYPE,PROBLEM_SETUP_FINISH_ACTION,ERR,ERROR,*999)
          !Finish problem control creation
          PROBLEM%CONTROL_LOOP%CONTROL_LOOP_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a control loop for a problem.
  !>The default values of the PROBLEM CONTROL LOOP attributes are:
  !>- LOOP_TYPE: PROBLEM_CONTROL_SIMPLE_TYPE
  !>- CONTROL_LOOP_LEVEL: 1
  !>- NUMBER_OF_SUB_LOOPS: 0
  SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of a control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_CONTROL_LOOP_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("The problem control loop is already associated.",ERR,ERROR,*999)        
      ELSE
        !Start the problem specific control setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_CONTROL_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the control loop for a problem.
  SUBROUTINE PROBLEM_CONTROL_LOOP_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_CONTROL_LOOP_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN        
        CALL CONTROL_LOOP_DESTROY(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_DESTROY

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control loop for a problem.
  SUBROUTINE PROBLEM_CONTROL_LOOP_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_CONTROL_LOOP_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_CONTROL_LOOP_GET_1(PROBLEM,(/CONTROL_LOOP_IDENTIFIER/),CONTROL_LOOP,ERR,ERROR,*999) 
       
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_GET_0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control_loop for a problem.
  SUBROUTINE PROBLEM_CONTROL_LOOP_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT
 
    CALL ENTERS("PROBLEM_CONTROL_LOOP_GET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("Control loop is already associated.",ERR,ERROR,*999)
      ELSE
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_1")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_LOOP_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_LOOP_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_LOOP_GET_1
  
  !
  !================================================================================================================================
  !

  !>Sets up the specifices for a problem.
  SUBROUTINE PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The problem setup type \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The problem setup action type \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%CLASS)
      CASE(PROBLEM_ELASTICITY_CLASS)
        CALL ELASTICITY_PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_FLUID_MECHANICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a problem solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,EQUATIONS_SET, &
    & EQUATIONS_SET_INDEX,ERR,ERROR,*)
    
    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to add the equations set to a problem solver
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On return, the index of the equations set that has been added.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0",ERR,ERROR,*999)
    
    CALL PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1(PROBLEM,(/CONTROL_LOOP_IDENTIFIER/),SOLVER_INDEX,EQUATIONS_SET, &
      & EQUATIONS_SET_INDEX,ERR,ERROR,*999)
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_0

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a problem solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,EQUATIONS_SET, &
    & EQUATIONS_SET_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to add the equations set to
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On return, the index of the equations set that has been added.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1",ERR,ERROR,*999)
    
    EQUATIONS_SET_INDEX=0
    IF(ASSOCIATED(PROBLEM)) THEN
      CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
      IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
        NULLIFY(CONTROL_LOOP)
        CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
        SOLVERS=>CONTROL_LOOP%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          IF(SOLVERS%SOLVERS_FINISHED) THEN
            IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
              SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
              IF(ASSOCIATED(SOLVER)) THEN
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
                    CALL PROBLEM_EQUATIONS_ADD_INITIALISE(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,EQUATIONS_SET, &
                      ERR,ERROR,*999)                    
                    !Add equation set via the problem specific setup
                    CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE,PROBLEM_SETUP_DO_ACTION,ERR,ERROR,*999)
                    EQUATIONS_SET_INDEX=PROBLEM%EQUATIONS_TO_ADD%EQUATIONS_SET_ADDED_INDEX
                    CALL PROBLEM_EQUATIONS_ADD_FINALISE(PROBLEM%EQUATIONS_TO_ADD,ERR,ERROR,*999)                   
                  ELSE
                    CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The index must be > 0 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers have not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)          
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD_1

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_SOLVER_EQUATIONS_GET_1(PROBLEM,(/CONTROL_LOOP_IDENTIFIER/),SOLVER_INDEX,SOLVER_EQUATIONS,ERR,ERROR,*999)
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in the solvers to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_GET_1",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        CALL FLAG_ERROR("The solver equations is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(SOLVER_EQUATIONS)
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
          SOLVERS=>CONTROL_LOOP%SOLVERS
          IF(ASSOCIATED(SOLVERS)) THEN            
            IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
              SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
              IF(ASSOCIATED(SOLVER)) THEN
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(.NOT.ASSOCIATED(SOLVER_EQUATIONS)) CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The index must be > 0 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)          
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_1")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_GET_1
  
  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear problem solver.
  SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE(SOLVER,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    
    CALL ENTERS("PROBLEM_SOLVER_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            !Copy the current solution vector to the depenent field
            CALL SOLVER_VARIABLES_STATIC_UPDATE(SOLVER,ERR,ERROR,*999)
            !Make sure the equations sets are up to date
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              !Assemble the equations 
              CALL EQUATIONS_SET_JACOBIAN_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
            ENDDO !equations_set_idx
            !Assemble the solver matrices
            CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_JACOBIAN_ONLY,ERR,ERROR,*999)          
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE
  
  !
  !================================================================================================================================
  !

  !>Evaluates the residual for a nonlinear problem solver.
  SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE(SOLVER,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    
    CALL ENTERS("PROBLEM_SOLVER_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            !Copy the current solution vector to the depenent field
            CALL SOLVER_VARIABLES_STATIC_UPDATE(SOLVER,ERR,ERROR,*999)
            !Make sure the equations sets are up to date
            DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              !Assemble the equations for linear problems
              CALL EQUATIONS_SET_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
            ENDDO !equations_set_idx
            !Assemble the solver matrices
!!TODO: need to work out wether to assemble rhs and residual or residual only
            CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Finish the creation of solvers for a problem.
  SUBROUTINE PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the creation of the solvers for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("PROBLEM_SOLVERS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN              
      !Finish the problem specific solvers setup.
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVERS_TYPE,PROBLEM_SETUP_FINISH_ACTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVERS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVERS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVERS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVERS_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a solvers for the problem
  SUBROUTINE PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to create the solvers for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVERS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN    
      !Start the problem specific solvers setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVERS_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVERS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVERS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVERS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVERS_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Solves a problem.
  SUBROUTINE PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    
    CALL ENTERS("PROBLEM_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          CALL PROBLEM_CONTROL_LOOP_SOLVE(CONTROL_LOOP,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves solver equations for a problem.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        SELECT CASE(SOLVER_EQUATIONS%TIME_DEPENDENCE)
        CASE(SOLVER_EQUATIONS_STATIC)
          SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
          CASE(SOLVER_EQUATIONS_LINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver equations linearity of "//TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(SOLVER_EQUATIONS_QUASISTATIC)
          SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
          CASE(SOLVER_EQUATIONS_LINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver equations linearity of "//TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
          SELECT CASE(SOLVER_EQUATIONS%LINEARITY)
          CASE(SOLVER_EQUATIONS_LINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            CALL PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver equations linearity of "//TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The solver equations time dependence type of "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Solver equations have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves dynamic linear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Make sure the equations sets are up to date
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx          
              !Get current control loop times
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Set the solver time
              CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the next time i.e., current time + time increment
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
              !Back-substitute to find flux values for linear problems
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx
            ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_LINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves dynamic nonlinear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVER)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Apply boundary conditition
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)    
              ENDDO !equations_set_idx          
              !Get current control loop times
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Set the solver time
              CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the next time i.e., current time + time increment
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DYNAMIC_NONLINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves quasistatic linear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
     
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Make sure the equations sets are up to date
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)    
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx          
              !Get current control loop times
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Set the solver time
              CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the next time i.e., current time + time increment
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
              !Back-substitute to find flux values for linear problems
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx
            ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_LINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves quasistatic nonlinear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
   
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVERS=>SOLVER%SOLVERS
        IF(ASSOCIATED(SOLVERS)) THEN
          CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
          IF(ASSOCIATED(CONTROL_LOOP)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
              !Make sure the equations sets are up to date
              DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)    
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx          
              !Get current control loop times
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Set the solver time
              CALL SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
              !Solve for the next time i.e., current time + time increment
              CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
             ELSE
              CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solvers is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_QUASISTATIC_NONLINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Solves static linear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    
    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          !Make sure the equations sets are up to date
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)    
            !Assemble the equations for linear problems
            CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
          ENDDO !equations_set_idx          
          !Solve
          CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
          !Back-substitute to find flux values for linear problems
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,ERR,ERROR,*999)
          ENDDO !equations_set_idx
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_LINEAR_SOLVE
  
  !
  !================================================================================================================================
  !

  !>Solves static nonlinear solver equations.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE(SOLVER_EQUATIONS,ERR,ERROR,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    
    CALL ENTERS("PROBLEM_SOLVER_EQUTIONS_STATIC_NONLINEAR_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER=>SOLVER_EQUATIONS%SOLVER
      IF(ASSOCIATED(SOLVER)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          !Apply boundary conditition
          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)    
          ENDDO !equations_set_idx          
          !Solve
          CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_STATIC_NONLINEAR_SOLVE

  !
  !================================================================================================================================
  !

  !>Destroy the solvers for a problem.
  SUBROUTINE PROBLEM_SOLVERS_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the solvers for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVERS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN        
        CALL CONTROL_LOOP_SOLVERS_DESTROY(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVERS_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVERS_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVERS_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVERS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finish the creation of the solver equations for the problem.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the solver equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN      
      !Finish problem specific startup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE,PROBLEM_SETUP_FINISH_ACTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a problem solver for a problem.
  !>The default values of the SOLVER attributes are:
  !>- SOLVE_TYPE: 1 (SOLVER_LINEAR_TYPE)
  !>- OUTPUT_TYPE: 0 (SOLVER_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (SOLVER_SPARSE_MATRICES)
  SUBROUTINE PROBLEM_SOLVER_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variablesg
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of the solver equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Start the problem specific control setup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the solver equations for a problem.
  SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the solver equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP

    CALL ENTERS("PROBLEM_SOLVER_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
      IF(ASSOCIATED(CONTROL_LOOP)) THEN
        CALL CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY(CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_EQUATIONS_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_EQUATIONS_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop.
  SUBROUTINE PROBLEM_SOLVER_GET_0(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER !<The control loop identifier
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLVER_GET_0",ERR,ERROR,*999)

    CALL PROBLEM_SOLVER_GET_1(PROBLEM,(/CONTROL_LOOP_IDENTIFIER/),SOLVER_INDEX,SOLVER,ERR,ERROR,*999) 
       
    CALL EXITS("PROBLEM_SOLVER_GET_0")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_GET_0",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_GET_0")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_GET_0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop.
  SUBROUTINE PROBLEM_SOLVER_GET_1(PROBLEM,CONTROL_LOOP_IDENTIFIER,SOLVER_INDEX,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: CONTROL_LOOP_IDENTIFIER(:) !<The control loop identifier to get the solver for.
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_SOLVER_GET_1",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        CALL FLAG_ERROR("Solver is already associated.",ERR,ERROR,*998)
      ELSE
        CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP_ROOT)) THEN
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
          SOLVERS=>CONTROL_LOOP%SOLVERS
          IF(ASSOCIATED(SOLVERS)) THEN
            IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
              SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
              IF(.NOT.ASSOCIATED(SOLVER)) CALL FLAG_ERROR("Solvers solver is not associated.",ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
                & " is invalid. The index must be > 0 and <= "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Control loop solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("PROBLEM_SOLVER_GET_1")
    RETURN
999 NULLIFY(SOLVER)
998 CALL ERRORS("PROBLEM_SOLVER_GET_1",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_GET_1")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_GET_1
  
  !
  !================================================================================================================================
  !

  !>Gets the problem specification i.e., problem class, type  and subtype for a problem identified by a user number.
  SUBROUTINE PROBLEM_SPECIFICATION_GET_NUMBER(USER_NUMBER,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the problem to set the specification for.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_CLASS !<On return, the problem class to get.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<On return, the problem equation to get.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<On return, the problem subtype.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("PROBLEM_SPECIFICATION_GET_NUMBER",ERR,ERROR,*999)

    CALL PROBLEM_USER_NUMBER_FIND(USER_NUMBER,PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SPECIFICATION_GET(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
           
    CALL EXITS("PROBLEM_SPECIFICATION_GET_NUMBER")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_GET_NUMBER",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_GET_NUMBER")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_GET_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Gets the problem specification i.e., problem class, type and subtype for a problem identified by a pointer.
  SUBROUTINE PROBLEM_SPECIFICATION_GET_PTR(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_CLASS !<On return, The problem class to set.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<On return, the problem equation type to set.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<On return, the problem subtype to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SPECIFICATION_GET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        PROBLEM_CLASS=PROBLEM%CLASS
        SELECT CASE(PROBLEM_CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_PROBLEM_CLASS_TYPE_GET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_FITTING_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_OPTIMISATION_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM_CLASS,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SPECIFICATION_GET_PTR")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_GET_PTR",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_GET_PTR")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_GET_PTR

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem specification i.e., problem class, type  and subtype for a problem identified by a user number.
  SUBROUTINE PROBLEM_SPECIFICATION_SET_NUMBER(USER_NUMBER,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_CLASS !<The problem class to set.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem equation to set.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM

    CALL ENTERS("PROBLEM_SPECIFICATION_SET_NUMBER",ERR,ERROR,*999)

    CALL PROBLEM_USER_NUMBER_FIND(USER_NUMBER,PROBLEM,ERR,ERROR,*999)
    CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
           
    CALL EXITS("PROBLEM_SPECIFICATION_SET_NUMBER")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_SET_NUMBER",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_SET_NUMBER")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_SET_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the problem specification i.e., problem class, type and subtype for a problem identified by a pointer.
  SUBROUTINE PROBLEM_SPECIFICATION_SET_PTR(PROBLEM,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_CLASS !<The problem class to set.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_EQUATION_TYPE !<The problem equation type to set.
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SPECIFICATION_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        CALL FLAG_ERROR("Problem has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(PROBLEM_CLASS)
        CASE(PROBLEM_ELASTICITY_CLASS)
          CALL ELASTICITY_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_FLUID_MECHANICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_PROBLEM_CLASS_TYPE_SET(PROBLEM,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*999)
        CASE(PROBLEM_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_FITTING_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(PROBLEM_OPTIMISATION_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(PROBLEM_CLASS,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("PROBLEM_SPECIFICATION_SET_PTR")
    RETURN
999 CALL ERRORS("PROBLEM_SPECIFICATION_SET_PTR",ERR,ERROR)
    CALL EXITS("PROBLEM_SPECIFICATION_SET_PTR")
    RETURN 1
  END SUBROUTINE PROBLEM_SPECIFICATION_SET_PTR
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in PROBLEM a pointer to the problem identified by USER_NUMBER. If no problem with that USER_NUMBER exists PROBLEM is left nullified.
  SUBROUTINE PROBLEM_USER_NUMBER_FIND(USER_NUMBER,PROBLEM,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<On return a pointer to the problem with the given user number. If no problem with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: problem_idx

    CALL ENTERS("PROBLEM_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CALL FLAG_ERROR("Problem is already associated.",ERR,ERROR,*999)
    ELSE
      problem_idx=1
      DO WHILE(problem_idx<=PROBLEMS%NUMBER_OF_PROBLEMS.AND..NOT.ASSOCIATED(PROBLEM))
        IF(PROBLEMS%PROBLEMS(problem_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
          PROBLEM=>PROBLEMS%PROBLEMS(problem_idx)%PTR
        ELSE
          problem_idx=problem_idx+1
        ENDIF
      ENDDO
    ENDIF
    
    CALL EXITS("PROBLEM_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("PROBLEM_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("PROBLEM_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE PROBLEM_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises all problems and deallocates all memory.
  SUBROUTINE PROBLEMS_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ::USER_NUMBER

    CALL ENTERS("PROBLEMS_FINALISE",ERR,ERROR,*999)

    DO WHILE(PROBLEMS%NUMBER_OF_PROBLEMS>0)
      USER_NUMBER=PROBLEMS%PROBLEMS(1)%PTR%USER_NUMBER
      CALL PROBLEM_DESTROY(USER_NUMBER,ERR,ERROR,*999)
    ENDDO !problem_idx
    
    CALL EXITS("PROBLEMS_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEMS_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEMS_FINALISE")
    RETURN 1   
  END SUBROUTINE PROBLEMS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all problems.
  SUBROUTINE PROBLEMS_INITIALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEMS_INITIALISE",ERR,ERROR,*999)

    PROBLEMS%NUMBER_OF_PROBLEMS=0
    NULLIFY(PROBLEMS%PROBLEMS)
    
    CALL EXITS("PROBLEMS_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEMS_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEMS_INITIALISE")
    RETURN 1   
  END SUBROUTINE PROBLEMS_INITIALISE
  
  !
  !================================================================================================================================
  !
  
END MODULE PROBLEM_ROUTINES

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the Jacobian for a Newton like nonlinear solver
SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE_PETSC(SNES,X,A,B,FLAG,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The PETSc B Mat
  INTEGER(INTG) :: FLAG !<The PETSC MatStructure flag
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOLVER_VECTOR
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    SOLVER_EQUATIONS=>CTX%SOLVER_EQUATIONS
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
      IF(ASSOCIATED(SOLVER_MATRICES)) THEN
        IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
          SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
          IF(ASSOCIATED(SOLVER_MATRIX)) THEN
            SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
            IF(ASSOCIATED(SOLVER_VECTOR)) THEN
              CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(SOLVER_VECTOR,X,ERR,ERROR,*999)
              
              CALL PROBLEM_SOLVER_JACOBIAN_EVALUATE(CTX,ERR,ERROR,*999)
              
              CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*998)              
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          LOCAL_ERROR="The number of solver matrices of "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
            & " is invalid. There should be 1 solver matrix."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*998)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",ERR,ERROR,*998)
  ENDIF
    
  RETURN
999 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL WRITE_ERROR(ERR,ERROR,*997)
997 CALL FLAG_WARNING("Error evaluating nonlinear Jacobian.",ERR,ERROR,*996)
996 RETURN 
END SUBROUTINE PROBLEM_SOLVER_JACOBIAN_EVALUATE_PETSC
        
!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the residual for a Newton like nonlinear solver
SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC(SNES,X,F,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES type
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec type
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: F !<The PETSc F Vec type
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL_VECTOR,SOLVER_VECTOR
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    SOLVER_EQUATIONS=>CTX%SOLVER_EQUATIONS
    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
      IF(ASSOCIATED(SOLVER_MATRICES)) THEN
        IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
          SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
          IF(ASSOCIATED(SOLVER_MATRIX)) THEN
            SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
            IF(ASSOCIATED(SOLVER_VECTOR)) THEN
              RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
              IF(ASSOCIATED(RESIDUAL_VECTOR)) THEN
                CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(SOLVER_VECTOR,X,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(RESIDUAL_VECTOR,F,ERR,ERROR,*999)                
                
                CALL PROBLEM_SOLVER_RESIDUAL_EVALUATE(CTX,ERR,ERROR,*999)
                
                CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(RESIDUAL_VECTOR,ERR,ERROR,*999)                
              ELSE
                CALL FLAG_ERROR("Residual vector is not associated.",ERR,ERROR,*997)                
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*997)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*997)
          ENDIF
        ELSE
          LOCAL_ERROR="The number of solver matrices of "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
            & " is invalid. There should be 1 solver matrix."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)          
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*997)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",ERR,ERROR,*997)
  ENDIF
  
  RETURN
999 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)  
998 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(RESIDUAL_VECTOR,DUMMY_ERR,DUMMY_ERROR,*997)
997 CALL WRITE_ERROR(ERR,ERROR,*996)
996 CALL FLAG_WARNING("Error evaluating nonlinear residual.",ERR,ERROR,*995)
995 RETURN    
END SUBROUTINE PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC
        
 
