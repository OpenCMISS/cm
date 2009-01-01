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
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLVER_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup PROBLEM_ROUTINES_LinearityTypes PROBLEM_ROUTINES::LinearityTypes
  !> \brief The linearity type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_PROBLEM_LINEARITIES=3 !<The number of problem linearity types defined. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR=1 !<The problemis linear. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR=2 !<The problem is non-linear. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_BCS=3 !<The problem has non-linear boundary conditions. \see PROBLEM_ROUTINES_LinearityTypes,PROBLEM_ROUTINES
  !>@}
  
  !> \addtogroup PROBLEM_ROUTINES_TimeDepedenceTypes PROBLEM_ROUTINES::TimeDepedenceTypes
  !> \brief The time dependence type parameters
  !> \see PROBLEM_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_PROBLEM_TIME_TYPES=3 !<The number of problem time dependence types defined. \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_STATIC=1 !<The problem is static and has no time dependence \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_DYNAMIC=2 !<The problem is dynamic \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_QUASISTATIC=3 !<The problem is quasi-static \see PROBLEM_ROUTINES_TimeDependenceTypes,PROBLEM_ROUTINES
  !>@}

  !Module types

  !Module variables

  TYPE(PROBLEMS_TYPE), TARGET :: PROBLEMS
  
  !Interfaces

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

  PUBLIC PROBLEM_LINEAR,PROBLEM_NONLINEAR,PROBLEM_NONLINEAR_BCS

  PUBLIC PROBLEM_STATIC,PROBLEM_DYNAMIC,PROBLEM_QUASISTATIC

  PUBLIC PROBLEMS_INITIALISE,PROBLEMS_FINALISE
  
  PUBLIC PROBLEM_CREATE_START,PROBLEM_CREATE_FINISH,PROBLEM_DESTROY

  PUBLIC PROBLEM_SPECIFICATION_SET

  PUBLIC PROBLEM_CONTROL_CREATE_START,PROBLEM_CONTROL_CREATE_FINISH,PROBLEM_CONTROL_DESTROY
  
  PUBLIC PROBLEM_SOLUTIONS_CREATE_START,PROBLEM_SOLUTIONS_CREATE_FINISH,PROBLEM_SOLUTION_EQUATIONS_SET_ADD

  PUBLIC PROBLEM_SOLUTION_JACOBIAN_EVALUATE,PROBLEM_SOLUTION_RESIDUAL_EVALUATE
  
  PUBLIC PROBLEM_SOLVER_CREATE_START,PROBLEM_SOLVER_CREATE_FINISH,PROBLEM_SOLVER_DESTROY,PROBLEM_SOLVER_GET
   
  PUBLIC PROBLEM_SOLVE

CONTAINS

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
      !Allocate and set up the problem solutions
      CALL PROBLEM_SOLUTIONS_INITIALISE(PROBLEM,ERR,ERROR,*999)      
      !Finish the problem creation
      PROBLEM%PROBLEM_FINISHED=.TRUE.
    ELSE        
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
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

  !>Finalise the problem and deallocate all memory.
  SUBROUTINE PROBLEM_FINALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      CALL CONTROL_LOOP_FINALISE(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      CALL PROBLEM_SOLUTIONS_FINALISE(PROBLEM,ERR,ERROR,*999)
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
      NULLIFY(PROBLEM%PROBLEMS)
      PROBLEM%CLASS=PROBLEM_NO_CLASS
      PROBLEM%TYPE=PROBLEM_NO_TYPE
      PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE
      NULLIFY(PROBLEM%CONTROL_LOOP)
      PROBLEM%NUMBER_OF_SOLUTIONS=0
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
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
  SUBROUTINE PROBLEM_CONTROL_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the control for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_CONTROL_CREATE_FINISH",ERR,ERROR,*999)

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
      
    CALL EXITS("PROBLEM_CONTROL_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a problem control for a problem.
  SUBROUTINE PROBLEM_CONTROL_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of a control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("PROBLEM_CONTROL_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN
        CALL FLAG_ERROR("The problem control loop is already associated.",ERR,ERROR,*998)        
      ELSE
        !Initialise the problem control loop
        CALL CONTROL_LOOP_INITIALISE(PROBLEM,ERR,ERROR,*999)
        !Start the problem specific control setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_CONTROL_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_CREATE_START")
    RETURN
999 CALL CONTROL_FINALISE(PROBLEM%CONTROL_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("PROBLEM_CONTROL_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the control for a problem.
  SUBROUTINE PROBLEM_CONTROL_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the control for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_CONTROL_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(PROBLEM%CONTROL_LOOP)) THEN        
        CALL CONTROL_LOOP_FINALISE(PROBLEM%CONTROL_LOOP,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem control loop is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_CONTROL_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_CONTROL_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_CONTROL_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_CONTROL_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear problem solution.
  SUBROUTINE PROBLEM_SOLUTION_JACOBIAN_EVALUATE(SOLUTION,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLUTION_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(SOLUTION%SOLUTION_FINISHED) THEN
        SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
          SOLVER=>SOLUTION%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            !Copy the current solution vector to the depenent field
            CALL SOLVER_VARIABLES_UPDATE(SOLVER,ERR,ERROR,*999)
            !Make sure the equations sets are up to date
            DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(EQUATIONS_SET%LINEARITY==EQUATIONS_SET_NONLINEAR) THEN
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_JACOBIAN_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Cannot evaluate the Jacobian for equations set index "// &
                  & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                  & " as it is not a nonlinear equations set."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !equations_set_idx
            !Assemble the solver matrices
            CALL SOLVER_MATRICES_ASSEMBLE(SOLVER,SOLVER_MATRICES_JACOBIAN_ONLY,ERR,ERROR,*999)          
          ELSE
            CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLUTION_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_JACOBIAN_EVALUATE
  
  !
  !================================================================================================================================
  !

  !>Evaluates the residual for a nonlinear problem solution.
  SUBROUTINE PROBLEM_SOLUTION_RESIDUAL_EVALUATE(SOLUTION,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("PROBLEM_SOLUTION_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(SOLUTION%SOLUTION_FINISHED) THEN
        SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
          SOLVER=>SOLUTION%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN            
            !Copy the current solution vector to the depenent field
            CALL SOLVER_VARIABLES_UPDATE(SOLVER,ERR,ERROR,*999)
            !Make sure the equations sets are up to date
            DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
              EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
              IF(SOLUTION%LINEARITY==PROBLEM_SOLUTION_NONLINEAR) THEN
                !Assemble the equations for linear problems
                CALL EQUATIONS_SET_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="Cannot evaluate a residual for equations set index "// &
                  & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                  & " as it is not a nonlinear equations set."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDDO !equations_set_idx
            !Assemble the solver matrices
!!TODO: need to work out wether to assemble rhs and residual or residual only
            CALL SOLVER_MATRICES_ASSEMBLE(SOLVER,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLUTION_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_RESIDUAL_EVALUATE

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

  !>Solves a problem.
  SUBROUTINE PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    
    CALL ENTERS("PROBLEM_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%PROBLEM_FINISHED) THEN
        CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
!!TODO SELECT ON CONTROL TYPES. Branch to separate subroutines.
          DO solution_idx=1,PROBLEM%NUMBER_OF_SOLUTIONS
            SOLUTION=>PROBLEM%SOLUTIONS(solution_idx)%PTR
            IF(ASSOCIATED(SOLUTION)) THEN
              CALL PROBLEM_SOLUTION_SOLVE(SOLUTION,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !solution_idx
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

  !>Solves a solution.
  SUBROUTINE PROBLEM_SOLUTION_SOLVE(SOLUTION,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    
    CALL ENTERS("PROBLEM_SOLUTION_SOLVE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLUTION)) THEN
      IF(SOLUTION%SOLUTION_FINISHED) THEN
        SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN          
          !Make sure the equations sets are up to date
          DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
            EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
            CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*999)	    
            IF(SOLUTION%LINEARITY==PROBLEM_SOLUTION_LINEAR) THEN
              !Assemble the equations for linear problems
              CALL EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*999)
            ENDIF	    
          ENDDO !equations_set_idx          
          SOLVER=>SOLUTION%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            !Solve
            CALL SOLVER_SOLVE(SOLVER,ERR,ERROR,*999)
            !Update depenent field with solution
            CALL SOLVER_VARIABLES_UPDATE(SOLVER,ERR,ERROR,*999)
            !Back-substitute to find flux values for linear problems
            IF(SOLUTION%LINEARITY==PROBLEM_SOLUTION_LINEAR) THEN
              DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                CALL EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,ERR,ERROR,*999)
              ENDDO !equations_set_idx
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solution has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("PROBLEM_SOLUTION_SOLVE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_SOLVE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_SOLVE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_SOLVE

  !
  !================================================================================================================================
  !

  !>Finish the creation of solutions for a problem.
  SUBROUTINE PROBLEM_SOLUTIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    
    CALL ENTERS("PROBLEM_SOLUTIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      DO solution_idx=1,PROBLEM%NUMBER_OF_SOLUTIONS
        SOLUTION=>PROBLEM%SOLUTIONS(solution_idx)%PTR
        IF(ASSOCIATED(SOLUTION)) THEN
          IF(SOLUTION%SOLUTION_FINISHED) THEN
            CALL FLAG_ERROR("Problem solution has been finished",ERR,ERROR,*999)
          ELSE
            !Finish the problem specific solution setup.
            CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLUTION_TYPE,PROBLEM_SETUP_FINISH_ACTION,ERR,ERROR,*999)
            !Finish the problem solution creation
            SOLUTION%SOLUTION_FINISHED=.TRUE.
          ENDIF          
        ELSE
          CALL FLAG_ERROR("The problem solution is not associated",ERR,ERROR,*999)
        ENDIF
      ENDDO !solution_idx
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a solution for the problem. \todo Should this return a pointer to the problem solution???
  SUBROUTINE PROBLEM_SOLUTIONS_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to create a solution for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLUTIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ALLOCATED(PROBLEM%SOLUTIONS)) THEN        
        !Start the problem specific solution setup
        CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLUTION_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("The problem solutions is not allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTIONS_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Adds an equations set to a problem solution.
  SUBROUTINE PROBLEM_SOLUTION_EQUATIONS_SET_ADD(PROBLEM,SOLUTION_INDEX,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to add the equations set to a solution
    INTEGER(INTG), INTENT(IN) :: SOLUTION_INDEX !<The solution index in the problem to add the equations set to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On return, the index of the equations set that has been added.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("PROBLEM_SOLUTION_EQUATIONS_SET_ADD",ERR,ERROR,*999)

    EQUATIONS_SET_INDEX=0
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SOLUTION_INDEX>0.AND.SOLUTION_INDEX<=PROBLEM%NUMBER_OF_SOLUTIONS) THEN
        SOLUTION=>PROBLEM%SOLUTIONS(SOLUTION_INDEX)%PTR
        IF(ASSOCIATED(SOLUTION)) THEN
          IF(SOLUTION%SOLUTION_FINISHED) THEN
            CALL FLAG_ERROR("Solution has already been finished.",ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(EQUATIONS_SET)) THEN
              IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
                IF(ASSOCIATED(SOLUTION%EQUATIONS_SET_TO_ADD)) THEN
                  CALL FLAG_ERROR("Solution equations set to add is already associated.",ERR,ERROR,*999)
                ELSE
                  SOLUTION%EQUATIONS_SET_TO_ADD=>EQUATIONS_SET
                  !Add equation set via the problem specific setup
                  CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLUTION_TYPE,PROBLEM_SETUP_DO_ACTION,ERR,ERROR,*999)
                  NULLIFY(SOLUTION%EQUATIONS_SET_TO_ADD)
                  EQUATIONS_SET_INDEX=SOLUTION%EQUATIONS_SET_ADDED_INDEX
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
            ENDIF            
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified solution index of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_INDEX,"*",ERR,ERROR))// &
          & " is invalid. The index must be > 0 and <= "//TRIM(NUMBER_TO_VSTRING(PROBLEM%NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_EQUATIONS_SET_ADD")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_EQUATIONS_SET_ADD",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_EQUATIONS_SET_ADD")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_EQUATIONS_SET_ADD

  !
  !================================================================================================================================
  !

  !>Finalises a solution and deallocates all memory.
  SUBROUTINE PROBLEM_SOLUTION_FINALISE(SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLUTION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(ASSOCIATED(SOLUTION%SOLUTION_MAPPING)) CALL SOLUTION_MAPPING_DESTROY(SOLUTION%SOLUTION_MAPPING,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLUTION%SOLVER)) CALL SOLVER_DESTROY(SOLUTION%SOLVER,ERR,ERROR,*999)
      DEALLOCATE(SOLUTION)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the solution for a problem.
  SUBROUTINE PROBLEM_SOLUTION_INITIALISE(SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("PROBLEM_SOLUTION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      SOLUTION%SOLUTION_NUMBER=0
      NULLIFY(SOLUTION%PROBLEM)
      SOLUTION%SOLUTION_FINISHED=.FALSE.
      SOLUTION%LINEARITY=PROBLEM_SOLUTION_LINEAR
      NULLIFY(SOLUTION%EQUATIONS_SET_TO_ADD)
      SOLUTION%EQUATIONS_SET_ADDED_INDEX=0
      NULLIFY(SOLUTION%SOLUTION_MAPPING)
      NULLIFY(SOLUTION%SOLVER)
    ELSE
      CALL FLAG_ERROR("Solution is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTION_INITIALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTION_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTION_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTION_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the solutions for a problem and deallocates all memory
  SUBROUTINE PROBLEM_SOLUTIONS_FINALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finalise the solutions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
 
    CALL ENTERS("PROBLEM_SOLUTIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ALLOCATED(PROBLEM%SOLUTIONS)) THEN
        DO solution_idx=1,SIZE(PROBLEM%SOLUTIONS,1)
          CALL PROBLEM_SOLUTION_FINALISE(PROBLEM%SOLUTIONS(solution_idx)%PTR,ERR,ERROR,*999)
        ENDDO ! solution_idx
        DEALLOCATE(PROBLEM%SOLUTIONS)
      ENDIF
      PROBLEM%NUMBER_OF_SOLUTIONS=0
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTIONS_FINALISE")
    RETURN
999 CALL ERRORS("PROBLEM_SOLUTIONS_FINALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTIONS_FINALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the solutions for a problem.
  SUBROUTINE PROBLEM_SOLUTIONS_INITIALISE(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to initialise the solutions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_SOLUTIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ALLOCATED(PROBLEM%SOLUTIONS)) THEN
        CALL FLAG_ERROR("Solutions is already allocated for this problem.",ERR,ERROR,*998)
      ELSE
        IF(PROBLEM%NUMBER_OF_SOLUTIONS>0) THEN
          ALLOCATE(PROBLEM%SOLUTIONS(PROBLEM%NUMBER_OF_SOLUTIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solutions.",ERR,ERROR,*999)
          DO solution_idx=1,PROBLEM%NUMBER_OF_SOLUTIONS
            NULLIFY(PROBLEM%SOLUTIONS(solution_idx)%PTR)
            ALLOCATE(PROBLEM%SOLUTIONS(solution_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution.",ERR,ERROR,*999)
            CALL PROBLEM_SOLUTION_INITIALISE(PROBLEM%SOLUTIONS(solution_idx)%PTR,ERR,ERROR,*999)
            PROBLEM%SOLUTIONS(solution_idx)%PTR%SOLUTION_NUMBER=solution_idx
            PROBLEM%SOLUTIONS(solution_idx)%PTR%PROBLEM=>PROBLEM
          ENDDO !solution_idx
        ELSE
          LOCAL_ERROR="The specified number of problem solutions of "// &
            & TRIM(NUMBER_TO_VSTRING(PROBLEM%NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))// &
            & " is invalid. The number of solutions must be > 0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLUTIONS_INITIALISE")
    RETURN
999 CALL PROBLEM_SOLUTIONS_FINALISE(PROBLEM,ERR,ERROR,*998)
998 CALL ERRORS("PROBLEM_SOLUTIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLUTIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLUTIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of the solver for the problem solutions.
  SUBROUTINE PROBLEM_SOLVER_CREATE_FINISH(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to finish the solvers for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("PROBLEM_SOLVER_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN      
      !Finish problem specific startup
      CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVER_TYPE,PROBLEM_SETUP_FINISH_ACTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("PROBLEM_SOLVER_CREATE_FINISH")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_CREATE_FINISH
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a problem solver for a problem.
  SUBROUTINE PROBLEM_SOLVER_CREATE_START(PROBLEM,ERR,ERROR,*)

    !Argument variablesg
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to start the creation of a solution for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    LOGICAL :: SOLUTIONS_FINISHED

    CALL ENTERS("PROBLEM_SOLVER_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      !Check problem solution has been finished
      IF(ALLOCATED(PROBLEM%SOLUTIONS)) THEN
        SOLUTIONS_FINISHED=.TRUE.
        DO solution_idx=1,PROBLEM%NUMBER_OF_SOLUTIONS
          SOLUTIONS_FINISHED=SOLUTIONS_FINISHED.AND.PROBLEM%SOLUTIONS(solution_idx)%PTR%SOLUTION_FINISHED
        ENDDO !solution_idx
        IF(SOLUTIONS_FINISHED) THEN
          !Start the problem specific control setup
          CALL PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP_SOLVER_TYPE,PROBLEM_SETUP_START_ACTION,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("The problem solutions have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Problem solutions are not allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_CREATE_START")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_CREATE_START",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_CREATE_START")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the solvers for a problem.
  SUBROUTINE PROBLEM_SOLVER_DESTROY(PROBLEM,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to destroy the solver for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solution_idx
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION

    CALL ENTERS("PROBLEM_SOLVER_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ALLOCATED(PROBLEM%SOLUTIONS)) THEN
        DO solution_idx=1,PROBLEM%NUMBER_OF_SOLUTIONS
          SOLUTION=>PROBLEM%SOLUTIONS(solution_idx)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            IF(ASSOCIATED(SOLUTION%SOLVER)) CALL SOLVER_DESTROY(SOLUTION%SOLVER,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !solution_idx
      ELSE
        CALL FLAG_ERROR("Problem solutions is not allocated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("PROBLEM_SOLVER_DESTROY")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_DESTROY",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_DESTROY")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_DESTROY

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a solution on a problme.
  SUBROUTINE PROBLEM_SOLVER_GET(PROBLEM,SOLUTION_INDEX,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: SOLUTION_INDEX !<The solution index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("PROBLEM_SOLVER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SOLUTION_INDEX>0.AND.SOLUTION_INDEX<=PROBLEM%NUMBER_OF_SOLUTIONS) THEN
        SOLUTION=>PROBLEM%SOLUTIONS(SOLUTION_INDEX)%PTR
        IF(ASSOCIATED(SOLUTION)) THEN
          IF(ASSOCIATED(SOLVER)) THEN
            CALL FLAG_ERROR("Solver is already associated.",ERR,ERROR,*999)
          ELSE
            SOLVER=>SOLUTION%SOLVER
            IF(.NOT.ASSOCIATED(SOLVER)) CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The specified solution index of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_INDEX,"*",ERR,ERROR))// &
          & " is invalid. The index must be > 0 and <= "//TRIM(NUMBER_TO_VSTRING(PROBLEM%NUMBER_OF_SOLUTIONS,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF        
       
    CALL EXITS("PROBLEM_SOLVER_GET")
    RETURN
999 CALL ERRORS("PROBLEM_SOLVER_GET",ERR,ERROR)
    CALL EXITS("PROBLEM_SOLVER_GET")
    RETURN 1
  END SUBROUTINE PROBLEM_SOLVER_GET
  
  !
  !================================================================================================================================
  !

  !>Gets the problem specification i.e., problem class, type  and subtype for a problem identified by a user number.
  SUBROUTINE PROBLEM_SPECIFICATION_GET_NUMBER(USER_NUMBER,PROBLEM_CLASS,PROBLEM_EQUATION_TYPE,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the problem to set the specification for.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_CLASS !<The problem class to get.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<The problem equation to get.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<The problem subtype.
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
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_CLASS !<The problem class to set.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_EQUATION_TYPE !<The problem equation type to set.
    INTEGER(INTG), INTENT(OUT) :: PROBLEM_SUBTYPE !<The problem subtype to set.
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
SUBROUTINE PROBLEM_SOLUTION_JACOBIAN_EVALUATE_PETSC(SNES,X,A,B,FLAG,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TYPES
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PETSC_MAT_TYPE), INTENT(INOUT) :: B !<The PETSc B Mat
  INTEGER(INTG) :: FLAG !<The PETSC MatStructure flag
  TYPE(SOLUTION_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOLVER_VECTOR
  TYPE(SOLVER_TYPE), POINTER ::SOLVER
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

    IF(ASSOCIATED(CTX)) THEN
    SOLVER=>CTX%SOLVER
    IF(ASSOCIATED(SOLVER)) THEN
      SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
      IF(ASSOCIATED(SOLVER_MATRICES)) THEN
        IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
          SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
          IF(ASSOCIATED(SOLVER_MATRIX)) THEN
            SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
            IF(ASSOCIATED(SOLVER_VECTOR)) THEN
              CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_ON(SOLVER_VECTOR,X,ERR,ERROR,*999)
              
              CALL PROBLEM_SOLUTION_JACOBIAN_EVALUATE(CTX,ERR,ERROR,*999)

              CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRICES,ERR,ERROR,*999)
              
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
        CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*998)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solution context is not associated.",ERR,ERROR,*998)
  ENDIF
    
  RETURN
999 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL WRITE_ERROR(ERR,ERROR,*997)
997 CALL FLAG_WARNING("Error evaluating nonlinear Jacobian.",ERR,ERROR,*996)
996 RETURN 
END SUBROUTINE PROBLEM_SOLUTION_JACOBIAN_EVALUATE_PETSC
        
!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the residual for a Newton like nonlinear solver
SUBROUTINE PROBLEM_SOLUTION_RESIDUAL_EVALUATE_PETSC(SNES,X,F,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_ROUTINES
  USE STRINGS
  USE TYPES
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES type
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The PETSc X Vec type
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: F !<The PETSc F Vec type
  TYPE(SOLUTION_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  INTEGER(INTG) :: DUMMY_ERR
  TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL_VECTOR,SOLVER_VECTOR
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
  TYPE(VARYING_STRING) :: DUMMY_ERROR,ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    SOLVER=>CTX%SOLVER
    IF(ASSOCIATED(SOLVER)) THEN
      SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
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
                
                CALL PROBLEM_SOLUTION_RESIDUAL_EVALUATE(CTX,ERR,ERROR,*999)
                
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
        CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*997)
    ENDIF
  ELSE
    CALL FLAG_ERROR("Solution context is not associated.",ERR,ERROR,*997)
  ENDIF
  
  RETURN
999 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(SOLVER_VECTOR,DUMMY_ERR,DUMMY_ERROR,*998)  
998 CALL DISTRIBUTED_VECTOR_OVERRIDE_SET_OFF(RESIDUAL_VECTOR,DUMMY_ERR,DUMMY_ERROR,*997)
997 CALL WRITE_ERROR(ERR,ERROR,*996)
996 CALL FLAG_WARNING("Error evaluating nonlinear residual.",ERR,ERROR,*995)
995 RETURN    
END SUBROUTINE PROBLEM_SOLUTION_RESIDUAL_EVALUATE_PETSC
        
 
