!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all the low-level base routines e.g., all debug, control, and low-level communication routines.
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

!> This module contains all the low-level base routines e.g., all debug, control, and low-level communication routines.
MODULE BASE_ROUTINES

  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  USE MACHINE_CONSTANTS 

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(INTG), PARAMETER :: MAX_OUTPUT_LINES=500 !<Maximum number of lines that can be output \see BASE_ROUTINES::WRITE_STR

  !> \addtogroup BASE_ROUTINES_OutputType BASE_ROUTINES::OutputType
  !> \brief Output type parameter
  !> \see BASE_ROUTINES
  !>@{  
  INTEGER(INTG), PARAMETER :: GENERAL_OUTPUT_TYPE=1 !<General output type \see BASE_ROUTINES_OutputType,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: DIAGNOSTIC_OUTPUT_TYPE=2 !<Diagnostic output type \see BASE_ROUTINES_OutputType,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: TIMING_OUTPUT_TYPE=3 !<Timing output type \see BASE_ROUTINES_OutputType,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: ERROR_OUTPUT_TYPE=4 !<Error output type \see BASE_ROUTINES_OutputType,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: WARNING_OUTPUT_TYPE=5 !<Warning output type \see BASE_ROUTINES_OutputType,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: HELP_OUTPUT_TYPE=6 !<Help output type \see BASE_ROUTINES_OutputType,BASE_ROUTINES
  !>@}

  !> \addtogroup BASE_ROUTINES_FileUnits BASE_ROUTINES::FileUnits
  !> \brief File unit parameters
  !> \see BASE_ROUTINES
  !>@{  
  INTEGER(INTG), PARAMETER :: ECHO_FILE_UNIT=10 !<File unit for echo files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: DIAGNOSTICS_FILE_UNIT=11 !<File unit for diagnostic files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: TIMING_FILE_UNIT=12 !<File unit for timing files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: LEARN_FILE_UNIT=13 !<File unit for learn files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IO1_FILE_UNIT=21 !<File unit for general IO 1 files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IO2_FILE_UNIT=22 !<File unit for general IO 2 files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IO3_FILE_UNIT=23 !<File unit for general IO 3 files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IO4_FILE_UNIT=24 !<File unit for general IO 4 files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IO5_FILE_UNIT=25 !<File unit for general IO 5 files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: TEMPORARY_FILE_UNIT=80 !<File unit for temporary files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: OPEN_COMFILE_UNIT=90 !<File unit for open command files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: START_READ_COMFILE_UNIT=90 !<First file unit for read command files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: STOP_READ_COMFILE_UNIT=99 !<Last file unit for read command files \see BASE_ROUTINES_FileUnits,BASE_ROUTINES
  !>@}

  !> \addtogroup BASE_ROUTINES_DiagnosticTypes BASE_ROUTINES::DiagnosticTypes
  !> \brief Diganostic type parameters
  !> \see BASE_ROUTINES
  !>@{  
  INTEGER(INTG), PARAMETER :: ALL_DIAG_TYPE=1 !<Type for setting diagnostic output in all routines \see BASE_ROUTINES_DiagnosticTypes,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IN_DIAG_TYPE=2 !<Type for setting diagnostic output in one routine \see BASE_ROUTINES_DiagnosticTypes,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: FROM_DIAG_TYPE=3 !<Type for setting diagnostic output from one routine downwards \see BASE_ROUTINES_DiagnosticTypes,BASE_ROUTINES
  !>@}

  !> \addtogroup BASE_ROUTINES_TimingTypes BASE_ROUTINES::TimingTypes
  !> \brief Timing type parameters
  !> \see BASE_ROUTINES
  !>@{  
  INTEGER(INTG), PARAMETER :: ALL_TIMING_TYPE=1 !<Type for setting timing output in all routines \see BASE_ROUTINES_TimingTypes,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: IN_TIMING_TYPE=2 !<Type for setting timing output in one routine \see BASE_ROUTINES_TimingTypes,BASE_ROUTINES
  INTEGER(INTG), PARAMETER :: FROM_TIMING_TYPE=3 !<Type for setting timing output from one routine downwards \see BASE_ROUTINES_TimingTypes,BASE_ROUTINES
  !>@}

  !Module types

  !>Contains information for an item in the routine list for diagnostics or timing
  TYPE ROUTINE_LIST_ITEM_TYPE
    TYPE(VARYING_STRING) :: NAME !<Name of the routine
    INTEGER(INTG) :: NUMBER_OF_INVOCATIONS !<Number of times the routine has been invocted
    REAL(SP) :: TOTAL_INCLUSIVE_CPU_TIME !<Total User CPU time spent in the routine inclusive of calls
    REAL(SP) :: TOTAL_INCLUSIVE_SYSTEM_TIME !<Total System CPU time spent in the routine inclusive of calls
    REAL(SP) :: TOTAL_EXCLUSIVE_CPU_TIME !<Total User CPU time spent in the routine exclusive of calls
    REAL(SP) :: TOTAL_EXCLUSIVE_SYSTEM_TIME !<Total System CPU time spent in the routine exclusive of calls
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: NEXT_ROUTINE !<Pointer to the next routine item in the routine list
  END TYPE ROUTINE_LIST_ITEM_TYPE

  !>Contains information for the routine list for diagnostics or timing
  TYPE ROUTINE_LIST_TYPE
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: HEAD !<A pointer to the head of the routine list.
  END TYPE ROUTINE_LIST_TYPE

  !>Contains information for an item in the routine invocation stack
  TYPE ROUTINE_STACK_ITEM_TYPE
    TYPE(VARYING_STRING) :: NAME !<The name of the routine
    REAL(SP) :: INCLUSIVE_CPU_TIME !<User CPU time spent in the routine inclusive of calls 
    REAL(SP) :: INCLUSIVE_SYSTEM_TIME !<System CPU time spent in the routine inclusive of calls 
    REAL(SP) :: EXCLUSIVE_CPU_TIME !<User CPU time spent in the routine exclusive of calls 
    REAL(SP) :: EXCLUSIVE_SYSTEM_TIME !<System CPU time spent in the routine exclusive of calls 
    LOGICAL :: DIAGNOSTICS !<.TRUE. if diagnostics are active in the routine
    LOGICAL :: TIMING !<.TRUE. if timing is active in the routine
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: ROUTINE_LIST_ITEM !<Pointer to the routine list item for diagnostics or timing
    TYPE(ROUTINE_STACK_ITEM_TYPE), POINTER :: PREVIOUS_ROUTINE !<Pointer to the previous routine in the routine stack
  END TYPE ROUTINE_STACK_ITEM_TYPE

  !>Contains information for the routine invocation stack
  TYPE ROUTINE_STACK_TYPE
    TYPE(ROUTINE_STACK_ITEM_TYPE), POINTER :: STACK_POINTER !<Pointer to the top of the stack
  END TYPE ROUTINE_STACK_TYPE

  !Module variables

  LOGICAL, SAVE :: DIAGNOSTICS !<.TRUE. if diagnostic output is required in any routines.
  LOGICAL, SAVE :: DIAGNOSTICS1 !<.TRUE. if level 1 diagnostic output is active in the current routine
  LOGICAL, SAVE :: DIAGNOSTICS2 !<.TRUE. if level 2 diagnostic output is active in the current routine
  LOGICAL, SAVE :: DIAGNOSTICS3 !<.TRUE. if level 3 diagnostic output is active in the current routine
  LOGICAL, SAVE :: DIAGNOSTICS4 !<.TRUE. if level 4 diagnostic output is active in the current routine
  LOGICAL, SAVE :: DIAGNOSTICS5 !<.TRUE. if level 5 diagnostic output is active in the current routine
  LOGICAL, SAVE :: DIAGNOSTICS_LEVEL1 !<.TRUE. if the user has requested level 1 diagnostic output to be active
  LOGICAL, SAVE :: DIAGNOSTICS_LEVEL2 !<.TRUE. if the user has requested level 2 diagnostic output to be active
  LOGICAL, SAVE :: DIAGNOSTICS_LEVEL3 !<.TRUE. if the user has requested level 3 diagnostic output to be active
  LOGICAL, SAVE :: DIAGNOSTICS_LEVEL4 !<.TRUE. if the user has requested level 4 diagnostic output to be active
  LOGICAL, SAVE :: DIAGNOSTICS_LEVEL5 !<.TRUE. if the user has requested level 5 diagnostic output to be active
  LOGICAL, SAVE :: DIAG_ALL_SUBROUTINES !<.TRUE. if diagnostic output is required in all routines
  LOGICAL, SAVE :: DIAG_FROM_SUBROUTINE !<.TRUE. if diagnostic output is required from a particular routine
  LOGICAL, SAVE :: DIAG_FILE_OPEN !<.TRUE. if the diagnostic output file is open
  LOGICAL, SAVE :: DIAG_OR_TIMING !<.TRUE. if diagnostics or time is .TRUE.
  LOGICAL, SAVE :: ECHO_OUTPUT !<.TRUE. if all output is to be echoed to the echo file
  LOGICAL, SAVE :: TIMING !<.TRUE. if timing output is required in any routines.
  LOGICAL, SAVE :: TIMING_SUMMARY !<.TRUE. if timing output will be summary form via a TIMING_SUMMARY_OUTPUT call otherwise timing will be output for routines when the routine exits \see BASE_ROUTINES::TIMING_SUMMARY_OUTPUT
  LOGICAL, SAVE :: TIMING_ALL_SUBROUTINES !<.TRUE. if timing output is required in all routines
  LOGICAL, SAVE :: TIMING_FROM_SUBROUTINE !<.TRUE. if timing output is required from a particular routine
  LOGICAL, SAVE :: TIMING_FILE_OPEN !<.TRUE. if the timing output file is open
  CHARACTER(LEN=MAXSTRLEN), SAVE :: OP_STRING(MAX_OUTPUT_LINES) !<The array of lines to output
  TYPE(ROUTINE_LIST_TYPE), SAVE :: DIAG_ROUTINE_LIST !<The list of routines for which diagnostic output is required
  TYPE(ROUTINE_LIST_TYPE), SAVE :: TIMING_ROUTINE_LIST !<The list of routines for which timing output is required
  TYPE(ROUTINE_STACK_TYPE), SAVE :: ROUTINE_STACK !<The routime invocation stack

  !Interfaces

  INTERFACE

!!!!NOTE: This module needs to call the c cputime function directly in order to avoid a circular module loop when timer uses
!!!!      base_routines.

    SUBROUTINE CPUTIMER(RETURN_TIME, TIME_TYPE, ERR, CERROR)
      !DEC$ ATTRIBUTES C, reference, alias:'_cputimer' :: cputimer
      USE KINDS
      REAL(DP), INTENT(OUT) :: RETURN_TIME
      INTEGER(INTG), INTENT(IN) :: TIME_TYPE
      INTEGER(INTG), INTENT(OUT) :: ERR
      INTEGER(INTG), INTENT(OUT) :: CERROR(*)
    END SUBROUTINE CPUTIMER

  END INTERFACE

  !>Flags an error condition \see BASE_ROUTINES
  INTERFACE FLAG_ERROR
    MODULE PROCEDURE FLAG_ERROR_C 
    MODULE PROCEDURE FLAG_ERROR_VS
  END INTERFACE !FLAG_ERROR
  
  !>Flags a warning to the user \see BASE_ROUTINES
  INTERFACE FLAG_WARNING
    MODULE PROCEDURE FLAG_WARNING_C
    MODULE PROCEDURE FLAG_WARNING_VS
  END INTERFACE !FLAG_WARNING
  
  PUBLIC GENERAL_OUTPUT_TYPE,DIAGNOSTIC_OUTPUT_TYPE,TIMING_OUTPUT_TYPE,ERROR_OUTPUT_TYPE,HELP_OUTPUT_TYPE,DIAGNOSTICS1, &
    & DIAGNOSTICS2,DIAGNOSTICS3,DIAGNOSTICS4,DIAGNOSTICS5,ALL_DIAG_TYPE,IN_DIAG_TYPE,FROM_DIAG_TYPE,OPEN_COMFILE_UNIT, &
    & START_READ_COMFILE_UNIT,STOP_READ_COMFILE_UNIT,TEMPORARY_FILE_UNIT,ALL_TIMING_TYPE,IN_TIMING_TYPE,FROM_TIMING_TYPE, &
    & LEARN_FILE_UNIT

  PUBLIC ENTERS,ERRORS,EXITS,EXTRACT_ERROR_MESSAGE,FLAG_ERROR,FLAG_WARNING,BASE_ROUTINES_FINALISE,BASE_ROUTINES_INITIALISE
  
  PUBLIC OP_STRING,DIAGNOSTICS_SET_ON,DIAGNOSTICS_SET_OFF,OUTPUT_SET_OFF,OUTPUT_SET_ON,TIMING_SET_ON,TIMING_SET_OFF, &
    & TIMING_SUMMARY_OUTPUT,WRITE_ERROR,WRITE_STR

CONTAINS

  !
  !================================================================================================================================
  !

  !>Records the entry into the named procedure and initialises the error code \see BASE_ROUTINES::EXITS
  SUBROUTINE ENTERS(NAME,ERR,ERROR,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: NAME !<The name of the routine being entered
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    REAL(DP) :: ENTERS_CPU_TIME,ENTERS_SYSTEM_TIME
    LOGICAL :: FINISHED
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: LIST_ROUTINE_PTR
    TYPE(ROUTINE_STACK_ITEM_TYPE), POINTER :: NEW_ROUTINE_PTR,ROUTINE_PTR

    IF(DIAG_OR_TIMING) THEN
      !$OMP CRITICAL(ENTERS_1)
      ALLOCATE(NEW_ROUTINE_PTR,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new routine stack item",ERR,ERROR,*999)
      NEW_ROUTINE_PTR%DIAGNOSTICS=.FALSE.
      NEW_ROUTINE_PTR%TIMING=.FALSE.
      NEW_ROUTINE_PTR%NAME=NAME(1:LEN_TRIM(NAME))
      IF(ASSOCIATED(ROUTINE_STACK%STACK_POINTER)) THEN
        NEW_ROUTINE_PTR%PREVIOUS_ROUTINE=>ROUTINE_STACK%STACK_POINTER
        ROUTINE_STACK%STACK_POINTER=>NEW_ROUTINE_PTR
      ELSE
        NULLIFY(NEW_ROUTINE_PTR%PREVIOUS_ROUTINE)
        ROUTINE_STACK%STACK_POINTER=>NEW_ROUTINE_PTR
      ENDIF
      ROUTINE_PTR=>ROUTINE_STACK%STACK_POINTER
      NULLIFY(ROUTINE_PTR%ROUTINE_LIST_ITEM)
      IF(DIAGNOSTICS) THEN
        IF(DIAG_ALL_SUBROUTINES) THEN !turn diagnostics on in all subroutines
          ROUTINE_PTR%DIAGNOSTICS=.TRUE.
        ELSE !diagnostics on in selected subroutines
          FINISHED=.FALSE.
          LIST_ROUTINE_PTR=>DIAG_ROUTINE_LIST%HEAD
          DO WHILE(ASSOCIATED(LIST_ROUTINE_PTR).AND..NOT.FINISHED)
            IF(LIST_ROUTINE_PTR%NAME==ROUTINE_PTR%NAME) THEN
              ROUTINE_PTR%DIAGNOSTICS=.TRUE.
              ROUTINE_PTR%ROUTINE_LIST_ITEM=>LIST_ROUTINE_PTR
              FINISHED=.TRUE.
            ELSE
              LIST_ROUTINE_PTR=>LIST_ROUTINE_PTR%NEXT_ROUTINE
            ENDIF
          ENDDO
          IF(DIAG_FROM_SUBROUTINE) THEN
            IF(ASSOCIATED(ROUTINE_PTR%PREVIOUS_ROUTINE)) THEN
              IF(ROUTINE_PTR%PREVIOUS_ROUTINE%DIAGNOSTICS) ROUTINE_PTR%DIAGNOSTICS=.TRUE.
            ENDIF
          ENDIF
        ENDIF
        IF(ROUTINE_PTR%DIAGNOSTICS) THEN
          DIAGNOSTICS1=DIAGNOSTICS_LEVEL1
          DIAGNOSTICS2=DIAGNOSTICS_LEVEL2
          DIAGNOSTICS3=DIAGNOSTICS_LEVEL3
          DIAGNOSTICS4=DIAGNOSTICS_LEVEL4
          DIAGNOSTICS5=DIAGNOSTICS_LEVEL5
        ELSE
          DIAGNOSTICS1=.FALSE.
          DIAGNOSTICS2=.FALSE.
          DIAGNOSTICS3=.FALSE.
          DIAGNOSTICS4=.FALSE.
          DIAGNOSTICS5=.FALSE.
        ENDIF
        IF(ROUTINE_PTR%DIAGNOSTICS) THEN
          WRITE(OP_STRING,'("*** Enters: ",A)') NAME(1:LEN_TRIM(NAME))
          CALL WRITE_STR(DIAGNOSTIC_OUTPUT_TYPE,ERR,ERROR,*999)
        ELSE IF(ASSOCIATED(ROUTINE_PTR%PREVIOUS_ROUTINE)) THEN
          !CPB 16/05/2007 Only show the calls if we have level 3 diagnostics or higher
          IF(DIAGNOSTICS3) THEN
            IF(ROUTINE_PTR%PREVIOUS_ROUTINE%DIAGNOSTICS) THEN
              WRITE(OP_STRING,'("*** Calls : ",A)') NAME(1:LEN_TRIM(NAME))
              CALL WRITE_STR(DIAGNOSTIC_OUTPUT_TYPE,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(TIMING) THEN
        CALL CPUTIMER(ENTERS_CPU_TIME,1,ERR,CERROR)
        CALL CPUTIMER(ENTERS_SYSTEM_TIME,2,ERR,CERROR)
        ROUTINE_PTR%INCLUSIVE_CPU_TIME=REAL(ENTERS_CPU_TIME,SP)
        ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME=REAL(ENTERS_SYSTEM_TIME,SP)
        ROUTINE_PTR%EXCLUSIVE_CPU_TIME=0.0_SP
        ROUTINE_PTR%EXCLUSIVE_SYSTEM_TIME=0.0_SP
        IF(TIMING_ALL_SUBROUTINES) THEN
          ROUTINE_PTR%TIMING=.TRUE.
        ELSE
          FINISHED=.FALSE.
          LIST_ROUTINE_PTR=>TIMING_ROUTINE_LIST%HEAD
          DO WHILE(ASSOCIATED(LIST_ROUTINE_PTR).AND..NOT.FINISHED)
            IF(LIST_ROUTINE_PTR%NAME==ROUTINE_PTR%NAME) THEN
              ROUTINE_PTR%TIMING=.TRUE.
              ROUTINE_PTR%ROUTINE_LIST_ITEM=>LIST_ROUTINE_PTR
              FINISHED=.TRUE.
            ELSE
              LIST_ROUTINE_PTR=>LIST_ROUTINE_PTR%NEXT_ROUTINE
            ENDIF
          ENDDO
          IF(TIMING_FROM_SUBROUTINE) THEN
            IF(ASSOCIATED(ROUTINE_PTR%PREVIOUS_ROUTINE)) THEN
              IF(ROUTINE_PTR%PREVIOUS_ROUTINE%TIMING) ROUTINE_PTR%TIMING=.TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      !$OMP END CRITICAL(ENTERS_1)
    ENDIF

    RETURN
999 RETURN 1
  END SUBROUTINE ENTERS

  !
  !================================================================================================================================
  !

  !>Records the exiting error of the subroutine 
  SUBROUTINE ERRORS(NAME,ERR,ERROR)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: NAME !<The name of the routine with an error condition
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    IF(ERR==0) ERR=1
    !CPB 20/02/07 aix compiler does not like varying strings so split the concatenate statement up into two statements
    LOCAL_ERROR=ERROR//ERROR_SEPARATOR_CONSTANT
    ERROR=LOCAL_ERROR//NAME(1:LEN_TRIM(NAME))

    RETURN
  END SUBROUTINE ERRORS

  !
  !================================================================================================================================
  !

  !>Records the exit out of the named procedure \see BASE_ROUTINES::ENTERS
  SUBROUTINE EXITS(NAME)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: NAME !<The name of the routine exiting
    !Local variables
    INTEGER(INTG) :: CERROR(100),ERR
    REAL(DP) :: EXITS_CPU_TIME,EXITS_SYSTEM_TIME
    TYPE(VARYING_STRING) :: ERROR
    TYPE(ROUTINE_STACK_ITEM_TYPE), POINTER :: PREVIOUS_ROUTINE_PTR,ROUTINE_PTR

    IF(DIAG_OR_TIMING) THEN
      !$OMP CRITICAL(EXITS_1)
      ROUTINE_PTR=>ROUTINE_STACK%STACK_POINTER
      IF(ASSOCIATED(ROUTINE_PTR)) THEN
        PREVIOUS_ROUTINE_PTR=>ROUTINE_PTR%PREVIOUS_ROUTINE
        IF(DIAGNOSTICS) THEN
          IF(ROUTINE_PTR%DIAGNOSTICS) THEN
            WRITE(OP_STRING,'("*** Exits : ",A)') NAME(1:LEN_TRIM(NAME))
            CALL WRITE_STR(DIAGNOSTIC_OUTPUT_TYPE,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(PREVIOUS_ROUTINE_PTR)) THEN
            IF(PREVIOUS_ROUTINE_PTR%DIAGNOSTICS) THEN
              DIAGNOSTICS1=DIAGNOSTICS_LEVEL1
              DIAGNOSTICS2=DIAGNOSTICS_LEVEL2
              DIAGNOSTICS3=DIAGNOSTICS_LEVEL3
              DIAGNOSTICS4=DIAGNOSTICS_LEVEL4
              DIAGNOSTICS5=DIAGNOSTICS_LEVEL5
            ELSE
              DIAGNOSTICS1=.FALSE.
              DIAGNOSTICS2=.FALSE.
              DIAGNOSTICS3=.FALSE.
              DIAGNOSTICS4=.FALSE.
              DIAGNOSTICS5=.FALSE.
            ENDIF
          ENDIF
        ENDIF

        IF(TIMING) THEN
          CALL CPUTIMER(EXITS_CPU_TIME,1,ERR,CERROR)
          CALL CPUTIMER(EXITS_SYSTEM_TIME,2,ERR,CERROR)
          ROUTINE_PTR%INCLUSIVE_CPU_TIME=ABS(REAL(EXITS_CPU_TIME,SP)-ROUTINE_PTR%INCLUSIVE_CPU_TIME)
          ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME=ABS(REAL(EXITS_SYSTEM_TIME,SP)-ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME)
          IF(ASSOCIATED(PREVIOUS_ROUTINE_PTR)) THEN
            PREVIOUS_ROUTINE_PTR%EXCLUSIVE_CPU_TIME=PREVIOUS_ROUTINE_PTR%EXCLUSIVE_CPU_TIME+ROUTINE_PTR%INCLUSIVE_CPU_TIME
            PREVIOUS_ROUTINE_PTR%EXCLUSIVE_SYSTEM_TIME=PREVIOUS_ROUTINE_PTR%EXCLUSIVE_SYSTEM_TIME+ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME
          ENDIF
          IF(ASSOCIATED(ROUTINE_PTR%ROUTINE_LIST_ITEM)) THEN
            ROUTINE_PTR%ROUTINE_LIST_ITEM%NUMBER_OF_INVOCATIONS=ROUTINE_PTR%ROUTINE_LIST_ITEM%NUMBER_OF_INVOCATIONS+1
            ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_INCLUSIVE_CPU_TIME=ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_INCLUSIVE_CPU_TIME+ &
              & ROUTINE_PTR%INCLUSIVE_CPU_TIME
            ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_INCLUSIVE_SYSTEM_TIME=ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_INCLUSIVE_SYSTEM_TIME+ &
              & ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME
            IF(ASSOCIATED(PREVIOUS_ROUTINE_PTR)) THEN
              IF(ASSOCIATED(PREVIOUS_ROUTINE_PTR%ROUTINE_LIST_ITEM)) THEN
                PREVIOUS_ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_EXCLUSIVE_CPU_TIME=PREVIOUS_ROUTINE_PTR%ROUTINE_LIST_ITEM% &
                  & TOTAL_EXCLUSIVE_CPU_TIME+PREVIOUS_ROUTINE_PTR%EXCLUSIVE_CPU_TIME
                PREVIOUS_ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_EXCLUSIVE_SYSTEM_TIME=PREVIOUS_ROUTINE_PTR%ROUTINE_LIST_ITEM% &
                  & TOTAL_EXCLUSIVE_SYSTEM_TIME+PREVIOUS_ROUTINE_PTR%EXCLUSIVE_SYSTEM_TIME
              ENDIF
            ENDIF
          ENDIF
          IF(ROUTINE_PTR%TIMING) THEN
            IF(.NOT.TIMING_SUMMARY) THEN
              WRITE(OP_STRING,'("*** Timing : ",A)') NAME(1:LEN_TRIM(NAME))
              CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
              IF(ASSOCIATED(ROUTINE_PTR%ROUTINE_LIST_ITEM)) THEN
                WRITE(OP_STRING,'("***    Number of invocations: ",I10)') ROUTINE_PTR%ROUTINE_LIST_ITEM%NUMBER_OF_INVOCATIONS
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
                WRITE(OP_STRING,'("***    Routine times:  Call Inclusive   Call Exclusive   Total Inclusive   Average Inclusive")')
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
                WRITE(OP_STRING,'("***    CPU       (s):  ",E14.6,"   ",E14.6,"   ",E15.6,"   ",E17.6)')  &
                  & ROUTINE_PTR%INCLUSIVE_CPU_TIME,ROUTINE_PTR%INCLUSIVE_CPU_TIME-ROUTINE_PTR%EXCLUSIVE_CPU_TIME, &
                  & ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_INCLUSIVE_CPU_TIME,ROUTINE_PTR%ROUTINE_LIST_ITEM% &
                  & TOTAL_INCLUSIVE_CPU_TIME/REAL(ROUTINE_PTR%ROUTINE_LIST_ITEM%NUMBER_OF_INVOCATIONS,SP)
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
                WRITE(OP_STRING,'("***    System    (s):  ",E14.6,"   ",E14.6,"   ",E15.6,"   ",E17.6)')  &
                  & ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME,ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME-ROUTINE_PTR%EXCLUSIVE_SYSTEM_TIME, &
                  & ROUTINE_PTR%ROUTINE_LIST_ITEM%TOTAL_INCLUSIVE_SYSTEM_TIME,ROUTINE_PTR%ROUTINE_LIST_ITEM% &
                  & TOTAL_INCLUSIVE_SYSTEM_TIME/REAL(ROUTINE_PTR%ROUTINE_LIST_ITEM%NUMBER_OF_INVOCATIONS,SP)
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
              ELSE
                WRITE(OP_STRING,'("***    Routine times:  Call Inclusive   Call Exclusive")')
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
                WRITE(OP_STRING,'("***    CPU       (s):  ",E14.6,"   ",E14.6)')  &
                  & ROUTINE_PTR%INCLUSIVE_CPU_TIME,ROUTINE_PTR%INCLUSIVE_CPU_TIME-ROUTINE_PTR%EXCLUSIVE_CPU_TIME
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
                WRITE(OP_STRING,'("***    System    (s):  ",E14.6,"   ",E14.6)')  &
                  & ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME,ROUTINE_PTR%INCLUSIVE_SYSTEM_TIME-ROUTINE_PTR%EXCLUSIVE_SYSTEM_TIME
                CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        IF(ASSOCIATED(PREVIOUS_ROUTINE_PTR)) THEN
          ROUTINE_STACK%STACK_POINTER=>PREVIOUS_ROUTINE_PTR
        ELSE
          NULLIFY(ROUTINE_STACK%STACK_POINTER)
        ENDIF
        DEALLOCATE(ROUTINE_PTR)

        !ELSE ERROR????
      ENDIF
      !$OMP END CRITICAL(EXITS_1)
    ENDIF

999 RETURN
  END SUBROUTINE EXITS

  !
  !================================================================================================================================
  !

  !>Extracts the error message from a CMISS error string and returns it as varying string
  FUNCTION EXTRACT_ERROR_MESSAGE(ERR,ERROR)

    !Argument variables    
    INTEGER(INTG), INTENT(IN) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(IN) :: ERROR !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: EXTRACT_ERROR_MESSAGE !<On exit, the error message contained in the error string
    !Local Variables
    INTEGER(INTG) :: POSITION
    
    POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
    EXTRACT_ERROR_MESSAGE=EXTRACT(ERROR,1,POSITION-1)
    
    RETURN    
  END FUNCTION EXTRACT_ERROR_MESSAGE
    
  !
  !================================================================================================================================
  !
 
  !>Sets the error string specified by a character string and flags an error 
  SUBROUTINE FLAG_ERROR_C(STRING,ERR,ERROR,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The error condition string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: STRING_LENGTH

    IF(ERR==0) ERR=1
    STRING_LENGTH=LEN_TRIM(STRING)
    ERROR=STRING(1:STRING_LENGTH)

    RETURN 1
  END SUBROUTINE FLAG_ERROR_C

  !
  !================================================================================================================================
  !

  !>Sets the error string specified by a varying string and flags an error.
  SUBROUTINE FLAG_ERROR_VS(STRING,ERR,ERROR,*)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The error condition string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    IF(ERR==0) ERR=1
    ERROR=STRING

    RETURN 1
  END SUBROUTINE FLAG_ERROR_VS

  !
  !================================================================================================================================
  !

  !>Writes a warning message specified by a character string to the user.
  SUBROUTINE FLAG_WARNING_C(STRING,ERR,ERROR,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The warning string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    WRITE(OP_STRING,'(">>WARNING: ",A)') STRING
    CALL WRITE_STR(WARNING_OUTPUT_TYPE,ERR,ERROR,*999)

    RETURN 
999 CALL ERRORS("FLAG_WARNING_C",ERR,ERROR)
    RETURN 1
  END SUBROUTINE FLAG_WARNING_C

  !
  !================================================================================================================================
  !

  !>Writes a warning message specified by a varying string to the user.
  SUBROUTINE FLAG_WARNING_VS(STRING,ERR,ERROR,*)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The warning string
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    WRITE(OP_STRING,'(">>WARNING: ",A)') CHAR(STRING)
    CALL WRITE_STR(WARNING_OUTPUT_TYPE,ERR,ERROR,*999)

    RETURN 
999 CALL ERRORS("FLAG_WARNING_VS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE FLAG_WARNING_VS

  !
  !================================================================================================================================
  !

  !>Finalises the base_routines module and deallocates all memory. \todo Finish this routine and deallocate memory.
  SUBROUTINE BASE_ROUTINES_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    ERR=0
    ERROR=""
    
    RETURN 
999 RETURN 1
  END SUBROUTINE BASE_ROUTINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the variables required for the base_routines module.
  SUBROUTINE BASE_ROUTINES_INITIALISE(ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j

    ERR=0
    ERROR=""
    DIAGNOSTICS=.FALSE.
    DIAGNOSTICS1=.FALSE.
    DIAGNOSTICS2=.FALSE.
    DIAGNOSTICS3=.FALSE.
    DIAGNOSTICS4=.FALSE.
    DIAGNOSTICS5=.FALSE.
    DIAGNOSTICS_LEVEL1=.FALSE.
    DIAGNOSTICS_LEVEL2=.FALSE.
    DIAGNOSTICS_LEVEL3=.FALSE.
    DIAGNOSTICS_LEVEL4=.FALSE.
    DIAGNOSTICS_LEVEL5=.FALSE.
    DIAG_ALL_SUBROUTINES=.TRUE.
    DIAG_FROM_SUBROUTINE=.FALSE.
    DIAG_FILE_OPEN=.FALSE.
    DIAG_OR_TIMING=.FALSE.
    ECHO_OUTPUT=.FALSE.
    TIMING=.FALSE.
    TIMING_SUMMARY=.FALSE.
    TIMING_ALL_SUBROUTINES=.TRUE.
    TIMING_FROM_SUBROUTINE=.FALSE.
    TIMING_FILE_OPEN=.FALSE.
    !Initialise loose tolerance here rather than in constants.f90
    LOOSE_TOLERANCE=SQRT(EPSILON(1.0_DP))
    LOOSE_TOLERANCE_SP=SQRT(EPSILON(1.0_SP))

    !Initialise OP_STRING
    SELECT CASE(MACHINE_OS)
    CASE(VMS_OS)
      DO i=1,MAX_OUTPUT_LINES
        OP_STRING(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE(IRIX_OS,LINUX_OS,AIX_OS)
      DO i=1,MAX_OUTPUT_LINES
        DO j=1,MAXSTRLEN
          OP_STRING(i)(j:j)=' '
        ENDDO !j
      ENDDO !i
    CASE(WINDOWS_OS)
      DO i=1,MAX_OUTPUT_LINES
        OP_STRING(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE DEFAULT
      CALL FLAG_ERROR("Operating system not implemented",ERR,ERROR,*999)
    END SELECT

    !Initialise diagnostics and tracing
    NULLIFY(ROUTINE_STACK%STACK_POINTER)
    NULLIFY(DIAG_ROUTINE_LIST%HEAD)
    NULLIFY(TIMING_ROUTINE_LIST%HEAD)

    RETURN 
999 RETURN 1
  END SUBROUTINE BASE_ROUTINES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets diagnositics off \see BASE_ROUTINES::DIAGNOSTICS_SET_ON
  SUBROUTINE DIAGNOSTICS_SET_OFF(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: NEXT_ROUTINE,ROUTINE

    CALL ENTERS("DIAGNOSTICS_SET_OFF",ERR,ERROR,*999)

    IF(DIAGNOSTICS) THEN
      IF(DIAG_FILE_OPEN) THEN
        DIAG_FILE_OPEN=.FALSE.
        CLOSE(UNIT=DIAGNOSTICS_FILE_UNIT)
      ENDIF
      IF(DIAG_ALL_SUBROUTINES) THEN
        DIAG_ALL_SUBROUTINES=.FALSE.
      ELSE
        ROUTINE=>DIAG_ROUTINE_LIST%HEAD
        DO WHILE(ASSOCIATED(ROUTINE))
          NEXT_ROUTINE=>ROUTINE%NEXT_ROUTINE
          DEALLOCATE(ROUTINE)
          ROUTINE=>NEXT_ROUTINE
        ENDDO
        NULLIFY(DIAG_ROUTINE_LIST%HEAD)
        DIAG_FROM_SUBROUTINE=.FALSE.
      ENDIF
      DIAGNOSTICS_LEVEL1=.FALSE.
      DIAGNOSTICS_LEVEL2=.FALSE.
      DIAGNOSTICS_LEVEL3=.FALSE.
      DIAGNOSTICS_LEVEL4=.FALSE.
      DIAGNOSTICS_LEVEL5=.FALSE.
      DIAGNOSTICS1=.FALSE.
      DIAGNOSTICS2=.FALSE.
      DIAGNOSTICS3=.FALSE.
      DIAGNOSTICS4=.FALSE.
      DIAGNOSTICS5=.FALSE.
      DIAGNOSTICS=.FALSE.
      DIAG_OR_TIMING=TIMING
    ELSE
      CALL FLAG_ERROR("Diagnositics is not on.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIAGNOSTICS_SET_OFF")
    RETURN
999 CALL ERRORS("DIAGNOSTICS_SET_OFF",ERR,ERROR)
    CALL EXITS("DIAGNOSTICS_SET_OFF")
    RETURN 1
  END SUBROUTINE DIAGNOSTICS_SET_OFF

  !
  !================================================================================================================================
  !

  !>Sets diagnositics on \see BASE_ROUTINES::DIAGNOSTICS_SET_OFF
  SUBROUTINE DIAGNOSTICS_SET_ON(DIAG_TYPE,LEVEL_LIST,DIAG_FILENAME,ROUTINE_LIST,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: DIAG_TYPE !<The type of diagnostics to set on \see BASE_ROUTINES_DiagnosticTypes
    INTEGER(INTG), INTENT(IN) :: LEVEL_LIST(:) !<The list of diagnostic levels to set on
    CHARACTER(LEN=*), INTENT(IN) :: DIAG_FILENAME !<If present the name of the file to output diagnostic information to. If omitted the diagnostic output is sent to the screen
    CHARACTER(LEN=*), INTENT(IN) :: ROUTINE_LIST(:) !<The list of routines to set diagnostics on in.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,LEVEL
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("A")-ICHAR("a"))
    CHARACTER(LEN=1) :: CHARAC1, CHARAC2
    CHARACTER(LEN=MAXSTRLEN) :: FILENAME
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: NEXT_ROUTINE,PREVIOUS_ROUTINE,ROUTINE
    TYPE(VARYING_STRING) :: LOCAL_STRING

    NULLIFY(ROUTINE)
    
    CALL ENTERS("DIAGNOSTICS_SET_ON",ERR,ERROR,*999)

    IF(LEN_TRIM(DIAG_FILENAME)>=1) THEN
      IF(DIAG_FILE_OPEN) CLOSE(UNIT=DIAGNOSTICS_FILE_UNIT)
      FILENAME=DIAG_FILENAME(1:LEN_TRIM(DIAG_FILENAME))//".diag"
      OPEN(UNIT=DIAGNOSTICS_FILE_UNIT,FILE=FILENAME(1:LEN_TRIM(FILENAME)),STATUS="UNKNOWN",IOSTAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not open diagnostics file",ERR,ERROR,*999)
      DIAG_FILE_OPEN=.TRUE.
    ENDIF
    SELECT CASE(DIAG_TYPE)
    CASE(ALL_DIAG_TYPE)
      DIAG_ALL_SUBROUTINES=.TRUE.
    CASE(IN_DIAG_TYPE,FROM_DIAG_TYPE)
      DIAG_ALL_SUBROUTINES=.FALSE.
      IF(ASSOCIATED(DIAG_ROUTINE_LIST%HEAD)) THEN
        ROUTINE=>DIAG_ROUTINE_LIST%HEAD
        DO WHILE(ASSOCIATED(ROUTINE))
          NEXT_ROUTINE=>ROUTINE%NEXT_ROUTINE
          DEALLOCATE(ROUTINE)
          ROUTINE=>NEXT_ROUTINE
        ENDDO
        NULLIFY(DIAG_ROUTINE_LIST%HEAD)
      ENDIF
      ALLOCATE(ROUTINE,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate routine list item",ERR,ERROR,*999)
      ROUTINE%NAME=""
      DO i=1,LEN(ROUTINE_LIST(1))
        CHARAC1=ROUTINE_LIST(1)(i:i)
        IF(LGE(CHARAC1,"a").AND.LLE(CHARAC1,"z")) THEN
          CHARAC2=CHAR(ICHAR(CHARAC1)+OFFSET)
        ELSE
          CHARAC2=CHARAC1
        ENDIF
        !CPB 26/9/05 Aix compiler doesn't like the vstring so split the statement and put a char() around it
        !ROUTINE%NAME=ROUTINE%NAME//CHARAC2
        LOCAL_STRING=ROUTINE%NAME//CHARAC2
        ROUTINE%NAME=CHAR(LOCAL_STRING)
      ENDDO !i
      PREVIOUS_ROUTINE=>ROUTINE
      NULLIFY(ROUTINE%NEXT_ROUTINE)
      DIAG_ROUTINE_LIST%HEAD=>ROUTINE
      DO i=2,SIZE(ROUTINE_LIST,1)
        ALLOCATE(ROUTINE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate routine list item",ERR,ERROR,*999)
        ROUTINE%NAME=""
        DO j=1,LEN(ROUTINE_LIST(i))
          CHARAC1=ROUTINE_LIST(i)(j:j)
          IF(LGE(CHARAC1,"a").AND.LLE(CHARAC1,"z")) THEN
            CHARAC2=CHAR(ICHAR(CHARAC1)+OFFSET)
          ELSE
            CHARAC2=CHARAC1
          ENDIF
          !CPB 26/9/05 Aix compiler doesn't like the vstring so split the statement and put a char() around it
          !ROUTINE%NAME=ROUTINE%NAME//CHARAC2
          LOCAL_STRING=ROUTINE%NAME//CHARAC2
          ROUTINE%NAME=CHAR(LOCAL_STRING)
        ENDDO !i
        NULLIFY(ROUTINE%NEXT_ROUTINE)
        PREVIOUS_ROUTINE%NEXT_ROUTINE=>ROUTINE
        PREVIOUS_ROUTINE=>ROUTINE
      ENDDO !i
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid diagnostic type",ERR,ERROR,*999)
    END SELECT
    DO i=1,SIZE(LEVEL_LIST,1)
      LEVEL=LEVEL_LIST(i)
      SELECT CASE(LEVEL)
      CASE(1)
        DIAGNOSTICS_LEVEL1=.TRUE.
      CASE(2)
        DIAGNOSTICS_LEVEL2=.TRUE.
      CASE(3)
        DIAGNOSTICS_LEVEL3=.TRUE.
      CASE(4)
        DIAGNOSTICS_LEVEL4=.TRUE.
      CASE(5)
        DIAGNOSTICS_LEVEL5=.TRUE.
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid diagnostic level",ERR,ERROR,*999)
      END SELECT
    ENDDO !i
    DIAGNOSTICS=.TRUE.
    DIAG_OR_TIMING=.TRUE.

    CALL EXITS("DIAGNOSTICS_SET_ON")
    RETURN
999 IF(DIAG_FILE_OPEN) THEN
      CLOSE(UNIT=DIAGNOSTICS_FILE_UNIT)
      DIAG_FILE_OPEN=.FALSE.
    ENDIF
    ROUTINE=>DIAG_ROUTINE_LIST%HEAD
    DO WHILE(ASSOCIATED(ROUTINE))
      NEXT_ROUTINE=>ROUTINE%NEXT_ROUTINE
      DEALLOCATE(ROUTINE)
      ROUTINE=>NEXT_ROUTINE
    ENDDO
    NULLIFY(DIAG_ROUTINE_LIST%HEAD)
    DIAG_ALL_SUBROUTINES=.FALSE.
    DIAG_FROM_SUBROUTINE=.FALSE.
    DIAGNOSTICS_LEVEL1=.FALSE.
    DIAGNOSTICS_LEVEL2=.FALSE.
    DIAGNOSTICS_LEVEL3=.FALSE.
    DIAGNOSTICS_LEVEL4=.FALSE.
    DIAGNOSTICS_LEVEL5=.FALSE.
    DIAGNOSTICS=.FALSE.
    DIAG_OR_TIMING=TIMING
    CALL ERRORS("DIAGNOSTICS_SET_ON",ERR,ERROR)
    CALL EXITS("DIAGNOSTICS_SET_ON")
    RETURN 1
  END SUBROUTINE DIAGNOSTICS_SET_ON

  !
  !================================================================================================================================
  !

  !>Sets writes file echo output off.
  SUBROUTINE OUTPUT_SET_OFF(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    CALL ENTERS("OUTPUT_SET_OFF",ERR,ERROR,*999)

    IF(ECHO_OUTPUT) THEN
      ECHO_OUTPUT=.FALSE.
      CLOSE(UNIT=ECHO_FILE_UNIT)
    ELSE
      CALL FLAG_ERROR("Write output is not on",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("OUTPUT_SET_OFF")
    RETURN
999 CALL ERRORS("OUTPUT_SET_OFF",ERR,ERROR)
    CALL EXITS("OUTPUT_SET_OFF")
    RETURN 1
  END SUBROUTINE OUTPUT_SET_OFF

  !
  !================================================================================================================================
  !

  !>Sets writes file echo output on.
  SUBROUTINE OUTPUT_SET_ON(ECHO_FILENAME,ERR,ERROR,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: ECHO_FILENAME !<The filename of the file to echo output to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: FILENAME

    CALL ENTERS("OUTPUT_SET_ON",ERR,ERROR,*999)

    IF(ECHO_OUTPUT) THEN
      CALL FLAG_ERROR("Write output is already on",ERR,ERROR,*999)
    ELSE
      FILENAME=ECHO_FILENAME(1:LEN_TRIM(ECHO_FILENAME))//".out"
      OPEN(UNIT=ECHO_FILE_UNIT,FILE=FILENAME(1:LEN_TRIM(FILENAME)),STATUS="UNKNOWN",IOSTAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not open write output file",ERR,ERROR,*999)
      ECHO_OUTPUT=.TRUE.
    ENDIF

    CALL EXITS("OUTPUT_SET_ON")
    RETURN
999 CALL ERRORS("OUTPUT_SET_ON",ERR,ERROR)
    CALL EXITS("OUTPUT_SET_ON")
    RETURN 1
  END SUBROUTINE OUTPUT_SET_ON

  !
  !================================================================================================================================
  !

  !>Sets timing off.
  SUBROUTINE TIMING_SET_OFF(ERR,ERROR,*)

   !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: NEXT_ROUTINE,ROUTINE

    CALL ENTERS("TIMING_SET_OFF",ERR,ERROR,*999)

    IF(TIMING) THEN
      IF(TIMING_FILE_OPEN) THEN
        TIMING_FILE_OPEN=.FALSE.
        CLOSE(UNIT=TIMING_FILE_UNIT)
      ENDIF
      IF(TIMING_ALL_SUBROUTINES) THEN
        TIMING_ALL_SUBROUTINES=.FALSE.
      ELSE
        ROUTINE=>TIMING_ROUTINE_LIST%HEAD
        DO WHILE(ASSOCIATED(ROUTINE))
          NEXT_ROUTINE=>ROUTINE%NEXT_ROUTINE
          DEALLOCATE(ROUTINE)
          ROUTINE=>NEXT_ROUTINE
        ENDDO
        NULLIFY(TIMING_ROUTINE_LIST%HEAD)
        TIMING_FROM_SUBROUTINE=.FALSE.
      ENDIF
      TIMING_SUMMARY=.FALSE.
      TIMING=.FALSE.
      DIAG_OR_TIMING=DIAGNOSTICS
    ELSE
      CALL FLAG_ERROR("Timing is not on",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TIMING_SET_OFF")
    RETURN
999 CALL ERRORS("TIMING_SET_OFF",ERR,ERROR)
    CALL EXITS("TIMING_SET_OFF")
    RETURN 1
  END SUBROUTINE TIMING_SET_OFF

  !
  !================================================================================================================================
  !

  !>Sets timing on.
  SUBROUTINE TIMING_SET_ON(TIMING_TYPE,TIMING_SUMMARY_FLAG,TIMING_FILENAME,ROUTINE_LIST,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TIMING_TYPE !<The type of timing to set on \see BASE_ROUTINES_TimingTypes
    LOGICAL, INTENT(IN) :: TIMING_SUMMARY_FLAG !<.TRUE. if the timing information will be output with subsequent TIMING_SUMMARY_OUTPUT calls, .FALSE. if the timing information will be output every time the routine exits 
    CHARACTER(LEN=*), INTENT(IN) :: TIMING_FILENAME !<If present the name of the file to output timing information to. If omitted the timing output is sent to the screen
    CHARACTER(LEN=*), INTENT(IN) :: ROUTINE_LIST(:) !<The list of routines to set timing on in.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("A")-ICHAR("a"))
    CHARACTER(LEN=1) :: CHARAC1,CHARAC2
    CHARACTER(LEN=MAXSTRLEN) :: FILENAME
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: NEXT_ROUTINE,PREVIOUS_ROUTINE,ROUTINE
    TYPE(VARYING_STRING) :: LOCAL_STRING

    CALL ENTERS("TIMING_SET_ON",ERR,ERROR,*999)

    NULLIFY(ROUTINE)
    IF(LEN_TRIM(TIMING_FILENAME)>=1) THEN
      IF(TIMING_FILE_OPEN) CLOSE(UNIT=TIMING_FILE_UNIT)
      FILENAME=TIMING_FILENAME(1:LEN_TRIM(TIMING_FILENAME))//".timing"
      OPEN(UNIT=TIMING_FILE_UNIT,FILE=FILENAME(1:LEN_TRIM(FILENAME)),STATUS="UNKNOWN",IOSTAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not open timing file",ERR,ERROR,*999)
      TIMING_FILE_OPEN=.TRUE.
    ENDIF
    SELECT CASE(TIMING_TYPE)
    CASE(ALL_TIMING_TYPE)
      TIMING_ALL_SUBROUTINES=.TRUE.
    CASE(IN_TIMING_TYPE,FROM_TIMING_TYPE)
      TIMING_ALL_SUBROUTINES=.FALSE.
      IF(ASSOCIATED(TIMING_ROUTINE_LIST%HEAD)) THEN
        ROUTINE=>TIMING_ROUTINE_LIST%HEAD
        DO WHILE(ASSOCIATED(ROUTINE))
          NEXT_ROUTINE=>ROUTINE%NEXT_ROUTINE
          DEALLOCATE(ROUTINE)
          ROUTINE=>NEXT_ROUTINE
        ENDDO
        NULLIFY(TIMING_ROUTINE_LIST%HEAD)
      ENDIF
      ALLOCATE(ROUTINE,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate routine list item",ERR,ERROR,*999)
      ROUTINE%NAME=""
      DO i=1,LEN(ROUTINE_LIST(1))
        CHARAC1=ROUTINE_LIST(1)(i:i)
        IF(LGE(CHARAC1,"a").AND.LLE(CHARAC1,"z")) THEN
          CHARAC2=CHAR(ICHAR(CHARAC1)+OFFSET)
        ELSE
          CHARAC2=CHARAC1
        ENDIF
        !CPB 26/9/05 Aix compiler doesn't like the vstring so split the statement and put a char() around it
        !ROUTINE%NAME=ROUTINE%NAME//CHARAC2
        LOCAL_STRING=ROUTINE%NAME//CHARAC2
        ROUTINE%NAME=CHAR(LOCAL_STRING)
      ENDDO !i
      PREVIOUS_ROUTINE=>ROUTINE
      NULLIFY(ROUTINE%NEXT_ROUTINE)
      TIMING_ROUTINE_LIST%HEAD=>ROUTINE
      ROUTINE%NUMBER_OF_INVOCATIONS=0
      ROUTINE%TOTAL_INCLUSIVE_CPU_TIME=0.0_SP
      ROUTINE%TOTAL_INCLUSIVE_SYSTEM_TIME=0.0_SP
      ROUTINE%TOTAL_EXCLUSIVE_CPU_TIME=0.0_SP
      ROUTINE%TOTAL_EXCLUSIVE_SYSTEM_TIME=0.0_SP
      DO i=2,SIZE(ROUTINE_LIST,1)
        ALLOCATE(ROUTINE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate routine list item",ERR,ERROR,*999)
        ROUTINE%NAME=""
        DO j=1,LEN(ROUTINE_LIST(i))
          CHARAC1=ROUTINE_LIST(i)(j:j)
          IF(LGE(CHARAC1,"a").AND.LLE(CHARAC1,"z")) THEN
            CHARAC2=CHAR(ICHAR(CHARAC1)+OFFSET)
          ELSE
            CHARAC2=CHARAC1
          ENDIF
          !CPB 26/9/05 Aix compiler doesn't like the vstring so split the statement and put a char() around it
          !ROUTINE%NAME=ROUTINE%NAME//CHARAC2
          LOCAL_STRING=ROUTINE%NAME//CHARAC2
          ROUTINE%NAME=CHAR(LOCAL_STRING)
        ENDDO !i
        NULLIFY(ROUTINE%NEXT_ROUTINE)
        PREVIOUS_ROUTINE%NEXT_ROUTINE=>ROUTINE
        PREVIOUS_ROUTINE=>ROUTINE
        ROUTINE%NUMBER_OF_INVOCATIONS=0
        ROUTINE%TOTAL_INCLUSIVE_CPU_TIME=0.0_SP
        ROUTINE%TOTAL_INCLUSIVE_CPU_TIME=0.0_SP
        ROUTINE%TOTAL_EXCLUSIVE_CPU_TIME=0.0_SP
        ROUTINE%TOTAL_EXCLUSIVE_CPU_TIME=0.0_SP
      ENDDO !i
    CASE DEFAULT
      CALL FLAG_ERROR("Invalid timing type",ERR,ERROR,*999)
    END SELECT
    TIMING_SUMMARY=TIMING_SUMMARY_FLAG
    TIMING=.TRUE.
    DIAG_OR_TIMING=.TRUE.

    CALL EXITS("TIMING_SET_ON")
    RETURN
999 IF(TIMING_FILE_OPEN) THEN
      CLOSE(UNIT=TIMING_FILE_UNIT)
      TIMING_FILE_OPEN=.FALSE.
    ENDIF
    ROUTINE=>TIMING_ROUTINE_LIST%HEAD
    DO WHILE(ASSOCIATED(ROUTINE))
      NEXT_ROUTINE=>ROUTINE%NEXT_ROUTINE
      DEALLOCATE(ROUTINE)
      ROUTINE=>NEXT_ROUTINE
    ENDDO
    NULLIFY(TIMING_ROUTINE_LIST%HEAD)
    TIMING_ALL_SUBROUTINES=.FALSE.
    TIMING_FROM_SUBROUTINE=.FALSE.
    TIMING=.FALSE.
    DIAG_OR_TIMING=DIAGNOSTICS
    CALL ERRORS("TIMING_SET_ON",ERR,ERROR)
    CALL EXITS("TIMING_SET_ON")
    RETURN 1
  END SUBROUTINE TIMING_SET_ON

  !
  !================================================================================================================================
  !

  !>Outputs the timing summary.
  SUBROUTINE TIMING_SUMMARY_OUTPUT(ERR,ERROR,*)    

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(ROUTINE_LIST_ITEM_TYPE), POINTER :: ROUTINE_PTR

    NULLIFY(ROUTINE_PTR)
    
    CALL ENTERS("TIMING_SUMMARY_OUTPUT",ERR,ERROR,*999)

    IF(TIMING) THEN
      WRITE(OP_STRING,'("*** Timing Summary: ")') 
      CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
      ROUTINE_PTR=>TIMING_ROUTINE_LIST%HEAD
      DO WHILE(ASSOCIATED(ROUTINE_PTR))
        WRITE(OP_STRING,'("*** Routine : ",A)') CHAR(TRIM(ROUTINE_PTR%NAME))
        CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
        WRITE(OP_STRING,'("***    Number of invocations: ",I10)') ROUTINE_PTR%NUMBER_OF_INVOCATIONS
        CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
        WRITE(OP_STRING,'("***    Routine times: Total Exclusive  Total Inclusive  Average Exclusive  Average Inclusive")')
        CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
        IF(ROUTINE_PTR%NUMBER_OF_INVOCATIONS==0) THEN
          WRITE(OP_STRING,'("***    CPU       (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & ROUTINE_PTR%TOTAL_EXCLUSIVE_CPU_TIME,ROUTINE_PTR%TOTAL_INCLUSIVE_CPU_TIME, &
            & REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP),REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP)
          CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
          WRITE(OP_STRING,'("***    System    (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & ROUTINE_PTR%TOTAL_EXCLUSIVE_SYSTEM_TIME,ROUTINE_PTR%TOTAL_INCLUSIVE_SYSTEM_TIME, &
            & REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP),REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP)
          CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
        ELSE
          WRITE(OP_STRING,'("***    CPU       (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & ROUTINE_PTR%TOTAL_EXCLUSIVE_CPU_TIME,ROUTINE_PTR%TOTAL_INCLUSIVE_CPU_TIME, &
            & ROUTINE_PTR%TOTAL_EXCLUSIVE_CPU_TIME/REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP), &
            & ROUTINE_PTR%TOTAL_INCLUSIVE_CPU_TIME/REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP)
          CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
          WRITE(OP_STRING,'("***    System    (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & ROUTINE_PTR%TOTAL_EXCLUSIVE_SYSTEM_TIME,ROUTINE_PTR%TOTAL_INCLUSIVE_SYSTEM_TIME, &
            & ROUTINE_PTR%TOTAL_EXCLUSIVE_SYSTEM_TIME/REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP), &
            & ROUTINE_PTR%TOTAL_INCLUSIVE_SYSTEM_TIME/REAL(ROUTINE_PTR%NUMBER_OF_INVOCATIONS,SP)
          CALL WRITE_STR(TIMING_OUTPUT_TYPE,ERR,ERROR,*999)
        ENDIF
        ROUTINE_PTR=>ROUTINE_PTR%NEXT_ROUTINE
      ENDDO
    ELSE
      CALL FLAG_ERROR("Timing is not on",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TIMING_SUMMARY_OUTPUT")
    RETURN
999 CALL ERRORS("TIMING_SUMMARY_OUTPUT",ERR,ERROR)
    CALL EXITS("TIMING_SUMMARY_OUTPUT")
    RETURN 1
  END SUBROUTINE TIMING_SUMMARY_OUTPUT

 !
  !================================================================================================================================
  !

  !>Writes the error string.
  SUBROUTINE WRITE_ERROR(ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: INDENT,LOCAL_ERR,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: INDENT_STRING=">>"
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_ERROR2
    
    POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
    WRITE(OP_STRING,'(">>ERROR: ",I5,": ",A)') ERR,CHAR(EXTRACT(ERROR,1,POSITION-1))
    CALL WRITE_STR(ERROR_OUTPUT_TYPE,LOCAL_ERR,LOCAL_ERROR2,*999)
    !CPB 20/02/07 aix compiler does not like varying strings so split the remove statement up into two statements
    LOCAL_ERROR=REMOVE(ERROR,1,POSITION)
    ERROR=LOCAL_ERROR
    POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
    INDENT=4
    DO WHILE(POSITION/=0)
      WRITE(OP_STRING,'(A)') INDENT_STRING(1:INDENT)//CHAR(EXTRACT(ERROR,1,POSITION-1))
      CALL WRITE_STR(ERROR_OUTPUT_TYPE,LOCAL_ERR,LOCAL_ERROR2,*999)
      !CPB 20/02/07 aix compiler does not like varying strings so split the remove statement up into two statements
      LOCAL_ERROR=REMOVE(ERROR,1,POSITION)
      ERROR=LOCAL_ERROR
      POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
      INDENT=INDENT+2
    ENDDO
    WRITE(OP_STRING,'(A)') INDENT_STRING(1:INDENT)//CHAR(ERROR)
    CALL WRITE_STR(ERROR_OUTPUT_TYPE,LOCAL_ERR,LOCAL_ERROR2,*999)
    ERR=0
    ERROR=""

    RETURN
    !Don't return an error code here otherwise we will get into a circular loop
999 RETURN 
  END SUBROUTINE WRITE_ERROR

  !
  !================================================================================================================================
  !

  !>Writes the output string to a specified output stream.
  SUBROUTINE WRITE_STR(ID,ERR,ERROR,*)


!!!!NOTE: No enters or exits is used here to avoid an infinite loop
!!!!NOTE: This routine is, in general, OS dependent but needs to be defined here so that an module loop is avoided when MACHINE
!!!!      module routines need to use this module.

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: END_LINE(MAX_OUTPUT_LINES),i,j,LENGTH,NUMBER_BLANKS,NUMBER_RECORDS

    !Calculate number of records in OP_STRING
    SELECT CASE(MACHINE_OS)
    CASE(VMS_OS)
      i=1 
      DO WHILE(OP_STRING(i)(1:1)/=CHAR(0).AND.i<MAX_OUTPUT_LINES)
        i=i+1
      ENDDO
      NUMBER_RECORDS=i-1
    CASE(IRIX_OS,LINUX_OS,AIX_OS)
      i=1
      NUMBER_BLANKS=0
      DO WHILE(i<MAX_OUTPUT_LINES.AND.NUMBER_BLANKS<2)
        i=i+1
        LENGTH=LEN_TRIM(OP_STRING(i))
        IF(LENGTH==0) THEN
          NUMBER_BLANKS=NUMBER_BLANKS+1
        ELSE
          NUMBER_BLANKS=0
        ENDIF
      ENDDO
      NUMBER_RECORDS=i-NUMBER_BLANKS
    CASE(WINDOWS_OS)
      i=1 
      DO WHILE(OP_STRING(i)(1:1)/=CHAR(0).AND.i<MAX_OUTPUT_LINES)
        i=i+1
      ENDDO
      NUMBER_RECORDS=i-1
    CASE DEFAULT
      CALL FLAG_ERROR("Operating system not implemented",ERR,ERROR,*999)
    END SELECT

    DO i=1,NUMBER_RECORDS
      END_LINE(i)=LEN_TRIM(OP_STRING(i))
    ENDDO !i

    IF(DIAG_FILE_OPEN.AND.ID==DIAGNOSTIC_OUTPUT_TYPE) THEN
      DO i=1,NUMBER_RECORDS
        IF(END_LINE(i)<=132) THEN
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') OP_STRING(i)(1:END_LINE(i))
        ELSE IF(END_LINE(i)>132.AND.END_LINE(i)<=MAXSTRLEN) THEN
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') OP_STRING(i)(1:132)
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') OP_STRING(i)(133:END_LINE(i))
        ELSE
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') OP_STRING(i)(1:132)
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') OP_STRING(i)(133:MAXSTRLEN)
        ENDIF
      ENDDO !i
    ELSE IF(TIMING_FILE_OPEN.AND.ID==TIMING_OUTPUT_TYPE) THEN
      DO i=1,NUMBER_RECORDS
        IF(END_LINE(i)<=132) THEN
          WRITE(TIMING_FILE_UNIT,'(A)') OP_STRING(i)(1:END_LINE(i))
        ELSE IF(END_LINE(i)>132.AND.END_LINE(i)<=MAXSTRLEN) THEN
          WRITE(TIMING_FILE_UNIT,'(A)') OP_STRING(i)(1:132)
          WRITE(TIMING_FILE_UNIT,'(A)') OP_STRING(i)(133:END_LINE(i))
        ELSE
          WRITE(TIMING_FILE_UNIT,'(A)') OP_STRING(i)(1:132)
          WRITE(TIMING_FILE_UNIT,'(A)') OP_STRING(i)(133:MAXSTRLEN)
        ENDIF
      ENDDO !i
    ELSE
      IF(ID<=9) THEN !not file output
        DO i=1,NUMBER_RECORDS
          IF(END_LINE(i)<=132) THEN
            WRITE(*,'(A)') OP_STRING(i)(1:END_LINE(i))
          ELSE IF(END_LINE(i)>132.AND.END_LINE(i)<=MAXSTRLEN) THEN
            WRITE(*,'(A)') OP_STRING(i)(1:132)
            WRITE(*,'(A)') OP_STRING(i)(133:END_LINE(i))
          ELSE
            WRITE(*,'(A)') OP_STRING(i)(1:132)
            WRITE(*,'(A)') OP_STRING(i)(133:MAXSTRLEN)
          ENDIF
        ENDDO !i
      ELSE !file output
        DO i=1,NUMBER_RECORDS
          WRITE(ID,'(A)') OP_STRING(i)(1:END_LINE(i))
        ENDDO !i
      ENDIF

      !Echo strings to output file if required

      IF(ECHO_OUTPUT) THEN
        DO i=1,NUMBER_RECORDS
          IF(END_LINE(i)<=132) THEN
            WRITE(ECHO_FILE_UNIT,'(A)') OP_STRING(i)(1:END_LINE(i))
          ELSE IF(END_LINE(i)>132.AND.END_LINE(i)<=MAXSTRLEN) THEN
            WRITE(ECHO_FILE_UNIT,'(A)') OP_STRING(i)(1:132)
            WRITE(ECHO_FILE_UNIT,'(A)') OP_STRING(i)(133:END_LINE(i))
          ELSE
            WRITE(ECHO_FILE_UNIT,'(A)') OP_STRING(i)(1:132)
            WRITE(ECHO_FILE_UNIT,'(A)') OP_STRING(i)(133:MAXSTRLEN)
          ENDIF
        ENDDO !i
      ENDIF
    ENDIF

    !Reset OP_STRING
    SELECT CASE(MACHINE_OS)
    CASE(VMS_OS)
      DO i=1,NUMBER_RECORDS
        OP_STRING(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE(IRIX_OS,LINUX_OS,AIX_OS)
      DO i=1,NUMBER_RECORDS
        DO j=1,MAXSTRLEN
          OP_STRING(i)(j:j)=' '
        ENDDO !j
      ENDDO !i
    CASE(WINDOWS_OS)
      DO i=1,NUMBER_RECORDS
        OP_STRING(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE DEFAULT
      CALL FLAG_ERROR("Operating system not implemented",ERR,ERROR,*999)
    END SELECT

    RETURN
999 CALL ERRORS("WRITE_STR",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_STR

  !
  !================================================================================================================================
  !

END MODULE BASE_ROUTINES
