!> \file
!> \author Chris Bradley
!> \brief This module handles all history file routines.
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
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
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

!> This module handles all history file routines.
MODULE HISTORY_ROUTINES
  
  USE BASE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

#include "macros.h"  
  
  IMPLICIT NONE
  
  PRIVATE
  
  !Module parameters

  !> \addtogroup HISTORY_ROUTINES_FileFormatTypes HISTORY_ROUTINES::FileFormatTypes
  !> \brief The types of a file format for a history file.
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: HISTORY_ASCII_FILE_FORMAT=1 !<ASCII history file format \see HISTORY_ROUTINES_FileFormatTypes,HISTORY_ROUTINES
  INTEGER(INTG), PARAMETER :: HISTORY_BINARY_FILE_FORMAT=2 !Binary history file format \see HISTORY_ROUTINES_FileFormatTypes,HISTORY_ROUTINES
  !>@}

  
  !Module types

  !Module variables
  
  !Interfaces

  !>Sets/changes the filename for a history file 
  INTERFACE HISTORY_FILENAME_SET
    MODULE PROCEDURE HISTORY_FILENAME_SET_C
    MODULE PROCEDURE HISTORY_FILENAME_SET_VS
  END INTERFACE !HISTORY_FILENAME_SET

  PUBLIC HISTORY_CREATE_FINISH,HISTORY_CREATE_START,HISTORY_DESTROY,HISTORY_FILE_FORMAT_SET,HISTORY_FILENAME_SET
  
CONTAINS
  
  !
  !================================================================================================================================
  !
  
  !>Closes a history file. 
  SUBROUTINE HISTORY_CLOSE(HISTORY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to the history file to close
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_CLOSE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      HISTORY%UNIT_NUMBER=0
    ELSE
      CALL FlagError("History is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_CLOSE")
    RETURN
999 ERRORSEXITS("HISTORY_CLOSE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_CLOSE

  !
  !================================================================================================================================
  !
  
  !>Finishes the process of creating a history file.
  SUBROUTINE HISTORY_CREATE_FINISH(HISTORY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<On exit, a pointer to the created history file.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_CREATE_FINISH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      HISTORY%HISTORY_FINISHED=.TRUE.      
    ELSE
      CALL FlagError("History is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("HISTORY_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_CREATE_FINISH
  
  !
  !================================================================================================================================
  !
  
  !>Starts the process of creating a history file.
  SUBROUTINE HISTORY_CREATE_START(CONTROL_LOOP,HISTORY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to setup
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<On exit, a pointer to the created history file. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(HISTORY)) THEN
        CALL FlagError("History is already associated.",ERR,ERROR,*999)
      ELSE
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_CREATE_START")
    RETURN
999 ERRORSEXITS("HISTORY_CREATE_START",ERR,ERROR)
    RETURN 1    
  END SUBROUTINE HISTORY_CREATE_START
  
  !
  !================================================================================================================================
  !
  
  !>Destroys a history file.
  SUBROUTINE HISTORY_DESTROY(HISTORY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to the history file to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    ENTERS("HISTORY_DESTROY",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      CALL HISTORY_FINALISE(HISTORY,ERR,ERROR,*999)
    ELSE
      CALL FlagError("History is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_DESTROY")
    RETURN
999 ERRORSEXITS("HISTORY_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_DESTROY
  
  !
  !================================================================================================================================
  !
  
  !>Finalises a history file and deallocates all memory. 
  SUBROUTINE HISTORY_FINALISE(HISTORY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to the history file to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_FINALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      IF(HISTORY%UNIT_NUMBER/=0) CALL HISTORY_CLOSE(HISTORY,ERR,ERROR,*999)
      DEALLOCATE(HISTORY)
    ENDIF
    
    EXITS("HISTORY_FINALISE")
    RETURN
999 ERRORSEXITS("HISTORY_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_FINALISE

  !
  !================================================================================================================================
  !
  
  !>Initialises a history file. 
  SUBROUTINE HISTORY_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control_loop to initialise the history file for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%HISTORY)) THEN
        CALL FlagError("Control loop history is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(CONTROL_LOOP%HISTORY,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate control loop history.",ERR,ERROR,*999)
        CONTROL_LOOP%HISTORY%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%HISTORY%HISTORY_FINISHED=.FALSE.
        CONTROL_LOOP%HISTORY%FILE_FORMAT=HISTORY_BINARY_FILE_FORMAT
        CONTROL_LOOP%HISTORY%FILENAME="History"
        CONTROL_LOOP%HISTORY%UNIT_NUMBER=0
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_INITIALISE")
    RETURN
999 ERRORSEXITS("HISTORY_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the file type for a history file 
  SUBROUTINE HISTORY_FILE_FORMAT_SET(HISTORY,FILE_FORMAT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to history to set the file format for
    INTEGER(INTG), INTENT(IN) :: FILE_FORMAT !<The file format to set \see HISTORY_ROUTINES_FileFormatTypes,HISTORY_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HISTORY_FILE_FORMAT_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      IF(HISTORY%HISTORY_FINISHED) THEN
        CALL FlagError("History has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(FILE_FORMAT)
        CASE(HISTORY_ASCII_FILE_FORMAT)
          HISTORY%FILE_FORMAT=HISTORY_ASCII_FILE_FORMAT
        CASE(HISTORY_BINARY_FILE_FORMAT)
          HISTORY%FILE_FORMAT=HISTORY_BINARY_FILE_FORMAT
        CASE DEFAULT
          LOCAL_ERROR="The supplied file format of "//TRIM(NUMBER_TO_VSTRING(FILE_FORMAT,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("History is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_FILE_FORMAT_SET")
    RETURN
999 ERRORSEXITS("HISTORY_FILE_FORMAT_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_FILE_FORMAT_SET

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the character string filename for a history file 
  SUBROUTINE HISTORY_FILENAME_SET_C(HISTORY,FILENAME,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to history to set the filename for
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME !<The filename to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_FILENAME_SET_C",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      IF(HISTORY%HISTORY_FINISHED) THEN
        CALL FlagError("History has been finished.",ERR,ERROR,*999)
      ELSE
!!TODO: Check filename???
        HISTORY%FILENAME=FILENAME
      ENDIF
    ELSE
      CALL FlagError("History is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_FILENAME_SET_C")
    RETURN
999 ERRORSEXITS("HISTORY_FILENAME_SET_C",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_FILENAME_SET_C

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the varying string filename for a history file 
  SUBROUTINE HISTORY_FILENAME_SET_VS(HISTORY,FILENAME,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to history to set the filename for
    TYPE(VARYING_STRING), INTENT(IN) :: FILENAME !<The filename to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_FILENAME_SET_VS",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      IF(HISTORY%HISTORY_FINISHED) THEN
        CALL FlagError("History has been finished.",ERR,ERROR,*999)
      ELSE
!!TODO: Check filename???
        HISTORY%FILENAME=FILENAME
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_FILENAME_SET_VS")
    RETURN
999 ERRORSEXITS("HISTORY_FILENAME_SET_VS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_FILENAME_SET_VS

  !
  !================================================================================================================================
  !
  
  !>Opens a history file. 
  SUBROUTINE HISTORY_OPEN(HISTORY,ERR,ERROR,*)
    
    !Argument variables
    TYPE(HISTORY_TYPE), POINTER :: HISTORY !<A pointer to the history file to open
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("HISTORY_OPEN",ERR,ERROR,*999)
    
    IF(ASSOCIATED(HISTORY)) THEN
      IF(HISTORY%HISTORY_FINISHED) THEN
      ELSE
        CALL FlagError("History has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("History is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("HISTORY_OPEN")
    RETURN
999 ERRORSEXITS("HISTORY_OPEN",ERR,ERROR)
    RETURN 1
  END SUBROUTINE HISTORY_OPEN

  !
  !================================================================================================================================
  !

END MODULE HISTORY_ROUTINES

