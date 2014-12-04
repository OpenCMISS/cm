!> \file
!> \author Chris Bradley
!> \brief The top level cmiss module.
!>
!> \mainpage OpenCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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
!>
!> The top level cmiss module.
MODULE CMISS

  USE ISO_C_BINDING
  
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MACHINE_CONSTANTS
  USE MPI
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(INTG), PARAMETER :: CMISS_MAJOR_VERSION = 0
  INTEGER(INTG), PARAMETER :: CMISS_MINOR_VERSION = 3
  INTEGER(INTG), PARAMETER :: CMISS_REVISION_VERSION = 0

  CHARACTER(LEN=MAXSTRLEN), PARAMETER :: CMISS_BUILD_VERSION = "$Rev"

  !> \addtogroup CMISS_ErrorHandlingModes CMISS::ErrorHandlingModes
  !> \brief Error handling mode parameters
  !> \see CMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMISS_RETURN_ERROR_CODE = 0 !<Just return the error code \see CMISS_ErrorHandlingModes,CMISS
  INTEGER(INTG), PARAMETER :: CMISS_OUTPUT_ERROR = 1 !<Output the error traceback and return the error code \see CMISS_ErrorHandlingModes,CMISS
  INTEGER(INTG), PARAMETER :: CMISS_TRAP_ERROR = 2 !<Trap the error by outputing the error traceback and stopping the program \see CMISS_ErrorHandlingModes,CMISS
  !>@}
  
  !Module types

  !Module variables

  INTEGER(INTG), SAVE :: CMISS_ERROR_HANDLING_MODE !<The current error handling mode for OpenCMISS \see CMISS_ErrorHandlingModes
 
  !Interfaces

  INTERFACE

    SUBROUTINE CMISSInitFatalHandler() BIND(C,NAME="CMISSInitFatalHandler")
    END SUBROUTINE CMISSINITFATALHANDLER

    SUBROUTINE CMISSResetFatalHandler() BIND(C,NAME="CMISSResetFatalHandler")
    END SUBROUTINE CMISSRESETFATALHANDLER
    
    SUBROUTINE CMISSSetFatalHandler() BIND(C,NAME="CMISSSetFatalHandler")
    END SUBROUTINE CMISSSETFATALHANDLER

  END INTERFACE

  !Allow using new code style
  INTERFACE CMISSHandleError
    MODULE PROCEDURE CMISS_HANDLE_ERROR
  END INTERFACE !CMISSHandleError

  PUBLIC CMISS_MAJOR_VERSION,CMISS_MINOR_VERSION,CMISS_REVISION_VERSION,CMISS_BUILD_VERSION

  PUBLIC CMISS_RETURN_ERROR_CODE,CMISS_OUTPUT_ERROR,CMISS_TRAP_ERROR

  PUBLIC CMISS_ERROR_HANDLING_MODE_GET,CMISS_ERROR_HANDLING_MODE_SET
  
  PUBLIC CMISS_HANDLE_ERROR,CMISSHandleError
  
  PUBLIC CMISS_WRITE_ERROR
  
  PUBLIC CMISS_FINALISE,CMISS_INITIALISE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns the error handling mode for CMISS \see OPENCMISS::CMISSErrorHandlingModeGet
  SUBROUTINE CMISS_ERROR_HANDLING_MODE_GET(ERROR_HANDLING_MODE,ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERROR_HANDLING_MODE !<On return, the error handling mode. \see CMISS_ErrorHandlingModes,CMISS
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error code
    !Local Variables

!    CALL ENTERS("CMISS_ERROR_HANDLING_MODE_GET",ERR,ERROR,*999)

    ERROR_HANDLING_MODE=CMISS_ERROR_HANDLING_MODE
    
!    CALL EXITS("CMISS_ERROR_HANDLING_MODE_GET")
    RETURN
999 CALL ERRORS("CMISS_ERROR_HANDLING_MODE_GET",ERR,ERROR)
!    CALL EXITS("CMISS_ERROR_HANDLING_MODE_GET")
    RETURN 1
  END SUBROUTINE CMISS_ERROR_HANDLING_MODE_GET

  !
  !================================================================================================================================
  !

  !>Sets the error handling mode for cmiss \see OPENCMISS::CMISSErrorHandlingModeSet
  SUBROUTINE CMISS_ERROR_HANDLING_MODE_SET(ERROR_HANDLING_MODE,ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ERROR_HANDLING_MODE !<The error handling mode to set. \see CMISS_ErrorHandlingModes,CMISS
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error code
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
!    CALL ENTERS("CMISS_ERROR_HANDLING_MODE_SET",ERR,ERROR,*999)

    SELECT CASE(ERROR_HANDLING_MODE)
    CASE(CMISS_RETURN_ERROR_CODE)
      CMISS_ERROR_HANDLING_MODE=CMISS_RETURN_ERROR_CODE
    CASE(CMISS_OUTPUT_ERROR)
      CMISS_ERROR_HANDLING_MODE=CMISS_OUTPUT_ERROR
    CASE(CMISS_TRAP_ERROR)
      CMISS_ERROR_HANDLING_MODE=CMISS_TRAP_ERROR
    CASE DEFAULT
      LOCAL_ERROR="The supplied error handling mode of "//TRIM(NUMBER_TO_VSTRING(ERROR_HANDLING_MODE,"*",ERR,ERROR))// &
        & " is invalid."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

 !    CALL EXITS("CMISS_ERROR_HANDLING_MODE_SET")
    RETURN
999 CALL ERRORS("CMISS_ERROR_HANDLING_MODE_SET",ERR,ERROR)
!    CALL EXITS("CMISS_ERROR_HANDLING_MODE_SET")
    RETURN 1
  END SUBROUTINE CMISS_ERROR_HANDLING_MODE_SET

  !
  !================================================================================================================================
  !

  !>Finalises CMISS. \see OPENCMISS::CMISSFinalise
  SUBROUTINE CMISS_FINALISE(ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error code
    !Local Variables

    !Finalise the problems
    CALL PROBLEMS_FINALISE(ERR,ERROR,*999)
    !Finalise the regions
    CALL REGIONS_FINALISE(ERR,ERROR,*999)
    !Finalise the coordinate systems
    CALL COORDINATE_SYSTEMS_FINALISE(ERR,ERROR,*999)
    !Finalise bases
    CALL BASES_FINALISE(ERR,ERROR,*999)
    !Reset the signal handler
    CALL CMISSResetFatalHandler()
    !Finalise computational enviroment
    CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(ERR,ERROR,*999)
    !Finalise the base routines
    CALL BASE_ROUTINES_FINALISE(ERR,ERROR,*999)
     
    RETURN
999 RETURN 1
  END SUBROUTINE CMISS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises CMISS. \see OPENCMISS::CMISSInitialise
  SUBROUTINE CMISS_INITIALISE(WORLD_REGION,ERR,ERROR,*)
  
    !Argument variables
    !TYPE(REGION_TYPE), POINTER, INTENT(OUT) :: WORLD_REGION !<On exit, a pointer to the world region. Must not be associated on entry.
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION !<On exit, a pointer to the world region. Must not be associated on entry.
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: VERSION_STRING

    !Initialise error mode
    CMISS_ERROR_HANDLING_MODE = CMISS_OUTPUT_ERROR !Default for now, maybe make CMISS_RETURN_ERROR_CODE the default
    !Initialise the base routines
    CALL BASE_ROUTINES_INITIALISE(ERR,ERROR,*999)
    !Intialise the computational environment
    CALL COMPUTATIONAL_ENVIRONMENT_INITIALISE(ERR,ERROR,*999)
    !Setup signal handling
    CALL CMISSInitFatalHandler()
    CALL CMISSSetFatalHandler()
    IF(ASSOCIATED(WORLD_REGION)) THEN
      CALL FLAG_ERROR("World region is already associated.",ERR,ERROR,*999)
    ELSE
      !Intialise the bases
      CALL BASES_INITIALISE(ERR,ERROR,*999)
      !Initialise the coordinate systems
      CALL COORDINATE_SYSTEMS_INITIALISE(ERR,ERROR,*999)
      !Initialise the regions 
      CALL REGIONS_INITIALISE(WORLD_REGION,ERR,ERROR,*999)
      !Initialise the problems
      CALL PROBLEMS_INITIALISE(ERR,ERROR,*999)
      
      !Write out the CMISS version
      IF(COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
        VERSION_STRING="OpenCMISS(cm) version "//TRIM(NUMBER_TO_VSTRING(CMISS_MAJOR_VERSION,"*",ERR,ERROR))
        VERSION_STRING=VERSION_STRING//"."
        VERSION_STRING=VERSION_STRING//TRIM(NUMBER_TO_VSTRING(CMISS_MINOR_VERSION,"*",ERR,ERROR))
        VERSION_STRING=VERSION_STRING//"."
        VERSION_STRING=VERSION_STRING//TRIM(NUMBER_TO_VSTRING(CMISS_REVISION_VERSION,"*",ERR,ERROR))
        !VERSION_STRING=VERSION_STRING//" ("
        !VERSION_STRING=VERSION_STRING//TRIM(CMISS_BUILD_VERSION(6:))
        !VERSION_STRING=VERSION_STRING//" )"
        
        WRITE(*,'(A)') CHAR(VERSION_STRING)

      ENDIF
    ENDIF
    
    RETURN
999 RETURN 1
  END SUBROUTINE CMISS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Writes the error string to screen. \todo replace with CMISS_HANDLE_ERROR.
  SUBROUTINE CMISS_WRITE_ERROR(ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local Variables
    
    CALL WRITE_ERROR(ERR,ERROR,*999)

    RETURN
999 RETURN

  END SUBROUTINE CMISS_WRITE_ERROR

  !
  !================================================================================================================================
  !

  !>Handle an error condition
  SUBROUTINE CMISS_HANDLE_ERROR(ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR
    
    SELECT CASE(CMISS_ERROR_HANDLING_MODE)
    CASE(CMISS_RETURN_ERROR_CODE)
      !Do nothing
    CASE(CMISS_OUTPUT_ERROR)
      CALL WRITE_ERROR(ERR,ERROR,*999)
    CASE(CMISS_TRAP_ERROR)
      CALL WRITE_ERROR(ERR,ERROR,*999)
      CALL MPI_ABORT(MPI_COMM_WORLD,ERR,MPI_IERROR)
      STOP
    CASE DEFAULT
      !Do nothing
    END SELECT

    RETURN
999 RETURN

  END SUBROUTINE CMISS_HANDLE_ERROR
  
  !
  !================================================================================================================================
  !

END MODULE CMISS
