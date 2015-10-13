!> \file
!> \author Chris Bradley
!> \brief The top level OpenCMISS Iron module.
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
MODULE Cmiss

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

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(INTG), PARAMETER :: CMFE_MAJOR_VERSION = 0
  INTEGER(INTG), PARAMETER :: CMFE_MINOR_VERSION = 3
  INTEGER(INTG), PARAMETER :: CMFE_REVISION_VERSION = 0

  CHARACTER(LEN=MAXSTRLEN), PARAMETER :: CMFE_BUILD_VERSION = "$Rev"

  !> \addtogroup CMFE_ErrorHandlingModes CMISS::ErrorHandlingModes
  !> \brief Error handling mode parameters
  !> \see CMISS
  !>@{
  INTEGER(INTG), PARAMETER :: CMFE_RETURN_ERROR_CODE = 0 !<Just return the error code \see CMFE_ErrorHandlingModes,CMISS
  INTEGER(INTG), PARAMETER :: CMFE_OUTPUT_ERROR = 1 !<Output the error traceback and return the error code \see CMFE_ErrorHandlingModes,CMISS
  INTEGER(INTG), PARAMETER :: CMFE_TRAP_ERROR = 2 !<Trap the error by outputing the error traceback and stopping the program \see CMFE_ErrorHandlingModes,CMISS
  !>@}
  
  !Module types

  !Module variables

  INTEGER(INTG), SAVE :: cmfe_ErrorHandlingMode !<The current error handling mode for OpenCMISS \see CMFE_ErrorHandlingModes
 
  !Interfaces

  INTERFACE

    SUBROUTINE cmfe_InitFatalHandler() BIND(C,NAME="cmfe_InitFatalHandler")
    END SUBROUTINE cmfe_InitFatalHandler

    SUBROUTINE cmfe_ResetFatalHandler() BIND(C,NAME="cmfe_ResetFatalHandler")
    END SUBROUTINE cmfe_ResetFatalHandler
    
    SUBROUTINE cmfe_SetFatalHandler() BIND(C,NAME="cmfe_SetFatalHandler")
    END SUBROUTINE cmfe_SetFatalHandler

  END INTERFACE

  PUBLIC CMFE_MAJOR_VERSION,CMFE_MINOR_VERSION,CMFE_REVISION_VERSION,CMFE_BUILD_VERSION

  PUBLIC CMFE_RETURN_ERROR_CODE,CMFE_OUTPUT_ERROR,CMFE_TRAP_ERROR

  PUBLIC cmfe_ErrorHandlingModeGet_,cmfe_ErrorHandlingModeSet_
  
  PUBLIC cmfe_HandleError
  
  PUBLIC cmfe_Finalise_,cmfe_Initialise_

CONTAINS

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Returns the error handling mode for CMISS \see OPENCMISS::CMISSErrorHandlingModeGet
  SUBROUTINE cmfe_ErrorHandlingModeGet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: errorHandlingMode !<On return, the error handling mode. \see CMFE_ErrorHandlingModes,CMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    ENTERS("cmfe_ErrorHandlingModeGet_",err,error,*999)

    errorHandlingMode=cmfe_ErrorHandlingMode
    
    EXITS("cmfe_ErrorHandlingModeGet_")
    RETURN
999 ERRORSEXITS("",err,error)
    RETURN 1
    
  END SUBROUTINE cmfe_ErrorHandlingModeGet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Sets the error handling mode for cmiss \see OPENCMISS::CMISSErrorHandlingModeSet
  SUBROUTINE cmfe_ErrorHandlingModeSet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: errorHandlingMode !<The error handling mode to set. \see CMFE_ErrorHandlingModes,CMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("cmfe_ErrorHandlingModeSet",err,error,*999)

    SELECT CASE(errorHandlingMode)
    CASE(CMFE_RETURN_ERROR_CODE)
      cmfe_ErrorHandlingMode=CMFE_RETURN_ERROR_CODE
    CASE(CMFE_OUTPUT_ERROR)
      cmfe_ErrorHandlingMode=CMFE_OUTPUT_ERROR
    CASE(CMFE_TRAP_ERROR)
      cmfe_ErrorHandlingMode=CMFE_TRAP_ERROR
    CASE DEFAULT
      localError="The supplied error handling mode of "//TRIM(NumberToVString(errorHandlingMode,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("cmfe_ErrorHandlingModeSet_")
    RETURN
999 ERRORSEXITS("cmfe_ErrorHandlingModeSet_",err,error)
    RETURN 1
    
  END SUBROUTINE cmfe_ErrorHandlingModeSet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.
  
  !>Finalises CMISS. \see OPENCMISS::CMISSFinalise
  SUBROUTINE cmfe_Finalise_(err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    !Finalise the problems
    CALL PROBLEMS_FINALISE(err,error,*999)
    !Finalise the regions
    CALL REGIONS_FINALISE(err,error,*999)
    !Finalise the coordinate systems
    CALL COORDINATE_SYSTEMS_FINALISE(err,error,*999)
    !Finalise bases
    CALL BASES_FINALISE(err,error,*999)
    !Reset the signal handler
    CALL cmfe_ResetFatalHandler()
    !Finalise computational enviroment
    CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(err,error,*999)
    !Finalise the base routines
    CALL BASE_ROUTINES_FINALISE(err,error,*999)
     
    RETURN
999 RETURN 1
    
  END SUBROUTINE cmfe_Finalise_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Initialises CMISS. \see OPENCMISS::CMISSInitialise
  SUBROUTINE cmfe_Initialise_(worldRegion,err,error,*)
  
    !Argument variables
    TYPE(REGION_TYPE), POINTER :: worldRegion !<On exit, a pointer to the world region. Must not be associated on entry.
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: versionString

    !Initialise error mode
    cmfe_ErrorHandlingMode = CMFE_OUTPUT_ERROR !Default for now, maybe make CMFE_RETURN_ERROR_CODE the default
    !Initialise the base routines
    CALL BASE_ROUTINES_INITIALISE(err,error,*999)
    !Intialise the computational environment
    CALL COMPUTATIONAL_ENVIRONMENT_INITIALISE(err,error,*999)
    !Setup signal handling
    CALL cmfe_InitFatalHandler()
    CALL cmfe_SetFatalHandler()
    IF(ASSOCIATED(worldRegion)) THEN
      CALL FlagError("World region is already associated.",err,error,*999)
    ELSE
      !Intialise the bases
      CALL BASES_INITIALISE(err,error,*999)
      !Initialise the coordinate systems
      CALL COORDINATE_SYSTEMS_INITIALISE(err,error,*999)
      !Initialise the regions 
      CALL REGIONS_INITIALISE(worldRegion,err,error,*999)
      !Initialise the problems
      CALL PROBLEMS_INITIALISE(err,error,*999)
      
      !Write out the CMISS version
      IF(COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
        versionString="OpenCMISS(Iron) version "//TRIM(NumberToVString(CMFE_MAJOR_VERSION,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(CMFE_MINOR_VERSION,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(CMFE_REVISION_VERSION,"*",err,error))
        !versionString=versionString//" ("
        !versionString=versionString//TRIM(CMFE_BUILD_VERSION(6:))
        !versionString=versionString//" )"
        
        WRITE(*,'(A)') CHAR(versionString)

      ENDIF
    ENDIF
    
    RETURN
999 RETURN 1
    
  END SUBROUTINE cmfe_Initialise_

  !
  !================================================================================================================================
  !

  !>Handle an error condition
  SUBROUTINE cmfe_HandleError(err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiError
    
    SELECT CASE(cmfe_ErrorHandlingMode)
    CASE(CMFE_RETURN_ERROR_CODE)
      !Do nothing
    CASE(CMFE_OUTPUT_ERROR)
      CALL WriteError(err,error,*999)
    CASE(CMFE_TRAP_ERROR)
      CALL WriteError(err,error,*999)
      CALL MPI_ABORT(MPI_COMM_WORLD,err,mpiError)
      STOP
    CASE DEFAULT
      !Do nothing
    END SELECT

    RETURN
999 RETURN

  END SUBROUTINE cmfe_HandleError
  
  !
  !================================================================================================================================
  !

END MODULE Cmiss
