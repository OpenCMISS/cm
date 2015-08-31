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
#ifndef NOMPIMOD
  USE MPI
#endif
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE STRINGS
  USE TYPES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

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

  INTEGER(INTG), SAVE :: CmissErrorHandlingMode !<The current error handling mode for OpenCMISS \see CMISS_ErrorHandlingModes
 
  !Interfaces

  INTERFACE

    SUBROUTINE CMISSInitFatalHandler() BIND(C,NAME="CMISSInitFatalHandler")
    END SUBROUTINE CMISSInitFatalHandler

    SUBROUTINE CMISSResetFatalHandler() BIND(C,NAME="CMISSResetFatalHandler")
    END SUBROUTINE CMISSResetFatalHandler
    
    SUBROUTINE CMISSSetFatalHandler() BIND(C,NAME="CMISSSetFatalHandler")
    END SUBROUTINE CMISSSetFatalHandler

  END INTERFACE

  PUBLIC CMISS_MAJOR_VERSION,CMISS_MINOR_VERSION,CMISS_REVISION_VERSION,CMISS_BUILD_VERSION

  PUBLIC CMISS_RETURN_ERROR_CODE,CMISS_OUTPUT_ERROR,CMISS_TRAP_ERROR

  PUBLIC CmissErrorHandlingModeGet_,CmissErrorHandlingModeSet_
  
  PUBLIC CmissHandleError
  
  PUBLIC CmissFinalise_,CmissInitialise_

CONTAINS

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Returns the error handling mode for CMISS \see OPENCMISS::CMISSErrorHandlingModeGet
  SUBROUTINE CmissErrorHandlingModeGet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: errorHandlingMode !<On return, the error handling mode. \see CMISS_ErrorHandlingModes,CMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    ENTERS("CmissErrorHandlingModeGet_",err,error,*999)

    errorHandlingMode=CmissErrorHandlingMode
    
    EXITS("CmissErrorHandlingModeGet_")
    RETURN
999 ERRORSEXITS("",err,error)
    RETURN 1
    
  END SUBROUTINE CmissErrorHandlingModeGet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Sets the error handling mode for cmiss \see OPENCMISS::CMISSErrorHandlingModeSet
  SUBROUTINE CmissErrorHandlingModeSet_(errorHandlingMode,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: errorHandlingMode !<The error handling mode to set. \see CMISS_ErrorHandlingModes,CMISS
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CmissErrorHandlingModeSet_",err,error,*999)

    SELECT CASE(errorHandlingMode)
    CASE(CMISS_RETURN_ERROR_CODE)
      CmissErrorHandlingMode=CMISS_RETURN_ERROR_CODE
    CASE(CMISS_OUTPUT_ERROR)
      CmissErrorHandlingMode=CMISS_OUTPUT_ERROR
    CASE(CMISS_TRAP_ERROR)
      CmissErrorHandlingMode=CMISS_TRAP_ERROR
    CASE DEFAULT
      localError="The supplied error handling mode of "//TRIM(NumberToVString(errorHandlingMode,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("CmissErrorHandlingModeSet_")
    RETURN
999 ERRORSEXITS("CmissErrorHandlingModeSet_",err,error)
    RETURN 1
    
  END SUBROUTINE CmissErrorHandlingModeSet_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.
  
  !>Finalises CMISS. \see OPENCMISS::CMISSFinalise
  SUBROUTINE CmissFinalise_(err,error,*)
  
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
    CALL CMISSResetFatalHandler()
    !Finalise computational enviroment
    CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(err,error,*999)
    !Finalise the base routines
    CALL BASE_ROUTINES_FINALISE(err,error,*999)
     
    RETURN
999 RETURN 1
    
  END SUBROUTINE CmissFinalise_

  !
  !================================================================================================================================
  !

!!TODO Underscore to avoid name clash. Can be removed upon prefix rename.

  !>Initialises CMISS. \see OPENCMISS::CMISSInitialise
  SUBROUTINE CmissInitialise_(worldRegion,err,error,*)
  
    !Argument variables
    TYPE(REGION_TYPE), POINTER :: worldRegion !<On exit, a pointer to the world region. Must not be associated on entry.
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: versionString

    !Initialise error mode
    CmissErrorHandlingMode = CMISS_OUTPUT_ERROR !Default for now, maybe make CMISS_RETURN_ERROR_CODE the default
    !Initialise the base routines
    CALL BASE_ROUTINES_INITIALISE(err,error,*999)
    !Intialise the computational environment
    CALL COMPUTATIONAL_ENVIRONMENT_INITIALISE(err,error,*999)
    !Setup signal handling
    CALL CMISSInitFatalHandler()
    CALL CMISSSetFatalHandler()
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
        versionString="OpenCMISS(cm) version "//TRIM(NumberToVString(CMISS_MAJOR_VERSION,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(CMISS_MINOR_VERSION,"*",err,error))
        versionString=versionString//"."
        versionString=versionString//TRIM(NumberToVString(CMISS_REVISION_VERSION,"*",err,error))
        !versionString=versionString//" ("
        !versionString=versionString//TRIM(CMISS_BUILD_VERSION(6:))
        !versionString=versionString//" )"
        
        WRITE(*,'(A)') CHAR(versionString)

      ENDIF
    ENDIF
    
    RETURN
999 RETURN 1
    
  END SUBROUTINE CmissInitialise_

  !
  !================================================================================================================================
  !

  !>Handle an error condition
  SUBROUTINE CmissHandleError(err,error)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: mpiError
    
    SELECT CASE(CmissErrorHandlingMode)
    CASE(CMISS_RETURN_ERROR_CODE)
      !Do nothing
    CASE(CMISS_OUTPUT_ERROR)
      CALL WriteError(err,error,*999)
    CASE(CMISS_TRAP_ERROR)
      CALL WriteError(err,error,*999)
      CALL MPI_ABORT(MPI_COMM_WORLD,err,mpiError)
      STOP
    CASE DEFAULT
      !Do nothing
    END SELECT

    RETURN
999 RETURN

  END SUBROUTINE CmissHandleError
  
  !
  !================================================================================================================================
  !

END MODULE Cmiss
