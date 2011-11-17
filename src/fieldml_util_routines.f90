!> \file
!> \author Caton Little
!> \brief This module handles non-IO FieldML logic.
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

!> Utility routines for FieldML

MODULE FIELDML_UTIL_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE FIELDML_API
  USE FIELDML_TYPES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE REGION_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Interfaces
  
  INTERFACE FIELDML_UTIL_CHECK_FIELDML_ERROR
    MODULE PROCEDURE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS
    MODULE PROCEDURE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC
  END INTERFACE FIELDML_UTIL_CHECK_FIELDML_ERROR

  PUBLIC :: FIELDML_UTIL_CHECK_FIELDML_ERROR

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS( ERROR_DESCRIPTION, FML_HANDLE, ERR, ERROR, * )
    TYPE(VARYING_STRING), INTENT(IN) :: ERROR_DESCRIPTION
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FML_ERR
    
    CALL ENTERS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS", ERR, ERROR, *999 )
    FML_ERR = Fieldml_GetLastError( FML_HANDLE )

    IF( FML_ERR == FML_ERR_NO_ERROR ) THEN
      CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS" )
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( ERROR_DESCRIPTION // " (error number " // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 CALL ERRORS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORVS
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC( ERROR_DESCRIPTION, FML_HANDLE, ERR, ERROR, * )
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_DESCRIPTION
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FML_ERR
    
    CALL ENTERS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC", ERR, ERROR, *999 )
    FML_ERR = Fieldml_GetLastError( FML_HANDLE )

    IF( FML_ERR == FML_ERR_NO_ERROR ) THEN
      CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC" )
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( ERROR_DESCRIPTION // " (error number " // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 CALL ERRORS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC
  
  !
  !================================================================================================================================
  !

END MODULE FIELDML_UTIL_ROUTINES
