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
    MODULE PROCEDURE FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS
    MODULE PROCEDURE FIELDML_UTIL_CHECK_FIELDML_SESSION_ERRORC
    MODULE PROCEDURE FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC
  END INTERFACE FIELDML_UTIL_CHECK_FIELDML_ERROR

  PUBLIC :: FIELDML_INFO_TYPE

  PUBLIC :: FIELDML_UTIL_INITIALISE_INFO, &
    & FIELDML_UTIL_FINALISE_INFO, FIELDML_UTIL_GET_COLLAPSE_SUFFIX, FIELDML_UTIL_CHECK_FIELDML_ERROR, &
    & FIELDML_UTIL_CHECK_ERROR_NUMBER

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS( ERROR_DESCRIPTION, FIELDML_INFO, ERR, ERROR, * )
    TYPE(VARYING_STRING), INTENT(IN) :: ERROR_DESCRIPTION
    TYPE(FIELDML_INFO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FML_ERR
    
    CALL ENTERS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS", ERR, ERROR, *999 )

    FML_ERR = Fieldml_GetLastError( FIELDML_INFO%FML_HANDLE )

    IF( FML_ERR == FML_ERR_NO_ERROR ) THEN
      CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS" )
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( ERROR_DESCRIPTION // " (error number " // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 CALL ERRORS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORVS
  
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
  
  SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC( ERROR_DESCRIPTION, FIELDML_INFO, ERR, ERROR, * )
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_DESCRIPTION
    TYPE(FIELDML_INFO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FML_ERR
    
    CALL ENTERS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC", ERR, ERROR, *999 )

    FML_ERR = Fieldml_GetLastError( FIELDML_INFO%FML_HANDLE )

    IF( FML_ERR == FML_ERR_NO_ERROR ) THEN
      CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC" )
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( ERROR_DESCRIPTION // " (error number " // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 CALL ERRORS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_CHECK_FIELDML_INFO_ERRORC
  
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
  
  SUBROUTINE FIELDML_UTIL_CHECK_ERROR_NUMBER( ERROR_DESCRIPTION, FML_ERROR, ERR, ERROR, * )
    CHARACTER(LEN=*), INTENT(IN) :: ERROR_DESCRIPTION
    INTEGER(INTG), INTENT(IN) :: FML_ERROR
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS( "FIELDML_UTIL_CHECK_ERROR_NUMBER", ERR, ERROR, *999 )
    IF( FML_ERROR == FML_ERR_NO_ERROR ) THEN
      CALL EXITS( "FIELDML_UTIL_CHECK_ERROR_NUMBER" )
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( ERROR_DESCRIPTION // " (error number " // TRIM(NUMBER_TO_VSTRING(FML_ERROR,"*",ERR,ERROR)) // ")", &
      & ERR, ERROR, *999 )

999 CALL ERRORS( "FIELDML_UTIL_CHECK_ERROR_NUMBER", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_CHECK_ERROR_NUMBER" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_CHECK_ERROR_NUMBER
  
  !
  !================================================================================================================================
  !
  
  !<Get the text suffix corresponding to the given array of collapse constants.
  SUBROUTINE FIELDML_UTIL_GET_COLLAPSE_SUFFIX( COLLAPSE_INFO, SUFFIX, ERR, ERROR )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: COLLAPSE_INFO(:) !<The collapse into from which to generate the suffix.
    TYPE(VARYING_STRING), INTENT(INOUT) :: SUFFIX !<The suffix string encoding the collapse info.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    INTEGER(INTG) :: I
    
    SUFFIX = ""
    DO I = 1, SIZE( COLLAPSE_INFO )
      IF( COLLAPSE_INFO( I ) == BASIS_XI_COLLAPSED ) THEN
        SUFFIX = SUFFIX // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//"C"
      ELSEIF( COLLAPSE_INFO( I ) == BASIS_COLLAPSED_AT_XI0 ) THEN
        SUFFIX = SUFFIX // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//"0"
      ELSEIF( COLLAPSE_INFO( I ) == BASIS_COLLAPSED_AT_XI1 ) THEN
        SUFFIX = SUFFIX // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//"1"
      ENDIF
    ENDDO
  
  END SUBROUTINE
  
  !
  !================================================================================================================================
  !
  
  !<Initialise up the given FieldML parsing state.
  SUBROUTINE FIELDML_UTIL_INITIALISE_INFO( FIELDML_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_INFO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state to clean up.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS( "FIELDML_UTIL_INITIALISE_INFO", ERR, ERROR, *999 )

    FIELDML_INFO%FML_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODES_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%MESH_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%ELEMENTS_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%XI_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODE_DOFS_HANDLE = FML_INVALID_HANDLE
    !fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    
    NULLIFY( FIELDML_INFO%COMPONENT_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%COMPONENT_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_MUTABLE_SET( FIELDML_INFO%COMPONENT_HANDLES, .TRUE., ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
    
    NULLIFY( FIELDML_INFO%BASIS_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    
    NULLIFY( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    
    NULLIFY( FIELDML_INFO%BASIS_LAYOUT_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_LAYOUT_HANDLES, LIST_INTG_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_UTIL_INITIALISE_INFO" )
    RETURN
999 CALL ERRORS( "FIELDML_UTIL_INITIALISE_INFO", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_INITIALISE_INFO" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_INITIALISE_INFO

  !
  !================================================================================================================================
  !

  !<Clean up the given FieldML parsing state.
  SUBROUTINE FIELDML_UTIL_FINALISE_INFO( FIELDML_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_INFO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state to clean up.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: FML_ERR

    CALL ENTERS( "FIELDML_UTIL_FINALISE_INFO", ERR, ERROR, *999 )

    FML_ERR = Fieldml_Destroy( FIELDML_INFO%FML_HANDLE )
    
    FIELDML_INFO%FML_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODES_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODES_ARGUMENT_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%MESH_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%ELEMENTS_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%ELEMENTS_ARGUMENT_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%XI_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%XI_ARGUMENT_HANDLE = FML_INVALID_HANDLE
    FIELDML_INFO%NODE_DOFS_HANDLE = FML_INVALID_HANDLE
!    fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
!    fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    
    CALL LIST_DESTROY( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DESTROY( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DESTROY( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DESTROY( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_UTIL_FINALISE_INFO" )
    RETURN
999 CALL ERRORS( "FIELDML_UTIL_FINALISE_INFO", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_FINALISE_INFO" )
    RETURN 1
    
  END SUBROUTINE FIELDML_UTIL_FINALISE_INFO

  !
  !================================================================================================================================
  !

END MODULE FIELDML_UTIL_ROUTINES
