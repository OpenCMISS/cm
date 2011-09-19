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

  PUBLIC :: FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE, FIELDML_UTIL_GET_GENERIC_TYPE, FIELDML_UTIL_INITIALISE_INFO, &
    & FIELDML_UTIL_GET_XI_TYPE, FIELDML_UTIL_GET_VALUE_TYPE, FIELDML_UTIL_FINALISE_INFO, FIELDML_UTIL_IMPORT_HANDLE, &
    & FIELDML_UTIL_GET_COLLAPSE_SUFFIX, FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE, FIELDML_UTIL_CHECK_FIELDML_ERROR, &
    & FIELDML_UTIL_CHECK_ERROR_NUMBER, FIELDML_UTIL_IMPORT

CONTAINS

  !
  !================================================================================================================================
  !

  !<Import the named object from the built-in library into the current FieldML document. The local name will be the same as the remote name.
  FUNCTION FIELDML_UTIL_IMPORT( FML_HANDLE, REMOTE_NAME )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    TYPE(VARYING_STRING), INTENT(IN) :: REMOTE_NAME !<The name of the object to import.

    INTEGER(C_INT) :: FIELDML_UTIL_IMPORT
    
    !Local variables
    INTEGER(C_INT) :: IMPORT_INDEX
    
    FIELDML_UTIL_IMPORT = Fieldml_GetObjectByName( FML_HANDLE, cchar(REMOTE_NAME) )
    IF( FIELDML_UTIL_IMPORT == FML_INVALID_HANDLE ) THEN
      IMPORT_INDEX = Fieldml_AddImportSource( FML_HANDLE, &
        & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//C_NULL_CHAR, "library"//C_NULL_CHAR )
      FIELDML_UTIL_IMPORT = Fieldml_AddImport( FML_HANDLE, IMPORT_INDEX, cchar(REMOTE_NAME), cchar(REMOTE_NAME) )
    ENDIF

  END FUNCTION FIELDML_UTIL_IMPORT
  
  !
  !================================================================================================================================
  !

  !<Import the given FieldML object if it is not already imported or local.
  FUNCTION FIELDML_UTIL_IMPORT_HANDLE( FML_HANDLE, HANDLE )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: HANDLE !<The FieldML object to import.

    INTEGER(C_INT) :: FIELDML_UTIL_IMPORT_HANDLE
    
    !Local variables
    INTEGER(C_INT) :: IMPORT_INDEX, LOCAL_HANDLE
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: NAME
    INTEGER(INTG) :: LENGTH
    
    FIELDML_UTIL_IMPORT_HANDLE = FML_INVALID_HANDLE
    LENGTH = Fieldml_CopyObjectDeclaredName( FML_HANDLE, HANDLE, NAME, MAXSTRLEN )
    
    IF( Fieldml_IsObjectLocal( FML_HANDLE, HANDLE ) /= 1 ) THEN
      IF( LENGTH > 0 ) THEN
        LOCAL_HANDLE = Fieldml_GetObjectByName( FML_HANDLE, NAME(1:LENGTH)//C_NULL_CHAR )
        IF( LOCAL_HANDLE == FML_INVALID_HANDLE ) THEN
          IMPORT_INDEX = Fieldml_AddImportSource( FML_HANDLE, &
            & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//C_NULL_CHAR, "library"//C_NULL_CHAR )
          FIELDML_UTIL_IMPORT_HANDLE = Fieldml_AddImport( FML_HANDLE, &
            & IMPORT_INDEX, NAME(1:LENGTH)//C_NULL_CHAR, NAME(1:LENGTH)//C_NULL_CHAR )
        ELSE IF( LOCAL_HANDLE == HANDLE ) THEN
          FIELDML_UTIL_IMPORT_HANDLE = HANDLE
        ENDIF
      ENDIF
    ENDIF

  END FUNCTION FIELDML_UTIL_IMPORT_HANDLE
  
  !
  !================================================================================================================================
  !
  
  !<Get the argument corresponding to the given type (named *.argument), importing it if needed.
  FUNCTION FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, TYPE_HANDLE, DO_IMPORT )
    !Argument variables
    TYPE(FIELDML_INFO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the argument.
    INTEGER(C_INT), INTENT(IN) :: TYPE_HANDLE !<The type out of whose name the argument name is built.
    
    INTEGER(C_INT) :: FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE

    !Local variables
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: NAME
    INTEGER(INTG) :: LENGTH
    INTEGER(C_INT) :: HANDLE, FML_ERR
    TYPE(VARYING_STRING) :: FULL_NAME
    
    LENGTH = Fieldml_CopyObjectName( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE, NAME, MAXSTRLEN )
    IF( LENGTH < 1 ) THEN
      LENGTH = Fieldml_CopyObjectDeclaredName( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE, NAME, MAXSTRLEN )
      FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE = FML_INVALID_HANDLE
      RETURN
    ENDIF

    IF( DO_IMPORT ) THEN
      FULL_NAME = NAME(1:LENGTH)//".argument"
      !Note: Don't need to check the result here, as Fieldml_GetObjectByName will fail if the import didn't work.
      FML_ERR = FIELDML_UTIL_IMPORT( FIELDML_INFO%FML_HANDLE, FULL_NAME )
    ENDIF
    
    HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, NAME(1:LENGTH)//".argument"//C_NULL_CHAR )
    IF( HANDLE == FML_INVALID_HANDLE ) THEN
      FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE = FML_INVALID_HANDLE
      RETURN
    ENDIF
    
    FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE = HANDLE
    
  END FUNCTION FIELDML_UTIL_GET_TYPE_ARGUMENT_HANDLE
  
  
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
    INTEGER(C_INT), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
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
    INTEGER(C_INT), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
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
  
  !<Get the FieldML built-in library type corresponding to the given OpenCMISS coordinate system type.
  SUBROUTINE FieldmlUtilGetCoordinatesType( FIELDML_HANDLE, COORDS_TYPE, DIMENSIONS, DO_IMPORT, TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: COORDS_TYPE !<The OpenCMISS coordinates type.
    INTEGER(INTG), INTENT(IN) :: DIMENSIONS !<The coordinate system's number of dimensions.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the FieldML type.
    INTEGER(C_INT), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type corresponding to the given OpenCMISS coordinate system type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: TYPE_NAME
    INTEGER(C_INT) :: TEMP

    CALL ENTERS( "FieldmlUtilGetCoordinatesType", ERR, ERROR, *999 )
    
    IF( COORDS_TYPE == COORDINATE_RECTANGULAR_CARTESIAN_TYPE ) THEN
      IF( DIMENSIONS == 1 ) THEN
        TYPE_NAME = "coordinates.rc.1d"
      ELSE IF( DIMENSIONS == 2 ) THEN
        TYPE_NAME = "coordinates.rc.2d"
      ELSE IF( DIMENSIONS == 3 ) THEN
        TYPE_NAME = "coordinates.rc.3d"
      ELSE
        TYPE_HANDLE = FML_INVALID_HANDLE
        CALL FLAG_ERROR( var_str("Cannot get FieldML RC coordinates type of dimension ")//DIMENSIONS//".", ERR, ERROR, *999)
      ENDIF
    ELSE
      TYPE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( var_str("Cannot get FieldML coordinates for OpenCMISS type ")//COORDS_TYPE//".", ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_UTIL_IMPORT( FIELDML_HANDLE, TYPE_NAME )
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(TYPE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get FieldML coordinates type "//char(TYPE_NAME)//".", FIELDML_HANDLE, &
      & ERR, ERROR, *999 )

    CALL EXITS( "FieldmlUtilGetCoordinatesType" )
    RETURN
999 CALL ERRORS( "FieldmlUtilGetCoordinatesType", ERR, ERROR )
    CALL EXITS( "FieldmlUtilGetCoordinatesType" )
    RETURN 1

  END SUBROUTINE FieldmlUtilGetCoordinatesType
  
  !
  !================================================================================================================================
  !
  
  !<Returns a generic n-dimensional real type from the built-in library.
  SUBROUTINE FIELDML_UTIL_GET_GENERIC_TYPE( FIELDML_HANDLE, DIMENSIONS, TYPE_HANDLE, DO_IMPORT, ERR, ERROR, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: DIMENSIONS !<The number of dimensions of the type.
    INTEGER(C_INT), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: TYPE_NAME
    INTEGER(C_INT) :: TEMP

    CALL ENTERS( "FIELDML_UTIL_GET_GENERIC_TYPE", ERR, ERROR, *999 )

    IF( DIMENSIONS == 1 ) THEN
      TYPE_NAME = "real.1d"
    ELSE IF( DIMENSIONS == 2 ) THEN
      TYPE_NAME = "real.2d"
    ELSE IF( DIMENSIONS == 3 ) THEN
      TYPE_NAME = "real.3d"
    ELSE
      TYPE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( var_str("Cannot get FieldML generic type of dimensionality ")//DIMENSIONS//".", ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_UTIL_IMPORT( FIELDML_HANDLE, TYPE_NAME )
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(TYPE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get generic type "//TYPE_NAME//".", FIELDML_HANDLE, ERR, ERROR, *999 )
    
    CALL EXITS( "FIELDML_UTIL_GET_GENERIC_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_UTIL_GET_GENERIC_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_GET_GENERIC_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_UTIL_GET_GENERIC_TYPE
  
  !
  !================================================================================================================================
  !
  
  !<Returns a type in the built-in library corresponding to a chart of the given dimensionality.
  SUBROUTINE FIELDML_UTIL_GET_XI_TYPE( FIELDML_HANDLE, DIMENSIONS, DO_IMPORT, TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: DIMENSIONS !<The number of dimensions of the chart.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the type.
    INTEGER(C_INT), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    INTEGER(C_INT) :: TEMP
    TYPE(VARYING_STRING) :: TYPE_NAME

    CALL ENTERS( "FIELDML_UTIL_GET_XI_TYPE", ERR, ERROR, *999 )
    
    IF( DIMENSIONS == 1 ) THEN
      TYPE_NAME = "chart.1d"
    ELSE IF( DIMENSIONS == 2 ) THEN
      TYPE_NAME = "chart.2d"
    ELSE IF( DIMENSIONS == 3 ) THEN
      TYPE_NAME = "chart.3d"
    ELSE
      TYPE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( var_str("Chart dimensionality ")//DIMENSIONS//" not supported.", ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_UTIL_IMPORT( FIELDML_HANDLE, TYPE_NAME )
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(TYPE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get xi type "//TYPE_NAME//".", FIELDML_HANDLE, ERR, ERROR, *999 )
    
    CALL EXITS( "FIELDML_UTIL_GET_XI_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_UTIL_GET_XI_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_GET_XI_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_UTIL_GET_XI_TYPE
  
  !
  !================================================================================================================================
  !
  
  !<Get the text suffix corresponding to the given array of collapse constants.
  SUBROUTINE FIELDML_UTIL_GET_COLLAPSE_SUFFIX( COLLAPSE_INFO, SUFFIX, ERR, ERROR )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: COLLAPSE_INFO(:) !<The collapse into from which to generate the suffix.
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
  
  !<Return the FieldML connectivity ensemble corresponding to the given tensor-product basis info.
  SUBROUTINE FieldmlUtilGetTPConnectivityEnsemble( FIELDML_HANDLE, XI_INTERPOLATIONS, COLLAPSE_INFO, DO_IMPORT, TYPE_HANDLE, &
    & ERR, ERROR, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: XI_INTERPOLATIONS(:) !<The per-xi interpolation of the given TP basis.
    INTEGER(C_INT), INTENT(IN) :: COLLAPSE_INFO(:) !<The collapse-constant for the given basis.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the connectivity ensemble.
    INTEGER(C_INT), INTENT(OUT) :: TYPE_HANDLE !<The FieldML ensemble handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    INTEGER(C_INT) :: XI_COUNT, FIRST_INTERPOLATION, I, IMPORT_INDEX, TEMP
    TYPE(VARYING_STRING) :: SUFFIX, LAYOUT_NAME
    
    CALL ENTERS( "FieldmlUtilGetTPConnectivityEnsemble", ERR, ERROR, *999 )

    XI_COUNT = SIZE( XI_INTERPOLATIONS )
  
    IMPORT_INDEX = Fieldml_AddImportSource( FIELDML_HANDLE, &
      & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//C_NULL_CHAR, "library"//C_NULL_CHAR )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot access built-in FieldML library.", FIELDML_HANDLE, ERR, ERROR, *999 )

    FIRST_INTERPOLATION = XI_INTERPOLATIONS(1)
    DO I = 2, XI_COUNT
      IF( XI_INTERPOLATIONS(I) /= FIRST_INTERPOLATION ) THEN
        !Do not yet support inhomogeneous TP bases
        CALL FLAG_ERROR( "FieldML translation of inhomogeneous tensor-product bases are not yet supported.", &
          & ERR, ERROR, *999 )
      ENDIF
    ENDDO

    CALL FIELDML_UTIL_GET_COLLAPSE_SUFFIX( COLLAPSE_INFO, SUFFIX, ERR, ERROR )
      
    IF( FIRST_INTERPOLATION == BASIS_QUADRATIC_LAGRANGE_INTERPOLATION ) THEN
      IF( XI_COUNT == 1 ) THEN
        LAYOUT_NAME = "localNodes.1d.line3"
      ELSE IF( XI_COUNT == 2 ) THEN
        LAYOUT_NAME = "localNodes.2d.square3x3"//SUFFIX
      ELSE IF( XI_COUNT == 3 ) THEN
        LAYOUT_NAME = "localNodes.3d.cube3x3x3"//SUFFIX
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( var_str("Quadratic Lagrangian interpolation not supported for ")//XI_COUNT//" dimensions.", &
          & ERR, ERROR, *999 )
      ENDIF
    ELSE IF( FIRST_INTERPOLATION == BASIS_LINEAR_LAGRANGE_INTERPOLATION ) THEN
      IF( XI_COUNT == 1 ) THEN
        LAYOUT_NAME = "localNodes.1d.line2"
      ELSE IF( XI_COUNT == 2 ) THEN
        LAYOUT_NAME = "localNodes.2d.square2x2"//SUFFIX
      ELSE IF( XI_COUNT == 3 ) THEN
        LAYOUT_NAME = "localNodes.3d.cube2x2x2"//SUFFIX
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( var_str("Linear Lagrangian interpolation not supported for ")//XI_COUNT//" dimensions.", &
          & ERR, ERROR, *999 )
      ENDIF
    ELSE
      CALL FLAG_ERROR( var_str("FieldML translation not yet supported for interpolation type ")//FIRST_INTERPOLATION//".", &
        & ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_UTIL_IMPORT( FIELDML_HANDLE, LAYOUT_NAME )
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(LAYOUT_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get local nodes type "//LAYOUT_NAME//".", FIELDML_HANDLE, ERR, ERROR, *999 )
    
    CALL EXITS( "FieldmlUtilGetTPConnectivityEnsemble" )
    RETURN
999 CALL ERRORS( "FieldmlUtilGetTPConnectivityEnsemble", ERR, ERROR )
    CALL EXITS( "FieldmlUtilGetTPConnectivityEnsemble" )
    RETURN 1

  END SUBROUTINE FieldmlUtilGetTPConnectivityEnsemble

  !
  !================================================================================================================================
  !

  !<Get the connectivity ensemble for the given basis. Currently, only tensor-product bases are supported.
  SUBROUTINE FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE( FIELDML_HANDLE, BASIS, TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The basis for which to return the connectivity.
    INTEGER(C_INT), INTENT(OUT) :: TYPE_HANDLE !<The FieldML connectivity ensemble.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    INTEGER(C_INT) :: BASISTYPE, XI_COUNT
    INTEGER(C_INT), ALLOCATABLE :: XI_INTERPOLATIONS(:), COLLAPSE_INFO(:)
    
    CALL ENTERS( "FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE", ERR, ERROR, *999 )

    TYPE_HANDLE = FML_INVALID_HANDLE

    CALL BASIS_TYPE_GET( BASIS, BASISTYPE, ERR, ERROR, *999 )
    CALL BASIS_NUMBER_OF_XI_GET( BASIS, XI_COUNT, ERR, ERROR, *999 )
    
    IF( BASISTYPE == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      ALLOCATE( XI_INTERPOLATIONS( XI_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate xi interpolations array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE_INFO( XI_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse info array.", ERR, ERROR, *999 )
      CALL BASIS_INTERPOLATION_XI_GET( BASIS, XI_INTERPOLATIONS, ERR, ERROR, *999 )
      CALL BASIS_COLLAPSED_XI_GET( BASIS, COLLAPSE_INFO, ERR, ERROR, *999 )
      
      CALL FieldmlUtilGetTPConnectivityEnsemble( FIELDML_HANDLE, XI_INTERPOLATIONS, COLLAPSE_INFO, .TRUE., TYPE_HANDLE, &
        & ERR, ERROR, *999 )
      
      DEALLOCATE( XI_INTERPOLATIONS )
      DEALLOCATE( COLLAPSE_INFO )
    ELSE
      CALL FLAG_ERROR( "Only translation of tensor product bases are currently supported", ERR, ERROR, *999 )
    ENDIF
    
    CALL EXITS( "FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE" )
    RETURN
999 CALL ERRORS( "FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE" )
    RETURN 1

  END SUBROUTINE FIELDML_UTIL_GET_CONNECTIVITY_ENSEMBLE

  !
  !================================================================================================================================
  !
  
  !<Returns a FieldML type appropriate for the given OpenCMISS field.
  SUBROUTINE FIELDML_UTIL_GET_VALUE_TYPE( FML_HANDLE, FIELD, TYPE_HANDLE, DO_IMPORT, ERR, ERROR, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field whose type is to be obtained.
    INTEGER(C_INT), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type handle.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    INTEGER(INTG) :: FIELDTYPE, SUB_TYPE, COUNT
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(REGION_TYPE), POINTER :: REGION
    
    CALL ENTERS( "FIELDML_UTIL_GET_VALUE_TYPE", ERR, ERROR, *999 )
    
    REGION => FIELD%REGION
    
    CALL FIELD_TYPE_GET( FIELD, FIELDTYPE, ERR, ERROR, *999 )
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( FIELD, FIELD_U_VARIABLE_TYPE, COUNT, ERR, ERROR, *999 )

    SELECT CASE( FIELDTYPE )
    CASE( FIELD_GEOMETRIC_TYPE )
      NULLIFY( COORDINATE_SYSTEM )
      CALL REGION_COORDINATE_SYSTEM_GET( REGION, COORDINATE_SYSTEM, ERR, ERROR, *999 )
      CALL COORDINATE_SYSTEM_TYPE_GET( COORDINATE_SYSTEM, SUB_TYPE, ERR, ERROR, *999 )
      CALL FieldmlUtilGetCoordinatesType( FML_HANDLE, SUB_TYPE, COUNT, DO_IMPORT, TYPE_HANDLE, ERR, ERROR, *999 )
    
    !CASE( CMISSFieldFibreType )

    !CASE( CMISSFieldGeneralType )

    !CASE( CMISSFieldMaterialType )

    CASE DEFAULT
      CALL FIELDML_UTIL_GET_GENERIC_TYPE( FML_HANDLE, COUNT, TYPE_HANDLE, DO_IMPORT, ERR, ERROR, *999 )
    END SELECT
  
    CALL EXITS( "FIELDML_UTIL_GET_VALUE_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_UTIL_GET_VALUE_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_UTIL_GET_VALUE_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_UTIL_GET_VALUE_TYPE
    
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
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%COMPONENT_HANDLES, LIST_C_INT_TYPE, ERR, ERROR, *999 )
    CALL LIST_MUTABLE_SET( FIELDML_INFO%COMPONENT_HANDLES, .TRUE., ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%COMPONENT_HANDLES, ERR, ERROR, *999 )
    
    NULLIFY( FIELDML_INFO%BASIS_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_HANDLES, LIST_C_INT_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_HANDLES, ERR, ERROR, *999 )
    
    NULLIFY( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, LIST_C_INT_TYPE, ERR, ERROR, *999 )
    CALL LIST_CREATE_FINISH( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, ERR, ERROR, *999 )
    
    NULLIFY( FIELDML_INFO%BASIS_LAYOUT_HANDLES )
    CALL LIST_CREATE_START( FIELDML_INFO%BASIS_LAYOUT_HANDLES, ERR, ERROR, *999 )
    CALL LIST_DATA_TYPE_SET( FIELDML_INFO%BASIS_LAYOUT_HANDLES, LIST_C_INT_TYPE, ERR, ERROR, *999 )
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
    INTEGER(C_INT) :: FML_ERR

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
