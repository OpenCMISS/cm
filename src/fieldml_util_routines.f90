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

  !Module parameters
  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  !Interfaces
  
  INTERFACE FieldmlUtil_CheckError
    MODULE PROCEDURE FieldmlUtil_CheckLastError
    MODULE PROCEDURE FieldmlUtil_CheckLastInfoError
  END INTERFACE

  PUBLIC :: FieldmlInfoType

  PUBLIC :: FieldmlUtil_GetConnectivityEnsemble, FieldmlUtil_GetGenericType, FieldmlUtil_InitialiseInfo, &
    & FieldmlUtil_GetXiType, FieldmlUtil_GetValueType, FieldmlUtil_FinaliseInfo, FieldmlUtil_ImportHandle, &
    & FieldmlUtil_GetCollapseSuffix, FieldmlUtil_GetTypeArgumentHandle, FieldmlUtil_CheckError, FieldmlUtil_CheckErrorNumber, &
    & FieldmlUtil_Import

CONTAINS

  !
  !================================================================================================================================
  !

  !<Import the named object from the built-in library into the current FieldML document. The local name will be the same as the remote name.
  FUNCTION FieldmlUtil_Import( fmlHandle, remoteName )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle.
    TYPE(VARYING_STRING), INTENT(IN) :: remoteName !<The name of the object to import.

    INTEGER(C_INT) :: FieldmlUtil_Import
    
    !Local variables
    INTEGER(C_INT) :: importIndex
    
    FieldmlUtil_Import = Fieldml_GetObjectByName( fmlHandle, char(remoteName)//NUL )
    IF( FieldmlUtil_Import == FML_INVALID_HANDLE ) THEN
      importIndex = Fieldml_AddImportSource( fmlHandle, &
        & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
      FieldmlUtil_Import = Fieldml_AddImport( fmlHandle, importIndex, char(remoteName)//NUL, char(remoteName)//NUL )
    ENDIF

  END FUNCTION FieldmlUtil_Import
  
  !
  !================================================================================================================================
  !

  !<Import the given FieldML object if it is not already imported or local.
  FUNCTION FieldmlUtil_ImportHandle( fmlHandle, handle )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: handle !<The FieldML object to import.

    INTEGER(C_INT) :: FieldmlUtil_ImportHandle
    
    !Local variables
    INTEGER(C_INT) :: importIndex, localHandle
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: name
    INTEGER(INTG) :: length
    
    FieldmlUtil_ImportHandle = FML_INVALID_HANDLE
    length = Fieldml_CopyObjectDeclaredName( fmlHandle, handle, name, MAXSTRLEN )
    
    IF( Fieldml_IsObjectLocal( fmlHandle, handle ) /= 1 ) THEN
      IF( length > 0 ) THEN
        localHandle = Fieldml_GetObjectByName( fmlHandle, name(1:length)//NUL )
        IF( localHandle == FML_INVALID_HANDLE ) THEN
          importIndex = Fieldml_AddImportSource( fmlHandle, &
            & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
          FieldmlUtil_ImportHandle = Fieldml_AddImport( fmlHandle, importIndex, name(1:length)//NUL, name(1:length)//NUL )
        ELSE IF( localHandle == handle ) THEN
          FieldmlUtil_ImportHandle = handle
        ENDIF
      ENDIF
    ENDIF

  END FUNCTION FieldmlUtil_ImportHandle
  
  !
  !================================================================================================================================
  !
  
  !<Get the argument corresponding to the given type (named *.argument), importing it if needed.
  FUNCTION FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, typeHandle, doImport )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    LOGICAL, INTENT(IN) :: doImport !<If true, import the argument.
    INTEGER(C_INT), INTENT(IN) :: typeHandle !<The type out of whose name the argument name is built.
    
    INTEGER(C_INT) :: FieldmlUtil_GetTypeArgumentHandle

    !Local variables
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: name
    INTEGER(INTG) :: length
    INTEGER(C_INT) :: handle, fmlErr
    TYPE(VARYING_STRING) :: fullName
    
    length = Fieldml_CopyObjectName( fieldmlInfo%fmlHandle, typeHandle, name, MAXSTRLEN )
    IF( length < 1 ) THEN
      length = Fieldml_CopyObjectDeclaredName( fieldmlInfo%fmlHandle, typeHandle, name, MAXSTRLEN )
      FieldmlUtil_GetTypeArgumentHandle = FML_INVALID_HANDLE
      RETURN
    ENDIF

    IF( doImport ) THEN
      fullName = name(1:length)//".argument"
      !Note: Don't need to check the result here, as Fieldml_GetObjectByName will fail if the import didn't work.
      fmlErr = FieldmlUtil_Import( fieldmlInfo%fmlHandle, fullName )
    ENDIF
    
    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, name(1:length)//".argument"//NUL )
    IF( handle == FML_INVALID_HANDLE ) THEN
      FieldmlUtil_GetTypeArgumentHandle = FML_INVALID_HANDLE
      RETURN
    ENDIF
    
    FieldmlUtil_GetTypeArgumentHandle = handle
    
  END FUNCTION FieldmlUtil_GetTypeArgumentHandle
  
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckLastInfoError( errorDescription, fieldmlInfo, err, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    INTEGER(INTG) :: fmlErr
    
    CALL ENTERS( "FieldmlUtil_CheckLastInfoError", err, errorString, *999 )

    fmlErr = Fieldml_GetLastError( fieldmlInfo%fmlHandle )

    IF( fmlErr == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( errorDescription // " (error number " // TRIM(NUMBER_TO_VSTRING(fmlErr,"*",err,errorString)) // ")", &
      & err, errorString, *999 )

999 CALL ERRORS( "FieldmlUtil_CheckLastInfoError", err, errorString )
    CALL EXITS( "FieldmlUtil_CheckLastInfoError" )
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckLastInfoError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckLastError( errorDescription, fmlHandle, err, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    INTEGER(INTG) :: fmlErr
    
    CALL ENTERS( "FieldmlUtil_CheckLastError", err, errorString, *999 )
    fmlErr = Fieldml_GetLastError( fmlHandle )

    IF( fmlErr == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( errorDescription // " (error number " // TRIM(NUMBER_TO_VSTRING(fmlErr,"*",err,errorString)) // ")", &
      & err, errorString, *999 )

999 CALL ERRORS( "FieldmlUtil_CheckLastError", err, errorString )
    CALL EXITS( "FieldmlUtil_CheckLastError" )
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckLastError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckErrorNumber( errorDescription, error, err, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    INTEGER(INTG), INTENT(IN) :: error
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    CALL ENTERS( "FieldmlUtil_CheckErrorNumber", err, errorString, *999 )
    IF( error == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    CALL FLAG_ERROR( errorDescription // " (error number " // TRIM(NUMBER_TO_VSTRING(error,"*",err,errorString)) // ")", &
      & err, errorString, *999 )

999 CALL ERRORS( "FieldmlUtil_CheckErrorNumber", err, errorString )
    CALL EXITS( "FieldmlUtil_CheckErrorNumber" )
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckErrorNumber
  
  !
  !================================================================================================================================
  !
  
  !<Get the FieldML built-in library type corresponding to the given OpenCMISS coordinate system type.
  SUBROUTINE FieldmlUtil_GetCoordinatesType( fieldmlHandle, coordsType, dimensions, doImport, typeHandle, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: coordsType !<The OpenCMISS coordinates type.
    INTEGER(INTG), INTENT(IN) :: dimensions !<The coordinate system's number of dimensions.
    LOGICAL, INTENT(IN) :: doImport !<If true, import the FieldML type.
    INTEGER(C_INT), INTENT(OUT) :: typeHandle !<The FieldML type corresponding to the given OpenCMISS coordinate system type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: typeName
    INTEGER(C_INT) :: temp

    CALL ENTERS( "FieldmlUtil_GetCoordinatesType", err, errorString, *999 )
    
    IF( coordsType == COORDINATE_RECTANGULAR_CARTESIAN_TYPE ) THEN
      IF( dimensions == 1 ) THEN
        typeName = "coordinates.rc.1d"
      ELSE IF( dimensions == 2 ) THEN
        typeName = "coordinates.rc.2d"
      ELSE IF( dimensions == 3 ) THEN
        typeName = "coordinates.rc.3d"
      ELSE
        typeHandle = FML_INVALID_HANDLE
        CALL FLAG_ERROR( "Cannot get RC coordinates type", err, errorString, *999 )
      ENDIF
    ELSE
      typeHandle = FML_INVALID_HANDLE
      CALL FLAG_ERROR( "Cannot get coordinates type", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, typeName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(typeName)//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get coordinates type", fieldmlHandle, err, errorString, *999 )

    CALL EXITS( "FieldmlUtil_GetCoordinatesType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetCoordinatesType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetCoordinatesType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetCoordinatesType
  
  !
  !================================================================================================================================
  !
  
  !<Returns a generic n-dimensional real type from the built-in library.
  SUBROUTINE FieldmlUtil_GetGenericType( fieldmlHandle, dimensions, typeHandle, doImport, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: dimensions !<The number of dimensions of the type.
    INTEGER(C_INT), INTENT(OUT) :: typeHandle !<The FieldML type.
    LOGICAL, INTENT(IN) :: doImport !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: typeName
    INTEGER(C_INT) :: temp

    CALL ENTERS( "FieldmlUtil_GetGenericType", err, errorString, *999 )

    IF( dimensions == 1 ) THEN
      typeName = "real.1d"
    ELSE IF( dimensions == 2 ) THEN
      typeName = "real.2d"
    ELSE IF( dimensions == 3 ) THEN
      typeName = "real.3d"
    ELSE
      typeHandle = FML_INVALID_HANDLE
      CALL FLAG_ERROR( "Cannot get generic type", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, typeName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(typeName)//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get generic type", fieldmlHandle, err, errorString, *999 )
    
    CALL EXITS( "FieldmlUtil_GetGenericType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetGenericType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetGenericType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetGenericType
  
  !
  !================================================================================================================================
  !
  
  !<Returns a type in the built-in library corresponding to a chart of the given dimensionality.
  SUBROUTINE FieldmlUtil_GetXiType( fieldmlHandle, dimensions, doImport, typeHandle, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: dimensions !<The number of dimensions of the chart.
    LOGICAL, INTENT(IN) :: doImport !<If true, import the type.
    INTEGER(C_INT), INTENT(OUT) :: typeHandle !<The FieldML type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    INTEGER(C_INT) :: temp
    TYPE(VARYING_STRING) :: typeName

    CALL ENTERS( "FieldmlUtil_GetXiType", err, errorString, *999 )
    
    IF( dimensions == 1 ) THEN
      typeName = "chart.1d"
    ELSE IF( dimensions == 2 ) THEN
      typeName = "chart.2d"
    ELSE IF( dimensions == 3 ) THEN
      typeName = "chart.3d"
    ELSE
      typeHandle = FML_INVALID_HANDLE
      CALL FLAG_ERROR( "Invalid chart dimensions", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, typeName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(typeName)//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get xi type", fieldmlHandle, err, errorString, *999 )
    
    CALL EXITS( "FieldmlUtil_GetXiType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetXiType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetXiType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetXiType
  
  !
  !================================================================================================================================
  !
  
  !<Get the text suffix corresponding to the given array of collapse constants.
  SUBROUTINE FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err, errorString )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:) !<The collapse into from which to generate the suffix.
    TYPE(VARYING_STRING), INTENT(INOUT) :: suffix !<The suffix string encoding the collapse info.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    INTEGER(INTG) :: i
    
    suffix = ""
    DO i = 1, SIZE( collapseInfo )
      IF( collapseInfo( i ) == BASIS_XI_COLLAPSED ) THEN
        suffix = suffix // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",err,errorString))//"C"
      ELSEIF( collapseInfo( i ) == BASIS_COLLAPSED_AT_XI0 ) THEN
        suffix = suffix // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",err,errorString))//"0"
      ELSEIF( collapseInfo( i ) == BASIS_COLLAPSED_AT_XI1 ) THEN
        suffix = suffix // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",err,errorString))//"1"
      ENDIF
    ENDDO
  
  END SUBROUTINE
  
  !
  !================================================================================================================================
  !
  
  !<Return the FieldML connectivity ensemble corresponding to the given tensor-product basis info.
  SUBROUTINE FieldmlUtil_GetTPConnectivityEnsemble( fieldmlHandle, xiInterpolations, collapseInfo, doImport, typeHandle, &
    & err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle !<The FieldML session handle.
    INTEGER(C_INT), INTENT(IN) :: xiInterpolations(:) !<The per-xi interpolation of the given TP basis.
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:) !<The collapse-constant for the given basis.
    LOGICAL, INTENT(IN) :: doImport !<If true, import the connectivity ensemble.
    INTEGER(C_INT), INTENT(OUT) :: typeHandle !<The FieldML ensemble handle.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    !Locals
    INTEGER(C_INT) :: xiCount, firstInterpolation, i, importIndex, temp
    TYPE(VARYING_STRING) :: suffix, layoutName
    
    CALL ENTERS( "FieldmlUtil_GetTPConnectivityEnsemble", err, errorString, *999 )

    xiCount = SIZE( xiInterpolations )
  
    importIndex = Fieldml_AddImportSource( fieldmlHandle, &
      & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot access library", fieldmlHandle, err, errorString, *999 )

    firstInterpolation = xiInterpolations(1)
    DO i = 2, xiCount
      IF( xiInterpolations(i) /= firstInterpolation ) THEN
        !Do not yet support inhomogeneous TP bases
        CALL FLAG_ERROR( "Inhomogeneous tensor-product bases are not yet supported", err, errorString, *999 )
      ENDIF
    ENDDO

    CALL FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err, errorString )
      
    IF( firstInterpolation == BASIS_QUADRATIC_LAGRANGE_INTERPOLATION ) THEN
      IF( xiCount == 1 ) THEN
        layoutName = "localNodes.1d.line3"
      ELSE IF( xiCount == 2 ) THEN
        layoutName = "localNodes.2d.square3x3"//char(suffix)
      ELSE IF( xiCount == 3 ) THEN
        layoutName = "localNodes.3d.cube3x3x3"//char(suffix)
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( "Quadratic interpolation not support for more than 3 dimensions", &
          & err, errorString, *999 )
      ENDIF
    ELSE IF( firstInterpolation == BASIS_LINEAR_LAGRANGE_INTERPOLATION ) THEN
      IF( xiCount == 1 ) THEN
        layoutName = "localNodes.1d.line2"
      ELSE IF( xiCount == 2 ) THEN
        layoutName = "localNodes.2d.square2x2"//char(suffix)
      ELSE IF( xiCount == 3 ) THEN
        layoutName = "localNodes.3d.cube2x2x2"//char(suffix)
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( "Linear interpolation not support for more than 3 dimensions", &
          & err, errorString, *999 )
      ENDIF
    ELSE
      CALL FLAG_ERROR( "Interpolation not yet supported", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, layoutName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(layoutName)//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get local nodes type", fieldmlHandle, err, errorString, *999 )
    
    CALL EXITS( "FieldmlUtil_GetTPConnectivityEnsemble" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetTPConnectivityEnsemble", err, errorString )
    CALL EXITS( "FieldmlUtil_GetTPConnectivityEnsemble" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetTPConnectivityEnsemble

  !
  !================================================================================================================================
  !

  !<Get the connectivity ensemble for the given basis. Currently, only tensor-product bases are supported.
  SUBROUTINE FieldmlUtil_GetConnectivityEnsemble( fieldmlHandle, basis, typeHandle, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle !<The FieldML session handle.
    TYPE(BASIS_TYPE), POINTER :: basis !<The basis for which to return the connectivity.
    INTEGER(C_INT), INTENT(OUT) :: typeHandle !<The FieldML connectivity ensemble.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    !Locals
    INTEGER(C_INT) :: basisType, xiCount
    INTEGER(C_INT), ALLOCATABLE :: xiInterpolations(:), collapseInfo(:)
    
    CALL ENTERS( "FieldmlUtil_GetConnectivityEnsemble", err, errorString, *999 )

    typeHandle = FML_INVALID_HANDLE

    CALL BASIS_TYPE_GET( basis, basisType, err, errorString, *999 )
    CALL BASIS_NUMBER_OF_XI_GET( basis, xiCount, err, errorString, *999 )
    
    IF( basisType == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      ALLOCATE( xiInterpolations( xiCount ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate xi interpolations array.", err, errorString, *999 )
      ALLOCATE( collapseInfo( xiCount ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse info array.", err, errorString, *999 )
      CALL BASIS_INTERPOLATION_XI_GET( basis, xiInterpolations, err, errorString, *999 )
      CALL BASIS_COLLAPSED_XI_GET( basis, collapseInfo, err, errorString, *999 )
      
      CALL FieldmlUtil_GetTPConnectivityEnsemble( fieldmlHandle, xiInterpolations, collapseInfo, .TRUE., typeHandle, &
        & err, errorString, *999 )
      
      DEALLOCATE( xiInterpolations )
      DEALLOCATE( collapseInfo )
    ELSE
      CALL FLAG_ERROR( "Only tensor product bases are currently supported", err, errorString, *999 )
    ENDIF
    
    CALL EXITS( "FieldmlUtil_GetConnectivityEnsemble" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetConnectivityEnsemble", err, errorString )
    CALL EXITS( "FieldmlUtil_GetConnectivityEnsemble" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetConnectivityEnsemble

  !
  !================================================================================================================================
  !
  
  !<Returns a FieldML type appropriate for the given OpenCMISS field.
  SUBROUTINE FieldmlUtil_GetValueType( fmlHandle, field, typeHandle, doImport, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field whose type is to be obtained.
    INTEGER(C_INT), INTENT(OUT) :: typeHandle !<The FieldML type handle.
    LOGICAL, INTENT(IN) :: doImport !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    !Locals
    INTEGER(INTG) :: fieldType, subType, count
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(REGION_TYPE), POINTER :: region
    
    CALL ENTERS( "FieldmlUtil_GetValueType", err, errorString, *999 )
    
    region => field%REGION
    
    CALL FIELD_TYPE_GET( field, fieldType, err, errorString, *999 )
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( field, FIELD_U_VARIABLE_TYPE, count, err, errorString, *999 )

    SELECT CASE( fieldType )
    CASE( FIELD_GEOMETRIC_TYPE )
      NULLIFY( coordinateSystem )
      CALL REGION_COORDINATE_SYSTEM_GET( region, coordinateSystem, err, errorString, *999 )
      CALL COORDINATE_SYSTEM_TYPE_GET( coordinateSystem, subType, err, errorString, *999 )
      CALL FieldmlUtil_GetCoordinatesType( fmlHandle, subType, count, doImport, typeHandle, err, errorString, *999 )
    
    !CASE( CMISSFieldFibreType )

    !CASE( CMISSFieldGeneralType )

    !CASE( CMISSFieldMaterialType )

    CASE DEFAULT
      CALL FieldmlUtil_GetGenericType( fmlHandle, count, typeHandle, doImport, err, errorString, *999 )
    END SELECT
  
    CALL EXITS( "FieldmlUtil_GetValueType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetValueType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetValueType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetValueType
    
  !
  !================================================================================================================================
  !
  
  !<Initialise up the given FieldML parsing state.
  SUBROUTINE FieldmlUtil_InitialiseInfo( fieldmlInfo, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state to clean up.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    CALL ENTERS( "FieldmlUtil_InitialiseInfo", err, errorString, *999 )

    fieldmlInfo%fmlHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodesHandle = FML_INVALID_HANDLE
    fieldmlInfo%meshHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodeDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    
    NULLIFY( fieldmlInfo%componentHandles )
    CALL LIST_CREATE_START( fieldmlInfo%componentHandles, err, errorString, *999 )
    CALL LIST_DATA_TYPE_SET( fieldmlInfo%componentHandles, LIST_C_INT_TYPE, err, errorString, *999 )
    CALL LIST_MUTABLE_SET( fieldmlInfo%componentHandles, .TRUE., err, errorString, *999 )
    CALL LIST_CREATE_FINISH( fieldmlInfo%componentHandles, err, errorString, *999 )
    
    NULLIFY( fieldmlInfo%basisHandles )
    CALL LIST_CREATE_START( fieldmlInfo%basisHandles, err, errorString, *999 )
    CALL LIST_DATA_TYPE_SET( fieldmlInfo%basisHandles, LIST_C_INT_TYPE, err, errorString, *999 )
    CALL LIST_CREATE_FINISH( fieldmlInfo%basisHandles, err, errorString, *999 )
    
    NULLIFY( fieldmlInfo%basisConnectivityHandles )
    CALL LIST_CREATE_START( fieldmlInfo%basisConnectivityHandles, err, errorString, *999 )
    CALL LIST_DATA_TYPE_SET( fieldmlInfo%basisConnectivityHandles, LIST_C_INT_TYPE, err, errorString, *999 )
    CALL LIST_CREATE_FINISH( fieldmlInfo%basisConnectivityHandles, err, errorString, *999 )
    
    NULLIFY( fieldmlInfo%basisLayoutHandles )
    CALL LIST_CREATE_START( fieldmlInfo%basisLayoutHandles, err, errorString, *999 )
    CALL LIST_DATA_TYPE_SET( fieldmlInfo%basisLayoutHandles, LIST_C_INT_TYPE, err, errorString, *999 )
    CALL LIST_CREATE_FINISH( fieldmlInfo%basisLayoutHandles, err, errorString, *999 )

    CALL EXITS( "FieldmlUtil_InitialiseInfo" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_InitialiseInfo", err, errorString )
    CALL EXITS( "FieldmlUtil_InitialiseInfo" )
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_InitialiseInfo

  !
  !================================================================================================================================
  !

  !<Clean up the given FieldML parsing state.
  SUBROUTINE FieldmlUtil_FinaliseInfo( fieldmlInfo, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state to clean up.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: fmlErr

    CALL ENTERS( "FieldmlUtil_FinaliseInfo", err, errorString, *999 )

    fmlErr = Fieldml_Destroy( fieldmlInfo%fmlHandle )
    
    fieldmlInfo%fmlHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodesHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodesArgumentHandle = FML_INVALID_HANDLE
    fieldmlInfo%meshHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsArgumentHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiArgumentHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodeDofsHandle = FML_INVALID_HANDLE
!    fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
!    fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    
    CALL LIST_DESTROY( fieldmlInfo%componentHandles, err, errorString, *999 )
    CALL LIST_DESTROY( fieldmlInfo%basisHandles, err, errorString, *999 )
    CALL LIST_DESTROY( fieldmlInfo%basisConnectivityHandles, err, errorString, *999 )
    CALL LIST_DESTROY( fieldmlInfo%basisLayoutHandles, err, errorString, *999 )

    CALL EXITS( "FieldmlUtil_FinaliseInfo" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_FinaliseInfo", err, errorString )
    CALL EXITS( "FieldmlUtil_FinaliseInfo" )
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_FinaliseInfo

  !
  !================================================================================================================================
  !

END MODULE FIELDML_UTIL_ROUTINES
