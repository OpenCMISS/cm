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

  USE KINDS
  USE FIELDML_API
  USE FIELDML_TYPES
  USE ISO_VARYING_STRING
  USE STRINGS
  USE BASE_ROUTINES
  USE REGION_ROUTINES
  USE COORDINATE_ROUTINES
  USE BASIS_ROUTINES
  USE FIELD_ROUTINES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: BUFFER_SIZE = 1024

  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  !Interfaces
  
  INTERFACE FieldmlUtil_CheckError
    MODULE PROCEDURE FieldmlUtil_CheckLastError
    MODULE PROCEDURE FieldmlUtil_CheckLastInfoError
  END INTERFACE

  PUBLIC :: FieldmlInfoType

  PUBLIC :: FieldmlUtil_GetConnectivityEnsemble, FieldmlUtil_GetGenericType, &
    & FieldmlUtil_GetXiType, FieldmlUtil_GetValueType, FieldmlUtil_FinaliseInfo, FieldmlUtil_ImportHandle, &
    & FieldmlUtil_GetCollapseSuffix, FieldmlUtil_GetTypeArgumentHandle, FieldmlUtil_CheckError, FieldmlUtil_CheckErrorNumber, &
    & FieldmlUtil_Import

CONTAINS

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlUtil_Import( fmlHandle, remoteName )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    TYPE(VARYING_STRING), INTENT(IN) :: remoteName

    INTEGER(C_INT) :: FieldmlUtil_Import
    
    !Local variables
    INTEGER(C_INT) :: importIndex
    
    FieldmlUtil_Import = Fieldml_GetObjectByName( fmlHandle, char(remoteName ) )
    IF( FieldmlUtil_Import == FML_INVALID_HANDLE ) THEN
      importIndex = Fieldml_AddImportSource( fmlHandle, &
        & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
      FieldmlUtil_Import = Fieldml_AddImport( fmlHandle, importIndex, char(remoteName), char(remoteName) )
    ENDIF

  END FUNCTION FieldmlUtil_Import
  
  !
  !================================================================================================================================
  !

  FUNCTION FieldmlUtil_ImportHandle( fmlHandle, handle )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: handle

    INTEGER(C_INT) :: FieldmlUtil_ImportHandle
    
    !Local variables
    INTEGER(C_INT) :: importIndex, localHandle
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length
    
    FieldmlUtil_ImportHandle = FML_INVALID_HANDLE
    length = Fieldml_CopyObjectDeclaredName( fmlHandle, handle, name, BUFFER_SIZE )
    
    IF( Fieldml_IsObjectLocal( fmlHandle, handle ) /= 1 ) THEN
      IF( length > 0 ) THEN
        localHandle = Fieldml_GetObjectByName( fmlHandle, name(1:length) )
        IF( localHandle == FML_INVALID_HANDLE ) THEN
          importIndex = Fieldml_AddImportSource( fmlHandle, &
            & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
          FieldmlUtil_ImportHandle = Fieldml_AddImport( fmlHandle, importIndex, name(1:length), name(1:length) )
        ELSE IF( localHandle == handle ) THEN
          FieldmlUtil_ImportHandle = handle
        ENDIF
      ENDIF
    ENDIF

  END FUNCTION FieldmlUtil_ImportHandle
  
  !
  !================================================================================================================================
  !
  
  FUNCTION FieldmlUtil_GetTypeArgumentHandle( fmlInfo, typeHandle, doImport )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    LOGICAL, INTENT(IN) :: doImport
    INTEGER(C_INT), INTENT(IN) :: typeHandle
    
    INTEGER(C_INT) :: FieldmlUtil_GetTypeArgumentHandle

    !Local variables
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length
    INTEGER(C_INT) :: handle, err
    TYPE(VARYING_STRING) :: fullName
    
    length = Fieldml_CopyObjectName( fmlInfo%fmlHandle, typeHandle, name, BUFFER_SIZE )
    IF( length < 1 ) THEN
      length = Fieldml_CopyObjectDeclaredName( fmlInfo%fmlHandle, typeHandle, name, BUFFER_SIZE )
      FieldmlUtil_GetTypeArgumentHandle = FML_INVALID_HANDLE
      RETURN
    ENDIF

    IF( doImport ) THEN
      fullName = name(1:length)//".argument"//NUL
      err = FieldmlUtil_Import( fmlInfo%fmlHandle, fullName )
    ENDIF
    
    handle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, name(1:length)//".argument"//NUL )
    IF( handle == FML_INVALID_HANDLE ) THEN
      FieldmlUtil_GetTypeArgumentHandle = FML_INVALID_HANDLE
      RETURN
    ENDIF
    
    FieldmlUtil_GetTypeArgumentHandle = handle
    
  END FUNCTION FieldmlUtil_GetTypeArgumentHandle
  
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckLastInfoError( errorDescription, fmlInfo, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    TYPE(VARYING_STRING), INTENT(INOUT) :: errorString
    
    INTEGER(INTG) :: err, fmlErr
    
    fmlErr = Fieldml_GetLastError( fmlInfo%fmlHandle )

    IF( fmlErr == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    errorString = errorDescription // " (error number " // TRIM(NUMBER_TO_VSTRING(fmlErr,"*",err,errorString)) // ")"
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckLastInfoError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckLastError( errorDescription, fmlHandle, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    TYPE(VARYING_STRING), INTENT(INOUT) :: errorString
    
    INTEGER(INTG) :: err, fmlErr
    
    fmlErr = Fieldml_GetLastError( fmlHandle )

    IF( fmlErr == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    errorString = errorDescription // " (error number " // TRIM(NUMBER_TO_VSTRING(fmlErr,"*",err,errorString)) // ")"
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckLastError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckErrorNumber( errorDescription, error, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    INTEGER(INTG), INTENT(IN) :: error
    TYPE(VARYING_STRING), INTENT(INOUT) :: errorString
    
    INTEGER(INTG) :: err

    IF( error == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    errorString = errorDescription // " (error number " // TRIM(NUMBER_TO_VSTRING(error,"*",err,errorString)) // ")"
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckErrorNumber
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetCoordinatesType( fieldmlHandle, coordsType, dimensions, doImport, typeHandle, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: coordsType
    INTEGER(C_INT), INTENT(IN) :: dimensions
    LOGICAL, INTENT(IN) :: doImport
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: typeName
    INTEGER(C_INT) :: temp

    CALL ENTERS( "FieldmlUtil_GetCoordinatesType", err, errorString, *999 )
    
    IF( coordsType == COORDINATE_RECTANGULAR_CARTESIAN_TYPE ) THEN
      IF( dimensions == 1 ) THEN
        typeName = "coordinates.rc.1d"//NUL
      ELSE IF( dimensions == 2 ) THEN
        typeName = "coordinates.rc.2d"//NUL
      ELSE IF( dimensions == 3 ) THEN
        typeName = "coordinates.rc.3d"//NUL
      ELSE
        typeHandle = FML_INVALID_HANDLE
        err = FML_ERR_UNSUPPORTED
        CALL FieldmlUtil_CheckErrorNumber( "Cannot get RC coordinates type", err, errorString, *999 )
      ENDIF
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckErrorNumber( "Cannot get coordinates type", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, typeName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(typeName) )
    CALL FieldmlUtil_CheckError( "Cannot get coordinates type", fieldmlHandle, errorString, *999 )

    CALL EXITS( "FieldmlUtil_GetCoordinatesType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetCoordinatesType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetCoordinatesType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetCoordinatesType
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetGenericType( fieldmlHandle, dimensions, typeHandle, doImport, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: dimensions
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    LOGICAL, INTENT(IN) :: doImport
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: typeName
    INTEGER(C_INT) :: temp

    CALL ENTERS( "FieldmlUtil_GetGenericType", err, errorString, *999 )

    IF( dimensions == 1 ) THEN
      typeName = "real.1d"//NUL
    ELSE IF( dimensions == 2 ) THEN
      typeName = "real.2d"//NUL
    ELSE IF( dimensions == 3 ) THEN
      typeName = "real.3d"//NUL
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckErrorNumber( "Cannot get generic type", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, typeName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(typeName) )
    CALL FieldmlUtil_CheckError( "Cannot get generic type", fieldmlHandle, errorString, *999 )
    
    CALL EXITS( "FieldmlUtil_GetGenericType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetGenericType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetGenericType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetGenericType
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetXiType( fieldmlHandle, dimensions, doImport, typeHandle, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: dimensions
    LOGICAL, INTENT(IN) :: doImport
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string
    
    !Locals
    INTEGER(C_INT) :: temp
    TYPE(VARYING_STRING) :: typeName

    CALL ENTERS( "FieldmlUtil_GetXiType", err, errorString, *999 )
    
    IF( dimensions == 1 ) THEN
      typeName = "chart.1d"//NUL
    ELSE IF( dimensions == 2 ) THEN
      typeName = "chart.2d"//NUL
    ELSE IF( dimensions == 3 ) THEN
      typeName = "chart.3d"//NUL
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get xi type", fieldmlHandle, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, typeName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(typeName) )
    CALL FieldmlUtil_CheckError( "Cannot get xi type", fieldmlHandle, errorString, *999 )
    
    CALL EXITS( "FieldmlUtil_GetXiType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetXiType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetXiType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetXiType
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err, errorString )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:)
    TYPE(VARYING_STRING), INTENT(OUT) :: suffix
    INTEGER(INTG), INTENT(OUT) :: err
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
  
  SUBROUTINE FieldmlUtil_GetTPConnectivityEnsemble( fieldmlHandle, xiInterpolations, collapseInfo, doImport, typeHandle, &
    & err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: xiInterpolations(:)
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:)
    LOGICAL, INTENT(IN) :: doImport
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    !Locals
    INTEGER(C_INT) :: xiCount, firstInterpolation, i, importIndex, temp
    TYPE(VARYING_STRING) :: suffix, layoutName
    
    CALL ENTERS( "FieldmlUtil_GetTPConnectivityEnsemble", err, errorString, *999 )

    xiCount = SIZE( xiInterpolations )
  
    importIndex = Fieldml_AddImportSource( fieldmlHandle, &
      & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot access library", fieldmlHandle, errorString, *999 )

    firstInterpolation = xiInterpolations(1)
    DO i = 2, xiCount
      IF( xiInterpolations(i) /= firstInterpolation ) THEN
        !Do not yet support inhomogeneous TP bases
        err = FML_ERR_INVALID_OBJECT
        CALL FieldmlUtil_CheckErrorNumber( "Inhomogeneous tensor-product bases are not yet supported", err, errorString, *999 )
      ENDIF
    ENDDO

    CALL FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err, errorString )
      
    IF( firstInterpolation == BASIS_QUADRATIC_LAGRANGE_INTERPOLATION ) THEN
      IF( xiCount == 1 ) THEN
        layoutName = "localNodes.1d.line3"//NUL
      ELSE IF( xiCount == 2 ) THEN
        layoutName = "localNodes.2d.square3x3"//char(suffix)//NUL
      ELSE IF( xiCount == 3 ) THEN
        layoutName = "localNodes.3d.cube3x3x3"//char(suffix)//NUL
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
        CALL FieldmlUtil_CheckErrorNumber( "Quadratic interpolation not support for more than 3 dimensions", &
          & err, errorString, *999 )
      ENDIF
    ELSE IF( firstInterpolation == BASIS_LINEAR_LAGRANGE_INTERPOLATION ) THEN
      IF( xiCount == 1 ) THEN
        layoutName = "localNodes.1d.line2"//NUL
      ELSE IF( xiCount == 2 ) THEN
        layoutName = "localNodes.2d.square2x2"//char(suffix)//NUL
      ELSE IF( xiCount == 3 ) THEN
        layoutName = "localNodes.3d.cube2x2x2"//char(suffix)//NUL
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
        CALL FieldmlUtil_CheckErrorNumber( "Linear interpolation not support for more than 3 dimensions", &
          & err, errorString, *999 )
      ENDIF
    ELSE
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckErrorNumber( "Interpolation not yet supported", err, errorString, *999 )
    ENDIF

    IF( doImport ) THEN
      temp = FieldmlUtil_Import( fieldmlHandle, layoutName )
    ENDIF
    typeHandle = Fieldml_GetObjectByName( fieldmlHandle, char(layoutName) )
    CALL FieldmlUtil_CheckError( "Cannot get local nodes type", fieldmlHandle, errorString, *999 )
    
    CALL EXITS( "FieldmlUtil_GetTPConnectivityEnsemble" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetTPConnectivityEnsemble", err, errorString )
    CALL EXITS( "FieldmlUtil_GetTPConnectivityEnsemble" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetTPConnectivityEnsemble

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlUtil_GetConnectivityEnsemble( fieldmlHandle, basis, typeHandle, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fieldmlHandle
    TYPE(BASIS_TYPE), POINTER :: basis
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    !Locals
    INTEGER(C_INT) :: basisType, xiCount
    INTEGER(C_INT), ALLOCATABLE :: xiInterpolations(:), collapseInfo(:)
    
    CALL ENTERS( "FieldmlUtil_GetConnectivityEnsemble", err, errorString, *999 )

    typeHandle = FML_INVALID_HANDLE

    CALL BASIS_TYPE_GET( basis, basisType, err, errorString, *999 )
    CALL BASIS_NUMBER_OF_XI_GET( basis, xiCount, err, errorString, *999 )
    
    IF( basisType == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      ALLOCATE( xiInterpolations( xiCount ) )
      ALLOCATE( collapseInfo( xiCount ) )
      CALL BASIS_INTERPOLATION_XI_GET( basis, xiInterpolations, err, errorString, *999 )
      CALL BASIS_COLLAPSED_XI_GET( basis, collapseInfo, err, errorString, *999 )
      
      CALL FieldmlUtil_GetTPConnectivityEnsemble( fieldmlHandle, xiInterpolations, collapseInfo, .TRUE., typeHandle, &
        & err, errorString, *999 )
      
      DEALLOCATE( xiInterpolations )
      DEALLOCATE( collapseInfo )
    ELSE
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckErrorNumber( "Only tensor product bases are currently supported", err, errorString, *999 )
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
  
  SUBROUTINE FieldmlUtil_GetValueType( fmlHandle, region, field, typeHandle, doImport, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(FIELD_TYPE), POINTER :: field
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    LOGICAL, INTENT(IN) :: doImport
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string

    !Locals
    INTEGER(INTG) :: fieldType, subType, count
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem

    CALL ENTERS( "FieldmlUtil_GetValueType", err, errorString, *999 )
    
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
  SUBROUTINE FieldmlUtil_FinaliseInfo( fieldmlInfo )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo

    !Locals
    INTEGER(INTG) :: err

    err = Fieldml_Destroy( fieldmlInfo%fmlHandle )
    
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
    
    IF( ALLOCATED( fieldmlInfo%componentHandles ) ) THEN
      DEALLOCATE( fieldmlInfo%componentHandles )
    ENDIF
    IF( ALLOCATED( fieldmlInfo%basisHandles ) ) THEN
      DEALLOCATE( fieldmlInfo%basisHandles )
    ENDIF
    IF( ALLOCATED( fieldmlInfo%basisConnectivityHandles ) ) THEN
      DEALLOCATE( fieldmlInfo%basisConnectivityHandles )
    ENDIF
    IF( ALLOCATED( fieldmlInfo%basisLayoutHandles ) ) THEN
      DEALLOCATE( fieldmlInfo%basisLayoutHandles )
    ENDIF
    
  END SUBROUTINE FieldmlUtil_FinaliseInfo

  !
  !================================================================================================================================
  !

END MODULE FIELDML_UTIL_ROUTINES
