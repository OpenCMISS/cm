!> \file
!> $Id$
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
  USE ISO_VARYING_STRING
  USE STRINGS
  USE OPENCMISS
  USE BASE_ROUTINES

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: BUFFER_SIZE = 1024

  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  TYPE(VARYING_STRING) :: errorString

  !Interfaces
  TYPE FieldmlInfoType
    TYPE(C_PTR) :: fmlHandle
    INTEGER(C_INT) :: nodesHandle
    INTEGER(C_INT) :: nodesVariableHandle
    INTEGER(C_INT) :: meshHandle
    INTEGER(C_INT) :: elementsHandle
    INTEGER(C_INT) :: elementsVariableHandle
    INTEGER(C_INT) :: xiHandle
    INTEGER(C_INT) :: xiVariableHandle
    INTEGER(C_INT) :: nodeDofsHandle
!    INTEGER(C_INT) :: elementDofsHandle
!    INTEGER(C_INT) :: constantDofsHandle
    INTEGER(C_INT), ALLOCATABLE :: componentHandles(:)
    INTEGER(C_INT), ALLOCATABLE :: basisHandles(:)
    INTEGER(C_INT), ALLOCATABLE :: basisConnectivityHandles(:)
  END TYPE FieldmlInfoType
  
  INTERFACE FieldmlUtil_CheckError
    MODULE PROCEDURE FieldmlUtil_CheckLastError
    MODULE PROCEDURE FieldmlUtil_CheckLastInfoError
    MODULE PROCEDURE FieldmlUtil_CheckThisError
  END INTERFACE

  PUBLIC :: FieldmlInfoType

  PUBLIC :: FieldmlUtil_GetConnectivityEnsemble, FieldmlUtil_GetGenericType, &
    & FieldmlUtil_GetXiEnsemble, FieldmlUtil_GetXiType, FieldmlUtil_GetValueType, FieldmlUtil_FinalizeInfo, &
    & FieldmlUtil_GetCollapseSuffix, FieldmlUtil_CheckError, FieldmlUtil_GetTypeVariableHandle

CONTAINS

  !
  !================================================================================================================================
  !
  
  FUNCTION FieldmlUtil_GetTypeVariableHandle( fmlInfo, typeHandle )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: typeHandle
    
    INTEGER(C_INT) :: FieldmlUtil_GetTypeVariableHandle

    !Local variables
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length, err
    INTEGER(C_INT) :: handle
    
    length = Fieldml_CopyObjectName( fmlInfo%fmlHandle, typeHandle, name, BUFFER_SIZE )
    IF( length < 1 ) THEN
      FieldmlUtil_GetTypeVariableHandle = FML_INVALID_HANDLE
      RETURN
    ENDIF
    
    handle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, name(1:length)//".variable"//NUL )
    IF( handle == FML_INVALID_HANDLE ) THEN
      FieldmlUtil_GetTypeVariableHandle = FML_INVALID_HANDLE
      RETURN
    ENDIF
    
    FieldmlUtil_GetTypeVariableHandle = handle
    
  END FUNCTION FieldmlUtil_GetTypeVariableHandle
  
  
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
    
    errorString = errorDescription // "(error number " // TRIM(NUMBER_TO_VSTRING(fmlErr,"*",err,errorString)) // ")"
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckLastInfoError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckLastError( errorDescription, fmlHandle, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    TYPE(VARYING_STRING), INTENT(INOUT) :: errorString
    
    INTEGER(INTG) :: err, fmlErr
    
    fmlErr = Fieldml_GetLastError( fmlHandle )

    IF( fmlErr == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    errorString = errorDescription // "(error number " // TRIM(NUMBER_TO_VSTRING(fmlErr,"*",err,errorString)) // ")"
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckLastError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CheckThisError( errorDescription, error, errorString, * )
    CHARACTER(LEN=*), INTENT(IN) :: errorDescription
    INTEGER(INTG), INTENT(IN) :: error
    TYPE(VARYING_STRING), INTENT(INOUT) :: errorString
    
    INTEGER(INTG) :: err

    IF( error == FML_ERR_NO_ERROR ) THEN
      RETURN
    ENDIF
    
    errorString = errorDescription // "(error number " // TRIM(NUMBER_TO_VSTRING(error,"*",err,errorString)) // ")"
    RETURN 1
    
  END SUBROUTINE FieldmlUtil_CheckThisError
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetCoordinatesType( fieldmlHandle, coordsType, dimensions, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: coordsType
    INTEGER(C_INT), INTENT(IN) :: dimensions
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlUtil_GetCoordinatesType", err, errorString, *999 )
    
    IF( coordsType == CMISSCoordinateRectangularCartesianType ) THEN
      IF( dimensions == 1 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.coordinates.rc.1d"//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get RC1 coordinates type", fieldmlHandle, errorString, *999 )
      ELSE IF( dimensions == 2 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.coordinates.rc.2d"//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get RC2 coordinates type", fieldmlHandle, errorString, *999 )
      ELSE IF( dimensions == 3 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.coordinates.rc.3d"//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get RC3 coordinates type", fieldmlHandle, errorString, *999 )
      ELSE
        typeHandle = FML_INVALID_HANDLE
        err = FML_ERR_UNSUPPORTED
        CALL FieldmlUtil_CheckError( "Cannot get RC coordinates type", err, errorString, *999 )
      ENDIF
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get coordinates type", err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlUtil_GetCoordinatesType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetCoordinatesType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetCoordinatesType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetCoordinatesType
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetGenericType( fieldmlHandle, dimensions, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: dimensions
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlUtil_GetGenericType", err, errorString, *999 )
    
    IF( dimensions == 1 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.real.1d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 1D generic type", fieldmlHandle, errorString, *999 )
    ELSE IF( dimensions == 2 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.real.2d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 2D generic type", fieldmlHandle, errorString, *999 )
    ELSE IF( dimensions == 3 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.real.3d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 3D generic type", fieldmlHandle, errorString, *999 )
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get generic type", err, errorString, *999 )
    ENDIF
    
    CALL EXITS( "FieldmlUtil_GetGenericType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetGenericType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetGenericType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetGenericType
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetXiEnsemble( fieldmlHandle, dimensions, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: dimensions
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlUtil_GetXiEnsemble", err, errorString, *999 )
    
    IF( dimensions == 1 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.ensemble.xi.1d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 1D xi ensemble", fieldmlHandle, errorString, *999 )
    ELSE IF( dimensions == 2 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.ensemble.xi.2d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 2D xi ensemble", fieldmlHandle, errorString, *999 )
    ELSE IF( dimensions == 3 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.ensemble.xi.3d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 3D xi ensemble", fieldmlHandle, errorString, *999 )
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get xi ensemble", err, errorString, *999 )
    ENDIF
    
    CALL EXITS( "FieldmlUtil_GetXiEnsemble" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetXiEnsemble", err, errorString )
    CALL EXITS( "FieldmlUtil_GetXiEnsemble" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetXiEnsemble
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetXiType( fieldmlHandle, dimensions, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: dimensions
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlUtil_GetXiType", err, errorString, *999 )
    
    IF( dimensions == 1 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.xi.1d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 1D xi type", fieldmlHandle, errorString, *999 )
    ELSE IF( dimensions == 2 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.xi.2d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 2D xi type", fieldmlHandle, errorString, *999 )
    ELSE IF( dimensions == 3 ) THEN
      typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.xi.3d"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot get 3D xi type", fieldmlHandle, errorString, *999 )
    ELSE
      typeHandle = FML_INVALID_HANDLE
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get xi type", fieldmlHandle, errorString, *999 )
    ENDIF
    
    CALL EXITS( "FieldmlUtil_GetXiType" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetXiType", err, errorString )
    CALL EXITS( "FieldmlUtil_GetXiType" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetXiType
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:)
    TYPE(VARYING_STRING), INTENT(OUT) :: suffix
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Locals
    INTEGER(INTG) :: i
    
    suffix = ""
    DO i = 1, SIZE( collapseInfo )
      IF( collapseInfo( i ) == CMISSBasisXiCollapsed ) THEN
        suffix = suffix // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",err,errorString))//"C"
      ELSEIF( collapseInfo( i ) == CMISSBasisCollapsedAtXi0 ) THEN
        suffix = suffix // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",err,errorString))//"0"
      ELSEIF( collapseInfo( i ) == CMISSBasisCollapsedAtXi1 ) THEN
        suffix = suffix // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",err,errorString))//"1"
      ENDIF
    ENDDO
  
  END SUBROUTINE
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetTPConnectivityEnsemble( fieldmlHandle, xiInterpolations, collapseInfo, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: xiInterpolations(:)
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:)
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: xiCount, firstInterpolation, i
    TYPE(VARYING_STRING) :: suffix
    
    CALL ENTERS( "FieldmlUtil_GetTPConnectivityEnsemble", err, errorString, *999 )

    xiCount = SIZE( xiInterpolations )
  
    firstInterpolation = xiInterpolations(1)
    DO i = 2, xiCount
      IF( xiInterpolations(i) /= firstInterpolation ) THEN
        !Do not yet support inhomogeneous TP bases
        err = FML_ERR_INVALID_OBJECT
        CALL FieldmlUtil_CheckError( "Inhomogeneous tensor-product bases are not yet supported", err, errorString, *999 )
      ENDIF
    ENDDO

    CALL FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err )
      
    IF( firstInterpolation == CMISSBasisQuadraticLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.local_nodes.line.3"//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get quadratic local nodes type", fieldmlHandle, errorString, *999 )
      ELSE IF( xiCount == 2 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.local_nodes.square.3x3"//char(suffix)//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get biquadratic local nodes type", fieldmlHandle, errorString, *999 )
      ELSE IF( xiCount == 3 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.local_nodes.cube.3x3x3"//char(suffix)//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get triquadratic local nodes type", fieldmlHandle, errorString, *999 )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
        CALL FieldmlUtil_CheckError( "Quadratic interpolation not support for more than 3 dimensions", &
          & err, errorString, *999 )
      ENDIF
    ELSE IF( firstInterpolation == CMISSBasisLinearLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.local_nodes.line.2"//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get linear local nodes type", fieldmlHandle, errorString, *999 )
      ELSE IF( xiCount == 2 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.local_nodes.square.2x2"//char(suffix)//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get bilinear local nodes type", fieldmlHandle, errorString, *999 )
      ELSE IF( xiCount == 3 ) THEN
        typeHandle = Fieldml_GetObjectByName( fieldmlHandle, "library.local_nodes.cube.2x2x2"//char(suffix)//NUL )
        CALL FieldmlUtil_CheckError( "Cannot get trilinear local nodes type", fieldmlHandle, errorString, *999 )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
        CALL FieldmlUtil_CheckError( "Linear interpolation not support for more than 3 dimensions", &
          & err, errorString, *999 )
      ENDIF
    ELSE
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckError( "Interpolation not yet supported", err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlUtil_GetTPConnectivityEnsemble" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_GetTPConnectivityEnsemble", err, errorString )
    CALL EXITS( "FieldmlUtil_GetTPConnectivityEnsemble" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_GetTPConnectivityEnsemble

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlUtil_GetConnectivityEnsemble( fieldmlHandle, basisNumber, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fieldmlHandle
    INTEGER(C_INT), INTENT(IN) :: basisNumber
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: basisType, xiCount
    INTEGER(C_INT), ALLOCATABLE :: xiInterpolations(:), collapseInfo(:)
    
    CALL ENTERS( "FieldmlUtil_GetConnectivityEnsemble", err, errorString, *999 )

    typeHandle = FML_INVALID_HANDLE
    
    CALL CMISSBasisTypeGet( basisNumber, basisType, err )
    CALL CMISSBasisNumberOfXiGet( basisNumber, xiCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get basis info for connectivity", err, errorString, *999 )
    
    IF( basisType == CMISSBasisLagrangeHermiteTPType ) THEN
      ALLOCATE( xiInterpolations( xiCount ) )
      ALLOCATE( collapseInfo( xiCount ) )
      CALL CMISSBasisInterpolationXiGet( basisNumber, xiInterpolations, err )
      CALL CMISSBasisCollapsedXiGet( basisNumber, collapseInfo, err )
      CALL FieldmlUtil_CheckError( "Cannot get basis interpolation info for connectivity", err, errorString, *999 )
      
      CALL FieldmlUtil_GetTPConnectivityEnsemble( fieldmlHandle, xiInterpolations, collapseInfo, typeHandle, err, *999 )
      
      DEALLOCATE( xiInterpolations )
      DEALLOCATE( collapseInfo )
    ELSE
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckError( "Only tensor product bases are currently supported", err, errorString, *999 )
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
  
  SUBROUTINE FieldmlUtil_GetValueType( fmlHandle, region, field, typeHandle, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    TYPE(CMISSRegionType), INTENT(IN) :: region
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(C_INT), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(INTG) :: fieldType, subType, count
    TYPE(CMISSCoordinateSystemType) coordinateSystem

    CALL ENTERS( "FieldmlUtil_GetValueType", err, errorString, *999 )
    
    CALL CMISSFieldTypeGet( field, fieldType, err )
    CALL CMISSFieldNumberOfComponentsGet( field, CMISSFieldUVariableType, count, err )
    CALL FieldmlUtil_CheckError( "Cannot get field info for value type", err, errorString, *999 )

    SELECT CASE( fieldType )
    CASE( CMISSFieldGeometricType )
      CALL CMISSCoordinateSystemTypeInitialise( coordinateSystem, err )
      CALL CMISSRegionCoordinateSystemGet( region, coordinateSystem, err )
      CALL CMISSCoordinateSystemTypeGet( coordinateSystem, subType, err )
      CALL FieldmlUtil_CheckError( "Cannot get coordinate system info for geometric field", err, errorString, *999 )
      CALL FieldmlUtil_GetCoordinatesType( fmlHandle, subType, count, typeHandle, err, *999 )
    
    !CASE( CMISSFieldFibreType )

    !CASE( CMISSFieldGeneralType )

    !CASE( CMISSFieldMaterialType )

    CASE DEFAULT
      CALL FieldmlUtil_GetGenericType( fmlHandle, count, typeHandle, err, *999 )
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
  SUBROUTINE FieldmlUtil_FinalizeInfo( fieldmlInfo )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo

    !Locals
    INTEGER(INTG) :: err

    err = Fieldml_Destroy( fieldmlInfo%fmlHandle )
    
    fieldmlInfo%fmlHandle = C_NULL_PTR
    fieldmlInfo%nodesHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodesVariableHandle = FML_INVALID_HANDLE
    fieldmlInfo%meshHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsVariableHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiVariableHandle = FML_INVALID_HANDLE
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
    
  END SUBROUTINE FieldmlUtil_FinalizeInfo

  !
  !================================================================================================================================
  !

END MODULE FIELDML_UTIL_ROUTINES
