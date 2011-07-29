!> \file
!> \author Caton Little
!> \brief This module handles reading in FieldML files.
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

!> Input routines for FieldML

MODULE FIELDML_INPUT_ROUTINES

  USE FIELDML_API
  USE FIELDML_UTIL_ROUTINES
  USE FIELDML_TYPES
  USE FIELD_ROUTINES
  USE BASIS_ROUTINES
  USE COORDINATE_ROUTINES
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE REGION_ROUTINES
  USE UTIL_ARRAY
  USE CMISS

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: BUFFER_SIZE = 1024

  INTEGER(INTG), PARAMETER :: FML_ERR_UNKNOWN_BASIS = 10001
  INTEGER(INTG), PARAMETER :: FML_ERR_INVALID_BASIS = 10002
  INTEGER(INTG), PARAMETER :: FML_ERR_UNKNOWN_MESH_XI = 10003
  INTEGER(INTG), PARAMETER :: FML_ERR_UNKNOWN_COORDINATE_TYPE = 10004
  INTEGER(INTG), PARAMETER :: FML_ERR_INVALID_PARAMETER = 10007
  INTEGER(INTG), PARAMETER :: FML_ERR_INVALID_MESH = 10008
  INTEGER(INTG), PARAMETER :: FML_ERR_INVALID_CONNECTIVITY = 10009
  INTEGER(INTG), PARAMETER :: FML_ERR_INVALID_READ = 10010

  CHARACTER(KIND=C_CHAR), PARAMETER :: NUL = C_NULL_CHAR

  !Shim type because Fortran can't handle arrays of arrays
  TYPE ArrayShimType
    INTEGER(C_INT), ALLOCATABLE :: array(:)
  END TYPE ArrayShimType

  !Interfaces

  INTERFACE

  END INTERFACE

  PUBLIC :: FieldmlInput_InitialiseFromFile, FieldmlInput_MeshCreateStart, &
    & FieldmlInput_CoordinateSystemCreateStart, FieldmlInput_BasisCreateStart, FieldmlInput_CreateMeshComponent, &
    & FieldmlInput_FieldCreateStart, FieldmlInput_FieldNodalParametersUpdate, FieldmlInput_NodesCreateStart

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_GetBasisConnectivityInfo( fmlInfo, basisHandle, paramArgHandle, connectivityHandle, layoutHandle, err, &
    & errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: basisHandle
    INTEGER(C_INT), INTENT(IN) :: paramArgHandle
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Local variables
    INTEGER(C_INT) :: count, bindNumber, paramsHandle, argHandle, layoutIndexHandle
    
    CALL ENTERS( "FieldmlInput_GetBasisConnectivityInfo", err, errorString, *999 )

    count = Fieldml_GetBindCount( fmlInfo%fmlHandle, basisHandle )
    IF( count /= 2 ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Library basis evaluators must have exactly two binds", fmlInfo, errorString, *999 )
    END IF
    
    paramsHandle = FML_INVALID_HANDLE
    DO bindNumber = 1, count
      argHandle = Fieldml_GetBindArgument( fmlInfo%fmlHandle, basisHandle, bindNumber )
      CALL FieldmlUtil_CheckError( "Cannot get bind for interpolator", fmlInfo%fmlHandle, errorString, *999 )
      IF( argHandle == paramArgHandle ) THEN
        paramsHandle = Fieldml_GetBindEvaluator( fmlInfo%fmlHandle, basisHandle, bindNumber )
      ENDIF
    ENDDO

    IF( paramsHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckErrorNumber( "Library interpolators must have a correct parameter bind", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( fmlInfo%fmlHandle, paramsHandle ) /= FHT_AGGREGATE_EVALUATOR ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckErrorNumber( "Parameter evaluator for interpolator must be an aggregate", err, errorString, *999 )
    ENDIF
    
    count = Fieldml_GetBindCount( fmlInfo%fmlHandle, paramsHandle )
    IF( count /= 1 ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckErrorNumber( "Nodal parameter evaluator must only have one bind", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetBindArgument( fmlInfo%fmlHandle, paramsHandle, 1 ) /= fmlInfo%nodesArgumentHandle ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckErrorNumber( "Nodal parameter evaluator must bind the nodes argument", err, errorString, *999 )
    ENDIF
    
    connectivityHandle = Fieldml_GetBindEvaluator( fmlInfo%fmlHandle, paramsHandle, 1 )
    CALL FieldmlUtil_CheckError( "Cannot get connectivity source for nodal parameters", fmlInfo%fmlHandle, &
      & errorString, *999 )
      
    layoutIndexHandle = Fieldml_GetIndexEvaluator( fmlInfo%fmlHandle, paramsHandle, 1 )
    layoutHandle = Fieldml_GetValueType( fmlInfo%fmlHandle, layoutIndexHandle )
    CALL FieldmlUtil_CheckError( "Cannot get connectivity source for nodal parameters", fmlInfo%fmlHandle, &
      & errorString, *999 )

    CALL EXITS( "FieldmlInput_GetBasisConnectivityInfo" )
    RETURN
999 CALL ERRORS( "FieldmlInput_GetBasisConnectivityInfo", err, errorString )
    CALL EXITS( "FieldmlInput_GetBasisConnectivityInfo" )
    RETURN 1

  END SUBROUTINE FieldmlInput_GetBasisConnectivityInfo

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_GetBasisCollapse( name, collapse )
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: collapse(:)

    collapse = BASIS_NOT_COLLAPSED
    
    IF( SIZE( collapse ) > 0 ) THEN
      IF( INDEX( name, "_xi1C" ) /= 0 ) THEN
        collapse(1) = BASIS_XI_COLLAPSED
      ELSE IF( INDEX( name, "_xi10" ) /= 0 ) THEN
        collapse(1) = BASIS_COLLAPSED_AT_XI0
      ELSE IF( INDEX( name, "_xi11" ) /= 0 ) THEN
        collapse(1) = BASIS_COLLAPSED_AT_XI1
      ENDIF
    ENDIF
  
    IF( SIZE( collapse ) > 1 ) THEN
      IF( INDEX( name, "_xi2C" ) /= 0 ) THEN
        collapse(2) = BASIS_XI_COLLAPSED
      ELSE IF( INDEX( name, "_xi20" ) /= 0 ) THEN
        collapse(2) = BASIS_COLLAPSED_AT_XI0
      ELSE IF( INDEX( name, "_xi21" ) /= 0 ) THEN
        collapse(2) = BASIS_COLLAPSED_AT_XI1
      ENDIF
    ENDIF
  
    IF( SIZE( collapse ) > 2 ) THEN
      IF( INDEX( name, "_xi3C" ) /= 0 ) THEN
        collapse(3) = BASIS_XI_COLLAPSED
      ELSE IF( INDEX( name, "_xi30" ) /= 0 ) THEN
        collapse(3) = BASIS_COLLAPSED_AT_XI0
      ELSE IF( INDEX( name, "_xi31" ) /= 0 ) THEN
        collapse(3) = BASIS_COLLAPSED_AT_XI1
      ENDIF
    ENDIF
  
  END SUBROUTINE

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_GetBasisInfo( fmlInfo, objectHandle, connectivityHandle, layoutHandle, basisType, basisInterpolations, &
    & collapse, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: objectHandle
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle
    INTEGER(INTG), INTENT(OUT) :: basisType
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: basisInterpolations(:)
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: collapse(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(C_INT) :: length, libraryBasisHandle, paramArgHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    CALL ENTERS( "FieldmlInput_GetBasisInfo", err, errorString, *999 )

    IF( .NOT. FieldmlInput_IsKnownBasis( fmlInfo, objectHandle, err ) ) THEN
      CALL FieldmlUtil_CheckErrorNumber( "Basis specified in FieldML file is not yet supported", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( fmlInfo%fmlHandle, objectHandle ) /= FHT_REFERENCE_EVALUATOR ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckErrorNumber( "Basis evaluator must be a continuous reference", err, errorString, *999 )
    ENDIF
    
    libraryBasisHandle = Fieldml_GetReferenceSourceEvaluator( fmlInfo%fmlHandle, objectHandle )
    CALL FieldmlUtil_CheckError( "Basis specified in FieldML is not a reference evaluator", fmlInfo, errorString, *999 )
    length = Fieldml_CopyObjectDeclaredName( fmlInfo%fmlHandle, libraryBasisHandle, name, BUFFER_SIZE )
    CALL FieldmlUtil_CheckError( "Cannot get name of basis evaluator", fmlInfo, errorString, *999 )

    IF( INDEX( name, 'interpolator.3d.unit.triquadraticLagrange') == 1 ) THEN
      paramArgHandle = Fieldml_GetObjectByDeclaredName( fmlInfo%fmlHandle, &
        & "parameters.3d.unit.triquadraticLagrange.argument"//NUL )
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      CALL REALLOCATE_INT( collapse, 3, "", err, errorString, *999 )
      basisInterpolations = BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
      basisType = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( name, 'interpolator.3d.unit.trilinearLagrange') == 1 ) THEN
      paramArgHandle = Fieldml_GetObjectByDeclaredName( fmlInfo%fmlHandle, "parameters.3d.unit.trilinearLagrange.argument"//NUL )
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      CALL REALLOCATE_INT( collapse, 3, "", err, errorString, *999 )
      basisInterpolations = BASIS_LINEAR_LAGRANGE_INTERPOLATION
      basisType = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE
      err = FML_ERR_UNKNOWN_BASIS
      CALL FieldmlUtil_CheckErrorNumber( "Basis cannot yet be interpreted", err, errorString, *999 )
    ENDIF
    
    IF( basisType == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      CALL FieldmlInput_GetBasisCollapse( name(1:length), collapse )
    ENDIF
    
    CALL FieldmlInput_GetBasisConnectivityInfo( fmlInfo, objectHandle, paramArgHandle, connectivityHandle, layoutHandle, &
      & err, errorString, *999 )

    CALL EXITS( "FieldmlInput_GetBasisInfo" )
    RETURN
999 CALL ERRORS( "FieldmlInput_GetBasisInfo", err, errorString )
    CALL EXITS( "FieldmlInput_GetBasisInfo" )
    RETURN 1

  END SUBROUTINE FieldmlInput_GetBasisInfo

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_IsKnownBasis( fmlInfo, objectHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: objectHandle
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Function
    LOGICAL :: FieldmlInput_IsKnownBasis

    !Locals
    INTEGER(C_INT) :: length, libraryBasisHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    FieldmlInput_IsKnownBasis = .FALSE.

    IF( Fieldml_GetObjectType( fmlInfo%fmlHandle, objectHandle ) /= FHT_REFERENCE_EVALUATOR ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF

    libraryBasisHandle = Fieldml_GetReferenceSourceEvaluator( fmlInfo%fmlHandle, objectHandle )
    length = Fieldml_CopyObjectDeclaredName( fmlInfo%fmlHandle, libraryBasisHandle, name, BUFFER_SIZE )

    IF( ( INDEX( name, 'interpolator.3d.unit.triquadraticLagrange') /= 1 ) .AND. &
      & ( INDEX( name, 'interpolator.3d.unit.trilinearLagrange') /= 1 ) ) THEN
      err = FML_ERR_UNKNOWN_BASIS
      RETURN
    ENDIF
    
    FieldmlInput_IsKnownBasis = .TRUE.
    RETURN

    err = FML_ERR_INVALID_BASIS
    
  END FUNCTION FieldmlInput_IsKnownBasis
  
  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_IsElementEvaluatorCompatible( fmlHandle, object, err )
    !Arguments
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: object
    INTEGER(INTG), INTENT(OUT) :: err

    LOGICAL :: FieldmlInput_IsElementEvaluatorCompatible

    INTEGER(C_INT) :: type, length, evaluatorHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name

    type = Fieldml_GetObjectType( fmlHandle, object )
    IF( type /= FHT_REFERENCE_EVALUATOR ) THEN
      FieldmlInput_IsElementEvaluatorCompatible = .FALSE.
      RETURN
    ENDIF

    evaluatorHandle = Fieldml_GetReferenceSourceEvaluator( fmlHandle, object )
    length = Fieldml_CopyObjectDeclaredName( fmlHandle, evaluatorHandle, name, BUFFER_SIZE )

    IF( INDEX( name, 'interpolator.3d.unit.trilinearLagrange' ) == 1 ) THEN
      FieldmlInput_IsElementEvaluatorCompatible = .TRUE.
    ELSE IF( INDEX( name, 'interpolator.3d.unit.triquadraticLagrange' ) == 1 ) THEN
      FieldmlInput_IsElementEvaluatorCompatible = .TRUE.
    ELSE
      FieldmlInput_IsElementEvaluatorCompatible = .FALSE.
    ENDIF

    err = FML_ERR_NO_ERROR

  END FUNCTION FieldmlInput_IsElementEvaluatorCompatible

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_IsTemplateCompatible( fmlHandle, object, elementType, err )
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: object
    INTEGER(C_INT), INTENT(IN) :: elementType
    INTEGER(INTG), INTENT(OUT) :: err

    LOGICAL :: FieldmlInput_IsTemplateCompatible

    INTEGER(C_INT) :: objectType, count, i, evaluator, type, firstEvaluator, evaluatorHandle, defaultEvaluator

    objectType = Fieldml_GetObjectType( fmlHandle, object )
    IF( objectType /= FHT_PIECEWISE_EVALUATOR ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    evaluatorHandle = Fieldml_GetIndexEvaluator( fmlHandle, object, 1 )
    type = Fieldml_GetValueType( fmlHandle, evaluatorHandle )
    IF( type /= elementType ) THEN
      FieldmlInput_IsTemplateCompatible = .TRUE.
      RETURN
    ENDIF

    count = Fieldml_GetEvaluatorCount( fmlHandle, object )
    defaultEvaluator = Fieldml_GetDefaultEvaluator( fmlHandle, object )
    
    IF( ( defaultEvaluator /= FML_INVALID_HANDLE ) .AND. .NOT. &
      & FieldmlInput_IsElementEvaluatorCompatible( fmlHandle, defaultEvaluator, err ) ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    IF( count == 0 ) THEN
      IF( defaultEvaluator == FML_INVALID_HANDLE ) THEN
          FieldmlInput_IsTemplateCompatible = .FALSE.
      ELSE
          FieldmlInput_IsTemplateCompatible = .TRUE.
      ENDIF
      RETURN
    ENDIF

    firstEvaluator = Fieldml_GetEvaluator( fmlHandle, object, 1 )
    IF( .NOT. FieldmlInput_IsElementEvaluatorCompatible( fmlHandle, firstEvaluator, err ) ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    !At the moment, the code does not support different evaluators per element.

    DO i = 2, count
      evaluator = Fieldml_GetEvaluator( fmlHandle, object, i )
      IF( evaluator /= firstEvaluator ) THEN
        FieldmlInput_IsTemplateCompatible = .FALSE.
        RETURN
      ENDIF
    ENDDO

    FieldmlInput_IsTemplateCompatible = .TRUE.

  END FUNCTION FieldmlInput_IsTemplateCompatible

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_IsFieldCompatible( fmlHandle, object, elementType, err )
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: object
    INTEGER(C_INT), INTENT(IN) :: elementType
    INTEGER(INTG), INTENT(OUT) :: err

    LOGICAL :: FieldmlInput_IsFieldCompatible

    INTEGER(C_INT) :: type, count, i, evaluator, defaultEvaluator

    type = Fieldml_GetObjectType( fmlHandle, object )

    IF( type /= FHT_AGGREGATE_EVALUATOR ) THEN
      FieldmlInput_IsFieldCompatible = .FALSE.
      RETURN
    ENDIF

    count = Fieldml_GetEvaluatorCount( fmlHandle, object )
    defaultEvaluator = Fieldml_GetDefaultEvaluator( fmlHandle, object )

    IF( ( defaultEvaluator /= FML_INVALID_HANDLE ) .AND. .NOT. &
      & FieldmlInput_IsElementEvaluatorCompatible( fmlHandle, defaultEvaluator, err ) ) THEN
      FieldmlInput_IsFieldCompatible = .FALSE.
      RETURN
    ENDIF

    IF( count == 0 ) THEN
      IF( defaultEvaluator == FML_INVALID_HANDLE ) THEN
          FieldmlInput_IsFieldCompatible = .FALSE.
      ELSE
          FieldmlInput_IsFieldCompatible = .TRUE.
      ENDIF
      RETURN
    ENDIF

    FieldmlInput_IsFieldCompatible = .TRUE.
    DO i = 1, count
      evaluator = Fieldml_GetEvaluator( fmlHandle, object, i )
      IF( .NOT. FieldmlInput_IsTemplateCompatible( fmlHandle, evaluator, elementType, err ) ) THEN
        FieldmlInput_IsFieldCompatible = .FALSE.
        RETURN
      ENDIF
    ENDDO

  END FUNCTION FieldmlInput_IsFieldCompatible

  !
  !================================================================================================================================
  !

  SUBROUTINE Fieldml_GetFieldHandles( fmlHandle, fieldHandles, meshHandle, err, errorString, * )
    INTEGER(C_INT), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), ALLOCATABLE :: fieldHandles(:)
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    INTEGER(C_INT) :: count, i, object, fieldCount, elementType

    CALL ENTERS( "Fieldml_GetFieldHandles", err, errorString, *999 )

    elementType = Fieldml_GetMeshElementsType( fmlHandle, meshHandle )
    CALL FieldmlUtil_CheckError( "Cannot get mesh element type", fmlHandle, errorString, *999 )

    fieldCount = 0
    count = Fieldml_GetObjectCount( fmlHandle, FHT_AGGREGATE_EVALUATOR )
    CALL FieldmlUtil_CheckError( "Cannot find any aggregate evaluators", fmlHandle, errorString, *999 )
    DO i = 1, count
      object = Fieldml_GetObject( fmlHandle, FHT_AGGREGATE_EVALUATOR, i )
      CALL FieldmlUtil_CheckError( "Cannot get aggregate evaluator", fmlHandle, errorString, *999 )

      IF( .NOT. FieldmlInput_IsFieldCompatible( fmlHandle, object, elementType, err ) ) THEN
        CYCLE
      ENDIF

      CALL GROW_ARRAY( fieldHandles, 1, "", err, errorString, *999 )
      fieldCount = fieldCount + 1
      fieldHandles( fieldCount ) = object
    ENDDO

    CALL EXITS( "Fieldml_GetFieldHandles" )
    RETURN
999 CALL ERRORS( "Fieldml_GetFieldHandles", err, errorString )
    CALL EXITS( "Fieldml_GetFieldHandles" )
    RETURN 1

  END SUBROUTINE Fieldml_GetFieldHandles


  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_CoordinateSystemCreateStart( fmlInfo, evaluatorName, coordinateSystem, userNumber, &
    & err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER, INTENT(IN) :: coordinateSystem
    INTEGER(INTG), INTENT(IN) :: userNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(C_INT) :: evaluatorHandle
    INTEGER(C_INT) :: typeHandle, length
    CHARACTER(LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: coordinateType
    INTEGER(INTG) :: coordinateCount

    CALL ENTERS( "FieldmlInput_CoordinateSystemCreateStart", err, errorString, *999 )

    coordinateType = 0 !There doesn't seem to be a COORDINATE_UNKNOWN_TYPE

    evaluatorHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, evaluatorName//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get coordinate evaluator for geometric field", fmlInfo%fmlHandle, errorString, *999 )

    typeHandle = Fieldml_GetValueType( fmlInfo%fmlHandle, evaluatorHandle )
    CALL FieldmlUtil_CheckError( "Cannot get value type for geometric field", fmlInfo%fmlHandle, errorString, *999 )

    length = Fieldml_CopyObjectDeclaredName( fmlInfo%fmlHandle, typeHandle, name, BUFFER_SIZE )

    IF( INDEX( name, 'coordinates.rc.3d' ) == 1 ) THEN
      coordinateType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      coordinateCount = 3
    ELSE IF( INDEX( name, 'coordinates.rc.2d' ) == 1 ) THEN
      coordinateType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      coordinateCount = 2
    ELSE
      err = FML_ERR_UNKNOWN_COORDINATE_TYPE
      CALL FieldmlUtil_CheckErrorNumber( "Coordinate system not yet supported", err, errorString, *999 )
    ENDIF

    CALL COORDINATE_SYSTEM_CREATE_START( userNumber, coordinateSystem, err, errorString, *999 )
    !Set the coordinate system dimension and type
    CALL COORDINATE_SYSTEM_DIMENSION_SET( coordinateSystem, coordinateCount, err, errorString, *999 )
    CALL COORDINATE_SYSTEM_TYPE_SET( coordinateSystem, coordinateType, err, errorString, *999 )

    CALL EXITS( "FieldmlInput_CoordinateSystemCreateStart" )
    RETURN
999 CALL ERRORS( "FieldmlInput_CoordinateSystemCreateStart", err, errorString )
    CALL EXITS( "FieldmlInput_CoordinateSystemCreateStart" )
    RETURN 1

  END SUBROUTINE FieldmlInput_CoordinateSystemCreateStart


  !
  !================================================================================================================================
  !
  
  
  SUBROUTINE FieldmlInput_NodesCreateStart( fmlInfo, nodesArgumentName, region, nodes, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: nodesArgumentName
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region
    TYPE(NODES_TYPE), POINTER, INTENT(INOUT) :: nodes
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(C_INT) :: nodesArgumentHandle, nodesHandle, nodeCount
    
    nodesArgumentHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, nodesArgumentName//NUL )
    IF( nodesArgumentHandle == FML_INVALID_HANDLE ) THEN
      CALL FieldmlUtil_CheckErrorNumber( "Nodes argument name is invalid", err, errorString, *999 )
    END IF
    
    nodesHandle = Fieldml_GetValueType( fmlInfo%fmlHandle, nodesArgumentHandle )
    IF( nodesHandle == FML_INVALID_HANDLE ) THEN
      CALL FieldmlUtil_CheckErrorNumber( "Nodes argument type is invalid", err, errorString, *999 )
    END IF

    fmlInfo%nodesArgumentHandle = nodesArgumentHandle
    fmlInfo%nodesHandle = nodesHandle

    nodeCount = Fieldml_GetMemberCount( fmlInfo%fmlHandle, fmlInfo%nodesHandle )
    NULLIFY( nodes )
    CALL NODES_CREATE_START( region, nodeCount, nodes, err, errorString, *999 )

    CALL EXITS( "FieldmlInput_NodesCreateStart" )
    RETURN
999 CALL ERRORS( "FieldmlInput_NodesCreateStart", err, errorString )
    CALL EXITS( "FieldmlInput_NodesCreateStart" )
    RETURN 1

  END SUBROUTINE FieldmlInput_NodesCreateStart


  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_MeshCreateStart( fmlInfo, meshArgumentName, mesh, meshNumber, region, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: meshArgumentName
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT) :: mesh
    INTEGER(INTG), INTENT(IN) :: meshNumber
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(INTG) :: count
    INTEGER(C_INT) :: meshArgument, xiDimensions, elementCount

    CALL ENTERS( "FieldmlInput_MeshCreateStart", err, errorString, *999 )
    
    meshArgument = Fieldml_GetObjectByName( fmlInfo%fmlHandle, meshArgumentName//NUL )
    IF( meshArgument == FML_INVALID_HANDLE ) THEN
      CALL FieldmlUtil_CheckError( "Named mesh argument not found", fmlInfo, errorString, *999 )
    ENDIF

    fmlInfo%meshHandle = Fieldml_GetValueType( fmlInfo%fmlHandle, meshArgument )
    IF( fmlInfo%meshHandle == FML_INVALID_HANDLE ) THEN
      CALL FieldmlUtil_CheckError( "Invalid mesh argument", fmlInfo, errorString, *999 )
    ENDIF
    
    fmlInfo%elementsHandle = Fieldml_GetMeshElementsType( fmlInfo%fmlHandle, fmlInfo%meshHandle )
    fmlInfo%elementsArgumentHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, meshArgumentName//".element"//NUL )

    fmlInfo%xiHandle = Fieldml_GetMeshChartType( fmlInfo%fmlHandle, fmlInfo%meshHandle )
    fmlInfo%xiArgumentHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, meshArgumentName//".xi"//NUL )

    count = Fieldml_GetTypeComponentCount( fmlInfo%fmlHandle, fmlInfo%xiHandle )
    IF( ( count < 1 ) .OR. ( count > 3 ) ) THEN
      err = FML_ERR_UNKNOWN_MESH_XI
      CALL FieldmlUtil_CheckErrorNumber( "Mesh dimension cannot be greater than 3, or less than 1", err, errorString, *999 )
    ENDIF

    xiDimensions = Fieldml_GetTypeComponentCount( fmlInfo%fmlHandle, fmlInfo%xiHandle )
    elementCount = Fieldml_GetMemberCount( fmlInfo%fmlHandle, fmlInfo%elementsHandle )
    NULLIFY( mesh )
    CALL MESH_CREATE_START( meshNumber, region, xiDimensions, mesh, err, errorString, *999 )
    CALL MESH_NUMBER_OF_ELEMENTS_SET( mesh, elementCount, err, errorString, *999 )
    
    CALL EXITS( "FieldmlInput_MeshCreateStart" )
    RETURN
999 CALL ERRORS( "FieldmlInput_MeshCreateStart", err, errorString )
    CALL EXITS( "FieldmlInput_MeshCreateStart" )
    RETURN 1

  END SUBROUTINE FieldmlInput_MeshCreateStart

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_BasisCreateStart( fieldmlInfo, evaluatorName, userNumber, basis, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    INTEGER(INTG), INTENT(IN) :: userNumber
    TYPE(BASIS_TYPE), POINTER, INTENT(INOUT) :: basis
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(INTG) :: count, i
    INTEGER(C_INT) :: handle, connectivityHandle, layoutHandle
    INTEGER(C_INT), ALLOCATABLE :: tempHandles(:)
    INTEGER(INTG) :: basisType
    INTEGER(INTG), ALLOCATABLE :: basisInterpolations(:)
    INTEGER(INTG), ALLOCATABLE :: collapse(:)
    
    CALL ENTERS( "FieldmlInput_BasisCreateStart", err, errorString, *999 )

    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, evaluatorName//NUL )
    CALL FieldmlUtil_CheckError( "Named basis not found", fieldmlInfo, errorString, *999 )
    CALL FieldmlInput_GetBasisInfo( fieldmlInfo, handle, connectivityHandle, layoutHandle, basisType, basisInterpolations, &
      & collapse, err, errorString, *999 )
    
    IF( ALLOCATED( fieldmlInfo%basisHandles ) ) THEN
      count = SIZE( fieldmlInfo%basisHandles )
      DO i = 1, count
        IF( fieldmlInfo%basisHandles( i ) == handle ) THEN
          CALL FieldmlUtil_CheckError( "Named basis already created", fieldmlInfo, errorString, *999 )
        ENDIF
      END DO
      ALLOCATE( tempHandles( count ) )

      tempHandles(1:count) = fieldmlInfo%basisHandles(1:count)
      DEALLOCATE( fieldmlInfo%basisHandles )
      ALLOCATE( fieldmlInfo%basisHandles( count + 1 ) )
      fieldmlInfo%basisHandles(1:count) = tempHandles(1:count)

      tempHandles(1:count) = fieldmlInfo%basisConnectivityHandles(1:count)
      DEALLOCATE( fieldmlInfo%basisConnectivityHandles )
      ALLOCATE( fieldmlInfo%basisConnectivityHandles( count + 1 ) )
      fieldmlInfo%basisConnectivityHandles(1:count) = tempHandles(1:count)

      tempHandles(1:count) = fieldmlInfo%basisLayoutHandles(1:count)
      DEALLOCATE( fieldmlInfo%basisLayoutHandles )
      ALLOCATE( fieldmlInfo%basisLayoutHandles( count + 1 ) )
      fieldmlInfo%basisLayoutHandles(1:count) = tempHandles(1:count)

      DEALLOCATE( tempHandles )
    ELSE
      count = 0
      ALLOCATE( fieldmlInfo%basisHandles( 1 ) )
      ALLOCATE( fieldmlInfo%basisConnectivityHandles( 1 ) )
      ALLOCATE( fieldmlInfo%basisLayoutHandles( 1 ) )
    ENDIF
    
    count = count + 1
    fieldmlInfo%basisHandles( count ) = handle
    fieldmlInfo%basisConnectivityHandles( count ) = connectivityHandle
    fieldmlInfo%basisLayoutHandles( count ) = layoutHandle
    err = Fieldml_SetObjectInt( fieldmlInfo%fmlHandle, handle, userNumber )
  
    NULLIFY(basis)
    CALL BASIS_CREATE_START( userNumber, basis, err, errorString, *999 )
    CALL BASIS_TYPE_SET( basis, basisType, err, errorString, *999 )
    CALL BASIS_NUMBER_OF_XI_SET( basis, size( basisInterpolations ), err, errorString, *999 )
    CALL BASIS_INTERPOLATION_XI_SET( basis, basisInterpolations, err, errorString, *999 )
    CALL BASIS_COLLAPSED_XI_SET( basis, collapse, err, errorString, *999 )
    
    IF( ALLOCATED( basisInterpolations ) ) THEN
      DEALLOCATE( basisInterpolations )
    ENDIF
    IF( ALLOCATED( collapse ) ) THEN
      DEALLOCATE( collapse )
    ENDIF
    
    CALL EXITS( "FieldmlInput_BasisCreateStart" )
    RETURN
999 CALL ERRORS( "FieldmlInput_BasisCreateStart", err, errorString )
    CALL EXITS( "FieldmlInput_BasisCreateStart" )
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_InitialiseFromFile( fieldmlInfo, filename, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Locals
    INTEGER(C_INT) :: length, count, i
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    CALL ENTERS( "FieldmlInput_InitialiseFromFile", err, errorString, *999 )

    fieldmlInfo%fmlHandle = Fieldml_CreateFromFile( filename//NUL )
    fieldmlInfo%nodesHandle = FML_INVALID_HANDLE
    fieldmlInfo%meshHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodeDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    
    err = Fieldml_GetLastError( fieldmlInfo%fmlHandle )
    IF( err /= FML_ERR_NO_ERROR ) THEN
      count = Fieldml_GetErrorCount( fieldmlInfo%fmlHandle )
      DO i = 1,count
        length = Fieldml_CopyError( fieldmlInfo%fmlHandle, i, name, BUFFER_SIZE )
        WRITE(*,*) "FieldML parse error: "//name(1:length)
      ENDDO
      CALL FieldmlUtil_CheckErrorNumber( "Cannot create FieldML handle from file", err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlInput_InitialiseFromFile" )
    RETURN
999 CALL ERRORS( "FieldmlInput_InitialiseFromFile", err, errorString )
    CALL EXITS( "FieldmlInput_InitialiseFromFile" )
    RETURN 1
    
  END SUBROUTINE FieldmlInput_InitialiseFromFile

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_ReadOrder( fieldmlInfo, orderHandle, order, count, err, errorString, * )
    !Argument
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    INTEGER(C_INT), INTENT(IN) :: orderHandle
    INTEGER(C_INT), ALLOCATABLE, INTENT(INOUT) :: order(:)
    INTEGER(C_INT), INTENT(IN) :: count
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Locals
    INTEGER(C_INT) :: readerHandle, i, readCount
    INTEGER(C_INT), TARGET :: buffer(1)
    
    CALL ENTERS( "FieldmlInput_ReadOrder", err, errorString, *999 )

    IF( orderHandle == FML_INVALID_HANDLE ) THEN
      !This is permitted, and indeed common.
      RETURN
    ENDIF
    
    readerHandle = Fieldml_OpenReader( fieldmlInfo%fmlHandle, orderHandle )
    CALL FieldmlUtil_CheckError( "Cannot open order reader", fieldmlInfo, errorString, *999 )
    
    ALLOCATE( order(count) )
    DO i = 1, count
      readCount = Fieldml_ReadIntValues( fieldmlInfo%fmlHandle, readerHandle, C_LOC(buffer), 1 )
      IF( readCount /= 1 ) THEN
        err = FML_ERR_IO_READ_ERR
        DEALLOCATE( order )
        CALL FieldmlUtil_CheckErrorNumber( "Cannot open order reader", err, errorString, *999 )
      ENDIF
      order(i) = buffer(1)
    END DO
    
    CALL EXITS( "FieldmlInput_ReadOrder" )
    RETURN
999 CALL ERRORS( "FieldmlInput_ReadOrder", err, errorString )
    CALL EXITS( "FieldmlInput_ReadOrder" )
    RETURN 1

  END SUBROUTINE FieldmlInput_ReadOrder

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldInput_Reorder( inputBuffer, order, count, outputBuffer, err, errorString, * )
    !Argument
    INTEGER(C_INT), INTENT(IN) :: inputBuffer(:)
    INTEGER(C_INT), ALLOCATABLE, INTENT(IN) :: order(:)
    INTEGER(C_INT), INTENT(IN) :: count
    INTEGER(C_INT), INTENT(INOUT) :: outputBuffer(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
  
    !Locals
    INTEGER(C_INT) :: i
    
    CALL ENTERS( "FieldInput_Reorder", err, errorString, *999 )
    
    IF( ALLOCATED( order ) ) THEN
      DO i = 1,count
        outputBuffer( i ) = inputBuffer( order( i ) )
      ENDDO
    ELSE
      outputBuffer = inputBuffer
    ENDIF

    CALL EXITS( "FieldInput_Reorder" )
    RETURN
999 CALL ERRORS( "FieldInput_Reorder", err, errorString )
    CALL EXITS( "FieldInput_Reorder" )
    RETURN 1
  
  END SUBROUTINE FieldInput_Reorder

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_CreateMeshComponent( fieldmlInfo, mesh, componentNumber, evaluatorName, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: mesh
    INTEGER(INTG), INTENT(IN) :: componentNumber
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Locals
    INTEGER(C_INT) :: handle, basisReferenceHandle, connectivityHandle, layoutHandle, basisNumber, lastBasisHandle, count
    INTEGER(C_INT), ALLOCATABLE, TARGET :: nodesBuffer(:), rawBuffer(:)
    INTEGER(INTG), ALLOCATABLE :: tempHandles(:)
    INTEGER(INTG) :: componentCount, elementCount, knownBasisCount, maxBasisNodesCount, basisNodesCount
    INTEGER(INTG) :: elementNumber, knownBasisNumber
    INTEGER(C_INT), ALLOCATABLE :: connectivityReaders(:), connectivityCounts(:)
    TYPE(ArrayShimType), ALLOCATABLE :: connectivityOrders(:)
    INTEGER(C_INT) :: tPtr, dataSource, orderHandle
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: meshElements

    CALL ENTERS( "FieldmlInput_CreateMeshComponent", err, errorString, *999 )
    
    err = FML_ERR_NO_ERROR
    NULLIFY( basis )
    NULLIFY( meshElements )    

    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, evaluatorName//NUL )
    IF( .NOT. FieldmlInput_IsTemplateCompatible( fieldmlInfo%fmlHandle, handle, fieldmlInfo%elementsHandle, err ) ) THEN
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckErrorNumber( "Mesh component cannot be created from this evaluator", err, errorString, *999 )
    ENDIF

    IF( .NOT. ALLOCATED( fieldmlInfo%componentHandles ) ) THEN
      ALLOCATE( fieldmlInfo%componentHandles( componentNumber ) )
    ELSE IF( SIZE( fieldmlInfo%componentHandles ) < componentNumber ) THEN
      componentCount = SIZE( fieldmlInfo%componentHandles )
      ALLOCATE( tempHandles( componentCount ) )
      tempHandles(1:componentCount) = fieldmlInfo%componentHandles(1:componentCount)
      DEALLOCATE( fieldmlInfo%componentHandles )
      ALLOCATE( fieldmlInfo%componentHandles( componentNumber ) )
      fieldmlInfo%componentHandles(1:componentCount) = tempHandles(1:componentCount)
      fieldmlInfo%componentHandles(componentCount+1:componentNumber) = FML_INVALID_HANDLE
      DEALLOCATE( tempHandles )
    ENDIF
    
    fieldmlInfo%componentHandles( componentNumber ) = handle
    
    knownBasisCount = SIZE( fieldmlInfo%basisHandles )
    ALLOCATE( connectivityReaders( knownBasisCount ) )
    ALLOCATE( connectivityCounts( knownBasisCount ) )
    ALLOCATE( connectivityOrders( knownBasisCount ) )
    
    maxBasisNodesCount = 0
    DO knownBasisNumber = 1, knownBasisCount
      layoutHandle = fieldmlInfo%basisLayoutHandles(knownBasisNumber)
      connectivityHandle = fieldmlInfo%basisConnectivityHandles(knownBasisNumber)
        
      basisNodesCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, layoutHandle )
      CALL FieldmlUtil_CheckError( "Cannot get local node count for layout", fieldmlInfo, errorString, *999 )

      IF( basisNodesCount > maxBasisNodesCount ) THEN
        maxBasisNodesCount = basisNodesCount
      ENDIF
      
      orderHandle = Fieldml_GetSemidenseIndexOrder( fieldmlInfo%fmlHandle, connectivityHandle, 1 )
      CALL FieldmlInput_ReadOrder( fieldmlInfo, orderHandle, connectivityOrders( knownBasisNumber )%array, &
        & basisNodesCount, err, errorString, *999 )
    
      dataSource = Fieldml_GetDataSource( fieldmlInfo%fmlHandle, connectivityHandle )
      connectivityReaders(knownBasisNumber) = Fieldml_OpenReader( fieldmlInfo%fmlHandle, dataSource )
      connectivityCounts(knownBasisNumber) = basisNodesCount
      CALL FieldmlUtil_CheckError( "Cannot open connectivity reader", fieldmlInfo, errorString, *999 )
      
    END DO

    ALLOCATE( nodesBuffer( maxBasisNodesCount ) )
    ALLOCATE( rawBuffer( maxBasisNodesCount ) )

    elementCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get element count for mesh", fieldmlInfo, errorString, *999 )
    
    lastBasisHandle = FML_INVALID_HANDLE
    
    DO elementNumber = 1, elementCount
      basisReferenceHandle = Fieldml_GetElementEvaluator( fieldmlInfo%fmlHandle, handle, elementNumber, 1 )
      CALL FieldmlUtil_CheckError( "Cannot get element evaluator from mesh component", fieldmlInfo, errorString, *999 )
      
      IF( basisReferenceHandle /= lastBasisHandle ) THEN
        basisNumber = Fieldml_GetObjectInt( fieldmlInfo%fmlHandle, basisReferenceHandle )
        CALL FieldmlUtil_CheckError( "Cannot get basis user number for element evaluator", fieldmlInfo, errorString, *999 )
        CALL BASIS_USER_NUMBER_FIND( basisNumber, basis, err, errorString, *999 )
      ENDIF

      IF( elementNumber == 1 ) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START( mesh, componentNumber, basis, meshElements, err, errorString, *999 )
      ENDIF
      
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET( elementNumber, meshElements, basis, err, errorString, *999 )
      
      DO knownBasisNumber = 1, knownBasisCount
        basisNodesCount = connectivityCounts( knownBasisNumber )
        !BUGFIX Intel compiler will explode if we don't use a temporary variable
        tPtr = connectivityReaders(knownBasisNumber)
        count = Fieldml_ReadIntValues( fieldmlInfo%fmlHandle, tPtr, C_LOC(rawBuffer), basisNodesCount )
        IF( count /= basisNodesCount ) THEN
          err = FML_ERR_INVALID_CONNECTIVITY
          CALL FieldmlUtil_CheckErrorNumber( "Error reading connectivity", err, errorString, *999 )
        ENDIF
        IF( fieldmlInfo%basisHandles(knownBasisNumber) == basisReferenceHandle ) THEN
          CALL FieldInput_Reorder( rawBuffer, connectivityOrders(knownBasisNumber)%array, basisNodesCount, &
            & nodesBuffer, err, errorString, *999 )
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET( elementNumber, meshElements, nodesBuffer(1:basisNodesCount), &
            & err, errorString, *999 )
        ENDIF
      ENDDO
  
    END DO
    
    DO knownBasisNumber = 1, knownBasisCount
      !BUGFIX Intel compiler will explode if we don't use a temporary variable
      tPtr = connectivityReaders(knownBasisNumber)
      err = Fieldml_CloseReader( fieldmlInfo%fmlHandle, tPtr )
      CALL FieldmlUtil_CheckErrorNumber( "Error closing connectivity reader", err, errorString, *999 )
      IF( ALLOCATED( connectivityOrders( knownBasisNumber )%array ) ) THEN
        DEALLOCATE( connectivityOrders( knownBasisNumber )%array )
      ENDIF
    ENDDO
    
    DEALLOCATE( nodesBuffer )
    DEALLOCATE( connectivityReaders )
    DEALLOCATE( connectivityCounts )
    DEALLOCATE( connectivityOrders )
    
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH( meshElements, err, errorString, *999 )
    
    err = Fieldml_SetObjectInt( fieldmlInfo%fmlHandle, handle, componentNumber )

    CALL EXITS( "FieldmlInput_CreateMeshComponent" )
    RETURN
999 CALL ERRORS( "FieldmlInput_CreateMeshComponent", err, errorString )
    IF( ALLOCATED( nodesBuffer ) ) THEN
      DEALLOCATE( nodesBuffer )
    ENDIF
    IF( ALLOCATED( connectivityReaders ) ) THEN
      DEALLOCATE( connectivityReaders )
    ENDIF
    IF( ALLOCATED( connectivityCounts ) ) THEN
      DEALLOCATE( connectivityCounts )
    ENDIF
    IF( ALLOCATED( connectivityOrders ) ) THEN
      DO knownBasisNumber = 1, knownBasisCount
        IF( ALLOCATED( connectivityOrders( knownBasisNumber )%array ) ) THEN
          DEALLOCATE( connectivityOrders( knownBasisNumber )%array )
        ENDIF
      ENDDO
    
      DEALLOCATE( connectivityOrders )
    ENDIF
    
    CALL EXITS( "FieldmlInput_CreateMeshComponent" )
    RETURN 1

  ENDSUBROUTINE

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_FieldCreateStart( fieldmlInfo, region, decomposition, fieldNumber, field, evaluatorName, &
    & err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: decomposition
    INTEGER(INTG), INTENT(IN) :: fieldNumber
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: field
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Locals
    INTEGER(C_INT) :: fieldHandle, templateHandle, typeHandle
    INTEGER(INTG) :: componentNumber, templateComponentNumber, fieldDimensions

    CALL ENTERS( "FieldmlInput_FieldCreateStart", err, errorString, *999 )

    fieldHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, evaluatorName//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get named field evaluator", fieldmlInfo, errorString, *999 )
    typeHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, fieldHandle )
    CALL FieldmlUtil_CheckError( "Cannot get named field evaluator's value type", fieldmlInfo, errorString, *999 )
    fieldDimensions = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    CALL FieldmlUtil_CheckError( "Cannot get named field evaluator's component count", fieldmlInfo, errorString, *999 )
    
    IF( .NOT. FieldmlInput_IsFieldCompatible( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%elementsHandle, err ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckError( "Cannot interpret given evaluator as a field", fieldmlInfo, errorString, *999 )
    ENDIF

    NULLIFY( field )
    CALL FIELD_CREATE_START( fieldNumber, region, field, err, errorString, *999 )
    CALL FIELD_TYPE_SET( field, FIELD_GEOMETRIC_TYPE, err, errorString, *999 )
    CALL FIELD_MESH_DECOMPOSITION_SET( field, decomposition, err, errorString, *999 )
    CALL FIELD_SCALING_TYPE_SET( field, FIELD_NO_SCALING, err, errorString, *999 )

    DO componentNumber = 1, fieldDimensions
      templateHandle = Fieldml_GetElementEvaluator( fieldmlInfo%fmlHandle, fieldHandle, componentNumber, 1 )
      CALL FieldmlUtil_CheckError( "Cannot get field component evaluator", fieldmlInfo, errorString, *999 )

      templateComponentNumber = Fieldml_GetObjectInt( fieldmlInfo%fmlHandle, templateHandle )
      CALL FieldmlUtil_CheckError( "Cannot get field component mesh component number", fieldmlInfo, errorString, *999 )

      CALL FIELD_COMPONENT_MESH_COMPONENT_SET( field, FIELD_U_VARIABLE_TYPE, componentNumber, templateComponentNumber, &
        & err, errorString, *999 )
    ENDDO

    CALL EXITS( "FieldmlInput_FieldCreateStart" )
    RETURN
999 CALL ERRORS( "FieldmlInput_FieldCreateStart", err, errorString )
    CALL EXITS( "FieldmlInput_FieldCreateStart" )
    RETURN 1
  
  END SUBROUTINE FieldmlInput_FieldCreateStart

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_FieldNodalParametersUpdate( fieldmlInfo, evaluatorName, mesh, field, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: mesh
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: field
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Locals
    TYPE(NODES_TYPE), POINTER :: nodes
    INTEGER(C_INT) :: nodalDofsHandle, count, dataSource
    INTEGER(INTG) :: versionNumber,componentNumber, nodeNumber, fieldDimensions, meshNodeCount, gNode
    INTEGER(INTG) :: meshComponent
    LOGICAL :: nodeExists
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: buffer(:)
    INTEGER(C_INT) :: reader

    nodalDofsHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, evaluatorName//NUL )
    CALL FieldmlUtil_CheckError( "Cannot get nodal field dofs", fieldmlInfo, errorString, *999 )
  
    dataSource = Fieldml_GetDataSource( fieldmlInfo%fmlHandle, nodalDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get nodal data source", fieldmlInfo, errorString, *999 )

    reader = Fieldml_OpenReader( fieldmlInfo%fmlHandle, dataSource )
    CALL FieldmlUtil_CheckError( "Cannot open nodal dofs reader", fieldmlInfo, errorString, *999 )
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( field, FIELD_U_VARIABLE_TYPE, fieldDimensions, err, errorString, *999 )
    
    IF( reader /= FML_INVALID_HANDLE ) THEN
      ALLOCATE( buffer( fieldDimensions ) )
      
      !TODO Code assumes that the data is dense in both node and component indexes.
      NULLIFY( nodes )
      CALL REGION_NODES_GET( mesh%REGION, nodes, err, errorString, *999 )
      CALL NODES_NUMBER_OF_NODES_GET( nodes, meshNodeCount, err, errorString, *999 )
      CALL FieldmlUtil_CheckError( "Cannot get mesh nodes count", fieldmlInfo, errorString, *999 )
      DO nodeNumber = 1, meshNodeCount
        count = Fieldml_ReadDoubleValues( fieldmlInfo%fmlHandle, reader, C_LOC(buffer), fieldDimensions )
        IF( count /= fieldDimensions ) THEN
          err = FML_ERR_INVALID_READ
          CALL FieldmlUtil_CheckErrorNumber( "Cannot read nodal dofs for field components", err, errorString, *999 )
        ENDIF

        DO componentNumber = 1, fieldDimensions
          CALL FIELD_COMPONENT_MESH_COMPONENT_GET( field, FIELD_U_VARIABLE_TYPE, componentNumber, meshComponent, &
            & err, errorString, *999 )
          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS( mesh, meshComponent, nodeNumber, nodeExists, gNode, &
            & err, errorString, *999 )
  
          IF( nodeExists ) THEN
            !Default to version 1 of each node derivative (value hardcoded in loop)
            versionNumber = 1
            CALL FIELD_PARAMETER_SET_UPDATE_NODE( field, FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE, versionNumber, &
              & NO_GLOBAL_DERIV, nodeNumber, componentNumber, buffer( componentNumber ), err, errorString, *999 )
          ENDIF
        ENDDO
      ENDDO
    
      DEALLOCATE( buffer )
  
      err = Fieldml_CloseReader( fieldmlInfo%fmlHandle, reader )
      CALL FieldmlUtil_CheckErrorNumber( "Cannot close nodal dofs reader", err, errorString, *999 )

      CALL FIELD_PARAMETER_SET_UPDATE_START( field, FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE, err, errorString, *999 )
      CALL FIELD_PARAMETER_SET_UPDATE_FINISH( field, FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE, err, errorString, *999 )
    ENDIF

    !TODO Set element and constant parameters
    
    CALL EXITS( "FieldmlInput_FieldNodalParametersUpdate" )
    RETURN
999 CALL ERRORS( "FieldmlInput_FieldNodalParametersUpdate", err, errorString )
    CALL EXITS( "FieldmlInput_FieldNodalParametersUpdate" )
    RETURN 1
  
  END SUBROUTINE FieldmlInput_FieldNodalParametersUpdate

  !
  !================================================================================================================================
  !

END MODULE FIELDML_INPUT_ROUTINES
