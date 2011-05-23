!> \file
!> $Id$
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


  !Interfaces

  INTERFACE

  END INTERFACE

  PUBLIC :: FieldmlInput_InitialiseFromFile, FieldmlInput_SetDofVariables, FieldmlInput_ReadMeshInfo, &
    & FieldmlInput_GetCoordinateSystemInfo, FieldmlInput_CreateBasis, FieldmlInput_CreateMeshComponent, &
    & FieldmlInput_CreateField

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_GetBasisConnectivityInfo( fmlInfo, basisHandle, connectivityHandle, layoutHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: basisHandle
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Local variables
    INTEGER(C_INT) :: count, paramsHandle, handle1, handle2
    
    CALL ENTERS( "FieldmlInput_GetBasisConnectivityInfo", err, errorString, *999 )

    count = Fieldml_GetBindCount( fmlInfo%fmlHandle, basisHandle )
    IF( count /= 2 ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Library basis evaluators must have exactly two binds", err, errorString, *999 )
    END IF
    
    handle1 = Fieldml_GetBindEvaluator( fmlInfo%fmlHandle, basisHandle, 1 )
    CALL FieldmlUtil_CheckError( "Cannot get first bind for FEM evaluator", fmlInfo%fmlHandle, errorString, *999 )
    handle2 = Fieldml_GetBindEvaluator( fmlInfo%fmlHandle, basisHandle, 2 )
    CALL FieldmlUtil_CheckError( "Cannot get second bind for FEM evaluator", fmlInfo%fmlHandle, errorString, *999 )

    IF( handle1 == fmlInfo%xiVariableHandle ) THEN
      paramsHandle = handle2
    ELSE IF( handle2 == fmlInfo%xiVariableHandle ) THEN
      paramsHandle = handle1
    ELSE
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Library FEM evaluators must a xi bind", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( fmlInfo%fmlHandle, paramsHandle ) /= FHT_AGGREGATE_EVALUATOR ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Parameter evaluator for interpolator must be a reference", err, errorString, *999 )
    ENDIF
    
    count = Fieldml_GetBindCount( fmlInfo%fmlHandle, paramsHandle )
    IF( count /= 1 ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Nodal parameter evaluator must only have one bind", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetBindVariable( fmlInfo%fmlHandle, paramsHandle, 1 ) /= fmlInfo%nodesVariableHandle ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Nodal parameter evaluator must bind the nodes variable", err, errorString, *999 )
    ENDIF
    
    connectivityHandle = Fieldml_GetBindEvaluator( fmlInfo%fmlHandle, paramsHandle, 1 )
    CALL FieldmlUtil_CheckError( "Cannot get connectivity source for nodal parameters", fmlInfo%fmlHandle, &
      & errorString, *999 )
      
    handle1 = Fieldml_GetIndexEvaluator( fmlInfo%fmlHandle, paramsHandle, 1 )
    layoutHandle = Fieldml_GetValueType( fmlInfo%fmlHandle, handle1 )
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

  SUBROUTINE FieldmlInput_GetBasisInfo( fmlInfo, objectHandle, connectivityHandle, basisType, basisInterpolations, collapse, &
    & err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: objectHandle
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle
    INTEGER(INTG), INTENT(OUT) :: basisType
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: basisInterpolations(:)
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: collapse(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(C_INT) :: length, layoutHandle, libraryBasisHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    CALL ENTERS( "FieldmlInput_GetBasisInfo", err, errorString, *999 )

    IF( .NOT. FieldmlInput_IsKnownBasis( fmlInfo, objectHandle, err ) ) THEN
      CALL FieldmlUtil_CheckError( "Basis specified in FieldML file is not yet supported", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( fmlInfo%fmlHandle, objectHandle ) /= FHT_REFERENCE_EVALUATOR ) THEN
      err = FML_ERR_INVALID_BASIS
      CALL FieldmlUtil_CheckError( "Basis evaluator must be a continuous reference", err, errorString, *999 )
    ENDIF
    
    libraryBasisHandle = Fieldml_GetReferenceRemoteEvaluator( fmlInfo%fmlHandle, objectHandle )
    CALL FieldmlUtil_CheckError( "Basis specified in FieldML is not a reference evaluator", fmlInfo, errorString, *999 )
    length = Fieldml_CopyObjectName( fmlInfo%fmlHandle, libraryBasisHandle, name, BUFFER_SIZE )
    CALL FieldmlUtil_CheckError( "Cannot get name of basis evaluator", fmlInfo, errorString, *999 )

    IF( INDEX( name, 'library.fem.triquadratic_lagrange') == 1 ) THEN
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      CALL REALLOCATE_INT( collapse, 3, "", err, errorString, *999 )
      basisInterpolations = BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
      basisType = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( name, 'library.fem.trilinear_lagrange') == 1 ) THEN
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      CALL REALLOCATE_INT( collapse, 3, "", err, errorString, *999 )
      basisInterpolations = BASIS_LINEAR_LAGRANGE_INTERPOLATION
      basisType = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE
      err = FML_ERR_UNKNOWN_BASIS
      CALL FieldmlUtil_CheckError( "Basis cannot yet be interpreted", err, errorString, *999 )
    ENDIF
    
    IF( basisType == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      CALL FieldmlInput_GetBasisCollapse( name(1:length), collapse )
    ENDIF
    
    CALL FieldmlInput_GetBasisConnectivityInfo( fmlInfo, objectHandle, connectivityHandle, layoutHandle, err, errorString, *999 )

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
    TYPE(VARYING_STRING) :: errorString

    !Locals
    INTEGER(C_INT) :: length, connectivityHandle, layoutHandle, libraryBasisHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    FieldmlInput_IsKnownBasis = .FALSE.

    IF( Fieldml_GetObjectType( fmlInfo%fmlHandle, objectHandle ) /= FHT_REFERENCE_EVALUATOR ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF

    libraryBasisHandle = Fieldml_GetReferenceRemoteEvaluator( fmlInfo%fmlHandle, objectHandle )
    length = Fieldml_CopyObjectName( fmlInfo%fmlHandle, libraryBasisHandle, name, BUFFER_SIZE )

    IF( ( INDEX( name, 'library.fem.triquadratic_lagrange') /= 1 ) .AND. &
      & ( INDEX( name, 'library.fem.trilinear_lagrange') /= 1 ) ) THEN
      err = FML_ERR_UNKNOWN_BASIS
      RETURN
    ENDIF
    
    CALL FieldmlInput_GetBasisConnectivityInfo( fmlInfo, objectHandle, connectivityHandle, layoutHandle, err, errorString, *999 )
    IF( connectivityHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF
    
    FieldmlInput_IsKnownBasis = .TRUE.
    RETURN

    err = FML_ERR_INVALID_BASIS
999 RETURN
    
  END FUNCTION FieldmlInput_IsKnownBasis
  
  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_IsElementEvaluatorCompatible( fmlHandle, object, err )
    !Arguments
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
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

    evaluatorHandle = Fieldml_GetReferenceRemoteEvaluator( fmlHandle, object )
    length = Fieldml_CopyObjectName( fmlHandle, evaluatorHandle, name, BUFFER_SIZE )

    IF( INDEX( name, 'library.fem.trilinear_lagrange' ) == 1 ) THEN
      FieldmlInput_IsElementEvaluatorCompatible = .TRUE.
    ELSE IF( INDEX( name, 'library.fem.triquadratic_lagrange' ) == 1 ) THEN
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
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
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
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
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
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), ALLOCATABLE :: fieldHandles(:)
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    INTEGER(C_INT) :: count, i, object, fieldCount, elementType

    CALL ENTERS( "Fieldml_GetFieldHandles", err, errorString, *999 )

    elementType = Fieldml_GetMeshElementType( fmlHandle, meshHandle )
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

  SUBROUTINE FieldmlInput_GetCoordinateSystemInfo( fmlInfo, evaluatorHandle, coordinateType, coordinateCount, &
    & err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: evaluatorHandle
    INTEGER(INTG), INTENT(OUT) :: coordinateType
    INTEGER(INTG), INTENT(OUT) :: coordinateCount
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(C_INT) :: typeHandle, length
    CHARACTER(LEN=BUFFER_SIZE) :: name

    CALL ENTERS( "FieldmlInput_GetCoordinateSystemInfo", err, errorString, *999 )

    coordinateType = 0 !There doesn't seem to be a COORDINATE_UNKNOWN_TYPE

    typeHandle = Fieldml_GetValueType( fmlInfo%fmlHandle, evaluatorHandle )
    CALL FieldmlUtil_CheckError( "Cannot get value type for geometric field", fmlInfo%fmlHandle, errorString, *999 )

    length = Fieldml_CopyObjectName( fmlInfo%fmlHandle, typeHandle, name, BUFFER_SIZE )

    IF( INDEX( name, 'library.coordinates.rc.3d' ) == 1 ) THEN
      coordinateType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      coordinateCount = 3
    ELSE IF( INDEX( name, 'library.coordinates.rc.2d' ) == 1 ) THEN
      coordinateType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      coordinateCount = 2
    ELSE
      err = FML_ERR_UNKNOWN_COORDINATE_TYPE
      CALL FieldmlUtil_CheckError( "Coordinate system not yet supported", err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlInput_GetCoordinateSystemInfo" )
    RETURN
999 CALL ERRORS( "FieldmlInput_GetCoordinateSystemInfo", err, errorString )
    CALL EXITS( "FieldmlInput_GetCoordinateSystemInfo" )
    RETURN 1

  END SUBROUTINE FieldmlInput_GetCoordinateSystemInfo


  !
  !================================================================================================================================
  !
  
  
  SUBROUTINE FieldmlInput_GetConnectivityInfo( fmlInfo, connectivityHandle, layoutHandle, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fmlInfo
    INTEGER(C_INT), INTENT(IN) :: connectivityHandle
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(C_INT) :: evaluatorHandle, type1, type2

    CALL ENTERS( "FieldmlInput_GetConnectivityInfo", err, errorString, *999 )

    IF( fmlInfo%nodesHandle /= Fieldml_GetValueType( fmlInfo%fmlHandle, connectivityHandle ) ) THEN
      err = FML_ERR_INVALID_CONNECTIVITY
      CALL FieldmlUtil_CheckError( "Connectivity evaluator must vary over the correct point ensemble", &
        & err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetIndexCount( fmlInfo%fmlHandle, connectivityHandle ) /= 2 ) THEN
      err = FML_ERR_INVALID_CONNECTIVITY
      CALL FieldmlUtil_CheckError( "Connectivity evaluator must only vary over two ensembles", err, errorString, *999 )
    END IF

    evaluatorHandle = Fieldml_GetIndexEvaluator( fmlInfo%fmlHandle, connectivityHandle, 1 )
    type1 = Fieldml_GetValueType( fmlInfo%fmlHandle, evaluatorHandle )

    evaluatorHandle = Fieldml_GetIndexEvaluator( fmlInfo%fmlHandle, connectivityHandle, 2 )
    type2 = Fieldml_GetValueType( fmlInfo%fmlHandle, evaluatorHandle )

    IF( ( type1 /= fmlInfo%elementsHandle ) .AND. &
      & ( type2 /= fmlInfo%elementsHandle ) ) THEN
    END IF
    
    IF( type1 == fmlInfo%elementsHandle ) THEN
      layoutHandle = type2
    ELSEIF( type2 == fmlInfo%elementsHandle ) THEN
      layoutHandle = type1
    ELSE
      err = FML_ERR_INVALID_CONNECTIVITY
      CALL FieldmlUtil_CheckError( "Connectivity evaluator must vary over mesh elements type", err, errorString, *999 )
    ENDIF
      
    CALL EXITS( "FieldmlInput_GetConnectivityInfo" )
    RETURN
999 CALL ERRORS( "FieldmlInput_GetConnectivityInfo", err, errorString )
    CALL EXITS( "FieldmlInput_GetConnectivityInfo" )
    RETURN 1

  END SUBROUTINE FieldmlInput_GetConnectivityInfo


  !
  !================================================================================================================================
  !


  SUBROUTINE FieldmlInput_ReadMeshInfo( fmlInfo, meshName, nodesEnsembleName, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: meshName
    CHARACTER(LEN=*), INTENT(IN) :: nodesEnsembleName
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(INTG) :: count

    CALL ENTERS( "FieldmlInput_ReadMeshInfo", err, errorString, *999 )

    fmlInfo%meshHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, meshName//NUL )
    IF( fmlInfo%meshHandle == FML_INVALID_HANDLE ) THEN
      err = Fieldml_GetLastError( fmlInfo%fmlHandle )
      CALL FieldmlUtil_CheckError( "Named mesh cannot be found", err, errorString, *999 )
    ENDIF
    
    fmlInfo%elementsHandle = Fieldml_GetMeshElementType( fmlInfo%fmlHandle, fmlInfo%meshHandle )
    fmlInfo%elementsVariableHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, meshName//".variable.element"//NUL )

    fmlInfo%xiHandle = Fieldml_GetMeshXiType( fmlInfo%fmlHandle, fmlInfo%meshHandle )
    fmlInfo%xiVariableHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, meshName//".variable.xi"//NUL )

    count = Fieldml_GetTypeComponentCount( fmlInfo%fmlHandle, fmlInfo%xiHandle )
    IF( ( count < 1 ) .OR. ( count > 3 ) ) THEN
      err = FML_ERR_UNKNOWN_MESH_XI
      CALL FieldmlUtil_CheckError( "Mesh dimension cannot be greater than 3, or less than 1", err, errorString, *999 )
    ENDIF
    
    fmlInfo%nodesHandle = Fieldml_GetObjectByName( fmlInfo%fmlHandle, nodesEnsembleName//NUL )
    IF( fmlInfo%nodesHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_CONNECTIVITY
      CALL FieldmlUtil_CheckError( "No valid point ensemble found for mesh connectivity", err, errorString, *999 )
    END IF

    fmlInfo%nodesVariableHandle = FieldmlUtil_GetTypeVariableHandle( fmlInfo, fmlInfo%nodesHandle )

    CALL EXITS( "FieldmlInput_ReadMeshInfo" )
    RETURN
999 CALL ERRORS( "FieldmlInput_ReadMeshInfo", err, errorString )
    CALL EXITS( "FieldmlInput_ReadMeshInfo" )
    RETURN 1

  END SUBROUTINE FieldmlInput_ReadMeshInfo

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_CreateBasis( fieldmlInfo, userNumber, evaluatorName, gaussQuadrature, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    INTEGER(INTG), INTENT(IN) :: userNumber
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    INTEGER(INTG), INTENT(IN) :: gaussQuadrature(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString

    !Locals
    INTEGER(INTG) :: count, i
    INTEGER(C_INT) :: handle, connectivityHandle
    INTEGER(C_INT), ALLOCATABLE :: tempHandles(:)
    INTEGER(INTG) :: basisType
    INTEGER(INTG), ALLOCATABLE :: basisInterpolations(:)
    INTEGER(INTG), ALLOCATABLE :: collapse(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    
    CALL ENTERS( "FieldmlInput_CreateBasis", err, errorString, *999 )

    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, evaluatorName//NUL )
    CALL FieldmlUtil_CheckError( "Named basis not found", fieldmlInfo, errorString, *999 )
    CALL FieldmlInput_GetBasisInfo( fieldmlInfo, handle, connectivityHandle, basisType, basisInterpolations, collapse, &
      & err, errorString, *999 )
    
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

      DEALLOCATE( tempHandles )
    ELSE
      count = 0
      ALLOCATE( fieldmlInfo%basisHandles( 1 ) )
      ALLOCATE( fieldmlInfo%basisConnectivityHandles( 1 ) )
    ENDIF
    
    count = count + 1
    fieldmlInfo%basisHandles( count ) = handle
    fieldmlInfo%basisConnectivityHandles( count ) = connectivityHandle
    err = Fieldml_SetObjectInt( fieldmlInfo%fmlHandle, handle, userNumber )
  
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START( userNumber, basis, err, errorString, *999 )
    CALL BASIS_TYPE_SET( basis, basisType, err, errorString, *999 )
    CALL BASIS_NUMBER_OF_XI_SET( basis, size( basisInterpolations ), err, errorString, *999 )
    CALL BASIS_INTERPOLATION_XI_SET( basis, basisInterpolations, err, errorString, *999 )
    CALL BASIS_COLLAPSED_XI_SET( basis, collapse, err, errorString, *999 )
    IF( SIZE( gaussQuadrature ) > 0 ) THEN
      CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET( basis, gaussQuadrature, err, errorString, *999 )
    ENDIF
    CALL BASIS_CREATE_FINISH( basis, err, errorString, *999 )
    
    IF( ALLOCATED( basisInterpolations ) ) THEN
      DEALLOCATE( basisInterpolations )
    ENDIF
    IF( ALLOCATED( collapse ) ) THEN
      DEALLOCATE( collapse )
    ENDIF
  
    CALL EXITS( "FieldmlInput_CreateBasis" )
    RETURN
999 CALL ERRORS( "FieldmlInput_CreateBasis", err, errorString )
    CALL EXITS( "FieldmlInput_CreateBasis" )
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
    
    CALL ENTERS( "FieldmlInput_InitialiseFromFile", err, errorString, *999 )

    fieldmlInfo%fmlHandle = Fieldml_CreateFromFile( filename//NUL )
    fieldmlInfo%nodesHandle = FML_INVALID_HANDLE
    fieldmlInfo%meshHandle = FML_INVALID_HANDLE
    fieldmlInfo%elementsHandle = FML_INVALID_HANDLE
    fieldmlInfo%xiHandle = FML_INVALID_HANDLE
    fieldmlInfo%nodeDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%elementDofsHandle = FML_INVALID_HANDLE
    !fieldmlInfo%constantDofsHandle = FML_INVALID_HANDLE
    
    CALL FieldmlUtil_CheckError( "Cannot create FieldML handle from file", fieldmlInfo, errorString, *999 )

    CALL EXITS( "FieldmlInput_InitialiseFromFile" )
    RETURN
999 CALL ERRORS( "FieldmlInput_InitialiseFromFile", err, errorString )
    CALL EXITS( "FieldmlInput_InitialiseFromFile" )
    RETURN 1
    
  END SUBROUTINE FieldmlInput_InitialiseFromFile

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
    INTEGER(C_INT), ALLOCATABLE, TARGET :: nodesBuffer(:)
    INTEGER(INTG), ALLOCATABLE :: tempHandles(:)
    INTEGER(INTG) :: componentCount, elementCount, knownBasisCount, maxBasisNodesCount, basisNodesCount
    INTEGER(INTG) :: elementNumber, knownBasisNumber
    TYPE(C_PTR), ALLOCATABLE :: connectivityReaders(:)
    INTEGER(C_INT), ALLOCATABLE :: connectivityCounts(:)
    TYPE(C_PTR) :: tPtr
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: meshElements
    
    CALL ENTERS( "FieldmlInput_CreateMeshComponent", err, errorString, *999 )
    err = FML_ERR_NO_ERROR
    NULLIFY( basis )
    NULLIFY( meshElements )    

    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, evaluatorName//NUL )
    IF( .NOT. FieldmlInput_IsTemplateCompatible( fieldmlInfo%fmlHandle, handle, fieldmlInfo%elementsHandle, err ) ) THEN
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Mesh component cannot be created from this evaluator", err, errorString, *999 )
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
    
    maxBasisNodesCount = 0
    DO knownBasisNumber = 1, knownBasisCount
      CALL FieldmlInput_GetBasisConnectivityInfo( fieldmlInfo, &
        & fieldmlInfo%basisHandles(knownBasisNumber), connectivityHandle, layoutHandle, err, errorString, *999 )
        
      basisNodesCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, layoutHandle )
      CALL FieldmlUtil_CheckError( "Cannot get local node count for layout", fieldmlInfo, errorString, *999 )

      IF( basisNodesCount > maxBasisNodesCount ) THEN
        maxBasisNodesCount = basisNodesCount
      ENDIF
      
      connectivityReaders(knownBasisNumber) = Fieldml_OpenReader( fieldmlInfo%fmlHandle, connectivityHandle )
      connectivityCounts(knownBasisNumber) = basisNodesCount
      CALL FieldmlUtil_CheckError( "Cannot open connectivity reader", fieldmlInfo, errorString, *999 )
      
    END DO

    ALLOCATE( nodesBuffer( maxBasisNodesCount ) )

    elementCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get element count for mesh", fieldmlInfo, errorString, *999 )
    
    lastBasisHandle = FML_INVALID_HANDLE
    
    DO elementNumber = 1, elementCount
      basisReferenceHandle = Fieldml_GetElementEvaluator( fieldmlInfo%fmlHandle, handle, elementNumber, 1 )
      CALL FieldmlUtil_CheckError( "Cannot get element evaluator from mesh component", fieldmlInfo, errorString, *999 )
      
      IF( basisReferenceHandle /= lastBasisHandle ) THEN
        CALL FieldmlInput_GetBasisConnectivityInfo( fieldmlInfo, basisReferenceHandle, &
          & connectivityHandle, layoutHandle, err, errorString, *999 )
    
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
        count = Fieldml_ReadIntValues( fieldmlInfo%fmlHandle, tPtr, C_LOC(nodesBuffer), basisNodesCount )
        IF( count /= basisNodesCount ) THEN
          err = FML_ERR_INVALID_CONNECTIVITY
          CALL FieldmlUtil_CheckError( "Error reading connectivity", err, errorString, *999 )
        ENDIF
        IF( fieldmlInfo%basisHandles(knownBasisNumber) == basisReferenceHandle ) THEN
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET( elementNumber, meshElements, nodesBuffer(1:basisNodesCount), &
            & err, errorString, *999 )
        ENDIF
      ENDDO
  
    END DO
    
    DO knownBasisNumber = 1, knownBasisCount
      !BUGFIX Intel compiler will explode if we don't use a temporary variable
      tPtr = connectivityReaders(knownBasisNumber)
      err = Fieldml_CloseReader( fieldmlInfo%fmlHandle, tPtr )
      CALL FieldmlUtil_CheckError( "Error closing connectivity reader", err, errorString, *999 )
    ENDDO
    
    DEALLOCATE( nodesBuffer )
    DEALLOCATE( connectivityReaders )
    DEALLOCATE( connectivityCounts )
    
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
    CALL EXITS( "FieldmlInput_CreateMeshComponent" )
    RETURN 1

  ENDSUBROUTINE

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_SetDofVariables( fieldmlInfo, nodeDofsName, elementDofsName, constantDofsName, err, errorString, * )
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    CHARACTER(LEN=*), INTENT(IN) :: nodeDofsName
    CHARACTER(LEN=*), INTENT(IN) :: elementDofsName
    CHARACTER(LEN=*), INTENT(IN) :: constantDofsName
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    CALL ENTERS( "FieldmlInput_SetDofVariables", err, errorString, *999 )

    !Some of these may not actually exist, but that's OK because that means they're not used.
    fieldmlInfo%nodeDofsHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, nodeDofsName//NUL )
    !fieldmlInfo%elementDofsHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, elementDofsName//NUL )
    !fieldmlInfo%constantDofsHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, constantDofsName//NUL )
    
    CALL EXITS( "FieldmlInput_SetDofVariables" )
    RETURN
999 CALL ERRORS( "FieldmlInput_SetDofVariables", err, errorString )
    CALL EXITS( "FieldmlInput_SetDofVariables" )
    RETURN 1
  
  END SUBROUTINE

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_CreateField( fieldmlInfo, region, mesh, decomposition, fieldNumber, field, evaluatorName, &
    & err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: mesh
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: decomposition
    INTEGER(INTG), INTENT(IN) :: fieldNumber
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: field
    CHARACTER(LEN=*), INTENT(IN) :: evaluatorName
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString
    
    !Locals
    INTEGER(C_INT) :: fieldHandle, templateHandle, nodalDofsHandle, typeHandle, count
    INTEGER(INTG) :: versionNumber,componentNumber, templateComponentNumber, nodeNumber, fieldDimensions, meshNodeCount, gNode
    INTEGER(INTG), ALLOCATABLE :: componentNumbers(:)
    LOGICAL :: nodeExists
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: reader

    CALL ENTERS( "FieldmlInput_CreateField", err, errorString, *999 )

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

    ALLOCATE( componentNumbers( fieldDimensions ) )
    DO componentNumber = 1, fieldDimensions
      templateHandle = Fieldml_GetElementEvaluator( fieldmlInfo%fmlHandle, fieldHandle, componentNumber, 1 )
      CALL FieldmlUtil_CheckError( "Cannot get field component evaluator", fieldmlInfo, errorString, *999 )

      templateComponentNumber = Fieldml_GetObjectInt( fieldmlInfo%fmlHandle, templateHandle )
      CALL FieldmlUtil_CheckError( "Cannot get field component mesh component number", fieldmlInfo, errorString, *999 )

      CALL FIELD_COMPONENT_MESH_COMPONENT_SET( field, FIELD_U_VARIABLE_TYPE, componentNumber, templateComponentNumber, &
        & err, errorString, *999 )

      componentNumbers( componentNumber ) = templateComponentNumber
    ENDDO

    CALL FIELD_CREATE_FINISH( field, err, errorString, *999 )

    nodalDofsHandle = Fieldml_GetBindByVariable( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%nodeDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get nodal field dofs", fieldmlInfo, errorString, *999 )
  
    reader = Fieldml_OpenReader( fieldmlInfo%fmlHandle, nodalDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot open nodal dofs reader", fieldmlInfo, errorString, *999 )
    IF( C_ASSOCIATED( reader ) ) THEN
      ALLOCATE( buffer( fieldDimensions ) )
      
      meshNodeCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle )
      CALL FieldmlUtil_CheckError( "Cannot get mesh nodes count", fieldmlInfo, errorString, *999 )
      DO nodeNumber = 1, meshNodeCount
        count = Fieldml_ReadDoubleValues( fieldmlInfo%fmlHandle, reader, C_LOC(buffer), fieldDimensions )
        IF( count /= fieldDimensions ) THEN
          err = FML_ERR_INVALID_READ
          CALL FieldmlUtil_CheckError( "Cannot read nodal dofs for field components", err, errorString, *999 )
        ENDIF

        DO componentNumber = 1, fieldDimensions
          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS( mesh, componentNumbers( componentNumber ), nodeNumber, nodeExists, gNode, &
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
      CALL FieldmlUtil_CheckError( "Cannot close nodal dofs reader", err, errorString, *999 )

      CALL FIELD_PARAMETER_SET_UPDATE_START( field, FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE, err, errorString, *999 )
      CALL FIELD_PARAMETER_SET_UPDATE_FINISH( field, FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE, err, errorString, *999 )
    ENDIF

    !TODO Set element and constant parameters
    
    DEALLOCATE( componentNumbers )

    CALL EXITS( "FieldmlInput_CreateField" )
    RETURN
999 CALL ERRORS( "FieldmlInput_CreateField", err, errorString )
    CALL EXITS( "FieldmlInput_CreateField" )
    RETURN 1
  
  END SUBROUTINE

  !
  !================================================================================================================================
  !

END MODULE FIELDML_INPUT_ROUTINES
