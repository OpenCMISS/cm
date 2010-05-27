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
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
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
  USE OPENCMISS
  USE UTIL_ARRAY

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

  TYPE(VARYING_STRING) :: errorString

  !Interfaces

  INTERFACE FieldmlInput_ReadRawData
    MODULE PROCEDURE FieldmlInput_ReadRawData_Int
    MODULE PROCEDURE FieldmlInput_ReadRawData_Real
  END INTERFACE FieldmlInput_ReadRawData

  INTERFACE

  END INTERFACE

  PUBLIC :: FieldmlInput_GetMeshInfo, FieldmlInput_GetCoordinateSystemInfo, FieldmlInput_GetBasisInfo, &
    & Fieldml_GetFieldHandles, FieldmlInput_GetComponentBasis, FieldmlInput_GetBasisConnectivityInfo, &
    & FieldmlInput_ReadRawData, FieldmlInput_GetBasisHandles

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlInput_ReadRawData_Int( fmlHandle, parametersHandle, array, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    INTEGER(C_INT), INTENT(INOUT) :: array(:,:)
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: i, indexCount, count1, count2, handle1, handle2
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(C_INT), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: reader
    
    indexCount = Fieldml_GetSemidenseIndexCount( fmlHandle, parametersHandle, 1 )
    IF( indexCount /= 0 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    indexCount = Fieldml_GetSemidenseIndexCount( fmlHandle, parametersHandle, 0 )
    IF( indexCount /= 2 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    handle1 = Fieldml_GetSemidenseIndex( fmlHandle, parametersHandle, 1, 0 )
    handle2 = Fieldml_GetSemidenseIndex( fmlHandle, parametersHandle, 2, 0 )
    IF( ( handle1 == FML_INVALID_HANDLE ) .OR. ( handle2 == FML_INVALID_HANDLE ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    count1 = Fieldml_GetEnsembleDomainElementCount( fmlHandle, handle1 )
    count2 = Fieldml_GetEnsembleDomainElementCount( fmlHandle, handle2 )

    ALLOCATE( buffer( count1 ) )

    reader = Fieldml_OpenReader( fmlHandle, parametersHandle )
    IF( .NOT. C_ASSOCIATED( reader ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    DO i = 1, count2
      err = Fieldml_ReadIntSlice( fmlHandle, reader, C_LOC(dummy), C_LOC(buffer) )
      array( i, 1:count1 ) = buffer( 1:count1 )
    ENDDO
    
    err = Fieldml_CloseReader( fmlHandle, reader )

    DEALLOCATE( buffer )
    
  END SUBROUTINE FieldmlInput_ReadRawData_Int

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_ReadRawData_Real( fmlHandle, parametersHandle, array, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    REAL(C_DOUBLE), INTENT(INOUT) :: array(:,:)
    INTEGER(INTG), INTENT(OUT) :: err

    INTEGER(C_INT) :: i, indexCount, count1, count2, handle1, handle2
    INTEGER(C_INT), TARGET :: dummy(0)
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: reader
    
    indexCount = Fieldml_GetSemidenseIndexCount( fmlHandle, parametersHandle, 1 )
    IF( indexCount /= 0 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    indexCount = Fieldml_GetSemidenseIndexCount( fmlHandle, parametersHandle, 0 )
    IF( indexCount /= 2 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    handle1 = Fieldml_GetSemidenseIndex( fmlHandle, parametersHandle, 1, 0 )
    handle2 = Fieldml_GetSemidenseIndex( fmlHandle, parametersHandle, 2, 0 )
    IF( ( handle1 == FML_INVALID_HANDLE ) .OR. ( handle2 == FML_INVALID_HANDLE ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    count1 = Fieldml_GetEnsembleDomainElementCount( fmlHandle, handle1 )
    count2 = Fieldml_GetEnsembleDomainElementCount( fmlHandle, handle2 )

    ALLOCATE( buffer( count1 ) )

    reader = Fieldml_OpenReader( fmlHandle, parametersHandle )
    IF( .NOT. C_ASSOCIATED( reader ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    DO i = 1, count2
      err = Fieldml_ReadDoubleSlice( fmlHandle, reader, C_LOC(dummy), C_LOC(buffer) )
      array( i, 1:count1 ) = buffer( 1:count1 )
    ENDDO
    
    err = Fieldml_CloseReader( fmlHandle, reader )

    DEALLOCATE( buffer )
    
  END SUBROUTINE FieldmlInput_ReadRawData_Real

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_GetBasisConnectivityInfo( fmlHandle, meshHandle, basisHandle, connectivityHandle, layoutHandle, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle !<The FieldML handle
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(C_INT), INTENT(IN) :: basisHandle
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    
    !Local variables
    INTEGER(C_INT) :: count, i, xiHandle, paramsHandle, handle1, handle2
    CHARACTER(LEN=BUFFER_SIZE) :: name
    INTEGER(C_INT) :: length
    
    count = Fieldml_GetAliasCount( fmlHandle, basisHandle )
    IF( count /= 2 ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    END IF
    
    xiHandle = Fieldml_GetMeshXiDomain( fmlHandle, meshHandle )

    handle1 = Fieldml_GetAliasLocal( fmlHandle, basisHandle, 1 )
    handle2 = Fieldml_GetAliasLocal( fmlHandle, basisHandle, 2 )

    IF( handle1 == xiHandle ) THEN
      paramsHandle = handle2
    ELSE IF( handle2 == xiHandle ) THEN
      paramsHandle = handle1
    ELSE
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF

    IF( Fieldml_GetObjectType( fmlHandle, paramsHandle ) /= FHT_CONTINUOUS_IMPORT ) THEN
      err = FML_ERR_INVALID_BASIS
      length = Fieldml_CopyObjectName( fmlHandle, paramsHandle, name, BUFFER_SIZE )
      RETURN
    ENDIF
    
    count = Fieldml_GetAliasCount( fmlHandle, paramsHandle )
    IF( count /= 1 ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF

    handle1 = Fieldml_GetAliasLocal( fmlHandle, paramsHandle, 1 )
    
    count = Fieldml_GetMeshConnectivityCount( fmlHandle, meshHandle )
    DO i = 1, count
      IF( Fieldml_GetMeshConnectivitySource( fmlHandle, meshHandle, i ) == handle1 ) THEN
        connectivityHandle = handle1
        layoutHandle = Fieldml_GetMeshConnectivityDomain( fmlHandle, meshHandle, i )
        RETURN
      ENDIF
    ENDDO

    err = FML_ERR_INVALID_BASIS
    RETURN
  END SUBROUTINE FieldmlInput_GetBasisConnectivityInfo

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_GetBasisInfo( fmlHandle, meshHandle, objectHandle, basisType, basisInterpolations, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(C_INT), INTENT(IN) :: objectHandle
    INTEGER(INTG), INTENT(OUT) :: basisType
    INTEGER(C_INT), ALLOCATABLE, INTENT(OUT) :: basisInterpolations(:)
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: length, connectivityHandle, layoutHandle, libraryBasisHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    IF( .NOT. FieldmlInput_IsKnownBasis( fmlHandle, meshHandle, objectHandle, err ) ) THEN
      RETURN
    ENDIF

    libraryBasisHandle = Fieldml_GetImportRemoteEvaluator( fmlHandle, objectHandle )
    length = Fieldml_CopyObjectName( fmlHandle, libraryBasisHandle, name, BUFFER_SIZE )

    IF( Fieldml_GetObjectType( fmlHandle, objectHandle ) /= FHT_CONTINUOUS_IMPORT ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF
    
    IF( INDEX( name, 'library.fem.triquadratic_lagrange') == 1 ) THEN
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      basisInterpolations = CMISSBasisQuadraticLagrangeInterpolation
      basisType = CMISSBasisLagrangeHermiteTPType
    ELSE IF( INDEX( name, 'library.fem.trilinear_lagrange') == 1 ) THEN
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      basisInterpolations = CMISSBasisLinearLagrangeInterpolation
      basisType = CMISSBasisLagrangeHermiteTPType
    ELSE
      err = FML_ERR_UNKNOWN_BASIS
      RETURN
    ENDIF
    
    CALL FieldmlInput_GetBasisConnectivityInfo( fmlHandle, meshHandle, objectHandle, connectivityHandle, layoutHandle, err )
    IF( connectivityHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF

999 RETURN
    !Deliberately not finalized, so the user can make OpenCMISS-specific tweaks.
  END SUBROUTINE FieldmlInput_GetBasisInfo

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_IsKnownBasis( fmlHandle, meshHandle, objectHandle, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(C_INT), INTENT(IN) :: objectHandle
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Function
    LOGICAL :: FieldmlInput_IsKnownBasis

    !Locals
    INTEGER(C_INT) :: length, connectivityHandle, layoutHandle, libraryBasisHandle
    CHARACTER(LEN=BUFFER_SIZE) :: name
    
    FieldmlInput_IsKnownBasis = .FALSE.

    libraryBasisHandle = Fieldml_GetImportRemoteEvaluator( fmlHandle, objectHandle )
    length = Fieldml_CopyObjectName( fmlHandle, libraryBasisHandle, name, BUFFER_SIZE )

    IF( Fieldml_GetObjectType( fmlHandle, objectHandle ) /= FHT_CONTINUOUS_IMPORT ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF

    IF( ( INDEX( name, 'library.fem.triquadratic_lagrange') /= 1 ) .AND. &
      & ( INDEX( name, 'library.fem.trilinear_lagrange') /= 1 ) ) THEN
      err = FML_ERR_UNKNOWN_BASIS
      RETURN
    ENDIF
    
    CALL FieldmlInput_GetBasisConnectivityInfo( fmlHandle, meshHandle, objectHandle, connectivityHandle, layoutHandle, err )
    IF( connectivityHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_BASIS
      RETURN
    ENDIF
    
    FieldmlInput_IsKnownBasis = .TRUE.
    
  END FUNCTION FieldmlInput_IsKnownBasis
  
  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_GetBasisHandles( fmlHandle, meshHandle, bases, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle !<The FieldML handle
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(C_INT), ALLOCATABLE, INTENT(INOUT) :: bases(:) !<An array to hold the identified bases
    INTEGER(INTG), INTENT(OUT) :: err !<The error code

    !Local variables
    INTEGER(INTG) :: basisCount
    INTEGER(C_INT) :: objectHandle, i, count

    count = Fieldml_GetObjectCount( fmlHandle, FHT_CONTINUOUS_IMPORT )
    basisCount = 0

    DO i = 1, count
      objectHandle = Fieldml_GetObject( fmlHandle, FHT_CONTINUOUS_IMPORT, i )

      IF( .NOT. FieldmlInput_IsKnownBasis( fmlHandle, meshHandle, objectHandle, err ) ) THEN
        CYCLE
      ENDIF
      
      basisCount = basisCount + 1
      CALL GROW_ARRAY( bases, 1, "", err, errorString, *999 )
      bases( basisCount ) = objectHandle

999   CYCLE
    ENDDO

  END SUBROUTINE FieldmlInput_GetBasisHandles

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_HasMarkup( fmlHandle, object, attribute, value, err )
    !Arguments
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: object
    CHARACTER(LEN=*), INTENT(IN) :: attribute
    CHARACTER(LEN=*), INTENT(IN) :: value
    INTEGER(INTG), INTENT(OUT) :: err

    LOGICAL :: FieldmlInput_HasMarkup

    !Locals
    INTEGER(C_INT) :: length
    CHARACTER(LEN=BUFFER_SIZE) :: buffer

    length = Fieldml_CopyMarkupAttributeValue( fmlHandle, object, attribute//C_NULL_CHAR, buffer, BUFFER_SIZE )

    FieldmlInput_HasMarkup = ( INDEX( buffer, value ) == 1 )

    err = FML_ERR_NO_ERROR

  END FUNCTION FieldmlInput_HasMarkup

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
    IF( type /= FHT_CONTINUOUS_IMPORT ) THEN
      FieldmlInput_IsElementEvaluatorCompatible = .FALSE.
      RETURN
    ENDIF

    evaluatorHandle = Fieldml_GetImportRemoteEvaluator( fmlHandle, object )
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

  FUNCTION FieldmlInput_IsTemplateCompatible( fmlHandle, object, elementDomain, err )
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: object
    INTEGER(C_INT), INTENT(IN) :: elementDomain
    INTEGER(INTG), INTENT(OUT) :: err

    LOGICAL :: FieldmlInput_IsTemplateCompatible

    INTEGER(C_INT) :: type, count, i, evaluator, domain, firstEvaluator

    type = Fieldml_GetObjectType( fmlHandle, object )
    IF( type /= FHT_CONTINUOUS_PIECEWISE ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    domain = Fieldml_GetIndexDomain( fmlHandle, object, 1 )
    IF( domain /= elementDomain ) THEN
      FieldmlInput_IsTemplateCompatible = .TRUE.
      RETURN
    ENDIF

    count = Fieldml_GetEvaluatorCount( fmlHandle, object )

    IF( count == 0 ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    firstEvaluator = Fieldml_GetEvaluator( fmlHandle, object, 1 )
    IF( .NOT. FieldmlInput_IsElementEvaluatorCompatible( fmlHandle, firstEvaluator, err ) ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    !At the moment, OpenCMISS does not support different evaluators per element.

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

  FUNCTION FieldmlInput_IsFieldCompatible( fmlHandle, object, elementDomain, err )
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: object
    INTEGER(C_INT), INTENT(IN) :: elementDomain
    INTEGER(INTG), INTENT(OUT) :: err

    LOGICAL :: FieldmlInput_IsFieldCompatible

    INTEGER(C_INT) :: type, count, i, evaluator

    type = Fieldml_GetObjectType( fmlHandle, object )

    IF( type /= FHT_CONTINUOUS_AGGREGATE ) THEN
      FieldmlInput_IsFieldCompatible = .FALSE.
      RETURN
    ENDIF

    count = Fieldml_GetEvaluatorCount( fmlHandle, object )
    IF( count < 1 ) THEN
      FieldmlInput_IsFieldCompatible = .FALSE.
      RETURN
    ENDIF

    FieldmlInput_IsFieldCompatible = .TRUE.
    DO i = 1, count
      evaluator = Fieldml_GetEvaluator( fmlHandle, object, i )
      IF( .NOT. FieldmlInput_IsTemplateCompatible( fmlHandle, evaluator, elementDomain, err ) ) THEN
        FieldmlInput_IsFieldCompatible = .FALSE.
        RETURN
      ENDIF
    ENDDO

  END FUNCTION FieldmlInput_IsFieldCompatible

  !
  !================================================================================================================================
  !

  SUBROUTINE Fieldml_GetFieldHandles( fmlHandle, fieldHandles, meshHandle, err )
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), ALLOCATABLE :: fieldHandles(:)
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(INTG), INTENT(OUT) :: err

    INTEGER(C_INT) :: count, i, object, fieldCount, elementDomain

    elementDomain = Fieldml_GetMeshElementDomain( fmlHandle, meshHandle )

    fieldCount = 0
    count = Fieldml_GetObjectCount( fmlHandle, FHT_CONTINUOUS_AGGREGATE )
    DO i = 1, count
      object = Fieldml_GetObject( fmlHandle, FHT_CONTINUOUS_AGGREGATE, i )
      IF( .NOT. FieldmlInput_HasMarkup( fmlHandle, object, 'field', 'true', err ) ) THEN
        CYCLE
      ENDIF

      IF( .NOT. FieldmlInput_IsFieldCompatible( fmlHandle, object, elementDomain, err ) ) THEN
        CYCLE
      ENDIF

      CALL GROW_ARRAY( fieldHandles, 1, "", err, errorString, *999 )
      fieldCount = fieldCount + 1
      fieldHandles( fieldCount ) = object
    ENDDO

999 RETURN
  END SUBROUTINE Fieldml_GetFieldHandles


  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlInput_GetCoordinateSystemInfo( fmlHandle, evaluatorHandle, coordinateType, coordinateCount, err )
    !Arguments
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: evaluatorHandle
    INTEGER(INTG), INTENT(OUT) :: coordinateType
    INTEGER(INTG), INTENT(OUT) :: coordinateCount
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: domainHandle, length
    CHARACTER(LEN=BUFFER_SIZE) :: name

    domainHandle = Fieldml_GetValueDomain( fmlHandle, evaluatorHandle )

    IF( domainHandle == FML_INVALID_HANDLE ) THEN
      coordinateType = 0 !There doesn't seem to be a COORDINATE_UNKNOWN_TYPE
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    length = Fieldml_CopyObjectName( fmlHandle, domainHandle, name, BUFFER_SIZE )

    IF( INDEX( name, 'library.coordinates.rc.3d' ) == 1 ) THEN
      coordinateType = CMISSCoordinateRectangularCartesianType
      coordinateCount = 3
    ELSE IF( INDEX( name, 'library.coordinates.rc.2d' ) == 1 ) THEN
      coordinateType = CMISSCoordinateRectangularCartesianType
      coordinateCount = 2
    ELSE
      coordinateType = 0 !There doesn't seem to be a COORDINATE_UNKNOWN_TYPE
      err = FML_ERR_UNKNOWN_COORDINATE_TYPE
    ENDIF

  END SUBROUTINE FieldmlInput_GetCoordinateSystemInfo


  !
  !================================================================================================================================
  !


  SUBROUTINE FieldmlInput_GetMeshInfo( fmlHandle, meshHandle, dimensions, elementCount, nodeDomain, err )
    !Arguments
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(INTG), INTENT(OUT) :: dimensions
    INTEGER(INTG), INTENT(OUT) :: elementCount
    INTEGER(INTG), INTENT(OUT) :: nodeDomain
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(INTG) :: componentHandle, elementHandle, count, i, handle

    nodeDomain = FML_INVALID_HANDLE
    elementHandle = Fieldml_GetMeshElementDomain( fmlHandle, meshHandle )
    elementCount = Fieldml_GetEnsembleDomainElementCount( fmlHandle, elementHandle )
    handle = Fieldml_GetMeshXiDomain( fmlHandle, meshHandle )
    componentHandle = Fieldml_GetDomainComponentEnsemble( fmlHandle, handle )

    count = Fieldml_GetMeshConnectivityCount( fmlHandle, meshHandle )

    IF( count == 0 ) THEN
      err = FML_ERR_INVALID_MESH
      RETURN
    END IF

    DO i = 1, count
      handle = Fieldml_GetMeshConnectivitySource( fmlHandle, meshHandle, i )
      IF( Fieldml_GetObjectType( fmlHandle, handle ) /= FHT_ENSEMBLE_PARAMETERS ) THEN
        err = FML_ERR_INVALID_CONNECTIVITY
        RETURN
      END IF

      IF( Fieldml_GetIndexCount( fmlHandle, handle ) /= 2 ) THEN
        err = FML_ERR_INVALID_CONNECTIVITY
        RETURN
      END IF

      IF( ( Fieldml_GetIndexDomain( fmlHandle, handle, 1 ) /= elementHandle ) .AND. &
        & ( Fieldml_GetIndexDomain( fmlHandle, handle, 2 ) /= elementHandle ) ) THEN
        err = FML_ERR_INVALID_CONNECTIVITY
        RETURN
      END IF
      
      IF( i == 1 ) THEN
        nodeDomain = Fieldml_GetValueDomain( fmlHandle, handle )
        IF( .NOT. FieldmlInput_HasMarkup( fmlHandle, nodeDomain, "geometric", "point", err ) ) THEN
          err = FML_ERR_INVALID_CONNECTIVITY
          RETURN
        END IF      
      ELSE IF( nodeDomain /= Fieldml_GetValueDomain( fmlHandle, handle ) ) THEN
        err = FML_ERR_INVALID_CONNECTIVITY
        RETURN
      ENDIF

    END DO

    IF( nodeDomain == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_CONNECTIVITY
      RETURN
    END IF

    !At the moment, library domains do not actually exist, so we can't directly
    !ask them what their cardinality is.
    dimensions = Fieldml_GetEnsembleDomainElementCount( fmlHandle, componentHandle )
    IF( ( dimensions < 1 ) .OR. ( dimensions > 3 ) ) THEN
      dimensions = 0
      err = FML_ERR_UNKNOWN_MESH_XI
      RETURN
    ENDIF

  END SUBROUTINE FieldmlInput_GetMeshInfo

  !
  !================================================================================================================================
  !

  FUNCTION FieldmlInput_GetComponentBasis( fmlHandle, fieldHandle, componentNumber, err )
    !Arguments
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: fieldHandle
    INTEGER(INTG), INTENT(IN) :: componentNumber
    INTEGER(INTG), INTENT(OUT) :: err

    !Function
    INTEGER(C_INT) :: FieldmlInput_GetComponentBasis

    !Locals
    INTEGER(INTG) :: templateHandle

    FieldmlInput_GetComponentBasis = FML_INVALID_HANDLE

    IF( Fieldml_GetObjectType( fmlHandle, fieldHandle ) /= FHT_CONTINUOUS_AGGREGATE ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    END IF

    templateHandle = Fieldml_GetEvaluator( fmlHandle, fieldHandle, componentNumber )
    IF( templateHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_PARAMETER
      RETURN
    END IF

    IF( Fieldml_GetObjectType( fmlHandle, templateHandle ) /= FHT_CONTINUOUS_PIECEWISE ) THEN
      err = FML_ERR_INVALID_PARAMETER
      RETURN
    END IF

    !At the moment, we don't need an element number, as the evaluator for all elements is the same
    FieldmlInput_GetComponentBasis = Fieldml_GetEvaluator( fmlHandle, templateHandle, 1 )

  END FUNCTION FieldmlInput_GetComponentBasis

  !
  !================================================================================================================================
  !

END MODULE FIELDML_INPUT_ROUTINES
