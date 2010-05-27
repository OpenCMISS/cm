!> \file
!> $Id$
!> \author Caton Little
!> \brief This module handles writing out FieldML files.
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

!> Output routines for FieldML

MODULE FIELDML_OUTPUT_ROUTINES

  USE KINDS
  USE FIELDML_API
  USE FIELDML_UTIL_ROUTINES
  USE ISO_VARYING_STRING
  USE OPENCMISS
  USE STRINGS

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: BUFFER_SIZE = 1024
  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  !INTEGER(INTG), PARAMETER ::
  !INTEGER(INTG), PARAMETER ::
  !INTEGER(INTG), PARAMETER ::

  TYPE(VARYING_STRING) :: errorString

  !Interfaces
  INTERFACE FieldmlOutput_WriteRawData
    MODULE PROCEDURE FieldmlOutput_WriteRawData_Int2
    MODULE PROCEDURE FieldmlOutput_WriteRawData_Real1
    MODULE PROCEDURE FieldmlOutput_WriteRawData_Real2
  END INTERFACE FieldmlOutput_WriteRawData

  TYPE ConnectivityInfoType
    INTEGER(C_INT) :: connectivityHandle
    INTEGER(C_INT) :: layoutHandle
  END TYPE ConnectivityInfoType

  TYPE BasisInfoType
    INTEGER(INTG) :: basisNumber
    INTEGER(C_INT) :: connectivityHandle
    INTEGER(C_INT) :: importHandle
  END TYPE BasisInfoType

  PUBLIC :: FieldmlOutput_WriteRawData, FieldmlOutput_AddConnectivity

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_GetTPBasisEvaluator( fmlhandle, xiInterpolations, evaluatorHandle, parametersHandle, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlhandle
    INTEGER(C_INT), INTENT(IN) :: xiInterpolations(:)
    INTEGER(C_INT), INTENT(OUT) :: evaluatorHandle
    INTEGER(C_INT), INTENT(OUT) :: parametersHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: xiCount, firstInterpolation, i
    
    xiCount = SIZE( xiInterpolations )
  
    firstInterpolation = xiInterpolations(1)
    DO i = 2, xiCount
      IF( xiInterpolations(i) /= firstInterpolation ) THEN
        !Do not yet support inhomogeneous TP bases
        err = FML_ERR_INVALID_OBJECT
        RETURN
      ENDIF
    ENDDO

    evaluatorHandle = FML_INVALID_HANDLE
    parametersHandle = FML_INVALID_HANDLE
      
    IF( firstInterpolation == CMISSBasisQuadraticLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlhandle, "library.fem.quadratic_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlhandle, "library.parameters.quadratic_lagrange"//NUL )
      ELSE IF( xiCount == 2 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlhandle, "library.fem.biquadratic_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlhandle, "library.parameters.biquadratic_lagrange"//NUL )
      ELSE IF( xiCount == 3 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlhandle, "library.fem.triquadratic_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlhandle, "library.parameters.triquadratic_lagrange"//NUL )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
      ENDIF
    ELSE IF( firstInterpolation == CMISSBasisLinearLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlhandle, "library.fem.linear_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlhandle, "library.parameters.linear_lagrange"//NUL )
      ELSE IF( xiCount == 2 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlhandle, "library.fem.bilinear_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlhandle, "library.parameters.bilinear_lagrange"//NUL )
      ELSE IF( xiCount == 3 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlhandle, "library.fem.trilinear_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlhandle, "library.parameters.trilinear_lagrange"//NUL )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
      ENDIF
    ELSE
      err = FML_ERR_INVALID_OBJECT
    ENDIF

  END SUBROUTINE FieldmlUtil_GetTPBasisEvaluator

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_WriteRawData_Int2( fmlHandle, parametersHandle, array, append, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    INTEGER(C_INT), INTENT(INOUT) :: array(:,:)
    INTEGER(C_INT), INTENT(IN) :: append
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: i, indexCount, count1, count2, handle1, handle2
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(C_INT), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: writer
    
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

    writer = Fieldml_OpenWriter( fmlHandle, parametersHandle, append )
    IF( .NOT. C_ASSOCIATED( writer ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    DO i = 1, count2
      buffer( 1:count1 ) = array( i, 1:count1 )
      err = Fieldml_WriteIntSlice( fmlHandle, writer, C_LOC(dummy), C_LOC(buffer) )
    ENDDO
    
    err = Fieldml_CloseWriter( fmlHandle, writer )

    DEALLOCATE( buffer )
    
  END SUBROUTINE FieldmlOutput_WriteRawData_Int2

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_WriteRawData_Real1( fmlHandle, parametersHandle, array, append, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    REAL(C_DOUBLE), INTENT(IN), ALLOCATABLE, TARGET :: array(:)
    INTEGER(C_INT), INTENT(IN) :: append
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: indexCount, count1, handle1
    INTEGER(C_INT), TARGET :: dummy(0)
    TYPE(C_PTR) :: writer
    
    indexCount = Fieldml_GetSemidenseIndexCount( fmlHandle, parametersHandle, 1 )
    IF( indexCount /= 0 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    indexCount = Fieldml_GetSemidenseIndexCount( fmlHandle, parametersHandle, 0 )
    IF( indexCount /= 1 ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    handle1 = Fieldml_GetSemidenseIndex( fmlHandle, parametersHandle, 1, 0 )
    IF( handle1 == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF
    
    count1 = Fieldml_GetEnsembleDomainElementCount( fmlHandle, handle1 )

    writer = Fieldml_OpenWriter( fmlHandle, parametersHandle, append )
    IF( .NOT. C_ASSOCIATED( writer ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    err = Fieldml_WriteDoubleSlice( fmlHandle, writer, C_LOC(dummy), C_LOC(array) )
    
    err = Fieldml_CloseWriter( fmlHandle, writer )
    
  END SUBROUTINE FieldmlOutput_WriteRawData_Real1

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_WriteRawData_Real2( fmlHandle, parametersHandle, array, append, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: parametersHandle
    REAL(C_DOUBLE), INTENT(INOUT), TARGET :: array(:,:)
    INTEGER(C_INT), INTENT(IN) :: append
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: i, indexCount, count1, count2, handle1, handle2
    INTEGER(C_INT), TARGET :: dummy(0)
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: buffer(:)
    TYPE(C_PTR) :: writer
    
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

    writer = Fieldml_OpenWriter( fmlHandle, parametersHandle, append )
    IF( .NOT. C_ASSOCIATED( writer ) ) THEN
      err = FML_ERR_INVALID_OBJECT
      RETURN
    ENDIF

    DO i = 1, count2
      buffer( 1:count1 ) = array( i, 1:count1 )
      err = Fieldml_WriteDoubleSlice( fmlHandle, writer, C_LOC(dummy), C_LOC(buffer) )
    ENDDO
    
    err = Fieldml_CloseWriter( fmlHandle, writer )

    DEALLOCATE( buffer )
    
  END SUBROUTINE FieldmlOutput_WriteRawData_Real2

  !
  !================================================================================================================================
  !
  
  FUNCTION FieldmlOutput_FindLayout( connectivityInfo, layoutHandle )
    !Argument variables
    TYPE(ConnectivityInfoType), INTENT(IN) :: connectivityInfo(:)
    INTEGER(C_INT), INTENT(IN) :: layoutHandle
    
    !Function
    INTEGER(INTG) :: FieldmlOutput_FindLayout
    
    !Locals
    INTEGER(INTG) :: i
    
    FieldmlOutput_FindLayout = -1
    DO i = 1, SIZE( connectivityInfo )
      IF( connectivityInfo(i)%layoutHandle == layoutHandle ) THEN
        FieldmlOutput_FindLayout = i
      ENDIF
    ENDDO
  
  END FUNCTION FieldmlOutput_FindLayout
  
  !
  !================================================================================================================================
  !
  
  FUNCTION FieldmlOutput_FindBasis( basisInfo, basisNumber )
    !Argument variables
    TYPE(BasisInfoType), INTENT(IN) :: basisInfo(:)
    INTEGER(INTG), INTENT(IN) :: basisNumber
    
    !Function
    INTEGER(INTG) :: FieldmlOutput_FindBasis
    
    !Locals
    INTEGER(INTG) :: i
    
    FieldmlOutput_FindBasis = -1
    DO i = 1, SIZE( basisInfo )
      IF( basisInfo(i)%basisNumber == basisNumber ) THEN
        FieldmlOutput_FindBasis = i
      ENDIF
    ENDDO
  
  END FUNCTION FieldmlOutput_FindBasis
  
  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_GetSimpleLayoutName( fmlHandle, layoutHandle, name, length, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: layoutHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: name
    INTEGER(C_INT) :: length
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: fullName
    
    length = Fieldml_CopyObjectName( fmlHandle, layoutHandle, fullName, BUFFER_SIZE )
    
    err = Fieldml_GetLastError( fmlHandle )
    IF( INDEX( fullName, 'library.local_nodes.') /= 1 ) THEN
      IF( INDEX( fullName, 'library.' ) /= 1 ) THEN
        name(1:length) = fullName(1:length)
      ELSE
        name(1:length - 7) = fullName(8:length)
        length = length - 7
      ENDIF
    ELSE
      name(1:length - 19) = fullName(20:length)
      length = length - 19
    ENDIF

  END SUBROUTINE

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_GetSimpleBasisName( fmlHandle, basisHandle, name, length, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: basisHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: name
    INTEGER(C_INT) :: length
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: fullName
    
    length = Fieldml_CopyObjectName( fmlHandle, basisHandle, fullName, BUFFER_SIZE )
    
    err = Fieldml_GetLastError( fmlHandle )
    IF( INDEX( fullName, 'library.fem.') /= 1 ) THEN
      IF( INDEX( fullName, 'library.' ) /= 1 ) THEN
        name(1:length) = fullName(1:length)
      ELSE
        name(1:length - 7) = fullName(8:length)
        length = length - 7
      ENDIF
    ELSE
      name(1:length - 11) = fullName(12:length)
      length = length - 11
    ENDIF

  END SUBROUTINE
    
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CreateBasisImport( fmlhandle, meshHandle, nodalDofsHandle, baseName, basisInfo, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlhandle
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(C_INT), INTENT(IN) :: nodalDofsHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(BasisInfoType), INTENT(INOUT) :: basisInfo
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: basisType, xiCount, dofsImportHandle, interpolationParametersHandle, handle, evaluatorHandle
    INTEGER(C_INT), ALLOCATABLE :: xiInterpolations(:)
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length
    TYPE(VARYING_STRING) :: importName
    
    CALL CMISSBasisTypeGet( basisInfo%basisNumber, basisType, err )
    
    CALL CMISSBasisNumberOfXiGet( basisInfo%basisNumber, xiCount, err )
    
    IF( basisType == CMISSBasisLagrangeHermiteTPType ) THEN
      ALLOCATE( xiInterpolations( xiCount ) )
      CALL CMISSBasisInterpolationXiGet( basisInfo%basisNumber, xiInterpolations, err )
      CALL FieldmlUtil_GetTPBasisEvaluator( fmlhandle, xiInterpolations, evaluatorHandle, interpolationParametersHandle, err )
      DEALLOCATE( xiInterpolations )

      CALL FieldmlOutput_GetSimpleBasisName( fmlHandle, evaluatorHandle, name, length, err )
      
      importName = baseName//name(1:length)//".parameters"
      
      handle = Fieldml_GetValueDomain( fmlHandle, basisInfo%connectivityHandle )
      dofsImportHandle = Fieldml_CreateContinuousImport( fmlHandle, char(importName//NUL), nodalDofsHandle, &
        & Fieldml_GetValueDomain( fmlHandle, nodalDofsHandle ) )
      err = Fieldml_SetAlias( fmlHandle, dofsImportHandle, handle, basisInfo%connectivityHandle )
      
      importName = baseName//name(1:length)//".evaluator"

      basisInfo%importHandle = Fieldml_CreateContinuousImport( fmlHandle, char(importName//NUL), evaluatorHandle, &
        & Fieldml_GetValueDomain( fmlHandle, nodalDofsHandle ) )
      CALL FieldmlUtil_GetXiDomain( fmlhandle, xiCount, handle, err )
      err = Fieldml_SetAlias( fmlHandle, basisInfo%importHandle, handle, Fieldml_GetMeshXiDomain( fmlhandle, meshHandle ) )
      err = Fieldml_SetAlias( fmlHandle, basisInfo%importHandle, interpolationParametersHandle, dofsImportHandle )
    ELSE
      basisInfo%importHandle = FML_INVALID_HANDLE
      err = FML_ERR_INVALID_OBJECT
    ENDIF
    
    IF( evaluatorHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_UNKNOWN_OBJECT
    ENDIF

  END SUBROUTINE FieldmlUtil_CreateBasisImport

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CreateLayoutParameters( fmlhandle, elementsHandle, nodesHandle, layoutHandle, componentName, &
    & connectivityInfo, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlhandle
    INTEGER(C_INT), INTENT(IN) :: elementsHandle
    INTEGER(C_INT), INTENT(IN) :: nodesHandle
    INTEGER(C_INT), INTENT(IN) :: layoutHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: componentName
    TYPE(ConnectivityInfoType), INTENT(INOUT) :: connectivityInfo
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length
    TYPE(VARYING_STRING) :: connectivityName

    CALL FieldmlOutput_GetSimpleLayoutName( fmlHandle, layoutHandle, name, length, err )

    connectivityName = componentName//name(1:length)

    connectivityInfo%layoutHandle = layoutHandle
    connectivityInfo%connectivityHandle = Fieldml_CreateEnsembleParameters( fmlHandle, &
      & char(connectivityName//NUL), nodesHandle )
    err = Fieldml_SetParameterDataDescription( fmlHandle, connectivityInfo%connectivityHandle, DESCRIPTION_SEMIDENSE )
    err = Fieldml_AddSemidenseIndex( fmlHandle, connectivityInfo%connectivityHandle, layoutHandle, 0 )
    err = Fieldml_AddSemidenseIndex( fmlHandle, connectivityInfo%connectivityHandle, elementsHandle, 0 )

  END SUBROUTINE FieldmlUtil_CreateLayoutParameters

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_AddConnectivity( fmlHandle, baseName, meshHandle, nodesHandle, nodalDofsHandle, component, &
    & meshElements, templateHandle, err )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: meshHandle
    INTEGER(C_INT), INTENT(IN) :: nodesHandle
    INTEGER(C_INT), INTENT(IN) :: nodalDofsHandle
    INTEGER(INTG), INTENT(IN) :: component
    TYPE(CMISSMeshElementsType), INTENT(IN) :: meshElements
    INTEGER(C_INT), INTENT(OUT) :: templateHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: layoutHandle, connectivityHandle, elementsHandle, elementCount, defaultHandle
    INTEGER(INTG) :: connectivityCount, basisCount, i, j, layoutNodeCount, basisNumber, idx
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(C_INT), ALLOCATABLE, TARGET :: iBuffer(:)
    TYPE(CMISSBasisType) :: basis
    TYPE(C_PTR) :: writer
    TYPE(ConnectivityInfoType), ALLOCATABLE :: connectivityInfo(:), tempConnectivityInfo(:)
    TYPE(BasisInfoType), ALLOCATABLE :: basisInfo(:), tempBasisInfo(:)
    TYPE(VARYING_STRING) :: componentName

    elementsHandle = Fieldml_GetMeshElementDomain( fmlHandle, meshHandle )
    elementCount = Fieldml_GetEnsembleDomainElementCount( fmlHandle, elementsHandle )
    
    connectivityCount = 0
    basisCount = 0
    
    err = FML_ERR_NO_ERROR
    
    componentName = baseName//".component"//TRIM(NUMBER_TO_VSTRING(component,"*",err,errorString))
    
    templateHandle = Fieldml_CreateContinuousPiecewise( fmlHandle, char(componentName//".template"//NUL), &
      & elementsHandle, Fieldml_GetValueDomain( fmlHandle, nodalDofsHandle ) )

    DO i = 1, elementCount
      CALL CMISSMeshElementsBasisGet( meshElements, i, basis, err )
      CALL CMISSUserNumberGet( basis, basisNumber, err )
      CALL FieldmlUtil_GetConnectivityEnsemble( fmlHandle, basisNumber, layoutHandle, err )
      
      idx = FieldmlOutput_FindLayout( connectivityInfo, layoutHandle )
      IF( idx == -1 ) THEN
        IF( connectivityCount == 0 ) THEN
          ALLOCATE( connectivityInfo( connectivityCount + 1 ) )
        ELSE
          ALLOCATE( tempConnectivityInfo( connectivityCount ) )
          tempConnectivityInfo(:) = connectivityInfo(:)
          DEALLOCATE( connectivityInfo )
          ALLOCATE( connectivityInfo( connectivityCount + 1 ) )
          connectivityInfo( 1:connectivityCount ) = tempConnectivityInfo( 1:connectivityCount )
        ENDIF
        
        CALL FieldmlUtil_CreateLayoutParameters( fmlHandle, elementsHandle, nodesHandle, layoutHandle, char(componentName), &
          & connectivityInfo(connectivityCount+1), err )

        err = Fieldml_SetParameterDataLocation( fmlHandle, connectivityInfo(connectivityCount+1)%connectivityHandle, &
          & LOCATION_FILE )
        err = Fieldml_SetParameterFileData( fmlHandle, connectivityInfo(connectivityCount+1)%connectivityHandle, &
          char(componentName//".connectivity"//NUL), TYPE_LINES, connectivityCount * elementCount )
        connectivityCount = connectivityCount + 1
        
        idx = connectivityCount
      ENDIF
      connectivityHandle = connectivityInfo(idx)%connectivityHandle

      idx = FieldmlOutput_FindBasis( basisInfo, basisNumber )
      IF( idx == -1 ) THEN
        IF( basisCount == 0 ) THEN
          ALLOCATE( basisInfo( basisCount + 1 ) )
        ELSE
          ALLOCATE( tempBasisInfo( basisCount ) )
          tempBasisInfo(:) = basisInfo(:)
          DEALLOCATE( basisInfo )
          ALLOCATE( basisInfo( basisCount + 1 ) )
          basisInfo( 1:basisCount ) = tempBasisInfo( 1:basisCount )
        ENDIF

        basisCount = basisCount + 1
        basisInfo( basisCount )%basisNumber = basisNumber
        basisInfo( basisCount )%connectivityHandle = connectivityHandle
        CALL FieldmlUtil_CreateBasisImport( fmlHandle, meshHandle, nodalDofsHandle, char(componentName), &
          & basisInfo(basisCount), err )
        idx = basisCount
      ENDIF

      IF( i == 1 ) THEN
        defaultHandle = basisInfo( idx )%importHandle
        err = Fieldml_SetDefaultEvaluator( fmlHandle, templateHandle, defaultHandle )
      ELSEIF( basisInfo( idx )%importHandle /= defaultHandle ) THEN
        err = Fieldml_SetEvaluator( fmlHandle, templateHandle, i, basisInfo( idx )%importHandle )
      ENDIF
      
    ENDDO

    DO i = 1, connectivityCount
      layoutNodeCount = Fieldml_GetEnsembleDomainElementCount( fmlHandle, connectivityInfo(i)%layoutHandle )
      writer = Fieldml_OpenWriter( fmlHandle, connectivityInfo(i)%connectivityHandle, 0 )
      ALLOCATE( iBuffer( layoutNodeCount ) )
      DO j = 1, elementCount
        CALL CMISSMeshElementsBasisGet( meshElements, j, basis, err )
        CALL CMISSUserNumberGet( basis, basisNumber, err )
        CALL FieldmlUtil_GetConnectivityEnsemble( fmlHandle, basisNumber, layoutHandle, err )
        IF( layoutHandle == connectivityInfo(i)%layoutHandle ) THEN
          CALL CMISSMeshElementsNodesGet( meshElements, j, iBuffer, err )
        ELSE
          iBuffer = 0
        ENDIF
        err = Fieldml_WriteIntSlice( fmlHandle, writer, C_LOC(dummy), C_LOC(iBuffer) )
      ENDDO
      DEALLOCATE( iBuffer )
      err = Fieldml_CloseWriter( fmlHandle, writer )
      err = Fieldml_SetMeshConnectivity( fmlHandle, meshHandle, connectivityInfo(i)%connectivityHandle, &
        & connectivityInfo(i)%layoutHandle )
    ENDDO
    
  END SUBROUTINE FieldmlOutput_AddConnectivity

  !
  !================================================================================================================================
  !

END MODULE FIELDML_OUTPUT_ROUTINES
