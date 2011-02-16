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

!> Output routines for FieldML

MODULE FIELDML_OUTPUT_ROUTINES

  USE KINDS
  USE FIELDML_API
  USE FIELDML_UTIL_ROUTINES
  USE ISO_VARYING_STRING
  USE OPENCMISS
  USE STRINGS
  USE CMISS
  USE BASE_ROUTINES

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: BUFFER_SIZE = 1024
  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  TYPE(VARYING_STRING) :: errorString

  !Interfaces
  TYPE ConnectivityInfoType
    INTEGER(C_INT) :: connectivityHandle
    INTEGER(C_INT) :: layoutHandle
  END TYPE ConnectivityInfoType

  TYPE BasisInfoType
    INTEGER(INTG) :: basisNumber
    INTEGER(C_INT) :: connectivityHandle
    INTEGER(C_INT) :: referenceHandle
  END TYPE BasisInfoType

  INTERFACE FieldmlOutput_AddField
    MODULE PROCEDURE FieldmlOutput_AddField_NoType
    MODULE PROCEDURE FieldmlOutput_AddField_WithType
  END INTERFACE
 
  PUBLIC :: FieldmlOutput_Write, FieldmlOutput_AddField, FieldmlOutput_InitializeInfo, &
    & FieldmlOutput_AddFieldComponents, FieldmlOutput_CreateEnsembleType, FieldmlOutput_CreateContinuousType

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlOutput_GetTPBasisEvaluator( fmlHandle, xiInterpolations, collapseInfo, evaluatorHandle, parametersHandle, &
    & err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: xiInterpolations(:)
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:)
    INTEGER(C_INT), INTENT(OUT) :: evaluatorHandle
    INTEGER(C_INT), INTENT(OUT) :: parametersHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: xiCount, firstInterpolation, i
    TYPE(VARYING_STRING) :: suffix

    CALL ENTERS( "FieldmlOutput_GetTPBasisEvaluator", err, errorString, *999 )
    
    xiCount = SIZE( xiInterpolations )
  
    DO i = 1, xiCount
      IF( i == 1 ) THEN
        firstInterpolation = xiInterpolations(i)
      ELSE IF( xiInterpolations(i) /= firstInterpolation ) THEN
        !Do not yet support inhomogeneous TP bases
        err = FML_ERR_INVALID_OBJECT
      ENDIF
    ENDDO
    CALL FieldmlUtil_CheckError( "Cannot handle inhomogeneous tensor-product basis", err, errorString, *999 )
    
    CALL FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err )

    evaluatorHandle = FML_INVALID_HANDLE
    parametersHandle = FML_INVALID_HANDLE
      
    IF( firstInterpolation == CMISSBasisQuadraticLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        evaluatorHandle = Fieldml_GetObjectByName( fmlHandle, "library.fem.quadratic_lagrange"//NUL )
        parametersHandle = Fieldml_GetObjectByName( fmlHandle, "library.parameters.quadratic_lagrange"//NUL )
      ELSE IF( xiCount == 2 ) THEN
        evaluatorHandle = Fieldml_GetObjectByName( fmlHandle, "library.fem.biquadratic_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetObjectByName( fmlHandle, "library.parameters.biquadratic_lagrange"//char(suffix)//NUL )
      ELSE IF( xiCount == 3 ) THEN
        evaluatorHandle = Fieldml_GetObjectByName( fmlHandle, "library.fem.triquadratic_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetObjectByName( fmlHandle, "library.parameters.triquadratic_lagrange"//char(suffix)//NUL )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
      ENDIF
    ELSE IF( firstInterpolation == CMISSBasisLinearLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        evaluatorHandle = Fieldml_GetObjectByName( fmlHandle, "library.fem.linear_lagrange"//NUL )
        parametersHandle = Fieldml_GetObjectByName( fmlHandle, "library.parameters.linear_lagrange"//NUL )
      ELSE IF( xiCount == 2 ) THEN
        evaluatorHandle = Fieldml_GetObjectByName( fmlHandle, "library.fem.bilinear_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetObjectByName( fmlHandle, "library.parameters.bilinear_lagrange"//char(suffix)//NUL )
      ELSE IF( xiCount == 3 ) THEN
        evaluatorHandle = Fieldml_GetObjectByName( fmlHandle, "library.fem.trilinear_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetObjectByName( fmlHandle, "library.parameters.trilinear_lagrange"//char(suffix)//NUL )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
      ENDIF
    ELSE
      err = FML_ERR_INVALID_OBJECT
    ENDIF
    
    IF( ( evaluatorHandle == FML_INVALID_HANDLE ) .OR. ( parametersHandle == FML_INVALID_HANDLE ) ) THEN
      err = FML_ERR_UNKNOWN_OBJECT
    ENDIF

    CALL FieldmlUtil_CheckError( "Cannot find an evaluator for the given basis", err, errorString, *999 )
    
    CALL EXITS( "FieldmlOutput_GetTPBasisEvaluator" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_GetTPBasisEvaluator", err, errorString )
    CALL EXITS( "FieldmlOutput_GetTPBasisEvaluator" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_GetTPBasisEvaluator

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

  SUBROUTINE FieldmlOutput_GetSimpleLayoutName( fmlHandle, layoutHandle, name, length, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: layoutHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: name
    INTEGER(C_INT) :: length
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: fullName

    CALL ENTERS( "FieldmlOutput_GetSimpleLayoutName", err, errorString, *999 )
    
    length = Fieldml_CopyObjectName( fmlHandle, layoutHandle, fullName, BUFFER_SIZE )
    CALL FieldmlUtil_CheckError("Cannot get name of layout ensemble",fmlHandle,errorString,*999 )
    
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

    CALL EXITS( "FieldmlOutput_GetSimpleLayoutName" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_GetSimpleLayoutName", err, errorString )
    CALL EXITS( "FieldmlOutput_GetSimpleLayoutName" )
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_GetSimpleBasisName( fmlHandle, basisHandle, name, length, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
    INTEGER(C_INT), INTENT(IN) :: basisHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: name
    INTEGER(C_INT) :: length
    INTEGER(INTG), INTENT(OUT) :: err
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: fullName
    
    CALL ENTERS( "FieldmlOutput_GetSimpleBasisName", err, errorString, *999 )

    length = Fieldml_CopyObjectName( fmlHandle, basisHandle, fullName, BUFFER_SIZE )
    CALL FieldmlUtil_CheckError("Cannot get name of basis evaluator",fmlHandle,errorString,*999 )
    
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

    CALL EXITS( "FieldmlOutput_GetSimpleBasisName" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_GetSimpleBasisName", err, errorString )
    CALL EXITS( "FieldmlOutput_GetSimpleBasisName" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_GetSimpleBasisName
    
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlOutput_CreateBasisReference( fieldmlInfo, baseName, basisInfo, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(BasisInfoType), INTENT(INOUT) :: basisInfo
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: basisType, xiCount, dofsReferenceHandle, interpolationParametersHandle, handle, evaluatorHandle
    INTEGER(C_INT) :: variableHandle, aggregateHandle
    INTEGER(C_INT), ALLOCATABLE :: xiInterpolations(:), collapseInfo(:)
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length
    TYPE(VARYING_STRING) :: referenceName

    CALL ENTERS( "FieldmlOutput_CreateBasisReference", err, errorString, *999 )
    
    CALL CMISSBasisTypeGet( basisInfo%basisNumber, basisType, err )
    CALL CMISSBasisNumberOfXiGet( basisInfo%basisNumber, xiCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get basis information", err, errorString, *999 )
    
    IF( basisType == CMISSBasisLagrangeHermiteTPType ) THEN
      ALLOCATE( xiInterpolations( xiCount ) )
      ALLOCATE( collapseInfo( xiCount ) )
      CALL CMISSBasisInterpolationXiGet( basisInfo%basisNumber, xiInterpolations, err )
      CALL CMISSBasisCollapsedXiGet( basisInfo%basisNumber, collapseInfo, err )
      CALL FieldmlUtil_CheckError( "Cannot get tensor-product basis information", err, errorString, *999 )
      
      CALL FieldmlOutput_GetTPBasisEvaluator( fieldmlInfo%fmlHandle, xiInterpolations, collapseInfo, evaluatorHandle, &
        & interpolationParametersHandle, err, *999 )
      DEALLOCATE( xiInterpolations )

      CALL FieldmlOutput_GetSimpleBasisName( fieldmlInfo%fmlHandle, evaluatorHandle, name, length, err, *999 )
      
      referenceName = baseName//name(1:length)//"_"//TRIM(NUMBER_TO_VSTRING(basisInfo%basisNumber,"*",err,errorString))// &
        & ".parameters"
      
      CALL FieldmlUtil_CheckError("Cannot get value type for basis connectivity",fieldmlInfo,errorString,*999 )
      dofsReferenceHandle = Fieldml_CreateAggregateEvaluator( fieldmlInfo%fmlHandle, char(referenceName//NUL), &
        & interpolationParametersHandle )
      CALL FieldmlUtil_CheckError( "Cannot create dofs for basis connectivity", fieldmlInfo, errorString, *999 )

      handle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, interpolationParametersHandle )
      aggregateHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, handle )

      err = Fieldml_SetIndexEvaluator( fieldmlInfo%fmlHandle, dofsReferenceHandle, 1, aggregateHandle )
      err = Fieldml_SetDefaultEvaluator( fieldmlInfo%fmlHandle, dofsReferenceHandle, fieldmlInfo%nodeDofsHandle )

      handle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, basisInfo%connectivityHandle )
      variableHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, handle )
      err = Fieldml_SetBind( fieldmlInfo%fmlHandle, dofsReferenceHandle, variableHandle, basisInfo%connectivityHandle )
      CALL FieldmlUtil_CheckError( "Cannot set bind for basis dofs", err, errorString, *999 )
      
      referenceName = baseName//name(1:length)//"_"//TRIM(NUMBER_TO_VSTRING(basisInfo%basisNumber,"*",err,errorString))// &
        & ".evaluator"

      basisInfo%referenceHandle = Fieldml_CreateReferenceEvaluator( fieldmlInfo%fmlHandle, char(referenceName//NUL), &
        & evaluatorHandle )

      CALL FieldmlUtil_GetXiType( fieldmlInfo%fmlHandle, xiCount, handle, err, *999 )
      variableHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, handle )
      err = Fieldml_SetBind( fieldmlInfo%fmlHandle, basisInfo%referenceHandle, variableHandle, fieldmlInfo%xiVariableHandle )

      variableHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, interpolationParametersHandle )
      err = Fieldml_SetBind( fieldmlInfo%fmlHandle, basisInfo%referenceHandle, variableHandle, &
        & dofsReferenceHandle )
    ELSE
      basisInfo%referenceHandle = FML_INVALID_HANDLE
      err = FML_ERR_INVALID_OBJECT
    ENDIF
    
    IF( evaluatorHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_UNKNOWN_OBJECT
    ENDIF

    CALL EXITS( "FieldmlOutput_CreateBasisReference" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateBasisReference", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateBasisReference" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_CreateBasisReference

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlUtil_CreateLayoutParameters( fieldmlInfo, layoutHandle, componentName, &
    & connectivityInfo, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    INTEGER(C_INT), INTENT(IN) :: layoutHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: componentName
    TYPE(ConnectivityInfoType), INTENT(INOUT) :: connectivityInfo
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    CHARACTER(KIND=C_CHAR,LEN=BUFFER_SIZE) :: name
    INTEGER(INTG) :: length
    INTEGER(C_INT) :: indexHandle
    TYPE(VARYING_STRING) :: connectivityName

    CALL ENTERS( "FieldmlUtil_CreateLayoutParameters", err, errorString, *999 )

    CALL FieldmlOutput_GetSimpleLayoutName( fieldmlInfo%fmlHandle, layoutHandle, name, length, err, *999 )
    connectivityName = componentName//name(1:length)

    connectivityInfo%layoutHandle = layoutHandle
    connectivityInfo%connectivityHandle = Fieldml_CreateParametersEvaluator( fieldmlInfo%fmlHandle, &
      & char(connectivityName//NUL), fieldmlInfo%nodesHandle )
    CALL FieldmlUtil_CheckError("Cannot create nodal parameters",fieldmlInfo%fmlHandle,errorString,*999 )

    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, connectivityInfo%connectivityHandle, &
      & DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError("Cannot set nodal parameters description",err,errorString,*999 )

    indexHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, layoutHandle )
    err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, connectivityInfo%connectivityHandle, indexHandle, &
      & FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError("Add layout index to nodal parameters",err,errorString,*999 )

    err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, connectivityInfo%connectivityHandle, &
      & fieldmlInfo%elementsVariableHandle, FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError("Add element index to nodal parameters",err,errorString,*999 )

    CALL EXITS( "FieldmlUtil_CreateLayoutParameters" )
    RETURN
999 CALL ERRORS( "FieldmlUtil_CreateLayoutParameters", err, errorString )
    CALL EXITS( "FieldmlUtil_CreateLayoutParameters" )
    RETURN 1

  END SUBROUTINE FieldmlUtil_CreateLayoutParameters

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_AddMeshComponent( fieldmlInfo, baseName, componentNumber, meshElements, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(INTG), INTENT(IN) :: componentNumber
    TYPE(CMISSMeshElementsType), INTENT(IN) :: meshElements
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: layoutHandle, connectivityHandle, elementCount, defaultHandle, templateHandle, typeHandle
    INTEGER(INTG) :: connectivityCount, basisCount, i, j, layoutNodeCount, basisNumber, idx
    INTEGER(C_INT), ALLOCATABLE, TARGET :: iBuffer(:)
    TYPE(CMISSBasisType) :: basis
    TYPE(C_PTR) :: writer
    TYPE(ConnectivityInfoType), ALLOCATABLE :: connectivityInfo(:), tempConnectivityInfo(:)
    TYPE(BasisInfoType), ALLOCATABLE :: basisInfo(:), tempBasisInfo(:)
    TYPE(VARYING_STRING) :: componentName
    
    CALL ENTERS( "FieldmlOutput_AddMeshComponent", err, errorString, *999 )

    elementCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    CALL FieldmlUtil_CheckError( "Cannot handle inhomogeneous tensor-product basis", fieldmlInfo, errorString, *999 )
    
    connectivityCount = 0
    basisCount = 0
    
    err = FML_ERR_NO_ERROR
    
    componentName = baseName//".component"//TRIM(NUMBER_TO_VSTRING(componentNumber,"*",err,errorString))
    
    typeHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, fieldmlInfo%nodeDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get node dofs type", fieldmlInfo, errorString, *999 )

    templateHandle = Fieldml_CreatePiecewiseEvaluator( fieldmlInfo%fmlHandle, char(componentName//".template"//NUL), &
      &  typeHandle )
    err = Fieldml_SetIndexEvaluator( fieldmlInfo%fmlHandle, templateHandle, 1, fieldmlInfo%elementsVariableHandle )
      
    CALL FieldmlUtil_CheckError( "Cannot create mesh component template", fieldmlInfo, errorString, *999 )

    DO i = 1, elementCount
      CALL CMISSMeshElementsBasisGet( meshElements, i, basis, err )
      CALL CMISSUserNumberGet( basis, basisNumber, err )
      CALL FieldmlUtil_CheckError( "Cannot OpenCMISS basis info", err, errorString, *999 )

      CALL FieldmlUtil_GetConnectivityEnsemble( fieldmlInfo%fmlHandle, basisNumber, layoutHandle, err, *999 )
      
      idx = -1
      IF( connectivityCount > 0 ) THEN
        idx = FieldmlOutput_FindLayout( connectivityInfo, layoutHandle )
      ENDIF

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
        
        CALL FieldmlUtil_CreateLayoutParameters( fieldmlInfo, &
          & layoutHandle, char(componentName), connectivityInfo(connectivityCount+1), err, *999 )

        err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, connectivityInfo(connectivityCount+1)%connectivityHandle, &
          & LOCATION_FILE )
        CALL FieldmlUtil_CheckError( "Cannot set connectivity data location", err, errorString, *999 )
  
        err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, connectivityInfo(connectivityCount+1)%connectivityHandle, &
          char(componentName//".connectivity"//NUL), TYPE_LINES, connectivityCount * elementCount )
        CALL FieldmlUtil_CheckError( "Cannot set connectivity data file", err, errorString, *999 )
  
        connectivityCount = connectivityCount + 1
        
        idx = connectivityCount
      ENDIF
      connectivityHandle = connectivityInfo(idx)%connectivityHandle

      IF( basisCount == 0 ) THEN
        idx = -1
      ELSE
        idx = FieldmlOutput_FindBasis( basisInfo, basisNumber )
      ENDIF
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
        CALL FieldmlOutput_CreateBasisReference( fieldmlInfo, char(componentName), basisInfo(basisCount), err, *999 )
        idx = basisCount
      ENDIF

      IF( i == 1 ) THEN
        defaultHandle = basisInfo( idx )%referenceHandle
        err = Fieldml_SetDefaultEvaluator( fieldmlInfo%fmlHandle, templateHandle, defaultHandle )
      ELSEIF( basisInfo( idx )%referenceHandle /= defaultHandle ) THEN
        err = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, templateHandle, i, basisInfo( idx )%referenceHandle )
      ENDIF
      CALL FieldmlUtil_CheckError( "Cannot set mesh connectivity evaluator", err, errorString, *999 )
      
    ENDDO

    DO i = 1, connectivityCount
      layoutNodeCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, connectivityInfo(i)%layoutHandle )
      CALL FieldmlUtil_CheckError( "Cannot get layout node count", err, errorString, *999 )
      IF( i == 1 ) THEN
        writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, connectivityInfo(i)%connectivityHandle, 0 )
      ELSE
        writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, connectivityInfo(i)%connectivityHandle, 1 )
      ENDIF
      CALL FieldmlUtil_CheckError( "Cannot open connectivity data writer", fieldmlInfo, errorString, *999 )
      
      ALLOCATE( iBuffer( layoutNodeCount ) )
      DO j = 1, elementCount
        CALL CMISSMeshElementsBasisGet( meshElements, j, basis, err )
        CALL CMISSUserNumberGet( basis, basisNumber, err )
        CALL FieldmlUtil_CheckError( "Cannot OpenCMISS basis info", err, errorString, *999 )
  
        CALL FieldmlUtil_GetConnectivityEnsemble( fieldmlInfo%fmlHandle, basisNumber, layoutHandle, err, *999 )
        IF( layoutHandle == connectivityInfo(i)%layoutHandle ) THEN
          CALL CMISSMeshElementsNodesGet( meshElements, j, iBuffer, err )
        ELSE
          iBuffer = 0
        ENDIF
        err = Fieldml_WriteIntValues( fieldmlInfo%fmlHandle, writer, C_LOC(iBuffer), layoutNodeCount )
        IF( err /= layoutNodeCount ) THEN
            CALL FieldmlUtil_CheckError( "Cannot write connectivity data", err, errorString, *999 )
        ENDIF
      ENDDO
      DEALLOCATE( iBuffer )
      err = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close connectivity data writer", err, errorString, *999 )
    ENDDO
    
    IF( ALLOCATED( basisInfo ) ) THEN
      DEALLOCATE( basisInfo )
    ENDIF
    IF( ALLOCATED( connectivityInfo ) ) THEN
      DEALLOCATE( connectivityInfo )
    ENDIF
    
    fieldmlInfo%componentHandles( componentNumber ) = templateHandle

    CALL EXITS( "FieldmlOutput_AddMeshComponent" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddMeshComponent", err, errorString )
    CALL EXITS( "FieldmlOutput_AddMeshComponent" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_AddMeshComponent

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlOutput_AddFieldNodeDofs( fieldmlInfo, baseName, typeHandle, mesh, field, fieldComponentNumbers, &
    & variableType, nodeDofsHandle, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: typeHandle
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(C_INT), INTENT(INOUT) :: nodeDofsHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: typeComponentHandle, real1DHandle, nodeCount, indexHandle
    INTEGER(INTG) :: versionNumber,componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    TYPE(C_PTR) :: writer
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    REAL(C_DOUBLE) :: dValue
    LOGICAL :: nodeExists
    LOGICAL, ALLOCATABLE :: isNodeBased(:)

    CALL ENTERS( "FieldmlOutput_AddFieldNodeDofs", err, errorString, *999 )
    
    CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, real1DHandle, err, *999 )

    componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    typeComponentHandle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
    nodeCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle )
    
    ALLOCATE( meshComponentNumbers( componentCount ) )
    ALLOCATE( isNodeBased( componentCount ) )

    DO i = 1, componentCount
      CALL CMISSFieldComponentMeshComponentGet( field, variableType, fieldComponentNumbers(i), &
        & meshComponentNumbers(i), Err)

      CALL CMISSFieldComponentInterpolationGet( field, variableType, fieldComponentNumbers(i), &
        & interpolationType, err )
      CALL FieldmlUtil_CheckError( "Cannot get mesh component interpolation info", err, errorString, *999 )
        
      isNodeBased( i ) = ( interpolationType == CMISSFieldNodeBasedInterpolation )
    ENDDO

    nodeDofsHandle = Fieldml_CreateParametersEvaluator( fieldmlInfo%fmlHandle, baseName//".dofs.node"//NUL, real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs parameter set", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, nodeDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter description", err, errorString, *999 )
    err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, nodeDofsHandle, LOCATION_FILE )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter location", err, errorString, *999 )
    err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, nodeDofsHandle, baseName//".dofs.node"//NUL, TYPE_LINES, 0 )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter information", err, errorString, *999 )

    IF( typeComponentHandle /= FML_INVALID_HANDLE ) THEN
      indexHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, typeComponentHandle )
      err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, nodeDofsHandle, indexHandle, FML_INVALID_HANDLE )
      CALL FieldmlUtil_CheckError( "Cannot add component index for nodal dofs parameter set", err, errorString, *999 )
    ENDIF
    err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, nodeDofsHandle, fieldmlInfo%nodesVariableHandle, &
      & FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError( "Cannot add layout index for nodal dofs parameter set", err, errorString, *999 )

    ALLOCATE( dBuffer( componentCount ) )
    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, nodeDofsHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open nodal parameter writer", fieldmlInfo, errorString, *999 )
    DO i = 1, nodeCount
      DO j = 1, componentCount
        dValue = 0
        err = 0
        IF( isNodeBased(j) ) THEN
          CALL CMISSMeshNodeExists( mesh, meshComponentNumbers(j), i, nodeExists, err )
          IF( nodeExists ) THEN
            !Default to version 1 of each node derivative (value hardcoded in loop)
            versionNumber = 1
            CALL CMISSFieldParameterSetGetNode( field, variableType, CMISSFieldValuesSetType, & 
              & versionNumber,CMISSNoGlobalDerivative, i, fieldComponentNumbers(j), dValue, err )
            CALL FieldmlUtil_CheckError( "Cannot get nodal dof value", err, errorString, *999 )
          ENDIF
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      err = Fieldml_WriteDoubleValues( fieldmlInfo%fmlHandle, writer, C_LOC(dBuffer), componentCount )
      IF( err /= componentCount ) THEN
        CALL FieldmlUtil_CheckError( "Cannot write nodal parameter values", err, errorString, *999 )
      ENDIF
    ENDDO
    err = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
    CALL FieldmlUtil_CheckError( "Cannot close nodal parameter writer", err, errorString, *999 )
    DEALLOCATE( dBuffer )
    
    DEALLOCATE( meshComponentNumbers )
    DEALLOCATE( isNodeBased )

    CALL EXITS( "FieldmlOutput_AddFieldNodeDofs" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddFieldNodeDofs", err, errorString )
    CALL EXITS( "FieldmlOutput_AddFieldNodeDofs" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_AddFieldNodeDofs
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlOutput_AddFieldElementDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
    & variableType, elementDofsHandle, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: typeHandle
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(C_INT), INTENT(INOUT) :: elementDofsHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: typeComponentHandle, real1DHandle, elementCount, indexHandle
    INTEGER(INTG) :: componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    TYPE(C_PTR) :: writer
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    REAL(C_DOUBLE) :: dValue
    LOGICAL, ALLOCATABLE :: isElementBased(:)

    CALL EXITS( "FieldmlOutput_AddFieldElementDofs" )
    
    CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, real1DHandle, err, *999 )

    componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    typeComponentHandle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
    
    elementCount = Fieldml_GetEnsembleTypeElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    
    ALLOCATE( meshComponentNumbers( componentCount ) )
    ALLOCATE( isElementBased( componentCount ) )

    DO i = 1, componentCount
      CALL CMISSFieldComponentMeshComponentGet( field, variableType, fieldComponentNumbers(i), &
        & meshComponentNumbers(i), Err)

      CALL CMISSFieldComponentInterpolationGet( field, variableType, fieldComponentNumbers(i), &
        & interpolationType, err )
        
      CALL FieldmlUtil_CheckError( "Cannot get mesh component interpolation info", err, errorString, *999 )

      isElementBased( i ) = ( interpolationType == CMISSFieldElementBasedInterpolation )
    ENDDO

    elementDofsHandle = Fieldml_CreateParametersEvaluator( fieldmlInfo%fmlHandle, baseName//".dofs.element"//NUL, real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs parameter set", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, elementDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter description", err, errorString, *999 )
    err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, elementDofsHandle, LOCATION_FILE )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter location", err, errorString, *999 )
    err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, elementDofsHandle, baseName//".dofs.element"//NUL, TYPE_LINES, 0 )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter information", err, errorString, *999 )

    IF( typeComponentHandle /= FML_INVALID_HANDLE ) THEN
      indexHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, typeComponentHandle )
      err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, elementDofsHandle, typeComponentHandle, FML_INVALID_HANDLE )
      CALL FieldmlUtil_CheckError( "Cannot add component index for element dofs parameter set", err, errorString, *999 )
    ENDIF
    err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, elementDofsHandle, fieldmlInfo%elementsVariableHandle, &
      & FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError( "Cannot add element index for element dofs parameter set", err, errorString, *999 )

    ALLOCATE( dBuffer( componentCount ) )
    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, elementDofsHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open element parameter writer", fieldmlInfo, errorString, *999 )
    DO i = 1, elementCount
      DO j = 1, componentCount
        dValue = 0
        IF( isElementBased(j) ) THEN
          CALL CMISSFieldParameterSetGetElement( field, variableType, CMISSFieldValuesSetType, & 
            & i, fieldComponentNumbers(j), dValue, err )
          CALL FieldmlUtil_CheckError( "Cannot get element dof value", err, errorString, *999 )
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      err = Fieldml_WriteDoubleValues( fieldmlInfo%fmlHandle, writer, C_LOC(dBuffer), componentCount )
      IF( err /= componentCount ) THEN
        CALL FieldmlUtil_CheckError( "Cannot write element parameter values", err, errorString, *999 )
      ENDIF
    ENDDO
    err = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
    CALL FieldmlUtil_CheckError( "Cannot close element parameter writer", err, errorString, *999 )
    DEALLOCATE( dBuffer )
    
    DEALLOCATE( meshComponentNumbers )
    DEALLOCATE( isElementBased )

    CALL EXITS( "FieldmlOutput_AddFieldElementDofs" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddFieldElementDofs", err, errorString )
    CALL EXITS( "FieldmlOutput_AddFieldElementDofs" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_AddFieldElementDofs
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlOutput_AddFieldConstantDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
    & variableType, constantDofsHandle, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: typeHandle
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(C_INT), INTENT(INOUT) :: constantDofsHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: doftypeHandle, typeType, componentType, dataType, indexHandle
    INTEGER(INTG) :: componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    TYPE(C_PTR) :: writer
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    INTEGER(C_INT), ALLOCATABLE, TARGET :: iBuffer(:)
    REAL(C_DOUBLE) :: dValue
    INTEGER(C_INT) :: iValue
    LOGICAL :: isReal
    LOGICAL, ALLOCATABLE :: isConstant(:)

    CALL ENTERS( "FieldmlOutput_AddFieldConstantDofs", err, errorString, *999 )
    
    typeType = Fieldml_GetObjectType(  fieldmlInfo%fmlHandle, typeHandle )

    IF( typeType == FHT_ENSEMBLE_TYPE ) THEN
      doftypeHandle = typeHandle
      componentCount = 1
      componentType = FML_INVALID_HANDLE
      isReal = .FALSE.
    ELSE
      CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, doftypeHandle, err, *999 )
      componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
      componentType = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
      isReal = .TRUE.
    ENDIF
    
    ALLOCATE( meshComponentNumbers( componentCount ) )
    ALLOCATE( isConstant( componentCount ) )

    DO i = 1, componentCount
      CALL CMISSFieldComponentMeshComponentGet( field, variableType, fieldComponentNumbers(i), &
        & meshComponentNumbers(i), Err)

      CALL CMISSFieldComponentInterpolationGet( field, variableType, fieldComponentNumbers(i), &
        & interpolationType, err )
        
      CALL FieldmlUtil_CheckError( "Cannot get mesh component interpolation info", err, errorString, *999 )

      isConstant( i ) = ( interpolationType == CMISSFieldConstantInterpolation )
    ENDDO

    constantDofsHandle = Fieldml_CreateParametersEvaluator( fieldmlInfo%fmlHandle, baseName//".dofs.constant"//NUL, &
      & doftypeHandle )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs parameter set", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, constantDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter description", err, errorString, *999 )
    err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, constantDofsHandle, LOCATION_FILE )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter location", err, errorString, *999 )
    err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, constantDofsHandle, baseName//".dofs.constant"//NUL, &
      & TYPE_LINES, 0 )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter information", err, errorString, *999 )


    IF( componentType /= FML_INVALID_HANDLE ) THEN
      indexHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, componentType )
      err = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, constantDofsHandle, indexHandle, FML_INVALID_HANDLE )
      CALL FieldmlUtil_CheckError( "Cannot add component index for constant dofs parameter set", err, errorString, *999 )
    ENDIF

    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, constantDofsHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open element parameter writer", fieldmlInfo, errorString, *999 )

    CALL CMISSFieldDataTypeGet( field, variableType, dataType, err )
    IF(dataType==CMISSFieldIntgType) THEN
      isReal = .false.
    ELSEIF(dataType==CMISSFieldDPType) THEN
      isReal = .true.
    ENDIF
    
    IF( isReal ) THEN
      ALLOCATE( dBuffer( componentCount ) )
      DO j = 1, componentCount
        dValue = 0
        IF( isConstant(j) ) THEN
          CALL CMISSFieldParameterSetGetConstant( field, variableType, CMISSFieldValuesSetType, & 
            & fieldComponentNumbers(j), dValue, err )
          CALL FieldmlUtil_CheckError( "Cannot get constant dof value", err, errorString, *999 )
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      err = Fieldml_WriteDoubleValues( fieldmlInfo%fmlHandle, writer, C_LOC(dBuffer), componentCount )
      IF( err /= componentCount ) THEN
        CALL FieldmlUtil_CheckError( "Cannot write constant parameter values", err, errorString, *999 )
      ENDIF
      err = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close constant parameter writer", err, errorString, *999 )
      DEALLOCATE( dBuffer )
    ELSE
      ALLOCATE( iBuffer( componentCount ) )
      DO j = 1, componentCount
        iValue = 0
        IF( isConstant(j) ) THEN
          CALL CMISSFieldParameterSetGetConstant( field, variableType, CMISSFieldValuesSetType, & 
            & fieldComponentNumbers(j), iValue, err )
          CALL FieldmlUtil_CheckError( "Cannot get constant dof value", err, errorString, *999 )
        ENDIF
        iBuffer( j ) = iValue
      ENDDO
      err = Fieldml_WriteIntValues( fieldmlInfo%fmlHandle, writer, C_LOC(iBuffer), componentCount )
      IF( err /= componentCount ) THEN
        CALL FieldmlUtil_CheckError( "Cannot write constant parameter values", err, errorString, *999 )
      ENDIF
      err = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close constant parameter writer", err, errorString, *999 )
      DEALLOCATE( iBuffer )
    ENDIF
    
    DEALLOCATE( meshComponentNumbers )
    DEALLOCATE( isConstant )

    CALL EXITS( "FieldmlOutput_AddFieldConstantDofs" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddFieldConstantDofs", err, errorString )
    CALL EXITS( "FieldmlOutput_AddFieldConstantDofs" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_AddFieldConstantDofs
  
  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_InitializeInfo( region, mesh, dimensions, location, baseName, fieldmlInfo, err )
    !Argument variables
    TYPE(CMISSRegionType), INTENT(IN) :: region
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    INTEGER(INTG), INTENT(IN) :: dimensions
    CHARACTER(KIND=C_CHAR,LEN=*) :: location
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(FieldmlInfoType), INTENT(OUT) :: fieldmlInfo
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(INTG) :: componentCount, i, nodeCount, elementCount
    INTEGER(C_INT) :: real1DHandle, xiHandle
    TYPE(CMISSMeshElementsType) :: meshElements
    TYPE(CMISSNodesType) :: Nodes

    CALL ENTERS( "FieldmlOutput_InitializeInfo", err, errorString, *999 )
    fieldmlInfo%fmlHandle = Fieldml_Create( location//NUL, baseName//NUL )
    IF( C_ASSOCIATED( fieldmlInfo%fmlHandle ) ) THEN
      CALL FieldmlUtil_CheckError( "Cannot create fieldml handle", err, errorString, *999 )
    ENDIF

    CALL CMISSNodesTypeInitialise( Nodes, err )
    CALL CMISSRegionNodesGet( Region, Nodes, err )
    CALL CMISSNodesNumberOfNodesGet( Nodes, nodeCount, err )
    CALL FieldmlUtil_CheckError( "Region does not have any nodes", err, errorString, *999 )

    fieldmlInfo%nodesHandle = Fieldml_CreateEnsembleType( fieldmlInfo%fmlHandle, baseName//".nodes"//NUL, 0 )
    CALL FieldmlUtil_CheckError( "Cannot create mesh nodes ensemble", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle, nodeCount )
    CALL FieldmlUtil_CheckError( "Cannot set mesh nodes ensemble bounds", err, errorString, *999 )
    
    fieldmlInfo%nodesVariableHandle = Fieldml_CreateAbstractEvaluator( fieldmlInfo%fmlHandle, baseName//".nodes.variable"//NUL, &
      & fieldmlInfo%nodesHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh nodes variable", fieldmlInfo, errorString, *999 )
    
    CALL CMISSMeshNumberOfElementsGet( Mesh, elementCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get mesh element count", fieldmlInfo, errorString, *999 )

    CALL FieldmlUtil_GetXiEnsemble( fieldmlInfo%fmlHandle, dimensions, xiHandle, err, *999 )
    fieldmlInfo%meshHandle = Fieldml_CreateMeshType( fieldmlInfo%fmlHandle, baseName//".mesh"//NUL, xiHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh type", fieldmlInfo, errorString, *999 )
    err = Fieldml_CreateAbstractEvaluator( fieldmlInfo%fmlHandle, baseName//".mesh.variable"//NUL, fieldmlInfo%meshHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh variable", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, elementCount )
    CALL FieldmlUtil_CheckError( "Cannot set mesh type element count", err, errorString, *999 )

    fieldmlInfo%xiHandle = Fieldml_GetMeshXiType( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle )
    fieldmlInfo%elementsHandle = Fieldml_GetMeshElementType( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle )
    
    fieldmlInfo%xiVariableHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, baseName//".mesh.variable.xi"//NUL )
    fieldmlInfo%elementsVariableHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle,baseName//".mesh.variable.element"//NUL )
    
    CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, real1DHandle, err, *999 )
    
    !TODO Some of these may end up being unused. Should use deferred assignment.
    fieldmlInfo%nodeDofsHandle = Fieldml_CreateAbstractEvaluator( fieldmlInfo%fmlHandle, baseName//".dofs.node"//NUL, &
      & real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs variable", fieldmlInfo, errorString, *999 )
!    fieldmlInfo%elementDofsHandle = Fieldml_CreateAbstractEvaluator( fieldmlInfo%fmlHandle, baseName//".dofs.element"//NUL, & 
!      & real1DHandle )
!    CALL FieldmlUtil_CheckError( "Cannot create element dofs variable", fieldmlInfo, errorString, *999 )
!    fieldmlInfo%constantDofsHandle = Fieldml_CreateAbstractEvaluator( fieldmlInfo%fmlHandle, baseName//".dofs.constant"//NUL, & 
!      & real1DHandle )
!    CALL FieldmlUtil_CheckError( "Cannot create constant dofs variable", fieldmlInfo, errorString, *999 )

    CALL CMISSMeshNumberOfComponentsGet( mesh, componentCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get mesh component count", err, errorString, *999 )
    ALLOCATE( fieldmlInfo%componentHandles( componentCount ) )
    DO i = 1, componentCount
      CALL CMISSMeshElementsTypeInitialise( meshElements, err )
      CALL CMISSMeshElementsGet( mesh, i, meshElements, err )
      CALL FieldmlUtil_CheckError( "Cannot get mesh component", err, errorString, *999 )
      CALL FieldmlOutput_AddMeshComponent( fieldmlInfo, baseName, i, meshElements, err, *999 )
    ENDDO
    
    !TODO Proper shape assignment.
    IF( dimensions == 2 ) THEN
      err = Fieldml_SetMeshDefaultShape( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "library.shape.square"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot set 2D mesh type element shape", err, errorString, *999 )
    ELSE
      err = Fieldml_SetMeshDefaultShape( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "library.shape.cube"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot set 3D mesh type element shape", err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlOutput_InitializeInfo" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_InitializeInfo", err, errorString )
    CALL EXITS( "FieldmlOutput_InitializeInfo" )
    CALL CMISS_HANDLE_ERROR( err, errorString )
    
  END SUBROUTINE FieldmlOutput_InitializeInfo

  !
  !================================================================================================================================
  !
  
  SUBROUTINE FieldmlOutput_AddFieldComponents( fieldmlInfo, typeHandle, baseName, mesh, field, fieldComponentNumbers, &
    & variableType, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    INTEGER(C_INT), INTENT(IN) :: typeHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: fieldHandle, componentHandle, nodalDofsHandle, elementDofsHandle, constantDofsHandle, indexHandle
    INTEGER(INTG) :: componentCount, i, meshComponentNumber, interpolationType
    INTEGER(C_INT), ALLOCATABLE, TARGET :: componentEvaluators(:)
  
    CALL ENTERS( "FieldmlOutput_AddFieldComponents", err, errorString, *999 )

    componentHandle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
    componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    ALLOCATE( componentEvaluators( componentCount ) )

    IF( SIZE( fieldComponentNumbers ) /= componentCount ) THEN
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckError( "Fieldml Component count must match value type component count", &
        & fieldmlInfo, errorString, *999 )
    ENDIF

    nodalDofsHandle = FML_INVALID_HANDLE
    elementDofsHandle = FML_INVALID_HANDLE
    constantDofsHandle = FML_INVALID_HANDLE
    !TODO Other types of interpolation not yet supported.
    DO i = 1, componentCount
      CALL CMISSFieldComponentInterpolationGet( field, variableType, fieldComponentNumbers(i), &
        interpolationType, err )
      CALL FieldmlUtil_CheckError( "Cannot get field component interpolation type", err, errorString, *999 )
        
      IF( interpolationType == CMISSFieldNodeBasedInterpolation ) THEN
        IF( nodalDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldNodeDofs( fieldmlInfo, baseName, typeHandle, mesh, field, fieldComponentNumbers, &
          & variableType, nodalDofsHandle, err, *999 )
        ENDIF
        CALL CMISSFieldComponentMeshComponentGet( field, variableType, fieldComponentNumbers(i), &
          & meshComponentNumber, err )
        CALL FieldmlUtil_CheckError( "Cannot get mesh component for field component", err, errorString, *999 )
        componentEvaluators( i ) = fieldmlInfo%componentHandles(meshComponentNumber)
      ELSEIF( interpolationType == CMISSFieldElementBasedInterpolation ) THEN
        IF( elementDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldElementDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
            & variableType, elementDofsHandle, err, *999 )
        ENDIF
        componentEvaluators( i ) = elementDofsHandle
      ELSEIF( interpolationType == CMISSFieldConstantInterpolation ) THEN
        IF( constantDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldConstantDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
            & variableType, constantDofsHandle, err, *999 )
        ENDIF
        componentEvaluators( i ) = constantDofsHandle
      ENDIF
    ENDDO
    

    IF( componentHandle /= FML_INVALID_HANDLE ) THEN
      fieldHandle = Fieldml_CreateAggregateEvaluator( fieldmlInfo%fmlHandle, baseName//NUL, typeHandle )
      CALL FieldmlUtil_CheckError( "Cannot create aggregate evaluator for field", fieldmlInfo, errorString, *999 )
      indexHandle = FieldmlUtil_GetTypeVariableHandle( fieldmlInfo, componentHandle )
      err = Fieldml_SetIndexEvaluator( fieldmlInfo%fmlHandle, fieldHandle, 1, indexHandle )

      DO i = 1, componentCount
        err = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, fieldHandle, i, componentEvaluators( i ) )
        CALL FieldmlUtil_CheckError( "Cannot set nodal evaluator for aggregate field component", err, errorString, *999 )
      ENDDO
    ELSE
      fieldHandle = Fieldml_CreateReferenceEvaluator( fieldmlInfo%fmlHandle, baseName//NUL, componentEvaluators( 1 ) )
      CALL FieldmlUtil_CheckError( "Cannot create aggregate evaluator for field", fieldmlInfo, errorString, *999 )
    ENDIF

    IF( nodalDofsHandle /= FML_INVALID_HANDLE ) THEN
      err = Fieldml_SetBind( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%nodeDofsHandle, nodalDofsHandle )
      CALL FieldmlUtil_CheckError( "Cannot set nodal dofs bind for field with interpolated elements", err, errorString, *999 )
    ENDIF
!    IF( elementDofsHandle /= FML_INVALID_HANDLE ) THEN
!      err = Fieldml_SetBind( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%elementDofsHandle, elementDofsHandle )
!      CALL FieldmlUtil_CheckError( "Cannot set element dofs bind for field with interpolated elements", err, errorString, *999 )
!    ENDIF
!    IF( constantDofsHandle /= FML_INVALID_HANDLE ) THEN
!      err = Fieldml_SetBind( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%constantDofsHandle, constantDofsHandle )
!      CALL FieldmlUtil_CheckError( "Cannot set constant dofs bind for field with interpolated elements", err, errorString, *999 )
!    ENDIF


    DEALLOCATE( componentEvaluators )
    CALL EXITS( "FieldmlOutput_AddFieldComponents" )
    RETURN
999 DEALLOCATE( componentEvaluators )
    CALL ERRORS( "FieldmlOutput_AddFieldComponents", err, errorString )
    CALL EXITS( "FieldmlOutput_AddFieldComponents" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_AddFieldComponents

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_AddField_NoType( fieldmlInfo, baseName, region, mesh, field, variableType, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(CMISSRegionType), INTENT(IN) :: region
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: typeHandle
    
    CALL ENTERS( "FieldmlOutput_AddField", err, errorString, *999 )

    CALL FieldmlUtil_GetValueType( fieldmlInfo%fmlHandle, region, field, typeHandle, err, *999 )

    CALL FieldmlOutput_AddField_WithType( fieldmlInfo, baseName, mesh, field, variableType, typeHandle, err )

    CALL EXITS( "FieldmlOutput_AddField_NoType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddField_NoType", err, errorString )
    CALL EXITS( "FieldmlOutput_AddField_NoType" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_AddField_NoType

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_AddField_WithType( fieldmlInfo, baseName, mesh, field, variableType, typeHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(INTG), INTENT(IN) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(INTG) :: i, componentCount
    INTEGER(INTG), ALLOCATABLE :: fieldComponentNumbers(:)

    IF( typeHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get value type for field", err, errorString, *999 )
    ENDIF
    
    CALL CMISSFieldNumberOfComponentsGet( field, CMISSFieldUVariableType, componentCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get component count for field", err, errorString, *999 )
    
    ALLOCATE( fieldComponentNumbers( componentCount ) )
    DO i = 1, componentCount
      fieldComponentNumbers(i) = i
    ENDDO

    CALL FieldmlOutput_AddFieldComponents( fieldmlInfo, typeHandle, baseName, mesh, field, fieldComponentNumbers, &
      & variableType, err )
    
    DEALLOCATE( fieldComponentNumbers )

    CALL EXITS( "FieldmlOutput_AddField_WithType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddField_WithType", err, errorString )
    CALL EXITS( "FieldmlOutput_AddField_WithType" )
    CALL CMISS_HANDLE_ERROR( err, errorString )
    
  END SUBROUTINE FieldmlOutput_AddField_WithType

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_CreateEnsembleType( fieldmlInfo, typeName, elementCount, typeHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: typeName
    INTEGER(INTG), INTENT(IN) :: elementCount
    INTEGER(INTG), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlOutput_CreateEnsembleType", err, errorString, *999 )

    typeHandle = Fieldml_CreateEnsembleType( fieldmlInfo%fmlHandle, typeName//NUL, 0 )
    CALL FieldmlUtil_CheckError( "Error creating ensemble type", err, errorString, *999 )
    
    err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, typeHandle, elementCount )
    CALL FieldmlUtil_CheckError( "Error setting ensemble type bounds", err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_CreateEnsembleType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateEnsembleType", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateEnsembleType" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_CreateEnsembleType

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_CreateContinuousType( fieldmlInfo, typeName, componentCount, typeHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: typeName
    INTEGER(INTG), INTENT(IN) :: componentCount
    INTEGER(INTG), INTENT(OUT) :: typeHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Local variables
    INTEGER(INTG) :: componentHandle

    CALL ENTERS( "FieldmlOutput_CreateContinuousType", err, errorString, *999 )

    componentHandle = FML_INVALID_HANDLE

    IF( componentCount > 1 ) THEN
      componentHandle = Fieldml_CreateEnsembleType( fieldmlInfo%fmlHandle, typeName//".components"//NUL, 1 )
      CALL FieldmlUtil_CheckError( "Error creating component type", err, errorString, *999 )
  
      err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, componentHandle, componentCount )
      CALL FieldmlUtil_CheckError( "Error setting component type bounds", err, errorString, *999 )
    ENDIF

    typeHandle = Fieldml_CreateContinuousType( fieldmlInfo%fmlHandle, typeName//NUL, componentHandle )
    CALL FieldmlUtil_CheckError( "Error creating continuous type", err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_CreateContinuousType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateContinuousType", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateContinuousType" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_CreateContinuousType

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_Write( fieldmlInfo, filename, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: filename
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlOutput_Write", err, errorString, *999 )

    err = Fieldml_WriteFile( fieldmlInfo%fmlHandle, filename//NUL )
    CALL FieldmlUtil_CheckError( "Error writing fieldml file", err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_Write" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_Write", err, errorString )
    CALL EXITS( "FieldmlOutput_Write" )
    CALL CMISS_HANDLE_ERROR( err, errorString )
  
  END SUBROUTINE

  !
  !================================================================================================================================
  !

END MODULE FIELDML_OUTPUT_ROUTINES
