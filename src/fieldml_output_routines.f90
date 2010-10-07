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
    MODULE PROCEDURE FieldmlOutput_AddField_NoDomain
    MODULE PROCEDURE FieldmlOutput_AddField_WithDomain
  END INTERFACE
 
  PUBLIC :: FieldmlOutput_Write, FieldmlOutput_AddField, FieldmlOutput_InitializeInfo, &
    & FieldmlOutput_AddFieldComponents, FieldmlOutput_CreateEnsembleDomain, FieldmlOutput_CreateContinuousDomain

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
        evaluatorHandle = Fieldml_GetNamedObject( fmlHandle, "library.fem.quadratic_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlHandle, "library.parameters.quadratic_lagrange"//NUL )
      ELSE IF( xiCount == 2 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlHandle, "library.fem.biquadratic_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlHandle, "library.parameters.biquadratic_lagrange"//char(suffix)//NUL )
      ELSE IF( xiCount == 3 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlHandle, "library.fem.triquadratic_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlHandle, "library.parameters.triquadratic_lagrange"//char(suffix)//NUL )
      ELSE
        !Do not yet support dimensions higher than 3.
        err = FML_ERR_INVALID_OBJECT
      ENDIF
    ELSE IF( firstInterpolation == CMISSBasisLinearLagrangeInterpolation ) THEN
      IF( xiCount == 1 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlHandle, "library.fem.linear_lagrange"//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlHandle, "library.parameters.linear_lagrange"//NUL )
      ELSE IF( xiCount == 2 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlHandle, "library.fem.bilinear_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlHandle, "library.parameters.bilinear_lagrange"//char(suffix)//NUL )
      ELSE IF( xiCount == 3 ) THEN
        evaluatorHandle = Fieldml_GetNamedObject( fmlHandle, "library.fem.trilinear_lagrange"//char(suffix)//NUL )
        parametersHandle = Fieldml_GetNamedObject( fmlHandle, "library.parameters.trilinear_lagrange"//char(suffix)//NUL )
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
      
      handle = Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, basisInfo%connectivityHandle )
      CALL FieldmlUtil_CheckError("Cannot get value domain for basis connectivity",fieldmlInfo,errorString,*999 )
      dofsReferenceHandle = Fieldml_CreateContinuousReference( fieldmlInfo%fmlHandle, char(referenceName//NUL), &
        &fieldmlInfo%nodeDofsHandle, Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, fieldmlInfo%nodeDofsHandle ) )
      CALL FieldmlUtil_CheckError( "Cannot create dofs for basis connectivity", fieldmlInfo, errorString, *999 )
      err = Fieldml_SetAlias( fieldmlInfo%fmlHandle, dofsReferenceHandle, handle, basisInfo%connectivityHandle )
      CALL FieldmlUtil_CheckError( "Cannot set alias for basis dofs", err, errorString, *999 )
      
      referenceName = baseName//name(1:length)//"_"//TRIM(NUMBER_TO_VSTRING(basisInfo%basisNumber,"*",err,errorString))// &
        & ".evaluator"

      basisInfo%referenceHandle = Fieldml_CreateContinuousReference( fieldmlInfo%fmlHandle, char(referenceName//NUL), &
        & evaluatorHandle, Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, fieldmlInfo%nodeDofsHandle ) )
      CALL FieldmlUtil_GetXiDomain( fieldmlInfo%fmlHandle, xiCount, handle, err, *999 )
      err = Fieldml_SetAlias( fieldmlInfo%fmlHandle, basisInfo%referenceHandle, handle, fieldmlInfo%xihandle )
      err = Fieldml_SetAlias( fieldmlInfo%fmlHandle, basisInfo%referenceHandle, interpolationParametersHandle, &
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
  
  SUBROUTINE FieldmlUtil_CreateLayoutParameters( fmlHandle, elementsHandle, nodesHandle, layoutHandle, componentName, &
    & connectivityInfo, err, * )
    !Argument variables
    TYPE(C_PTR), INTENT(IN) :: fmlHandle
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

    CALL ENTERS( "FieldmlUtil_CreateLayoutParameters", err, errorString, *999 )

    CALL FieldmlOutput_GetSimpleLayoutName( fmlHandle, layoutHandle, name, length, err, *999 )

    connectivityName = componentName//name(1:length)

    connectivityInfo%layoutHandle = layoutHandle
    connectivityInfo%connectivityHandle = Fieldml_CreateEnsembleParameters( fmlHandle, &
      & char(connectivityName//NUL), nodesHandle )
    CALL FieldmlUtil_CheckError("Cannot create nodal parameters",fmlHandle,errorString,*999 )

    err = Fieldml_SetParameterDataDescription( fmlHandle, connectivityInfo%connectivityHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError("Cannot set nodal parameters description",err,errorString,*999 )

    err = Fieldml_AddSemidenseIndex( fmlHandle, connectivityInfo%connectivityHandle, layoutHandle, 0 )
    CALL FieldmlUtil_CheckError("Add layout index to nodal parameters",err,errorString,*999 )

    err = Fieldml_AddSemidenseIndex( fmlHandle, connectivityInfo%connectivityHandle, elementsHandle, 0 )
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
    INTEGER(C_INT) :: layoutHandle, connectivityHandle, elementCount, defaultHandle, templateHandle
    INTEGER(INTG) :: connectivityCount, basisCount, i, j, layoutNodeCount, basisNumber, idx
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(C_INT), ALLOCATABLE, TARGET :: iBuffer(:)
    TYPE(CMISSBasisType) :: basis
    TYPE(C_PTR) :: writer
    TYPE(ConnectivityInfoType), ALLOCATABLE :: connectivityInfo(:), tempConnectivityInfo(:)
    TYPE(BasisInfoType), ALLOCATABLE :: basisInfo(:), tempBasisInfo(:)
    TYPE(VARYING_STRING) :: componentName
    
    CALL ENTERS( "FieldmlOutput_AddMeshComponent", err, errorString, *999 )

    elementCount = Fieldml_GetEnsembleDomainElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    CALL FieldmlUtil_CheckError( "Cannot handle inhomogeneous tensor-product basis", fieldmlInfo, errorString, *999 )
    
    connectivityCount = 0
    basisCount = 0
    
    err = FML_ERR_NO_ERROR
    
    componentName = baseName//".component"//TRIM(NUMBER_TO_VSTRING(componentNumber,"*",err,errorString))
    
    templateHandle = Fieldml_CreateContinuousPiecewise( fieldmlInfo%fmlHandle, char(componentName//".template"//NUL), &
      & fieldmlInfo%elementsHandle, Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, fieldmlInfo%nodeDofsHandle ) )
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
        
        CALL FieldmlUtil_CreateLayoutParameters( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle, fieldmlInfo%nodesHandle, &
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
      layoutNodeCount = Fieldml_GetEnsembleDomainElementCount( fieldmlInfo%fmlHandle, connectivityInfo(i)%layoutHandle )
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
        err = Fieldml_WriteIntSlice( fieldmlInfo%fmlHandle, writer, C_LOC(dummy), C_LOC(iBuffer) )
        CALL FieldmlUtil_CheckError( "Cannot write connectivity data", err, errorString, *999 )
      ENDDO
      DEALLOCATE( iBuffer )
      err = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close connectivity data writer", err, errorString, *999 )
      err = Fieldml_SetMeshConnectivity( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, connectivityInfo(i)%connectivityHandle, &
        & connectivityInfo(i)%layoutHandle )
      CALL FieldmlUtil_CheckError( "Cannot set connectivity data", err, errorString, *999 )
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
  
  SUBROUTINE FieldmlOutput_AddFieldNodeDofs( fieldmlInfo, baseName, fieldHandle, mesh, field, fieldComponentNumbers, &
    & variableType, nodeDofsHandle, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: fieldHandle
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(C_INT), INTENT(INOUT) :: nodeDofsHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: domainHandle, real1DHandle, nodeCount
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(INTG) :: componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    TYPE(C_PTR) :: writer
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    REAL(C_DOUBLE) :: dValue
    LOGICAL :: nodeExists
    LOGICAL, ALLOCATABLE :: isNodeBased(:)

    CALL ENTERS( "FieldmlOutput_AddFieldNodeDofs", err, errorString, *999 )
    
    CALL FieldmlUtil_GetGenericDomain( fieldmlInfo%fmlHandle, 1, real1DHandle, err, *999 )

    domainHandle = Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, fieldHandle )
    componentCount = Fieldml_GetDomainComponentCount( fieldmlInfo%fmlHandle, domainHandle )
    domainHandle = Fieldml_GetDomainComponentEnsemble( fieldmlInfo%fmlHandle, domainHandle )
    nodeCount = Fieldml_GetEnsembleDomainElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle )
    
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

    nodeDofsHandle = Fieldml_CreateContinuousParameters( fieldmlInfo%fmlHandle, baseName//".dofs.node"//NUL, real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs parameter set", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, nodeDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter description", err, errorString, *999 )
    err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, nodeDofsHandle, LOCATION_FILE )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter location", err, errorString, *999 )
    err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, nodeDofsHandle, baseName//".dofs.node"//NUL, TYPE_LINES, 0 )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter information", err, errorString, *999 )

    IF( domainHandle /= FML_INVALID_HANDLE ) THEN
      err = Fieldml_AddSemidenseIndex( fieldmlInfo%fmlHandle, nodeDofsHandle, domainHandle, 0 )
      CALL FieldmlUtil_CheckError( "Cannot add component index for nodal dofs parameter set", err, errorString, *999 )
    ENDIF
    err = Fieldml_AddSemidenseIndex( fieldmlInfo%fmlHandle, nodeDofsHandle, fieldmlInfo%nodesHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot add layout index for nodal dofs parameter set", err, errorString, *999 )
    err = Fieldml_SetAlias( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%nodeDofsHandle, nodeDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs alias for field with interpolated elements", err, errorString, *999 )

    ALLOCATE( dBuffer( componentCount ) )
    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, nodeDofsHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open nodal parameter writer", fieldmlInfo, errorString, *999 )
    DO i = 1, nodeCount
      DO j = 1, componentCount
        dValue = 0
        IF( isNodeBased(j) ) THEN
          CALL CMISSMeshNodeExists( mesh, meshComponentNumbers(j), i, nodeExists, err )
          IF( nodeExists ) THEN
            CALL CMISSFieldParameterSetGetNode( field, variableType, CMISSFieldValuesSetType, & 
              & CMISSNoGlobalDerivative, i, fieldComponentNumbers(j), dValue, err )
          ENDIF
          CALL FieldmlUtil_CheckError( "Cannot get nodal dof value", err, errorString, *999 )
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      err = Fieldml_WriteDoubleSlice( fieldmlInfo%fmlHandle, writer, C_LOC(dummy), C_LOC(dBuffer) )
      CALL FieldmlUtil_CheckError( "Cannot write nodal parameter values", err, errorString, *999 )
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
  
  SUBROUTINE FieldmlOutput_AddFieldElementDofs( fieldmlInfo, baseName, fieldHandle, field, fieldComponentNumbers, &
    & variableType, elementDofsHandle, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: fieldHandle
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(C_INT), INTENT(INOUT) :: elementDofsHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: domainHandle, real1DHandle, elementCount
    INTEGER(C_INT), TARGET :: dummy(0)
    INTEGER(INTG) :: componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    TYPE(C_PTR) :: writer
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    REAL(C_DOUBLE) :: dValue
    LOGICAL, ALLOCATABLE :: isElementBased(:)

    CALL EXITS( "FieldmlOutput_AddFieldElementDofs" )
    
    CALL FieldmlUtil_GetGenericDomain( fieldmlInfo%fmlHandle, 1, real1DHandle, err, *999 )

    domainHandle = Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, fieldHandle )
    componentCount = Fieldml_GetDomainComponentCount( fieldmlInfo%fmlHandle, domainHandle )
    domainHandle = Fieldml_GetDomainComponentEnsemble( fieldmlInfo%fmlHandle, domainHandle )
    
    elementCount = Fieldml_GetEnsembleDomainElementCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    
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

    elementDofsHandle = Fieldml_CreateContinuousParameters( fieldmlInfo%fmlHandle, baseName//".dofs.element"//NUL, real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs parameter set", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, elementDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter description", err, errorString, *999 )
    err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, elementDofsHandle, LOCATION_FILE )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter location", err, errorString, *999 )
    err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, elementDofsHandle, baseName//".dofs.element"//NUL, TYPE_LINES, 0 )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter information", err, errorString, *999 )

    IF( domainHandle /= FML_INVALID_HANDLE ) THEN
      err = Fieldml_AddSemidenseIndex( fieldmlInfo%fmlHandle, elementDofsHandle, domainHandle, 0 )
      CALL FieldmlUtil_CheckError( "Cannot add component index for element dofs parameter set", err, errorString, *999 )
    ENDIF
    err = Fieldml_AddSemidenseIndex( fieldmlInfo%fmlHandle, elementDofsHandle, fieldmlInfo%elementsHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot add element index for element dofs parameter set", err, errorString, *999 )
!    err = Fieldml_SetAlias( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%elementDofsHandle, elementDofsHandle )
!    CALL FieldmlUtil_CheckError( "Cannot set element dofs alias for field with constant elements", err, errorString, *999 )

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
      err = Fieldml_WriteDoubleSlice( fieldmlInfo%fmlHandle, writer, C_LOC(dummy), C_LOC(dBuffer) )
      CALL FieldmlUtil_CheckError( "Cannot write element parameter values", err, errorString, *999 )
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
  
  SUBROUTINE FieldmlOutput_AddFieldConstantDofs( fieldmlInfo, baseName, fieldHandle, field, fieldComponentNumbers, &
    & variableType, constantDofsHandle, err, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    INTEGER(C_INT), INTENT(IN) :: fieldHandle
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(C_INT), INTENT(INOUT) :: constantDofsHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: domainHandle, dofDomainHandle, domainType, componentDomain, dataType
    INTEGER(C_INT), TARGET :: dummy(0)
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
    
    domainHandle = Fieldml_GetValueDomain( fieldmlInfo%fmlHandle, fieldHandle )
    domainType = Fieldml_GetObjectType(  fieldmlInfo%fmlHandle, domainHandle )

    IF( domainType == FHT_ENSEMBLE_DOMAIN ) THEN
      dofDomainHandle = domainHandle
      componentCount = 1
      componentDomain = FML_INVALID_HANDLE
      isReal = .FALSE.
    ELSE
      CALL FieldmlUtil_GetGenericDomain( fieldmlInfo%fmlHandle, 1, dofDomainHandle, err, *999 )
      componentCount = Fieldml_GetDomainComponentCount( fieldmlInfo%fmlHandle, domainHandle )
      componentDomain = Fieldml_GetDomainComponentEnsemble( fieldmlInfo%fmlHandle, domainHandle )
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

    constantDofsHandle = Fieldml_CreateContinuousParameters( fieldmlInfo%fmlHandle, baseName//".dofs.constant"//NUL, &
      & dofDomainHandle )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs parameter set", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, constantDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter description", err, errorString, *999 )
    err = Fieldml_SetParameterDataLocation( fieldmlInfo%fmlHandle, constantDofsHandle, LOCATION_FILE )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter location", err, errorString, *999 )
    err = Fieldml_SetParameterFileData( fieldmlInfo%fmlHandle, constantDofsHandle, baseName//".dofs.constant"//NUL, &
      & TYPE_LINES, 0 )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter information", err, errorString, *999 )


    IF( componentDomain /= FML_INVALID_HANDLE ) THEN
      err = Fieldml_AddSemidenseIndex( fieldmlInfo%fmlHandle, constantDofsHandle, componentDomain, 0 )
      CALL FieldmlUtil_CheckError( "Cannot add component index for constant dofs parameter set", err, errorString, *999 )
    ENDIF
!    err = Fieldml_SetAlias( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%constantDofsHandle, constantDofsHandle )
!    CALL FieldmlUtil_CheckError( "Cannot set constant dofs alias for field with constant components", err, errorString, *999 )

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
      err = Fieldml_WriteDoubleSlice( fieldmlInfo%fmlHandle, writer, C_LOC(dummy), C_LOC(dBuffer) )
      CALL FieldmlUtil_CheckError( "Cannot write constant parameter values", err, errorString, *999 )
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
      err = Fieldml_WriteIntSlice( fieldmlInfo%fmlHandle, writer, C_LOC(dummy), C_LOC(iBuffer) )
      CALL FieldmlUtil_CheckError( "Cannot write constant parameter values", err, errorString, *999 )
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

    fieldmlInfo%nodesHandle = Fieldml_CreateEnsembleDomain( fieldmlInfo%fmlHandle, baseName//".nodes"//NUL, FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError( "Cannot create mesh nodes ensemble", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle, nodeCount )
    CALL FieldmlUtil_CheckError( "Cannot set mesh nodes ensemble bounds", err, errorString, *999 )
    err = Fieldml_SetMarkup( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle, "geometric"//NUL, "point"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot set mesh nodes ensemble markup", err, errorString, *999 )
    
    CALL CMISSMeshNumberOfElementsGet( Mesh, elementCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get mesh element count", fieldmlInfo, errorString, *999 )

    CALL FieldmlUtil_GetXiEnsemble( fieldmlInfo%fmlHandle, dimensions, xiHandle, err, *999 )
    fieldmlInfo%meshHandle = Fieldml_CreateMeshDomain( fieldmlInfo%fmlHandle, baseName//".mesh"//NUL, xiHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh domain", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, elementCount )
    CALL FieldmlUtil_CheckError( "Cannot set mesh domain element count", err, errorString, *999 )

    fieldmlInfo%xiHandle = Fieldml_GetMeshXiDomain( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle )
    fieldmlInfo%elementsHandle = Fieldml_GetMeshElementDomain( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle )
    
    CALL FieldmlUtil_GetGenericDomain( fieldmlInfo%fmlHandle, 1, real1DHandle, err, *999 )
    
    !TODO Some of these may end up being unused. Should use deferred assignment.
    fieldmlInfo%nodeDofsHandle = Fieldml_CreateContinuousVariable( fieldmlInfo%fmlHandle, baseName//".dofs.node"//NUL, &
      & real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs variable", fieldmlInfo, errorString, *999 )
!    fieldmlInfo%elementDofsHandle = Fieldml_CreateContinuousVariable( fieldmlInfo%fmlHandle, baseName//".dofs.element"//NUL, & 
!      & real1DHandle )
!    CALL FieldmlUtil_CheckError( "Cannot create element dofs variable", fieldmlInfo, errorString, *999 )
!    fieldmlInfo%constantDofsHandle = Fieldml_CreateContinuousVariable( fieldmlInfo%fmlHandle, baseName//".dofs.constant"//NUL, & 
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
      CALL FieldmlUtil_CheckError( "Cannot set 2D mesh domain element shape", err, errorString, *999 )
    ELSE
      err = Fieldml_SetMeshDefaultShape( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "library.shape.cube"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot set 3D mesh domain element shape", err, errorString, *999 )
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
  
  SUBROUTINE FieldmlOutput_AddFieldComponents( fieldmlInfo, domainHandle, baseName, mesh, field, fieldComponentNumbers, &
    & variableType, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    INTEGER(C_INT), INTENT(IN) :: domainHandle
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:)
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: fieldHandle, componentHandle, nodalDofsHandle, elementDofsHandle, constantDofsHandle
    INTEGER(INTG) :: componentCount, i, meshComponentNumber, interpolationType
  
    CALL ENTERS( "FieldmlOutput_AddFieldComponents", err, errorString, *999 )

    componentHandle = Fieldml_GetDomainComponentEnsemble( fieldmlInfo%fmlHandle, domainHandle )
    componentCount = Fieldml_GetDomainComponentCount( fieldmlInfo%fmlHandle, domainHandle )

    IF( SIZE( fieldComponentNumbers ) /= componentCount ) THEN
      err = FML_ERR_INVALID_OBJECT
      CALL FieldmlUtil_CheckError( "Fieldml Component count must match value domain component count", &
        & fieldmlInfo, errorString, *999 )
    ENDIF

    fieldHandle = Fieldml_CreateContinuousAggregate( fieldmlInfo%fmlHandle, baseName//NUL, domainHandle )
    CALL FieldmlUtil_CheckError( "Cannot create aggregate evaluator for field", fieldmlInfo, errorString, *999 )
    err = Fieldml_SetMarkup( fieldmlInfo%fmlHandle, fieldHandle, "field"//NUL, "true"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot set aggregate evaluator markup", fieldmlInfo, errorString, *999 )

    nodalDofsHandle = FML_INVALID_HANDLE
    elementDofsHandle = FML_INVALID_HANDLE
    constantDofsHandle = FML_INVALID_HANDLE
    !TODO Other types or interpolation not yet supported.
    DO i = 1, componentCount
      CALL CMISSFieldComponentInterpolationGet( field, variableType, fieldComponentNumbers(i), &
        interpolationType, err )
      CALL FieldmlUtil_CheckError( "Cannot get field component interpolation type", err, errorString, *999 )
        
      IF( interpolationType == CMISSFieldNodeBasedInterpolation ) THEN
        IF( nodalDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldNodeDofs( fieldmlInfo, baseName, fieldHandle, mesh, field, fieldComponentNumbers, &
          & variableType, nodalDofsHandle, err, *999 )
        ENDIF
        CALL CMISSFieldComponentMeshComponentGet( field, variableType, fieldComponentNumbers(i), &
          & meshComponentNumber, err )
        CALL FieldmlUtil_CheckError( "Cannot get mesh component for field component", err, errorString, *999 )
        err = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, fieldHandle, i, fieldmlInfo%componentHandles(meshComponentNumber) )
        CALL FieldmlUtil_CheckError( "Cannot set nodal evaluator for aggregate field component", err, errorString, *999 )
      ELSEIF( interpolationType == CMISSFieldElementBasedInterpolation ) THEN
        IF( elementDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldElementDofs( fieldmlInfo, baseName, fieldHandle, field, fieldComponentNumbers, &
            & variableType, elementDofsHandle, err, *999 )
        ENDIF
        err = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, fieldHandle, i, elementDofsHandle )
        CALL FieldmlUtil_CheckError( "Cannot set element evaluator for aggregate field component", err, errorString, *999 )
      ELSEIF( interpolationType == CMISSFieldConstantInterpolation ) THEN
        IF( constantDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldConstantDofs( fieldmlInfo, baseName, fieldHandle, field, fieldComponentNumbers, &
            & variableType, constantDofsHandle, err, *999 )
        ENDIF
        err = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, fieldHandle, i, constantDofsHandle )
        CALL FieldmlUtil_CheckError( "Cannot set constant evaluator for aggregate field component", err, errorString, *999 )
      ENDIF
    ENDDO
    
    CALL EXITS( "FieldmlOutput_AddFieldComponents" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddFieldComponents", err, errorString )
    CALL EXITS( "FieldmlOutput_AddFieldComponents" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_AddFieldComponents

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_AddField_NoDomain( fieldmlInfo, baseName, region, mesh, field, variableType, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(CMISSRegionType), INTENT(IN) :: region
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(C_INT) :: domainHandle
    
    CALL ENTERS( "FieldmlOutput_AddField", err, errorString, *999 )

    CALL FieldmlUtil_GetValueDomain( fieldmlInfo%fmlHandle, region, field, domainHandle, err, *999 )

    CALL FieldmlOutput_AddField_WithDomain( fieldmlInfo, baseName, mesh, field, variableType, domainHandle, err )

    CALL EXITS( "FieldmlOutput_AddField_NoDomain" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddField_NoDomain", err, errorString )
    CALL EXITS( "FieldmlOutput_AddField_NoDomain" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_AddField_NoDomain

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_AddField_WithDomain( fieldmlInfo, baseName, mesh, field, variableType, domainHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: baseName
    TYPE(CMISSMeshType), INTENT(IN) :: mesh
    TYPE(CMISSFieldType), INTENT(IN) :: field
    INTEGER(INTG), INTENT(IN) :: variableType
    INTEGER(INTG), INTENT(IN) :: domainHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Locals
    INTEGER(INTG) :: i, componentCount
    INTEGER(INTG), ALLOCATABLE :: fieldComponentNumbers(:)

    IF( domainHandle == FML_INVALID_HANDLE ) THEN
      err = FML_ERR_UNSUPPORTED
      CALL FieldmlUtil_CheckError( "Cannot get value domain for field", err, errorString, *999 )
    ENDIF
    
    CALL CMISSFieldNumberOfComponentsGet( field, CMISSFieldUVariableType, componentCount, err )
    CALL FieldmlUtil_CheckError( "Cannot get component count for field", err, errorString, *999 )
    
    ALLOCATE( fieldComponentNumbers( componentCount ) )
    DO i = 1, componentCount
      fieldComponentNumbers(i) = i
    ENDDO
    
    CALL FieldmlOutput_AddFieldComponents( fieldmlInfo, domainHandle, baseName, mesh, field, fieldComponentNumbers, &
      & variableType, err )
    
    DEALLOCATE( fieldComponentNumbers )

    CALL EXITS( "FieldmlOutput_AddField_WithDomain" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddField_WithDomain", err, errorString )
    CALL EXITS( "FieldmlOutput_AddField_WithDomain" )
    CALL CMISS_HANDLE_ERROR( err, errorString )
    
  END SUBROUTINE FieldmlOutput_AddField_WithDomain

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_CreateEnsembleDomain( fieldmlInfo, domainName, elementCount, domainHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: domainName
    INTEGER(INTG), INTENT(IN) :: elementCount
    INTEGER(INTG), INTENT(OUT) :: domainHandle
    INTEGER(INTG), INTENT(OUT) :: err

    CALL ENTERS( "FieldmlOutput_CreateEnsembleDomain", err, errorString, *999 )

    domainHandle = Fieldml_CreateEnsembleDomain( fieldmlInfo%fmlHandle, domainName//NUL, FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError( "Error creating ensemble domain", err, errorString, *999 )
    
    err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, domainHandle, elementCount )
    CALL FieldmlUtil_CheckError( "Error setting ensemble domain bounds", err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_CreateEnsembleDomain" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateEnsembleDomain", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateEnsembleDomain" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_CreateEnsembleDomain

  !
  !================================================================================================================================
  !

  SUBROUTINE FieldmlOutput_CreateContinuousDomain( fieldmlInfo, domainName, componentCount, domainHandle, err )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo
    CHARACTER(KIND=C_CHAR,LEN=*) :: domainName
    INTEGER(INTG), INTENT(IN) :: componentCount
    INTEGER(INTG), INTENT(OUT) :: domainHandle
    INTEGER(INTG), INTENT(OUT) :: err

    !Local variables
    INTEGER(INTG) :: componentHandle

    CALL ENTERS( "FieldmlOutput_CreateContinuousDomain", err, errorString, *999 )

    componentHandle = FML_INVALID_HANDLE

    IF( componentCount > 1 ) THEN
      componentHandle = Fieldml_CreateComponentEnsembleDomain( fieldmlInfo%fmlHandle, domainName//".components"//NUL )
      CALL FieldmlUtil_CheckError( "Error creating component domain", err, errorString, *999 )
  
      err = Fieldml_SetContiguousBoundsCount( fieldmlInfo%fmlHandle, componentHandle, componentCount )
      CALL FieldmlUtil_CheckError( "Error setting component domain bounds", err, errorString, *999 )
    ENDIF

    domainHandle = Fieldml_CreateContinuousDomain( fieldmlInfo%fmlHandle, domainName//NUL, componentHandle )
    CALL FieldmlUtil_CheckError( "Error creating continuous domain", err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_CreateContinuousDomain" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateContinuousDomain", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateContinuousDomain" )
    CALL CMISS_HANDLE_ERROR( err, errorString )

  END SUBROUTINE FieldmlOutput_CreateContinuousDomain

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
