!> \file
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

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE FIELD_ROUTINES
  USE FIELDML_API
  USE FIELDML_TYPES
  USE FIELDML_UTIL_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE REGION_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  CHARACTER(C_CHAR), PARAMETER :: NUL=C_NULL_CHAR

  !Interfaces
  TYPE ConnectivityInfoType
    INTEGER(C_INT) :: connectivityHandle !<The basis connectivity evaluator handle.
    INTEGER(C_INT) :: layoutHandle !<The local node layout.
  END TYPE ConnectivityInfoType

  TYPE BasisInfoType
    TYPE(BASIS_TYPE), POINTER :: basis
    INTEGER(C_INT) :: connectivityHandle !<The basis connectivity evaluator handle.
    INTEGER(C_INT) :: referenceHandle !<The reference evaluator representing the basis.
    INTEGER(C_INT) :: layoutHandle !<The local node layout.
  END TYPE BasisInfoType

  INTERFACE FieldmlOutput_AddField
    MODULE PROCEDURE FieldmlOutput_AddField_NoType
    MODULE PROCEDURE FieldmlOutput_AddField_WithType
  END INTERFACE
 
  PUBLIC :: FieldmlOutput_Write, FieldmlOutput_AddField, FieldmlOutput_InitialiseInfo, &
    & FieldmlOutput_AddFieldComponents, FieldmlOutput_CreateEnsembleType, FieldmlOutput_CreateContinuousType

CONTAINS

  !
  !================================================================================================================================
  !
  
  !<Get an evaluator from the built-in library that corresponds to the given OpenCMISS tensor-product basis.
  SUBROUTINE FieldmlOutput_GetTPBasisEvaluator( fmlHandle, xiInterpolations, collapseInfo, evaluatorHandle, parametersHandle, &
    & err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle
    INTEGER(C_INT), INTENT(IN) :: xiInterpolations(:) !<The per-xi interpolations used by the basis.
    INTEGER(C_INT), INTENT(IN) :: collapseInfo(:) !<The basis collapse info.
    INTEGER(C_INT), INTENT(OUT) :: evaluatorHandle !<The evaluator handle for the basis.
    INTEGER(C_INT), INTENT(OUT) :: parametersHandle !<The basis parameters argument evaluator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: xiCount, firstInterpolation, i, importIndex
    TYPE(VARYING_STRING) :: suffix, interpolatorName, parameterName

    CALL ENTERS( "FieldmlOutput_GetTPBasisEvaluator", err, errorString, *999 )
    
    xiCount = SIZE( xiInterpolations )
  
    DO i = 1, xiCount
      IF( i == 1 ) THEN
        firstInterpolation = xiInterpolations(i)
      ELSE IF( xiInterpolations(i) /= firstInterpolation ) THEN
        !Do not yet support inhomogeneous TP bases
        CALL FLAG_ERROR( "Cannot yet handle inhomogeneous tensor-product basis", err, errorString, *999 )
      ENDIF
    ENDDO
   
    CALL FieldmlUtil_GetCollapseSuffix( collapseInfo, suffix, err, errorString )

    evaluatorHandle = FML_INVALID_HANDLE
    parametersHandle = FML_INVALID_HANDLE
      
    IF( firstInterpolation == BASIS_QUADRATIC_LAGRANGE_INTERPOLATION ) THEN
      IF( xiCount == 1 ) THEN
        interpolatorName = "interpolator.1d.unit.quadraticLagrange"
        parameterName = "parameters.1d.unit.quadraticLagrange"
      ELSE IF( xiCount == 2 ) THEN
        interpolatorName = "interpolator.2d.unit.biquadraticLagrange"//suffix
        parameterName = "parameters.2d.unit.biquadraticLagrange"//suffix
      ELSE IF( xiCount == 3 ) THEN
        interpolatorName = "interpolator.3d.unit.triquadraticLagrange"//suffix
        parameterName = "parameters.3d.unit.triquadraticLagrange"//suffix
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( "Cannot find an evaluator for the given basis", err, errorString, *999 )
      ENDIF
    ELSE IF( firstInterpolation == BASIS_LINEAR_LAGRANGE_INTERPOLATION ) THEN
      IF( xiCount == 1 ) THEN
        interpolatorName = "interpolator.1d.unit.linearLagrange"
        parameterName = "parameters.1d.unit.linearLagrange"
      ELSE IF( xiCount == 2 ) THEN
        interpolatorName = "interpolator.2d.unit.bilinearLagrange"//suffix
        parameterName = "parameters.2d.unit.bilinearLagrange"//suffix
      ELSE IF( xiCount == 3 ) THEN
        interpolatorName = "interpolator.3d.unit.trilinearLagrange"//suffix
        parameterName = "parameters.3d.unit.trilinearLagrange"//suffix
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( "Cannot find an evaluator for the given basis", err, errorString, *999 )
      ENDIF
    ELSE
      CALL FLAG_ERROR( "Cannot find an evaluator for the given basis", err, errorString, *999 )
    ENDIF

    importIndex = Fieldml_AddImportSource( fmlHandle, &
      & "http://www.fieldml.org/resources/xml/0.4/FieldML_Library_0.4.xml"//NUL, "library"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot access library", fmlHandle, err, errorString, *999 )

    evaluatorHandle = Fieldml_AddImport( fmlHandle, importIndex, cchar(interpolatorName), cchar(interpolatorName) )
    parametersHandle = Fieldml_AddImport( fmlHandle, importIndex, cchar(parameterName), cchar(parameterName) )
    
    IF( ( evaluatorHandle == FML_INVALID_HANDLE ) .OR. ( parametersHandle == FML_INVALID_HANDLE ) ) THEN
      CALL FLAG_ERROR( "Cannot find an evaluator for the given basis", err, errorString, *999 )
    ENDIF
    
    CALL EXITS( "FieldmlOutput_GetTPBasisEvaluator" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_GetTPBasisEvaluator", err, errorString )
    CALL EXITS( "FieldmlOutput_GetTPBasisEvaluator" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_GetTPBasisEvaluator

  !
  !================================================================================================================================
  !
  
  !<Returns the index of the layout handle used by the given connectivity info array, or -1 if none can be found.
  FUNCTION FieldmlOutput_FindLayout( connectivityInfo, layoutHandle )
    !Argument variables
    TYPE(ConnectivityInfoType), INTENT(IN) :: connectivityInfo(:) !<The connectivity info array to search.
    INTEGER(C_INT), INTENT(IN) :: layoutHandle !<The local node layout handle to search for.
    
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
  
  !<Returns the index of the basis handle used by the given basis info array, or -1 if none can be found.
  FUNCTION FieldmlOutput_FindBasis( basisInfo, basis )
    !Argument variables
    TYPE(BasisInfoType), INTENT(IN) :: basisInfo(:) !<The basis info array to search.
    TYPE(BASIS_TYPE), POINTER, INTENT(IN) :: basis !<The basis handle to search for. 
    
    !Function
    INTEGER(INTG) :: FieldmlOutput_FindBasis
    
    !Locals
    INTEGER(INTG) :: i
    
    FieldmlOutput_FindBasis = -1
    DO i = 1, SIZE( basisInfo )
      IF( ASSOCIATED( basisInfo(i)%basis, TARGET = basis ) ) THEN
        FieldmlOutput_FindBasis = i
      ENDIF
    ENDDO
  
  END FUNCTION FieldmlOutput_FindBasis
  
  !
  !================================================================================================================================
  !

  !<Returns the simplified name of the given layout. This is used for naming associated connectivity evaluators.
  SUBROUTINE FieldmlOutput_GetSimpleLayoutName( fmlHandle, layoutHandle, name, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle
    INTEGER(C_INT), INTENT(IN) :: layoutHandle !<The local node layout.
    TYPE(VARYING_STRING), INTENT(INOUT) :: name !<The simplified name.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: fullName
    INTEGER(INTG) :: length

    CALL ENTERS( "FieldmlOutput_GetSimpleLayoutName", err, errorString, *999 )
    
    length = Fieldml_CopyObjectDeclaredName( fmlHandle, layoutHandle, fullName, MAXSTRLEN )
    CALL FieldmlUtil_CheckError("Cannot get name of layout ensemble", fmlHandle, err, errorString, *999 )
    
    IF( INDEX( fullName, 'localNodes.') /= 1 ) THEN
      name = fullName(1:length)
    ELSE
      name = fullName(12:length)
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

  !<Returns the simplified name of the given basis. This is used for naming associated reference evaluators.
  SUBROUTINE FieldmlOutput_GetSimpleBasisName( fmlHandle, basisHandle, name, err, errorString, * )
    !Argument variables
    INTEGER(C_INT), INTENT(IN) :: fmlHandle !<The FieldML session handle
    INTEGER(C_INT), INTENT(IN) :: basisHandle !<The basis handle.
    TYPE(VARYING_STRING), INTENT(INOUT) :: name !<The simplified name.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: fullName
    INTEGER(INTG) :: length
    
    CALL ENTERS( "FieldmlOutput_GetSimpleBasisName", err, errorString, *999 )

    length = Fieldml_CopyObjectDeclaredName( fmlHandle, basisHandle, fullName, MAXSTRLEN )
    CALL FieldmlUtil_CheckError("Cannot get name of basis evaluator", fmlHandle, err, errorString, *999 )
    
    IF( INDEX( fullName, 'interpolator.1d.unit.') == 1 ) THEN
      name = fullName(22:length)
    ELSEIF( INDEX( fullName, 'interpolator.2d.unit.') == 1 ) THEN
      name = fullName(22:length)
    ELSEIF( INDEX( fullName, 'interpolator.3d.unit.') == 1 ) THEN
      name = fullName(22:length)
    ELSE
      name = fullName(1:length)
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
  
  !<Create a basis evaluator from the given basis info.
  SUBROUTINE FieldmlOutput_CreateBasisReference( fieldmlInfo, baseName, basisInfo, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    TYPE(BasisInfoType), INTENT(INOUT) :: basisInfo !<The basis info describing the basis to create.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: basisType, xiCount, interpolationParametersHandle, handle, evaluatorHandle
    INTEGER(C_INT) :: variableHandle, aggregateHandle, indexEvaluatorHandle, fmlErr
    INTEGER(C_INT), ALLOCATABLE :: xiInterpolations(:), collapseInfo(:)
    TYPE(VARYING_STRING) :: referenceName, name

    CALL ENTERS( "FieldmlOutput_CreateBasisReference", err, errorString, *999 )
    
    CALL BASIS_TYPE_GET( basisInfo%basis, basisType, err, errorString, *999 )
    CALL BASIS_NUMBER_OF_XI_GET( basisInfo%basis, xiCount, err, errorString, *999 )
    
    IF( basisType == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      ALLOCATE( xiInterpolations( xiCount ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate xi interpolation array.", err, errorString, *999 )
      ALLOCATE( collapseInfo( xiCount ), STAT = err )
      CALL BASIS_INTERPOLATION_XI_GET( basisInfo%basis, xiInterpolations, err, errorString, *999 )
      CALL BASIS_COLLAPSED_XI_GET( basisInfo%basis, collapseInfo, err, errorString, *999 )
      
      CALL FieldmlOutput_GetTPBasisEvaluator( fieldmlInfo%fmlHandle, xiInterpolations, collapseInfo, evaluatorHandle, &
        & interpolationParametersHandle, err, errorString, *999 )
      DEALLOCATE( xiInterpolations )
      DEALLOCATE( collapseInfo )

      CALL FieldmlOutput_GetSimpleBasisName( fieldmlInfo%fmlHandle, evaluatorHandle, name, err, errorString, *999 )
      
      referenceName = baseName//name//"_"//TRIM(NUMBER_TO_VSTRING(basisInfo%basis%USER_NUMBER,"*",err,errorString))// &
        & ".parameters"
      
      aggregateHandle = Fieldml_CreateAggregateEvaluator( fieldmlInfo%fmlHandle, cchar(referenceName), &
        & interpolationParametersHandle )
      CALL FieldmlUtil_CheckError( "Cannot create dofs for basis connectivity", fieldmlInfo, err, errorString, *999 )

      indexEvaluatorHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, basisInfo%layoutHandle, .TRUE. )

      fmlErr = Fieldml_SetIndexEvaluator( fieldmlInfo%fmlHandle, aggregateHandle, 1, indexEvaluatorHandle )
      CALL FieldmlUtil_CheckError( "Cannot set field component index evaluator.", fieldmlInfo, err, errorString, *999 )
      
      fmlErr = Fieldml_SetDefaultEvaluator( fieldmlInfo%fmlHandle, aggregateHandle, fieldmlInfo%nodeDofsHandle )
      CALL FieldmlUtil_CheckError( "Cannot set nodal field dofs.", fieldmlInfo, err, errorString, *999 )

      handle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, basisInfo%connectivityHandle )
      variableHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, handle, .FALSE. )
      fmlErr = Fieldml_SetBind( fieldmlInfo%fmlHandle, aggregateHandle, variableHandle, basisInfo%connectivityHandle )
      CALL FieldmlUtil_CheckError( "Cannot set bind for basis dofs", fieldmlInfo, err, errorString, *999 )
      
      referenceName = baseName//name//"_"//TRIM(NUMBER_TO_VSTRING(basisInfo%basis%USER_NUMBER,"*",err,errorString))// &
        & ".evaluator"

      basisInfo%referenceHandle = Fieldml_CreateReferenceEvaluator( fieldmlInfo%fmlHandle, cchar(referenceName), &
        & evaluatorHandle )

      CALL FieldmlUtil_GetXiType( fieldmlInfo%fmlHandle, xiCount, .TRUE., handle, err, errorString, *999 )
      variableHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, handle, .TRUE. )
      fmlErr = Fieldml_SetBind( fieldmlInfo%fmlHandle, basisInfo%referenceHandle, variableHandle, fieldmlInfo%xiArgumentHandle )
      CALL FieldmlUtil_CheckError( "Cannot bind xi to basis evaluator.", fieldmlInfo, err, errorString, *999 )

      variableHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, interpolationParametersHandle, .TRUE. )
      fmlErr = Fieldml_SetBind( fieldmlInfo%fmlHandle, basisInfo%referenceHandle, variableHandle, &
        & aggregateHandle )
      CALL FieldmlUtil_CheckError( "Cannot bind parameters to basis evaluator.", fieldmlInfo, err, errorString, *999 )
    ELSE
      basisInfo%referenceHandle = FML_INVALID_HANDLE
      CALL FLAG_ERROR( "FieldML export code can currently only handle tensor-product bases.", err, errorString, *999 )
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
  
  !<Create a parameter evaluator for the given local node layout.
  SUBROUTINE FieldmlOutput_CreateLayoutParameters( fieldmlInfo, layoutHandle, componentName, &
    & connectivityInfo, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: layoutHandle !<The local node layout.
    TYPE(VARYING_STRING), INTENT(IN) :: componentName !<The component name.
    TYPE(ConnectivityInfoType), INTENT(INOUT) :: connectivityInfo !<The connectivity info for the local node layout.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    TYPE(VARYING_STRING) :: name
    INTEGER(C_INT) :: indexHandle, fmlErr
    TYPE(VARYING_STRING) :: connectivityName

    CALL ENTERS( "FieldmlOutput_CreateLayoutParameters", err, errorString, *999 )

    CALL FieldmlOutput_GetSimpleLayoutName( fieldmlInfo%fmlHandle, layoutHandle, name, err, errorString, *999 )
    connectivityName = componentName//name

    connectivityInfo%layoutHandle = layoutHandle
    connectivityInfo%connectivityHandle = Fieldml_CreateParameterEvaluator( fieldmlInfo%fmlHandle, &
      & cchar(connectivityName), fieldmlInfo%nodesHandle )
    CALL FieldmlUtil_CheckError("Cannot create nodal parameters", fieldmlInfo, err, errorString, *999 )

    fmlErr = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, connectivityInfo%connectivityHandle, &
      & DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError("Cannot set nodal parameters description",fieldmlInfo, err, errorString, *999 )

    indexHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, layoutHandle, .TRUE. )
    fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, connectivityInfo%connectivityHandle, indexHandle, &
      & FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError("Cannot add layout index to nodal parameters", fieldmlInfo, err, errorString, *999 )

    fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, connectivityInfo%connectivityHandle, &
      & fieldmlInfo%elementsArgumentHandle, FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError("Cannot add element index to nodal parameters", fieldmlInfo, err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_CreateLayoutParameters" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateLayoutParameters", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateLayoutParameters" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_CreateLayoutParameters

  !
  !================================================================================================================================
  !

  !<Add an evaluator corresponding to the given component of the given OpenCMISS mesh.
  SUBROUTINE FieldmlOutput_AddMeshComponent( fieldmlInfo, baseName, componentNumber, meshElements, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The mesh component number for which an evaluator should be constructed.
    TYPE(MESH_ELEMENTS_TYPE), POINTER, INTENT(IN) :: meshElements !<The mesh element from which to obtain topology info.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: layoutHandle, connectivityHandle, elementCount, defaultHandle, templateHandle, typeHandle, resourceHandle
    INTEGER(INTG) :: connectivityCount, basisCount, i, j, layoutNodeCount, idx
    INTEGER(C_INT), ALLOCATABLE, TARGET :: iBuffer(:)
    TYPE(BASIS_TYPE), POINTER :: basis
    INTEGER(C_INT) :: writer, sourceHandle, writeCount, fmlErr
    TYPE(ConnectivityInfoType), ALLOCATABLE :: connectivityInfo(:), tempConnectivityInfo(:)
    TYPE(BasisInfoType), ALLOCATABLE :: basisInfo(:), tempBasisInfo(:)
    TYPE(VARYING_STRING) :: componentName !<The component name.
    
    CALL ENTERS( "FieldmlOutput_AddMeshComponent", err, errorString, *999 )

    elementCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    CALL FieldmlUtil_CheckError( "Cannot handle inhomogeneous tensor-product basis", fieldmlInfo, err, errorString, *999 )
    
    connectivityCount = 0
    basisCount = 0
    
    componentName = baseName//".component"//TRIM(NUMBER_TO_VSTRING(componentNumber,"*",err,errorString))
    
    typeHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, fieldmlInfo%nodeDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get node dofs type", fieldmlInfo, err, errorString, *999 )

    templateHandle = Fieldml_CreatePiecewiseEvaluator( fieldmlInfo%fmlHandle, cchar(componentName//".template"), &
      &  typeHandle )
    fmlErr = Fieldml_SetIndexEvaluator( fieldmlInfo%fmlHandle, templateHandle, 1, fieldmlInfo%elementsArgumentHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh component template", fieldmlInfo, err, errorString, *999 )
    
    resourceHandle = Fieldml_CreateTextFileDataResource( fieldmlInfo%fmlHandle, &
      & cchar(componentName//".connectivity.resource"), cchar(componentName//".connectivity") )
    CALL FieldmlUtil_CheckError( "Cannot create mesh component connectivity resource", fieldmlInfo, err, errorString, *999 )

    DO i = 1, elementCount
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET( i, meshElements, basis, err, errorString, *999 )

      CALL FieldmlUtil_GetConnectivityEnsemble( fieldmlInfo%fmlHandle, basis, layoutHandle, err, errorString, *999 )
      
      idx = -1
      IF( connectivityCount > 0 ) THEN
        idx = FieldmlOutput_FindLayout( connectivityInfo, layoutHandle )
      ENDIF

      IF( idx == -1 ) THEN
        IF( connectivityCount == 0 ) THEN
          ALLOCATE( connectivityInfo( connectivityCount + 1 ), STAT = err )
          IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity info array.", err, errorString, *999 )
        ELSE
          ALLOCATE( tempConnectivityInfo( connectivityCount ), STAT = err )
          IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate temporary connectivity array.", err, errorString, *999 )
          tempConnectivityInfo(:) = connectivityInfo(:)
          DEALLOCATE( connectivityInfo )
          ALLOCATE( connectivityInfo( connectivityCount + 1 ), STAT = err )
          IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate new connectivity info array.", err, errorString, *999 )
          connectivityInfo( 1:connectivityCount ) = tempConnectivityInfo( 1:connectivityCount )
        ENDIF
        
        CALL FieldmlOutput_CreateLayoutParameters( fieldmlInfo, layoutHandle, componentName, &
          & connectivityInfo(connectivityCount+1), err, errorString, *999 )
          
        layoutNodeCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, connectivityInfo(connectivityCount+1)%layoutHandle )
        sourceHandle = Fieldml_CreateTextDataSource( fieldmlInfo%fmlHandle, cchar(componentName//".connectivity"), &
          & resourceHandle, ( connectivityCount * elementCount ) + 1, elementCount, layoutNodeCount, 0, 0 )
        CALL FieldmlUtil_CheckError( "Cannot create connectivity data source", fieldmlInfo, err, errorString, *999 )

        fmlErr = Fieldml_SetDataSource( fieldmlInfo%fmlHandle, connectivityInfo(connectivityCount+1)%connectivityHandle, &
          & sourceHandle )
        CALL FieldmlUtil_CheckError( "Cannot set connectivity data source", fieldmlInfo, err, errorString, *999 )
  
        connectivityCount = connectivityCount + 1
        
        idx = connectivityCount
      ENDIF
      connectivityHandle = connectivityInfo(idx)%connectivityHandle

      IF( basisCount == 0 ) THEN
        idx = -1
      ELSE
        idx = FieldmlOutput_FindBasis( basisInfo, basis )
      ENDIF
      IF( idx == -1 ) THEN
        IF( basisCount == 0 ) THEN
          ALLOCATE( basisInfo( basisCount + 1 ), STAT = err )
          IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate basis info array.", err, errorString, *999 )
        ELSE
          ALLOCATE( tempBasisInfo( basisCount ), STAT = err )
          IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate temporary basis info array.", err, errorString, *999 )
          tempBasisInfo(:) = basisInfo(:)
          DEALLOCATE( basisInfo )
          ALLOCATE( basisInfo( basisCount + 1 ), STAT = err )
          IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate new basis info array.", err, errorString, *999 )
          basisInfo( 1:basisCount ) = tempBasisInfo( 1:basisCount )
        ENDIF

        basisCount = basisCount + 1
        basisInfo( basisCount )%basis => basis
        basisInfo( basisCount )%connectivityHandle = connectivityHandle
        basisInfo( basisCount )%layoutHandle = layoutHandle
        CALL FieldmlOutput_CreateBasisReference( fieldmlInfo, componentName, basisInfo(basisCount), err, errorString, *999 )
        idx = basisCount
      ENDIF

      IF( i == 1 ) THEN
        defaultHandle = basisInfo( idx )%referenceHandle
        fmlErr = Fieldml_SetDefaultEvaluator( fieldmlInfo%fmlHandle, templateHandle, defaultHandle )
      ELSEIF( basisInfo( idx )%referenceHandle /= defaultHandle ) THEN
        fmlErr = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, templateHandle, i, basisInfo( idx )%referenceHandle )
      ENDIF
      CALL FieldmlUtil_CheckError( "Cannot set mesh connectivity evaluator", fieldmlInfo, err, errorString, *999 )
      
    ENDDO

    DO i = 1, connectivityCount
      layoutNodeCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, connectivityInfo(i)%layoutHandle )
      CALL FieldmlUtil_CheckError( "Cannot get layout node count", fieldmlInfo, err, errorString, *999 )
      
      sourceHandle = Fieldml_GetDataSource( fieldmlInfo%fmlHandle, connectivityInfo(i)%connectivityHandle )
      IF( i == 1 ) THEN
        writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, sourceHandle, 0 )
      ELSE
        writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, sourceHandle, 1 )
      ENDIF
      CALL FieldmlUtil_CheckError( "Cannot open connectivity data writer", fieldmlInfo, err, errorString, *999 )
      
      ALLOCATE( iBuffer( layoutNodeCount ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate layout buffer.", err, errorString, *999 )
      DO j = 1, elementCount
        CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET( j, meshElements, basis, err, errorString, *999 )
  
        CALL FieldmlUtil_GetConnectivityEnsemble( fieldmlInfo%fmlHandle, basis, layoutHandle, err, errorString, *999 )
        IF( layoutHandle == connectivityInfo(i)%layoutHandle ) THEN
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET( j, meshElements, iBuffer, err, errorString, *999 )
        ELSE
          iBuffer = 0
        ENDIF
        writeCount = Fieldml_WriteIntValues( fieldmlInfo%fmlHandle, writer, C_LOC(iBuffer), layoutNodeCount )
        IF( writeCount /= layoutNodeCount ) THEN
          CALL FLAG_ERROR( "Cannot write connectivity data", err, errorString, *999 )
        ENDIF
      ENDDO
      DEALLOCATE( iBuffer )
      fmlErr = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close connectivity data writer", fieldmlInfo, err, errorString, *999 )
    ENDDO
    
    IF( ALLOCATED( basisInfo ) ) THEN
      DEALLOCATE( basisInfo )
    ENDIF
    IF( ALLOCATED( connectivityInfo ) ) THEN
      DEALLOCATE( connectivityInfo )
    ENDIF
    
    CALL LIST_ITEM_SET_C_INT( fieldmlInfo%componentHandles, componentNumber, templateHandle, err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_AddMeshComponent" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddMeshComponent", err, errorString )
    CALL EXITS( "FieldmlOutput_AddMeshComponent" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_AddMeshComponent

  !
  !================================================================================================================================
  !
  
  !<Create a parameter evaluator and associated data source containing the nodal dofs for the given field components.
  SUBROUTINE FieldmlOutput_AddFieldNodeDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
    & variableType, nodeDofsHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    INTEGER(C_INT), INTENT(IN) :: typeHandle !<The FieldML type handle for the field.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:) !<The field component numbers for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: variableType !<The OpenCMISS variable type to generate dofs for.
    INTEGER(C_INT), INTENT(INOUT) :: nodeDofsHandle !<The handle of the nodal dofs parameter evaluator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    TYPE(MESH_TYPE), POINTER :: mesh
    INTEGER(C_INT) :: typeComponentHandle, real1DHandle, nodeCount, indexHandle, resourceHandle, sourceHandle
    INTEGER(INTG) :: versionNumber,componentCount, i, j, interpolationType, globalNodeNumber
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    INTEGER(C_INT) :: writer, writeCount, fmlErr
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    REAL(C_DOUBLE) :: dValue
    LOGICAL :: nodeExists
    LOGICAL, ALLOCATABLE :: isNodeBased(:)

    CALL ENTERS( "FieldmlOutput_AddFieldNodeDofs", err, errorString, *999 )
    
    mesh => field%DECOMPOSITION%MESH
    
    CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, real1DHandle, .TRUE., err, errorString, *999 )

    componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    typeComponentHandle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
    nodeCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle )
    
    ALLOCATE( meshComponentNumbers( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate mesh component array.", err, errorString, *999 )
    ALLOCATE( isNodeBased( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate nodal component array.", err, errorString, *999 )

    DO i = 1, componentCount
      CALL FIELD_COMPONENT_MESH_COMPONENT_GET( field, variableType, fieldComponentNumbers(i), &
        & meshComponentNumbers(i), err, errorString, *999 )
      CALL FIELD_COMPONENT_INTERPOLATION_GET( field, variableType, fieldComponentNumbers(i), interpolationType, &
        & err, errorString, *999 )
        
      isNodeBased( i ) = ( interpolationType == FIELD_NODE_BASED_INTERPOLATION )
    ENDDO

    resourceHandle = Fieldml_CreateTextFileDataResource( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.node.resource"), &
      & cchar(baseName//".dofs.node") )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs data resource", fieldmlInfo, err, errorString, *999 )
    
    nodeDofsHandle = Fieldml_CreateParameterEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.node"), real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs parameter set", fieldmlInfo, err, errorString, *999 )
    fmlErr = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, nodeDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs parameter description", fieldmlInfo, err, errorString, *999 )
    
    sourceHandle = Fieldml_CreateTextDataSource( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.node.data"), resourceHandle, &
      & 1, nodeCount, componentCount, 0, 0 )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs data source", fieldmlInfo, err, errorString, *999 )
    
    fmlErr = Fieldml_SetDataSource( fieldmlInfo%fmlHandle, nodeDofsHandle, sourceHandle )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs data source", fieldmlInfo, err, errorString, *999 )

    IF( typeComponentHandle /= FML_INVALID_HANDLE ) THEN
      typeComponentHandle = FieldmlUtil_ImportHandle( fieldmlInfo%fmlHandle, typeComponentHandle )
      indexHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, typeComponentHandle, .TRUE. )
      fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, nodeDofsHandle, indexHandle, FML_INVALID_HANDLE )
      CALL FieldmlUtil_CheckError( "Cannot add component index for nodal dofs parameter set", fieldmlInfo, &
        & err, errorString, *999 )
    ENDIF
    fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, nodeDofsHandle, fieldmlInfo%nodesArgumentHandle, &
      & FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError( "Cannot add layout index for nodal dofs parameter set", fieldmlInfo, err, errorString, *999 )

    ALLOCATE( dBuffer( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate nodal dofs array.", err, errorString, *999 )
    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, sourceHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open nodal parameter writer", fieldmlInfo, err, errorString, *999 )
    DO i = 1, nodeCount
      DO j = 1, componentCount
        dValue = 0
        IF( isNodeBased(j) ) THEN
          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS( mesh, meshComponentNumbers(j), i, nodeExists, globalNodeNumber, &
            & err, errorString, *999 )
          IF( nodeExists ) THEN
            !Default to version 1 of each node derivative (value hardcoded in loop)
            versionNumber = 1
            CALL FIELD_PARAMETER_SET_GET_NODE( field, variableType, FIELD_VALUES_SET_TYPE, versionNumber, &
              & NO_GLOBAL_DERIV, i, fieldComponentNumbers(j), dValue, err, errorString, *999 )
          ENDIF
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      writeCount = Fieldml_WriteDoubleValues( fieldmlInfo%fmlHandle, writer, C_LOC(dBuffer), componentCount )
      IF( writeCount /= componentCount ) THEN
        CALL FLAG_ERROR( "Cannot write nodal parameter values", err, errorString, *999 )
      ENDIF
    ENDDO
    fmlErr = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
    CALL FieldmlUtil_CheckError( "Cannot close nodal parameter writer", fieldmlInfo, err, errorString, *999 )
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
  
  !<Create a parameter evaluator and associated data source containing the element dofs for the given field components.
  SUBROUTINE FieldmlOutput_AddFieldElementDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
    & variableType, elementDofsHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    INTEGER(C_INT), INTENT(IN) :: typeHandle !<The FieldML type handle for the field.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:) !<The field component numbers for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: variableType !<The OpenCMISS variable type to generate dofs for.
    INTEGER(C_INT), INTENT(INOUT) :: elementDofsHandle !<The handle of the element dofs parameter evaluator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: typeComponentHandle, real1DHandle, elementCount, indexHandle, resourceHandle, sourceHandle
    INTEGER(INTG) :: componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    INTEGER(C_INT) :: writer, writeCount, fmlErr
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: dBuffer(:)
    REAL(C_DOUBLE) :: dValue
    LOGICAL, ALLOCATABLE :: isElementBased(:)

    CALL EXITS( "FieldmlOutput_AddFieldElementDofs" )
    
    CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, real1DHandle, .TRUE., err, errorString, *999 )

    componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    typeComponentHandle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
    
    elementCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    
    ALLOCATE( meshComponentNumbers( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate mesh component number array.", err, errorString, *999 )
    ALLOCATE( isElementBased( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate element component array.", err, errorString, *999 )

    DO i = 1, componentCount
      CALL FIELD_COMPONENT_MESH_COMPONENT_GET( field, variableType, fieldComponentNumbers(i), &
        & meshComponentNumbers(i), err, errorString, *999 )
      CALL FIELD_COMPONENT_INTERPOLATION_GET( field, variableType, fieldComponentNumbers(i), interpolationType, &
        & err, errorString, *999 )

      isElementBased( i ) = ( interpolationType == FIELD_ELEMENT_BASED_INTERPOLATION )
    ENDDO

    resourceHandle = Fieldml_CreateTextFileDataResource( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.element.resource"), &
      & cchar(baseName//".dofs.element") )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs data resource", fieldmlInfo, err, errorString, *999 )
    
    elementDofsHandle = Fieldml_CreateParameterEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.element"), &
      & real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs parameter set", fieldmlInfo, err, errorString, *999 )
    fmlErr = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, elementDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set element dofs parameter description", fieldmlInfo, err, errorString, *999 )

    sourceHandle = Fieldml_CreateTextDataSource( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.element.data"), &
      & resourceHandle, 1, elementCount, componentCount, 0, 0 )
    CALL FieldmlUtil_CheckError( "Cannot create element dofs data source", fieldmlInfo, err, errorString, *999 )
    
    fmlErr = Fieldml_SetDataSource( fieldmlInfo%fmlHandle, elementDofsHandle, sourceHandle )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs data source", fieldmlInfo, err, errorString, *999 )

    IF( typeComponentHandle /= FML_INVALID_HANDLE ) THEN
      typeComponentHandle = FieldmlUtil_ImportHandle( fieldmlInfo%fmlHandle, typeComponentHandle )
      indexHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, typeComponentHandle, .TRUE. )
      fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, elementDofsHandle, typeComponentHandle, FML_INVALID_HANDLE )
      CALL FieldmlUtil_CheckError( "Cannot add component index for element dofs parameter set", fieldmlInfo, &
        & err, errorString, *999 )
    ENDIF
    fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, elementDofsHandle, fieldmlInfo%elementsArgumentHandle, &
      & FML_INVALID_HANDLE )
    CALL FieldmlUtil_CheckError( "Cannot add element index for element dofs parameter set", fieldmlInfo, &
      & err, errorString, *999 )

    ALLOCATE( dBuffer( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate element dofs buffer.", err, errorString, *999 )
    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, sourceHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open element parameter writer", fieldmlInfo, err, errorString, *999 )
    DO i = 1, elementCount
      DO j = 1, componentCount
        dValue = 0
        IF( isElementBased(j) ) THEN
          CALL FIELD_PARAMETER_SET_GET_ELEMENT( field, variableType, FIELD_VALUES_SET_TYPE, i, &
            & fieldComponentNumbers(j), dValue, err, errorString, *999 )
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      writeCount = Fieldml_WriteDoubleValues( fieldmlInfo%fmlHandle, writer, C_LOC(dBuffer), componentCount )
      IF( writeCount /= componentCount ) THEN
        CALL FLAG_ERROR( "Cannot write element parameter values", err, errorString, *999 )
      ENDIF
    ENDDO
    fmlErr = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
    CALL FieldmlUtil_CheckError( "Cannot close element parameter writer", fieldmlInfo, err, errorString, *999 )
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
  
  !<Create a parameter evaluator and associated data source containing the globally constant dofs for the given field components.
  SUBROUTINE FieldmlOutput_AddFieldConstantDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
    & variableType, constantDofsHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    INTEGER(C_INT), INTENT(IN) :: typeHandle !<The FieldML type handle for the field.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:) !<The field component numbers for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: variableType !<The OpenCMISS variable type to generate dofs for.
    INTEGER(C_INT), INTENT(INOUT) :: constantDofsHandle !<The handle of the constant dofs parameter evaluator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: doftypeHandle, typeType, componentType, dataType, indexHandle, resourceHandle, sourceHandle
    INTEGER(INTG) :: componentCount, i, j, interpolationType
    INTEGER(INTG), ALLOCATABLE :: meshComponentNumbers(:)
    INTEGER(C_INT) :: writer, writeCount, fmlErr
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
      CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, doftypeHandle, .TRUE., err, errorString, *999 )
      componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
      componentType = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
      isReal = .TRUE.
    ENDIF
    
    ALLOCATE( meshComponentNumbers( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate mesh component array.", err, errorString, *999 )
    ALLOCATE( isConstant( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate constant component array.", err, errorString, *999 )

    DO i = 1, componentCount
      CALL FIELD_COMPONENT_MESH_COMPONENT_GET( field, variableType, fieldComponentNumbers(i), &
        & meshComponentNumbers(i), err, errorString, *999 )
      CALL FIELD_COMPONENT_INTERPOLATION_GET( field, variableType, fieldComponentNumbers(i), interpolationType, &
        & err, errorString, *999 )

      isConstant( i ) = ( interpolationType == FIELD_CONSTANT_INTERPOLATION )
    ENDDO

    resourceHandle = Fieldml_CreateTextFileDataResource( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.constant.resource"), &
      & cchar(baseName//".dofs.constant") )
    CALL FieldmlUtil_CheckError( "Cannot create constant dofs data resource", fieldmlInfo, err, errorString, *999 )
    
    constantDofsHandle = Fieldml_CreateParameterEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.constant"), &
      & doftypeHandle )
    CALL FieldmlUtil_CheckError( "Cannot create constant dofs parameter set", fieldmlInfo, err, errorString, *999 )
    fmlErr = Fieldml_SetParameterDataDescription( fieldmlInfo%fmlHandle, constantDofsHandle, DESCRIPTION_SEMIDENSE )
    CALL FieldmlUtil_CheckError( "Cannot set constant dofs parameter description", fieldmlInfo, err, errorString, *999 )

    sourceHandle = Fieldml_CreateTextDataSource( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.constant.data"), &
      & resourceHandle, 1, 1, componentCount, 0, 0 )
    CALL FieldmlUtil_CheckError( "Cannot create constant dofs data source", fieldmlInfo, err, errorString, *999 )
    
    fmlErr = Fieldml_SetDataSource( fieldmlInfo%fmlHandle, constantDofsHandle, sourceHandle )
    CALL FieldmlUtil_CheckError( "Cannot set nodal dofs data source", fieldmlInfo, err, errorString, *999 )

    IF( componentType /= FML_INVALID_HANDLE ) THEN
      componentType = FieldmlUtil_ImportHandle( fieldmlInfo%fmlHandle, componentType )
      indexHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, componentType, .TRUE. )
      fmlErr = Fieldml_AddDenseIndexEvaluator( fieldmlInfo%fmlHandle, constantDofsHandle, indexHandle, FML_INVALID_HANDLE )
      CALL FieldmlUtil_CheckError( "Cannot add component index for constant dofs parameter set", fieldmlInfo, &
        & err, errorString, *999 )
    ENDIF

    writer = Fieldml_OpenWriter( fieldmlInfo%fmlHandle, sourceHandle, 0 )
    CALL FieldmlUtil_CheckError( "Cannot open constant parameter writer", fieldmlInfo, err, errorString, *999 )

    CALL FIELD_DATA_TYPE_GET( field, variableType, dataType, err, errorString, *999 )
    IF( dataType == FIELD_INTG_TYPE ) THEN
      isReal = .false.
    ELSEIF( dataType == FIELD_DP_TYPE ) THEN
      isReal = .true.
    ENDIF
    
    IF( isReal ) THEN
      ALLOCATE( dBuffer( componentCount ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate constant dofs buffer.", err, errorString, *999 )
      DO j = 1, componentCount
        dValue = 0
        IF( isConstant(j) ) THEN
          CALL FIELD_PARAMETER_SET_GET_CONSTANT( field, variableType, FIELD_VALUES_SET_TYPE, &
            & fieldComponentNumbers(j), dValue, err, errorString, *999 )
        ENDIF
        dBuffer( j ) = dValue
      ENDDO
      writeCount = Fieldml_WriteDoubleValues( fieldmlInfo%fmlHandle, writer, C_LOC(dBuffer), componentCount )
      IF( writeCount /= componentCount ) THEN
        CALL FLAG_ERROR( "Cannot write constant parameter values", err, errorString, *999 )
      ENDIF
      fmlErr = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close constant parameter writer", fieldmlInfo, err, errorString, *999 )
      DEALLOCATE( dBuffer )
    ELSE
      ALLOCATE( iBuffer( componentCount ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate constant dofs buffer.", err, errorString, *999 )
      DO j = 1, componentCount
        iValue = 0
        IF( isConstant(j) ) THEN
          CALL FIELD_PARAMETER_SET_GET_CONSTANT( field, variableType, FIELD_VALUES_SET_TYPE, &
            & fieldComponentNumbers(j), iValue, err, errorString, *999 )
        ENDIF
        iBuffer( j ) = iValue
      ENDDO
      writeCount = Fieldml_WriteIntValues( fieldmlInfo%fmlHandle, writer, C_LOC(iBuffer), componentCount )
      IF( writeCount /= componentCount ) THEN
        CALL FLAG_ERROR( "Cannot write constant parameter values", err, errorString, *999 )
      ENDIF
      fmlErr = Fieldml_CloseWriter( fieldmlInfo%fmlHandle, writer )
      CALL FieldmlUtil_CheckError( "Cannot close constant parameter writer", fieldmlInfo, err, errorString, *999 )
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

  !<Initialize the given FieldML parsing state for use with the given mesh
  SUBROUTINE FieldmlOutput_InitialiseInfo( mesh, location, baseName, fieldmlInfo, err, errorString, * )
    !Argument variables
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: mesh !<The mesh with which the FieldML document is associated.
    TYPE(VARYING_STRING), INTENT(IN) :: location !<The location of the FieldML file. Data resources will be created here.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    TYPE(FieldmlInfoType), INTENT(OUT) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    TYPE(REGION_TYPE), POINTER :: region
    INTEGER(INTG) :: componentCount, i, nodeCount, elementCount, dimensions
    INTEGER(C_INT) :: real1DHandle, xiComponentHandle, fmlErr
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: meshElements
    TYPE(NODES_TYPE), POINTER :: nodes

    CALL ENTERS( "FieldmlOutput_InitialiseInfo", err, errorString, *999 )
    
    region => mesh%REGION
    
    dimensions = mesh%NUMBER_OF_DIMENSIONS
    
    CALL FieldmlUtil_InitialiseInfo( fieldmlInfo, err, errorString, *999 )
    
    fieldmlInfo%fmlHandle = Fieldml_Create( cchar(location), cchar(baseName) )
    CALL FieldmlUtil_CheckError( "Cannot create fieldml handle", fieldmlInfo, err, errorString, *999 )

    NULLIFY( nodes )
    CALL REGION_NODES_GET( region, nodes, err, errorString, *999 )
    CALL NODES_NUMBER_OF_NODES_GET( nodes, nodeCount, err, errorString, *999 )

    fieldmlInfo%nodesHandle = Fieldml_CreateEnsembleType( fieldmlInfo%fmlHandle, cchar(baseName//".nodes") )
    CALL FieldmlUtil_CheckError( "Cannot create mesh nodes ensemble", fieldmlInfo, err, errorString, *999 )
    fmlErr = Fieldml_SetEnsembleMembersRange( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle, 1, nodeCount, 1 )
    CALL FieldmlUtil_CheckError( "Cannot set mesh nodes ensemble bounds", fieldmlInfo, err, errorString, *999 )
    
    fieldmlInfo%nodesArgumentHandle = Fieldml_CreateArgumentEvaluator( fieldmlInfo%fmlHandle, &
      & cchar(baseName//".nodes.argument"), fieldmlInfo%nodesHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh nodes variable", fieldmlInfo, err, errorString, *999 )
    
    CALL MESH_NUMBER_OF_ELEMENTS_GET( mesh, elementCount, err, errorString, *999 )

    fieldmlInfo%meshHandle = Fieldml_CreateMeshType( fieldmlInfo%fmlHandle, cchar(baseName//".mesh") )
    CALL FieldmlUtil_CheckError( "Cannot create mesh type", fieldmlInfo, err, errorString, *999 )

    fieldmlInfo%elementsHandle = Fieldml_CreateMeshElementsType( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "element"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot create mesh elements type", fieldmlInfo, err, errorString, *999 )
    fmlErr = Fieldml_SetEnsembleMembersRange( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle, 1, elementCount, 1 )
    CALL FieldmlUtil_CheckError( "Cannot set mesh type element count", fieldmlInfo, err, errorString, *999 )

    fieldmlInfo%xiHandle = Fieldml_CreateMeshChartType( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "xi"//NUL )
    CALL FieldmlUtil_CheckError( "Cannot create mesh chart type", fieldmlInfo, err, errorString, *999 )
    xiComponentHandle = Fieldml_CreateContinuousTypeComponents( fieldmlInfo%fmlHandle, fieldmlInfo%xiHandle, &
      & cchar(baseName//".mesh.xi.component"), dimensions )
    CALL FieldmlUtil_CheckError( "Cannot create mesh chart type", fieldmlInfo, err, errorString, *999 )
    
    fmlErr = Fieldml_CreateArgumentEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".mesh.argument"), &
      & fieldmlInfo%meshHandle )
    CALL FieldmlUtil_CheckError( "Cannot create mesh variable", fieldmlInfo, err, errorString, *999 )

    fieldmlInfo%xiArgumentHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(baseName//".mesh.argument.xi") )
    CALL FieldmlUtil_CheckError( "Cannot get mesh xi variable", fieldmlInfo, err, errorString, *999 )
    fieldmlInfo%elementsArgumentHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, &
      & cchar(baseName//".mesh.argument.element") )
    CALL FieldmlUtil_CheckError( "Cannot get mesh element variable", fieldmlInfo, err, errorString, *999 )
    
    CALL FieldmlUtil_GetGenericType( fieldmlInfo%fmlHandle, 1, real1DHandle, .TRUE., err, errorString, *999 )
    
    !TODO Some of these may end up being unused. Should use deferred assignment.
    fieldmlInfo%nodeDofsHandle = Fieldml_CreateArgumentEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.node"), &
      & real1DHandle )
    CALL FieldmlUtil_CheckError( "Cannot create nodal dofs variable", fieldmlInfo, err, errorString, *999 )
!    fieldmlInfo%elementDofsHandle = Fieldml_CreateArgumentEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.element"), & 
!      & real1DHandle )
!    CALL FieldmlUtil_CheckError( "Cannot create element dofs variable", fieldmlInfo, err, errorString, *999 )
!    fieldmlInfo%constantDofsHandle = Fieldml_CreateArgumentEvaluator( fieldmlInfo%fmlHandle, cchar(baseName//".dofs.constant"), & 
!      & real1DHandle )
!    CALL FieldmlUtil_CheckError( "Cannot create constant dofs variable", fieldmlInfo, err, errorString, *999 )

    CALL MESH_NUMBER_OF_COMPONENTS_GET( mesh, componentCount, err, errorString, *999 )
    DO i = 1, componentCount
      NULLIFY( meshElements )
      CALL LIST_ITEM_ADD_C_INT( fieldmlInfo%componentHandles, FML_INVALID_HANDLE, err, errorString, *999 )
      CALL MESH_TOPOLOGY_ELEMENTS_GET( mesh, i, meshElements, err, errorString, *999 )
      CALL FieldmlOutput_AddMeshComponent( fieldmlInfo, baseName, i, meshElements, err, errorString, *999 )
    ENDDO
    
    !TODO Proper shape assignment.
    IF( dimensions == 2 ) THEN
      fmlErr = Fieldml_SetMeshDefaultShape( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "shape.square"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot set 2D mesh type element shape", fieldmlInfo, err, errorString, *999 )
    ELSE
      fmlErr = Fieldml_SetMeshDefaultShape( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle, "shape.cube"//NUL )
      CALL FieldmlUtil_CheckError( "Cannot set 3D mesh type element shape", fieldmlInfo, err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlOutput_InitialiseInfo" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_InitialiseInfo", err, errorString )
    CALL EXITS( "FieldmlOutput_InitialiseInfo" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_InitialiseInfo

  !
  !================================================================================================================================
  !
  
  !< Add the components of the given field to the given FieldML evaluator, creating component templates as needed.
  SUBROUTINE FieldmlOutput_AddFieldComponents( fieldmlInfo, typeHandle, baseName, field, fieldComponentNumbers, &
    & variableType, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: typeHandle !<The FieldML type handle for the field.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field for which dof components are to be created.
    INTEGER(INTG), INTENT(IN) :: fieldComponentNumbers(:) !<The field component numbers for which evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: variableType !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    TYPE(MESH_TYPE), POINTER :: mesh
    INTEGER(C_INT) :: fieldHandle, componentHandle, nodalDofsHandle, elementDofsHandle, constantDofsHandle, indexHandle, fmlErr
    INTEGER(INTG) :: componentCount, i, meshComponentNumber, interpolationType
    INTEGER(C_INT), ALLOCATABLE, TARGET :: componentEvaluators(:)
  
    CALL ENTERS( "FieldmlOutput_AddFieldComponents", err, errorString, *999 )
    
    mesh => field%DECOMPOSITION%MESH

    componentHandle = Fieldml_GetTypeComponentEnsemble( fieldmlInfo%fmlHandle, typeHandle )
    componentCount = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    ALLOCATE( componentEvaluators( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate component evaluators array.", err, errorString, *999 )

    IF( SIZE( fieldComponentNumbers ) /= componentCount ) THEN
      CALL FLAG_ERROR( "Fieldml Component count must match value type component count", err, errorString, *999 )
    ENDIF

    nodalDofsHandle = FML_INVALID_HANDLE
    elementDofsHandle = FML_INVALID_HANDLE
    constantDofsHandle = FML_INVALID_HANDLE
    !TODO Other types of interpolation not yet supported.
    DO i = 1, componentCount
      CALL FIELD_COMPONENT_INTERPOLATION_GET( field, variableType, fieldComponentNumbers(i), interpolationType, &
        & err, errorString, *999 )
        
      IF( interpolationType == FIELD_NODE_BASED_INTERPOLATION ) THEN
        IF( nodalDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldNodeDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
          & variableType, nodalDofsHandle, err, errorString, *999 )
        ENDIF
        CALL FIELD_COMPONENT_MESH_COMPONENT_GET( field, variableType, fieldComponentNumbers(i), &
          & meshComponentNumber, err, errorString, *999 )
        CALL LIST_ITEM_GET_C_INT( fieldmlInfo%componentHandles, meshComponentNumber, componentEvaluators( i ), &
          & err, errorString, *999 )
      ELSEIF( interpolationType == FIELD_ELEMENT_BASED_INTERPOLATION ) THEN
        IF( elementDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldElementDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
            & variableType, elementDofsHandle, err, errorString, *999 )
        ENDIF
        componentEvaluators( i ) = elementDofsHandle
      ELSEIF( interpolationType == FIELD_CONSTANT_INTERPOLATION ) THEN
        IF( constantDofsHandle == FML_INVALID_HANDLE ) THEN
          CALL FieldmlOutput_AddFieldConstantDofs( fieldmlInfo, baseName, typeHandle, field, fieldComponentNumbers, &
            & variableType, constantDofsHandle, err, errorString, *999 )
        ENDIF
        componentEvaluators( i ) = constantDofsHandle
      ENDIF
    ENDDO
    

    IF( componentHandle /= FML_INVALID_HANDLE ) THEN
      fieldHandle = Fieldml_CreateAggregateEvaluator( fieldmlInfo%fmlHandle, cchar(baseName), typeHandle )
      CALL FieldmlUtil_CheckError( "Cannot create aggregate evaluator for field", fieldmlInfo, err, errorString, *999 )
      indexHandle = FieldmlUtil_GetTypeArgumentHandle( fieldmlInfo, componentHandle, .TRUE. )
      fmlErr = Fieldml_SetIndexEvaluator( fieldmlInfo%fmlHandle, fieldHandle, 1, indexHandle )
      CALL FieldmlUtil_CheckError( "Cannot set index evaluator for aggregate field component", fieldmlInfo, &
        & err, errorString, *999 )

      DO i = 1, componentCount
        fmlErr = Fieldml_SetEvaluator( fieldmlInfo%fmlHandle, fieldHandle, i, componentEvaluators( i ) )
        CALL FieldmlUtil_CheckError( "Cannot set nodal evaluator for aggregate field component", fieldmlInfo, &
          & err, errorString, *999 )
      ENDDO
    ELSE
      fieldHandle = Fieldml_CreateReferenceEvaluator( fieldmlInfo%fmlHandle, cchar(baseName), componentEvaluators( 1 ) )
      CALL FieldmlUtil_CheckError( "Cannot create aggregate evaluator for field", fieldmlInfo, err, errorString, *999 )
    ENDIF

    IF( nodalDofsHandle /= FML_INVALID_HANDLE ) THEN
      fmlErr = Fieldml_SetBind( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%nodeDofsHandle, nodalDofsHandle )
      CALL FieldmlUtil_CheckError( "Cannot set nodal dofs bind for field with interpolated elements", fieldmlInfo, err, &
        & errorString, *999 )
    ENDIF
!    IF( elementDofsHandle /= FML_INVALID_HANDLE ) THEN
!      fmlErr = Fieldml_SetBind( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%elementDofsHandle, elementDofsHandle )
!      CALL FieldmlUtil_CheckError( "Cannot set element dofs bind for field with interpolated elements", fieldmlInfo, &
!  &err, errorString, *999 )
!    ENDIF
!    IF( constantDofsHandle /= FML_INVALID_HANDLE ) THEN
!      fmlErr = Fieldml_SetBind( fieldmlInfo%fmlHandle, fieldHandle, fieldmlInfo%constantDofsHandle, constantDofsHandle )
!      CALL FieldmlUtil_CheckError( "Cannot set constant dofs bind for field with interpolated elements", fieldmlInfo, &
!  &err, errorString, *999 )
!    ENDIF


    DEALLOCATE( componentEvaluators )
    CALL EXITS( "FieldmlOutput_AddFieldComponents" )
    RETURN
999 DEALLOCATE( componentEvaluators )
    CALL ERRORS( "FieldmlOutput_AddFieldComponents", err, errorString )
    CALL EXITS( "FieldmlOutput_AddFieldComponents" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_AddFieldComponents

  !
  !================================================================================================================================
  !

  !<Add the given field to the given FieldML document. The field's type will be determined by FieldmlUtil_GetValueType. \see Fieldml_Util_Routines::FieldmlUtil_GetValueType
  SUBROUTINE FieldmlOutput_AddField_NoType( fieldmlInfo, baseName, field, variableType, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field for which evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: variableType !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: typeHandle
    
    CALL ENTERS( "FieldmlOutput_AddField", err, errorString, *999 )

    CALL FieldmlUtil_GetValueType( fieldmlInfo%fmlHandle, field, typeHandle, .TRUE., err, errorString, *999 )

    CALL FieldmlOutput_AddField_WithType( fieldmlInfo, baseName, field, variableType, typeHandle, err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_AddField_NoType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddField_NoType", err, errorString )
    CALL EXITS( "FieldmlOutput_AddField_NoType" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_AddField_NoType

  !
  !================================================================================================================================
  !

  !<Add the given field to the given FieldML document using the given FieldML type.
  SUBROUTINE FieldmlOutput_AddField_WithType( fieldmlInfo, baseName, field, variableType, typeHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: baseName !<The root name of the basis evaluator.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: field !<The field for which evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: variableType !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: typeHandle !<The FieldML type handle for the field.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(INTG) :: i, componentCount
    INTEGER(INTG), ALLOCATABLE :: fieldComponentNumbers(:)
    TYPE(MESH_TYPE), POINTER :: mesh
    
    mesh => field%DECOMPOSITION%MESH

    IF( typeHandle == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Cannot get value type for field.", err, errorString, *999 )
    ENDIF
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( field, FIELD_U_VARIABLE_TYPE, componentCount, err, errorString, *999 )
    
    ALLOCATE( fieldComponentNumbers( componentCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate component numbers array.", err, errorString, *999 )
    DO i = 1, componentCount
      fieldComponentNumbers(i) = i
    ENDDO

    CALL FieldmlOutput_AddFieldComponents( fieldmlInfo, typeHandle, baseName, field, fieldComponentNumbers, &
      & variableType, err, errorString, *999 )
    
    DEALLOCATE( fieldComponentNumbers )

    CALL EXITS( "FieldmlOutput_AddField_WithType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_AddField_WithType", err, errorString )
    CALL EXITS( "FieldmlOutput_AddField_WithType" )
    RETURN 1
    
  END SUBROUTINE FieldmlOutput_AddField_WithType

  !
  !================================================================================================================================
  !

  !<Create a FieldML ensemble with the given name and number of members.
  SUBROUTINE FieldmlOutput_CreateEnsembleType( fieldmlInfo, typeName, elementCount, typeHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: typeName !<The name of the ensemble type.
    INTEGER(INTG), INTENT(IN) :: elementCount !<The number of members in the ensemble.
    INTEGER(INTG), INTENT(OUT) :: typeHandle !<The FieldML type handle for the ensemble.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: fmlErr

    CALL ENTERS( "FieldmlOutput_CreateEnsembleType", err, errorString, *999 )

    typeHandle = Fieldml_CreateEnsembleType( fieldmlInfo%fmlHandle, cchar(typeName) )
    CALL FieldmlUtil_CheckError( "Error creating ensemble type", fieldmlInfo, err, errorString, *999 )
    
    fmlErr = Fieldml_SetEnsembleMembersRange( fieldmlInfo%fmlHandle, typeHandle, 1, elementCount, 1 )
    CALL FieldmlUtil_CheckError( "Error setting ensemble type bounds", fieldmlInfo, err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_CreateEnsembleType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateEnsembleType", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateEnsembleType" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_CreateEnsembleType

  !
  !================================================================================================================================
  !

  !<Create a FieldML continuous type with the given name and number of components.
  SUBROUTINE FieldmlOutput_CreateContinuousType( fieldmlInfo, typeName, componentCount, typeHandle, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: typeName !<The name of the continuous type.
    INTEGER(INTG), INTENT(IN) :: componentCount !<The number of components in the continuous type.
    INTEGER(INTG), INTENT(OUT) :: typeHandle !<The FieldML type handle for the continuous type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Local variables
    INTEGER(INTG) :: componentHandle

    CALL ENTERS( "FieldmlOutput_CreateContinuousType", err, errorString, *999 )

    componentHandle = FML_INVALID_HANDLE

    typeHandle = Fieldml_CreateContinuousType( fieldmlInfo%fmlHandle, cchar(typeName) )
    CALL FieldmlUtil_CheckError( "Error creating continuous type", fieldmlInfo, err, errorString, *999 )

    IF( componentCount > 1 ) THEN
      componentHandle = Fieldml_CreateContinuousTypeComponents( fieldmlInfo%fmlHandle, typeHandle, &
        & cchar(typeName//".component"), componentCount )
      CALL FieldmlUtil_CheckError( "Error creating component type", fieldmlInfo, err, errorString, *999 )
    ENDIF

    CALL EXITS( "FieldmlOutput_CreateContinuousType" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_CreateContinuousType", err, errorString )
    CALL EXITS( "FieldmlOutput_CreateContinuousType" )
    RETURN 1

  END SUBROUTINE FieldmlOutput_CreateContinuousType

  !
  !================================================================================================================================
  !

  !<Write the given FieldML document to the given file.
  SUBROUTINE FieldmlOutput_Write( fieldmlInfo, filename, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: filename !<The file to write the FieldML document to.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: fmlErr

    CALL ENTERS( "FieldmlOutput_Write", err, errorString, *999 )

    fmlErr = Fieldml_WriteFile( fieldmlInfo%fmlHandle, cchar(filename) )
    CALL FieldmlUtil_CheckError( "Error writing fieldml file", fieldmlInfo, err, errorString, *999 )

    CALL EXITS( "FieldmlOutput_Write" )
    RETURN
999 CALL ERRORS( "FieldmlOutput_Write", err, errorString )
    CALL EXITS( "FieldmlOutput_Write" )
    RETURN 1
  
  END SUBROUTINE

  !
  !================================================================================================================================
  !

END MODULE FIELDML_OUTPUT_ROUTINES
