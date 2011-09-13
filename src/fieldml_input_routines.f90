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

  USE BASIS_ROUTINES
  USE CMISS
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE FIELDML_API
  USE FIELDML_TYPES
  USE FIELDML_UTIL_ROUTINES
  USE LISTS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE REGION_ROUTINES
  USE UTIL_ARRAY

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  CHARACTER(KIND=C_CHAR), PARAMETER :: NUL = C_NULL_CHAR

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
  
  !>Determines the connectivity evaluator and layout argument for the given basis.
  SUBROUTINE FieldmlInput_GetBasisConnectivityInfo( fieldmlInfo, basisHandle, paramArgHandle, connectivityHandle, layoutHandle, &
    & err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: basisHandle !<The basis handle.
    INTEGER(C_INT), INTENT(IN) :: paramArgHandle !<The basis parameters argument handle.
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle !<The basis connectivity evaluator handle.
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle !<The local node layout.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Local variables
    INTEGER(C_INT) :: count, bindNumber, paramsHandle, argHandle, layoutIndexHandle
    
    CALL ENTERS( "FieldmlInput_GetBasisConnectivityInfo", err, errorString, *999 )

    count = Fieldml_GetBindCount( fieldmlInfo%fmlHandle, basisHandle )
    IF( count /= 2 ) THEN
      CALL FLAG_ERROR( "Library basis evaluators must have exactly two binds", err, errorString, *999 )
    END IF
    
    paramsHandle = FML_INVALID_HANDLE
    DO bindNumber = 1, count
      argHandle = Fieldml_GetBindArgument( fieldmlInfo%fmlHandle, basisHandle, bindNumber )
      CALL FieldmlUtil_CheckError( "Cannot get bind for interpolator", fieldmlInfo%fmlHandle, err, errorString, *999 )
      IF( argHandle == paramArgHandle ) THEN
        paramsHandle = Fieldml_GetBindEvaluator( fieldmlInfo%fmlHandle, basisHandle, bindNumber )
      ENDIF
    ENDDO

    IF( paramsHandle == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Library interpolators must have a correct parameter bind", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( fieldmlInfo%fmlHandle, paramsHandle ) /= FHT_AGGREGATE_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Parameter evaluator for interpolator must be an aggregate", err, errorString, *999 )
    ENDIF
    
    count = Fieldml_GetBindCount( fieldmlInfo%fmlHandle, paramsHandle )
    IF( count /= 1 ) THEN
      CALL FLAG_ERROR( "Nodal parameter evaluator must only have one bind", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetBindArgument( fieldmlInfo%fmlHandle, paramsHandle, 1 ) /= fieldmlInfo%nodesArgumentHandle ) THEN
      CALL FLAG_ERROR( "Nodal parameter evaluator must bind the nodes argument", err, errorString, *999 )
    ENDIF
    
    connectivityHandle = Fieldml_GetBindEvaluator( fieldmlInfo%fmlHandle, paramsHandle, 1 )
    CALL FieldmlUtil_CheckError( "Cannot get connectivity source for nodal parameters", fieldmlInfo%fmlHandle, &
      & err, errorString, *999 )
      
    layoutIndexHandle = Fieldml_GetIndexEvaluator( fieldmlInfo%fmlHandle, paramsHandle, 1 )
    layoutHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, layoutIndexHandle )
    CALL FieldmlUtil_CheckError( "Cannot get connectivity source for nodal parameters", fieldmlInfo%fmlHandle, &
      & err, errorString, *999 )

    CALL EXITS( "FieldmlInput_GetBasisConnectivityInfo" )
    RETURN
999 CALL ERRORS( "FieldmlInput_GetBasisConnectivityInfo", err, errorString )
    CALL EXITS( "FieldmlInput_GetBasisConnectivityInfo" )
    RETURN 1

  END SUBROUTINE FieldmlInput_GetBasisConnectivityInfo

  !
  !================================================================================================================================
  !
  
  !<Determine the basis collapse parameters from the given evaluator's name.
  SUBROUTINE FieldmlInput_GetBasisCollapse( name, collapse )
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: name !<The basis evaluator name.
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: collapse(:) !<The array of OpenCMISS basis collapse constants.

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

  !<Determines the basis configuration from the given basis evaluator.
  SUBROUTINE FieldmlInput_GetBasisInfo( fieldmlInfo, basisHandle, connectivityHandle, layoutHandle, basisType, &
    & basisInterpolations, collapse, err, errorString, * )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: basisHandle !<The basis evaluator handle.
    INTEGER(C_INT), INTENT(OUT) :: connectivityHandle !<The basis connectivity evaluator handle.
    INTEGER(C_INT), INTENT(OUT) :: layoutHandle !<The local node layout.
    INTEGER(INTG), INTENT(OUT) :: basisType !<The OpenCMISS basis type.
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: basisInterpolations(:) !<The per-xi basis interpolations (for TP bases).
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: collapse(:) !<The collapse constants for the basis.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: length, libraryBasisHandle, paramArgHandle
    CHARACTER(LEN=MAXSTRLEN) :: name
    TYPE(VARYING_STRING) :: collapseName
    
    CALL ENTERS( "FieldmlInput_GetBasisInfo", err, errorString, *999 )

    IF( .NOT. FieldmlInput_IsKnownBasis( fieldmlInfo, basisHandle ) ) THEN
      CALL FLAG_ERROR( "Basis specified in FieldML file is not yet supported", err, errorString, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( fieldmlInfo%fmlHandle, basisHandle ) /= FHT_REFERENCE_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Basis evaluator must be a continuous reference", err, errorString, *999 )
    ENDIF
    
    libraryBasisHandle = Fieldml_GetReferenceSourceEvaluator( fieldmlInfo%fmlHandle, basisHandle )
    CALL FieldmlUtil_CheckError( "Basis specified in FieldML is not a reference evaluator", fieldmlInfo, err, errorString, *999 )
    length = Fieldml_CopyObjectDeclaredName( fieldmlInfo%fmlHandle, libraryBasisHandle, name, MAXSTRLEN )
    CALL FieldmlUtil_CheckError( "Cannot get name of basis evaluator", fieldmlInfo, err, errorString, *999 )

    IF( INDEX( name, 'interpolator.3d.unit.triquadraticLagrange') == 1 ) THEN
      paramArgHandle = Fieldml_GetObjectByDeclaredName( fieldmlInfo%fmlHandle, &
        & "parameters.3d.unit.triquadraticLagrange.argument"//NUL )
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      CALL REALLOCATE_INT( collapse, 3, "", err, errorString, *999 )
      basisInterpolations = BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
      basisType = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( name, 'interpolator.3d.unit.trilinearLagrange') == 1 ) THEN
      paramArgHandle = Fieldml_GetObjectByDeclaredName( fieldmlInfo%fmlHandle, &
        & "parameters.3d.unit.trilinearLagrange.argument"//NUL )
      CALL REALLOCATE_INT( basisInterpolations, 3, "", err, errorString, *999 )
      CALL REALLOCATE_INT( collapse, 3, "", err, errorString, *999 )
      basisInterpolations = BASIS_LINEAR_LAGRANGE_INTERPOLATION
      basisType = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE
      CALL FLAG_ERROR( "Basis cannot yet be interpreted", err, errorString, *999 )
    ENDIF
    
    IF( basisType == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      collapseName = name(1:length)
      CALL FieldmlInput_GetBasisCollapse( collapseName, collapse )
    ENDIF
    
    CALL FieldmlInput_GetBasisConnectivityInfo( fieldmlInfo, basisHandle, paramArgHandle, connectivityHandle, layoutHandle, &
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

  !<Determines whether or not the given basis evaluator is known to OpenCMISS.
  FUNCTION FieldmlInput_IsKnownBasis( fieldmlInfo, basisHandle )
    !Argument variables
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: basisHandle !<The basis handle.
    
    !Function
    LOGICAL :: FieldmlInput_IsKnownBasis

    !Locals
    INTEGER(C_INT) :: length, libraryBasisHandle
    CHARACTER(LEN=MAXSTRLEN) :: name
    
    IF( Fieldml_GetObjectType( fieldmlInfo%fmlHandle, basisHandle ) /= FHT_REFERENCE_EVALUATOR ) THEN
      FieldmlInput_IsKnownBasis = .FALSE.
      RETURN
    ENDIF

    libraryBasisHandle = Fieldml_GetReferenceSourceEvaluator( fieldmlInfo%fmlHandle, basisHandle )
    length = Fieldml_CopyObjectDeclaredName( fieldmlInfo%fmlHandle, libraryBasisHandle, name, MAXSTRLEN )

    IF( ( INDEX( name, 'interpolator.3d.unit.triquadraticLagrange') /= 1 ) .AND. &
      & ( INDEX( name, 'interpolator.3d.unit.trilinearLagrange') /= 1 ) ) THEN
      FieldmlInput_IsKnownBasis = .FALSE.
      RETURN
    ENDIF
    
    FieldmlInput_IsKnownBasis = .TRUE.
    RETURN
    
  END FUNCTION FieldmlInput_IsKnownBasis
  
  !
  !================================================================================================================================
  !

  !<Determines whether or not the given evaluator is a recognisable mesh component evaluator.
  FUNCTION FieldmlInput_IsTemplateCompatible( fieldmlInfo, componentHandle, elementType )
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: componentHandle !<The mesh component evaluator.
    INTEGER(C_INT), INTENT(IN) :: elementType !<The element ensemble type.

    LOGICAL :: FieldmlInput_IsTemplateCompatible

    INTEGER(C_INT) :: objectType, count, i, evaluator, type, firstEvaluator, evaluatorHandle, defaultEvaluator

    objectType = Fieldml_GetObjectType( fieldmlInfo%fmlHandle, componentHandle )
    IF( objectType /= FHT_PIECEWISE_EVALUATOR ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    evaluatorHandle = Fieldml_GetIndexEvaluator( fieldmlInfo%fmlHandle, componentHandle, 1 )
    type = Fieldml_GetValueType( fieldmlInfo%fmlHandle, evaluatorHandle )
    IF( type /= elementType ) THEN
      FieldmlInput_IsTemplateCompatible = .TRUE.
      RETURN
    ENDIF

    count = Fieldml_GetEvaluatorCount( fieldmlInfo%fmlHandle, componentHandle )
    defaultEvaluator = Fieldml_GetDefaultEvaluator( fieldmlInfo%fmlHandle, componentHandle )
    
    IF( ( defaultEvaluator /= FML_INVALID_HANDLE ) .AND. .NOT. &
      & FieldmlInput_IsKnownBasis( fieldmlInfo, defaultEvaluator ) ) THEN
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

    firstEvaluator = Fieldml_GetEvaluator( fieldmlInfo%fmlHandle, componentHandle, 1 )
    IF( .NOT. FieldmlInput_IsKnownBasis( fieldmlInfo, firstEvaluator ) ) THEN
      FieldmlInput_IsTemplateCompatible = .FALSE.
      RETURN
    ENDIF

    !At the moment, the code does not support different evaluators per element.

    DO i = 2, count
      evaluator = Fieldml_GetEvaluator( fieldmlInfo%fmlHandle, componentHandle, i )
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

  !<Determines whether or not the given field evaluator can be parsed as an OpenCMISS field. 
  SUBROUTINE FieldmlInput_CheckFieldCompatible( fieldmlInfo, fieldHandle, elementType, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: fieldHandle !<The field evaluator handle.
    INTEGER(C_INT), INTENT(IN) :: elementType !<The element ensemble type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: type, count, i, evaluator, defaultEvaluator

    CALL ENTERS( "FieldmlInput_CheckFieldCompatible", err, errorString, *999 )

    type = Fieldml_GetObjectType( fieldmlInfo%fmlHandle, fieldHandle )

    IF( type /= FHT_AGGREGATE_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Field evaluator must be an aggregate evaluator.", err, errorString, *999 )
    ENDIF

    count = Fieldml_GetEvaluatorCount( fieldmlInfo%fmlHandle, fieldHandle )
    defaultEvaluator = Fieldml_GetDefaultEvaluator( fieldmlInfo%fmlHandle, fieldHandle )

    IF( ( defaultEvaluator /= FML_INVALID_HANDLE ) .AND. .NOT. &
      & FieldmlInput_IsTemplateCompatible( fieldmlInfo, defaultEvaluator, elementType ) ) THEN
      CALL FLAG_ERROR( "Field evaluator must be use a compatible default.", err, errorString, *999 )
      RETURN
    ENDIF

    IF( count == 0 ) THEN
      IF( defaultEvaluator == FML_INVALID_HANDLE ) THEN
        CALL FLAG_ERROR( "Field evaluator must be able to evaluator all field components.", err, errorString, *999 )
      ENDIF
      RETURN
    ENDIF

    DO i = 1, count
      evaluator = Fieldml_GetEvaluator( fieldmlInfo%fmlHandle, fieldHandle, i )
      IF( .NOT. FieldmlInput_IsTemplateCompatible( fieldmlInfo, evaluator, elementType ) ) THEN
        CALL FLAG_ERROR( "Field evaluator must use a compatible component evaluator.", err, errorString, *999 )
        RETURN
      ENDIF
    ENDDO

    CALL EXITS( "FieldmlInput_CheckFieldCompatible" )
    RETURN
999 CALL ERRORS( "FieldmlInput_CheckFieldCompatible", err, errorString )
    CALL EXITS( "FieldmlInput_CheckFieldCompatible" )
    RETURN 1

  END SUBROUTINE FieldmlInput_CheckFieldCompatible

  !
  !================================================================================================================================
  !

  !<Creates an OpenCMISS coordinate system using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FieldmlInput_CoordinateSystemCreateStart( fieldmlInfo, evaluatorName, coordinateSystem, userNumber, &
    & err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: evaluatorName !<The name of the coordinate system evaluator.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER, INTENT(IN) :: coordinateSystem !<The OpenCMISS coordinate system to initialize.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to assign to the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: evaluatorHandle
    INTEGER(C_INT) :: typeHandle, length
    CHARACTER(LEN=MAXSTRLEN) :: name
    INTEGER(INTG) :: coordinateType
    INTEGER(INTG) :: coordinateCount

    CALL ENTERS( "FieldmlInput_CoordinateSystemCreateStart", err, errorString, *999 )

    coordinateType = 0 !There doesn't seem to be a COORDINATE_UNKNOWN_TYPE

    evaluatorHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(evaluatorName) )
    CALL FieldmlUtil_CheckError( "Cannot get coordinate evaluator for geometric field", fieldmlInfo%fmlHandle, &
      & err, errorString, *999 )

    typeHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, evaluatorHandle )
    CALL FieldmlUtil_CheckError( "Cannot get value type for geometric field", fieldmlInfo%fmlHandle, &
      & err, errorString, *999 )

    length = Fieldml_CopyObjectDeclaredName( fieldmlInfo%fmlHandle, typeHandle, name, MAXSTRLEN )

    IF( INDEX( name, 'coordinates.rc.3d' ) == 1 ) THEN
      coordinateType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      coordinateCount = 3
    ELSE IF( INDEX( name, 'coordinates.rc.2d' ) == 1 ) THEN
      coordinateType = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      coordinateCount = 2
    ELSE
      CALL FLAG_ERROR( "Coordinate system not yet supported", err, errorString, *999 )
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
  
  !<Creates an OpenCMISS nodes object using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FieldmlInput_NodesCreateStart( fieldmlInfo, nodesArgumentName, region, nodes, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: nodesArgumentName !<The argument evaluator used as the node index in relevant evaluators.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region !<The region in which to create the nodes.
    TYPE(NODES_TYPE), POINTER, INTENT(INOUT) :: nodes !<The OpenCMISS nodes object to create.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(C_INT) :: nodesArgumentHandle, nodesHandle, nodeCount
    
    nodesArgumentHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(nodesArgumentName) )
    IF( nodesArgumentHandle == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Nodes argument name is invalid", err, errorString, *999 )
    END IF
    
    nodesHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, nodesArgumentHandle )
    IF( nodesHandle == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Nodes argument type is invalid", err, errorString, *999 )
    END IF

    fieldmlInfo%nodesArgumentHandle = nodesArgumentHandle
    fieldmlInfo%nodesHandle = nodesHandle

    nodeCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%nodesHandle )
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

  !<Creates an OpenCMISS mesh using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FieldmlInput_MeshCreateStart( fieldmlInfo, meshArgumentName, mesh, meshNumber, region, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: meshArgumentName !<The argument evaluator used as the mesh location in relevant evaluators.
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT) :: mesh !<The OpenCMISS mesh object to create.
    INTEGER(INTG), INTENT(IN) :: meshNumber !<The user number to assign to the mesh.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region !<The region in which to create the mesh.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(INTG) :: count
    INTEGER(C_INT) :: meshArgument, xiDimensions, elementCount

    CALL ENTERS( "FieldmlInput_MeshCreateStart", err, errorString, *999 )
    
    meshArgument = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(meshArgumentName) )
    IF( meshArgument == FML_INVALID_HANDLE ) THEN
      CALL FieldmlUtil_CheckError( "Named mesh argument not found", fieldmlInfo, err, errorString, *999 )
    ENDIF

    fieldmlInfo%meshHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, meshArgument )
    IF( fieldmlInfo%meshHandle == FML_INVALID_HANDLE ) THEN
      CALL FieldmlUtil_CheckError( "Invalid mesh argument", fieldmlInfo, err, errorString, *999 )
    ENDIF
    
    fieldmlInfo%elementsHandle = Fieldml_GetMeshElementsType( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle )
    fieldmlInfo%elementsArgumentHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(meshArgumentName//".element"))

    fieldmlInfo%xiHandle = Fieldml_GetMeshChartType( fieldmlInfo%fmlHandle, fieldmlInfo%meshHandle )
    fieldmlInfo%xiArgumentHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(meshArgumentName//".xi") )

    count = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, fieldmlInfo%xiHandle )
    IF( ( count < 1 ) .OR. ( count > 3 ) ) THEN
      CALL FLAG_ERROR( "Mesh dimension cannot be greater than 3, or less than 1", err, errorString, *999 )
    ENDIF

    xiDimensions = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, fieldmlInfo%xiHandle )
    elementCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
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
  
  !<Creates an OpenCMISS basis object using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FieldmlInput_BasisCreateStart( fieldmlInfo, evaluatorName, userNumber, basis, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: evaluatorName !<The name of the basis evaluator.
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to assign to the basis.
    TYPE(BASIS_TYPE), POINTER, INTENT(INOUT) :: basis !<The OpenCMISS basis object to create.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.

    !Locals
    INTEGER(INTG) :: listIndex
    INTEGER(C_INT) :: handle, connectivityHandle, layoutHandle, fmlErr
    INTEGER(INTG) :: basisType
    INTEGER(INTG), ALLOCATABLE :: basisInterpolations(:)
    INTEGER(INTG), ALLOCATABLE :: collapse(:)
    
    CALL ENTERS( "FieldmlInput_BasisCreateStart", err, errorString, *999 )

    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(evaluatorName) )
    CALL FieldmlUtil_CheckError( "Named basis not found", fieldmlInfo, err, errorString, *999 )
    CALL LIST_ITEM_IN_LIST_C_INT( fieldmlInfo%basisHandles, handle, listIndex, err, errorString, *999 )
    IF( listIndex /= 0 ) THEN
      CALL FLAG_ERROR( "Named basis already created", err, errorString, *999 )
    ENDIF
    
    CALL FieldmlInput_GetBasisInfo( fieldmlInfo, handle, connectivityHandle, layoutHandle, basisType, basisInterpolations, &
      & collapse, err, errorString, *999 )
    
    CALL LIST_ITEM_ADD_C_INT( fieldmlInfo%basisHandles, handle, err, errorString, *999 )
    CALL LIST_ITEM_ADD_C_INT( fieldmlInfo%basisConnectivityHandles, connectivityHandle, err, errorString, *999 )
    CALL LIST_ITEM_ADD_C_INT( fieldmlInfo%basisLayoutHandles, layoutHandle, err, errorString, *999 )
    fmlErr = Fieldml_SetObjectInt( fieldmlInfo%fmlHandle, handle, userNumber )
    CALL FieldmlUtil_CheckError( "Cannot set basis object's user number", fieldmlInfo, err, errorString, *999 )
  
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

  !<Initialize the given FieldML parsing state from the given FieldML file.
  SUBROUTINE FieldmlInput_InitialiseFromFile( fieldmlInfo, filename, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: filename !<The name of the FieldML file to parse.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: length, count, i, fmlErr
    CHARACTER(LEN=MAXSTRLEN) :: name
    
    CALL ENTERS( "FieldmlInput_InitialiseFromFile", err, errorString, *999 )

    CALL FieldmlUtil_InitialiseInfo( fieldmlInfo, err, errorString, *999 )

    fieldmlInfo%fmlHandle = Fieldml_CreateFromFile( cchar(filename) )
    
    fmlErr = Fieldml_GetLastError( fieldmlInfo%fmlHandle )
    IF( fmlErr /= FML_ERR_NO_ERROR ) THEN
      count = Fieldml_GetErrorCount( fieldmlInfo%fmlHandle )
      DO i = 1,count
        length = Fieldml_CopyError( fieldmlInfo%fmlHandle, i, name, MAXSTRLEN )
        !TODO Use diagnostic output routines.
        WRITE(*,*) "FieldML parse error: "//name(1:length)
      ENDDO
      CALL FLAG_ERROR( "Cannot create FieldML handle from file", err, errorString, *999 )
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
  
  !<Reads an ensemble ordering using the given data source.
  SUBROUTINE FieldmlInput_ReadOrder( fieldmlInfo, orderHandle, order, count, err, errorString, * )
    !Argument
    TYPE(FieldmlInfoType), INTENT(IN) :: fieldmlInfo !<The FieldML parsing state.
    INTEGER(C_INT), INTENT(IN) :: orderHandle !<The data source containing the ordering.
    INTEGER(C_INT), ALLOCATABLE, INTENT(INOUT) :: order(:) !<The array in which the order is stored.
    INTEGER(C_INT), INTENT(IN) :: count !<The number of entries in the ordering.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: readerHandle, i, readCount
    INTEGER(C_INT), TARGET :: buffer(1)
    
    CALL ENTERS( "FieldmlInput_ReadOrder", err, errorString, *999 )

    IF( orderHandle == FML_INVALID_HANDLE ) THEN
      !This is permitted, and indeed common.
      RETURN
    ENDIF
    
    readerHandle = Fieldml_OpenReader( fieldmlInfo%fmlHandle, orderHandle )
    CALL FieldmlUtil_CheckError( "Cannot open order reader", fieldmlInfo, err, errorString, *999 )
    
    ALLOCATE( order(count), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not order array.", err, errorString, *999 )
    DO i = 1, count
      readCount = Fieldml_ReadIntValues( fieldmlInfo%fmlHandle, readerHandle, C_LOC(buffer), 1 )
      IF( readCount /= 1 ) THEN
        DEALLOCATE( order )
        CALL FLAG_ERROR( "Cannot open order reader", err, errorString, *999 )
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

  !<Reorder the given values according to the given ordering.
  SUBROUTINE FieldInput_Reorder( inputBuffer, order, count, outputBuffer, err, errorString, * )
    !Argument
    INTEGER(C_INT), INTENT(IN) :: inputBuffer(:) !<The values to reorder.
    INTEGER(C_INT), ALLOCATABLE, INTENT(IN) :: order(:) !<The ordering to apply.
    INTEGER(C_INT), INTENT(IN) :: count !<The number of values to reorder.
    INTEGER(C_INT), INTENT(INOUT) :: outputBuffer(:) !<The reordered values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
  
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

  !<Creates an OpenCMISS mesh component using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FieldmlInput_CreateMeshComponent( fieldmlInfo, mesh, componentNumber, evaluatorName, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: mesh !<The OpenCMISS mesh in which to create the component.
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number to create.
    TYPE(VARYING_STRING), INTENT(IN) :: evaluatorName !<The name of the mesh component evaluator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: handle, basisReferenceHandle, connectivityHandle, layoutHandle, basisNumber, lastBasisHandle, readCount
    INTEGER(C_INT), ALLOCATABLE, TARGET :: nodesBuffer(:), rawBuffer(:)
    INTEGER(INTG) :: componentCount, elementCount, knownBasisCount, maxBasisNodesCount, basisNodesCount
    INTEGER(INTG) :: elementNumber, knownBasisNumber, count
    INTEGER(C_INT), ALLOCATABLE :: connectivityReaders(:), connectivityCounts(:)
    TYPE(INTEGER_CINT_ALLOC_TYPE), ALLOCATABLE :: connectivityOrders(:)
    INTEGER(C_INT) :: tPtr, dataSource, orderHandle, tmpBasisHandle, fmlErr
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: meshElements

    CALL ENTERS( "FieldmlInput_CreateMeshComponent", err, errorString, *999 )
    
    NULLIFY( basis )
    NULLIFY( meshElements )    

    handle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(evaluatorName) )
    IF( .NOT. FieldmlInput_IsTemplateCompatible( fieldmlInfo, handle, fieldmlInfo%elementsHandle ) ) THEN
      CALL FLAG_ERROR( "Mesh component cannot be created from this evaluator", err, errorString, *999 )
    ENDIF
    
    CALL LIST_NUMBER_OF_ITEMS_GET( fieldmlInfo%componentHandles, count, err, errorString, *999 )
    IF( count < componentNumber ) THEN
      DO componentCount = count + 1, componentNumber
        CALL LIST_ITEM_ADD_C_INT( fieldmlInfo%componentHandles, FML_INVALID_HANDLE, err, errorString, *999 )
      ENDDO
    ENDIF
    
    CALL LIST_ITEM_SET_C_INT( fieldmlInfo%componentHandles, componentNumber, handle, err, errorString, *999 )
    
    CALL LIST_NUMBER_OF_ITEMS_GET( fieldmlInfo%basisHandles, knownBasisCount, err, errorString, *999 )
    ALLOCATE( connectivityReaders( knownBasisCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity readers.", err, errorString, *999 )
    ALLOCATE( connectivityCounts( knownBasisCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity counts.", err, errorString, *999 )
    ALLOCATE( connectivityOrders( knownBasisCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity orders.", err, errorString, *999 )
    
    maxBasisNodesCount = 0
    DO knownBasisNumber = 1, knownBasisCount
      CALL LIST_ITEM_GET_C_INT( fieldmlInfo%basisLayoutHandles, knownBasisNumber, layoutHandle, err, errorString, *999 )
      CALL LIST_ITEM_GET_C_INT( fieldmlInfo%basisConnectivityHandles, knownBasisNumber, connectivityHandle, &
        & err, errorString, *999 )
        
      basisNodesCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, layoutHandle )
      CALL FieldmlUtil_CheckError( "Cannot get local node count for layout", fieldmlInfo, err, errorString, *999 )

      IF( basisNodesCount > maxBasisNodesCount ) THEN
        maxBasisNodesCount = basisNodesCount
      ENDIF
      
      orderHandle = Fieldml_GetSemidenseIndexOrder( fieldmlInfo%fmlHandle, connectivityHandle, 1 )
      CALL FieldmlInput_ReadOrder( fieldmlInfo, orderHandle, connectivityOrders( knownBasisNumber )%ARRAY, &
        & basisNodesCount, err, errorString, *999 )
    
      dataSource = Fieldml_GetDataSource( fieldmlInfo%fmlHandle, connectivityHandle )
      connectivityReaders(knownBasisNumber) = Fieldml_OpenReader( fieldmlInfo%fmlHandle, dataSource )
      connectivityCounts(knownBasisNumber) = basisNodesCount
      CALL FieldmlUtil_CheckError( "Cannot open connectivity reader", fieldmlInfo, err, errorString, *999 )
      
    END DO

    ALLOCATE( nodesBuffer( maxBasisNodesCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate nodes buffer.", err, errorString, *999 )
    ALLOCATE( rawBuffer( maxBasisNodesCount ), STAT = err )
    IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate raw nodes buffer.", err, errorString, *999 )

    elementCount = Fieldml_GetMemberCount( fieldmlInfo%fmlHandle, fieldmlInfo%elementsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get element count for mesh", fieldmlInfo, err, errorString, *999 )
    
    lastBasisHandle = FML_INVALID_HANDLE
    
    DO elementNumber = 1, elementCount
      basisReferenceHandle = Fieldml_GetElementEvaluator( fieldmlInfo%fmlHandle, handle, elementNumber, 1 )
      CALL FieldmlUtil_CheckError( "Cannot get element evaluator from mesh component", fieldmlInfo, err, errorString, *999 )
      
      IF( basisReferenceHandle /= lastBasisHandle ) THEN
        basisNumber = Fieldml_GetObjectInt( fieldmlInfo%fmlHandle, basisReferenceHandle )
        CALL FieldmlUtil_CheckError( "Cannot get basis user number for element evaluator", fieldmlInfo, err, errorString, *999 )
        CALL BASIS_USER_NUMBER_FIND( basisNumber, basis, err, errorString, *999 )
        IF( .NOT. ASSOCIATED( basis ) ) THEN
          CALL FLAG_ERROR( "Basis not found.", err, errorString, *999 ) 
        ENDIF
        lastBasisHandle = basisReferenceHandle
      ENDIF

      IF( elementNumber == 1 ) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START( mesh, componentNumber, basis, meshElements, err, errorString, *999 )
      ENDIF
      
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET( elementNumber, meshElements, basis, err, errorString, *999 )
      
      DO knownBasisNumber = 1, knownBasisCount
        basisNodesCount = connectivityCounts( knownBasisNumber )
        !BUGFIX Intel compiler will explode if we don't use a temporary variable
        tPtr = connectivityReaders(knownBasisNumber)
        readCount = Fieldml_ReadIntValues( fieldmlInfo%fmlHandle, tPtr, C_LOC(rawBuffer), basisNodesCount )
        IF( readCount /= basisNodesCount ) THEN
          CALL FLAG_ERROR( "Error reading connectivity", err, errorString, *999 )
        ENDIF
        CALL LIST_ITEM_GET_C_INT( fieldmlInfo%basisHandles, knownBasisNumber, tmpBasisHandle, err, errorString, *999 )
        IF( tmpBasisHandle == basisReferenceHandle ) THEN
          CALL FieldInput_Reorder( rawBuffer, connectivityOrders(knownBasisNumber)%ARRAY, basisNodesCount, &
            & nodesBuffer, err, errorString, *999 )
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET( elementNumber, meshElements, nodesBuffer(1:basisNodesCount), &
            & err, errorString, *999 )
        ENDIF
      ENDDO
  
    END DO
    
    DO knownBasisNumber = 1, knownBasisCount
      !BUGFIX Intel compiler will explode if we don't use a temporary variable
      tPtr = connectivityReaders(knownBasisNumber)
      fmlErr = Fieldml_CloseReader( fieldmlInfo%fmlHandle, tPtr )
      IF( fmlErr /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( "Error closing connectivity reader", err, errorString, *999 )
      ENDIF
      IF( ALLOCATED( connectivityOrders( knownBasisNumber )%ARRAY ) ) THEN
        DEALLOCATE( connectivityOrders( knownBasisNumber )%ARRAY )
      ENDIF
    ENDDO
    
    DEALLOCATE( nodesBuffer )
    DEALLOCATE( connectivityReaders )
    DEALLOCATE( connectivityCounts )
    DEALLOCATE( connectivityOrders )
    
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH( meshElements, err, errorString, *999 )
    
    fmlErr = Fieldml_SetObjectInt( fieldmlInfo%fmlHandle, handle, componentNumber )

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
        IF( ALLOCATED( connectivityOrders( knownBasisNumber )%ARRAY ) ) THEN
          DEALLOCATE( connectivityOrders( knownBasisNumber )%ARRAY )
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
  
  !<Creates an OpenCMISS field using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FieldmlInput_FieldCreateStart( fieldmlInfo, region, decomposition, fieldNumber, field, evaluatorName, &
    & err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: region !<The region in which to create the field.
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: decomposition !<The decomposition to use when creating the field.
    INTEGER(INTG), INTENT(IN) :: fieldNumber !<The user number to assign to the created field.
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: field !<The OpenCMISS field object to create.
    TYPE(VARYING_STRING), INTENT(IN) :: evaluatorName !<The name of the field evaluator.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    INTEGER(C_INT) :: fieldHandle, templateHandle, typeHandle
    INTEGER(INTG) :: componentNumber, templateComponentNumber, fieldDimensions

    CALL ENTERS( "FieldmlInput_FieldCreateStart", err, errorString, *999 )

    fieldHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(evaluatorName) )
    CALL FieldmlUtil_CheckError( "Cannot get named field evaluator", fieldmlInfo, err, errorString, *999 )
    typeHandle = Fieldml_GetValueType( fieldmlInfo%fmlHandle, fieldHandle )
    CALL FieldmlUtil_CheckError( "Cannot get named field evaluator's value type", fieldmlInfo, err, errorString, *999 )
    fieldDimensions = Fieldml_GetTypeComponentCount( fieldmlInfo%fmlHandle, typeHandle )
    CALL FieldmlUtil_CheckError( "Cannot get named field evaluator's component count", fieldmlInfo, err, errorString, *999 )
    
    CALL FieldmlInput_CheckFieldCompatible( fieldmlInfo, fieldHandle, fieldmlInfo%elementsHandle, err, errorString, *999 )

    NULLIFY( field )
    CALL FIELD_CREATE_START( fieldNumber, region, field, err, errorString, *999 )
    CALL FIELD_TYPE_SET( field, FIELD_GEOMETRIC_TYPE, err, errorString, *999 )
    CALL FIELD_MESH_DECOMPOSITION_SET( field, decomposition, err, errorString, *999 )
    CALL FIELD_SCALING_TYPE_SET( field, FIELD_NO_SCALING, err, errorString, *999 )

    DO componentNumber = 1, fieldDimensions
      templateHandle = Fieldml_GetElementEvaluator( fieldmlInfo%fmlHandle, fieldHandle, componentNumber, 1 )
      CALL FieldmlUtil_CheckError( "Cannot get field component evaluator", fieldmlInfo, err, errorString, *999 )

      templateComponentNumber = Fieldml_GetObjectInt( fieldmlInfo%fmlHandle, templateHandle )
      CALL FieldmlUtil_CheckError( "Cannot get field component mesh component number", fieldmlInfo, err, errorString, *999 )

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

  !<Update the given field's nodal parameters using the given parameter evaluator.
  SUBROUTINE FieldmlInput_FieldNodalParametersUpdate( fieldmlInfo, evaluatorName, field, err, errorString, * )
    !Arguments
    TYPE(FieldmlInfoType), INTENT(INOUT) :: fieldmlInfo !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: evaluatorName !<The name of the nodal dofs evaluator.
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: field !<The field whose parameters are to be updated.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: errorString !<The error string.
    
    !Locals
    TYPE(MESH_TYPE), POINTER :: mesh
    TYPE(NODES_TYPE), POINTER :: nodes
    INTEGER(C_INT) :: nodalDofsHandle, count, dataSource, fmlErr
    INTEGER(INTG) :: versionNumber,componentNumber, nodeNumber, fieldDimensions, meshNodeCount, gNode
    INTEGER(INTG) :: meshComponent
    LOGICAL :: nodeExists
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: buffer(:)
    INTEGER(C_INT) :: reader
    
    mesh => field%DECOMPOSITION%MESH

    nodalDofsHandle = Fieldml_GetObjectByName( fieldmlInfo%fmlHandle, cchar(evaluatorName) )
    CALL FieldmlUtil_CheckError( "Cannot get nodal field dofs", fieldmlInfo, err, errorString, *999 )
  
    dataSource = Fieldml_GetDataSource( fieldmlInfo%fmlHandle, nodalDofsHandle )
    CALL FieldmlUtil_CheckError( "Cannot get nodal data source", fieldmlInfo, err, errorString, *999 )

    reader = Fieldml_OpenReader( fieldmlInfo%fmlHandle, dataSource )
    CALL FieldmlUtil_CheckError( "Cannot open nodal dofs reader", fieldmlInfo, err, errorString, *999 )
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( field, FIELD_U_VARIABLE_TYPE, fieldDimensions, err, errorString, *999 )
    
    IF( reader /= FML_INVALID_HANDLE ) THEN
      ALLOCATE( buffer( fieldDimensions ), STAT = err )
      IF( err /= 0 ) CALL FLAG_ERROR( "Could not allocate raw nodes buffer.", err, errorString, *999 )
      
      !TODO Code assumes that the data is dense in both node and component indexes.
      NULLIFY( nodes )
      CALL REGION_NODES_GET( mesh%REGION, nodes, err, errorString, *999 )
      CALL NODES_NUMBER_OF_NODES_GET( nodes, meshNodeCount, err, errorString, *999 )
      CALL FieldmlUtil_CheckError( "Cannot get mesh nodes count", fieldmlInfo, err, errorString, *999 )
      DO nodeNumber = 1, meshNodeCount
        count = Fieldml_ReadDoubleValues( fieldmlInfo%fmlHandle, reader, C_LOC(buffer), fieldDimensions )
        IF( count /= fieldDimensions ) THEN
          CALL FLAG_ERROR( "Cannot read nodal dofs for field components", err, errorString, *999 )
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
  
      fmlErr = Fieldml_CloseReader( fieldmlInfo%fmlHandle, reader )
      IF( fmlErr /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( "Error closing nodal dofs reader", err, errorString, *999 )
      ENDIF

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
