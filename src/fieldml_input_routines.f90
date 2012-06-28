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

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS
  USE CONSTANTS
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE FIELD_ROUTINES
  USE FIELDML_API
  USE FIELDML_TYPES
  USE FIELDML_UTIL_ROUTINES
  USE INPUT_OUTPUT
  USE LISTS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE REGION_ROUTINES
  USE STRINGS

  IMPLICIT NONE

  PRIVATE

  !Interfaces

  INTERFACE

  END INTERFACE

  PUBLIC :: FIELDML_INPUT_INITIALISE_FROM_FILE, FIELDML_INPUT_MESH_CREATE_START, &
    & FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START, FIELDML_INPUT_BASIS_CREATE_START, FIELDML_INPUT_CREATE_MESH_COMPONENT, &
    & FIELDML_INPUT_FIELD_CREATE_START, FIELDML_INPUT_FIELD_PARAMETERS_UPDATE, FIELDML_INPUT_NODES_CREATE_START

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    CALL ENTERS( "FIELDML_ASSERT_IS_IN", ERR, ERROR, *999 )

    IF( FIELDML_INFO%IS_OUT ) THEN
      CALL FLAG_ERROR( "Outbound FieldML handle used four an input-only operation.", ERR, ERROR, *999 )
    ENDIF
    
    CALL EXITS( "FIELDML_ASSERT_IS_IN" )
    RETURN
999 CALL ERRORS( "FIELDML_ASSERT_IS_IN", ERR, ERROR )
    CALL EXITS( "FIELDML_ASSERT_IS_IN" )
    RETURN 1
    
  END SUBROUTINE FIELDML_ASSERT_IS_IN
  
  !
  !================================================================================================================================
  !
  
  !>Determines the connectivity evaluator and layout argument for the given basis.
  SUBROUTINE FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO( FIELDML_INFO, BASIS_HANDLE, PARAM_ARG_HANDLE, CONNECTIVITY_HANDLE, &
    & LAYOUT_HANDLE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: BASIS_HANDLE !<The basis handle.
    INTEGER(INTG), INTENT(IN) :: PARAM_ARG_HANDLE !<The basis parameters argument handle.
    INTEGER(INTG), INTENT(OUT) :: CONNECTIVITY_HANDLE !<The basis connectivity evaluator handle.
    INTEGER(INTG), INTENT(OUT) :: LAYOUT_HANDLE !<The local node layout.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Local variables
    INTEGER(INTG) :: COUNT, BIND_NUMBER, PARAMS_HANDLE, ARG_HANDLE, LAYOUT_INDEX_HANDLE
    
    CALL ENTERS( "FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO", ERR, ERROR, *999 )

    COUNT = Fieldml_GetBindCount( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE )
    IF( COUNT /= 2 ) THEN
      CALL FLAG_ERROR( "Library basis evaluators must have exactly two binds.", ERR, ERROR, *999 )
    END IF
    
    PARAMS_HANDLE = FML_INVALID_HANDLE
    DO BIND_NUMBER = 1, COUNT
      ARG_HANDLE = Fieldml_GetBindArgument( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE, BIND_NUMBER )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get bind for interpolator.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      IF( ARG_HANDLE == PARAM_ARG_HANDLE ) THEN
        PARAMS_HANDLE = Fieldml_GetBindEvaluator( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE, BIND_NUMBER )
      ENDIF
    ENDDO

    IF( PARAMS_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Library interpolators must have a correct parameter bind.", ERR, ERROR, *999 )
    ENDIF

    IF( Fieldml_GetObjectType( FIELDML_INFO%FML_HANDLE, PARAMS_HANDLE ) /= FHT_AGGREGATE_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Parameter evaluator for interpolator must be an aggregate.", ERR, ERROR, *999 )
    ENDIF
    
    COUNT = Fieldml_GetBindCount( FIELDML_INFO%FML_HANDLE, PARAMS_HANDLE )
    IF( COUNT /= 1 ) THEN
      CALL FLAG_ERROR( "Nodal parameter evaluator must only have one bind.", ERR, ERROR, *999 )
    ENDIF

    IF( Fieldml_GetBindArgument( FIELDML_INFO%FML_HANDLE, PARAMS_HANDLE, 1 ) /= FIELDML_INFO%NODES_ARGUMENT_HANDLE ) THEN
      CALL FLAG_ERROR( "Nodal parameter evaluator must bind the nodes argument.", ERR, ERROR, *999 )
    ENDIF
    
    CONNECTIVITY_HANDLE = Fieldml_GetBindEvaluator( FIELDML_INFO%FML_HANDLE, PARAMS_HANDLE, 1 )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get connectivity source for nodal parameters.", FIELDML_INFO%FML_HANDLE, &
      & ERR, ERROR, *999 )
      
    LAYOUT_INDEX_HANDLE = Fieldml_GetIndexEvaluator( FIELDML_INFO%FML_HANDLE, PARAMS_HANDLE, 1 )
    LAYOUT_HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, LAYOUT_INDEX_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get connectivity source for nodal parameters.", FIELDML_INFO%FML_HANDLE, &
      & ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO

  !
  !================================================================================================================================
  !
  
  !>Determine the basis collapse parameters from the given evaluator's name.
  SUBROUTINE FIELDML_INPUT_GET_BASIS_COLLAPSE( NAME, COLLAPSE, ERR, ERROR, * )
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: NAME !<The basis evaluator name.
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: COLLAPSE(:) !<The array of OpenCMISS basis collapse constants.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    CALL ENTERS( "FIELDML_INPUT_GET_BASIS_COLLAPSE", ERR, ERROR, *999 )

    COLLAPSE = BASIS_NOT_COLLAPSED
    
    IF( SIZE( COLLAPSE ) > 0 ) THEN
      IF( INDEX( NAME, "_xi1C" ) /= 0 ) THEN
        COLLAPSE(1) = BASIS_XI_COLLAPSED
      ELSE IF( INDEX( NAME, "_xi10" ) /= 0 ) THEN
        COLLAPSE(1) = BASIS_COLLAPSED_AT_XI0
      ELSE IF( INDEX( NAME, "_xi11" ) /= 0 ) THEN
        COLLAPSE(1) = BASIS_COLLAPSED_AT_XI1
      ENDIF
    ENDIF
  
    IF( SIZE( COLLAPSE ) > 1 ) THEN
      IF( INDEX( NAME, "_xi2C" ) /= 0 ) THEN
        COLLAPSE(2) = BASIS_XI_COLLAPSED
      ELSE IF( INDEX( NAME, "_xi20" ) /= 0 ) THEN
        COLLAPSE(2) = BASIS_COLLAPSED_AT_XI0
      ELSE IF( INDEX( NAME, "_xi21" ) /= 0 ) THEN
        COLLAPSE(2) = BASIS_COLLAPSED_AT_XI1
      ENDIF
    ENDIF
  
    IF( SIZE( COLLAPSE ) > 2 ) THEN
      IF( INDEX( NAME, "_xi3C" ) /= 0 ) THEN
        COLLAPSE(3) = BASIS_XI_COLLAPSED
      ELSE IF( INDEX( NAME, "_xi30" ) /= 0 ) THEN
        COLLAPSE(3) = BASIS_COLLAPSED_AT_XI0
      ELSE IF( INDEX( NAME, "_xi31" ) /= 0 ) THEN
        COLLAPSE(3) = BASIS_COLLAPSED_AT_XI1
      ENDIF
    ENDIF

    CALL EXITS( "FIELDML_INPUT_GET_BASIS_COLLAPSE" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_GET_BASIS_COLLAPSE", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_GET_BASIS_COLLAPSE" )
    RETURN 1
  
  END SUBROUTINE FIELDML_INPUT_GET_BASIS_COLLAPSE

  !
  !================================================================================================================================
  !

  !>Determines the basis configuration from the given basis evaluator.
  SUBROUTINE FIELDML_INPUT_GET_BASIS_INFO( FIELDML_INFO, BASIS_HANDLE, CONNECTIVITY_HANDLE, LAYOUT_HANDLE, BASISTYPE, &
    & BASIS_INTERPOLATIONS, COLLAPSE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: BASIS_HANDLE !<The basis evaluator handle.
    INTEGER(INTG), INTENT(OUT) :: CONNECTIVITY_HANDLE !<The basis connectivity evaluator handle.
    INTEGER(INTG), INTENT(OUT) :: LAYOUT_HANDLE !<The local node layout.
    INTEGER(INTG), INTENT(OUT) :: BASISTYPE !<The OpenCMISS basis type.
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: BASIS_INTERPOLATIONS(:) !<The per-xi basis interpolations (for TP bases).
    INTEGER(INTG), ALLOCATABLE, INTENT(OUT) :: COLLAPSE(:) !<The collapse constants for the basis.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: LENGTH, LIBRARY_BASIS_HANDLE, PARAM_ARG_HANDLE
    CHARACTER(LEN=MAXSTRLEN) :: NAME
    TYPE(VARYING_STRING) :: COLLAPSE_NAME
    
    CALL ENTERS( "FIELDML_INPUT_GET_BASIS_INFO", ERR, ERROR, *999 )

    IF( .NOT. FIELDML_INPUT_IS_KNOWN_BASIS( FIELDML_INFO, BASIS_HANDLE, ERR, ERROR ) ) THEN
      CALL FLAG_ERROR( "Basis specified in FieldML file is not yet supported.", ERR, ERROR, *999 )
    ENDIF
    IF(ERR/=0) GOTO 999

    IF( Fieldml_GetObjectType( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE ) /= FHT_REFERENCE_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Basis evaluator must be a continuous reference.", ERR, ERROR, *999 )
    ENDIF
    
    LIBRARY_BASIS_HANDLE = Fieldml_GetReferenceSourceEvaluator( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Basis specified in FieldML is not a reference evaluator.", FIELDML_INFO%FML_HANDLE, &
      & ERR, ERROR, *999 )
    LENGTH = Fieldml_CopyObjectDeclaredName( FIELDML_INFO%FML_HANDLE, LIBRARY_BASIS_HANDLE, NAME, MAXSTRLEN )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get name of basis evaluator.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    IF( INDEX( NAME, 'interpolator.3d.unit.triquadraticLagrange') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.3d.unit.triquadraticLagrange.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(3), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE(3), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
      BASISTYPE = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.3d.unit.trilinearLagrange') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.3d.unit.trilinearLagrange.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(3), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE(3), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_LINEAR_LAGRANGE_INTERPOLATION
      BASISTYPE = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.2d.unit.biquadraticLagrange') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.2d.unit.biquadraticLagrange.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(2), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE(2), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
      BASISTYPE = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.2d.unit.bilinearLagrange') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.2d.unit.bilinearLagrange.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(2), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE(2), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_LINEAR_LAGRANGE_INTERPOLATION
      BASISTYPE = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.1d.unit.linearLagrange') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.1d.unit.linearLagrange.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(1), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE(1), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_LINEAR_LAGRANGE_INTERPOLATION
      BASISTYPE = BASIS_LAGRANGE_HERMITE_TP_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.2d.unit.bilinearSimplex') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.2d.unit.bilinearSimplex.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(2), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_LINEAR_SIMPLEX_INTERPOLATION
      BASISTYPE = BASIS_SIMPLEX_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.2d.unit.biquadraticSimplex') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.2d.unit.biquadraticSimplex.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(2), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_QUADRATIC_SIMPLEX_INTERPOLATION
      BASISTYPE = BASIS_SIMPLEX_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.3d.unit.trilinearSimplex') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.3d.unit.trilinearSimplex.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(3), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_LINEAR_SIMPLEX_INTERPOLATION
      BASISTYPE = BASIS_SIMPLEX_TYPE
    ELSE IF( INDEX( NAME, 'interpolator.3d.unit.triquadraticSimplex') == 1 ) THEN
      PARAM_ARG_HANDLE = Fieldml_GetObjectByDeclaredName( FIELDML_INFO%FML_HANDLE, &
        & "parameters.3d.unit.triquadraticSimplex.argument"//C_NULL_CHAR )
      ALLOCATE( BASIS_INTERPOLATIONS(3), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate interpolation array.", ERR, ERROR, *999 )
      BASIS_INTERPOLATIONS = BASIS_QUADRATIC_SIMPLEX_INTERPOLATION
      BASISTYPE = BASIS_SIMPLEX_TYPE
    ELSE
      CALL FLAG_ERROR( "Basis "//NAME(1:LENGTH)//" cannot yet be interpreted.", ERR, ERROR, *999 )
    ENDIF
    
    IF( BASISTYPE == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      COLLAPSE_NAME = NAME(1:LENGTH)
      CALL FIELDML_INPUT_GET_BASIS_COLLAPSE( COLLAPSE_NAME, COLLAPSE, ERR, ERROR, *999 )
    ENDIF
    
    CALL FIELDML_INPUT_GET_BASIS_CONNECTIVITY_INFO( FIELDML_INFO, BASIS_HANDLE, PARAM_ARG_HANDLE, CONNECTIVITY_HANDLE, &
      & LAYOUT_HANDLE, ERR, ERROR, *999 )

    CALL ENTERS( "FIELDML_INPUT_GET_BASIS_INFO", ERR, ERROR, *999 )
    CALL EXITS( "FIELDML_INPUT_GET_BASIS_INFO" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_GET_BASIS_INFO", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_GET_BASIS_INFO" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_GET_BASIS_INFO

  !
  !================================================================================================================================
  !

  !>Determines whether or not the given basis evaluator is known to OpenCMISS.
  FUNCTION FIELDML_INPUT_IS_KNOWN_BASIS( FIELDML_INFO, BASIS_HANDLE, ERR, ERROR )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: BASIS_HANDLE !<The basis handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Function
    LOGICAL :: FIELDML_INPUT_IS_KNOWN_BASIS

    !Locals
    INTEGER(INTG) :: LENGTH, LIBRARY_BASIS_HANDLE
    CHARACTER(LEN=MAXSTRLEN) :: NAME
    
    CALL ENTERS( "FIELDML_INPUT_IS_KNOWN_BASIS", ERR, ERROR, *999 )

    IF( Fieldml_GetObjectType( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE ) /= FHT_REFERENCE_EVALUATOR ) THEN
      FIELDML_INPUT_IS_KNOWN_BASIS = .FALSE.
      CALL EXITS( "FIELDML_INPUT_IS_KNOWN_BASIS" )
      RETURN
    ENDIF

    LIBRARY_BASIS_HANDLE = Fieldml_GetReferenceSourceEvaluator( FIELDML_INFO%FML_HANDLE, BASIS_HANDLE )
    LENGTH = Fieldml_CopyObjectDeclaredName( FIELDML_INFO%FML_HANDLE, LIBRARY_BASIS_HANDLE, NAME, MAXSTRLEN )

    IF( ( INDEX( NAME, 'interpolator.3d.unit.triquadraticLagrange') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.1d.unit.linearLagrange') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.2d.unit.biquadraticLagrange') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.2d.unit.bilinearLagrange') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.3d.unit.trilinearLagrange') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.2d.unit.bilinearSimplex') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.2d.unit.biquadraticSimplex') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.3d.unit.trilinearSimplex') /= 1 ) .AND. &
      & ( INDEX( NAME, 'interpolator.3d.unit.triquadraticSimplex') /= 1 ) ) THEN
      FIELDML_INPUT_IS_KNOWN_BASIS = .FALSE.
    ELSE
      FIELDML_INPUT_IS_KNOWN_BASIS = .TRUE.
    ENDIF

    CALL EXITS( "FIELDML_INPUT_IS_KNOWN_BASIS" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_IS_KNOWN_BASIS", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_IS_KNOWN_BASIS" )
    
  END FUNCTION FIELDML_INPUT_IS_KNOWN_BASIS
  
  !
  !================================================================================================================================
  !

  !>Determines whether or not the given evaluator is a recognisable mesh component evaluator.
  FUNCTION FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE( FIELDML_INFO, COMPONENT_HANDLE, ELEMENT_TYPE, ERR, ERROR )
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_HANDLE !<The mesh component evaluator.
    INTEGER(INTG), INTENT(IN) :: ELEMENT_TYPE !<The element ensemble type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    LOGICAL :: FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE

    INTEGER(INTG) :: OBJECT_TYPE, COUNT, I, EVALUATOR, TYPE, FIRST_EVALUATOR, EVALUATOR_HANDLE, DEFAULT_EVALUATOR

    CALL ENTERS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE", ERR, ERROR, *999 )
  
    OBJECT_TYPE = Fieldml_GetObjectType( FIELDML_INFO%FML_HANDLE, COMPONENT_HANDLE )
    IF( OBJECT_TYPE /= FHT_PIECEWISE_EVALUATOR ) THEN
      FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .FALSE.
      CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
      RETURN
    ENDIF

    EVALUATOR_HANDLE = Fieldml_GetIndexEvaluator( FIELDML_INFO%FML_HANDLE, COMPONENT_HANDLE, 1 )
    TYPE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, EVALUATOR_HANDLE )
    IF( TYPE /= ELEMENT_TYPE ) THEN
      FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .TRUE.
      CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
      RETURN
    ENDIF

    COUNT = Fieldml_GetEvaluatorCount( FIELDML_INFO%FML_HANDLE, COMPONENT_HANDLE )
    DEFAULT_EVALUATOR = Fieldml_GetDefaultEvaluator( FIELDML_INFO%FML_HANDLE, COMPONENT_HANDLE )
    
    IF( ( DEFAULT_EVALUATOR /= FML_INVALID_HANDLE ) .AND. .NOT. &
      & FIELDML_INPUT_IS_KNOWN_BASIS( FIELDML_INFO, DEFAULT_EVALUATOR, ERR, ERROR ) ) THEN
      FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .FALSE.
      CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
      RETURN
    ENDIF
    IF(ERR/=0) GOTO 999

    IF( COUNT == 0 ) THEN
      IF( DEFAULT_EVALUATOR == FML_INVALID_HANDLE ) THEN
          FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .FALSE.
      ELSE
          FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .TRUE.
      ENDIF
      CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
      RETURN
    ENDIF

    FIRST_EVALUATOR = Fieldml_GetEvaluator( FIELDML_INFO%FML_HANDLE, COMPONENT_HANDLE, 1 )
    IF( .NOT. FIELDML_INPUT_IS_KNOWN_BASIS( FIELDML_INFO, FIRST_EVALUATOR, ERR, ERROR ) ) THEN
      FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .FALSE.
      CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
      RETURN
    ENDIF
    IF(ERR/=0) GOTO 999

    !At the moment, the code does not support different evaluators per element.

    DO I = 2, COUNT
      EVALUATOR = Fieldml_GetEvaluator( FIELDML_INFO%FML_HANDLE, COMPONENT_HANDLE, I )
      IF( EVALUATOR /= FIRST_EVALUATOR ) THEN
        FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .FALSE.
        CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
        RETURN
      ENDIF
    ENDDO

    FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE = .TRUE.

    CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE" )

  END FUNCTION FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE

  !
  !================================================================================================================================
  !

  !>Determines whether or not the given field evaluator can be parsed as an OpenCMISS field. 
  SUBROUTINE FIELDML_INPUT_CHECK_FIELD_COMPATIBLE( FIELDML_INFO, FIELD_HANDLE, ELEMENT_TYPE, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: FIELD_HANDLE !<The field evaluator handle.
    INTEGER(INTG), INTENT(IN) :: ELEMENT_TYPE !<The element ensemble type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: TYPE, COUNT, I, EVALUATOR, DEFAULT_EVALUATOR

    CALL ENTERS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE", ERR, ERROR, *999 )

    TYPE = Fieldml_GetObjectType( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE )

    IF( TYPE /= FHT_AGGREGATE_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Field evaluator must be an aggregate evaluator.", ERR, ERROR, *999 )
    ENDIF

    COUNT = Fieldml_GetEvaluatorCount( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE )
    DEFAULT_EVALUATOR = Fieldml_GetDefaultEvaluator( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE )

    IF( ( DEFAULT_EVALUATOR /= FML_INVALID_HANDLE ) .AND. .NOT. &
      & FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE( FIELDML_INFO, DEFAULT_EVALUATOR, ELEMENT_TYPE, ERR, ERROR ) ) THEN
      CALL FLAG_ERROR( "Field evaluator must be use a compatible default.", ERR, ERROR, *999 )
      CALL EXITS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE" )
      RETURN
    ENDIF
    IF(ERR/=0) GOTO 999

    IF( COUNT == 0 ) THEN
      IF( DEFAULT_EVALUATOR == FML_INVALID_HANDLE ) THEN
        CALL FLAG_ERROR( "Field evaluator must be able to evaluator all field components.", ERR, ERROR, *999 )
      ENDIF
      CALL EXITS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE" )
      RETURN
    ENDIF

    DO I = 1, COUNT
      EVALUATOR = Fieldml_GetEvaluator( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE, I )
      IF( .NOT. FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE( FIELDML_INFO, EVALUATOR, ELEMENT_TYPE, ERR, ERROR ) ) THEN
        CALL FLAG_ERROR( "Field evaluator must use a compatible component evaluator.", ERR, ERROR, *999 )
        CALL EXITS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE" )
        RETURN
      ENDIF
      IF(ERR/=0) GOTO 999
    ENDDO

    CALL EXITS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_CHECK_FIELD_COMPATIBLE" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_CHECK_FIELD_COMPATIBLE

  !
  !================================================================================================================================
  !

  !>Creates an OpenCMISS coordinate system using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START( FIELDML_INFO, EVALUATOR_NAME, COORDINATE_SYSTEM, USER_NUMBER, &
    & ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: EVALUATOR_NAME !<The name of the coordinate system evaluator.
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER, INTENT(IN) :: COORDINATE_SYSTEM !<The OpenCMISS coordinate system to initialize.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to assign to the coordinate system.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: EVALUATOR_HANDLE
    INTEGER(INTG) :: TYPE_HANDLE, LENGTH
    CHARACTER(LEN=MAXSTRLEN) :: NAME
    INTEGER(INTG) :: COORDINATE_TYPE
    INTEGER(INTG) :: COORDINATE_COUNT

    CALL ENTERS( "FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )

    COORDINATE_TYPE = 0 !There doesn't seem to be a COORDINATE_UNKNOWN_TYPE

    EVALUATOR_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(EVALUATOR_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get coordinate evaluator for geometric field "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    TYPE_HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, EVALUATOR_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get value type for geometric field "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    LENGTH = Fieldml_CopyObjectDeclaredName( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE, NAME, MAXSTRLEN )

    IF( INDEX( NAME, 'coordinates.rc.3d' ) == 1 ) THEN
      COORDINATE_TYPE = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      COORDINATE_COUNT = 3
    ELSE IF( INDEX( NAME, 'coordinates.rc.2d' ) == 1 ) THEN
      COORDINATE_TYPE = COORDINATE_RECTANGULAR_CARTESIAN_TYPE
      COORDINATE_COUNT = 2
    ELSE
      CALL FLAG_ERROR( "Coordinate system "//NAME(1:LENGTH)//" not yet supported.", ERR, ERROR, *999 )
    ENDIF

    CALL COORDINATE_SYSTEM_CREATE_START( USER_NUMBER, COORDINATE_SYSTEM, ERR, ERROR, *999 )
    !Set the coordinate system dimension and type
    CALL COORDINATE_SYSTEM_DIMENSION_SET( COORDINATE_SYSTEM, COORDINATE_COUNT, ERR, ERROR, *999 )
    CALL COORDINATE_SYSTEM_TYPE_SET( COORDINATE_SYSTEM, COORDINATE_TYPE, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_COORDINATE_SYSTEM_CREATE_START


  !
  !================================================================================================================================
  !
  
  !>Creates an OpenCMISS nodes object using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FIELDML_INPUT_NODES_CREATE_START( FIELDML_INFO, NODES_ARGUMENT_NAME, REGION, NODES, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: NODES_ARGUMENT_NAME !<The argument evaluator used as the node index in relevant evaluators.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: REGION !<The region in which to create the nodes.
    TYPE(NODES_TYPE), POINTER, INTENT(INOUT) :: NODES !<The OpenCMISS nodes object to create.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: NODES_ARGUMENT_HANDLE, NODES_HANDLE, NODE_COUNT
    
    CALL ENTERS( "FIELDML_INPUT_NODES_CREATE_START", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )

    NODES_ARGUMENT_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(NODES_ARGUMENT_NAME) )
    IF( NODES_ARGUMENT_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Nodes argument name "//NODES_ARGUMENT_NAME//" is invalid.", ERR, ERROR, *999 )
    END IF

    IF( Fieldml_GetObjectType( FIELDML_INFO%FML_HANDLE, NODES_ARGUMENT_HANDLE ) /= FHT_ARGUMENT_EVALUATOR ) THEN
      CALL FLAG_ERROR( "Nodes argument "//NODES_ARGUMENT_NAME//" type is not an argument evaluator.", ERR, ERROR, *999 )
    ENDIF

    NODES_HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, NODES_ARGUMENT_HANDLE )
    IF( NODES_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Nodes argument "//NODES_ARGUMENT_NAME//" type is invalid.", ERR, ERROR, *999 )
    END IF

    FIELDML_INFO%NODES_ARGUMENT_HANDLE = NODES_ARGUMENT_HANDLE
    FIELDML_INFO%NODES_HANDLE = NODES_HANDLE

    NODE_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%NODES_HANDLE )
    NULLIFY( NODES )
    CALL NODES_CREATE_START( REGION, NODE_COUNT, NODES, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_INPUT_NODES_CREATE_START" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_NODES_CREATE_START", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_NODES_CREATE_START" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_NODES_CREATE_START


  !
  !================================================================================================================================
  !

  !>Creates an OpenCMISS mesh using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FIELDML_INPUT_MESH_CREATE_START( FIELDML_INFO, MESH_ARGUMENT_NAME, MESH, MESH_NUMBER, REGION, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: MESH_ARGUMENT_NAME !<The argument evaluator used as the mesh location in relevant evaluators.
    TYPE(MESH_TYPE), POINTER, INTENT(INOUT) :: MESH !<The OpenCMISS mesh object to create.
    INTEGER(INTG), INTENT(IN) :: MESH_NUMBER !<The user number to assign to the mesh.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: REGION !<The region in which to create the mesh.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: COUNT
    INTEGER(INTG) :: MESH_ARGUMENT, XI_DIMENSIONS, ELEMENT_COUNT

    CALL ENTERS( "FIELDML_INPUT_MESH_CREATE_START", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )
    
    MESH_ARGUMENT = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(MESH_ARGUMENT_NAME) )
    IF( MESH_ARGUMENT == FML_INVALID_HANDLE ) THEN
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Named mesh argument "//MESH_ARGUMENT_NAME//" not found.", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF

    FIELDML_INFO%MESH_HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, MESH_ARGUMENT )
    IF( FIELDML_INFO%MESH_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Invalid mesh argument type for "//MESH_ARGUMENT_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF
    
    FIELDML_INFO%ELEMENTS_HANDLE = Fieldml_GetMeshElementsType( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%MESH_HANDLE )
    FIELDML_INFO%ELEMENTS_ARGUMENT_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, &
      & cchar(MESH_ARGUMENT_NAME//".element"))

    FIELDML_INFO%XI_HANDLE = Fieldml_GetMeshChartType( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%MESH_HANDLE )
    FIELDML_INFO%XI_ARGUMENT_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(MESH_ARGUMENT_NAME//".xi") )

    COUNT = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%XI_HANDLE )
    IF( ( COUNT < 1 ) .OR. ( COUNT > 3 ) ) THEN
      CALL FLAG_ERROR( "Mesh "//MESH_ARGUMENT_NAME//" dimension cannot be greater than 3, or less than 1.", &
        & ERR, ERROR, *999 )
    ENDIF

    XI_DIMENSIONS = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%XI_HANDLE )
    ELEMENT_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%ELEMENTS_HANDLE )
    NULLIFY( MESH )
    CALL MESH_CREATE_START( MESH_NUMBER, REGION, XI_DIMENSIONS, MESH, ERR, ERROR, *999 )
    CALL MESH_NUMBER_OF_ELEMENTS_SET( MESH, ELEMENT_COUNT, ERR, ERROR, *999 )
    
    CALL EXITS( "FIELDML_INPUT_MESH_CREATE_START" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_MESH_CREATE_START", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_MESH_CREATE_START" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_MESH_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Creates an OpenCMISS basis object using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FIELDML_INPUT_BASIS_CREATE_START( FIELDML_INFO, EVALUATOR_NAME, USER_NUMBER, BASIS, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: EVALUATOR_NAME !<The name of the basis evaluator.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to assign to the basis.
    TYPE(BASIS_TYPE), POINTER, INTENT(INOUT) :: BASIS !<The OpenCMISS basis object to create.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: LIST_INDEX
    INTEGER(INTG) :: HANDLE, CONNECTIVITY_HANDLE, LAYOUT_HANDLE, FML_ERR
    INTEGER(INTG) :: BASISTYPE
    INTEGER(INTG), ALLOCATABLE :: BASIS_INTERPOLATIONS(:)
    INTEGER(INTG), ALLOCATABLE :: COLLAPSE(:)
    
    CALL ENTERS( "FIELDML_INPUT_BASIS_CREATE_START", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )

    HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(EVALUATOR_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot find basis evaluator "//EVALUATOR_NAME//".", FIELDML_INFO%FML_HANDLE, &
      & ERR, ERROR, *999 )
    CALL LIST_ITEM_IN_LIST( FIELDML_INFO%BASIS_HANDLES, HANDLE, LIST_INDEX, ERR, ERROR, *999 )
    IF( LIST_INDEX /= 0 ) THEN
      CALL FLAG_ERROR( "Named basis "//EVALUATOR_NAME//" already created", ERR, ERROR, *999 )
    ENDIF
    
    CALL FIELDML_INPUT_GET_BASIS_INFO( FIELDML_INFO, HANDLE, CONNECTIVITY_HANDLE, LAYOUT_HANDLE, BASISTYPE, &
      & BASIS_INTERPOLATIONS, COLLAPSE, ERR, ERROR, *999 )
    
    CALL LIST_ITEM_ADD( FIELDML_INFO%BASIS_HANDLES, HANDLE, ERR, ERROR, *999 )
    CALL LIST_ITEM_ADD( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, CONNECTIVITY_HANDLE, ERR, ERROR, *999 )
    CALL LIST_ITEM_ADD( FIELDML_INFO%BASIS_LAYOUT_HANDLES, LAYOUT_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetObjectInt( FIELDML_INFO%FML_HANDLE, HANDLE, USER_NUMBER )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set user number for basis "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
  
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START( USER_NUMBER, BASIS, ERR, ERROR, *999 )
    CALL BASIS_TYPE_SET( BASIS, BASISTYPE, ERR, ERROR, *999 )
    CALL BASIS_NUMBER_OF_XI_SET( BASIS, size( BASIS_INTERPOLATIONS ), ERR, ERROR, *999 )
    CALL BASIS_INTERPOLATION_XI_SET( BASIS, BASIS_INTERPOLATIONS, ERR, ERROR, *999 )
    !Note: collapse bases currently only supported for BASIS_LAGRANGE_HERMITE_TP_TYPE
    IF( size( BASIS_INTERPOLATIONS ) > 1 .AND. ALLOCATED(COLLAPSE)) THEN
      CALL BASIS_COLLAPSED_XI_SET( BASIS, COLLAPSE, ERR, ERROR, *999 )
    ENDIF
    
    IF( ALLOCATED( BASIS_INTERPOLATIONS ) ) THEN
      DEALLOCATE( BASIS_INTERPOLATIONS )
    ENDIF
    IF( ALLOCATED( COLLAPSE ) ) THEN
      DEALLOCATE( COLLAPSE )
    ENDIF
    
    CALL EXITS( "FIELDML_INPUT_BASIS_CREATE_START" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_BASIS_CREATE_START", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_BASIS_CREATE_START" )
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !>Initialize the given FieldML parsing state from the given FieldML file.
  SUBROUTINE FIELDML_INPUT_INITIALISE_FROM_FILE( FIELDML_INFO, FILENAME, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: FILENAME !<The name of the FieldML file to parse.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: LENGTH, COUNT, I, FML_ERR
    CHARACTER(LEN=MAXSTRLEN) :: NAME
    
    CALL ENTERS( "FIELDML_INPUT_INITIALISE_FROM_FILE", ERR, ERROR, *999 )

    CALL FIELDML_IO_INITIALISE( FIELDML_INFO, .FALSE., ERR, ERROR, *999 )

    FIELDML_INFO%FML_HANDLE = Fieldml_CreateFromFile( cchar(FILENAME) )
    
    FML_ERR = Fieldml_GetLastError( FIELDML_INFO%FML_HANDLE )
    IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
      COUNT = Fieldml_GetErrorCount( FIELDML_INFO%FML_HANDLE )
      DO I = 1,COUNT
        LENGTH = Fieldml_CopyError( FIELDML_INFO%FML_HANDLE, I, NAME, MAXSTRLEN )
         CALL WRITE_STRING_VALUE(ERROR_OUTPUT_TYPE,"FieldML parse error: ",NAME(1:LENGTH),ERR,ERROR,*999)
      ENDDO
      CALL FLAG_ERROR( "Cannot create FieldML handle from file "//FILENAME//".", ERR, ERROR, *999 )
    ENDIF

    CALL EXITS( "FIELDML_INPUT_INITIALISE_FROM_FILE" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_INITIALISE_FROM_FILE", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_INITIALISE_FROM_FILE" )
    RETURN 1
    
  END SUBROUTINE FIELDML_INPUT_INITIALISE_FROM_FILE

  !
  !================================================================================================================================
  !
  
  !>Reads an ensemble ordering using the given data source.
  SUBROUTINE FIELDML_INPUT_READ_ORDER( FIELDML_INFO, ORDER_HANDLE, ORDER, COUNT, ERR, ERROR, * )
    !Argument
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: ORDER_HANDLE !<The data source containing the ordering.
    INTEGER(INTG), ALLOCATABLE, TARGET, INTENT(INOUT) :: ORDER(:) !<The array in which the order is stored.
    INTEGER(INTG), INTENT(IN) :: COUNT !<The number of entries in the ordering.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: READER_HANDLE, RANK, FML_ERR
    INTEGER(INTG), TARGET :: OFFSETS(1), SIZES(1)
    
    CALL ENTERS( "FIELDML_INPUT_READ_ORDER", ERR, ERROR, *999 )

    IF( ORDER_HANDLE == FML_INVALID_HANDLE ) THEN
      !This is permitted, and indeed common.
      CALL EXITS( "FIELDML_INPUT_READ_ORDER" )
      RETURN
    ENDIF
    
    RANK = Fieldml_GetArrayDataSourceRank( FIELDML_INFO%FML_HANDLE, ORDER_HANDLE )
    IF( RANK /= 1 ) THEN
      CALL FLAG_ERROR( "Invalid rank for ensemble order.", err, ERROR, *999 )
    ENDIF
    
    READER_HANDLE = Fieldml_OpenReader( FIELDML_INFO%FML_HANDLE, ORDER_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open order reader.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    ALLOCATE( ORDER(COUNT), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate order array.", ERR, ERROR, *999 )
    OFFSETS(:) = 0
    SIZES(1) = COUNT

    FML_ERR = Fieldml_ReadIntSlab( READER_HANDLE, C_LOC(OFFSETS), C_LOC(SIZES), C_LOC(ORDER) )
    IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
      CALL FLAG_ERROR( "Error reading order data"//"("// TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", &
        & ERR, ERROR, *999 )
    ENDIF

    FML_ERR = Fieldml_CloseReader( READER_HANDLE )
    
    CALL EXITS( "FIELDML_INPUT_READ_ORDER" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_READ_ORDER", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_READ_ORDER" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_READ_ORDER

  !
  !================================================================================================================================
  !

  !>Reorder the given values according to the given ordering.
  SUBROUTINE FIELDML_INPUT_REORDER( INPUT_BUFFER, ORDER, COUNT, OUTPUT_BUFFER, ERR, ERROR, * )
    !Argument
    INTEGER(INTG), INTENT(IN) :: INPUT_BUFFER(:) !<The values to reorder.
    INTEGER(INTG), ALLOCATABLE, INTENT(IN) :: ORDER(:) !<The ordering to apply.
    INTEGER(INTG), INTENT(IN) :: COUNT !<The number of values to reorder.
    INTEGER(INTG), INTENT(INOUT) :: OUTPUT_BUFFER(:) !<The reordered values.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
  
    !Locals
    INTEGER(INTG) :: I
    
    CALL ENTERS( "FIELDML_INPUT_REORDER", ERR, ERROR, *999 )
    
    IF( ALLOCATED( ORDER ) ) THEN
      DO I = 1,COUNT
        OUTPUT_BUFFER( I ) = INPUT_BUFFER( ORDER( I ) )
      ENDDO
    ELSE
      OUTPUT_BUFFER = INPUT_BUFFER
    ENDIF

    CALL EXITS( "FIELDML_INPUT_REORDER" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_REORDER", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_REORDER" )
    RETURN 1
  
  END SUBROUTINE FIELDML_INPUT_REORDER

  !
  !================================================================================================================================
  !

  !>Creates an OpenCMISS mesh component using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FIELDML_INPUT_CREATE_MESH_COMPONENT( FIELDML_INFO, MESH, COMPONENT_NUMBER, EVALUATOR_NAME, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: MESH !<The OpenCMISS mesh in which to create the component.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to create.
    TYPE(VARYING_STRING), INTENT(IN) :: EVALUATOR_NAME !<The name of the mesh component evaluator.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: HANDLE, BASIS_REFERENCE_HANDLE, CONNECTIVITY_HANDLE, LAYOUT_HANDLE, BASIS_NUMBER, LAST_BASIS_HANDLE
    INTEGER(INTG), ALLOCATABLE, TARGET :: NODES_BUFFER(:), RAW_BUFFER(:)
    INTEGER(INTG) :: COMPONENT_COUNT, ELEMENT_COUNT, KNOWN_BASIS_COUNT, MAX_BASIS_NODES_COUNT, BASIS_NODES_COUNT
    INTEGER(INTG) :: ELEMENT_NUMBER, KNOWN_BASIS_NUMBER, COUNT
    INTEGER(INTG), TARGET :: OFFSETS(2), SIZES(2)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_READERS(:), CONNECTIVITY_COUNTS(:)
    TYPE(INTEGER_CINT_ALLOC_TYPE), ALLOCATABLE :: CONNECTIVITY_ORDERS(:)
    INTEGER(INTG) :: TEMP_POINTER, DATA_SOURCE, ORDER_HANDLE, TEMP_BASIS_HANDLE, FML_ERR
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS

    CALL ENTERS( "FIELDML_INPUT_CREATE_MESH_COMPONENT", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )
    
    NULLIFY( BASIS )
    NULLIFY( MESH_ELEMENTS )    

    HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(EVALUATOR_NAME) )
    IF( .NOT. FIELDML_INPUT_IS_TEMPLATE_COMPATIBLE( FIELDML_INFO, HANDLE, FIELDML_INFO%ELEMENTS_HANDLE, ERR, ERROR ) ) THEN
      CALL FLAG_ERROR( "Mesh component cannot be created from evaluator "//EVALUATOR_NAME//".", ERR, ERROR, *999 )
    ENDIF
    IF(ERR/=0) GOTO 999
    
    CALL LIST_NUMBER_OF_ITEMS_GET( FIELDML_INFO%COMPONENT_HANDLES, COUNT, ERR, ERROR, *999 )
    IF( COUNT < COMPONENT_NUMBER ) THEN
      DO COMPONENT_COUNT = COUNT + 1, COMPONENT_NUMBER
        CALL LIST_ITEM_ADD( FIELDML_INFO%COMPONENT_HANDLES, FML_INVALID_HANDLE, ERR, ERROR, *999 )
      ENDDO
    ENDIF
    
    CALL LIST_ITEM_SET( FIELDML_INFO%COMPONENT_HANDLES, COMPONENT_NUMBER, HANDLE, ERR, ERROR, *999 )
    
    CALL LIST_NUMBER_OF_ITEMS_GET( FIELDML_INFO%BASIS_HANDLES, KNOWN_BASIS_COUNT, ERR, ERROR, *999 )
    ALLOCATE( CONNECTIVITY_READERS( KNOWN_BASIS_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity readers for "//EVALUATOR_NAME//".", ERR, ERROR, *999 )
    ALLOCATE( CONNECTIVITY_COUNTS( KNOWN_BASIS_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity counts for "//EVALUATOR_NAME//".", ERR, ERROR, *999 )
    ALLOCATE( CONNECTIVITY_ORDERS( KNOWN_BASIS_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity orders for "//EVALUATOR_NAME//".", ERR, ERROR, *999 )
    
    MAX_BASIS_NODES_COUNT = 0
    DO KNOWN_BASIS_NUMBER = 1, KNOWN_BASIS_COUNT
      CALL LIST_ITEM_GET( FIELDML_INFO%BASIS_LAYOUT_HANDLES, KNOWN_BASIS_NUMBER, LAYOUT_HANDLE, ERR, ERROR, *999 )
      CALL LIST_ITEM_GET( FIELDML_INFO%BASIS_CONNECTIVITY_HANDLES, KNOWN_BASIS_NUMBER, CONNECTIVITY_HANDLE, &
        & ERR, ERROR, *999 )
        
      BASIS_NODES_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, LAYOUT_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get local node count for layout for mesh component "//EVALUATOR_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      IF( BASIS_NODES_COUNT > MAX_BASIS_NODES_COUNT ) THEN
        MAX_BASIS_NODES_COUNT = BASIS_NODES_COUNT
      ENDIF
      
      ORDER_HANDLE = Fieldml_GetParameterIndexOrder( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_HANDLE, 1 )
      CALL FIELDML_INPUT_READ_ORDER( FIELDML_INFO, ORDER_HANDLE, CONNECTIVITY_ORDERS( KNOWN_BASIS_NUMBER )%ARRAY, &
        & BASIS_NODES_COUNT, ERR, ERROR, *999 )
    
      DATA_SOURCE = Fieldml_GetDataSource( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_HANDLE )
      CONNECTIVITY_READERS(KNOWN_BASIS_NUMBER) = Fieldml_OpenReader( FIELDML_INFO%FML_HANDLE, DATA_SOURCE )
      CONNECTIVITY_COUNTS(KNOWN_BASIS_NUMBER) = BASIS_NODES_COUNT
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open connectivity reader for mesh component "//EVALUATOR_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
    END DO

    ALLOCATE( NODES_BUFFER( MAX_BASIS_NODES_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate nodes buffer for "//EVALUATOR_NAME//".", ERR, ERROR, *999 )
    ALLOCATE( RAW_BUFFER( MAX_BASIS_NODES_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate raw nodes buffer for "//EVALUATOR_NAME//".", ERR, ERROR, *999 )

    ELEMENT_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%ELEMENTS_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get element count for mesh with component "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    LAST_BASIS_HANDLE = FML_INVALID_HANDLE
    
    OFFSETS(:) = 0
    SIZES(1) = 1
    SIZES(2) = 0
    
    DO ELEMENT_NUMBER = 1, ELEMENT_COUNT
      BASIS_REFERENCE_HANDLE = Fieldml_GetElementEvaluator( FIELDML_INFO%FML_HANDLE, HANDLE, ELEMENT_NUMBER, 1 )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get element evaluator from mesh component "//EVALUATOR_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
      IF( BASIS_REFERENCE_HANDLE /= LAST_BASIS_HANDLE ) THEN
        BASIS_NUMBER = Fieldml_GetObjectInt( FIELDML_INFO%FML_HANDLE, BASIS_REFERENCE_HANDLE )
        CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get basis user number for element evaluator for mesh component "//&
          & EVALUATOR_NAME//".", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
        CALL BASIS_USER_NUMBER_FIND( BASIS_NUMBER, BASIS, ERR, ERROR, *999 )
        IF( .NOT. ASSOCIATED( BASIS ) ) THEN
          CALL FLAG_ERROR( "Basis not found for component "//EVALUATOR_NAME//".", ERR, ERROR, *999 ) 
        ENDIF
        LAST_BASIS_HANDLE = BASIS_REFERENCE_HANDLE
      ENDIF

      IF( ELEMENT_NUMBER == 1 ) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START( MESH, COMPONENT_NUMBER, BASIS, MESH_ELEMENTS, ERR, ERROR, *999 )
      ENDIF
      
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET( ELEMENT_NUMBER, MESH_ELEMENTS, BASIS, ERR, ERROR, *999 )
      
      DO KNOWN_BASIS_NUMBER = 1, KNOWN_BASIS_COUNT
        BASIS_NODES_COUNT = CONNECTIVITY_COUNTS( KNOWN_BASIS_NUMBER )
        !BUGFIX Intel compiler will explode if we don't use a temporary variable
        TEMP_POINTER = CONNECTIVITY_READERS(KNOWN_BASIS_NUMBER)
        SIZES(2) = BASIS_NODES_COUNT
        FML_ERR = Fieldml_ReadIntSlab( TEMP_POINTER, &
          & C_LOC(OFFSETS), C_LOC(SIZES), C_LOC(RAW_BUFFER) )
        IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
          CALL FLAG_ERROR( "Error reading connectivity for "//EVALUATOR_NAME//"("// &
            & TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", ERR, ERROR, *999 )
        ENDIF
        CALL LIST_ITEM_GET( FIELDML_INFO%BASIS_HANDLES, KNOWN_BASIS_NUMBER, TEMP_BASIS_HANDLE, ERR, ERROR, *999 )
        IF( TEMP_BASIS_HANDLE == BASIS_REFERENCE_HANDLE ) THEN
          CALL FIELDML_INPUT_REORDER( RAW_BUFFER, CONNECTIVITY_ORDERS(KNOWN_BASIS_NUMBER)%ARRAY, BASIS_NODES_COUNT, &
            & NODES_BUFFER, ERR, ERROR, *999 )
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET( ELEMENT_NUMBER, MESH_ELEMENTS, NODES_BUFFER(1:BASIS_NODES_COUNT), &
            & ERR, ERROR, *999 )
        ENDIF
      ENDDO
      
      OFFSETS(1) = OFFSETS(1) + 1
  
    END DO
    
    DO KNOWN_BASIS_NUMBER = 1, KNOWN_BASIS_COUNT
      !BUGFIX Intel compiler will explode if we don't use a temporary variable
      TEMP_POINTER = CONNECTIVITY_READERS(KNOWN_BASIS_NUMBER)
      FML_ERR = Fieldml_CloseReader( TEMP_POINTER )
      IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( "Error closing connectivity reader for "//EVALUATOR_NAME//"("// &
          & TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", ERR, ERROR, *999 )
      ENDIF
      IF( ALLOCATED( CONNECTIVITY_ORDERS( KNOWN_BASIS_NUMBER )%ARRAY ) ) THEN
        DEALLOCATE( CONNECTIVITY_ORDERS( KNOWN_BASIS_NUMBER )%ARRAY )
      ENDIF
    ENDDO
    
    DEALLOCATE( NODES_BUFFER )
    DEALLOCATE( CONNECTIVITY_READERS )
    DEALLOCATE( CONNECTIVITY_COUNTS )
    DEALLOCATE( CONNECTIVITY_ORDERS )
    
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH( MESH_ELEMENTS, ERR, ERROR, *999 )
    
    FML_ERR = Fieldml_SetObjectInt( FIELDML_INFO%FML_HANDLE, HANDLE, COMPONENT_NUMBER )

    CALL EXITS( "FIELDML_INPUT_CREATE_MESH_COMPONENT" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_CREATE_MESH_COMPONENT", ERR, ERROR )
    IF( ALLOCATED( NODES_BUFFER ) ) THEN
      DEALLOCATE( NODES_BUFFER )
    ENDIF
    IF( ALLOCATED( CONNECTIVITY_READERS ) ) THEN
      DEALLOCATE( CONNECTIVITY_READERS )
    ENDIF
    IF( ALLOCATED( CONNECTIVITY_COUNTS ) ) THEN
      DEALLOCATE( CONNECTIVITY_COUNTS )
    ENDIF
    IF( ALLOCATED( CONNECTIVITY_ORDERS ) ) THEN
      DO KNOWN_BASIS_NUMBER = 1, KNOWN_BASIS_COUNT
        IF( ALLOCATED( CONNECTIVITY_ORDERS( KNOWN_BASIS_NUMBER )%ARRAY ) ) THEN
          DEALLOCATE( CONNECTIVITY_ORDERS( KNOWN_BASIS_NUMBER )%ARRAY )
        ENDIF
      ENDDO
    
      DEALLOCATE( CONNECTIVITY_ORDERS )
    ENDIF
    
    CALL EXITS( "FIELDML_INPUT_CREATE_MESH_COMPONENT" )
    RETURN 1

  END SUBROUTINE FIELDML_INPUT_CREATE_MESH_COMPONENT

  !
  !================================================================================================================================
  !
  
  !>Creates an OpenCMISS field using relevant parameters from FieldML. Does not call CreateFinish.
  SUBROUTINE FIELDML_INPUT_FIELD_CREATE_START( FIELDML_INFO, REGION, DECOMPOSITION, FIELD_NUMBER, FIELD, VARIABLE_TYPE, &
    & EVALUATOR_NAME, ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(REGION_TYPE), POINTER, INTENT(IN) :: REGION !<The region in which to create the field.
    TYPE(DECOMPOSITION_TYPE), POINTER, INTENT(IN) :: DECOMPOSITION !<The decomposition to use when creating the field.
    INTEGER(INTG), INTENT(IN) :: FIELD_NUMBER !<The user number to assign to the created field.
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: FIELD !<The OpenCMISS field object to create.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type.
    TYPE(VARYING_STRING), INTENT(IN) :: EVALUATOR_NAME !<The name of the field evaluator.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: FIELD_HANDLE, TEMPLATE_HANDLE, TYPE_HANDLE
    INTEGER(INTG) :: COMPONENT_NUMBER, TEMPLATE_COMPONENT_NUMBER, FIELD_DIMENSIONS

    CALL ENTERS( "FIELDML_INPUT_FIELD_CREATE_START", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )

    FIELD_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(EVALUATOR_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get named field evaluator "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    TYPE_HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get named field evaluator's value type for "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FIELD_DIMENSIONS = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get named field evaluator's component count for "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    CALL FIELDML_INPUT_CHECK_FIELD_COMPATIBLE( FIELDML_INFO, FIELD_HANDLE, FIELDML_INFO%ELEMENTS_HANDLE, ERR, ERROR, *999 )

    NULLIFY( FIELD )
    CALL FIELD_CREATE_START( FIELD_NUMBER, REGION, FIELD, ERR, ERROR, *999 )
    CALL FIELD_TYPE_SET( FIELD, FIELD_GEOMETRIC_TYPE, ERR, ERROR, *999 )
    CALL FIELD_MESH_DECOMPOSITION_SET( FIELD, DECOMPOSITION, ERR, ERROR, *999 )
    CALL FIELD_SCALING_TYPE_SET( FIELD, FIELD_NO_SCALING, ERR, ERROR, *999 )

    DO COMPONENT_NUMBER = 1, FIELD_DIMENSIONS
      TEMPLATE_HANDLE = Fieldml_GetElementEvaluator( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE, COMPONENT_NUMBER, 1 )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( var_str("Cannot get field component ")//COMPONENT_NUMBER//" evaluator for "//&
        & EVALUATOR_NAME//".", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      TEMPLATE_COMPONENT_NUMBER = Fieldml_GetObjectInt( FIELDML_INFO%FML_HANDLE, TEMPLATE_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( var_str("Cannot get mesh component number for field component ")//COMPONENT_NUMBER//&
        & " of "//EVALUATOR_NAME//".", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      CALL FIELD_COMPONENT_MESH_COMPONENT_SET( FIELD, VARIABLE_TYPE, COMPONENT_NUMBER, TEMPLATE_COMPONENT_NUMBER, &
        & ERR, ERROR, *999 )
    ENDDO

    CALL EXITS( "FIELDML_INPUT_FIELD_CREATE_START" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_FIELD_CREATE_START", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_FIELD_CREATE_START" )
    RETURN 1
  
  END SUBROUTINE FIELDML_INPUT_FIELD_CREATE_START

  !
  !================================================================================================================================
  !

  !>Inputs from a FieldML file the parameters for a field variable parameter set.
  SUBROUTINE FIELDML_INPUT_FIELD_PARAMETERS_UPDATE(FIELDML_INFO, EVALUATOR_NAME, FIELD, VARIABLE_TYPE, SET_TYPE, &
    & ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: EVALUATOR_NAME !<The name of the nodal dofs evaluator.
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: FIELD !<The field whose parameters are to be updated.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    INTEGER(INTG) :: component_idx,INTERPOLATION_TYPE,MESH_COMPONENT1,MESH_COMPONENT2,NUMBER_OF_COMPONENTS
    LOGICAL :: IS_ALL_NODAL_INTERPOLATION,IS_SAME_MESH_COMPONENTS

    CALL ENTERS("FIELDML_INPUT_FIELD_PARAMETERS_UPDATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_NUMBER_OF_COMPONENTS_GET(FIELD,VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      IF(NUMBER_OF_COMPONENTS>0) THEN
        CALL FIELD_COMPONENT_INTERPOLATION_GET(FIELD,VARIABLE_TYPE,1,INTERPOLATION_TYPE,ERR,ERROR,*999)
        CALL FIELD_COMPONENT_MESH_COMPONENT_GET(FIELD,VARIABLE_TYPE,1,MESH_COMPONENT1,ERR,ERROR,*999)
        IS_ALL_NODAL_INTERPOLATION=INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION
        IS_SAME_MESH_COMPONENTS=.TRUE.
        DO component_idx=2,NUMBER_OF_COMPONENTS
          CALL FIELD_COMPONENT_INTERPOLATION_GET(FIELD,VARIABLE_TYPE,component_idx,INTERPOLATION_TYPE,ERR,ERROR,*999)
          CALL FIELD_COMPONENT_MESH_COMPONENT_GET(FIELD,VARIABLE_TYPE,component_idx,MESH_COMPONENT2,ERR,ERROR,*999)
          IS_ALL_NODAL_INTERPOLATION=IS_ALL_NODAL_INTERPOLATION.AND.INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION
          IS_SAME_MESH_COMPONENTS=IS_SAME_MESH_COMPONENTS.AND.MESH_COMPONENT2==MESH_COMPONENT1
        ENDDO !component_idx
        IF(IS_ALL_NODAL_INTERPOLATION) THEN
          IF(IS_SAME_MESH_COMPONENTS) THEN
            CALL FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE(FIELDML_INFO,EVALUATOR_NAME,FIELD,VARIABLE_TYPE,SET_TYPE, &
              & ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR( &
              & "FieldML input parameters only implemented for fields where all components have the same mesh component.", &
              & ERR,ERROR,*999)            
          ENDIF
        ELSE
          CALL FLAG_ERROR("FieldML input parameters only implemented for fields where all components are nodally interpolated.", &
            & ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field does not have any components.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELDML_INPUT_FIELD_PARAMETERS_UPDATE")
    RETURN
999 CALL ERRORS("FIELDML_INPUT_FIELD_PARAMETERS_UPDATE",ERR,ERROR)
    CALL EXITS("FIELDML_INPUT_FIELD_PARAMETERS_UPDATE")
    RETURN 1
  END SUBROUTINE FIELDML_INPUT_FIELD_PARAMETERS_UPDATE

  !
  !================================================================================================================================
  !

  !>Update the given field's nodal parameters using the given parameter evaluator.
  SUBROUTINE FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE( FIELDML_INFO, EVALUATOR_NAME, FIELD, VARIABLE_TYPE, SET_TYPE, &
    & ERR, ERROR, * )
    !Arguments
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: EVALUATOR_NAME !<The name of the nodal dofs evaluator.
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: FIELD !<The field whose parameters are to be updated.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG) :: NODAL_DOFS_HANDLE, DATA_SOURCE, FML_ERR, RANK
    INTEGER(INTG) :: VERSION_NUMBER,COMPONENT_NUMBER, NODE_NUMBER, FIELD_DIMENSIONS, MESH_NODE_COUNT
    INTEGER(INTG), TARGET :: OFFSETS(2), SIZES(2)
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: BUFFER(:)
    INTEGER(INTG) :: READER
    INTEGER(INTG) :: myComputationalNodeNumber,nodeDomain,meshComponentNumber
    
    CALL ENTERS( "FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_IN( FIELDML_INFO, ERR, ERROR, *999 )
    
    MESH => FIELD%DECOMPOSITION%MESH

    NODAL_DOFS_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(EVALUATOR_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get nodal field dofs evaluator "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
  
    DATA_SOURCE = Fieldml_GetDataSource( FIELDML_INFO%FML_HANDLE, NODAL_DOFS_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get nodal data source for "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
    RANK = Fieldml_GetArrayDataSourceRank( FIELDML_INFO%FML_HANDLE, DATA_SOURCE )
    IF( RANK /= 2 ) THEN
      CALL FLAG_ERROR( "Invalid rank for nodal dofs.", err, ERROR, *999 )
    ENDIF

    READER = Fieldml_OpenReader( FIELDML_INFO%FML_HANDLE, DATA_SOURCE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open nodal dofs reader for "//EVALUATOR_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( FIELD, VARIABLE_TYPE, FIELD_DIMENSIONS, ERR, ERROR, *999 )
    
    ALLOCATE( BUFFER( FIELD_DIMENSIONS ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate raw nodes buffer for "//EVALUATOR_NAME//".", ERR, ERROR, *999 )
      
    !TODO Code assumes that the data is dense in both node and component indexes.
    NULLIFY( NODES )
    CALL REGION_NODES_GET( MESH%REGION, NODES, ERR, ERROR, *999 )
    CALL NODES_NUMBER_OF_NODES_GET( NODES, MESH_NODE_COUNT, ERR, ERROR, *999 )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( var_str("Cannot get mesh nodes count for mesh ")//mesh%USER_NUMBER//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    OFFSETS(:) = 0
    SIZES(1) = 1
    SIZES(2) = FIELD_DIMENSIONS
    
    DO NODE_NUMBER = 1, MESH_NODE_COUNT
      FML_ERR = Fieldml_ReadDoubleSlab( READER, C_LOC(OFFSETS), C_LOC(SIZES), C_LOC(BUFFER) )
      OFFSETS(1) = OFFSETS(1) + 1
      IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( "Cannot read nodal dofs from "//EVALUATOR_NAME//"("&
          & // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", ERR, ERROR, *999 )
      ENDIF

      DO COMPONENT_NUMBER = 1, FIELD_DIMENSIONS
        !Default to version 1 of each node derivative (value hardcoded in loop)
        VERSION_NUMBER = 1

        myComputationalNodeNumber = COMPUTATIONAL_NODE_NUMBER_GET(err,error)
        CALL DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(FIELD%DECOMPOSITION,meshComponentNumber,err,error,*999)
        CALL DECOMPOSITION_NODE_DOMAIN_GET(FIELD%DECOMPOSITION,NODE_NUMBER,meshComponentNumber,nodeDomain,err,error,*999)
        IF(nodeDomain==myComputationalNodeNumber) THEN
          CALL FIELD_PARAMETER_SET_UPDATE_NODE( FIELD, VARIABLE_TYPE, SET_TYPE, VERSION_NUMBER, &
            & NO_GLOBAL_DERIV, NODE_NUMBER, COMPONENT_NUMBER, BUFFER( COMPONENT_NUMBER ), ERR, ERROR, *999 )
        ENDIF

      ENDDO
    ENDDO
    
    DEALLOCATE( BUFFER )
  
    FML_ERR = Fieldml_CloseReader( READER )
    IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
      CALL FLAG_ERROR( "Error closing nodal dofs reader for "//EVALUATOR_NAME//"("&
        & // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", ERR, ERROR, *999 )
    ENDIF

    !TODO Set element and constant parameters
    
    CALL EXITS( "FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE" )
    RETURN
999 CALL ERRORS( "FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE", ERR, ERROR )
    CALL EXITS( "FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE" )
    RETURN 1
  
  END SUBROUTINE FIELDML_INPUT_FIELD_NODAL_PARAMETERS_UPDATE

  !
  !================================================================================================================================
  !

END MODULE FIELDML_INPUT_ROUTINES
