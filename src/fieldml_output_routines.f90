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
  USE COORDINATE_ROUTINES
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

  !Interfaces
  TYPE CONNECTIVITY_INFO_TYPE
    INTEGER(INTG) :: CONNECTIVITY_HANDLE !<The basis connectivity evaluator handle.
    INTEGER(INTG) :: LAYOUT_HANDLE !<The local node layout.
  END TYPE CONNECTIVITY_INFO_TYPE

  TYPE BASIS_INFO_TYPE
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The OpenCMISS basis.
    INTEGER(INTG) :: CONNECTIVITY_HANDLE !<The basis connectivity evaluator handle.
    INTEGER(INTG) :: REFERENCE_HANDLE !<The reference evaluator representing the basis.
    INTEGER(INTG) :: LAYOUT_HANDLE !<The local node layout.
  END TYPE BASIS_INFO_TYPE

  INTERFACE FIELDML_OUTPUT_ADD_FIELD
    MODULE PROCEDURE FIELDML_OUTPUT_ADD_FIELD_NO_TYPE
    MODULE PROCEDURE FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE
  END INTERFACE
 
  PUBLIC :: FIELDML_OUTPUT_WRITE, FIELDML_OUTPUT_ADD_FIELD, FIELDML_OUTPUT_INITIALISE_INFO, FIELDML_OUTPUT_IMPORT, &
    & FIELDML_OUTPUT_ADD_FIELD_COMPONENTS

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE FIELDML_ASSERT_IS_OUT( FIELDML_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    CALL ENTERS( "FIELDML_ASSERT_IS_OUT", ERR, ERROR, *999 )

    IF( .NOT. FIELDML_INFO%IS_OUT ) THEN
      CALL FLAG_ERROR( "Inbound FieldML handle used four an output-only operation.", ERR, ERROR, *999 )
    ENDIF
    
    CALL EXITS( "FIELDML_ASSERT_IS_OUT" )
    RETURN
999 CALL ERRORS( "FIELDML_ASSERT_IS_OUT", ERR, ERROR )
    CALL EXITS( "FIELDML_ASSERT_IS_OUT" )
    RETURN 1
    
  END SUBROUTINE FIELDML_ASSERT_IS_OUT
  
  !
  !================================================================================================================================
  !
  
  !>Get the text suffix corresponding to the given array of collapse constants.
  SUBROUTINE FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX( COLLAPSE_INFO, SUFFIX, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: COLLAPSE_INFO(:) !<The collapse into from which to generate the suffix.
    TYPE(VARYING_STRING), INTENT(INOUT) :: SUFFIX !<The suffix string encoding the collapse info.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    INTEGER(INTG) :: I
    
    CALL ENTERS( "FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX", ERR, ERROR, *999 )

    SUFFIX = ""
    DO I = 1, SIZE( COLLAPSE_INFO )
      IF( COLLAPSE_INFO( I ) == BASIS_XI_COLLAPSED ) THEN
        SUFFIX = SUFFIX // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//"C"
      ELSEIF( COLLAPSE_INFO( I ) == BASIS_COLLAPSED_AT_XI0 ) THEN
        SUFFIX = SUFFIX // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//"0"
      ELSEIF( COLLAPSE_INFO( I ) == BASIS_COLLAPSED_AT_XI1 ) THEN
        SUFFIX = SUFFIX // "_xi"//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//"1"
      ENDIF
    ENDDO
  
    CALL EXITS( "FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX
  
  !
  !================================================================================================================================
  !

  !>Import the named object from the built-in library into the current FieldML document. The local name will be the same as the remote name.
  FUNCTION FIELDML_OUTPUT_IMPORT_FML( FML_HANDLE, REMOTE_NAME, ERR, ERROR )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    TYPE(VARYING_STRING), INTENT(IN) :: REMOTE_NAME !<The name of the object to import.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    INTEGER(INTG) :: FIELDML_OUTPUT_IMPORT_FML
    
    !Local variables
    INTEGER(INTG) :: IMPORT_INDEX
    
    CALL ENTERS( "FIELDML_OUTPUT_IMPORT_FML", ERR, ERROR, *999 )

    FIELDML_OUTPUT_IMPORT_FML = Fieldml_GetObjectByName( FML_HANDLE, cchar(REMOTE_NAME) )
    IF( FIELDML_OUTPUT_IMPORT_FML == FML_INVALID_HANDLE ) THEN
      IMPORT_INDEX = Fieldml_AddImportSource( FML_HANDLE, &
        & "http://www.fieldml.org/resources/xml/0.5/FieldML_Library_0.5.xml"//C_NULL_CHAR, "library"//C_NULL_CHAR )
      FIELDML_OUTPUT_IMPORT_FML = Fieldml_AddImport( FML_HANDLE, IMPORT_INDEX, cchar(REMOTE_NAME), cchar(REMOTE_NAME) )
      IF( FIELDML_OUTPUT_IMPORT_FML == FML_INVALID_HANDLE ) ERR = 1
    ENDIF

    CALL EXITS( "FIELDML_OUTPUT_IMPORT_FML" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_IMPORT_FML", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_IMPORT_FML" )

  END FUNCTION FIELDML_OUTPUT_IMPORT_FML
  
  !
  !================================================================================================================================
  !

  !>Import the named object from the built-in library into the current FieldML document. The local name will be the same as the remote name.
  FUNCTION FIELDML_OUTPUT_IMPORT( FIELDML_INFO, REMOTE_NAME, ERR, ERROR )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: REMOTE_NAME !<The name of the object to import.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    INTEGER(INTG) :: FIELDML_OUTPUT_IMPORT
    
    CALL ENTERS( "FIELDML_OUTPUT_IMPORT", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_OUT( FIELDML_INFO, ERR, ERROR, *999 )

    FIELDML_OUTPUT_IMPORT = FIELDML_OUTPUT_IMPORT_FML( FIELDML_INFO%FML_HANDLE, REMOTE_NAME, ERR, ERROR )
 
    CALL EXITS( "FIELDML_OUTPUT_IMPORT" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_IMPORT", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_IMPORT" )

  END FUNCTION FIELDML_OUTPUT_IMPORT
  
  !
  !================================================================================================================================
  !

  !>Import the given FieldML object if it is not already imported or local.
  FUNCTION FIELDML_OUTPUT_IMPORT_HANDLE( FML_HANDLE, HANDLE, ERR, ERROR )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: HANDLE !<The FieldML object to import.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    INTEGER(INTG) :: FIELDML_OUTPUT_IMPORT_HANDLE
    
    !Local variables
    INTEGER(INTG) :: IMPORT_INDEX, LOCAL_HANDLE
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: NAME
    INTEGER(INTG) :: LENGTH
    
    CALL ENTERS( "FIELDML_OUTPUT_IMPORT_HANDLE", ERR, ERROR, *999 )

    FIELDML_OUTPUT_IMPORT_HANDLE = FML_INVALID_HANDLE
    LENGTH = Fieldml_CopyObjectDeclaredName( FML_HANDLE, HANDLE, NAME, MAXSTRLEN )
    
    IF( Fieldml_IsObjectLocal( FML_HANDLE, HANDLE , 1 ) /= 1 ) THEN
      IF( LENGTH > 0 ) THEN
        LOCAL_HANDLE = Fieldml_GetObjectByName( FML_HANDLE, NAME(1:LENGTH)//C_NULL_CHAR )
        IF( LOCAL_HANDLE == FML_INVALID_HANDLE ) THEN
          IMPORT_INDEX = Fieldml_AddImportSource( FML_HANDLE, &
            & "http://www.fieldml.org/resources/xml/0.5/FieldML_Library_0.5.xml"//C_NULL_CHAR, "library"//C_NULL_CHAR )
          FIELDML_OUTPUT_IMPORT_HANDLE = Fieldml_AddImport( FML_HANDLE, &
            & IMPORT_INDEX, NAME(1:LENGTH)//C_NULL_CHAR, NAME(1:LENGTH)//C_NULL_CHAR )
        ELSE IF( LOCAL_HANDLE == HANDLE ) THEN
          FIELDML_OUTPUT_IMPORT_HANDLE = HANDLE
        ENDIF
      ENDIF
    ENDIF

    CALL EXITS( "FIELDML_OUTPUT_IMPORT_HANDLE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_IMPORT_HANDLE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_IMPORT_HANDLE" )

  END FUNCTION FIELDML_OUTPUT_IMPORT_HANDLE
  
  !
  !================================================================================================================================
  !
  
  !>Get the argument corresponding to the given type (named *.argument), importing it if needed.
  FUNCTION FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, TYPE_HANDLE, DO_IMPORT, ERR, ERROR )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the argument.
    INTEGER(INTG), INTENT(IN) :: TYPE_HANDLE !<The type out of whose name the argument name is built.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG) :: FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE

    !Local variables
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: NAME
    INTEGER(INTG) :: LENGTH
    INTEGER(INTG) :: HANDLE, FML_ERR
    TYPE(VARYING_STRING) :: FULL_NAME
    
    CALL ENTERS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE", ERR, ERROR, *999 )

    LENGTH = Fieldml_CopyObjectName( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE, NAME, MAXSTRLEN )
    IF( LENGTH < 1 ) THEN
      LENGTH = Fieldml_CopyObjectDeclaredName( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE, NAME, MAXSTRLEN )
      FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE = FML_INVALID_HANDLE
      CALL EXITS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE" )
      RETURN
    ENDIF

    IF( DO_IMPORT ) THEN
      FULL_NAME = NAME(1:LENGTH)//".argument"
      FML_ERR = FIELDML_OUTPUT_IMPORT( FIELDML_INFO, FULL_NAME, ERR, ERROR )
      IF(ERR/=0) GOTO 999
    ENDIF
    
    HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, NAME(1:LENGTH)//".argument"//C_NULL_CHAR )
    IF( HANDLE == FML_INVALID_HANDLE ) THEN
      FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE = FML_INVALID_HANDLE
      CALL EXITS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE" )
      RETURN
    ENDIF
    
    FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE = HANDLE

    CALL ENTERS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE", ERR, ERROR, *999 )
    CALL EXITS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE" )
    
  END FUNCTION FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE
  
  !
  !================================================================================================================================
  !
  
  !>Get the FieldML built-in library type corresponding to the given OpenCMISS coordinate system type.
  SUBROUTINE FIELDML_OUTPUT_GET_COORDINATES_TYPE( FIELDML_HANDLE, COORDS_TYPE, DIMENSIONS, DO_IMPORT, TYPE_HANDLE, &
    & ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: COORDS_TYPE !<The OpenCMISS coordinates type.
    INTEGER(INTG), INTENT(IN) :: DIMENSIONS !<The coordinate system's number of dimensions.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the FieldML type.
    INTEGER(INTG), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type corresponding to the given OpenCMISS coordinate system type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: TYPE_NAME
    INTEGER(INTG) :: TEMP

    CALL ENTERS( "FIELDML_OUTPUT_GET_COORDINATES_TYPE", ERR, ERROR, *999 )
    
    IF( COORDS_TYPE == COORDINATE_RECTANGULAR_CARTESIAN_TYPE ) THEN
      IF( DIMENSIONS == 1 ) THEN
        TYPE_NAME = "coordinates.rc.1d"
      ELSE IF( DIMENSIONS == 2 ) THEN
        TYPE_NAME = "coordinates.rc.2d"
      ELSE IF( DIMENSIONS == 3 ) THEN
        TYPE_NAME = "coordinates.rc.3d"
      ELSE
        TYPE_HANDLE = FML_INVALID_HANDLE
        CALL FLAG_ERROR( var_str("Cannot get FieldML RC coordinates type of dimension ")//DIMENSIONS//".", ERR, ERROR, *999)
      ENDIF
    ELSE
      TYPE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( var_str("Cannot get FieldML coordinates for OpenCMISS type ")//COORDS_TYPE//".", ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_OUTPUT_IMPORT_FML( FIELDML_HANDLE, TYPE_NAME, ERR, ERROR )
      IF(ERR/=0) GOTO 999
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(TYPE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get FieldML coordinates type "//char(TYPE_NAME)//".", FIELDML_HANDLE, &
      & ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_OUTPUT_GET_COORDINATES_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_COORDINATES_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_COORDINATES_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_COORDINATES_TYPE
  
  !
  !================================================================================================================================
  !
  
  !>Returns a generic n-dimensional real type from the built-in library.
  SUBROUTINE FIELDML_OUTPUT_GET_GENERIC_TYPE( FIELDML_HANDLE, DIMENSIONS, TYPE_HANDLE, DO_IMPORT, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: DIMENSIONS !<The number of dimensions of the type.
    INTEGER(INTG), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    TYPE(VARYING_STRING) :: TYPE_NAME
    INTEGER(INTG) :: TEMP

    CALL ENTERS( "FIELDML_OUTPUT_GET_GENERIC_TYPE", ERR, ERROR, *999 )

    IF( DIMENSIONS == 1 ) THEN
      TYPE_NAME = "real.1d"
    ELSE IF( DIMENSIONS == 2 ) THEN
      TYPE_NAME = "real.2d"
    ELSE IF( DIMENSIONS == 3 ) THEN
      TYPE_NAME = "real.3d"
    ELSE
      TYPE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( var_str("Cannot get FieldML generic type of dimensionality ")//DIMENSIONS//".", ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_OUTPUT_IMPORT_FML( FIELDML_HANDLE, TYPE_NAME, ERR, ERROR )
      IF(ERR/=0) GOTO 999
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(TYPE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get generic type "//TYPE_NAME//".", FIELDML_HANDLE, ERR, ERROR, *999 )
    
    CALL EXITS( "FIELDML_OUTPUT_GET_GENERIC_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_GENERIC_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_GENERIC_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_GENERIC_TYPE
  
  !
  !================================================================================================================================
  !
  
  !>Returns a type in the built-in library corresponding to a chart of the given dimensionality.
  SUBROUTINE FIELDML_OUTPUT_GET_XI_TYPE( FIELDML_HANDLE, DIMENSIONS, DO_IMPORT, TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: DIMENSIONS !<The number of dimensions of the chart.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Locals
    INTEGER(INTG) :: TEMP
    TYPE(VARYING_STRING) :: TYPE_NAME

    CALL ENTERS( "FIELDML_OUTPUT_GET_XI_TYPE", ERR, ERROR, *999 )
    
    IF( DIMENSIONS == 1 ) THEN
      TYPE_NAME = "chart.1d"
    ELSE IF( DIMENSIONS == 2 ) THEN
      TYPE_NAME = "chart.2d"
    ELSE IF( DIMENSIONS == 3 ) THEN
      TYPE_NAME = "chart.3d"
    ELSE
      TYPE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( var_str("Chart dimensionality ")//DIMENSIONS//" not supported.", ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_OUTPUT_IMPORT_FML( FIELDML_HANDLE, TYPE_NAME, ERR, ERROR )
      IF(ERR/=0) GOTO 999
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(TYPE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get xi type "//TYPE_NAME//".", FIELDML_HANDLE, ERR, ERROR, *999 )
    
    CALL EXITS( "FIELDML_OUTPUT_GET_XI_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_XI_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_XI_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_XI_TYPE
  
  !
  !================================================================================================================================
  !
  
  !>Returns a FieldML type appropriate for the given OpenCMISS field.
  SUBROUTINE FIELDML_OUTPUT_GET_VALUE_TYPE( FML_HANDLE, FIELD, VARIABLE_TYPE, DO_IMPORT, TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field whose type is to be obtained.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the type.
    INTEGER(INTG), INTENT(OUT) :: TYPE_HANDLE !<The FieldML type handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    INTEGER(INTG) :: FIELDTYPE, SUB_TYPE, COUNT
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(REGION_TYPE), POINTER :: REGION
    
    CALL ENTERS( "FIELDML_OUTPUT_GET_VALUE_TYPE", ERR, ERROR, *999 )
    
    REGION => FIELD%REGION
    
    CALL FIELD_TYPE_GET( FIELD, FIELDTYPE, ERR, ERROR, *999 )
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( FIELD, VARIABLE_TYPE, COUNT, ERR, ERROR, *999 )

    SELECT CASE( FIELDTYPE )
    CASE( FIELD_GEOMETRIC_TYPE )
      NULLIFY( COORDINATE_SYSTEM )
      CALL REGION_COORDINATE_SYSTEM_GET( REGION, COORDINATE_SYSTEM, ERR, ERROR, *999 )
      CALL COORDINATE_SYSTEM_TYPE_GET( COORDINATE_SYSTEM, SUB_TYPE, ERR, ERROR, *999 )
      CALL FIELDML_OUTPUT_GET_COORDINATES_TYPE( FML_HANDLE, SUB_TYPE, COUNT, DO_IMPORT, TYPE_HANDLE, ERR, ERROR, *999 )
    
    !CASE( CMISSFieldFibreType )

    !CASE( CMISSFieldGeneralType )

    !CASE( CMISSFieldMaterialType )

    CASE DEFAULT
      CALL FIELDML_OUTPUT_GET_GENERIC_TYPE( FML_HANDLE, COUNT, TYPE_HANDLE, DO_IMPORT, ERR, ERROR, *999 )
    END SELECT
  
    CALL EXITS( "FIELDML_OUTPUT_GET_VALUE_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_VALUE_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_VALUE_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_VALUE_TYPE
    
  !
  !================================================================================================================================
  !

  !>Get an evaluator from the built-in library that corresponds to the given OpenCMISS tensor-product basis.
  SUBROUTINE FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR( FML_HANDLE, XI_INTERPOLATIONS, COLLAPSE_INFO, EVALUATOR_HANDLE, &
    & PARAMETERS_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle
    INTEGER(INTG), INTENT(IN) :: XI_INTERPOLATIONS(:) !<The per-xi interpolations used by the basis.
    INTEGER(INTG), INTENT(IN) :: COLLAPSE_INFO(:) !<The basis collapse info.
    INTEGER(INTG), INTENT(OUT) :: EVALUATOR_HANDLE !<The evaluator handle for the basis.
    INTEGER(INTG), INTENT(OUT) :: PARAMETERS_HANDLE !<The basis parameters argument evaluator.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: XI_COUNT, FIRST_INTERPOLATION, I
    TYPE(VARYING_STRING) :: SUFFIX, INTERPOLATOR_NAME, PARAMETER_NAME

    CALL ENTERS( "FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR", ERR, ERROR, *999 )
    
    XI_COUNT = SIZE( XI_INTERPOLATIONS )
  
    DO I = 1, XI_COUNT
      IF( I == 1 ) THEN
        FIRST_INTERPOLATION = XI_INTERPOLATIONS(I)
      ELSE IF( XI_INTERPOLATIONS(I) /= FIRST_INTERPOLATION ) THEN
        !Do not yet support inhomogeneous TP bases
        CALL FLAG_ERROR( "Translation of inhomogeneous tensor-product basis not yet supported.", ERR, ERROR, *999 )
      ENDIF
    ENDDO
   
    CALL FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX( COLLAPSE_INFO, SUFFIX, ERR, ERROR, *999 )

    EVALUATOR_HANDLE = FML_INVALID_HANDLE
    PARAMETERS_HANDLE = FML_INVALID_HANDLE
      
    IF( FIRST_INTERPOLATION == BASIS_QUADRATIC_LAGRANGE_INTERPOLATION ) THEN
      IF( XI_COUNT == 1 ) THEN
        INTERPOLATOR_NAME = "interpolator.1d.unit.quadraticLagrange"
        PARAMETER_NAME = "parameters.1d.unit.quadraticLagrange"
      ELSE IF( XI_COUNT == 2 ) THEN
        INTERPOLATOR_NAME = "interpolator.2d.unit.biquadraticLagrange"//SUFFIX
        PARAMETER_NAME = "parameters.2d.unit.biquadraticLagrange"//SUFFIX
      ELSE IF( XI_COUNT == 3 ) THEN
        INTERPOLATOR_NAME = "interpolator.3d.unit.triquadraticLagrange"//SUFFIX
        PARAMETER_NAME = "parameters.3d.unit.triquadraticLagrange"//SUFFIX
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( var_str("Quadratic Lagrangian interpolation not supported for ")//XI_COUNT//" dimensions.", &
          & ERR, ERROR, *999 )
      ENDIF
    ELSE IF( FIRST_INTERPOLATION == BASIS_LINEAR_LAGRANGE_INTERPOLATION ) THEN
      IF( XI_COUNT == 1 ) THEN
        INTERPOLATOR_NAME = "interpolator.1d.unit.linearLagrange"
        PARAMETER_NAME = "parameters.1d.unit.linearLagrange"
      ELSE IF( XI_COUNT == 2 ) THEN
        INTERPOLATOR_NAME = "interpolator.2d.unit.bilinearLagrange"//SUFFIX
        PARAMETER_NAME = "parameters.2d.unit.bilinearLagrange"//SUFFIX
      ELSE IF( XI_COUNT == 3 ) THEN
        INTERPOLATOR_NAME = "interpolator.3d.unit.trilinearLagrange"//SUFFIX
        PARAMETER_NAME = "parameters.3d.unit.trilinearLagrange"//SUFFIX
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( var_str("Quadratic Lagrangian interpolation not supported for ")//XI_COUNT//" dimensions.", &
          & ERR, ERROR, *999 )
      ENDIF
    ELSE
      CALL FLAG_ERROR( var_str("FieldML translation not yet supported for interpolation type ")//FIRST_INTERPOLATION//".", &
        & ERR, ERROR, *999 )
    ENDIF

    EVALUATOR_HANDLE = FIELDML_OUTPUT_IMPORT_FML( FML_HANDLE, INTERPOLATOR_NAME, ERR, ERROR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not import interpolator "//char(INTERPOLATOR_NAME)//".", ERR, ERROR, *999 )

    PARAMETERS_HANDLE = FIELDML_OUTPUT_IMPORT_FML( FML_HANDLE, PARAMETER_NAME, ERR, ERROR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not import parameter type "//char(INTERPOLATOR_NAME)//".", ERR, ERROR, *999 )
    
    IF( EVALUATOR_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Cannot get a handle for basis evaluator "//char(INTERPOLATOR_NAME)//".", ERR, ERROR, *999 )
    ENDIF

    IF( PARAMETERS_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( "Cannot get a handle for basis parameters "//char(PARAMETER_NAME)//".", ERR, ERROR, *999 )
    ENDIF
    
    CALL EXITS( "FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR

  !
  !================================================================================================================================
  !
  
  !>Return the FieldML connectivity ensemble corresponding to the given tensor-product basis info.
  SUBROUTINE FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE( FIELDML_HANDLE, XI_INTERPOLATIONS, COLLAPSE_INFO, DO_IMPORT, TYPE_HANDLE, &
    & ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    INTEGER(INTG), INTENT(IN) :: XI_INTERPOLATIONS(:) !<The per-xi interpolation of the given TP basis.
    INTEGER(INTG), INTENT(IN) :: COLLAPSE_INFO(:) !<The collapse-constant for the given basis.
    LOGICAL, INTENT(IN) :: DO_IMPORT !<If true, import the connectivity ensemble.
    INTEGER(INTG), INTENT(OUT) :: TYPE_HANDLE !<The FieldML ensemble handle.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    INTEGER(INTG) :: XI_COUNT, FIRST_INTERPOLATION, I, IMPORT_INDEX, TEMP
    TYPE(VARYING_STRING) :: SUFFIX, LAYOUT_NAME
    
    CALL ENTERS( "FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE", ERR, ERROR, *999 )

    XI_COUNT = SIZE( XI_INTERPOLATIONS )
  
    IMPORT_INDEX = Fieldml_AddImportSource( FIELDML_HANDLE, &
      & "http://www.fieldml.org/resources/xml/0.5/FieldML_Library_0.5.xml"//C_NULL_CHAR, "library"//C_NULL_CHAR )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot access built-in FieldML library.", FIELDML_HANDLE, ERR, ERROR, *999 )

    FIRST_INTERPOLATION = XI_INTERPOLATIONS(1)
    DO I = 2, XI_COUNT
      IF( XI_INTERPOLATIONS(I) /= FIRST_INTERPOLATION ) THEN
        !Do not yet support inhomogeneous TP bases
        CALL FLAG_ERROR( "FieldML translation of inhomogeneous tensor-product bases are not yet supported.", &
          & ERR, ERROR, *999 )
      ENDIF
    ENDDO

    CALL FIELDML_OUTPUT_GET_COLLAPSE_SUFFIX( COLLAPSE_INFO, SUFFIX, ERR, ERROR, *999 )
      
    IF( FIRST_INTERPOLATION == BASIS_QUADRATIC_LAGRANGE_INTERPOLATION ) THEN
      IF( XI_COUNT == 1 ) THEN
        LAYOUT_NAME = "localNodes.1d.line3"
      ELSE IF( XI_COUNT == 2 ) THEN
        LAYOUT_NAME = "localNodes.2d.square3x3"//SUFFIX
      ELSE IF( XI_COUNT == 3 ) THEN
        LAYOUT_NAME = "localNodes.3d.cube3x3x3"//SUFFIX
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( var_str("Quadratic Lagrangian interpolation not supported for ")//XI_COUNT//" dimensions.", &
          & ERR, ERROR, *999 )
      ENDIF
    ELSE IF( FIRST_INTERPOLATION == BASIS_LINEAR_LAGRANGE_INTERPOLATION ) THEN
      IF( XI_COUNT == 1 ) THEN
        LAYOUT_NAME = "localNodes.1d.line2"
      ELSE IF( XI_COUNT == 2 ) THEN
        LAYOUT_NAME = "localNodes.2d.square2x2"//SUFFIX
      ELSE IF( XI_COUNT == 3 ) THEN
        LAYOUT_NAME = "localNodes.3d.cube2x2x2"//SUFFIX
      ELSE
        !Do not yet support dimensions higher than 3.
        CALL FLAG_ERROR( var_str("Linear Lagrangian interpolation not supported for ")//XI_COUNT//" dimensions.", &
          & ERR, ERROR, *999 )
      ENDIF
    ELSE
      CALL FLAG_ERROR( var_str("FieldML translation not yet supported for interpolation type ")//FIRST_INTERPOLATION//".", &
        & ERR, ERROR, *999 )
    ENDIF

    IF( DO_IMPORT ) THEN
      TEMP = FIELDML_OUTPUT_IMPORT_FML( FIELDML_HANDLE, LAYOUT_NAME, ERR, ERROR )
      IF(ERR/=0) GOTO 999
    ENDIF
    TYPE_HANDLE = Fieldml_GetObjectByName( FIELDML_HANDLE, cchar(LAYOUT_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get local nodes type "//LAYOUT_NAME//".", FIELDML_HANDLE, ERR, ERROR, *999 )
    
    CALL EXITS( "FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE

  !
  !================================================================================================================================
  !

  !>Get the connectivity ensemble for the given basis. Currently, only tensor-product bases are supported.
  SUBROUTINE FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE( FIELDML_HANDLE, BASIS, TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FIELDML_HANDLE !<The FieldML session handle.
    TYPE(BASIS_TYPE), POINTER :: BASIS !<The basis for which to return the connectivity.
    INTEGER(INTG), INTENT(OUT) :: TYPE_HANDLE !<The FieldML connectivity ensemble.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Locals
    INTEGER(INTG) :: BASISTYPE, XI_COUNT
    INTEGER(INTG), ALLOCATABLE :: XI_INTERPOLATIONS(:), COLLAPSE_INFO(:)
    
    CALL ENTERS( "FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE", ERR, ERROR, *999 )

    TYPE_HANDLE = FML_INVALID_HANDLE

    CALL BASIS_TYPE_GET( BASIS, BASISTYPE, ERR, ERROR, *999 )
    CALL BASIS_NUMBER_OF_XI_GET( BASIS, XI_COUNT, ERR, ERROR, *999 )
    
    IF( BASISTYPE == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      ALLOCATE( XI_INTERPOLATIONS( XI_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate xi interpolations array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE_INFO( XI_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate collapse info array.", ERR, ERROR, *999 )
      CALL BASIS_INTERPOLATION_XI_GET( BASIS, XI_INTERPOLATIONS, ERR, ERROR, *999 )
      CALL BASIS_COLLAPSED_XI_GET( BASIS, COLLAPSE_INFO, ERR, ERROR, *999 )
      
      CALL FIELDML_OUTPUT_GET_TP_CONNECTIVITY_TYPE( FIELDML_HANDLE, XI_INTERPOLATIONS, COLLAPSE_INFO, .TRUE., TYPE_HANDLE, &
        & ERR, ERROR, *999 )
      
      DEALLOCATE( XI_INTERPOLATIONS )
      DEALLOCATE( COLLAPSE_INFO )
    ELSE
      CALL FLAG_ERROR( "Only translation of tensor product bases are currently supported", ERR, ERROR, *999 )
    ENDIF
    
    CALL EXITS( "FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE

  !
  !================================================================================================================================
  !
  
  !>Returns the index of the layout handle used by the given connectivity info array, or -1 if none can be found.
  FUNCTION FIELDML_OUTPUT_FIND_LAYOUT( CONNECTIVITY_INFO, LAYOUT_HANDLE, ERR, ERROR )
    !Argument variables
    TYPE(CONNECTIVITY_INFO_TYPE), INTENT(IN) :: CONNECTIVITY_INFO(:) !<The connectivity info array to search.
    INTEGER(INTG), INTENT(IN) :: LAYOUT_HANDLE !<The local node layout handle to search for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Function
    INTEGER(INTG) :: FIELDML_OUTPUT_FIND_LAYOUT
    
    !Locals
    INTEGER(INTG) :: I
    
    CALL ENTERS( "FIELDML_OUTPUT_FIND_LAYOUT", ERR, ERROR, *999 )

    FIELDML_OUTPUT_FIND_LAYOUT = -1
    DO I = 1, SIZE( CONNECTIVITY_INFO )
      IF( CONNECTIVITY_INFO(I)%LAYOUT_HANDLE == LAYOUT_HANDLE ) THEN
        FIELDML_OUTPUT_FIND_LAYOUT = I
      ENDIF
    ENDDO
  
    CALL EXITS( "FIELDML_OUTPUT_FIND_LAYOUT" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_FIND_LAYOUT", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_FIND_LAYOUT" )

  END FUNCTION FIELDML_OUTPUT_FIND_LAYOUT
  
  !
  !================================================================================================================================
  !
  
  !>Returns the index of the basis handle used by the given basis info array, or -1 if none can be found.
  FUNCTION FIELDML_OUTPUT_FIND_BASIS( BASIS_INFO, BASIS, ERR, ERROR )
    !Argument variables
    TYPE(BASIS_INFO_TYPE), INTENT(IN) :: BASIS_INFO(:) !<The basis info array to search.
    TYPE(BASIS_TYPE), POINTER, INTENT(IN) :: BASIS !<The basis handle to search for. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !Function
    INTEGER(INTG) :: FIELDML_OUTPUT_FIND_BASIS
    
    !Locals
    INTEGER(INTG) :: I
    
    CALL ENTERS( "FIELDML_OUTPUT_FIND_BASIS", ERR, ERROR, *999 )

    FIELDML_OUTPUT_FIND_BASIS = -1
    DO I = 1, SIZE( BASIS_INFO )
      IF( ASSOCIATED( BASIS_INFO(I)%BASIS, TARGET = BASIS ) ) THEN
        FIELDML_OUTPUT_FIND_BASIS = I
      ENDIF
    ENDDO
  
    CALL EXITS( "FIELDML_OUTPUT_FIND_BASIS" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_FIND_BASIS", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_FIND_BASIS" )

  END FUNCTION FIELDML_OUTPUT_FIND_BASIS
  
  !
  !================================================================================================================================
  !

  !>Returns the simplified name of the given layout. This is used for naming associated connectivity evaluators.
  SUBROUTINE FIELDML_OUTPUT_GET_SIMPLE_LAYOUT_NAME( FML_HANDLE, LAYOUT_HANDLE, NAME, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle
    INTEGER(INTG), INTENT(IN) :: LAYOUT_HANDLE !<The local node layout.
    TYPE(VARYING_STRING), INTENT(INOUT) :: NAME !<The simplified name.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: FULL_NAME
    INTEGER(INTG) :: LENGTH

    CALL ENTERS( "FIELDML_OUTPUT_GET_SIMPLE_LAYOUT_NAME", ERR, ERROR, *999 )
    
    LENGTH = Fieldml_CopyObjectDeclaredName( FML_HANDLE, LAYOUT_HANDLE, FULL_NAME, MAXSTRLEN )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR("Cannot get name of layout ensemble.", FML_HANDLE, ERR, ERROR, *999 )
    
    IF( INDEX( FULL_NAME, 'localNodes.') /= 1 ) THEN
      NAME = FULL_NAME(1:LENGTH)
    ELSE
      NAME = FULL_NAME(12:LENGTH)
    ENDIF

    CALL EXITS( "FIELDML_OUTPUT_GET_SIMPLE_LAYOUT_NAME" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_SIMPLE_LAYOUT_NAME", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_SIMPLE_LAYOUT_NAME" )
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !>Returns the simplified name of the given basis. This is used for naming associated reference evaluators.
  SUBROUTINE FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME( FML_HANDLE, BASIS_HANDLE, NAME, ERR, ERROR, * )
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FML_HANDLE !<The FieldML session handle
    INTEGER(INTG), INTENT(IN) :: BASIS_HANDLE !<The basis handle.
    TYPE(VARYING_STRING), INTENT(INOUT) :: NAME !<The simplified name.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: FULL_NAME
    INTEGER(INTG) :: LENGTH
    
    CALL ENTERS( "FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME", ERR, ERROR, *999 )

    LENGTH = Fieldml_CopyObjectDeclaredName( FML_HANDLE, BASIS_HANDLE, FULL_NAME, MAXSTRLEN )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR("Cannot get name of basis evaluator.", FML_HANDLE, ERR, ERROR, *999 )
    
    IF( INDEX( FULL_NAME, 'interpolator.1d.unit.') == 1 ) THEN
      NAME = FULL_NAME(22:LENGTH)
    ELSEIF( INDEX( FULL_NAME, 'interpolator.2d.unit.') == 1 ) THEN
      NAME = FULL_NAME(22:LENGTH)
    ELSEIF( INDEX( FULL_NAME, 'interpolator.3d.unit.') == 1 ) THEN
      NAME = FULL_NAME(22:LENGTH)
    ELSE
      NAME = FULL_NAME(1:LENGTH)
    ENDIF

    CALL EXITS( "FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME
    
  !
  !================================================================================================================================
  !
  
  !>Create a basis evaluator from the given basis info.
  SUBROUTINE FIELDML_OUTPUT_CREATE_BASIS_REFERENCE( FIELDML_INFO, BASE_NAME, BASIS_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(BASIS_INFO_TYPE), INTENT(INOUT) :: BASIS_INFO !<The basis info describing the basis to create.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: BASIS_TYPE, XI_COUNT, INTERPOLATION_PARAMETERS_HANDLE, HANDLE, EVALUATOR_HANDLE
    INTEGER(INTG) :: VARIABLE_HANDLE, AGGREGATE_HANDLE, INDEX_EVALUATOR_HANDLE, FML_ERR
    INTEGER(INTG), ALLOCATABLE :: XI_INTERPOLATIONS(:), COLLAPSE_INFO(:)
    TYPE(VARYING_STRING) :: REFERENCE_NAME, NAME

    CALL ENTERS( "FIELDML_OUTPUT_CREATE_BASIS_REFERENCE", ERR, ERROR, *999 )
    
    CALL BASIS_TYPE_GET( BASIS_INFO%BASIS, BASIS_TYPE, ERR, ERROR, *999 )
    CALL BASIS_NUMBER_OF_XI_GET( BASIS_INFO%BASIS, XI_COUNT, ERR, ERROR, *999 )
    
    IF( BASIS_TYPE == BASIS_LAGRANGE_HERMITE_TP_TYPE ) THEN
      ALLOCATE( XI_INTERPOLATIONS( XI_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate xi interpolation array.", ERR, ERROR, *999 )
      ALLOCATE( COLLAPSE_INFO( XI_COUNT ), STAT = ERR )
      CALL BASIS_INTERPOLATION_XI_GET( BASIS_INFO%BASIS, XI_INTERPOLATIONS, ERR, ERROR, *999 )
      CALL BASIS_COLLAPSED_XI_GET( BASIS_INFO%BASIS, COLLAPSE_INFO, ERR, ERROR, *999 )
      
      CALL FIELDML_OUTPUT_GET_TP_BASIS_EVALUATOR( FIELDML_INFO%FML_HANDLE, XI_INTERPOLATIONS, COLLAPSE_INFO, EVALUATOR_HANDLE, &
        & INTERPOLATION_PARAMETERS_HANDLE, ERR, ERROR, *999 )
      DEALLOCATE( XI_INTERPOLATIONS )
      DEALLOCATE( COLLAPSE_INFO )

      CALL FIELDML_OUTPUT_GET_SIMPLE_BASIS_NAME( FIELDML_INFO%FML_HANDLE, EVALUATOR_HANDLE, NAME, ERR, ERROR, *999 )
      
      REFERENCE_NAME = BASE_NAME//NAME//"_"//TRIM(NUMBER_TO_VSTRING(BASIS_INFO%BASIS%USER_NUMBER,"*",ERR,ERROR))// &
        & ".parameters"
      
      AGGREGATE_HANDLE = Fieldml_CreateAggregateEvaluator( FIELDML_INFO%FML_HANDLE, cchar(REFERENCE_NAME), &
        & INTERPOLATION_PARAMETERS_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create dofs for basis connectivity for "//NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      INDEX_EVALUATOR_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, BASIS_INFO%LAYOUT_HANDLE, .TRUE., &
        & ERR, ERROR )
      IF(ERR/=0) GOTO 999

      FML_ERR = Fieldml_SetIndexEvaluator( FIELDML_INFO%FML_HANDLE, AGGREGATE_HANDLE, 1, INDEX_EVALUATOR_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set field component index evaluator for "//REFERENCE_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
      FML_ERR = Fieldml_SetDefaultEvaluator( FIELDML_INFO%FML_HANDLE, AGGREGATE_HANDLE, FIELDML_INFO%NODE_DOFS_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal field dofs for "//REFERENCE_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, BASIS_INFO%CONNECTIVITY_HANDLE )
      VARIABLE_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, HANDLE, .FALSE., ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_SetBind( FIELDML_INFO%FML_HANDLE, AGGREGATE_HANDLE, VARIABLE_HANDLE, BASIS_INFO%CONNECTIVITY_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set bind for basis dofs for"//REFERENCE_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
      REFERENCE_NAME = BASE_NAME//NAME//"_"//TRIM(NUMBER_TO_VSTRING(BASIS_INFO%BASIS%USER_NUMBER,"*",ERR,ERROR))// &
        & ".evaluator"

      BASIS_INFO%REFERENCE_HANDLE = Fieldml_CreateReferenceEvaluator( FIELDML_INFO%FML_HANDLE, cchar(REFERENCE_NAME), &
        & EVALUATOR_HANDLE )

      CALL FIELDML_OUTPUT_GET_XI_TYPE( FIELDML_INFO%FML_HANDLE, XI_COUNT, .TRUE., HANDLE, ERR, ERROR, *999 )
      VARIABLE_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, HANDLE, .TRUE., ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_SetBind( FIELDML_INFO%FML_HANDLE, BASIS_INFO%REFERENCE_HANDLE, VARIABLE_HANDLE, &
        & FIELDML_INFO%XI_ARGUMENT_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot bind xi to basis evaluator "//REFERENCE_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      VARIABLE_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, INTERPOLATION_PARAMETERS_HANDLE, .TRUE.,&
        & ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_SetBind( FIELDML_INFO%FML_HANDLE, BASIS_INFO%REFERENCE_HANDLE, VARIABLE_HANDLE, &
        & AGGREGATE_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot bind parameters to basis evaluator "//REFERENCE_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ELSE
      BASIS_INFO%REFERENCE_HANDLE = FML_INVALID_HANDLE
      CALL FLAG_ERROR( "FieldML export code can currently only translate tensor-product bases.", ERR, ERROR, *999 )
    ENDIF
    
    CALL EXITS( "FIELDML_OUTPUT_CREATE_BASIS_REFERENCE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_CREATE_BASIS_REFERENCE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_CREATE_BASIS_REFERENCE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_CREATE_BASIS_REFERENCE

  !
  !================================================================================================================================
  !
  
  !>Create a parameter evaluator for the given local node layout.
  SUBROUTINE FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS( FIELDML_INFO, LAYOUT_HANDLE, COMPONENT_NAME, &
    & CONNECTIVITY_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: LAYOUT_HANDLE !<The local node layout.
    TYPE(VARYING_STRING), INTENT(IN) :: COMPONENT_NAME !<The component name.
    TYPE(CONNECTIVITY_INFO_TYPE), INTENT(INOUT) :: CONNECTIVITY_INFO !<The connectivity info for the local node layout.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    TYPE(VARYING_STRING) :: NAME
    INTEGER(INTG) :: INDEX_HANDLE, FML_ERR
    TYPE(VARYING_STRING) :: CONNECTIVITY_NAME

    CALL ENTERS( "FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS", ERR, ERROR, *999 )

    CALL FIELDML_OUTPUT_GET_SIMPLE_LAYOUT_NAME( FIELDML_INFO%FML_HANDLE, LAYOUT_HANDLE, NAME, ERR, ERROR, *999 )
    CONNECTIVITY_NAME = COMPONENT_NAME//NAME

    CONNECTIVITY_INFO%LAYOUT_HANDLE = LAYOUT_HANDLE
    CONNECTIVITY_INFO%CONNECTIVITY_HANDLE = Fieldml_CreateParameterEvaluator( FIELDML_INFO%FML_HANDLE, &
      & cchar(CONNECTIVITY_NAME), FIELDML_INFO%NODES_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR("Cannot create nodal parameters for "//CONNECTIVITY_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FML_ERR = Fieldml_SetParameterDataDescription( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_INFO%CONNECTIVITY_HANDLE, &
      & FML_DATA_DESCRIPTION_DENSE_ARRAY )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR("Cannot set nodal parameters description for "//CONNECTIVITY_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_INFO%CONNECTIVITY_HANDLE, &
      & FIELDML_INFO%ELEMENTS_ARGUMENT_HANDLE, FML_INVALID_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR("Cannot add element index to nodal parameters "//CONNECTIVITY_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    INDEX_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, LAYOUT_HANDLE, .TRUE., ERR, ERROR )
    IF(ERR/=0) GOTO 999
    FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_INFO%CONNECTIVITY_HANDLE, INDEX_HANDLE, &
      & FML_INVALID_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR("Cannot add layout index to nodal parameters "//CONNECTIVITY_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS

  !
  !================================================================================================================================
  !

  !>Add an evaluator corresponding to the given component of the given OpenCMISS mesh.
  SUBROUTINE FIELDML_OUTPUT_ADD_MESH_COMPONENT( FIELDML_INFO, BASE_NAME, CONNECTIVITY_FORMAT, COMPONENT_NUMBER, &
    & MESH_ELEMENTS, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(INOUT) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: CONNECTIVITY_FORMAT !<The name of the format to use when writing connectivity data.
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The mesh component number for which an evaluator should be constructed.
    TYPE(MESH_ELEMENTS_TYPE), POINTER, INTENT(IN) :: MESH_ELEMENTS !<The mesh element from which to obtain topology info.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: LAYOUT_HANDLE, CONNECTIVITY_HANDLE, ELEMENT_COUNT, DEFAULT_HANDLE, TEMPLATE_HANDLE, TYPE_HANDLE
    INTEGER(INTG) :: CONNECTIVITY_COUNT, BASIS_COUNT, I, J, LAYOUT_NODE_COUNT, IDX
    INTEGER(INTG), ALLOCATABLE, TARGET :: IBUFFER(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG) :: WRITER, SOURCE_HANDLE, FML_ERR, RESOURCE_HANDLE
    TYPE(CONNECTIVITY_INFO_TYPE), ALLOCATABLE :: CONNECTIVITY_INFO(:), TEMP_CONNECTIVITY_INFO(:)
    TYPE(BASIS_INFO_TYPE), ALLOCATABLE :: BASIS_INFO(:), TEMP_BASIS_INFO(:)
    TYPE(VARYING_STRING) :: COMPONENT_NAME, ARRAY_LOCATION
    INTEGER(INTG), TARGET :: OFFSETS(2), SIZES(2)
    
    CALL ENTERS( "FIELDML_OUTPUT_ADD_MESH_COMPONENT", ERR, ERROR, *999 )

    ELEMENT_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%ELEMENTS_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get element count for mesh "//BASE_NAME//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    CONNECTIVITY_COUNT = 0
    BASIS_COUNT = 0
    
    COMPONENT_NAME = BASE_NAME//".component"//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))
    
    TYPE_HANDLE = Fieldml_GetValueType( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%NODE_DOFS_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get node dofs FieldML type.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    RESOURCE_HANDLE = Fieldml_CreateHrefDataResource( FIELDML_INFO%FML_HANDLE, &
      & cchar(COMPONENT_NAME//".connectivity.resource"), cchar( CONNECTIVITY_FORMAT ), &
      & cchar(COMPONENT_NAME//".connectivity") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh component connectivity resource "//COMPONENT_NAME//&
      & ".connectivity.resource", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    TEMPLATE_HANDLE = Fieldml_CreatePiecewiseEvaluator( FIELDML_INFO%FML_HANDLE, cchar(COMPONENT_NAME//".template"), &
      &  TYPE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh component template "//COMPONENT_NAME//".template.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetIndexEvaluator( FIELDML_INFO%FML_HANDLE, TEMPLATE_HANDLE, 1, FIELDML_INFO%ELEMENTS_ARGUMENT_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( &
      & "Cannot set index evaluator for mesh omponent template "//COMPONENT_NAME//".template.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    DO I = 1, ELEMENT_COUNT
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET( I, MESH_ELEMENTS, BASIS, ERR, ERROR, *999 )

      CALL FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE( FIELDML_INFO%FML_HANDLE, BASIS, LAYOUT_HANDLE, ERR, ERROR, *999 )
      
      IDX = -1
      IF( CONNECTIVITY_COUNT > 0 ) THEN
        IDX = FIELDML_OUTPUT_FIND_LAYOUT( CONNECTIVITY_INFO, LAYOUT_HANDLE, ERR, ERROR )
        IF(ERR/=0) GOTO 999
      ENDIF

      IF( IDX == -1 ) THEN
        IF( CONNECTIVITY_COUNT == 0 ) THEN
          ALLOCATE( CONNECTIVITY_INFO( CONNECTIVITY_COUNT + 1 ), STAT = ERR )
          IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate connectivity info array.", ERR, ERROR, *999 )
        ELSE
          ALLOCATE( TEMP_CONNECTIVITY_INFO( CONNECTIVITY_COUNT ), STAT = ERR )
          IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate temporary connectivity array.", ERR, ERROR, *999 )
          TEMP_CONNECTIVITY_INFO(:) = CONNECTIVITY_INFO(:)
          DEALLOCATE( CONNECTIVITY_INFO )
          ALLOCATE( CONNECTIVITY_INFO( CONNECTIVITY_COUNT + 1 ), STAT = ERR )
          IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate new connectivity info array.", ERR, ERROR, *999 )
          CONNECTIVITY_INFO( 1:CONNECTIVITY_COUNT ) = TEMP_CONNECTIVITY_INFO( 1:CONNECTIVITY_COUNT )
        ENDIF
        
        CALL FIELDML_OUTPUT_CREATE_LAYOUT_PARAMETERS( FIELDML_INFO, LAYOUT_HANDLE, COMPONENT_NAME, &
          & CONNECTIVITY_INFO(CONNECTIVITY_COUNT+1), ERR, ERROR, *999 )
          
        LAYOUT_NODE_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, &
          & CONNECTIVITY_INFO(CONNECTIVITY_COUNT+1)%LAYOUT_HANDLE )
          
        ARRAY_LOCATION = ""
        ARRAY_LOCATION = ARRAY_LOCATION//( CONNECTIVITY_COUNT + 1 )
        SIZES(1) = ELEMENT_COUNT
        SIZES(2) = LAYOUT_NODE_COUNT
        SOURCE_HANDLE = Fieldml_CreateArrayDataSource( FIELDML_INFO%FML_HANDLE, cchar(COMPONENT_NAME//".connectivity"), &
          & RESOURCE_HANDLE, cchar(ARRAY_LOCATION), 2 )
        FML_ERR = Fieldml_SetArrayDataSourceRawSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, C_LOC(SIZES) )
        FML_ERR = Fieldml_SetArrayDataSourceSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, C_LOC(SIZES) )
        CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create connectivity data source "//COMPONENT_NAME//".connectivity", &
          & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

        FML_ERR = Fieldml_SetDataSource( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_INFO(CONNECTIVITY_COUNT+1)%CONNECTIVITY_HANDLE, &
          & SOURCE_HANDLE )
        CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set connectivity data source to "//COMPONENT_NAME//".connectivity.",&
          & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
  
        CONNECTIVITY_COUNT = CONNECTIVITY_COUNT + 1
        
        IDX = CONNECTIVITY_COUNT
      ENDIF
      CONNECTIVITY_HANDLE = CONNECTIVITY_INFO(IDX)%CONNECTIVITY_HANDLE

      IF( BASIS_COUNT == 0 ) THEN
        IDX = -1
      ELSE
        IDX = FIELDML_OUTPUT_FIND_BASIS( BASIS_INFO, BASIS, ERR, ERROR )
        IF(ERR/=0) GOTO 999
      ENDIF
      IF( IDX == -1 ) THEN
        IF( BASIS_COUNT == 0 ) THEN
          ALLOCATE( BASIS_INFO( BASIS_COUNT + 1 ), STAT = ERR )
          IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate basis info array.", ERR, ERROR, *999 )
        ELSE
          ALLOCATE( TEMP_BASIS_INFO( BASIS_COUNT ), STAT = ERR )
          IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate temporary basis info array.", ERR, ERROR, *999 )
          TEMP_BASIS_INFO(:) = BASIS_INFO(:)
          DEALLOCATE( BASIS_INFO )
          ALLOCATE( BASIS_INFO( BASIS_COUNT + 1 ), STAT = ERR )
          IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate new basis info array.", ERR, ERROR, *999 )
          BASIS_INFO( 1:BASIS_COUNT ) = TEMP_BASIS_INFO( 1:BASIS_COUNT )
        ENDIF

        BASIS_COUNT = BASIS_COUNT + 1
        BASIS_INFO( BASIS_COUNT )%BASIS => BASIS
        BASIS_INFO( BASIS_COUNT )%CONNECTIVITY_HANDLE = CONNECTIVITY_HANDLE
        BASIS_INFO( BASIS_COUNT )%LAYOUT_HANDLE = LAYOUT_HANDLE
        CALL FIELDML_OUTPUT_CREATE_BASIS_REFERENCE( FIELDML_INFO, COMPONENT_NAME, BASIS_INFO(BASIS_COUNT), ERR, ERROR, *999 )
        IDX = BASIS_COUNT
      ENDIF

      IF( I == 1 ) THEN
        DEFAULT_HANDLE = BASIS_INFO( IDX )%REFERENCE_HANDLE
        FML_ERR = Fieldml_SetDefaultEvaluator( FIELDML_INFO%FML_HANDLE, TEMPLATE_HANDLE, DEFAULT_HANDLE )
      ELSEIF( BASIS_INFO( IDX )%REFERENCE_HANDLE /= DEFAULT_HANDLE ) THEN
        FML_ERR = Fieldml_SetEvaluator( FIELDML_INFO%FML_HANDLE, TEMPLATE_HANDLE, I, BASIS_INFO( IDX )%REFERENCE_HANDLE )
      ENDIF
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set mesh connectivity evaluator to "//COMPONENT_NAME//".template.", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
    ENDDO

    DO I = 1, CONNECTIVITY_COUNT
      LAYOUT_NODE_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_INFO(I)%LAYOUT_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get layout node count.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
      SOURCE_HANDLE = Fieldml_GetDataSource( FIELDML_INFO%FML_HANDLE, CONNECTIVITY_INFO(I)%CONNECTIVITY_HANDLE )

      SIZES(1) = ELEMENT_COUNT
      SIZES(2) = LAYOUT_NODE_COUNT
      IF( I == 1 ) THEN
        WRITER = Fieldml_OpenArrayWriter( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, FIELDML_INFO%NODES_HANDLE, 0, C_LOC(SIZES), 2)
      ELSE
        WRITER = Fieldml_OpenArrayWriter( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, FIELDML_INFO%NODES_HANDLE, 1, C_LOC(SIZES), 2)
      ENDIF
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open connectivity data writer.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
      ALLOCATE( IBUFFER( LAYOUT_NODE_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate layout buffer.", ERR, ERROR, *999 )
      SIZES(1) = 1
      SIZES(2) = LAYOUT_NODE_COUNT
      OFFSETS(:) = 0
      DO J = 1, ELEMENT_COUNT
        CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET( J, MESH_ELEMENTS, BASIS, ERR, ERROR, *999 )
  
        CALL FIELDML_OUTPUT_GET_CONNECTIVITY_ENSEMBLE( FIELDML_INFO%FML_HANDLE, BASIS, LAYOUT_HANDLE, ERR, ERROR, *999 )
        IF( LAYOUT_HANDLE == CONNECTIVITY_INFO(I)%LAYOUT_HANDLE ) THEN
          CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET( J, MESH_ELEMENTS, IBUFFER, ERR, ERROR, *999 )
        ELSE
          IBUFFER(:) = 0
        ENDIF
        FML_ERR = Fieldml_WriteIntSlab( WRITER, C_LOC(OFFSETS), C_LOC(SIZES), C_LOC(IBUFFER) )
        IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
          CALL FLAG_ERROR( var_str("I/O error while writing connectivity data for ")//BASE_NAME//"("&
            & // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", &
            & ERR, ERROR, *999 )
        ENDIF
        OFFSETS(1) = OFFSETS(1) + 1
      ENDDO
      DEALLOCATE( IBUFFER )
      FML_ERR = Fieldml_CloseWriter( WRITER )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot close connectivity data writer.", FIELDML_INFO%FML_HANDLE, &
        & ERR, ERROR, *999 )
    ENDDO
    
    IF( ALLOCATED( BASIS_INFO ) ) THEN
      DEALLOCATE( BASIS_INFO )
    ENDIF
    IF( ALLOCATED( CONNECTIVITY_INFO ) ) THEN
      DEALLOCATE( CONNECTIVITY_INFO )
    ENDIF
    
    CALL LIST_ITEM_SET( FIELDML_INFO%COMPONENT_HANDLES, COMPONENT_NUMBER, TEMPLATE_HANDLE, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_OUTPUT_ADD_MESH_COMPONENT" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_ADD_MESH_COMPONENT", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_MESH_COMPONENT" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_ADD_MESH_COMPONENT

  !
  !================================================================================================================================
  !
  
  !>Create a parameter evaluator and associated data source containing the nodal dofs for the given field components.
  SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS( FIELDML_INFO, BASE_NAME, DOF_FORMAT, TYPE_HANDLE, FIELD, &
    & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, NODE_DOFS_HANDLE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: DOF_FORMAT !<The name of the format to use when writing dof data.
    INTEGER(INTG), INTENT(IN) :: TYPE_HANDLE !<The FieldML type handle for the field.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: FIELD_COMPONENT_NUMBERS(:) !<The field component numbers for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(INOUT) :: NODE_DOFS_HANDLE !<The handle of the nodal dofs parameter evaluator.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG) :: TYPE_COMPONENT_HANDLE, REAL_1D_HANDLE, NODE_COUNT, INDEX_HANDLE, RESOURCE_HANDLE, SOURCE_HANDLE
    INTEGER(INTG) :: VERSION_NUMBER,COMPONENT_COUNT, I, J, INTERPOLATION_TYPE, GLOBAL_NODE_NUMBER, RANK
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENT_NUMBERS(:)
    INTEGER(INTG), TARGET :: SIZES(2), OFFSETS(2), SINGLE_SIZE
    INTEGER(INTG) :: WRITER, FML_ERR
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: DBUFFER(:)
    REAL(C_DOUBLE) :: DVALUE
    LOGICAL :: NODE_EXISTS
    LOGICAL, ALLOCATABLE :: IS_NODE_BASED(:)
    TYPE(C_PTR) :: SIZE_POINTER
    TYPE(VARYING_STRING) :: ARRAY_LOCATION

    CALL ENTERS( "FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS", ERR, ERROR, *999 )
    
    MESH => FIELD%DECOMPOSITION%MESH
    
    CALL FIELDML_OUTPUT_GET_GENERIC_TYPE( FIELDML_INFO%FML_HANDLE, 1, REAL_1D_HANDLE, .TRUE., ERR, ERROR, *999 )

    COMPONENT_COUNT = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    TYPE_COMPONENT_HANDLE = Fieldml_GetTypeComponentEnsemble( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    NODE_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%NODES_HANDLE )
    
    ALLOCATE( MESH_COMPONENT_NUMBERS( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate mesh component array.", ERR, ERROR, *999 )
    ALLOCATE( IS_NODE_BASED( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate nodal component array.", ERR, ERROR, *999 )

    DO I = 1, COMPONENT_COUNT
      CALL FIELD_COMPONENT_MESH_COMPONENT_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), &
        & MESH_COMPONENT_NUMBERS(I), ERR, ERROR, *999 )
      CALL FIELD_COMPONENT_INTERPOLATION_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), INTERPOLATION_TYPE, &
        & ERR, ERROR, *999 )
        
      IS_NODE_BASED( I ) = ( INTERPOLATION_TYPE == FIELD_NODE_BASED_INTERPOLATION )
    ENDDO

    RESOURCE_HANDLE = Fieldml_CreateHrefDataResource( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.node.resource"), &
      & cchar( DOF_FORMAT ), cchar(BASE_NAME//".dofs.node") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create nodal dofs data resource "//BASE_NAME//".dofs.node.resource", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    NODE_DOFS_HANDLE = Fieldml_CreateParameterEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.node"), REAL_1D_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create nodal dofs parameter set "//BASE_NAME//".dofs.node.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetParameterDataDescription( FIELDML_INFO%FML_HANDLE, NODE_DOFS_HANDLE, FML_DATA_DESCRIPTION_DENSE_ARRAY )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal dofs parameter description for "//BASE_NAME//".dofs.node.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    SIZES( 1 ) = NODE_COUNT
    SIZES( 2 ) = COMPONENT_COUNT
    SINGLE_SIZE = NODE_COUNT
 
    IF( COMPONENT_COUNT == 1 ) THEN
      RANK = 1
      SIZE_POINTER = C_LOC(SINGLE_SIZE)
    ELSE
       RANK = 2
       SIZE_POINTER = C_LOC(SIZES)
    ENDIF

    ARRAY_LOCATION = ARRAY_LOCATION//1
    SOURCE_HANDLE = Fieldml_CreateArrayDataSource( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.node.data"), &
      & RESOURCE_HANDLE, cchar(ARRAY_LOCATION), RANK )
    FML_ERR = Fieldml_SetArrayDataSourceRawSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, SIZE_POINTER )
    FML_ERR = Fieldml_SetArrayDataSourceSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, SIZE_POINTER )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create nodal dofs data source "//BASE_NAME//".dofs.node.data.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    FML_ERR = Fieldml_SetDataSource( FIELDML_INFO%FML_HANDLE, NODE_DOFS_HANDLE, SOURCE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal dofs data source to "//BASE_NAME//".dofs.node.data", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, NODE_DOFS_HANDLE, FIELDML_INFO%NODES_ARGUMENT_HANDLE, &
      & FML_INVALID_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot add layout index for nodal dofs parameter set "//BASE_NAME//".dofs.node.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    IF( TYPE_COMPONENT_HANDLE /= FML_INVALID_HANDLE ) THEN
      TYPE_COMPONENT_HANDLE = FIELDML_OUTPUT_IMPORT_HANDLE( FIELDML_INFO%FML_HANDLE, TYPE_COMPONENT_HANDLE, ERR, ERROR )
      IF(ERR/=0) GOTO 999
      INDEX_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, TYPE_COMPONENT_HANDLE, .TRUE., ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, NODE_DOFS_HANDLE, INDEX_HANDLE, FML_INVALID_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( &
        & "Cannot add component index for nodal dofs parameter set "//BASE_NAME//".dofs.node.", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF

    ALLOCATE( DBUFFER( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate nodal dofs array.", ERR, ERROR, *999 )
    WRITER = Fieldml_OpenArrayWriter( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, REAL_1D_HANDLE, 0, SIZE_POINTER, RANK )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open nodal parameter writer for "//BASE_NAME//".dofs.node.data.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
    OFFSETS(:) = 0
    SIZES(1) = 1
    SIZES(2) = COMPONENT_COUNT
    DO I = 1, NODE_COUNT
      DO J = 1, COMPONENT_COUNT
        DVALUE = 0
        IF( IS_NODE_BASED(J) ) THEN
          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS( MESH, MESH_COMPONENT_NUMBERS(J), I, NODE_EXISTS, GLOBAL_NODE_NUMBER, &
            & ERR, ERROR, *999 )
          IF( NODE_EXISTS ) THEN
            !Default to version 1 of each node derivative (value hardcoded in loop)
            VERSION_NUMBER = 1
            CALL FIELD_PARAMETER_SET_GET_NODE( FIELD, VARIABLE_TYPE, SET_TYPE, VERSION_NUMBER, &
              & NO_GLOBAL_DERIV, I, FIELD_COMPONENT_NUMBERS(J), DVALUE, ERR, ERROR, *999 )
          ENDIF
        ENDIF
        DBUFFER( J ) = DVALUE
      ENDDO
      FML_ERR = Fieldml_WriteDoubleSlab( WRITER, C_LOC(OFFSETS), C_LOC(SIZES), C_LOC(DBUFFER) )
      IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( var_str("I/O error while writing nodal parameter values for ")//BASE_NAME//"("// &
          & TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", err, ERROR, *999 )
      ENDIF
      OFFSETS(1) = OFFSETS(1) + 1
    ENDDO
    FML_ERR = Fieldml_CloseWriter( WRITER )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot close nodal parameter writer for "//BASE_NAME//".dofs.node.data.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    DEALLOCATE( DBUFFER )
    
    DEALLOCATE( MESH_COMPONENT_NUMBERS )
    DEALLOCATE( IS_NODE_BASED )

    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS
  
  !
  !================================================================================================================================
  !
  
  !>Create a parameter evaluator and associated data source containing the element dofs for the given field components.
  SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS( FIELDML_INFO, BASE_NAME, DOF_FORMAT, TYPE_HANDLE, FIELD, &
    & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, ELEMENT_DOFS_HANDLE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: DOF_FORMAT !<The name of the format to use when writing dof data.
    INTEGER(INTG), INTENT(IN) :: TYPE_HANDLE !<The FieldML type handle for the field.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: FIELD_COMPONENT_NUMBERS(:) !<The field component numbers for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(INOUT) :: ELEMENT_DOFS_HANDLE !<The handle of the element dofs parameter evaluator.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: TYPE_COMPONENT_HANDLE, REAL_1D_HANDLE, ELEMENT_COUNT, INDEX_HANDLE, RESOURCE_HANDLE, SOURCE_HANDLE
    INTEGER(INTG) :: COMPONENT_COUNT, I, J, INTERPOLATION_TYPE
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENT_NUMBERS(:)
    INTEGER(INTG) :: WRITER, FML_ERR
    INTEGER(INTG), TARGET :: SIZES(2), OFFSETS(2)
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: DBUFFER(:)
    REAL(C_DOUBLE) :: DVALUE
    LOGICAL, ALLOCATABLE :: IS_ELEMENT_BASED(:)
    TYPE(VARYING_STRING) :: ARRAY_LOCATION

    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS" )
    
    CALL FIELDML_OUTPUT_GET_GENERIC_TYPE( FIELDML_INFO%FML_HANDLE, 1, REAL_1D_HANDLE, .TRUE., ERR, ERROR, *999 )

    COMPONENT_COUNT = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    TYPE_COMPONENT_HANDLE = Fieldml_GetTypeComponentEnsemble( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    
    ELEMENT_COUNT = Fieldml_GetMemberCount( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%ELEMENTS_HANDLE )
    
    ALLOCATE( MESH_COMPONENT_NUMBERS( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate mesh component number array.", ERR, ERROR, *999 )
    ALLOCATE( IS_ELEMENT_BASED( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate element component array.", ERR, ERROR, *999 )

    DO I = 1, COMPONENT_COUNT
      CALL FIELD_COMPONENT_MESH_COMPONENT_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), &
        & MESH_COMPONENT_NUMBERS(I), ERR, ERROR, *999 )
      CALL FIELD_COMPONENT_INTERPOLATION_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), INTERPOLATION_TYPE, &
        & ERR, ERROR, *999 )

      IS_ELEMENT_BASED( I ) = ( INTERPOLATION_TYPE == FIELD_ELEMENT_BASED_INTERPOLATION )
    ENDDO

    RESOURCE_HANDLE = Fieldml_CreateHrefDataResource( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.element.resource"), &
      & cchar( DOF_FORMAT ), cchar(BASE_NAME//".dofs.element") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create element dofs data resource "//BASE_NAME//".dofs.element.resource.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    ELEMENT_DOFS_HANDLE = Fieldml_CreateParameterEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.element"), &
      & REAL_1D_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create element dofs parameter set "//BASE_NAME//".dofs.element.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetParameterDataDescription( FIELDML_INFO%FML_HANDLE, ELEMENT_DOFS_HANDLE, FML_DATA_DESCRIPTION_DENSE_ARRAY )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set element dofs parameter description for "//BASE_NAME//".dofs.element.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    ARRAY_LOCATION = ARRAY_LOCATION//1
    SOURCE_HANDLE = Fieldml_CreateArrayDataSource( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.element.data"), &
      & RESOURCE_HANDLE, cchar(ARRAY_LOCATION), 2 )
    SIZES( 1 ) = ELEMENT_COUNT
    SIZES( 2 ) = COMPONENT_COUNT
    FML_ERR = Fieldml_SetArrayDataSourceRawSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, C_LOC( SIZES ) )
    FML_ERR = Fieldml_SetArrayDataSourceSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, C_LOC( SIZES ) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create element dofs data source "//BASE_NAME//".dofs.element.data.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    FML_ERR = Fieldml_SetDataSource( FIELDML_INFO%FML_HANDLE, ELEMENT_DOFS_HANDLE, SOURCE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal dofs data source for "//BASE_NAME//".dofs.element.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, ELEMENT_DOFS_HANDLE, &
      & FIELDML_INFO%ELEMENTS_ARGUMENT_HANDLE, FML_INVALID_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot add element index for element dofs parameter set "//BASE_NAME//".dofs.element."&
      & , FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    IF( TYPE_COMPONENT_HANDLE /= FML_INVALID_HANDLE ) THEN
      TYPE_COMPONENT_HANDLE = FIELDML_OUTPUT_IMPORT_HANDLE( FIELDML_INFO%FML_HANDLE, TYPE_COMPONENT_HANDLE, ERR, ERROR )
      IF(ERR/=0) GOTO 999
      INDEX_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, TYPE_COMPONENT_HANDLE, .TRUE., ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, ELEMENT_DOFS_HANDLE, TYPE_COMPONENT_HANDLE, &
        & FML_INVALID_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot add component index for element dofs parameter set "//BASE_NAME//&
        & ".dofs.element.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF

    ALLOCATE( DBUFFER( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate element dofs buffer.", ERR, ERROR, *999 )
    WRITER = Fieldml_OpenArrayWriter( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, REAL_1D_HANDLE, 0, C_LOC(SIZES), 2 )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open element parameter writer for "//BASE_NAME//".dofs.element.data.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      
    OFFSETS(:) = 0
    SIZES(1) = 1
    SIZES(2) = COMPONENT_COUNT
    DO I = 1, ELEMENT_COUNT
      DO J = 1, COMPONENT_COUNT
        DVALUE = 0
        IF( IS_ELEMENT_BASED(J) ) THEN
          CALL FIELD_PARAMETER_SET_GET_ELEMENT( FIELD, VARIABLE_TYPE, SET_TYPE, I, &
            & FIELD_COMPONENT_NUMBERS(J), DVALUE, ERR, ERROR, *999 )
        ENDIF
        DBUFFER( J ) = DVALUE
      ENDDO
      FML_ERR = Fieldml_WriteDoubleSlab( WRITER, C_LOC(OFFSETS), C_LOC(SIZES), C_LOC(DBUFFER) )
      IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( var_str("I/O error while writing element parameter values for")//BASE_NAME//"("&
          & // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", err, ERROR, *999 )
      ENDIF
      OFFSETS(1) = OFFSETS(1) + 1
    ENDDO
    FML_ERR = Fieldml_CloseWriter( WRITER )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot close element parameter writer for "//BASE_NAME//".dofs.element.data", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    DEALLOCATE( DBUFFER )
    
    DEALLOCATE( MESH_COMPONENT_NUMBERS )
    DEALLOCATE( IS_ELEMENT_BASED )

    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS
  
  !
  !================================================================================================================================
  !
  
  !>Create a parameter evaluator and associated data source containing the globally constant dofs for the given field components.
  SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS( FIELDML_INFO, BASE_NAME, DOF_FORMAT, TYPE_HANDLE, FIELD, &
    & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, CONSTANT_DOFS_HANDLE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: DOF_FORMAT !<The name of the format to use when writing dof data.
    INTEGER(INTG), INTENT(IN) :: TYPE_HANDLE !<The FieldML type handle for the field.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: FIELD_COMPONENT_NUMBERS(:) !<The field component numbers for which dof evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(INOUT) :: CONSTANT_DOFS_HANDLE !<The handle of the constant dofs parameter evaluator.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: DOFTYPE_HANDLE, TYPE_TYPE, COMPONENT_TYPE, DATA_TYPE, INDEX_HANDLE, RESOURCE_HANDLE, SOURCE_HANDLE
    INTEGER(INTG) :: COMPONENT_COUNT, I, J, INTERPOLATION_TYPE
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENT_NUMBERS(:)
    INTEGER(INTG), TARGET :: OFFSETS(2), SINGLE_SIZE
    TYPE(VARYING_STRING) :: ARRAY_LOCATION
    INTEGER(INTG) :: WRITER, FML_ERR
    REAL(C_DOUBLE), ALLOCATABLE, TARGET :: DBUFFER(:)
    INTEGER(INTG), ALLOCATABLE, TARGET :: IBUFFER(:)
    REAL(C_DOUBLE) :: DVALUE
    INTEGER(INTG) :: IVALUE
    LOGICAL :: IS_REAL
    LOGICAL, ALLOCATABLE :: IS_CONSTANT(:)

    CALL ENTERS( "FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS", ERR, ERROR, *999 )
    
    TYPE_TYPE = Fieldml_GetObjectType(  FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )

    IF( TYPE_TYPE == FHT_ENSEMBLE_TYPE ) THEN
      DOFTYPE_HANDLE = TYPE_HANDLE
      COMPONENT_COUNT = 1
      COMPONENT_TYPE = FML_INVALID_HANDLE
      IS_REAL = .FALSE.
    ELSE
      CALL FIELDML_OUTPUT_GET_GENERIC_TYPE( FIELDML_INFO%FML_HANDLE, 1, DOFTYPE_HANDLE, .TRUE., ERR, ERROR, *999 )
      COMPONENT_COUNT = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
      COMPONENT_TYPE = Fieldml_GetTypeComponentEnsemble( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
      IS_REAL = .TRUE.
    ENDIF
    
    ALLOCATE( MESH_COMPONENT_NUMBERS( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate mesh component array.", ERR, ERROR, *999 )
    ALLOCATE( IS_CONSTANT( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate constant component array.", ERR, ERROR, *999 )

    DO I = 1, COMPONENT_COUNT
      CALL FIELD_COMPONENT_MESH_COMPONENT_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), &
        & MESH_COMPONENT_NUMBERS(I), ERR, ERROR, *999 )
      CALL FIELD_COMPONENT_INTERPOLATION_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), INTERPOLATION_TYPE, &
        & ERR, ERROR, *999 )

      IS_CONSTANT( I ) = ( INTERPOLATION_TYPE == FIELD_CONSTANT_INTERPOLATION )
    ENDDO

    RESOURCE_HANDLE = Fieldml_CreateHrefDataResource( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.constant.resource"), &
      & cchar( DOF_FORMAT ), cchar(BASE_NAME//".dofs.constant") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create constant dofs data resource "//BASE_NAME//".dofs.constant.resource.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    CONSTANT_DOFS_HANDLE = Fieldml_CreateParameterEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.constant"), &
      & DOFTYPE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create constant dofs parameter set "//BASE_NAME//".dofs.constant.", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetParameterDataDescription( FIELDML_INFO%FML_HANDLE, CONSTANT_DOFS_HANDLE, FML_DATA_DESCRIPTION_DENSE_ARRAY )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set constant dofs parameter description for "//BASE_NAME//".dofs.constant", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    ARRAY_LOCATION = ARRAY_LOCATION//1
    SOURCE_HANDLE = Fieldml_CreateArrayDataSource( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.element.data"), &
      & RESOURCE_HANDLE, cchar(ARRAY_LOCATION), 1 )
    SINGLE_SIZE = COMPONENT_COUNT
    FML_ERR = Fieldml_SetArrayDataSourceRawSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, C_LOC(SINGLE_SIZE) )
    FML_ERR = Fieldml_SetArrayDataSourceSizes( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, C_LOC(SINGLE_SIZE) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create constant dofs data source "//BASE_NAME//".dofs.constant.data", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    FML_ERR = Fieldml_SetDataSource( FIELDML_INFO%FML_HANDLE, CONSTANT_DOFS_HANDLE, SOURCE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal dofs data source for "//BASE_NAME//".dofs.constant", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    IF( COMPONENT_TYPE /= FML_INVALID_HANDLE ) THEN
      COMPONENT_TYPE = FIELDML_OUTPUT_IMPORT_HANDLE( FIELDML_INFO%FML_HANDLE, COMPONENT_TYPE, ERR, ERROR )
      IF(ERR/=0) GOTO 999
      INDEX_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, COMPONENT_TYPE, .TRUE., ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_AddDenseIndexEvaluator( FIELDML_INFO%FML_HANDLE, CONSTANT_DOFS_HANDLE, INDEX_HANDLE, FML_INVALID_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot add component index for constant dofs parameter set "//BASE_NAME//&
        & ".dofs.constant", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF

    WRITER = Fieldml_OpenArrayWriter( FIELDML_INFO%FML_HANDLE, SOURCE_HANDLE, DOFTYPE_HANDLE, 0, C_LOC(SINGLE_SIZE), 1 )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot open constant parameter writer for "//BASE_NAME//".dofs.constant.data", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    CALL FIELD_DATA_TYPE_GET( FIELD, VARIABLE_TYPE, DATA_TYPE, ERR, ERROR, *999 )
    IF( DATA_TYPE == FIELD_INTG_TYPE ) THEN
      IS_REAL = .false.
    ELSEIF( DATA_TYPE == FIELD_DP_TYPE ) THEN
      IS_REAL = .true.
    ENDIF
    
    OFFSETS(:) = 0
    SINGLE_SIZE = COMPONENT_COUNT
    IF( IS_REAL ) THEN
      ALLOCATE( DBUFFER( COMPONENT_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate constant dofs buffer.", ERR, ERROR, *999 )
      DO J = 1, COMPONENT_COUNT
        DVALUE = 0
        IF( IS_CONSTANT(J) ) THEN
          CALL FIELD_PARAMETER_SET_GET_CONSTANT( FIELD, VARIABLE_TYPE, SET_TYPE, &
            & FIELD_COMPONENT_NUMBERS(J), DVALUE, ERR, ERROR, *999 )
        ENDIF
        DBUFFER( J ) = DVALUE
      ENDDO
      FML_ERR = Fieldml_WriteDoubleSlab( WRITER, C_LOC(OFFSETS), C_LOC(SINGLE_SIZE), C_LOC(DBUFFER) )
      IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( var_str("I/O error while writing constant parameter values for ")//BASE_NAME//"(" &
        & // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", err, ERROR, *999)
      ENDIF
      FML_ERR = Fieldml_CloseWriter( WRITER )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot close constant parameter writer for "//BASE_NAME//".dofs.constant.data", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      DEALLOCATE( DBUFFER )
    ELSE
      ALLOCATE( IBUFFER( COMPONENT_COUNT ), STAT = ERR )
      IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate constant dofs buffer.", ERR, ERROR, *999 )
      DO J = 1, COMPONENT_COUNT
        IVALUE = 0
        IF( IS_CONSTANT(J) ) THEN
          CALL FIELD_PARAMETER_SET_GET_CONSTANT( FIELD, VARIABLE_TYPE, SET_TYPE, &
            & FIELD_COMPONENT_NUMBERS(J), IVALUE, ERR, ERROR, *999 )
        ENDIF
        IBUFFER( J ) = IVALUE
      ENDDO
      FML_ERR = Fieldml_WriteIntSlab( WRITER, C_LOC(OFFSETS), C_LOC(SINGLE_SIZE), C_LOC(IBUFFER) )
      IF( FML_ERR /= FML_ERR_NO_ERROR ) THEN
        CALL FLAG_ERROR( var_str("I/O while writing constant parameter values for ")//BASE_NAME//"(" &
          & // TRIM(NUMBER_TO_VSTRING(FML_ERR,"*",ERR,ERROR)) //").", ERR, ERROR, *999 )
      ENDIF
      FML_ERR = Fieldml_CloseWriter( WRITER )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot close constant parameter writer for "//BASE_NAME//".dofs.constant.data", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      DEALLOCATE( IBUFFER )
    ENDIF
    
    DEALLOCATE( MESH_COMPONENT_NUMBERS )
    DEALLOCATE( IS_CONSTANT )

    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS
  
  !
  !================================================================================================================================
  !

  !>Initialize the given FieldML parsing state for use with the given mesh
  SUBROUTINE FIELDML_OUTPUT_INITIALISE_INFO( MESH, LOCATION, BASE_NAME, CONNECTIVITY_FORMAT, FIELDML_INFO, ERR, ERROR, * )
    !Argument variables
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: MESH !<The mesh with which the FieldML document is associated.
    TYPE(VARYING_STRING), INTENT(IN) :: LOCATION !<The location of the FieldML file. Data resources will be created here.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: CONNECTIVITY_FORMAT !<The name of the format to use when writing connectivity data.
    TYPE(FIELDML_IO_TYPE), INTENT(OUT) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG) :: COMPONENT_COUNT, I, NODE_COUNT, ELEMENT_COUNT, DIMENSIONS
    INTEGER(INTG) :: REAL_1D_HANDLE, XI_COMPONENT_HANDLE, FML_ERR, SHAPE_HANDLE
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(VARYING_STRING) :: SHAPE_NAME

    CALL ENTERS( "FIELDML_OUTPUT_INITIALISE_INFO", ERR, ERROR, *999 )
    
    REGION => MESH%REGION
    
    DIMENSIONS = MESH%NUMBER_OF_DIMENSIONS
    
    CALL FIELDML_IO_INITIALISE( FIELDML_INFO, .TRUE., ERR, ERROR, *999 )
    
    FIELDML_INFO%FML_HANDLE = Fieldml_Create( cchar(LOCATION), cchar(BASE_NAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create fieldml handle for "//BASE_NAME//" at "//LOCATION//".", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    NULLIFY( NODES )
    CALL REGION_NODES_GET( REGION, NODES, ERR, ERROR, *999 )
    CALL NODES_NUMBER_OF_NODES_GET( NODES, NODE_COUNT, ERR, ERROR, *999 )

    FIELDML_INFO%NODES_HANDLE = Fieldml_CreateEnsembleType( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".nodes") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh nodes ensemble "//BASE_NAME//".nodes", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetEnsembleMembersRange( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%NODES_HANDLE, 1, NODE_COUNT, 1 )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set mesh nodes ensemble bounds for "//BASE_NAME//".nodes", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FIELDML_INFO%NODES_ARGUMENT_HANDLE = Fieldml_CreateArgumentEvaluator( FIELDML_INFO%FML_HANDLE, &
      & cchar(BASE_NAME//".nodes.argument"), FIELDML_INFO%NODES_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh nodes variable "//BASE_NAME//".nodes.argument", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    CALL MESH_NUMBER_OF_ELEMENTS_GET( MESH, ELEMENT_COUNT, ERR, ERROR, *999 )

    FIELDML_INFO%MESH_HANDLE = Fieldml_CreateMeshType( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".mesh") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh type "//BASE_NAME//".mesh", FIELDML_INFO%FML_HANDLE, &
      & ERR, ERROR, *999 )

    FIELDML_INFO%ELEMENTS_HANDLE = Fieldml_CreateMeshElementsType( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%MESH_HANDLE, &
      & "element"//C_NULL_CHAR )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh elements type for "//BASE_NAME//".mesh", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FML_ERR = Fieldml_SetEnsembleMembersRange( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%ELEMENTS_HANDLE, 1, ELEMENT_COUNT, 1 )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set mesh type element count for "//BASE_NAME//".mesh", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FIELDML_INFO%XI_HANDLE = Fieldml_CreateMeshChartType( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%MESH_HANDLE, "xi"//C_NULL_CHAR )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh chart type for "//BASE_NAME//".mesh", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    XI_COMPONENT_HANDLE = Fieldml_CreateContinuousTypeComponents( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%XI_HANDLE, &
      & cchar(BASE_NAME//".mesh.xi.component"), DIMENSIONS )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh chart components for "//BASE_NAME//".mesh", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FML_ERR = Fieldml_CreateArgumentEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".mesh.argument"), &
      & FIELDML_INFO%MESH_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create mesh variable "//BASE_NAME//".mesh.argument", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    FIELDML_INFO%XI_ARGUMENT_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".mesh.argument.xi") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get mesh xi variable for "//BASE_NAME//".mesh", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    FIELDML_INFO%ELEMENTS_ARGUMENT_HANDLE = Fieldml_GetObjectByName( FIELDML_INFO%FML_HANDLE, &
      & cchar(BASE_NAME//".mesh.argument.element") )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot get mesh element variable for "//BASE_NAME//".mesh", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    
    CALL FIELDML_OUTPUT_GET_GENERIC_TYPE( FIELDML_INFO%FML_HANDLE, 1, REAL_1D_HANDLE, .TRUE., ERR, ERROR, *999 )
 
    !TODO Some of these may end up being unused. Should use deferred assignment.
    FIELDML_INFO%NODE_DOFS_HANDLE = Fieldml_CreateArgumentEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME//".dofs.node"), &
      & REAL_1D_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create nodal dofs variable "//BASE_NAME//".dofs.node", &
      & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
!    fieldmlInfo%elementDofsHandle = Fieldml_CreateArgumentEvaluator( fieldmlInfo%FML_HANDLE, cchar(baseName//".dofs.element"), & 
!      & real1DHandle )
!    CALL FieldmlUtilCheckFieldmlError( "Cannot create element dofs variable "//".dofs.element", &
!      & fieldmlInfo, err, errorString, *999 )
!    fieldmlInfo%constantDofsHandle = Fieldml_CreateArgumentEvaluator( fieldmlInfo%FML_HANDLE, cchar(baseName//".dofs.constant"), & 
!      & real1DHandle )
!    CALL FieldmlUtilCheckFieldmlError( "Cannot create constant dofs variable "//".dofs.constant", &
!      & fieldmlInfo, err, errorString, *999 )

    CALL MESH_NUMBER_OF_COMPONENTS_GET( MESH, COMPONENT_COUNT, ERR, ERROR, *999 )
    DO I = 1, COMPONENT_COUNT
      NULLIFY( MESH_ELEMENTS )
      CALL LIST_ITEM_ADD( FIELDML_INFO%COMPONENT_HANDLES, FML_INVALID_HANDLE, ERR, ERROR, *999 )
      CALL MESH_TOPOLOGY_ELEMENTS_GET( MESH, I, MESH_ELEMENTS, ERR, ERROR, *999 )
      CALL FIELDML_OUTPUT_ADD_MESH_COMPONENT( FIELDML_INFO, BASE_NAME, CONNECTIVITY_FORMAT, I, MESH_ELEMENTS, &
        & ERR, ERROR, *999 )
    ENDDO

    
    !TODO Proper shape assignment.
    IF( DIMENSIONS == 1 ) THEN
      SHAPE_NAME = "shape.unit.line"
    ELSE IF( DIMENSIONS == 2 ) THEN
      SHAPE_NAME = "shape.unit.square"
    ELSE
      SHAPE_NAME = "shape.unit.cube"
    ENDIF

    SHAPE_HANDLE = FIELDML_OUTPUT_IMPORT( FIELDML_INFO, SHAPE_NAME, ERR, ERROR )

    IF(ERR/=0) GOTO 999
    FML_ERR = Fieldml_SetMeshShapes( FIELDML_INFO%FML_HANDLE, FIELDML_INFO%MESH_HANDLE, SHAPE_HANDLE )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set mesh type element shape.", FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_OUTPUT_INITIALISE_INFO" )
    RETURN

999 CALL ERRORS( "FIELDML_OUTPUT_INITIALISE_INFO", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_INITIALISE_INFO" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_INITIALISE_INFO

  !
  !================================================================================================================================
  !
  
  !> Add the components of the given field to the given FieldML evaluator, creating component templates as needed.
  SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_COMPONENTS( FIELDML_INFO, TYPE_HANDLE, BASE_NAME, DOF_FORMAT, FIELD, &
    & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    INTEGER(INTG), INTENT(IN) :: TYPE_HANDLE !<The FieldML type handle for the field.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: DOF_FORMAT !<The name of the format to use when writing dof data.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field for which dof components are to be created.
    INTEGER(INTG), INTENT(IN) :: FIELD_COMPONENT_NUMBERS(:) !<The field component numbers for which evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG) :: FIELD_HANDLE, COMPONENT_HANDLE, NODAL_DOFS_HANDLE, ELEMENT_DOFS_HANDLE, CONSTANT_DOFS_HANDLE, INDEX_HANDLE
    INTEGER(INTG) :: COMPONENT_COUNT, I, MESH_COMPONENT_NUMBER, INTERPOLATION_TYPE, FML_ERR
    INTEGER(INTG), ALLOCATABLE, TARGET :: COMPONENT_EVALUATORS(:)
  
    CALL ENTERS( "FIELDML_OUTPUT_ADD_FIELD_COMPONENTS", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_OUT( FIELDML_INFO, ERR, ERROR, *999 )
    
    MESH => FIELD%DECOMPOSITION%MESH

    COMPONENT_HANDLE = Fieldml_GetTypeComponentEnsemble( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    COMPONENT_COUNT = Fieldml_GetTypeComponentCount( FIELDML_INFO%FML_HANDLE, TYPE_HANDLE )
    ALLOCATE( COMPONENT_EVALUATORS( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate component evaluators array.", ERR, ERROR, *999 )

    IF( SIZE( FIELD_COMPONENT_NUMBERS ) /= COMPONENT_COUNT ) THEN
      CALL FLAG_ERROR( var_str("Fieldml Component count ")//SIZE( FIELD_COMPONENT_NUMBERS )//&
        & " must match value type component count "//COMPONENT_COUNT//".", ERR, ERROR, *999 )
    ENDIF

    NODAL_DOFS_HANDLE = FML_INVALID_HANDLE
    ELEMENT_DOFS_HANDLE = FML_INVALID_HANDLE
    CONSTANT_DOFS_HANDLE = FML_INVALID_HANDLE
    !TODO Other types of interpolation not yet supported.
    DO I = 1, COMPONENT_COUNT
      CALL FIELD_COMPONENT_INTERPOLATION_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), INTERPOLATION_TYPE, &
        & ERR, ERROR, *999 )
        
      IF( INTERPOLATION_TYPE == FIELD_NODE_BASED_INTERPOLATION ) THEN
        IF( NODAL_DOFS_HANDLE == FML_INVALID_HANDLE ) THEN
          CALL FIELDML_OUTPUT_ADD_FIELD_NODE_DOFS( FIELDML_INFO, BASE_NAME, DOF_FORMAT, TYPE_HANDLE, FIELD, &
          & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, NODAL_DOFS_HANDLE, ERR, ERROR, *999 )
        ENDIF
        CALL FIELD_COMPONENT_MESH_COMPONENT_GET( FIELD, VARIABLE_TYPE, FIELD_COMPONENT_NUMBERS(I), &
          & MESH_COMPONENT_NUMBER, ERR, ERROR, *999 )
        CALL LIST_ITEM_GET( FIELDML_INFO%COMPONENT_HANDLES, MESH_COMPONENT_NUMBER, COMPONENT_EVALUATORS( I ), &
          & ERR, ERROR, *999 )
      ELSEIF( INTERPOLATION_TYPE == FIELD_ELEMENT_BASED_INTERPOLATION ) THEN
        IF( ELEMENT_DOFS_HANDLE == FML_INVALID_HANDLE ) THEN
          CALL FIELDML_OUTPUT_ADD_FIELD_ELEMENT_DOFS( FIELDML_INFO, BASE_NAME, DOF_FORMAT, TYPE_HANDLE, FIELD, &
            & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, ELEMENT_DOFS_HANDLE, ERR, ERROR, *999 )
        ENDIF
        COMPONENT_EVALUATORS( I ) = ELEMENT_DOFS_HANDLE
      ELSEIF( INTERPOLATION_TYPE == FIELD_CONSTANT_INTERPOLATION ) THEN
        IF( CONSTANT_DOFS_HANDLE == FML_INVALID_HANDLE ) THEN
          CALL FIELDML_OUTPUT_ADD_FIELD_CONSTANT_DOFS( FIELDML_INFO, BASE_NAME, DOF_FORMAT, TYPE_HANDLE, FIELD, &
            & FIELD_COMPONENT_NUMBERS, VARIABLE_TYPE, SET_TYPE, CONSTANT_DOFS_HANDLE, ERR, ERROR, *999 )
        ENDIF
        COMPONENT_EVALUATORS( I ) = CONSTANT_DOFS_HANDLE
      ENDIF
    ENDDO
    
    IF( COMPONENT_HANDLE /= FML_INVALID_HANDLE ) THEN
      FIELD_HANDLE = Fieldml_CreateAggregateEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME), TYPE_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create field aggregate evaluator "//BASE_NAME, &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      INDEX_HANDLE = FIELDML_OUTPUT_GET_TYPE_ARGUMENT_HANDLE( FIELDML_INFO, COMPONENT_HANDLE, .TRUE., ERR, ERROR )
      IF(ERR/=0) GOTO 999
      FML_ERR = Fieldml_SetIndexEvaluator( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE, 1, INDEX_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set index evaluator for aggregate evaluator "//BASE_NAME//".", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )

      DO I = 1, COMPONENT_COUNT
        FML_ERR = Fieldml_SetEvaluator( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE, I, COMPONENT_EVALUATORS( I ) )
        CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal evaluator for aggregate evaluator "//BASE_NAME//".", &
          & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
      ENDDO
    ELSE
      FIELD_HANDLE = Fieldml_CreateReferenceEvaluator( FIELDML_INFO%FML_HANDLE, cchar(BASE_NAME), COMPONENT_EVALUATORS( 1 ) )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot create reference evaluator for field "//BASE_NAME, &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF

    IF( NODAL_DOFS_HANDLE /= FML_INVALID_HANDLE ) THEN
      FML_ERR = Fieldml_SetBind( FIELDML_INFO%FML_HANDLE, FIELD_HANDLE, FIELDML_INFO%NODE_DOFS_HANDLE, NODAL_DOFS_HANDLE )
      CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Cannot set nodal dofs bind for field "//BASE_NAME//" with interpolated elements", &
        & FIELDML_INFO%FML_HANDLE, ERR, ERROR, *999 )
    ENDIF
!    IF( elementDofsHandle /= FML_INVALID_HANDLE ) THEN
!      fmlErr = Fieldml_SetBind( fieldmlInfo%FML_HANDLE, fieldHandle, fieldmlInfo%elementDofsHandle, elementDofsHandle )
!      CALL FieldmlUtilCheckFieldmlError( "Cannot set element dofs bind for field with constant elements", fieldmlInfo, &
!  &err, errorString, *999 )
!    ENDIF
!    IF( constantDofsHandle /= FML_INVALID_HANDLE ) THEN
!      fmlErr = Fieldml_SetBind( fieldmlInfo%FML_HANDLE, fieldHandle, fieldmlInfo%constantDofsHandle, constantDofsHandle )
!      CALL FieldmlUtilCheckFieldmlError( "Cannot set constant dofs bind for field with constant value", fieldmlInfo, &
!  &err, errorString, *999 )
!    ENDIF


    DEALLOCATE( COMPONENT_EVALUATORS )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_COMPONENTS" )
    RETURN
999 DEALLOCATE( COMPONENT_EVALUATORS )
    CALL ERRORS( "FIELDML_OUTPUT_ADD_FIELD_COMPONENTS", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_COMPONENTS" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_COMPONENTS

  !
  !================================================================================================================================
  !

  !>Add the given field to the given FieldML document. The field's type will be determined by FieldmlUtilGetValueType. \see Fieldml_Util_Routines::FieldmlUtilGetValueType
  SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_NO_TYPE( FIELDML_INFO, BASE_NAME, DOF_FORMAT, FIELD, VARIABLE_TYPE, SET_TYPE, &
    & ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: DOF_FORMAT !<The name of the format to use when writing dof data.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field for which evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: TYPE_HANDLE
    
    CALL ENTERS( "FIELDML_OUTPUT_ADD_FIELD_NO_TYPE", ERR, ERROR, *999 )
    
    CALL FIELDML_ASSERT_IS_OUT( FIELDML_INFO, ERR, ERROR, *999 )

    CALL FIELDML_OUTPUT_GET_VALUE_TYPE( FIELDML_INFO%FML_HANDLE, FIELD, VARIABLE_TYPE, .TRUE., TYPE_HANDLE, ERR, ERROR, *999 )

    CALL FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE( FIELDML_INFO, BASE_NAME, DOF_FORMAT, FIELD, VARIABLE_TYPE, SET_TYPE, TYPE_HANDLE, &
      & ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_NO_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_ADD_FIELD_NO_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_NO_TYPE" )
    RETURN 1

  END SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_NO_TYPE

  !
  !================================================================================================================================
  !

  !>Add the given field to the given FieldML document using the given FieldML type.
  SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE( FIELDML_INFO, BASE_NAME, DOF_FORMAT, FIELD, VARIABLE_TYPE, SET_TYPE, &
    & TYPE_HANDLE, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: BASE_NAME !<The root name of the basis evaluator.
    TYPE(VARYING_STRING), INTENT(IN) :: DOF_FORMAT !<The name of the format to use when writing dof data.
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: FIELD !<The field for which evaluators are to be created.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The OpenCMISS variable type to generate dofs for.
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The parameter set type.
    INTEGER(INTG), INTENT(IN) :: TYPE_HANDLE !<The FieldML type handle for the field.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.

    !Locals
    INTEGER(INTG) :: I, COMPONENT_COUNT
    INTEGER(INTG), ALLOCATABLE :: FIELD_COMPONENT_NUMBERS(:)
    TYPE(MESH_TYPE), POINTER :: MESH
    
    CALL ENTERS( "FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE", ERR, ERROR, *999 )

    CALL FIELDML_ASSERT_IS_OUT( FIELDML_INFO, ERR, ERROR, *999 )

    MESH => FIELD%DECOMPOSITION%MESH

    IF( TYPE_HANDLE == FML_INVALID_HANDLE ) THEN
      CALL FLAG_ERROR( var_str("Cannot get value type for field ")//BASE_NAME//".", ERR, ERROR, *999 )
    ENDIF
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET( FIELD, VARIABLE_TYPE, COMPONENT_COUNT, ERR, ERROR, *999 )
    
    ALLOCATE( FIELD_COMPONENT_NUMBERS( COMPONENT_COUNT ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( "Could not allocate component numbers array.", ERR, ERROR, *999 )
    DO I = 1, COMPONENT_COUNT
      FIELD_COMPONENT_NUMBERS(I) = I
    ENDDO

    CALL FIELDML_OUTPUT_ADD_FIELD_COMPONENTS( FIELDML_INFO, TYPE_HANDLE, BASE_NAME, DOF_FORMAT, FIELD, FIELD_COMPONENT_NUMBERS, &
      & VARIABLE_TYPE, SET_TYPE, ERR, ERROR, *999 )
    
    DEALLOCATE( FIELD_COMPONENT_NUMBERS )

    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE" )
    RETURN 1
    
  END SUBROUTINE FIELDML_OUTPUT_ADD_FIELD_WITH_TYPE

  !
  !================================================================================================================================
  !

  !>Write the given FieldML document to the given file.
  SUBROUTINE FIELDML_OUTPUT_WRITE( FIELDML_INFO, FILENAME, ERR, ERROR, * )
    !Argument variables
    TYPE(FIELDML_IO_TYPE), INTENT(IN) :: FIELDML_INFO !<The FieldML parsing state.
    TYPE(VARYING_STRING), INTENT(IN) :: FILENAME !<The file to write the FieldML document to.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    
    !Locals
    INTEGER(INTG) :: FML_ERR

    CALL ENTERS( "FIELDML_OUTPUT_WRITE", ERR, ERROR, *999 )

    CALL FIELDML_ASSERT_IS_OUT( FIELDML_INFO, ERR, ERROR, *999 )

    FML_ERR = Fieldml_WriteFile( FIELDML_INFO%FML_HANDLE, cchar(FILENAME) )
    CALL FIELDML_UTIL_CHECK_FIELDML_ERROR( "Error writing fieldml file "//FILENAME//".", FIELDML_INFO%FML_HANDLE, &
      & ERR, ERROR, *999 )

    CALL EXITS( "FIELDML_OUTPUT_WRITE" )
    RETURN
999 CALL ERRORS( "FIELDML_OUTPUT_WRITE", ERR, ERROR )
    CALL EXITS( "FIELDML_OUTPUT_WRITE" )
    RETURN 1
  
  END SUBROUTINE

  !
  !================================================================================================================================
  !

END MODULE FIELDML_OUTPUT_ROUTINES
