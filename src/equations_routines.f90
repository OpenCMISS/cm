!> \file
!> \author Chris Bradley
!> \brief This module handles all equations routines.
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

!> This module handles all equations routines.
MODULE EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE


  !> \addtogroup EQUATIONS_ROUTINES_OutputTypes EQUATIONS_ROUTINES::OutputTypes
  !> \brief The equations output types
  !> \see EQUATIONS_ROUTINES,OPENCMISS_EquationsConstants
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_NO_OUTPUT=0 !<No output. \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_TIMING_OUTPUT=1 !<Timing information output. \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output. \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output. \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_NODAL_MATRIX_OUTPUT=4 !<All below and nodal matrices output. \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
  !>@}

  !> \addtogroup EQUATIONS_ROUTINES_SparsityTypes EQUATIONS_ROUTINES::SparsityTypes
  !> \brief Equations matrices sparsity types
  !> \see EQUATIONS_ROUTINES,OPENCMISS_EquationsSparsityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations. \see EQUATIONS_ROUTINES_SparsityTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. \see EQUATIONS_ROUTINES_SparsityTypes,EQUATIONS_ROUTINES
 !>@}
 
  !> \addtogroup EQUATIONS_ROUTINES_LumpingTypes EQUATIONS_ROUTINES::LumpingTypes
  !> \brief Equations matrices lumping types
  !> \see EQUATIONS_ROUTINES,OPENCMISS_EquationsLumpingTypes
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_UNLUMPED_MATRICES=1 !<The equations matrices are not lumped. \see EQUATIONS_ROUTINES_LumpingTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: EQUATIONS_LUMPED_MATRICES=2 !<The equations matrices are "mass" lumped. \see EQUATIONS_ROUTINES_LumpingTypes,EQUATIONS_ROUTINES
 !>@}
 
  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATIONS_NO_OUTPUT,EQUATIONS_TIMING_OUTPUT,EQUATIONS_MATRIX_OUTPUT,EQUATIONS_ELEMENT_MATRIX_OUTPUT

  PUBLIC EQUATIONS_NODAL_MATRIX_OUTPUT

  PUBLIC EQUATIONS_SPARSE_MATRICES,EQUATIONS_FULL_MATRICES

  PUBLIC EQUATIONS_UNLUMPED_MATRICES,EQUATIONS_LUMPED_MATRICES
  
  PUBLIC EQUATIONS_CREATE_START,EQUATIONS_CREATE_FINISH

  PUBLIC EQUATIONS_DESTROY

  PUBLIC EQUATIONS_INITIALISE,EQUATIONS_FINALISE

  PUBLIC EQUATIONS_LINEARITY_TYPE_GET,EQUATIONS_LINEARITY_TYPE_SET

  PUBLIC EQUATIONS_LUMPING_TYPE_GET,EQUATIONS_LUMPING_TYPE_SET

  PUBLIC EQUATIONS_OUTPUT_TYPE_GET,EQUATIONS_OUTPUT_TYPE_SET

  PUBLIC EQUATIONS_SPARSITY_TYPE_GET,EQUATIONS_SPARSITY_TYPE_SET

  PUBLIC EQUATIONS_TIME_DEPENDENCE_TYPE_GET,EQUATIONS_TIME_DEPENDENCE_TYPE_SET

  PUBLIC EQUATIONS_SET_EQUATIONS_GET

  PUBLIC Equations_DerivedVariableGet

  PUBLIC Equations_NumberOfLinearMatricesGet

  PUBLIC Equations_NumberOfJacobianMatricesGet

  PUBLIC Equations_NumberOfDynamicMatricesGet

  PUBLIC Equations_LinearMatrixGet

  PUBLIC Equations_JacobianMatrixGet

  PUBLIC Equations_DynamicMatrixGet

  PUBLIC Equations_DynamicMatrixGetByType

  PUBLIC Equations_DynamicMatrixTypeGet

  PUBLIC Equations_RhsVectorGet

  PUBLIC Equations_ResidualVectorGet

  PUBLIC Equations_ResidualNumberOfVariablesGet

  PUBLIC Equations_ResidualVariablesGet

  PUBLIC Equations_SourceVectorGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of equations.
  SUBROUTINE EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FlagError("Equations have already been finished.",ERR,ERROR,*999)        
      ELSE
        !Set the finished flag
        EQUATIONS%EQUATIONS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of equations for the equation set.
  SUBROUTINE EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create equations for
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) THEN
        CALL FlagError("Equations are already associated for the equations set.",ERR,ERROR,*999)        
      ELSE
        IF(ASSOCIATED(EQUATIONS)) THEN
          CALL FlagError("Equations is already associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the equations
          CALL EQUATIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          !Return the pointer
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_CREATE_START")
    RETURN
999 ERRORSEXITS("EQUATIONS_CREATE_START",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys equations
  SUBROUTINE EQUATIONS_DESTROY(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      CALL EQUATIONS_FINALISE(EQUATIONS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_DESTROY",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the equations and deallocate all memory.
  SUBROUTINE EQUATIONS_FINALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EQUATIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MAPPING)) CALL EQUATIONS_MAPPING_DESTROY(EQUATIONS%EQUATIONS_MAPPING,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MATRICES)) CALL EQUATIONS_MATRICES_DESTROY(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS)
    ENDIF
       
    EXITS("EQUATIONS_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the equations for an equations set.
  SUBROUTINE EQUATIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("EQUATIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) THEN
        CALL FlagError("Equations is already associated for this equations set.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%EQUATIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate equations.",ERR,ERROR,*999)
        EQUATIONS_SET%EQUATIONS%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%EQUATIONS%LINEARITY=EQUATIONS_LINEAR
        EQUATIONS_SET%EQUATIONS%TIME_DEPENDENCE=EQUATIONS_STATIC
        EQUATIONS_SET%EQUATIONS%OUTPUT_TYPE=EQUATIONS_NO_OUTPUT
        EQUATIONS_SET%EQUATIONS%SPARSITY_TYPE=EQUATIONS_SPARSE_MATRICES
        EQUATIONS_SET%EQUATIONS%LUMPING_TYPE=EQUATIONS_UNLUMPED_MATRICES
        NULLIFY(EQUATIONS_SET%EQUATIONS%INTERPOLATION)
        NULLIFY(EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING)
        NULLIFY(EQUATIONS_SET%EQUATIONS%EQUATIONS_MATRICES)
        EQUATIONS_SET%EQUATIONS%EQUATIONS_FINISHED=.FALSE.
        CALL EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS_SET%EQUATIONS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*998)
    ENDIF
       
    EXITS("EQUATIONS_INITIALISE")
    RETURN
999 CALL EQUATIONS_FINALISE(EQUATIONS_SET%EQUATIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_INITIALISE
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the interpolation information for equations and deallocates all memory
  SUBROUTINE EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS_INTERPOLATION,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_INTERPOLATION_TYPE), POINTER :: EQUATIONS_INTERPOLATION !<A pointer to the equations interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("EQUATIONS_INTERPOLATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_INTERPOLATION)) THEN
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%FIBRE_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%INDEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(EQUATIONS_INTERPOLATION%GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(EQUATIONS_INTERPOLATION%DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(EQUATIONS_INTERPOLATION%INDEPENDENT_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(EQUATIONS_INTERPOLATION%FIBRE_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(EQUATIONS_INTERPOLATION%MATERIALS_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(EQUATIONS_INTERPOLATION%SOURCE_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_PHYSICAL_POINTS_FINALISE(EQUATIONS_INTERPOLATION%DEPENDENT_PHYSICAL_POINT,ERR,ERROR,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(EQUATIONS_INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(EQUATIONS_INTERPOLATION%INDEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(EQUATIONS_INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
      CALL Field_InterpolatedPointsMetricsFinalise(EQUATIONS_INTERPOLATION%FIBRE_INTERP_POINT_METRICS,ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_INTERPOLATION)
    ENDIF
       
    EXITS("EQUATIONS_INTERPOLATION_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_INTERPOLATION_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_INTERPOLATION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interpolation information for equations
  SUBROUTINE EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<The pointer to the equations to initialise the interpolation for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    ENTERS("EQUATIONS_INTERPOLATION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        IF(ASSOCIATED(EQUATIONS%INTERPOLATION)) THEN
          CALL FlagError("Interpolation is already associated for these equations.",ERR,ERROR,*998)
        ELSE
          ALLOCATE(EQUATIONS%INTERPOLATION,STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate equations interpolation",ERR,ERROR,*999)
          EQUATIONS%INTERPOLATION%EQUATIONS=>EQUATIONS
          NULLIFY(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_PHYSICAL_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS)
          NULLIFY(EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT_METRICS)
          NULLIFY(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS)
          NULLIFY(EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT_METRICS)
          
          EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          EQUATIONS%INTERPOLATION%FIBRE_FIELD=>EQUATIONS_SET%GEOMETRY%FIBRE_FIELD
          EQUATIONS%INTERPOLATION%DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
            EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
          ELSE
            NULLIFY(EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD)
          ENDIF
          IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
            EQUATIONS%INTERPOLATION%MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          ELSE
            NULLIFY(EQUATIONS%INTERPOLATION%MATERIALS_FIELD)
          ENDIF
          IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
            EQUATIONS%INTERPOLATION%SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
          ELSE
            NULLIFY(EQUATIONS%INTERPOLATION%SOURCE_FIELD)
          ENDIF

          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD, &
            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS, &
            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
          CALL Field_InterpolatedPointsMetricsInitialise(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT, &
            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%DEPENDENT_FIELD, &
            & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS, &
            & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
!           CALL FIELD_PHYSICAL_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT, &
!             & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT,EQUATIONS%INTERPOLATION%DEPENDENT_PHYSICAL_POINT, &
!             & ERR,ERROR,*999)
          IF(EQUATIONS%INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
            & EQUATIONS%INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE.OR. &
            & EQUATIONS%INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
            CALL Field_InterpolatedPointsMetricsInitialise(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT, &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%FIBRE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%FIBRE_FIELD, &
              & EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS,  &
              & EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT,ERR,ERROR,*999)
            CALL Field_InterpolatedPointsMetricsInitialise(EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT,  &
              & EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD, &
              & EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_PARAMETERS,  &
              & EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT,ERR,ERROR,*999)
            IF(EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
              & EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
              CALL Field_InterpolatedPointsMetricsInitialise(EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT,  &
                &  EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
            END IF
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%MATERIALS_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%MATERIALS_FIELD, &
              & EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS,  &
              & EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%SOURCE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%SOURCE_FIELD, &
              & EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS, &
              & EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT,ERR,ERROR,*999)
          ENDIF
          
        ENDIF
      ELSE
        CALL FlagError("Equations equation set is not associated",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated",ERR,ERROR,*998)
    ENDIF
       
    EXITS("EQUATIONS_INTERPOLATION_INITIALISE")
    RETURN
999 CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("EQUATIONS_INTERPOLATION_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_INTERPOLATION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the linearity type for equations.
  SUBROUTINE EQUATIONS_LINEARITY_TYPE_GET(EQUATIONS,LINEARITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to get the linearity for
    INTEGER(INTG), INTENT(OUT) :: LINEARITY_TYPE !<On exit, the linearity type of the equations. \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_LINEARITY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        LINEARITY_TYPE=EQUATIONS%LINEARITY
      ELSE
        CALL FlagError("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_LINEARITY_TYPE_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_LINEARITY_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_LINEARITY_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for equations.
  SUBROUTINE EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,LINEARITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: LINEARITY_TYPE !<The linearity type to set \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("EQUATIONS_LINEARITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FlagError("Equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(LINEARITY_TYPE)
        CASE(EQUATIONS_LINEAR)
          EQUATIONS%LINEARITY=EQUATIONS_LINEAR
        CASE(EQUATIONS_NONLINEAR)
          EQUATIONS%LINEARITY=EQUATIONS_NONLINEAR
        CASE(EQUATIONS_NONLINEAR_BCS)
          EQUATIONS%LINEARITY=EQUATIONS_NONLINEAR_BCS
        CASE DEFAULT
          LOCAL_ERROR="The specified linearity type of "//TRIM(NUMBER_TO_VSTRING(LINEARITY_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_LINEARITY_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_LINEARITY_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_LINEARITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the lumping type for equations.
  SUBROUTINE EQUATIONS_LUMPING_TYPE_GET(EQUATIONS,LUMPING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to get the lumping type for
    INTEGER(INTG), INTENT(OUT) :: LUMPING_TYPE !<On exit, the lumping type of the equations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_LUMPING_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        LUMPING_TYPE=EQUATIONS%LUMPING_TYPE
      ELSE
        CALL FlagError("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_LUMPING_TYPE_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_LUMPING_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_LUMPING_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the matrix lumping for the equations.
  SUBROUTINE EQUATIONS_LUMPING_TYPE_SET(EQUATIONS,LUMPING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to set the lumping for
    INTEGER(INTG), INTENT(IN) :: LUMPING_TYPE !<The lumping type to set \see EQUATIONS_ROUTINES_LumpingTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("EQUATIONS_LUMPING_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FlagError("Equations has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(EQUATIONS%TIME_DEPENDENCE==EQUATIONS_FIRST_ORDER_DYNAMIC.OR. &
          & EQUATIONS%TIME_DEPENDENCE==EQUATIONS_SECOND_ORDER_DYNAMIC) THEN
          SELECT CASE(LUMPING_TYPE)
          CASE(EQUATIONS_UNLUMPED_MATRICES)
            EQUATIONS%LUMPING_TYPE=EQUATIONS_UNLUMPED_MATRICES
          CASE(EQUATIONS_LUMPED_MATRICES)
            EQUATIONS%LUMPING_TYPE=EQUATIONS_LUMPED_MATRICES
          CASE DEFAULT
            LOCAL_ERROR="The specified lumping type of "//TRIM(NUMBER_TO_VSTRING(LUMPING_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          LOCAL_ERROR="Invalid equations time dependence. The equations time dependence of "// &
            & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))// &
            & " does not correspond to dynamic equations. You can only set lumping for dynamic equations."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_LUMPING_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_LUMPING_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_LUMPING_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the output type for equations.
  SUBROUTINE EQUATIONS_OUTPUT_TYPE_GET(EQUATIONS,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: OUTPUT_TYPE !<On exit, the output type of the equations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_OUTPUT_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        OUTPUT_TYPE=EQUATIONS%OUTPUT_TYPE
      ELSE
        CALL FlagError("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_OUTPUT_TYPE_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_OUTPUT_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_OUTPUT_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the equations.
  SUBROUTINE EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to set the output type for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The output type to set \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FlagError("Equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(OUTPUT_TYPE)
        CASE(EQUATIONS_NO_OUTPUT)
          EQUATIONS%OUTPUT_TYPE=EQUATIONS_NO_OUTPUT
        CASE(EQUATIONS_TIMING_OUTPUT)
          EQUATIONS%OUTPUT_TYPE=EQUATIONS_TIMING_OUTPUT
        CASE(EQUATIONS_MATRIX_OUTPUT)
          EQUATIONS%OUTPUT_TYPE=EQUATIONS_MATRIX_OUTPUT
        CASE(EQUATIONS_ELEMENT_MATRIX_OUTPUT)
          EQUATIONS%OUTPUT_TYPE=EQUATIONS_ELEMENT_MATRIX_OUTPUT
        CASE(EQUATIONS_NODAL_MATRIX_OUTPUT)
          EQUATIONS%OUTPUT_TYPE=EQUATIONS_NODAL_MATRIX_OUTPUT
        CASE DEFAULT
          LOCAL_ERROR="The specified output type of "//TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_OUTPUT_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_OUTPUT_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Gets the sparsity type for equations.
  SUBROUTINE EQUATIONS_SPARSITY_TYPE_GET(EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: SPARSITY_TYPE !<On exit, the sparsity type of the equations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_SPARSITY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        SPARSITY_TYPE=EQUATIONS%SPARSITY_TYPE
      ELSE
        CALL FlagError("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_SPARSITY_TYPE_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SPARSITY_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_SPARSITY_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the equations.
  SUBROUTINE EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The sparsity type to set \see EQUATIONS_ROUTINES_SparsityTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FlagError("Equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(SPARSITY_TYPE)
        CASE(EQUATIONS_SPARSE_MATRICES)
          EQUATIONS%SPARSITY_TYPE=EQUATIONS_SPARSE_MATRICES
        CASE(EQUATIONS_FULL_MATRICES)
          EQUATIONS%SPARSITY_TYPE=EQUATIONS_FULL_MATRICES
        CASE DEFAULT
          LOCAL_ERROR="The specified sparsity type of "//TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_SPARSITY_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_SPARSITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for equations.
  SUBROUTINE EQUATIONS_TIME_DEPENDENCE_TYPE_GET(EQUATIONS,TIME_DEPENDENCE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: TIME_DEPENDENCE_TYPE !<On exit, the time dependence type of the equations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_TIME_DEPENDENCE_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        TIME_DEPENDENCE_TYPE=EQUATIONS%TIME_DEPENDENCE
      ELSE
        CALL FlagError("Equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_TIME_DEPENDENCE_TYPE_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_TIME_DEPENDENCE_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_TIME_DEPENDENCE_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for equations.
  SUBROUTINE EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,TIME_DEPENDENCE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: TIME_DEPENDENCE_TYPE !<The time dependence type to set \see EQUATIONS_ROUTINES_TimeDependenceTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("EQUATIONS_TIME_DEPENDENCE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(EQUATIONS%EQUATIONS_FINISHED) THEN
        CALL FlagError("Equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TIME_DEPENDENCE_TYPE)
        CASE(EQUATIONS_STATIC)
          EQUATIONS%TIME_DEPENDENCE=EQUATIONS_STATIC
        CASE(EQUATIONS_QUASISTATIC)
          EQUATIONS%TIME_DEPENDENCE=EQUATIONS_QUASISTATIC
        CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
          EQUATIONS%TIME_DEPENDENCE=EQUATIONS_FIRST_ORDER_DYNAMIC
        CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
          EQUATIONS%TIME_DEPENDENCE=EQUATIONS_SECOND_ORDER_DYNAMIC
        CASE DEFAULT
          LOCAL_ERROR="The specified time dependence type of "//TRIM(NUMBER_TO_VSTRING(TIME_DEPENDENCE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_TIME_DEPENDENCE_TYPE_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_TIME_DEPENDENCE_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE EQUATIONS_TIME_DEPENDENCE_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Gets the equations for an equations set.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to get the equations for
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<On exit, a pointer to the equations in the specified equations set. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_SET_EQUATIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        IF(ASSOCIATED(EQUATIONS)) THEN
          CALL FlagError("Equations is already associated.",ERR,ERROR,*999)
        ELSE
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
          IF(.NOT.ASSOCIATED(EQUATIONS)) CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_EQUATIONS_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_EQUATIONS_GET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_GET

  !
  !================================================================================================================================
  !

  !>Gets the field variable for the derived variable type
  SUBROUTINE Equations_DerivedVariableGet(equations,derivedType,fieldVariable,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<A pointer to the equations to get the field variable for.
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to get the field variable for. \see EQUATIONS_SET_CONSTANTS_DerivedTypes.
    TYPE(FIELD_VARIABLE_TYPE), POINTER, INTENT(INOUT) :: fieldVariable !<On return, the field variable for the derived variable type.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: derivedField
    INTEGER(INTG) :: fieldVariableType

    ENTERS("Equations_DerivedVariableGet",err,error,*999)

    NULLIFY(derivedField)

    !Check pointers
    IF(ASSOCIATED(equations)) THEN
      equationsSet=>equations%equations_set
      IF(ASSOCIATED(equationsSet)) THEN
        IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) THEN
          CALL FlagError("Equations set has not been finished.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Equations set is not associated.",err,error,*999)
      END IF
      IF(ASSOCIATED(fieldVariable)) THEN
        CALL FlagError("Derived field variable is already associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations are not associated.",err,error,*999)
    END IF

    IF(ASSOCIATED(equationsSet%derived)) THEN
      IF(equationsSet%derived%derivedFinished) THEN
        IF(ASSOCIATED(equationsSet%derived%derivedField)) THEN
          IF(derivedType>0.AND.derivedType<=EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES) THEN
            fieldVariableType=equationsSet%derived%variableTypes(derivedType)
            IF(fieldVariableType/=0) THEN
              CALL FIELD_VARIABLE_GET(equationsSet%derived%derivedField,fieldVariableType,fieldVariable,err,error,*999)
            ELSE
              CALL FlagError("The field variable type for the derived variable type of "// &
                & TRIM(NUMBER_TO_VSTRING(derivedType,"*",err,error))//" has not been set.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("The derived variable type of "//TRIM(NUMBER_TO_VSTRING(derivedType,"*",err,error))// &
              & " is invalid. It should be between 1 and "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES,"*", &
              & err,error))//" inclusive.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set derived field is not associated",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Equations set derived information is not finished",err,error,*999)
      END IF
    END IF

    EXITS("Equations_DerivedVariableGet")
    RETURN
999 ERRORSEXITS("Equations_DerivedVariableGet",err,error)
    RETURN 1
  END SUBROUTINE Equations_DerivedVariableGet

  !
  !================================================================================================================================
  !

  !>Get the number of linear matrices in the equations
  SUBROUTINE Equations_NumberOfLinearMatricesGet(equations,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the number of linear matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of linear matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices

    ENTERS("Equations_NumberOfLinearMatricesGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        linearMatrices=>equationsMatrices%linear_matrices
        IF(ASSOCIATED(linearMatrices)) THEN
          numberOfMatrices=linearMatrices%number_of_linear_matrices
        ELSE
          numberOfMatrices=0
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_NumberOfLinearMatricesGet")
    RETURN
999 ERRORSEXITS("Equations_NumberOfLinearMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_NumberOfLinearMatricesGet

  !
  !================================================================================================================================
  !

  !>Get the number of Jacobian matrices in the equations
  SUBROUTINE Equations_NumberOfJacobianMatricesGet(equations,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the number of Jacobian matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of Jacobian matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices

    ENTERS("Equations_NumberOfJacobianMatricesGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        nonlinearMatrices=>equationsMatrices%nonlinear_matrices
        IF(ASSOCIATED(nonlinearMatrices)) THEN
          numberOfMatrices=nonlinearMatrices%number_of_jacobians
        ELSE
          numberOfMatrices=0
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_NumberOfJacobianMatricesGet")
    RETURN
999 ERRORSEXITS("Equations_NumberOfJacobianMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_NumberOfJacobianMatricesGet

  !
  !================================================================================================================================
  !

  !>Get the number of dynamic matrices in the equations
  SUBROUTINE Equations_NumberOfDynamicMatricesGet(equations,numberOfMatrices,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the number of dynamic matrices for
    INTEGER(INTG), INTENT(OUT) :: numberOfMatrices !<On return, the number of dynamic matrices
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices

    ENTERS("Equations_NumberOfDynamicMatricesGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        dynamicMatrices=>equationsMatrices%dynamic_matrices
        IF(ASSOCIATED(dynamicMatrices)) THEN
          numberOfMatrices=dynamicMatrices%number_of_dynamic_matrices
        ELSE
          numberOfMatrices=0
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_NumberOfDynamicMatricesGet")
    RETURN
999 ERRORSEXITS("Equations_NumberOfDynamicMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_NumberOfDynamicMatricesGet

  !
  !================================================================================================================================
  !

  !>Get a linear equations matrix from equations
  SUBROUTINE Equations_LinearMatrixGet(equations,matrixIndex,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the linear matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The index of the linear matrix to get
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the linear matrix requested
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices

    ENTERS("Equations_LinearMatrixGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        linearMatrices=>equationsMatrices%linear_matrices
        IF(ASSOCIATED(linearMatrices)) THEN
          IF(matrixIndex>0.AND.matrixIndex<=linearMatrices%number_of_linear_matrices) THEN
            IF(.NOT.ASSOCIATED(matrix)) THEN
              equationsMatrix=>linearMatrices%matrices(matrixIndex)%ptr
              IF(ASSOCIATED(equationsMatrix)) THEN
                matrix=>equationsMatrix%matrix
              ELSE
                CALL FlagError("The equations matrix is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("The matrix is already associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Invalid matrix index. The matrix index must be greater than zero and less than or equal to "// &
              & TRIM(NumberToVstring(linearMatrices%number_of_linear_matrices,"*",err,error))//".",err,error,*999)
          END IF
        ELSE
          CALL FlagError("The equations linear matrices are not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_LinearMatrixGet")
    RETURN
999 ERRORSEXITS("Equations_LinearMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_LinearMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a Jacobian matrix from equations
  SUBROUTINE Equations_JacobianMatrixGet(equations,residualIndex,variableType,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the Jacobian matrix for
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type that the residual is differentiated with respect to for this Jacobian
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the requested Jacobian matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: matrixIndex,variableIndex
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: equationsJacobian

    ENTERS("Equations_JacobianMatrixGet",err,error,*999)

    !Check for pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        nonlinearMapping=>equationsMapping%nonlinear_mapping
        IF(.NOT.ASSOCIATED(nonlinearMapping)) THEN
          CALL FlagError("The equations nonlinear mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF
    IF(ASSOCIATED(matrix)) THEN
      CALL FlagError("The matrix is already associated.",err,error,*999)
    END IF

    IF(residualIndex/=1) THEN
      CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    END IF

    !Find Jacobian matrix index using the nonlinear equations mapping
    matrixIndex=0
    DO variableIndex=1,nonlinearMapping%number_of_residual_variables
      IF(nonlinearMapping%residual_variables(variableIndex)%ptr%variable_type==variableType) THEN
        matrixIndex=nonlinearMapping%var_to_jacobian_map(variableIndex)%jacobian_number
      END IF
    END DO
    IF(matrixIndex==0) THEN
      CALL FlagError("Equations do not have a Jacobian matrix for residual index "// &
        & TRIM(NumberToVstring(residualIndex,"*",err,error))//" and variable type "// &
        & TRIM(NumberToVstring(variableType,"*",err,error))//".",err,error,*999)
    END IF

    !Now get Jacobian matrix using the matrix index
    equationsJacobian=>nonlinearMapping%jacobian_to_var_map(matrixIndex)%jacobian
    IF(ASSOCIATED(equationsJacobian)) THEN
      matrix=>equationsJacobian%jacobian
    ELSE
      CALL FlagError("The equations Jacobian matrix is not associated.",err,error,*999)
    END IF

    EXITS("Equations_JacobianMatrixGet")
    RETURN
999 ERRORSEXITS("Equations_JacobianMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_JacobianMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a dynamic equations matrix from equations using the dynamic matrix index
  SUBROUTINE Equations_DynamicMatrixGet(equations,matrixIndex,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The number of the dynamic matrix to get
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the requested dynamic matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices

    ENTERS("Equations_DynamicMatrixGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        dynamicMatrices=>equationsMatrices%dynamic_matrices
        IF(ASSOCIATED(dynamicMatrices)) THEN
          IF(matrixIndex>0.AND.matrixIndex<=dynamicMatrices%number_of_dynamic_matrices) THEN
            IF(.NOT.ASSOCIATED(matrix)) THEN
              equationsMatrix=>dynamicMatrices%matrices(matrixIndex)%ptr
              IF(ASSOCIATED(equationsMatrix)) THEN
                matrix=>equationsMatrix%matrix
              ELSE
                CALL FlagError("The equations matrix is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("The matrix is already associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Invalid matrix index. The matrix index must be greater than zero and less than or equal to "// &
              & TRIM(NumberToVstring(dynamicMatrices%number_of_dynamic_matrices,"*",err,error))//".",err,error,*999)
          END IF
        ELSE
          CALL FlagError("The equations dynamic matrices are not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_DynamicMatrixGet")
    RETURN
999 ERRORSEXITS("Equations_DynamicMatrixGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixGet

  !
  !================================================================================================================================
  !

  !>Get a dynamic equations matrix from equations using the dynamic matrix type
  SUBROUTINE Equations_DynamicMatrixGetByType(equations,matrixType,matrix,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixType !<The type of the dynamic matrix to get. \see EQUATIONS_SET_CONSTANTS_DynamicMatrixTypes,EQUATIONS_SET_CONSTANTS
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER, INTENT(INOUT) :: matrix !<On return, the requested dynamic matrix
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: matrixIndex
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: equationsMatrix
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: dynamicMatrices
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping

    ENTERS("Equations_DynamicMatrixGetByType",err,error,*999)

    !Check all pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        dynamicMatrices=>equationsMatrices%dynamic_matrices
        IF(.NOT.ASSOCIATED(dynamicMatrices)) THEN
          CALL FlagError("The equations dynamic matrices are not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
        IF(.NOT.ASSOCIATED(dynamicMapping)) THEN
          CALL FlagError("The equations dynamic mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF
    IF(ASSOCIATED(matrix)) THEN
      CALL FlagError("The matrix is already associated.",err,error,*999)
    END IF

    !Now get the dynamic matrix
    !Find matrix index using the equations mapping
    SELECT CASE(matrixType)
    CASE(EQUATIONS_MATRIX_STIFFNESS)
      matrixIndex=dynamicMapping%stiffness_matrix_number
    CASE(EQUATIONS_MATRIX_DAMPING)
      matrixIndex=dynamicMapping%damping_matrix_number
    CASE(EQUATIONS_MATRIX_MASS)
      matrixIndex=dynamicMapping%mass_matrix_number
    CASE DEFAULT
      CALL FlagError("Invalid dynamic matrix type "//TRIM(NumberToVstring(matrixType,"*",err,error))// &
        & " specified.",err,error,*999)
    END SELECT
    IF(matrixIndex==0) THEN
      CALL FlagError("The equations dynamic matrices do not have a matrix with the specified type of "// &
        & TRIM(NumberToVstring(matrixType,"*",err,error))//".",err,error,*999)
    ELSE
      equationsMatrix=>dynamicMatrices%matrices(matrixIndex)%ptr
      IF(ASSOCIATED(equationsMatrix)) THEN
        matrix=>equationsMatrix%matrix
      ELSE
        CALL FlagError("The equations dynamic matrix for index "// &
          & TRIM(NumberToVstring(matrixIndex,"*",err,error))//" is not associated.",err,error,*999)
      END IF
    END IF

    EXITS("Equations_DynamicMatrixGetByType")
    RETURN
999 ERRORSEXITS("Equations_DynamicMatrixGetByType",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixGetByType

  !
  !================================================================================================================================
  !

  !>Get the type of a dynamic matrix, eg. stiffness, damping or mass
  SUBROUTINE Equations_DynamicMatrixTypeGet(equations,matrixIndex,matrixType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the dynamic matrix for
    INTEGER(INTG), INTENT(IN) :: matrixIndex !<The number of the dynamic matrix to get
    INTEGER(INTG), INTENT(INOUT) :: matrixType !<On return, the type of the dynamic matrix. \see EQUATIONS_MATRICES_ROUTINES_DynamicMatrixTypes,EQUATIONS_MATRICES_ROUTINES
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: dynamicMapping

    ENTERS("Equations_DynamicMatrixTypeGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        dynamicMapping=>equationsMapping%DYNAMIC_MAPPING
        IF(ASSOCIATED(dynamicMapping)) THEN
          IF(matrixIndex>0.AND.matrixIndex<=dynamicMapping%number_of_dynamic_equations_matrices) THEN
            IF(matrixIndex==dynamicMapping%stiffness_matrix_number) THEN
              matrixType=EQUATIONS_MATRIX_STIFFNESS
            ELSE IF(matrixIndex==dynamicMapping%damping_matrix_number) THEN
              matrixType=EQUATIONS_MATRIX_DAMPING
            ELSE IF(matrixIndex==dynamicMapping%mass_matrix_number) THEN
              matrixType=EQUATIONS_MATRIX_MASS
            ELSE
              CALL FlagError("Could not find dynamic matrix type.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Invalid matrix index. The matrix index must be greater than zero and less than or equal to "// &
              & TRIM(NumberToVstring(dynamicMapping%number_of_dynamic_equations_matrices,"*",err,error))//".",err,error,*999)
          END IF
        ELSE
          CALL FlagError("The equations dynamic mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_DynamicMatrixTypeGet")
    RETURN
999 ERRORSEXITS("Equations_DynamicMatrixTypeGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_DynamicMatrixTypeGet

  !
  !================================================================================================================================
  !

  !>Get the right hand side vector for equations
  SUBROUTINE Equations_RhsVectorGet(equations,vector,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the right hand side vector for
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER, INTENT(INOUT) :: vector !<On return, the right hand side vector for the equations
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices

    ENTERS("Equations_RhsVectorGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        rhsVector=>equationsMatrices%rhs_vector
        IF(ASSOCIATED(rhsVector)) THEN
          IF(.NOT.ASSOCIATED(vector)) THEN
            vector=>rhsVector%vector
          ELSE
            CALL FlagError("The vector is already associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("The equations matrices right hand side vector is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_RhsVectorGet")
    RETURN
999 ERRORSEXITS("Equations_RhsVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_RhsVectorGet

  !
  !================================================================================================================================
  !

  !>Get a residual vector for nonlinear equations
  SUBROUTINE Equations_ResidualVectorGet(equations,residualIndex,vector,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER, INTENT(INOUT) :: vector !<On return, the residual vector for the equations
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices

    ENTERS("Equations_ResidualVectorGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        nonlinearMatrices=>equationsMatrices%nonlinear_matrices
        IF(ASSOCIATED(nonlinearMatrices)) THEN
          IF(.NOT.ASSOCIATED(vector)) THEN
            IF(residualIndex==1) THEN
              vector=>nonlinearMatrices%residual
            ELSE
              CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("The vector is already associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("The equations matrices nonlinear matrices are not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_ResidualVectorGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualVectorGet

  !
  !================================================================================================================================
  !

  !>Get the number of field variables that contribute to the residual vector
  SUBROUTINE Equations_ResidualNumberOfVariablesGet(equations,residualIndex,numberOfVariables,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector number of variables for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: numberOfVariables !<On return, the number of variables that contribute to the residual vector
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping

    ENTERS("Equations_ResidualNumberOfVariablesGet",err,error,*999)

    !Check for pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        nonlinearMapping=>equationsMapping%nonlinear_mapping
        IF(.NOT.ASSOCIATED(nonlinearMapping)) THEN
          CALL FlagError("The equations nonlinear mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    IF(residualIndex==1) THEN
      numberOfVariables=nonlinearMapping%number_of_residual_variables
    ELSE
      CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    END IF

    EXITS("Equations_ResidualNumberOfVariablesGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualNumberOfVariablesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualNumberOfVariablesGet

  !
  !================================================================================================================================
  !

  !>Get the field variables that contribute to the residual vector
  SUBROUTINE Equations_ResidualVariablesGet(equations,residualIndex,residualVariables,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the residual vector variables for
    INTEGER(INTG), INTENT(IN) :: residualIndex !<The index of the residual vector to get the variables for
    INTEGER(INTG), INTENT(OUT) :: residualVariables(:) !<residualVariables(varIdx). On return, the field variable type for the varIdx'th residual variable
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    INTEGER(INTG) :: numberOfVariables,variableIdx
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping

    ENTERS("Equations_ResidualVariablesGet",err,error,*999)

    !Check for pointer associations
    IF(ASSOCIATED(equations)) THEN
      equationsMapping=>equations%equations_mapping
      IF(ASSOCIATED(equationsMapping)) THEN
        nonlinearMapping=>equationsMapping%nonlinear_mapping
        IF(.NOT.ASSOCIATED(nonlinearMapping)) THEN
          CALL FlagError("The equations nonlinear mapping is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations mapping is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    IF(residualIndex==1) THEN
      numberOfVariables=nonlinearMapping%number_of_residual_variables
      IF(SIZE(residualVariables,1)>=numberOfVariables) THEN
        DO variableIdx=1,numberOfVariables
          residualVariables(variableIdx)=nonlinearMapping%residual_variables(variableIdx)%ptr%variable_type
        END DO
      ELSE
        CALL FlagError("residualVariables array must have size of at least "// &
          & TRIM(numberToVstring(numberOfVariables,"*",err,error))//".",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Multiple residual vectors are not yet implemented.",err,error,*999)
    END IF

    EXITS("Equations_ResidualVariablesGet")
    RETURN
999 ERRORSEXITS("Equations_ResidualVariablesGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_ResidualVariablesGet

  !
  !================================================================================================================================
  !

  !>Get the source vector for equations
  SUBROUTINE Equations_SourceVectorGet(equations,vector,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER, INTENT(IN) :: equations !<The equations to get the source vector for
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER, INTENT(INOUT) :: vector !<On return, the source vector for the equations
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error message
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    !Local variables
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: matricesSource
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices

    ENTERS("Equations_SourceVectorGet",err,error,*999)

    IF(ASSOCIATED(equations)) THEN
      equationsMatrices=>equations%equations_matrices
      IF(ASSOCIATED(equationsMatrices)) THEN
        matricesSource=>equationsMatrices%source_vector
        IF(ASSOCIATED(matricesSource)) THEN
          IF(.NOT.ASSOCIATED(vector)) THEN
            vector=>matricesSource%vector
          ELSE
            CALL FlagError("The vector is already associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("The equations matrices source vector is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("The equations matrices are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("The equations are not associated.",err,error,*999)
    END IF

    EXITS("Equations_SourceVectorGet")
    RETURN
999 ERRORSEXITS("Equations_SourceVectorGet",err,error)
    RETURN 1

  END SUBROUTINE Equations_SourceVectorGet

  !
  !================================================================================================================================
  !

END MODULE EQUATIONS_ROUTINES
