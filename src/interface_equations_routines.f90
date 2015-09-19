!> \file
!> \author Chris Bradley
!> \brief This module handles all interface equations routines.
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

!>This module handles all interface equations routines.
MODULE INTERFACE_EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE EQUATIONS_ROUTINES
  USE FIELD_ROUTINES
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_MAPPING_ROUTINES
  USE INTERFACE_MATRICES_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup INTERFACE_EQUATIONS_ROUTINES_OutputTypes INTERFACE_EQUATIONS_ROUTINES::OutputTypes
  !> \brief The interface equations output types
  !> \see INTERFACE_EQUATIONS_ROUTINES,OPENCMISS_InterfaceEquationsConstants
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_NO_OUTPUT=0 !<No output. \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_TIMING_OUTPUT=1 !<Timing information output. \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output. \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output .\see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES 
  !>@}

  !> \addtogroup INTERFACE_EQUATIONS_ROUTINES_SparsityTypes INTERFACE_EQUATIONS_ROUTINES::SparsityTypes
  !> \brief Interface equations matrices sparsity types
  !> \see INTERFACE_EQUATIONS_ROUTINES,OPENCMISS_InterfaceEquationsSparsityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the interface equations. \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the interface equations. \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
 !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC INTERFACE_EQUATIONS_NO_OUTPUT,INTERFACE_EQUATIONS_TIMING_OUTPUT,INTERFACE_EQUATIONS_MATRIX_OUTPUT, &
    & INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT

  PUBLIC INTERFACE_EQUATIONS_SPARSE_MATRICES,INTERFACE_EQUATIONS_FULL_MATRICES

  PUBLIC INTERFACE_EQUATIONS_CREATE_FINISH,INTERFACE_EQUATIONS_CREATE_START

  PUBLIC INTERFACE_EQUATIONS_DESTROY

  PUBLIC InterfaceEquations_InterfaceInterpSetsNumberSet

  PUBLIC INTERFACE_EQUATIONS_OUTPUT_TYPE_GET,INTERFACE_EQUATIONS_OUTPUT_TYPE_SET
  
  PUBLIC INTERFACE_EQUATIONS_SPARSITY_TYPE_GET,INTERFACE_EQUATIONS_SPARSITY_TYPE_SET

  PUBLIC InterfaceEquations_VariableInterpSetsNumberSet

  PUBLIC INTERFACE_EQUATIONS_LINEARITY_TYPE_GET,INTERFACE_EQUATIONS_LINEARITY_TYPE_SET

  PUBLIC InterfaceEquations_TimeDependenceTypeGet,InterfaceEquationsTimeDependenceTypeSet

  PUBLIC INTERFACE_CONDITION_EQUATIONS_GET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_CREATE_FINISH(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,LAGRANGE_FIELD,PENALTY_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("INTERFACE_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations have already been finished.",ERR,ERROR,*999)
      ELSE
        !Create the interpolation sets
        INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
        IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERPOLATION)) THEN
                  GEOMETRIC_FIELD=>INTERFACE_CONDITION%GEOMETRY%GEOMETRIC_FIELD
                  LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
                  NULLIFY(PENALTY_FIELD)
                  IF(ASSOCIATED(INTERFACE_CONDITION%PENALTY)) THEN
                    PENALTY_FIELD=>INTERFACE_CONDITION%PENALTY%PENALTY_FIELD
                  ENDIF
                  !\todo Truncating subroutine name from INTERFACE_EQUATIONS_DOMAIN_INTERFACE_INTERPOLATION_SETUP until bug in gfortran 4.6 is fixed http://gcc.gnu.org/bugzilla/show_bug.cgi?id=46971
                  CALL INTERFACE_EQUATIONS_DOMAIN_INTERFACE_INTERPOLATION_(INTERFACE_EQUATIONS%INTERPOLATION% &
                    & INTERFACE_INTERPOLATION,GEOMETRIC_FIELD,LAGRANGE_FIELD,PENALTY_FIELD,ERR,ERROR,*999)
                  DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
                    DEPENDENT_VARIABLE=>INTERFACE_DEPENDENT%FIELD_VARIABLES(variable_idx)%PTR
                    IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                      DEPENDENT_FIELD=>DEPENDENT_VARIABLE%FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        GEOMETRIC_FIELD=>DEPENDENT_FIELD%GEOMETRIC_FIELD
                        CALL InterfaceEquations_DomainVariableInterpolationSetup(INTERFACE_EQUATIONS%INTERPOLATION% &
                          & VARIABLE_INTERPOLATION(variable_idx),GEOMETRIC_FIELD,DEPENDENT_FIELD,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("Dependent variable field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Dependent variable is not associated for variable index "// &
                        & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !variable_idx
                ELSE
                  CALL FlagError("Interface equations interpolation is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition dependent is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Set the finished flag
          INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED=.TRUE.
        ELSE
          CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of interface equations for an interface condition.
  SUBROUTINE INTERFACE_EQUATIONS_CREATE_START(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<On exit, a pointer to the created interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) THEN
        CALL FlagError("Interface equations are already associated for the interface condition.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          CALL FlagError("Interface equations is already associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the equations
          CALL INTERFACE_EQUATIONS_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_CREATE_START")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_CREATE_START",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_DESTROY(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      CALL INTERFACE_EQUATIONS_FINALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_DESTROY")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_DESTROY",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the interface equations domain interpolation and deallocates all memory.
  SUBROUTINE InterfaceEquations_DomainInterpolationFinalise(DOMAIN_INTERPOLATION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE) :: DOMAIN_INTERPOLATION !<The domain interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interpolation_set_idx
 
    ENTERS("InterfaceEquations_DomainInterpolationFinalise",ERR,ERROR,*999)

    NULLIFY(DOMAIN_INTERPOLATION%GEOMETRIC_FIELD)
    IF(ALLOCATED(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION)) THEN
      DO interpolation_set_idx=1,SIZE(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION,1)
        CALL InterfaceEquations_InterpolationSetFinalise(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx), &
          & ERR,ERROR,*999)
      ENDDO !interpolation_set_idx
      DEALLOCATE(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION)
    ENDIF
    DOMAIN_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS=0
    NULLIFY(DOMAIN_INTERPOLATION%DEPENDENT_FIELD)
    IF(ALLOCATED(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION)) THEN
     DO interpolation_set_idx=1,SIZE(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION,1)
        CALL InterfaceEquations_InterpolationSetFinalise(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(interpolation_set_idx), &
          & ERR,ERROR,*999)
      ENDDO !interpolation_set_idx
      DEALLOCATE(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION)
    ENDIF
    DOMAIN_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS=0
       
    EXITS("InterfaceEquations_DomainInterpolationFinalise")
    RETURN
999 ERRORS("InterfaceEquations_DomainInterpolationFinalise",ERR,ERROR)
    EXITS("InterfaceEquations_DomainInterpolationFinalise")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_DomainInterpolationFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations domain interpolation.
  SUBROUTINE InterfaceEquations_DomainInterpolationInitialise(DOMAIN_INTERPOLATION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE) :: DOMAIN_INTERPOLATION !<The domain interpolation to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_DomainInterpolationInitialise",ERR,ERROR,*999)

    NULLIFY(DOMAIN_INTERPOLATION%PENALTY_FIELD)
    DOMAIN_INTERPOLATION%NUMBER_OF_PENALTY_INTERPOLATION_SETS=1
    NULLIFY(DOMAIN_INTERPOLATION%GEOMETRIC_FIELD)
    DOMAIN_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS=1
    NULLIFY(DOMAIN_INTERPOLATION%DEPENDENT_FIELD)
    DOMAIN_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS=1
       
    EXITS("InterfaceEquations_DomainInterpolationInitialise")
    RETURN
999 ERRORS("InterfaceEquations_DomainInterpolationInitialise",ERR,ERROR)
    EXITS("InterfaceEquations_DomainInterpolationInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_DomainInterpolationInitialise

  !
  !================================================================================================================================
  !

  !>Sets up the interface equations domain interface interpolation. !\todo Truncating subroutine name from INTERFACE_EQUATIONS_DOMAIN_INTERFACE_INTERPOLATION_SETUP until bug in gfortran 4.6 is fixed http://gcc.gnu.org/bugzilla/show_bug.cgi?id=46971
  SUBROUTINE INTERFACE_EQUATIONS_DOMAIN_INTERFACE_INTERPOLATION_(DOMAIN_INTERPOLATION,GEOMETRIC_FIELD,LAGRANGE_FIELD, &
    & PENALTY_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE) :: DOMAIN_INTERPOLATION !<The domain interpolation to initialise
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field to set up the domain interpolation for
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD !<A pointer to the Lagrange field to set up the domain interpoaltion for
    TYPE(FIELD_TYPE), POINTER :: PENALTY_FIELD !<A pointer to the penalty field to set up the domain interpoaltion for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,interpolation_set_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("InterfaceEquations_DomainInterpolationSet",ERR,ERROR,*998)

    IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
      IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
        DOMAIN_INTERPOLATION%GEOMETRIC_FIELD=>GEOMETRIC_FIELD
        ALLOCATE(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(DOMAIN_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS), &
          & STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain interpolation geometric interpolation.",ERR,ERROR,*999)
        DO interpolation_set_idx=1,DOMAIN_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS
          CALL InterfaceEquations_InterpolationSetInitialise(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION( &
            & interpolation_set_idx),ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DOMAIN_INTERPOLATION%GEOMETRIC_FIELD,DOMAIN_INTERPOLATION% &
            & GEOMETRIC_INTERPOLATION(interpolation_set_idx)%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINTS_INITIALISE(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx)% &
            & INTERPOLATION_PARAMETERS,DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx)%INTERPOLATED_POINT, &
            & ERR,ERROR,*999)
          IF(DOMAIN_INTERPOLATION%GEOMETRIC_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
            & DOMAIN_INTERPOLATION%GEOMETRIC_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
            CALL Field_InterpolatedPointsMetricsInitialise(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION( &
              & interpolation_set_idx)%INTERPOLATED_POINT,DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx)% &
              & INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
        ENDDO !interpolation_set_idx
        DOMAIN_INTERPOLATION%DEPENDENT_FIELD=>LAGRANGE_FIELD
        ALLOCATE(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(DOMAIN_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS), &
          & STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain interpolation dependent interpolation.",ERR,ERROR,*999)
        DO interpolation_set_idx=1,DOMAIN_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS
          CALL InterfaceEquations_InterpolationSetInitialise(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION( &
            & interpolation_set_idx),ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DOMAIN_INTERPOLATION%DEPENDENT_FIELD,DOMAIN_INTERPOLATION% &
            & DEPENDENT_INTERPOLATION(interpolation_set_idx)%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINTS_INITIALISE(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(interpolation_set_idx)% &
            & INTERPOLATION_PARAMETERS,DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(interpolation_set_idx)%INTERPOLATED_POINT, &
            & ERR,ERROR,*999)
        ENDDO !interpolation_set_idx
        IF(ASSOCIATED(PENALTY_FIELD)) THEN
          DOMAIN_INTERPOLATION%PENALTY_FIELD=>PENALTY_FIELD
          ALLOCATE(DOMAIN_INTERPOLATION%PENALTY_INTERPOLATION(DOMAIN_INTERPOLATION%NUMBER_OF_PENALTY_INTERPOLATION_SETS), &
            & STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate domain interpolation dependent interpolation.",ERR,ERROR,*999)
          DO interpolation_set_idx=1,DOMAIN_INTERPOLATION%NUMBER_OF_PENALTY_INTERPOLATION_SETS
            CALL InterfaceEquations_InterpolationSetInitialise(DOMAIN_INTERPOLATION%PENALTY_INTERPOLATION( &
              & interpolation_set_idx),ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DOMAIN_INTERPOLATION%PENALTY_FIELD,DOMAIN_INTERPOLATION% &
              & PENALTY_INTERPOLATION(interpolation_set_idx)%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(DOMAIN_INTERPOLATION%PENALTY_INTERPOLATION(interpolation_set_idx)% &
              & INTERPOLATION_PARAMETERS,DOMAIN_INTERPOLATION%PENALTY_INTERPOLATION(interpolation_set_idx)%INTERPOLATED_POINT, &
              & ERR,ERROR,*999)
          ENDDO !interpolation_set_idx
        ENDIF
      ELSE
        CALL FlagError("Lagrange field is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("InterfaceEquations_DomainInterpolationSet")
    RETURN
999 CALL InterfaceEquations_DomainInterpolationFinalise(DOMAIN_INTERPOLATION,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("InterfaceEquations_DomainInterpolationSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_DOMAIN_INTERFACE_INTERPOLATION_

  !
  !================================================================================================================================
  !

  !>Sets up the interface equations domain variable interpolation. \todo Truncating subroutine name until bug in gfortran 4.6 is fixed http://gcc.gnu.org/bugzilla/show_bug.cgi?id=46971
  SUBROUTINE InterfaceEquations_DomainVariableInterpolationSetup(DOMAIN_INTERPOLATION,GEOMETRIC_FIELD,DEPENDENT_FIELD, &
    &  ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE) :: DOMAIN_INTERPOLATION !<The domain interpolation to initialise
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field to set up the domain interpolation for
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<A pointer to the depdendent field to set up the domain interpoaltion for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,interpolation_set_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("InterfaceEquations_DomainVariableInterpolationSetup",ERR,ERROR,*998)

    IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        DOMAIN_INTERPOLATION%GEOMETRIC_FIELD=>GEOMETRIC_FIELD
        ALLOCATE(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(DOMAIN_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS), &
          & STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain interpolation geometric interpolation.",ERR,ERROR,*999)
        DO interpolation_set_idx=1,DOMAIN_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS
          CALL InterfaceEquations_InterpolationSetInitialise(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION( &
            & interpolation_set_idx),ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DOMAIN_INTERPOLATION%GEOMETRIC_FIELD,DOMAIN_INTERPOLATION% &
            & GEOMETRIC_INTERPOLATION(interpolation_set_idx)%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINTS_INITIALISE(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx)% &
            & INTERPOLATION_PARAMETERS,DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx)%INTERPOLATED_POINT, &
            & ERR,ERROR,*999)
          IF(DOMAIN_INTERPOLATION%GEOMETRIC_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
            & DOMAIN_INTERPOLATION%GEOMETRIC_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
            CALL Field_InterpolatedPointsMetricsInitialise(DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION( &
              & interpolation_set_idx)%INTERPOLATED_POINT,DOMAIN_INTERPOLATION%GEOMETRIC_INTERPOLATION(interpolation_set_idx)% &
              & INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
        ENDDO !interpolation_set_idx
        DOMAIN_INTERPOLATION%DEPENDENT_FIELD=>DEPENDENT_FIELD
        ALLOCATE(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(DOMAIN_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS), &
          & STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate domain interpolation dependent interpolation.",ERR,ERROR,*999)
        DO interpolation_set_idx=1,DOMAIN_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS
          CALL InterfaceEquations_InterpolationSetInitialise(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION( &
            & interpolation_set_idx),ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DOMAIN_INTERPOLATION%DEPENDENT_FIELD,DOMAIN_INTERPOLATION% &
            & DEPENDENT_INTERPOLATION(interpolation_set_idx)%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINTS_INITIALISE(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(interpolation_set_idx)% &
            & INTERPOLATION_PARAMETERS,DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(interpolation_set_idx)%INTERPOLATED_POINT, &
            & ERR,ERROR,*999)
          IF(DOMAIN_INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
            & DOMAIN_INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
            CALL Field_InterpolatedPointsMetricsInitialise(DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION( &
              & interpolation_set_idx)%INTERPOLATED_POINT,DOMAIN_INTERPOLATION%DEPENDENT_INTERPOLATION(interpolation_set_idx)% &
              & INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
        ENDDO !interpolation_set_idx
      ELSE
        CALL FlagError("Dependent field is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Geometric field is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("InterfaceEquations_DomainVariableInterpolationSetup")
    RETURN
999 CALL InterfaceEquations_DomainInterpolationFinalise(DOMAIN_INTERPOLATION,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORS("InterfaceEquations_DomainVariableInterpolationSetup",ERR,ERROR)
    EXITS("InterfaceEquations_DomainVariableInterpolationSetup")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_DomainVariableInterpolationSetup

  !
  !================================================================================================================================
  !

  SUBROUTINE InterfaceEquations_InterfaceInterpSetsNumberSet(INTERFACE_EQUATIONS,NUMBER_OF_GEOMETRIC_SETS, &
     & NUMBER_OF_DEPENDENT_SETS,NUMBER_OF_PENALTY_SETS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<The interface equations to set the interface interpolation sets for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_GEOMETRIC_SETS !<The number of geometric interface interpolation sets to set
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DEPENDENT_SETS !<The number of dependent interface interpolation sets to set
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_PENALTY_SETS !<The number of penalty interface interpolation sets to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("InterfaceEquations_InterfaceInterpSetsNumberSet",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERPOLATION)) THEN
          IF(NUMBER_OF_GEOMETRIC_SETS>0) THEN
            IF(NUMBER_OF_DEPENDENT_SETS>0) THEN
              IF(NUMBER_OF_PENALTY_SETS>=0) THEN
                INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS= &
                  & NUMBER_OF_GEOMETRIC_SETS
                INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%NUMBER_OF_DEPENDENT_INTERPOLATION_SETS= &
                  & NUMBER_OF_DEPENDENT_SETS
                INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%NUMBER_OF_PENALTY_INTERPOLATION_SETS= &
                  & NUMBER_OF_PENALTY_SETS
              ELSE
                LOCAL_ERROR="The specified number of penalty sets of "// &
                  & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_PENALTY_SETS,"*",ERR,ERROR))// &
                  & " is invalid. The number of penalty sets must be > 0."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified number of dependent sets of "// &
                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DEPENDENT_SETS,"*",ERR,ERROR))// &
                & " is invalid. The number of dependent sets must be > 0."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified number of geometric sets of "// &
              & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GEOMETRIC_SETS,"*",ERR,ERROR))// &
              & " is invalid. The number of geometric sets must be > 0."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface equations interpolation is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("InterfaceEquations_InterfaceInterpSetsNumberSet")
    RETURN
999 ERRORS("InterfaceEquations_InterfaceInterpSetsNumberSet",ERR,ERROR)
    EXITS("InterfaceEquations_InterfaceInterpSetsNumberSet")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterfaceInterpSetsNumberSet

  !
  !================================================================================================================================
  !

  !>Finalise the interface equations and deallocate all memory.
  SUBROUTINE INTERFACE_EQUATIONS_FINALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("INTERFACE_EQUATIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      CALL INTERFACE_EQUATIONS_INTERPOLATION_FINALISE(INTERFACE_EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MAPPING)) &
        & CALL INTERFACE_MAPPING_DESTROY(INTERFACE_EQUATIONS%INTERFACE_MAPPING,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MATRICES)) &
        & CALL INTERFACE_MATRICES_DESTROY(INTERFACE_EQUATIONS%INTERFACE_MATRICES,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_EQUATIONS)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations for an interface condition
  SUBROUTINE INTERFACE_EQUATIONS_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to initialise the interface equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("INTERFACE_EQUATIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) THEN
        CALL FlagError("Interface equations is already associated for this interface condition.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_CONDITION%INTERFACE_EQUATIONS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate interface equations.",ERR,ERROR,*999)
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_CONDITION=>INTERFACE_CONDITION
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%LINEARITY=INTERFACE_CONDITION_LINEAR
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%TIME_DEPENDENCE=INTERFACE_CONDITION_STATIC
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_NO_OUTPUT
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%SPARSITY_TYPE=INTERFACE_EQUATIONS_SPARSE_MATRICES
        NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERPOLATION)
        NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_MAPPING)
        NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_MATRICES)
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED=.FALSE.
        CALL InterfaceEquations_InterpolationInitialise(INTERFACE_CONDITION%INTERFACE_EQUATIONS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface condition is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_INITIALISE")
    RETURN
999 CALL INTERFACE_EQUATIONS_FINALISE(INTERFACE_CONDITION%INTERFACE_EQUATIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("INTERFACE_EQUATIONS_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the interface equations interpolation and deallocates all memory.
  SUBROUTINE INTERFACE_EQUATIONS_INTERPOLATION_FINALISE(INTERFACE_EQUATIONS_INTERPOLATION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_INTERPOLATION_TYPE), POINTER :: INTERFACE_EQUATIONS_INTERPOLATION !<A pointer to the interface equations interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
 
    ENTERS("INTERFACE_EQUATIONS_INTERPOLATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS_INTERPOLATION)) THEN
      CALL InterfaceEquations_DomainInterpolationFinalise(INTERFACE_EQUATIONS_INTERPOLATION%INTERFACE_INTERPOLATION, &
        & ERR,ERROR,*999)
      IF(ALLOCATED(INTERFACE_EQUATIONS_INTERPOLATION%VARIABLE_INTERPOLATION)) THEN
        DO variable_idx=1,SIZE(INTERFACE_EQUATIONS_INTERPOLATION%VARIABLE_INTERPOLATION,1)
          CALL InterfaceEquations_DomainInterpolationFinalise(INTERFACE_EQUATIONS_INTERPOLATION% &
            & VARIABLE_INTERPOLATION(variable_idx),ERR,ERROR,*999)
        ENDDO !variable_idx
        DEALLOCATE(INTERFACE_EQUATIONS_INTERPOLATION%VARIABLE_INTERPOLATION)
      ENDIF
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_INTERPOLATION_FINALISE")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_INTERPOLATION_FINALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_INTERPOLATION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations interpolation.
  SUBROUTINE InterfaceEquations_InterpolationInitialise(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to initialise the interpolation for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_idx
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("InterfaceEquations_InterpolationInitialise",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
      IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
        IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERPOLATION)) THEN
          CALL FlagError("Interface equations interpolation is already associated.",ERR,ERROR,*998)
        ELSE
          INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
          IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
            ALLOCATE(INTERFACE_EQUATIONS%INTERPOLATION,STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface equations interpolation.",ERR,ERROR,*999)
            INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
            CALL InterfaceEquations_DomainInterpolationInitialise(INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION, &
              & ERR,ERROR,*999)
            INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%INTERPOLATION=>INTERFACE_EQUATIONS%INTERPOLATION
            ALLOCATE(INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES), &
              & STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate interface equations interpolation mesh interpolation.",ERR,ERROR,*999)
            DO variable_idx=1,INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES
              CALL InterfaceEquations_DomainInterpolationInitialise(INTERFACE_EQUATIONS%INTERPOLATION% &
                & VARIABLE_INTERPOLATION(variable_idx),ERR,ERROR,*999)
                INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(variable_idx)%INTERPOLATION=> &
                  & INTERFACE_EQUATIONS%INTERPOLATION
            ENDDO !variable_idx
          ELSE
            CALL FlagError("Interface condition dependent is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("InterfaceEquations_InterpolationInitialise")
    RETURN
999 CALL INTERFACE_EQUATIONS_INTERPOLATION_FINALISE(INTERFACE_EQUATIONS%INTERPOLATION,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("InterfaceEquations_InterpolationInitialise",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the interface equations interpolation set and deallocates all memory.
  SUBROUTINE InterfaceEquations_InterpolationSetFinalise(INTERPOLATION_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_INTERPOLATION_SET_TYPE) :: INTERPOLATION_SET !<The interpolation set to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterpolationSetFinalise",ERR,ERROR,*999)

    CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(INTERPOLATION_SET%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    CALL FIELD_INTERPOLATED_POINTS_FINALISE(INTERPOLATION_SET%INTERPOLATED_POINT,ERR,ERROR,*999)
    CALL Field_InterpolatedPointsMetricsFinalise(INTERPOLATION_SET%INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
       
    EXITS("InterfaceEquations_InterpolationSetFinalise")
    RETURN
999 ERRORSEXITS("InterfaceEquations_InterpolationSetFinalise",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationSetFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations interpolation set.
  SUBROUTINE InterfaceEquations_InterpolationSetInitialise(INTERPOLATION_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_INTERPOLATION_SET_TYPE) :: INTERPOLATION_SET !<The interpolation set to intialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_InterpolationSetInitialise",ERR,ERROR,*999)

    NULLIFY(INTERPOLATION_SET%INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATION_SET%INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_SET%INTERPOLATED_POINT_METRICS)
       
    EXITS("InterfaceEquations_InterpolationSetInitialise")
    RETURN
999 ERRORS("InterfaceEquations_InterpolationSetInitialise",ERR,ERROR)
    EXITS("InterfaceEquations_InterpolationSetInitialise")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_InterpolationSetInitialise

  !
  !================================================================================================================================
  !

  !>Gets the output type for interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_GET(INTERFACE_EQUATIONS,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: OUTPUT_TYPE !<On exit, the output type of the interface equations \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        OUTPUT_TYPE=INTERFACE_EQUATIONS%OUTPUT_TYPE
      ELSE
        CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_SET(INTERFACE_EQUATIONS,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to set the output type for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The output type to set \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(OUTPUT_TYPE)
        CASE(INTERFACE_EQUATIONS_NO_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_NO_OUTPUT
        CASE(INTERFACE_EQUATIONS_TIMING_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_TIMING_OUTPUT
        CASE(INTERFACE_EQUATIONS_MATRIX_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_MATRIX_OUTPUT
        CASE(INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
        CASE DEFAULT
          LOCAL_ERROR="The specified output type of "//TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Gets the sparsity type for interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_GET(INTERFACE_EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: SPARSITY_TYPE !<On exit, the sparsity type of the interface equations. \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        SPARSITY_TYPE=INTERFACE_EQUATIONS%SPARSITY_TYPE
      ELSE
        CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_SET(INTERFACE_EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The sparsity type to set \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(SPARSITY_TYPE)
        CASE(INTERFACE_EQUATIONS_SPARSE_MATRICES)
          INTERFACE_EQUATIONS%SPARSITY_TYPE=INTERFACE_EQUATIONS_SPARSE_MATRICES
        CASE(INTERFACE_EQUATIONS_FULL_MATRICES)
          INTERFACE_EQUATIONS%SPARSITY_TYPE=INTERFACE_EQUATIONS_FULL_MATRICES
        CASE DEFAULT
          LOCAL_ERROR="The specified sparsity type of "//TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the linearity type for interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_LINEARITY_TYPE_GET(INTERFACE_EQUATIONS,LINEARITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to get the linearity for
    INTEGER(INTG), INTENT(OUT) :: LINEARITY_TYPE !<On exit, the linearity type of the interface equations. \see INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes,INTERFACE_CONDITIONS_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("INTERFACE_EQUATIONS_LINEARITY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        LINEARITY_TYPE=INTERFACE_EQUATIONS%LINEARITY
      ELSE
        CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_LINEARITY_TYPE_GET")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_LINEARITY_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_LINEARITY_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_LINEARITY_TYPE_SET(INTERFACE_EQUATIONS,LINEARITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: LINEARITY_TYPE !<The linearity type to set \see INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes,INTERFACE_CONDITIONS_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("INTERFACE_EQUATIONS_LINEARITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(LINEARITY_TYPE)
        CASE(INTERFACE_CONDITION_LINEAR)
          INTERFACE_EQUATIONS%LINEARITY=INTERFACE_CONDITION_LINEAR
        CASE(INTERFACE_CONDITION_NONLINEAR)
          INTERFACE_EQUATIONS%LINEARITY=INTERFACE_CONDITION_NONLINEAR
        CASE(INTERFACE_CONDITION_NONLINEAR_BCS)
          INTERFACE_EQUATIONS%LINEARITY=INTERFACE_CONDITION_NONLINEAR_BCS
        CASE DEFAULT
          LOCAL_ERROR="The specified linearity type of "//TRIM(NUMBER_TO_VSTRING(LINEARITY_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_EQUATIONS_LINEARITY_TYPE_SET")
    RETURN
999 ERRORSEXITS("INTERFACE_EQUATIONS_LINEARITY_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_LINEARITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the time dependence type for interface equations.
  SUBROUTINE InterfaceEquations_TimeDependenceTypeGet(INTERFACE_EQUATIONS,TIME_DEPENDENCE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: TIME_DEPENDENCE_TYPE !<On exit, the time dependence type of the interface equations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("InterfaceEquations_TimeDependenceTypeGet",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        TIME_DEPENDENCE_TYPE=INTERFACE_EQUATIONS%TIME_DEPENDENCE
      ELSE
        CALL FlagError("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("InterfaceEquations_TimeDependenceTypeGet")
    RETURN
999 ERRORSEXITS("InterfaceEquations_TimeDependenceTypeGet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_TimeDependenceTypeGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for interface equations.
  SUBROUTINE InterfaceEquationsTimeDependenceTypeSet(INTERFACE_EQUATIONS,TIME_DEPENDENCE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to set the linearity for
    INTEGER(INTG), INTENT(IN) :: TIME_DEPENDENCE_TYPE !<The time dependence type to set \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,INTERFACE_CONDITIONS_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("InterfaceEquationsTimeDependenceTypeSet",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TIME_DEPENDENCE_TYPE)
        CASE(INTERFACE_CONDITION_STATIC)
          INTERFACE_EQUATIONS%TIME_DEPENDENCE=INTERFACE_CONDITION_STATIC
        CASE(INTERFACE_CONDITION_QUASISTATIC)
          INTERFACE_EQUATIONS%TIME_DEPENDENCE=INTERFACE_CONDITION_QUASISTATIC
        CASE(INTERFACE_CONDITION_FIRST_ORDER_DYNAMIC)
          INTERFACE_EQUATIONS%TIME_DEPENDENCE=INTERFACE_CONDITION_FIRST_ORDER_DYNAMIC
        CASE(INTERFACE_CONDITION_SECOND_ORDER_DYNAMIC)
          INTERFACE_EQUATIONS%TIME_DEPENDENCE=INTERFACE_CONDITION_SECOND_ORDER_DYNAMIC
        CASE DEFAULT
          LOCAL_ERROR="The specified time dependence type of "//TRIM(NUMBER_TO_VSTRING(TIME_DEPENDENCE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("InterfaceEquationsTimeDependenceTypeSet")
    RETURN
999 ERRORSEXITS("InterfaceEquationsTimeDependenceTypeSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE InterfaceEquationsTimeDependenceTypeSet

  !
  !================================================================================================================================
  !

  SUBROUTINE InterfaceEquations_VariableInterpSetsNumberSet(INTERFACE_EQUATIONS,VARIABLE_INDEX, &
    & NUMBER_OF_GEOMETRIC_SETS,NUMBER_OF_DEPENDENT_SETS,NUMBER_OF_PENALTY_SETS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<The interface equations to set the interface interpolation sets for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_INDEX !<The variable index number to set the number of interpolation sets for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_GEOMETRIC_SETS !<The number of geometric interface interpolation sets to set
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DEPENDENT_SETS !<The number of dependent interface interpolation sets to set
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_PENALTY_SETS !<The number of dependent interface interpolation sets to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_DEPENDENT_TYPE), POINTER :: INTERFACE_DEPENDENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("InterfaceEquations_VariableInterpSetsNumberSet",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FlagError("Interface equations have already been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERPOLATION)) THEN
          IF(ALLOCATED(INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION)) THEN
            INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
            IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
              INTERFACE_DEPENDENT=>INTERFACE_CONDITION%DEPENDENT
              IF(ASSOCIATED(INTERFACE_DEPENDENT)) THEN
                IF(VARIABLE_INDEX>0.AND.VARIABLE_INDEX<=INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES) THEN
                  IF(NUMBER_OF_GEOMETRIC_SETS>0) THEN
                    IF(NUMBER_OF_DEPENDENT_SETS>0) THEN
                      INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(VARIABLE_INDEX)% &
                        & NUMBER_OF_GEOMETRIC_INTERPOLATION_SETS=NUMBER_OF_GEOMETRIC_SETS
                      INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(VARIABLE_INDEX)% &
                        & NUMBER_OF_DEPENDENT_INTERPOLATION_SETS=NUMBER_OF_DEPENDENT_SETS
                      INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(VARIABLE_INDEX)% &
                        & NUMBER_OF_PENALTY_INTERPOLATION_SETS=NUMBER_OF_PENALTY_SETS
                    ELSE
                      LOCAL_ERROR="The specified number of dependent sets of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DEPENDENT_SETS,"*",ERR,ERROR))// &
                        & " is invalid. The number of dependent sets must be > 0."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The specified number of geometric sets of "// &
                      & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GEOMETRIC_SETS,"*",ERR,ERROR))// &
                      & " is invalid. The number of geometric sets must be > 0."
                    CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The specified variable index of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_INDEX,"*",ERR,ERROR))// &
                    & " is invalid. The index needs to be > 0 and <= "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERFACE_DEPENDENT%NUMBER_OF_DEPENDENT_VARIABLES,"*",ERR,ERROR))//"."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Interface condition dependent is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Interface equations interface condition is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Interface equations interpolation variable interpolation is not allocated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Interface equations interpolation is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("InterfaceEquations_VariableInterpSetsNumberSet")
    RETURN
999 ERRORS("InterfaceEquations_VariableInterpSetsNumberSet",ERR,ERROR)
    EXITS("InterfaceEquations_VariableInterpSetsNumberSet")
    RETURN 1
    
  END SUBROUTINE InterfaceEquations_VariableInterpSetsNumberSet

  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface conditions.
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_GET(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface conditions to get the interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<On exit, a pointer to the interface equations in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("INTERFACE_CONDITION_EQUATIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          CALL FlagError("Interface equations is already associated.",ERR,ERROR,*999)
        ELSE
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          IF(.NOT.ASSOCIATED(INTERFACE_EQUATIONS)) &
            & CALL FlagError("Interface equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("INTERFACE_CONDITION_EQUATIONS_GET")
    RETURN
999 ERRORSEXITS("INTERFACE_CONDITION_EQUATIONS_GET",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_GET





  !
  !================================================================================================================================
  !
  
END MODULE INTERFACE_EQUATIONS_ROUTINES
