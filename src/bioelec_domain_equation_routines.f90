!> \file
!> $Id: bioelec_domain_equation_routines.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief This module handles all bioelectric domain equation routines.
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
!> The Original Code is openCMISS
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

!>This module handles all bioelectric domain equation routines.
MODULE BIOELEC_DOMAIN_EQUATION_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLUTIONS_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP,BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET, &
    & BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE,BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric domain equation type of a bioelectric equations set class.
  SUBROUTINE BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a bioelectric domain equation on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
!!\todo sort out mono/bidomain
      SELECT CASE(SETUP_TYPE)
      CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FIXED_CONDITIONS_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
          & " is invalid for a bioelectric domain equation."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a bioelectric domain equation type of a bioelectric equations set class.
  SUBROUTINE BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%TYPE)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET_SUBTYPE)
        CASE(EQUATIONS_SET_NO_SUBTYPE)
          EQUATIONS_SET%CLASS=EQUATIONS_SET_BIOELECTRICS_CLASS
          EQUATIONS_SET%TYPE=EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE
          EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SUBTYPE
          CALL BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE, &
            & EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a monodomain equation type of a bioelectric equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        SELECT CASE(EQUATIONS_SET_SUBTYPE)
        CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
          EQUATIONS_SET%CLASS=EQUATIONS_SET_BIOELECTRICS_CLASS
          EQUATIONS_SET%TYPE=EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE
          EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE
          CALL BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE, &
            & EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
          EQUATIONS_SET%CLASS=EQUATIONS_SET_BIOELECTRICS_CLASS
          EQUATIONS_SET%TYPE=EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE
          EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE
          CALL BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE, &
            & EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a bidomain equation type of a bioelectric equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TYPE,"*",ERR,ERROR))// &
          & " is not valid for a bioelectric equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE BIOELEC_DOMAIN_EQUATION_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectric domain solution.
  SUBROUTINE BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the solutions set to setup a bioelectric domain equation on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
!!\todo sort out mono/bi domain
      SELECT CASE(SETUP_TYPE)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_DO_ACTION)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
        CASE(PROBLEM_SETUP_DO_ACTION)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric domain equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLUTION_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
        CASE(PROBLEM_SETUP_DO_ACTION)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_TYPE)
        SELECT CASE(ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
        CASE(PROBLEM_SETUP_DO_ACTION)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
            & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a bioelectric equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
          & " is invalid for a bioelectric domain equation."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE BIOELEC_DOMAIN_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectric domain equation finite element equations set.
  SUBROUTINE BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%TYPE)
        CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
          CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
          CASE DEFAULT
            LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The equations set type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE BIOELEC_DOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !
   
END MODULE BIOELEC_DOMAIN_EQUATION_ROUTINES
