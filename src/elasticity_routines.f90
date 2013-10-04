!> \file
!> \author Chris Bradley
!> \brief This module handles all elasticity routines.
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

!> This module handles all elasticity class routines.
MODULE ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE CONTROL_LOOP_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LINEAR_ELASTICITY_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Elasticity_EquationsSetSpecificationSet

  PUBLIC ELASTICITY_FINITE_ELEMENT_CALCULATE

  PUBLIC ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE,ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  PUBLIC ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE,ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE

  PUBLIC ELASTICITY_EQUATIONS_SET_SETUP

  PUBLIC ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  PUBLIC ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC
  
  PUBLIC Elasticity_ProblemSpecificationSet

  PUBLIC ELASTICITY_PROBLEM_SETUP

  PUBLIC ELASTICITY_PRE_SOLVE,ELASTICITY_POST_SOLVE

  PUBLIC ELASTICITY_CONTROL_LOOP_PRE_LOOP,Elasticity_ControlLoopPostLoop

  PUBLIC ELASTICITY_LOAD_INCREMENT_APPLY

CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Sets the problem specification for an elasticity equation set class.
  SUBROUTINE Elasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("Elasticity_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LinearElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE DEFAULT
        localError="Equations set equation type "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for an elasticity equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

    CALL Exits("Elasticity_EquationsSetSpecificationSet")
    RETURN
999 CALL Errors("Elasticity_EquationsSetSpecificationSet",err,error)
    CALL Exits("Elasticity_EquationsSetSpecificationSet")
    RETURN 1
  END SUBROUTINE Elasticity_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for an elasticity class finite element equation set.
  SUBROUTINE ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LINEAR_ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for the given element number for an elasticity class finite element equation set.
  SUBROUTINE ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual and rhs vector for the given element number for an elasticity class finite element equation set.
  SUBROUTINE ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the strain field for an elasticity class finite element equation set.
  SUBROUTINE ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_STRAIN_CALCULATE

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for an elasticity class finite element equation set.
  SUBROUTINE ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Cannot pre-evaluate the residual for a linear equations set.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for an elasticity class finite element equation set.
  SUBROUTINE ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FLAG_ERROR("Cannot post-evaluate the residual for a linear equations set.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for an elasticity equations set class.
  SUBROUTINE ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LINEAR_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the solution method for an elasticity equation set class.
  SUBROUTINE ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LINEAR_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for an elasticity equation set class.
  SUBROUTINE ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LINEAR_ELASTICITY_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN
999 CALL ERRORS("ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR)
    CALL EXITS("ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN 1
  END SUBROUTINE ELASTICITY_EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for an elasticity problem class.
  SUBROUTINE Elasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification array
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    CALL Enters("Elasticity_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=2) THEN
        problemType=problemSpecification(2)
        SELECT CASE(problemType)
        CASE(PROBLEM_LINEAR_ELASTICITY_TYPE)
          CALL LinearElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_FINITE_ELASTICITY_TYPE)
          CALL FiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
          CALL FLAG_ERROR("Not implemented yet.",err,error,*999)
        CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
          CALL FiniteElasticity_ContactProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE DEFAULT
          localError="Problem type "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for an elasticity problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Elasticity problem specification requires a type.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated",err,error,*999)
    END IF

    CALL Exits("Elasticity_ProblemSpecificationSet")
    RETURN
999 CALL Errors("Elasticity_ProblemSpecificationSet",err,error)
    CALL Exits("Elasticity_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE Elasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for an elasticity problem class.
  SUBROUTINE ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_LINEAR_ELASTICITY_TYPE)
        CALL LINEAR_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
        CALL FLAG_ERROR("Not implemented yet.",ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        CALL FiniteElasticity_ContactProblemSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("ELASTICITY_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
  
  !>Performs pre-solve actions for an elasticity problem class.
  SUBROUTINE ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        !Do Nothing
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
        !Do Nothing
      CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        !Do Nothing
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_PRE_SOLVE")
    RETURN
999 CALL ERRORS("ELASTICITY_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("ELASTICITY_PRE_SOLVE")
    RETURN 1
    
  END SUBROUTINE ELASTICITY_PRE_SOLVE

  !
  !================================================================================================================================
  !
  
  !>Sets up the output type for an elasticity problem class.
  SUBROUTINE ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        !Do Nothing
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
        !Do Nothing
      CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        !Do Nothing
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("ELASTICITY_POST_SOLVE")
    RETURN
999 CALL ERRORS("ELASTICITY_POST_SOLVE",ERR,ERROR)
    CALL EXITS("ELASTICITY_POST_SOLVE")
    RETURN 1
  END SUBROUTINE ELASTICITY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE ELASTICITY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

    CALL ENTERS("ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
      CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"====== Starting time step",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"CURRENT_TIME          = ",CURRENT_TIME,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"TIME_INCREMENT        = ",TIME_INCREMENT,ERR,ERROR,*999)
        IF(DIAGNOSTICS1) THEN
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"====== Starting time step",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"CURRENT_TIME          = ",CURRENT_TIME,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"TIME_INCREMENT        = ",TIME_INCREMENT,ERR,ERROR,*999)
        ENDIF
        IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
          CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
        END IF
        SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
        CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
            !do nothing for now
        CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
            CALL FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
            & " is not valid for an elasticity problem class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        !do nothing
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE ELASTICITY_CONTROL_LOOP_PRE_LOOP
  
  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop
  SUBROUTINE Elasticity_ControlLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Elasticity_ControlLoopPostLoop",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      problem=>controlLoop%PROBLEM
      IF(ASSOCIATED(problem)) THEN
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          IF(.NOT.ALLOCATED(problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(problem%specification,1)<2) THEN
            CALL FlagError("Problem specification must have at least two entries for an elasticity problem.",err,error,*999)
          END IF
          SELECT CASE(problem%specification(2))
          CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE)
            !Do nothing
          CASE(PROBLEM_FINITE_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
            CALL FiniteElasticity_ControlLoadIncrementLoopPostLoop(controlLoop,err,error,*999)
          CASE DEFAULT
            localError="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",err,error))// &
              & " is not valid for a elasticity problem class."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("Elasticity_ControlLoopPostLoop")
    RETURN
999 CALL ERRORS("Elasticity_ControlLoopPostLoop",err,error)
    CALL EXITS("Elasticity_ControlLoopPostLoop")
    RETURN 1
  END SUBROUTINE Elasticity_ControlLoopPostLoop

  !
  !================================================================================================================================
  !

  !> Apply load increments for equations sets
  SUBROUTINE ELASTICITY_LOAD_INCREMENT_APPLY(EQUATIONS_SET,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("ELASTICITY_LOAD_INCREMENT_APPLY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FINITE_ELASTICITY_LOAD_INCREMENT_APPLY(EQUATIONS_SET,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
      CASE DEFAULT
        !Do nothing
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("ELASTICITY_LOAD_INCREMENT_APPLY")
    RETURN
999 CALL ERRORS("ELASTICITY_LOAD_INCREMENT_APPLY",ERR,ERROR)
    CALL EXITS("ELASTICITY_LOAD_INCREMENT_APPLY")
    RETURN 1

  END SUBROUTINE ELASTICITY_LOAD_INCREMENT_APPLY

  !
  !================================================================================================================================
  !

END MODULE ELASTICITY_ROUTINES

