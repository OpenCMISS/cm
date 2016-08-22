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

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Elasticity_EquationsSetSpecificationSet

  PUBLIC ELASTICITY_FINITE_ELEMENT_CALCULATE

  PUBLIC ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE,ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  PUBLIC Elasticity_FiniteElementPreResidualEvaluate,Elasticity_FiniteElementPostResidualEvaluate

  PUBLIC ELASTICITY_EQUATIONS_SET_SETUP

  PUBLIC Elasticity_EquationsSetSolutionMethodSet

  PUBLIC Elasticity_EquationsSetDerivedVariableCalculate

  PUBLIC Elasticity_StrainInterpolateXi

  PUBLIC Elasticity_BoundaryConditionsAnalyticCalculate
  
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

    ENTERS("Elasticity_EquationsSetSpecificationSet",err,error,*999)

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
        localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for an elasticity equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Elasticity_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Elasticity_EquationsSetSpecificationSet",err,error)
    EXITS("Elasticity_EquationsSetSpecificationSet")
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
    
    ENTERS("ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

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
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("ELASTICITY_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
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
    
    ENTERS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 ERRORSEXITS("ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
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
    
    ENTERS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 ERRORSEXITS("ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementPreResidualEvaluate(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("Elasticity_FiniteElementPreResidualEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FlagError("Cannot pre-evaluate the residual for a linear equations set.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_FiniteElementPreResidualEvaluate(EQUATIONS_SET,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Elasticity_FiniteElementPreResidualEvaluate")
    RETURN
999 ERRORSEXITS("Elasticity_FiniteElementPreResidualEvaluate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementPreResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for an elasticity class finite element equation set.
  SUBROUTINE Elasticity_FiniteElementPostResidualEvaluate(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("Elasticity_FiniteElementPostResidualEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL FlagError("Cannot post-evaluate the residual for a linear equations set.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_FiniteElementPostResidualEvaluate(EQUATIONS_SET,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Elasticity_FiniteElementPostResidualEvaluate")
    RETURN
999 ERRORS("Elasticity_FiniteElementPostResidualEvaluate",ERR,ERROR)
    EXITS("Elasticity_FiniteElementPostResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE Elasticity_FiniteElementPostResidualEvaluate

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
    
    ENTERS("ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

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
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the solution method for an elasticity equation set class.
  SUBROUTINE Elasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Elasticity_EquationsSetSolutionMethodSet",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LinearElasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Elasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Elasticity_EquationsSetSolutionMethodSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Elasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Calculates a derived value for the elasticity equations set. \see OPENCMISS::CMISSEquationsSet_DerivedCalculate
  SUBROUTINE Elasticity_EquationsSetDerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived field type to calculate. \see EQUATIONS_SET_CONSTANTS_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Elasticity_EquationsSetDerivedVariableCalculate",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) THEN
        CALL FlagError("Equations set has not been finished.",err,error,*999)
      ELSE
        SELECT CASE(equationsSet%specification(2))
        CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
          CALL FiniteElasticityEquationsSet_DerivedVariableCalculate(equationsSet,derivedType, &
            & err,error,*999)
        CASE DEFAULT
          CALL FlagError("The second equations set specification of "// &
            & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(2),"*",err,error))// &
            & " is not valid for an elasticity equations set.",err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Elasticity_EquationsSetDerivedVariableCalculate")
    RETURN
999 ERRORS("Elasticity_EquationsSetDerivedVariableCalculate",err,error)
    EXITS("Elasticity_EquationsSetDerivedVariableCalculate")
    RETURN 1
    
  END SUBROUTINE Elasticity_EquationsSetDerivedVariableCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the strain tensor at a given element xi location.
  SUBROUTINE Elasticity_StrainInterpolateXi(equationsSet,userElementNumber,xi,values,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate strain for.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(IN) :: xi(:) !<The element xi to interpolate the field at.
    REAL(DP), INTENT(OUT) :: values(6) !<The interpolated strain tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("Elasticity_StrainInterpolateXi",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    SELECT CASE(equationsSet%specification(2))
    CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
      CALL FiniteElasticity_StrainInterpolateXi(equationsSet,userElementNumber,xi,values,err,error,*999)
    CASE DEFAULT
      CALL FlagError("The second equations set specification of "// &
        & TRIM(NumberToVstring(equationsSet%specification(2),"*",err,error))// &
        & " is not valid for an elasticity equation set.",err,error,*999)
    END SELECT

    EXITS("Elasticity_StrainInterpolateXi")
    RETURN
999 ERRORSEXITS("Elasticity_StrainInterpolateXi",err,error)
    RETURN 1
  END SUBROUTINE Elasticity_StrainInterpolateXi

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for an elasticity equation set class.
  SUBROUTINE Elasticity_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("Elasticity_BoundaryConditionsAnalyticCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for an elasticity class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TYPE)
        CALL LinearElasticity_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_TYPE)
        CALL FiniteElasticity_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("Elasticity_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("Elasticity_BoundaryConditionsAnalyticCalculate",ERR,ERROR)
    EXITS("Elasticity_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Elasticity_BoundaryConditionsAnalyticCalculate

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
          localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for an elasticity problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Elasticity problem specification requires a type.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
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
    
    ENTERS("ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

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
        CALL FlagError("Not implemented yet.",ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        CALL FiniteElasticity_ContactProblemSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for an elasticity problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
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
    
    ENTERS("ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

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
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("ELASTICITY_PRE_SOLVE",ERR,ERROR)
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
    
    ENTERS("ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

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
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_POST_SOLVE")
    RETURN
999 ERRORSEXITS("ELASTICITY_POST_SOLVE",ERR,ERROR)
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

    ENTERS("ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

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
            CALL FiniteElasticity_ControlTimeLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
            & " is not valid for an elasticity problem class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        !do nothing
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("ELASTICITY_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 ERRORSEXITS("ELASTICITY_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
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

    ENTERS("Elasticity_ControlLoopPostLoop",err,error,*999)

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
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Elasticity_ControlLoopPostLoop")
    RETURN
999 ERRORSEXITS("Elasticity_ControlLoopPostLoop",err,error)
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

    ENTERS("ELASTICITY_LOAD_INCREMENT_APPLY",ERR,ERROR,*999)

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
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("ELASTICITY_LOAD_INCREMENT_APPLY")
    RETURN
999 ERRORSEXITS("ELASTICITY_LOAD_INCREMENT_APPLY",ERR,ERROR)
    RETURN 1

  END SUBROUTINE ELASTICITY_LOAD_INCREMENT_APPLY

  !
  !================================================================================================================================
  !

END MODULE ELASTICITY_ROUTINES

