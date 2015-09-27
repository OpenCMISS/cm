!> \file
!> \author Chris Bradley
!> \brief This module handles all classical field routines.
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

!> This module handles all classical field class routines.
MODULE CLASSICAL_FIELD_ROUTINES

  USE ADVECTION_EQUATION_ROUTINES
  USE ADVECTION_DIFFUSION_EQUATION_ROUTINES
  USE BASE_ROUTINES
  USE DIFFUSION_EQUATION_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE HELMHOLTZ_EQUATIONS_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE HAMILTON_JACOBI_EQUATIONS_ROUTINES
  USE LAPLACE_EQUATIONS_ROUTINES
  USE POISSON_EQUATIONS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE REACTION_DIFFUSION_EQUATION_ROUTINES
  USE STRINGS
  USE TYPES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE
  
  PUBLIC CLASSICAL_FIELD_CONTROL_LOOP_POST_LOOP

  PUBLIC ClassicalField_FiniteElementJacobianEvaluate,ClassicalField_FiniteElementResidualEvaluate
  
  PUBLIC ClassicalField_EquationsSetSpecificationSet

  PUBLIC CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE
  
  PUBLIC CLASSICAL_FIELD_EQUATIONS_SET_SETUP

  PUBLIC ClassicalField_EquationsSetSolutionMethodSet

  PUBLIC ClassicalField_ProblemSpecificationSet

  PUBLIC CLASSICAL_FIELD_PROBLEM_SETUP

  PUBLIC CLASSICAL_FIELD_PRE_SOLVE,CLASSICAL_FIELD_POST_SOLVE

  PUBLIC ClassicalField_BoundaryConditionsAnalyticCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluate the analytic solution for a classical field equations set.
  SUBROUTINE CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_TYPE,ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS, &
    & NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the analytic for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_TYPE !<The type of equation to evaluate
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: POSITION(:) !<POSITION(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: TANGENTS(:,:) !<TANGENTS(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: NORMAL(:) !<NORMAL(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: TIME !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: GLOBAL_DERIVATIVE !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: ANALYTIC_PARAMETERS(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: MATERIALS_PARAMETERS(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the analtyic function value.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
     SELECT CASE(EQUATIONS_TYPE)
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,POSITION, &
          & TANGENTS,NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS, &
          & MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a classical field equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop for bioelectric problems, i.e., after each time step for a time loop
  SUBROUTINE CLASSICAL_FIELD_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("CLASSICAL_FIELD_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      PROBLEM=>CONTROL_LOOP%PROBLEM
      IF(ASSOCIATED(PROBLEM)) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          SELECT CASE(PROBLEM%specification(2))
          CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
            CALL DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The second problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(2),"*",ERR,ERROR))// &
              & " is not valid for a classical field problem."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("CLASSICAL_FIELD_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_LOOP_POST_LOOP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CLASSICAL_FIELD_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

  !>Sets the equations set specification for a classical field equation set.
  SUBROUTINE ClassicalField_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("ClassicalField_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL Laplace_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL HJEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL Poisson_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL Helmholtz_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL Advection_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL AdvectionDiffusion_EquationEquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL ReactionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for a classical field equations set."
       CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

    EXITS("ClassicalField_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("ClassicalField_EquationsSetSpecificationSet",err,error)
    EXITS("ClassicalField_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ClassicalField_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a clasical field class finite element equation set.
  SUBROUTINE CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL HJ_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL POISSON_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL AdvectionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL ReactionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a clasical field class finite element equation set.
  SUBROUTINE ClassicalField_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("ClassicalField_FiniteElementJacobianEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL Poisson_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ClassicalField_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("ClassicalField_FiniteElementJacobianEvaluate",ERR,ERROR)
    EXITS("ClassicalField_FiniteElementJacobianEvaluate")
    RETURN 1
    
  END SUBROUTINE ClassicalField_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a clasical field class finite element equation set.
  SUBROUTINE ClassicalField_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("ClassicalField_FiniteElementResidualEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL Poisson_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL ADVECTION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ClassicalField_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("ClassicalField_FiniteElementResidualEvaluate",ERR,ERROR)
    EXITS("ClassicalField_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE ClassicalField_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a classical field equations set class.
  SUBROUTINE CLASSICAL_FIELD_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CLASSICAL_FIELD_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL LAPLACE_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL HJ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL POISSON_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL AdvectionDiffusion_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL Advection_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL ReactionDiffusion_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CLASSICAL_FIELD_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLASSICAL_FIELD_EQUATIONS_SET_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a classical field equation set class.
  SUBROUTINE ClassicalField_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CLASSICAL_FIELD_EQUATIONS_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL Laplace_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL HJEquation_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL Poisson_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL Helmholtz_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL AdvectionDiffusion_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_EQUATION_TYPE)
        CALL Advection_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL ReactionDiffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ClassicalField_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("ClassicalField_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("ClassicalField_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE ClassicalField_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for a classical field equation set class.
  SUBROUTINE ClassicalField_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("ClassicalField_BoundaryConditionsAnalyticCalculate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a classical field class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
        CALL Laplace_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HJ_EQUATION_TYPE)
        CALL HJ_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISSON_EQUATION_TYPE)
        CALL Poisson_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE)
        CALL Helmholtz_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TYPE)
        CALL Diffusion_BoundaryConditionAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL AdvectionDiffusion_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    EXITS("ClassicalField_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("ClassicalField_BoundaryConditionsAnalyticCalculate",ERR,ERROR)
    EXITS("ClassicalField_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE ClassicalField_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a classical field problem class.
  SUBROUTINE ClassicalField_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    ENTERS("ClassicalField_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=2) THEN
        problemType=problemSpecification(2)
        SELECT CASE(problemType)
        CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
          CALL Laplace_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_HJ_EQUATION_TYPE)
          CALL HJEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_POISSON_EQUATION_TYPE)
          CALL Poisson_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
          CALL Helmholtz_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_WAVE_EQUATION_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
          CALL Diffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
          CALL Advection_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
          CALL AdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
          CALL ReactionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a classical field problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Classical field problem specification must have a type set.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("ClassicalField_ProblemSpecificationSet")
    RETURN
999 ERRORS("ClassicalField_ProblemSpecificationSet",err,error)
    EXITS("ClassicalField_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE ClassicalField_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a classical field problem class.
  SUBROUTINE CLASSICAL_FIELD_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CLASSICAL_FIELD_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a classical field problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
        CALL LAPLACE_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        CALL HJ_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_POISSON_EQUATION_TYPE)
        CALL POISSON_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
        CALL HELMHOLTZ_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_WAVE_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
        CALL ADVECTION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL REACTION_DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CLASSICAL_FIELD_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLASSICAL_FIELD_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a classical field problem class.
  SUBROUTINE CLASSICAL_FIELD_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CLASSICAL_FIELD_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a classical field problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
        !CALL LAPLACE_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        !CALL HJ_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_POISSON_EQUATION_TYPE)
        CALL POISSON_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
        !CALL HELMHOLTZ_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_WAVE_EQUATION_TYPE)
        !CALL WAVE_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL DIFFUSION_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
        CALL ADVECTION_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL ADVECTION_DIFFUSION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL REACTION_DIFFUSION_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
        !CALL BIHARMONIC_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CLASSICAL_FIELD_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLASSICAL_FIELD_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a classical field problem class.
  SUBROUTINE CLASSICAL_FIELD_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CLASSICAL_FIELD_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a classical field problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_LAPLACE_EQUATION_TYPE)
        !CALL LAPLACE_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_HJ_EQUATION_TYPE)
        !CALL HJ_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_POISSON_EQUATION_TYPE)
        CALL POISSON_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_HELMHOLTZ_EQUATION_TYPE)
        !CALL HELMHOLTZ_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_WAVE_EQUATION_TYPE)
        !CALL WAVE_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_EQUATION_TYPE)
        CALL DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_EQUATION_TYPE)
        !CALL ADVECTION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE)
        CALL ADVECTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
        CALL REACTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_BIHARMONIC_EQUATION_TYPE)
        !CALL BIHARMONIC_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a classical field problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CLASSICAL_FIELD_POST_SOLVE")
    RETURN
999 ERRORSEXITS("CLASSICAL_FIELD_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLASSICAL_FIELD_POST_SOLVE


  !
  !================================================================================================================================

END MODULE CLASSICAL_FIELD_ROUTINES

