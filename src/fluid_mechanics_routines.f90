!> \file
!> \author Sebastian Krittian
!> \brief This module handles all fluid mechanics routines.
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
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): David Ladd
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

!> This module handles all fluid mechanics class routines.
MODULE FLUID_MECHANICS_ROUTINES

  USE BASE_ROUTINES
  USE BURGERS_EQUATION_ROUTINES
  USE CHARACTERISTIC_EQUATION_ROUTINES
  USE CONTROL_LOOP_ROUTINES
  USE DARCY_EQUATIONS_ROUTINES
  USE DARCY_PRESSURE_EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STOKES_EQUATIONS_ROUTINES
  USE NAVIER_STOKES_EQUATIONS_ROUTINES
  USE POISEUILLE_EQUATIONS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE TYPES


  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE,FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE
  PUBLIC FluidMechanics_NodalJacobianEvaluate,FluidMechanics_NodalResidualEvaluate

  PUBLIC FluidMechanics_EquationsSetSpecificationSet,FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE, &
    & FLUID_MECHANICS_EQUATIONS_SET_SETUP,FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC, &
    & FluidMechanics_ProblemSpecificationSet,FLUID_MECHANICS_PROBLEM_SETUP, &
    & FLUID_MECHANICS_POST_SOLVE,FLUID_MECHANICS_PRE_SOLVE,FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP, &
    & FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP,FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluate the analytic solution for a fluid mechanics equations set.
  SUBROUTINE FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_TYPE,ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS, &
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

    CALL ENTERS("FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_TYPE)
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,POSITION, &
          & TANGENTS,NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS, &
          & MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_TYPE,"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_ANALYTIC_FUNCTIONS_EVALUATE

   !
  !================================================================================================================================
  !

  !>Sets the problem specification for a fluid mechanics equation set class.
  SUBROUTINE FluidMechanics_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL Enters("FluidMechanics_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL Stokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL NavierStokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL DarcyEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL DarcyPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL PoiseuilleEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BurgersEquation_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL Characteristic_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE DEFAULT
        localError="Equations set equation type "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

    CALL Exits("FluidMechanics_EquationsSetSpecificationSet")
    RETURN
999 CALL Errors("FluidMechanics_EquationsSetSpecificationSet",err,error)
    CALL Exits("FluidMechanics_EquationsSetSpecificationSet")
    RETURN 1
  END SUBROUTINE FluidMechanics_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a fluid mechanics class finite element equation set.
  SUBROUTINE FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL STOKES_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("There are no finite element matrices to be calculated for Navier-Stokes equation.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL FLAG_ERROR("There is no element stiffness matrix to be calculated for Darcy pressure.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL FLAG_ERROR("There are no finite element matrices to be calculated for Burgers equation.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL FLAG_ERROR("There are no finite element matrices to be calculated for Characteristic equations.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a fluid mechanics class finite element equation set.
  SUBROUTINE FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("There is no Jacobian to be evaluated for Stokes.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL FLAG_ERROR("There is no Jacobian to be evaluated for Poiseuille.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a fluid mechanics class finite element equation set.
  SUBROUTINE FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("There is no residual to be evaluated for Stokes.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL DARCY_PRESSURE_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL FLAG_ERROR("There is no residual to be evaluated for Poiseuille.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the nodal Jacobian matrix for the given node number for a fluid mechanics class nodal equation set.
  SUBROUTINE FluidMechanics_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("FluidMechanics_NodalJacobianEvaluate",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL Characteristic_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CASE DEFAULT
        localError="Equations set type "//TRIM(NUMBER_TO_VSTRING(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for a fluid mechanics equation set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",err,error,*999)
    ENDIF
       
    CALL EXITS("FluidMechanics_NodalJacobianEvaluate")
    RETURN
999 CALL ERRORS("FluidMechanics_NodalJacobianEvaluate",err,error)
    CALL EXITS("FluidMechanics_NodalJacobianEvaluate")
    RETURN 1
  END SUBROUTINE FluidMechanics_NodalJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the nodal residual and rhs vectors for the given node number for a fluid mechanics class nodal equation set.
  SUBROUTINE FluidMechanics_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("FluidMechanics_NodalResidualEvaluate",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL Characteristic_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CASE DEFAULT
        localError="Equations set type "//TRIM(NUMBER_TO_VSTRING(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for a fluid mechanics equation set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",err,error,*999)
    ENDIF
       
    CALL EXITS("FluidMechanics_NodalResidualEvaluate")
    RETURN
999 CALL ERRORS("FluidMechanics_NodalResidualEvaluate",err,error)
    CALL EXITS("FluidMechanics_NodalResidualEvaluate")
    RETURN 1
  END SUBROUTINE FluidMechanics_NodalResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a fluid mechanics equations set class.
  SUBROUTINE FLUID_MECHANICS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL STOKES_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL DARCY_PRESSURE_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL Characteristic_EquationsSet_Setup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_EQUATIONS_SET_SETUP
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a fluid mechanics equation set class.
  SUBROUTINE FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_EQUATIONS_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL STOKES_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL DARCY_PRESSURE_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL Characteristic_EquationsSet_SolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechancis equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets the analytic boundary conditions for a fluid mechanics equation set class.
  SUBROUTINE FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditionsn to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a fluid mechanics class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_STOKES_EQUATION_TYPE)
        CALL STOKES_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE)
        CALL NavierStokes_AnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_EQUATIONS_BOUNDARY_CONDITIONS_ANALYTIC

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a fluid mechanics problem class.
  SUBROUTINE FluidMechanics_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    CALL Enters("FluidMechanics_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=2) THEN
        problemType=problemSpecification(2)
        SELECT CASE(problemType)
        CASE(PROBLEM_STOKES_EQUATION_TYPE)
          CALL Stokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
          CALL NavierStokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_DARCY_EQUATION_TYPE)
          CALL DarcyEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
          CALL PoiseuilleEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE(PROBLEM_BURGERS_EQUATION_TYPE)
          CALL BurgersEquation_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
        CASE DEFAULT
          localError="Problem equation type "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a fluid mechanics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Fluid mechanics problem specification must have a type set.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated",err,error,*999)
    END IF

    CALL Exits("FluidMechanics_ProblemSpecificationSet")
    RETURN
999 CALL Errors("FluidMechanics_ProblemSpecificationSet",err,error)
    CALL Exits("FluidMechanics_ProblemSpecificationSet")
    RETURN 1
  END SUBROUTINE FluidMechanics_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a fluid mechanics problem class.
  SUBROUTINE FLUID_MECHANICS_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_STOKES_EQUATION_TYPE)
        CALL STOKES_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a fluid mechanics problem class.
  SUBROUTINE FLUID_MECHANICS_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(control_loop%problem%specification,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_STOKES_EQUATION_TYPE)
        CALL STOKES_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_POST_SOLVE")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_POST_SOLVE",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_POST_SOLVE")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_POST_SOLVE

  !
  !================================================================================================================================


  !>Sets up the output type for a fluid mechanics problem class.
  SUBROUTINE FLUID_MECHANICS_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FLUID_MECHANICS_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(control_loop%problem%specification,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_STOKES_EQUATION_TYPE)
        CALL STOKES_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
        CALL NAVIER_STOKES_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DARCY_EQUATION_TYPE)
        CALL DARCY_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
        CALL POISEUILLE_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_BURGERS_EQUATION_TYPE)
        CALL BURGERS_EQUATION_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a fluid mechanics problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("FLUID_MECHANICS_PRE_SOLVE")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
      CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
        IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(control_loop%problem%specification,1)<2) THEN
          CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
        END IF
        SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
        CASE(PROBLEM_STOKES_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_DARCY_EQUATION_TYPE)
          CALL DARCY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
        CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_BURGERS_EQUATION_TYPE)
          !do nothing
        CASE DEFAULT
          LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
            & " is not valid for a fluid mechanics problem class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        !do nothing
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
      CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
        IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(control_loop%problem%specification,1)<2) THEN
          CALL FlagError("Problem specification must have at least two entries for a fluid mechanics problem.",err,error,*999)
        END IF
        SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
        CASE(PROBLEM_STOKES_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_NAVIER_STOKES_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_DARCY_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_POISEUILLE_EQUATION_TYPE)
          !do nothing
        CASE(PROBLEM_BURGERS_EQUATION_TYPE)
          !do nothing
        CASE DEFAULT
          LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
            & " is not valid for a fluid mechanics problem class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        !do nothing
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP")
    RETURN
999 CALL ERRORS("FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    CALL EXITS("FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP")
    RETURN 1
  END SUBROUTINE FLUID_MECHANICS_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

END MODULE FLUID_MECHANICS_ROUTINES

