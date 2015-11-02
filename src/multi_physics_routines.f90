!> \file
!> $Id: multi_physics_routines.f90 177 2009-04-20 
!> \authors Christian Michler, Jack Lee
!> \brief This module handles all multi physics routines.
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

!> This module handles all multi physics class routines.
MODULE MULTI_PHYSICS_ROUTINES

  USE BASE_ROUTINES
  USE DIFFUSION_ADVECTION_DIFFUSION_ROUTINES
  USE DIFFUSION_DIFFUSION_ROUTINES
  USE FINITE_ELASTICITY_DARCY_ROUTINES
  USE FINITE_ELASTICITY_FLUID_PRESSURE_ROUTINES
  USE FSI_ROUTINES
  USE BIOELECTRIC_FINITE_ELASTICITY_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE MULTI_COMPARTMENT_TRANSPORT_ROUTINES
  USE NAVIER_STOKES_EQUATIONS_ROUTINES
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

  PUBLIC MultiPhysics_FiniteElementJacobianEvaluate,MultiPhysics_FiniteElementResidualEvaluate
  
  PUBLIC MultiPhysics_EquationsSetSpecificationSet,MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE, &
    & MULTI_PHYSICS_EQUATIONS_SET_SETUP,MultiPhysics_EquationsSetSolnMethodSet, &
    & MultiPhysics_ProblemSpecificationSet,MULTI_PHYSICS_PROBLEM_SETUP, &
    & MULTI_PHYSICS_POST_SOLVE,MULTI_PHYSICS_PRE_SOLVE,MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP, &
    & MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a multi physics equation set class.
  SUBROUTINE MultiPhysics_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiPhysics_EquationsSetSpecificationSet",err,error,*999)

    !Not that in general, this routine is never used as most multi-physics problems
    !use standard equations sets and couples them, rather than having a special
    !multi-physics problem equations set

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a multiphysics equations set.", &
          & err,error,*999)
      ENDIF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
        CALL FinElasticityFluidPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
        CALL DiffusionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DiffusionAdvectionDiffusion_EquationsSetSpecSet(equationsSet,specification,err,error,*999)
      CASE DEFAULT
        localError="The second equations set specification of "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for a multi physics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF

    EXITS("MultiPhysics_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("MultiPhysics_EquationsSetSpecificationSet",err,error)
    EXITS("MultiPhysics_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE MultiPhysics_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a multi physics class finite element equation set.
  SUBROUTINE MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a "// &
          & "multi-physics class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
        CALL FinElasticityFluidPressure_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
        CALL DiffusionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DiffusionAdvectionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a multi physics class finite element equation set.
  SUBROUTINE MultiPhysics_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MultiPhysics_FiniteElementJacobianEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a "// &
          & "multi-physics class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MultiPhysics_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("MultiPhysics_FiniteElementJacobianEvaluate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vectors for the given element number for a multi physics class finite element equation set.
  SUBROUTINE MultiPhysics_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MultiPhysics_FiniteElementResidualEvaluate",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a "// &
          & "multi-physics class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
!         CALL ELASTICITY_DARCY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MultiPhysics_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("MultiPhysics_FiniteElementResidualEvaluate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MultiPhysics_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the equations set for a multi physics equations set class.
  SUBROUTINE MULTI_PHYSICS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MULTI_PHYSICS_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a "// &
          & "multi-physics class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
        CALL ELASTICITY_DARCY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
        CALL DIFFUSION_DIFFUSION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DiffusionAdvectionDiffusion_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equation set type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MULTI_PHYSICS_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_EQUATIONS_SET_SETUP
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a multi physics equation set class.
  SUBROUTINE MultiPhysics_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MultiPhysics_EquationsSetSolnMethodSet",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a "// &
          & "multi-physics class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_DARCY_TYPE)
        CALL FinElasticityFluidPressure_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_DIFFUSION_TYPE)
        CALL DiffusionDiffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set equation type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MultiPhysics_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("MultiPhysics_EquationsSetSolnMethodSet",ERR,ERROR)
    EXITS("MultiPhysics_EquationsSetSolnMethodSet")
    RETURN 1
    
  END SUBROUTINE MultiPhysics_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem type and subtype for a multi physics problem class.
  SUBROUTINE MultiPhysics_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType

    ENTERS("MultiPhysics_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)<2) THEN
        CALL FlagError("Multi physics problem specification requires at least two entries.",err,error,*999)
      ENDIF
      problemType=problemSpecification(2)
      SELECT CASE(problemType)
      CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
        CALL FiniteElasticityDarcy_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
        CALL FinElasticityFluidPressure_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL BioelectricFiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FSI_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
        CALL DiffusionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DiffusionAdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
        CALL MultiCompartmentTransport_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
      CASE DEFAULT
        localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
          & " is not valid for a multi physics problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("MultiPhysics_ProblemSpecificationSet")
    RETURN
999 ERRORS("MultiPhysics_ProblemSpecificationSet",err,error)
    EXITS("MultiPhysics_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE MultiPhysics_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the problem for a multi physics problem class.
  SUBROUTINE MULTI_PHYSICS_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MULTI_PHYSICS_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a multi physics problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
        CALL ELASTICITY_DARCY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
        CALL ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL BIOELECTRIC_FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999) 
      CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        !CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CALL FSI_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
        CALL DIFFUSION_DIFFUSION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
        CALL MULTI_COMPARTMENT_TRANSPORT_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MULTI_PHYSICS_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a multi physics problem class.
  SUBROUTINE MULTI_PHYSICS_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MULTI_PHYSICS_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a multi physics problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
        CALL ELASTICITY_DARCY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
        CALL ELASTICITY_FLUID_PRESSURE_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL BIOELECTRIC_FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FSI_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
        CALL DIFFUSION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        !CALL DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
        CALL MULTI_COMPARTMENT_TRANSPORT_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MULTI_PHYSICS_POST_SOLVE")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a multi physics problem class.
  SUBROUTINE MULTI_PHYSICS_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MULTI_PHYSICS_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a multi physics problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
        CALL ELASTICITY_DARCY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
        CALL ELASTICITY_FLUID_PRESSURE_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL BIOELECTRIC_FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
        CALL FlagError("Not implemented.",ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FSI_PRE_SOLVE(CONTROL_LOOP,SOLVER,Err,Error,*999)
      CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
        CALL DIFFUSION_DIFFUSION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        CALL DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
        CALL MULTI_COMPARTMENT_TRANSPORT_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MULTI_PHYSICS_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a multi physics problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
        CALL ELASTICITY_DARCY_CONTROL_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
        !do nothing
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL BioelectricFiniteElasticity_ControlLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
        !do nothing
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        !TODO Store previous data?
      CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
        !do nothing
      CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        !do nothing
      CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
        !do nothing
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a multi physics problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_FINITE_ELASTICITY_DARCY_TYPE)
        CALL ELASTICITY_DARCY_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE)
        !do nothing
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL BioelectricFiniteElasticity_ControlLoopPostLoop(CONTROL_LOOP,ERR,ERROR,*999)
      CASE(PROBLEM_FINITE_ELASTICITY_STOKES_TYPE)
        !do nothing
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        CALL FSI_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*999)
      CASE(PROBLEM_DIFFUSION_DIFFUSION_TYPE)
        !do nothing
      CASE(PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE)
        !do nothing
      CASE(PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE)
        !do nothing
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a multi physics problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

END MODULE MULTI_PHYSICS_ROUTINES

