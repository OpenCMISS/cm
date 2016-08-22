!> \file  
!> \author Sebastian Krittian
!> \brief This module handles all Navier-Stokes fluid routines.
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
!> Contributor(s): David Ladd, Soroush Safaei
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

!>This module handles all Navier-Stokes fluid routines.
MODULE NAVIER_STOKES_EQUATIONS_ROUTINES

  USE ADVECTION_EQUATION_ROUTINES
  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CHARACTERISTIC_EQUATION_ROUTINES
  USE CMISS_MPI  
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LAPACK
  USE MATHS
  USE MATRIX_VECTOR
  USE MESH_ROUTINES
  USE MPI
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STREE_EQUATION_ROUTINES
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  PUBLIC NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE
  
  PUBLIC NavierStokes_EquationsSetSpecificationSet
  
  PUBLIC NavierStokes_EquationsSetSolutionMethodSet
  
  PUBLIC NAVIER_STOKES_EQUATIONS_SET_SETUP
  
  PUBLIC NavierStokes_PreSolveALEUpdateParameters
  
  PUBLIC NavierStokes_PreSolveUpdateBoundaryConditions

  PUBLIC NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH

  PUBLIC NAVIER_STOKES_PRE_SOLVE

  PUBLIC NAVIER_STOKES_POST_SOLVE

  PUBLIC NavierStokes_ProblemSpecificationSet
  
  PUBLIC NAVIER_STOKES_PROBLEM_SETUP
  
  PUBLIC NavierStokes_FiniteElementResidualEvaluate
  
  PUBLIC NavierStokes_FiniteElementJacobianEvaluate

  PUBLIC NavierStokes_BoundaryConditionsAnalyticCalculate
  
  PUBLIC NavierStokes_ResidualBasedStabilisation
  
  PUBLIC NavierStokes_Couple1D0D
  
  PUBLIC NavierStokes_CoupleCharacteristics
  
  PUBLIC NavierStokes_FiniteElementPreResidualEvaluate
  
  PUBLIC NavierStokes_ControlLoopPostLoop
  
  PUBLIC NavierStokes_UpdateMultiscaleBoundary

CONTAINS 

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Navier-Stokes flow equation type of an fluid mechanics equations set class.
  SUBROUTINE NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NavierStokes_EquationsSetSolutionMethodSet",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_NODAL_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Navier-Stokes flow equation type of a fluid mechanics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokes_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("NavierStokes_EquationsSetSolutionMethodSet",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NavierStokes_EquationsSetSolutionMethodSet

!
!================================================================================================================================
!

  !>Sets the equation specification for a Navier-Stokes fluid type of a fluid mechanics equations set class.
  SUBROUTINE NavierStokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: specification(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("NavierStokes_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      ENDIF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_Coupled1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
        !ok
      CASE(EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE)
        CALL FlagError("Not implemented yet.",err,error,*999)
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      ENDIF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokes_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("NavierStokes_EquationsSetSpecificationSet",err,error)
    EXITS("NavierStokes_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE NavierStokes_EquationsSetSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Navier-Stokes fluid setup.
  SUBROUTINE NAVIER_STOKES_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD,ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    INTEGER(INTG) :: GEOMETRIC_SCALING_TYPE,GEOMETRIC_MESH_COMPONENT,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_ANALYTIC_COMPONENTS,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,GEOMETRIC_COMPONENT_NUMBER,I,compIdx,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: MATERIAL_FIELD_NUMBER_OF_VARIABLES,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
    INTEGER(INTG) :: elementBasedComponents,nodeBasedComponents,constantBasedComponents
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NAVIER_STOKES_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS_EQUATIONS_SET_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_FIELD)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
              EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                & "*",ERR,ERROR))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP% &
                & SETUP_TYPE,"*",ERR,ERROR))// " is not implemented for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
              EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1 
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 1
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF (EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Penalty Coefficient",ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                !Default the penalty coefficient value to 1E4
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0E4_DP,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                & "*",ERR,ERROR))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP% &
                & SETUP_TYPE,"*",ERR,ERROR))// " is not implemented for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR,ERROR,*999)
              EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 3
              nodeBasedComponents = 1  ! boundary flux
              elementBasedComponents = 10  ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, C1
              constantBasedComponents = 4  ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                   & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                 CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],ERR,ERROR,*999)
                 CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & "BoundaryFlow",ERR,ERROR,*999)
                 CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_V_VARIABLE_TYPE, &
                   & "ElementMetrics",ERR,ERROR,*999)
                 CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_U1_VARIABLE_TYPE, &
                   & "EquationsConstants",ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_V_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U1_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & FIELD_U_VARIABLE_TYPE,nodeBasedComponents,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & FIELD_V_VARIABLE_TYPE,elementBasedComponents,ERR,ERROR,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & FIELD_U1_VARIABLE_TYPE,constantBasedComponents,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="User-specified fields are not yet implemented for an equations set field field &
                  & setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP% &
                  & SETUP_TYPE,"*",ERR,ERROR))// " for a Navier-Stokes fluid."
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                !Default the Element Metrics parameter values 0.0
                nodeBasedComponents = 1  ! boundary flux
                elementBasedComponents = 10  ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, C1
                constantBasedComponents = 4  ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
                ! Init boundary flux to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                ! Init Element Metrics to 0 (except C1)
                DO compIdx=1,elementBasedComponents-1
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,compIdx,0.0_DP,ERR,ERROR,*999)
                END DO
                ! Default C1 to -1 for now, will be calculated in ResidualBasedStabilisation if not specified by user
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elementBasedComponents,-1.0_DP,ERR,ERROR,*999)
                ! Boundary stabilisation scale factor (beta): default to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                ! Max Courant (CFL) number: default to 1.0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,1.0_DP,ERR,ERROR,*999)
                ! Init Time increment to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3,0.0_DP,ERR,ERROR,*999)
                ! Stabilisation type: default to 1 for RBS (0=none, 2=RBVM)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,1.0_DP,ERR,ERROR,*999)
              ELSE
                LOCAL_ERROR="User-specified fields are not yet implemented for an equations set field field &
                  & setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP% &
                  & SETUP_TYPE,"*",ERR,ERROR))// " for a Navier-Stokes fluid."
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                & "*",ERR,ERROR))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP% &
                & SETUP_TYPE,"*",ERR,ERROR))// " is not implemented for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            !Do nothing???
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 1 
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,& 
                  & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                DO compIdx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,compIdx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                END DO
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              nodeBasedComponents = 1  ! boundary flux
              elementBasedComponents = 10  ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, C1
              constantBasedComponents = 4  ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
              IF (EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                  & ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ! Element-based fields
                DO compIdx = 1, elementBasedComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,compIdx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,compIdx,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                END DO
                ! Constant fields: boundary stabilisation scale factor and max courant #
                DO compIdx = 1, constantBasedComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,compIdx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,compIdx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                END DO
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                !Do nothing
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF (EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "del U/del n",ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension and one for pressure
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                DO compIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !compIdx
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO compIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !compIdx
                  !Default geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                !Other solutions not defined yet
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE 
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension and one for pressure
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO compIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                      & compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                      & compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO !compIdx
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                    & "*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                 & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                 & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_PRESSURE_VALUES_SET_TYPE,ERR,ERROR,*999)
                  DO compIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_PRESSURE_VALUES_SET_TYPE,compIdx,0.0_DP,ERR,ERROR,*999)
                  END DO
                END IF
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE, & 
                & "*",ERR,ERROR))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
               EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !set number of variables to 5 (U,DELUDELN,V,U1,U2)
              DEPENDENT_FIELD_NUMBER_OF_VARIABLES=5
              !calculate number of components (Q,A) for U and dUdN
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
                !define new created field to be dependent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)                  
                !set number of variables to 6 (U,DELUDELN,V,U1,U2)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                !set data type
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                  
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                !2 component (W1,W2) for V
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & NUMBER_OF_DIMENSIONS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                END DO
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & FIELD_U1_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & FIELD_U2_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE 
                !Check the user specified field- Characteristic equations
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES, &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)

                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components (Q,A) for U and dUdN
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                !2 component (W1,W2) for V
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                    & "*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE, & 
                & "*",ERR,ERROR))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF (EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
                !define new created field to be dependent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                !set number of variables to 5 (U,DELUDELN,V,U1,U2)
                DEPENDENT_FIELD_NUMBER_OF_VARIABLES=5
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_U3_VARIABLE_TYPE],ERR,ERROR,*999)
                ELSE
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                    & ERR,ERROR,*999)
                END IF
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components (Q,A)
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,1,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U2_VARIABLE_TYPE,1,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & NUMBER_OF_DIMENSIONS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                END DO
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U2_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE 
                !set number of variables to 5 (U,DELUDELN,V,U1,U2)
                DEPENDENT_FIELD_NUMBER_OF_VARIABLES=5
                !Check the user specified field- Characteristic equations
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components (Q,A) for U and dUdN
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FlagError("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                    & "*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE, & 
                & "*",ERR,ERROR))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF (EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% & 
                  & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !compIdx
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !compIdx
                  !Default geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                !Other solutions not defined yet
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT 
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS                  
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO !compIdx
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                    &"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              END IF    
              !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !set number of variables to 1 (W)
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              !normalDirection for wave relative to node for W1,W2
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
              !Create the auto created independent field
              IF (EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !start field creation with name 'INDEPENDENT_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR, & 
                  & *999)
                !define new created field to be independent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% & 
                  & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                !set number of variables to 1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !compIdx
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !compIdx
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT 
              ELSE
                !Check the user specified field- Characteristic equation
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES, &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              END IF    
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
         CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !set number of variables to 1
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              !normalDirection for wave relative to element for W1,W2
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=1
              !Create the auto created independent field
              IF (EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                ! Do nothing? independent field should be set up by characteristic equation routines
              ELSE
                !Check the user specified field- Characteristic equation
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES, &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              END IF    
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% & 
                  & GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                END DO !compIdx
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO !compIdx
                  !Default geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                !Other solutions not defined yet
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT 
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   DO compIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS                  
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO !compIdx
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                    &"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              END IF    
              !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF (EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
              IF (ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
                IF (EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF (ASSOCIATED(DEPENDENT_FIELD)) THEN
                    EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                    IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                      IF (EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                        GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                        IF (ASSOCIATED(GEOMETRIC_FIELD)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE
                            !Check that domain is 2D
                            IF (NUMBER_OF_DIMENSIONS/=2) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 2 geometric dimensions."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END IF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            NUMBER_OF_ANALYTIC_COMPONENTS=4
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                             & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN, &
                             & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART, &
                             & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                            !Check that this is a 1D equations set
                            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                              !Set analytic function type
                              EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                              !Set numbrer of components- Q,A (same as N-S depenedent field)
                              NUMBER_OF_ANALYTIC_COMPONENTS=2
                            ELSE
                              LOCAL_ERROR="The third equations set specification must by a TRANSIENT1D or COUPLED1D0D "// &
                                & "to use an analytic function of type "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))//"."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END IF
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                            !Check that domain is 2D/3D
                            IF(NUMBER_OF_DIMENSIONS<2 .OR. NUMBER_OF_DIMENSIONS>3) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 2 or 3 geometric dimensions."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END IF
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            !Set numbrer of components
                            NUMBER_OF_ANALYTIC_COMPONENTS=10
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN
                            !Check that domain is 2D
                            IF (NUMBER_OF_DIMENSIONS/=2) THEN
                              LOCAL_ERROR="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                & " requires that there be 2 geometric dimensions."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            END IF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            NUMBER_OF_ANALYTIC_COMPONENTS=2
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1
                          CASE DEFAULT
                            LOCAL_ERROR="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                              & " is invalid for an analytic Navier-Stokes problem."
                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                          !Create analytic field if required
                          IF (NUMBER_OF_ANALYTIC_COMPONENTS>=1) THEN
                            IF (EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                              !Create the auto created analytic field
                              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                                & EQUATIONS_ANALYTIC%ANALYTIC_FIELD,ERR,ERROR,*999)
                              CALL FIELD_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,"Analytic Field",ERR,ERROR,*999)
                              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_INDEPENDENT_TYPE, &
                                & ERR,ERROR,*999)
                              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                                & ERR,ERROR,*999)
                              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                                & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,EQUATIONS_SET%GEOMETRY% &
                                & GEOMETRIC_FIELD,ERR,ERROR,*999)
                              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,1,ERR,ERROR,*999)
                              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,[FIELD_U_VARIABLE_TYPE], &
                                & ERR,ERROR,*999)
                              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & "Analytic",ERR,ERROR,*999)
                              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_DP_TYPE,ERR,ERROR,*999)
                              !Set the number of analytic components
                              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS,ERR,ERROR,*999)
                              !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                              DO compIdx=1,NUMBER_OF_ANALYTIC_COMPONENTS
                                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                  & compIdx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                                IF(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE == &
                                 & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                    & compIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                                ELSE
                                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                    & compIdx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                                END IF
                              END DO !compIdx
                              !Default the field scaling to that of the geometric field
                              CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                                & ERR,ERROR,*999)
                              CALL FIELD_SCALING_TYPE_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                            ELSE
                              !Check the user specified field
                              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                              IF (NUMBER_OF_ANALYTIC_COMPONENTS==1) THEN
                                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                              ELSE
                                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                              END IF
                              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                                & ERR,ERROR,*999)
                              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & NUMBER_OF_ANALYTIC_COMPONENTS,ERR,ERROR,*999)
                            END IF
                          END IF
                        ELSE
                          CALL FlagError("Equations set materials is not finished.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations set dependent field has not been finished.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations analytic is not associated.",ERR,ERROR,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
              IF (ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
                ANALYTIC_FIELD=>EQUATIONS_ANALYTIC%ANALYTIC_FIELD
                IF (ASSOCIATED(ANALYTIC_FIELD)) THEN
                  IF (EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                    !Finish creating the analytic field
                    CALL FIELD_CREATE_FINISH(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,ERR,ERROR,*999)
                    !Set the default values for the analytic field
                    SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
                      SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                        !Default the analytic parameter values (L, H, U_mean, Pout) to 0.0
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,2,0.0_DP,ERR,ERROR,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,3,0.0_DP,ERR,ERROR,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,4,0.0_DP,ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " is invalid for an analytical static Navier-Stokes equation."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                      SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
                        !Default the analytic parameter values (U_characteristic, L) to 0.0
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,2,0.0_DP,ERR,ERROR,*999)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                        !Default the analytic parameter values to 0
                        NUMBER_OF_ANALYTIC_COMPONENTS = 10
                        DO compIdx = 1,NUMBER_OF_ANALYTIC_COMPONENTS
                          CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,compIdx,0.0_DP,ERR,ERROR,*999)
                        END DO
                      CASE DEFAULT
                        LOCAL_ERROR="The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " is invalid for an analytical transient Navier-Stokes equation."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                      SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                        !Default the analytic parameter period values to 0
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " is invalid for a 1D Navier-Stokes equation."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    CASE DEFAULT
                      LOCAL_ERROR="The third equations set specification of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
                        & " is invalid for an analytical Navier-Stokes equation set."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  END IF
                END IF
              ELSE
                CALL FlagError("Equations set analytic is not associated.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for an analytic Navier-Stokes problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation set."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d 
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=1
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS1=2! viscosity, density
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF (EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, & 
                    & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, & 
                    & ERR,ERROR,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, & 
                    & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & "Materials",ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
              END IF
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF (EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  !Set the default values for the materials field
                  ! viscosity,density=1
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2,1.0_DP,ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", & 
                & ERR,ERROR))//" is invalid for Navier-Stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=2
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS1=2! U_var (constant)  : viscosity scale, density
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS2=2! V_var (gaussBased): viscosity, shear rate
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF (EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, & 
                    & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"MaterialsField",ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, & 
                    & ERR,ERROR,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, & 
                    & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    &[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
                  ! Set up U_VARIABLE (constants)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & "MaterialsConstants",ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  DO compIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & compIdx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  ! Set up V_VARIABLE (gauss-point based, CellML in/out parameters)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                     & "ConstitutiveValues",ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  DO compIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & compIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD, &
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
                  ! Check the U_VARIABLE
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,ERR,ERROR,*999)
                  ! Check the U_VARIABLE
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
              END IF
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  !Set the default values for the materials constants (viscosity scale, density)
                  DO compIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,compIdx,1.0_DP,ERR,ERROR,*999)
                  END DO
                  !Set the default values for the materials consitutive parameters (viscosity scale, density)
                  DO compIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,compIdx,1.0_DP,ERR,ERROR,*999)
                  END DO
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", & 
                & ERR,ERROR))//" is invalid for Navier-Stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! 1 variables for the 1D Navier-Stokes materials
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=2
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS1=9
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS2=8
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, & 
                    & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, & 
                    & ERR,ERROR,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, & 
                    & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                  ! 2 U,V materials field
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    &[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  ! Set up Navier-Stokes materials parameters
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  DO I=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1 !(MU,RHO,alpha,pressureExternal,LengthScale,TimeScale,MassScale)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & I,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  END DO  
                  DO I=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2 !(A0,E,H0)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  ! Set up coupling materials parameters
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE], &
                    & ERR,ERROR,*999)
                  ! Check N-S field variable
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                    & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
              END IF
              !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", & 
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", & 
                & ERR,ERROR))//" is invalid for Navier-Stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, & 
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            !\todo: Think about gravity
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !Do nothing
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
              !? Maybe set finished flag????
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & ERR,ERROR))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", & 
                & ERR,ERROR))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              &  " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              &  " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN              
                IF (EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                  CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                  CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations materials is not associated.",ERR,ERROR,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                !Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
              CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                !Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, & 
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)

            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN              
                IF (EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                  CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                  CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations materials is not associated.",ERR,ERROR,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EquationsMapping_ResidualVariableTypesSet(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR, & 
                  & ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                !Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE, &
                    & MATRIX_BLOCK_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, & 
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(EQUATIONS_MATRICES, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
               CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-Stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF (ASSOCIATED(EQUATIONS_MATERIALS)) THEN              
                IF (EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                  CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                  CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Equations materials is not associated.",ERR,ERROR,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                !Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(EQUATIONS_MATRICES, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(EQUATIONS_MATRICES, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-Stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a Navier-Stokes fluid subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    END IF

    EXITS("NAVIER_STOKES_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NAVIER_STOKES_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes problem pre solve.
  SUBROUTINE NAVIER_STOKES_PRE_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2,cellmlSolver
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    INTEGER(INTG) :: solver_matrix_idx,iteration
    REAL(DP) :: timeIncrement,currentTime
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(SOLVER2)
    
    ENTERS("NAVIER_STOKES_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          !Since we can have a fluid mechanics navier stokes equations set in a coupled problem setup we do not necessarily
          !have PROBLEM%SPECIFICATION(1)==FLUID_MECHANICS_CLASS
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF (ASSOCIATED(SOLVER_MAPPING)) THEN
                  ! TODO: Set up for multiple equations sets
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                  IF (ASSOCIATED(EQUATIONS_SET)) THEN
                    EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
                    IF (ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
                      !Update boundary conditions and any analytic values
                      CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
              !Update transient boundary conditions and any analytic values
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
              !Update transient boundary conditions
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              !Update transient boundary conditions
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
              !CALL NavierStokes_CalculateBoundaryFlux(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              nonlinearSolver=>SOLVER%DYNAMIC_SOLVER%NONLINEAR_SOLVER%NONLINEAR_SOLVER
              IF (ASSOCIATED(nonlinearSolver)) THEN
                !check for a linked CellML solver 
                cellmlSolver=>nonlinearSolver%NEWTON_SOLVER%CELLML_EVALUATOR_SOLVER
                IF (ASSOCIATED(cellmlSolver)) THEN
                  ! Calculate the CellML equations
                  CALL SOLVER_SOLVE(cellmlSolver,ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Nonlinear solver is not associated.",ERR,ERROR,*999)
              END IF

            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)

              SELECT CASE(SOLVER%SOLVE_TYPE)
              ! This switch takes advantage of the uniqueness of the solver types to do pre-solve operations
              ! for each of solvers in the various possible 1D subloops

              ! --- C h a r a c t e r i s t i c   S o l v e r ---
              CASE(SOLVER_NONLINEAR_TYPE)
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,currentTime,timeIncrement,ERR,ERROR,*999)
                iteration = CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
                EQUATIONS_SET=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                ! Characteristic solver effectively solves for the mass/momentum conserving fluxes at the
                ! *NEXT* timestep by extrapolating current field values and then solving a system of nonlinear
                ! equations: cons mass, continuity of pressure, and the characteristics.
                NULLIFY(fieldVariable)
                CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_INPUT_DATA1_SET_TYPE)%PTR)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_INPUT_DATA1_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_INPUT_DATA2_SET_TYPE,ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                 & FIELD_INPUT_DATA1_SET_TYPE,1.0_DP,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_RESIDUAL_SET_TYPE, &
                 & FIELD_INPUT_DATA2_SET_TYPE,1.0_DP,ERR,ERROR,*999)

                IF(iteration == 1) THEN
                  NULLIFY(fieldVariable)
                  CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                  END IF
                  ! Extrapolate new W from Q,A if this is the first timestep (otherwise will be calculated based on Navier-Stokes
                  ! values)
                  CALL Characteristic_Extrapolate(SOLVER,timeIncrement,ERR,ERROR,*999)
                END IF

              ! --- 1 D   N a v i e r - S t o k e s   S o l v e r ---
              CASE(SOLVER_DYNAMIC_TYPE)
                IF(SOLVER%global_number==2) THEN
                  ! update solver matrix
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                      SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
                      IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                        DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                          SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                          IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                            SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                          ELSE
                            CALL FlagError("Solver Matrix is not associated.",ERR,ERROR,*999)
                          END IF
                        END DO
                      ELSE
                        CALL FlagError("Solver Matrices is not associated.",ERR,ERROR,*999)
                      END IF
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF(ASSOCIATED(dependentField)) THEN
                          CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
                           & FIELD_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE, &
                           & FIELD_RESIDUAL_SET_TYPE,1.0_DP,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Dependent field is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  ! --- A d v e c t i o n   S o l v e r ---
                  CALL ADVECTION_PRE_SOLVE(SOLVER,ERR,ERROR,*999)
                END IF
                ! Update boundary conditions
                CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)

              ! --- C e l l M L    S o l v e r ---
              CASE(SOLVER_DAE_TYPE)
                ! DAE solver-set time
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,currentTime,timeIncrement,ERR,ERROR,*999)
                CALL SOLVER_DAE_TIMES_SET(SOLVER,currentTime,currentTime + timeIncrement,ERR,ERROR,*999)
                CALL SOLVER_DAE_TIME_STEP_SET(SOLVER,timeIncrement/1000.0_DP,ERR,ERROR,*999)

              ! --- S T R E E    S o l v e r ---
              CASE(SOLVER_LINEAR_TYPE)
                CALL Stree_PRE_SOLVE(SOLVER,ERR,ERROR,*999)

              CASE DEFAULT
                LOCAL_ERROR="The solve type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a 1D Navier-Stokes problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT   

            CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
              !do nothing ???
              !First update mesh and calculates boundary velocity values
              CALL NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH(SOLVER,ERR,ERROR,*999)
              !Then apply both normal and moving mesh boundary conditions
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              IF (SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement pre solve... ",ERR,ERROR,*999)
                !Update boundary conditions for mesh-movement
                CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,ERR,ERROR,*999)
                IF (ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                  SOLVER2%DYNAMIC_SOLVER%ALE=.FALSE.
                ELSE  
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",ERR,ERROR,*999)
                END IF
                !Update material properties for Laplace mesh movement
                CALL NavierStokes_PreSolveALEUpdateParameters(SOLVER,ERR,ERROR,*999)
                !Pre solve for the linear solver
              ELSE IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Navier-Stokes pre solve... ",ERR,ERROR,*999)
                IF (SOLVER%DYNAMIC_SOLVER%ALE) THEN
                  !First update mesh and calculates boundary velocity values
                  CALL NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH(SOLVER,ERR,ERROR,*999)
                  !Then apply both normal and moving mesh boundary conditions
                  CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                ELSE  
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
                END IF
              ELSE  
                CALL FlagError("Solver type is not associated for ALE problem.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
                !Pre solve for the linear solver
                IF (SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement pre solve... ",ERR,ERROR,*999)
                  !TODO if first time step smooth imported mesh with respect to absolute nodal position?

                  !Update boundary conditions for mesh-movement
                  CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER2,ERR,ERROR,*999)
                  IF (ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                    SOLVER2%DYNAMIC_SOLVER%ALE=.FALSE.
                  ELSE  
                    CALL FlagError("Dynamic solver is not associated for ALE problem.",ERR,ERROR,*999)
                  END IF
                  !Update material properties for Laplace mesh movement
                  CALL NavierStokes_PreSolveALEUpdateParameters(SOLVER,ERR,ERROR,*999)
                  !Pre solve for the dynamic solver which deals with the coupled FiniteElasticity-NavierStokes problem
                ELSE IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Navier-Stokes pre solve... ",ERR,ERROR,*999)
                  IF (SOLVER%DYNAMIC_SOLVER%ALE) THEN
                    !Apply both normal and moving mesh boundary conditions
                    CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                  ELSE  
                    CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
                  END IF
                ELSE  
                  CALL FlagError("Solver type is not associated for ALE problem.",ERR,ERROR,*999)
                END IF
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
                CALL FlagError(LOCAL_ERROR,Err,Error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
                & " is not valid for NAVIER_STOKES_PRE_SOLVE of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,Err,Error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(1),"*",ERR,ERROR))// &
              & " is not valid for Navier-Stokes fluid types."
            CALL FlagError(LOCAL_ERROR,Err,Error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Solvers are not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    END IF

    EXITS("NAVIER_STOKES_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_PRE_SOLVE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NAVIER_STOKES_PRE_SOLVE

! 
!================================================================================================================================
!

  !>Sets/changes the problem subtype for a Navier-Stokes fluid type.
  SUBROUTINE NavierStokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("NavierStokes_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
          !All ok
        CASE(PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE)
          CALL FlagError("Not implemented yet.",err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Navier-Stokes fluid mechanics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_NAVIER_STOKES_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Navier-Stokes problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_ProblemSpecificationSet")
    RETURN
999 ERRORS("NavierStokes_ProblemSpecificationSet",err,error)
    EXITS("NavierStokes_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE NavierStokes_ProblemSpecificationSet

! 
!================================================================================================================================
!
  
  !>Sets up the Navier-Stokes problem.
  SUBROUTINE NAVIER_STOKES_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: iterativeWhileLoop,iterativeWhileLoop2,simpleLoop
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,MESH_SOLVER_EQUATIONS,BIF_SOLVER_EQUATIONS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER, MESH_SOLVER,BIF_SOLVER,cellmlSolver
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NAVIER_STOKES_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(BIF_SOLVER)
    NULLIFY(BIF_SOLVER_EQUATIONS)
    NULLIFY(cellmlSolver)
    NULLIFY(CELLML_EQUATIONS)
    NULLIFY(CONTROL_LOOP)
    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(MESH_SOLVER)
    NULLIFY(MESH_SOLVER_EQUATIONS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
        !All steady state cases of Navier-Stokes
        CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing???
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a simple control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
                  !Set the solver to be a nonlinear solver
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                  !Set solver defaults
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !Transient cases and moving mesh
        CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
           & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing???
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a transient Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a time control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a transient Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
                  !Set the solver to be a first order dynamic solver 
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                  !Set solver defaults
                  CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                  !setup CellML evaluator
                  IF(PROBLEM%specification(3)==PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE) THEN
                    !Create the CellML evaluator solver
                    CALL SOLVER_NEWTON_CELLML_EVALUATOR_CREATE(SOLVER,cellmlSolver,ERR,ERROR,*999)
                    !Link the CellML evaluator solver to the solver
                    CALL SOLVER_LINKED_SOLVER_ADD(SOLVER,cellmlSolver,SOLVER_CELLML_EVALUATOR_TYPE,ERR,ERROR,*999)
                  END IF
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a transient Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                  & ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Get the control loop
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                !Get the solver
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                !Get the CellML evaluator solver
                CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,cellmlSolver,ERR,ERROR,*999)
                !Create the CellML equations
                CALL CELLML_EQUATIONS_CREATE_START(cellmlSolver,CELLML_EQUATIONS, &
                  & ERR,ERROR,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Get the control loop
                CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                !Get the solver
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                !Get the CellML evaluator solver
                CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,cellmlSolver,ERR,ERROR,*999)
                !Get the CellML equations for the CellML evaluator solver
                CALL SOLVER_CELLML_EQUATIONS_GET(cellmlSolver,CELLML_EQUATIONS,ERR,ERROR,*999)
                !Finish the CellML equations creation
                CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a CellML setup for a  transient Navier-Stokes equation."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a transient Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,     & !1D Navier-Stokes
           & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE,     & !with coupled 0D boundaries
           & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, & !with coupled advection
           & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, & !with coupled 0D boundaries and advection
           & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE,       & !with stree
           & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)     !with stree and advection

          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
          CASE(PROBLEM_SETUP_INITIAL_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              !Do nothing
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for Coupled1dDaeNavierStokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_SETUP_CONTROL_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              NULLIFY(CONTROL_LOOP_ROOT)
              !Time Loop
              CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
              NULLIFY(iterativeWhileLoop)
              IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                ! The 1D-0D boundary value iterative coupling loop
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",ERR,ERROR,*999)
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple CellML solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(simpleLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop2,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                 & ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                ! The 1D-0D boundary value iterative coupling loop
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",ERR,ERROR,*999)
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple CellML solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(simpleLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop2,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                 & ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple Advection solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(simpleLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"Advection",ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3) == PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E3_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3) == PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E6_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple Advection solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(simpleLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"Advection",ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3) == PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                ! The 1D-0D boundary value iterative coupling loop
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",ERR,ERROR,*999)
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple CellML solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(simpleLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop2,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                 & ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3) == PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                ! The 1D-0D boundary value iterative coupling loop
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",ERR,ERROR,*999)
                CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple CellML solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(simpleLoop,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,ERR,ERROR,*999)
                CALL CONTROL_LOOP_OUTPUT_TYPE_SET(iterativeWhileLoop2,CONTROL_LOOP_NO_OUTPUT,ERR,ERROR,*999)
                CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
                CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                 & ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                ! The simple Advection solver loop
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)                
                CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
                CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"Advection",ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a 1d transient Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          !Create the solvers
          CASE(PROBLEM_SETUP_SOLVERS_TYPE)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              SELECT CASE(PROBLEM%specification(3))
              CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(iterativeWhileLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,2,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !!!-- N A V I E R   S T O K E S --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(iterativeWhileLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(iterativeWhileLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,2,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !!!-- N A V I E R   S T O K E S --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                ! Simple loop 1 contains the Advection solver
                ! (this subloop holds 1 solver)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Advection Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)

                ! Simple loop 1 contains the 0D/CellML DAE solver
                ! (this subloop holds 1 solver)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- D A E --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(solvers,1,solver,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(solver,SOLVER_DAE_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"DAE Solver",ERR,ERROR,*999)

                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,2,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !!!-- N A V I E R   S T O K E S --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)

                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Advection Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                ! Simple loop 1 contains the 0D/CellML DAE solver
                ! (this subloop holds 1 solver)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- D A E --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(solvers,1,solver,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(solver,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Linear Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,2,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !!!-- N A V I E R   S T O K E S --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Advection Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)

                ! Simple loop 1 contains the 0D/CellML DAE solver
                ! (this subloop holds 1 solver)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- D A E --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(solvers,1,solver,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(solver,SOLVER_DAE_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"DAE Solver",ERR,ERROR,*999)

                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,2,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !!!-- N A V I E R   S T O K E S --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)

                ! Simple loop 1 contains the 0D/CellML DAE solver
                ! (this subloop holds 1 solver)
                NULLIFY(simpleLoop)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
                !!!-- D A E --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(solvers,1,solver,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(solver,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Linear Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)

                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(solvers,2,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                !!!-- N A V I E R   S T O K E S --!!!
                NULLIFY(SOLVER)
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
                  & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              IF(PROBLEM%specification(3)==PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3)==PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3)==PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3)==PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3)==PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              ELSE IF(PROBLEM%specification(3)==PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                NULLIFY(SOLVERS)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a 1d transient Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          !Create the solver equations
          CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              NULLIFY(CONTROL_LOOP)
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              SELECT CASE(PROBLEM%specification(3))
              CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(simpleLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(simpleLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- D A E --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(simpleLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- D A E --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
                  & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              NULLIFY(CONTROL_LOOP)
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              SELECT CASE(PROBLEM%specification(3))
              CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(simpleLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(simpleLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- D A E --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(simpleLoop)
                ! Iterative loop couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                !!!-- A D V E C T I O N --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
                NULLIFY(iterativeWhileLoop)
                ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
                ! (this subloop holds 2 subloops)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- D A E --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVERS)
                NULLIFY(SOLVER_EQUATIONS)
                NULLIFY(iterativeWhileLoop2)
                ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
                ! (this subloop holds 2 solvers)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,ERR,ERROR,*999)
                !!!-- C H A R A C T E R I S T I C --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
                NULLIFY(SOLVER)
                NULLIFY(SOLVER_EQUATIONS)
                !!!-- N A V I E R   S T O K E S --!!!
                CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
                  & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          !Create the CELLML solver equations
          CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP              
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
                & PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              ELSE
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
              END IF
              NULLIFY(SOLVER)
              NULLIFY(cellMLSolver)
              NULLIFY(CELLML_EQUATIONS)
              SELECT CASE(PROBLEM%specification(3))
              CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                 & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL CELLML_EQUATIONS_CREATE_START(solver,CELLML_EQUATIONS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
                  & " is not valid for cellML equations setup Navier-Stokes equation type of a fluid mechanics problem class."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
                & PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(iterativeWhileLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
                NULLIFY(simpleLoop)
                CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
                CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              ELSE
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
              END IF
              NULLIFY(SOLVER)
              SELECT CASE(PROBLEM%specification(3))
              CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                 & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_CELLML_EQUATIONS_GET(solver,CELLML_EQUATIONS,ERR,ERROR,*999)
                CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The third problem specification of "// &
                  & TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
                  & " is not valid for cellML equations setup Navier-Stokes fluid mechanics problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a CellML setup for a 1D Navier-Stokes equation."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a 1d transient Navier-Stokes fluid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
          !Quasi-static Navier-Stokes
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing???
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a quasistatic Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a time control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a quasistatic Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
              CASE(PROBLEM_SETUP_START_ACTION)
                !Start the solvers creation
                CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
                !Set the solver to be a nonlinear solver
                CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
                !Set solver defaults
                CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
              CASE(PROBLEM_SETUP_FINISH_ACTION)
                !Get the solvers
                CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                !Finish the solvers creation
                CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                  & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a quasistatic Navier-Stokes equation."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a quasistatic Navier-Stokes equation."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a quasistatic Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !Navier-Stokes ALE cases
        CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
            CASE(PROBLEM_SETUP_INITIAL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Do nothing????
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Do nothing????
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a ALE Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_CONTROL_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Set up a time control loop
                  CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Finish the control loops
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a ALE Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVERS_TYPE)
              !Get the control loop
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Start the solvers creation
                  CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_NUMBER_SET(SOLVERS,2,ERR,ERROR,*999)
                  !Set the first solver to be a linear solver for the Laplace mesh movement problem
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_TYPE_SET(MESH_SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
                  !Set solver defaults
                  CALL SOLVER_LIBRARY_TYPE_SET(MESH_SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                  !Set the solver to be a first order dynamic solver 
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
                  !Set solver defaults
                  CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
                  CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
                  CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the solvers
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  !Finish the solvers creation
                  CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a ALE Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
              SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
                CASE(PROBLEM_SETUP_START_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(MESH_SOLVER,MESH_SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                  !Create the solver equations
                  CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                  & ERR,ERROR,*999)
                  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                CASE(PROBLEM_SETUP_FINISH_ACTION)
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(MESH_SOLVER,MESH_SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(MESH_SOLVER_EQUATIONS,ERR,ERROR,*999)             
                  !Get the solver equations
                  CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Finish the solver equations creation
                  CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
                CASE DEFAULT
                  LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a Navier-Stokes fluid."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a ALE Navier-Stokes fluid."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The third problem specification of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%specification(3),"*",ERR,ERROR))// &
            & " is not valid for a Navier-Stokes fluid mechanics problem."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    END IF

    EXITS("NAVIER_STOKES_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NAVIER_STOKES_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Navier-Stokes equation finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS,INDEPENDENT_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: STIFFNESS_MATRIX,DAMPING_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,INDEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    INTEGER(INTG) :: ng,mh,mhs,mi,ms,nh,nhs,ni,ns,nhs_max,mhs_max,nhs_min,mhs_min,xv,out
    INTEGER(INTG) :: FIELD_VAR_TYPE,MESH_COMPONENT1,MESH_COMPONENT2,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: nodeIdx,derivIdx,versionIdx,compIdx,xiIdx,coordIdx
    INTEGER(INTG) :: numberOfVersions,nodeNumber,numberOfElementNodes,numberOfParameters,firstNode,lastNode
    REAL(DP) :: JGW,SUM,X(3),DXI_DX(3,3),DPHIMS_DXI(3),DPHINS_DXI(3),PHIMS,PHINS,momentum,mass,Qupwind,Aupwind,kp,k1,k2,k3,b1
    REAL(DP) :: U_VALUE(3),W_VALUE(3),U_DERIV(3,3),Q_VALUE,A_VALUE,Q_DERIV,A_DERIV,pressure,normalWave,lengthScale,timeScale
    REAL(DP) :: MU,RHO,Fr,A0,E,H,A0_DERIV,E_DERIV,H_DERIV,alpha,beta,G0,x1,y1,z1,x2,y2,z2,slope,muScale,massScale,Pext
    REAL(DP), POINTER :: dependentParameters(:),materialsParameters(:),materialsParameters1(:)
    LOGICAL :: UPDATE_STIFFNESS_MATRIX,UPDATE_DAMPING_MATRIX,UPDATE_RHS_VECTOR,UPDATE_NONLINEAR_RESIDUAL
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NavierStokes_FiniteElementResidualEvaluate",ERR,ERROR,*999)

    UPDATE_STIFFNESS_MATRIX=.FALSE.
    UPDATE_DAMPING_MATRIX=.FALSE.
    UPDATE_RHS_VECTOR=.FALSE.
    UPDATE_NONLINEAR_RESIDUAL=.FALSE.
    X=0.0_DP
    out=0

    NULLIFY(DEPENDENT_BASIS,GEOMETRIC_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(LINEAR_MAPPING)
    NULLIFY(NONLINEAR_MAPPING)
    NULLIFY(DYNAMIC_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(LINEAR_MATRICES)
    NULLIFY(NONLINEAR_MATRICES)
    NULLIFY(DYNAMIC_MATRICES)
    NULLIFY(RHS_VECTOR)
    NULLIFY(STIFFNESS_MATRIX,DAMPING_MATRIX)
    NULLIFY(DEPENDENT_FIELD,INDEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD)
    NULLIFY(dependentParameters,materialsParameters,materialsParameters1)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(QUADRATURE_SCHEME)
    NULLIFY(QUADRATURE_SCHEME1,QUADRATURE_SCHEME2)
    NULLIFY(DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          !Set general and specific pointers
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            STIFFNESS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            STIFFNESS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
          CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE)
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
            DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(DAMPING_MATRIX)) UPDATE_DAMPING_MATRIX=DAMPING_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, & 
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
            DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(DAMPING_MATRIX)) UPDATE_DAMPING_MATRIX=DAMPING_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
             & MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
            DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(DAMPING_MATRIX)) UPDATE_DAMPING_MATRIX=DAMPING_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
          CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            &  EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            DECOMPOSITION => DEPENDENT_FIELD%DECOMPOSITION
            MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
            DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(DAMPING_MATRIX)) UPDATE_DAMPING_MATRIX=DAMPING_MATRIX%UPDATE_MATRIX
            IF(ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF(ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
          CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
            INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
            INDEPENDENT_BASIS=>INDEPENDENT_FIELD%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)% & 
              & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
            DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF (ASSOCIATED(STIFFNESS_MATRIX)) UPDATE_STIFFNESS_MATRIX=STIFFNESS_MATRIX%UPDATE_MATRIX
            IF (ASSOCIATED(DAMPING_MATRIX)) UPDATE_DAMPING_MATRIX=DAMPING_MATRIX%UPDATE_MATRIX
            IF (ASSOCIATED(RHS_VECTOR)) UPDATE_RHS_VECTOR=RHS_VECTOR%UPDATE_VECTOR
            IF (ASSOCIATED(NONLINEAR_MATRICES)) UPDATE_NONLINEAR_RESIDUAL=NONLINEAR_MATRICES%UPDATE_RESIDUAL
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_MESH_VELOCITY_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & INDEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)

          !Calculate slope of vessels for gravitational force
          ELEMENTS_TOPOLOGY=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%components(1)%domain%topology%elements
          numberOfElementNodes=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_NODES
          firstNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(1)
          lastNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(numberOfElementNodes)
          CALL Field_ParameterSetGetLocalNode(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,firstNode, & 
            & 1,x1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,firstNode, &
            & 2,y1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,firstNode, & 
            & 3,z1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,lastNode, & 
            & 1,x2,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,lastNode, &
            & 2,y2,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,lastNode, &
            & 3,z2,err,error,*999)
          slope=(z2-z1)/(((x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0)**0.5)

          !Loop over Gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              W_VALUE(1)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              W_VALUE(2)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              IF (FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                W_VALUE(3)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
              END IF 
            ELSE
              W_VALUE=0.0_DP
            END IF

            ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
              ! Note the constant from the U_VARIABLE is a scale factor
              muScale = MU
              ! Get the gauss point based value returned from the CellML solver
              CALL Field_ParameterSetGetLocalGaussPoint(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ng,ELEMENT_NUMBER,1,MU,ERR,ERROR,*999)
              MU=MU*muScale
            END IF

            !Start with matrix calculations
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
              !Loop over field components
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                  & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)

                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  IF (UPDATE_STIFFNESS_MATRIX.OR.UPDATE_DAMPING_MATRIX) THEN
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                        &(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                      ! JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME2%&
                      ! &GAUSS_WEIGHTS(ng)                        
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        !Calculate some general values
                        DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                          DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DXI_DX(mi,ni)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                              & DXI_DX(mi,ni)
                          ENDDO
                          DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        ENDDO !ni
                        PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        !Laplace only matrix
                        IF (UPDATE_STIFFNESS_MATRIX) THEN
                          !LAPLACE TYPE 
                          IF (nh==mh) THEN 
                            SUM=0.0_DP
                            !Calculate SUM 
                            DO xv=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                  SUM=SUM+MU*DPHINS_DXI(ni)*DXI_DX(ni,xv)*DPHIMS_DXI(mi)*DXI_DX(mi,xv)
                                END DO !ni
                              END DO !mi
                            END DO !x 
                            !Calculate MATRIX  
                            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                          ENDIF
                        ENDIF
                        !General matrix
                        IF (UPDATE_STIFFNESS_MATRIX) THEN
                          !GRADIENT TRANSPOSE TYPE
                          IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE) THEN 
                            IF(nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN 
                              SUM=0.0_DP
                              !Calculate SUM 
                              DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                  !note mh/nh derivative in DXI_DX 
                                  SUM=SUM+MU*DPHINS_DXI(mi)*DXI_DX(mi,mh)*DPHIMS_DXI(ni)*DXI_DX(ni,nh)
                                END DO !ni
                              END DO !mi
                              !Calculate MATRIX
                              STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) &
                               & +SUM*JGW
                           END IF
                          END IF
                        END IF
                        !Contribution through ALE
                        IF(UPDATE_STIFFNESS_MATRIX) THEN
                          !GRADIENT TRANSPOSE TYPE
                          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. & 
                            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN 
                            IF(nh==mh) THEN 
                              SUM=0.0_DP
                              !Calculate SUM 
                              DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                  SUM=SUM-RHO*W_VALUE(mi)*DPHINS_DXI(ni)*DXI_DX(ni,mi)*PHIMS
                                END DO !ni
                              END DO !mi
                              !Calculate MATRIX
                              STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                               & SUM*JGW
                            END IF
                          END IF
                        END IF
                        !Pressure contribution (B transpose)
                        IF (UPDATE_STIFFNESS_MATRIX) THEN
                          !LAPLACE TYPE 
                          IF (nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN 
                            SUM=0.0_DP
                            !Calculate SUM 
                            DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              SUM=SUM-PHINS*DPHIMS_DXI(ni)*DXI_DX(ni,mh)
                            END DO !ni
                            !Calculate MATRIX
                            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                          ENDIF
                        ENDIF
                        !Damping matrix
                        IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE .OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                          IF(UPDATE_DAMPING_MATRIX) THEN
                            IF(nh==mh) THEN 
                              SUM=0.0_DP 
                              !Calculate SUM 
                              SUM=PHIMS*PHINS*RHO
                              !Calculate MATRIX
                              DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                            END IF
                          END IF
                        END IF
                      END DO !ns    
                    END DO !nh
                  END IF
                END DO !ms
              END DO !mh
              !Analytic RHS vector
              IF (RHS_VECTOR%FIRST_ASSEMBLY) THEN
                IF (UPDATE_RHS_VECTOR) THEN
                  IF (ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                    IF (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5) THEN
                      mhs=0
                      DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                        MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                        JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                          & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)

                        DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                          mhs=mhs+1
                          PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          !note mh value derivative 
                          SUM=0.0_DP 
                          X(1)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1)
                          X(2)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,1)
                          IF (DEPENDENT_BASIS1%NUMBER_OF_XI==3) THEN
                            X(3)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,1)
                          ENDIF
                          IF (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1) THEN
                            IF (mh==1) THEN 
                              !Calculate SUM 
                              SUM=0.0_DP                         
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-2.0_DP/3.0_DP*(X(1)**3*RHO+3.0_DP*MU*10.0_DP**2- &
                                & 3.0_DP*RHO*X(2)**2*X(1))/(10.0_DP**4))
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2) &
                            & THEN
                            IF (mh==1) THEN 
                              !Calculate SUM 
                              SUM=0.0_DP                               
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-4.0_DP*MU/10.0_DP/10.0_DP*EXP((X(1)-X(2))/10.0_DP))
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3) & 
                           & THEN
                            IF(mh==1) THEN 
                              !Calculate SUM 
                              SUM=0.0_DP         
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(16.0_DP*MU*PI**2/10.0_DP**2*COS(2.0_DP*PI*X(2)/10.0_DP)* &
                                & COS(2.0_DP*PI*X(1)/10.0_DP)- &
                                & 2.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)*RHO*PI/10.0_DP)
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4) & 
                           & THEN
                            IF(mh==1) THEN 
                              !Calculate SUM 
                              SUM=PHIMS*(2.0_DP*SIN(X(1))*COS(X(2)))*MU
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-2.0_DP*COS(X(1))*SIN(X(2)))*MU
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5) & 
                           & THEN
                            !do nothing
                          ELSE IF (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== & 
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                            IF (mh==1) THEN 
                              !Calculate SUM 
                              SUM=0.0_DP       
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-2.0_DP/3.0_DP*(RHO*X(1)**3+6.0_DP*RHO*X(1)*X(3)*X(2)+ &
                                & 6.0_DP*MU*10.0_DP**2- & 
                                & 3.0_DP*RHO*X(2)**2*X(1)-3.0_DP*RHO*X(3)*X(1)**2-3.0_DP*RHO*X(3)*X(2)**2)/ &
                                  & (10.0_DP**4))
                            ELSE IF(mh==3) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-2.0_DP/3.0_DP*(6.0_DP*RHO*X(1)*X(3)*X(2)+RHO*X(1)**3+ &
                                & 6.0_DP*MU*10.0_DP**2- & 
                                & 3.0_DP*RHO*X(1)*X(3)**2-3.0_DP*RHO*X(2)*X(1)**2-3.0_DP*RHO*X(2)*X(3)**2)/ & 
                                & (10.0_DP**4))
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== & 
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2) THEN
                            IF (mh==1) THEN 
                              !Calculate SUM 
                              SUM=0.0_DP         
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*((-4.0_DP*MU*EXP((X(1)-X(2))/10.0_DP)-2.0_DP*MU*EXP((X(2)-X(3))/10.0_DP)+ & 
                                & RHO*EXP((X(3)-X(2))/10.0_DP)*10.0_DP)/10.0_DP**2)
                            ELSE IF (mh==3) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-(4.0_DP*MU*EXP((X(3)-X(1))/10.0_DP)+2.0_DP*MU*EXP((X(2)-X(3))/10.0_DP)+ & 
                                & RHO*EXP((X(3)-X(2))/10.0_DP)*10.0_DP)/10.0_DP** 2)
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== & 
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3) THEN
                            IF (mh==1) THEN 
                              !Calculate SUM 
                              SUM=0.0_DP         
                            ELSE IF (mh==2) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(2.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*(18.0_DP*COS(2.0_DP*PI*X(1)/10.0_DP)* &
                                & MU*PI*SIN(2.0_DP*PI*X(3)/10.0_DP)-3.0_DP*RHO*COS(2.0_DP*PI*X(1)/10.0_DP)**2* &
                                & SIN(2.0_DP*PI*X(2)/10.0_DP)*10.0_DP-2.0_DP*RHO*SIN(2.0_DP*PI*X(2)/10.0_DP)*10.0_DP+ & 
                                & 2.0_DP*RHO*SIN(2.0_DP*PI*X(2)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*X(3)/10.0_DP)**2)*PI/ &
                                & 10.0_DP**2)
                            ELSE IF (mh==3) THEN
                              !Calculate SUM 
                              SUM=PHIMS*(-2.0_DP*PI*COS(2.0_DP*PI*X(3)/10.0_DP)*RHO*SIN(2.0_DP*PI*X(3)/10.0_DP)* & 
                                & (-1.0_DP+COS(2.0_DP*PI*X(2)/10.0_DP)**2)/10.0_DP)
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4) THEN
                            IF(mh==1) THEN 
                              !Calculate SUM 
                              !SUM=PHIMS*(2.0_DP*SIN(X(1))*COS(X(2)))*MU_PARAM
                            ELSE IF(mh==2) THEN
                              !Calculate SUM 
                              !SUM=PHIMS*(-2.0_DP*COS(X(1))*SIN(X(2)))*MU_PARAM
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5) THEN
                            !do nothing
                          END IF
                          !Calculate RH VECTOR
                           RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+SUM*JGW
                        END DO !ms
                      END DO !mh
                    ELSE
                      RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
                    END IF                 
                  END IF
                END IF                                                                     
              END IF

              !Calculate nonlinear vector
              IF(UPDATE_NONLINEAR_RESIDUAL) THEN
                ! Get interpolated velocity and velocity gradient values for nonlinear term
                U_VALUE(1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                U_VALUE(2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                U_DERIV(1,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S1)
                U_DERIV(1,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S2)
                U_DERIV(2,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,PART_DERIV_S1)
                U_DERIV(2,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,PART_DERIV_S2)
                IF (FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                  U_VALUE(3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                  U_DERIV(3,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,PART_DERIV_S1)
                  U_DERIV(3,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,PART_DERIV_S2)
                  U_DERIV(3,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,PART_DERIV_S3)
                  U_DERIV(1,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S3)
                  U_DERIV(2,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,PART_DERIV_S3) 
                ELSE
                  U_VALUE(3)=0.0_DP
                  U_DERIV(3,1)=0.0_DP
                  U_DERIV(3,2)=0.0_DP
                  U_DERIV(3,3)=0.0_DP
                  U_DERIV(1,3)=0.0_DP
                  U_DERIV(2,3)=0.0_DP
                ENDIF
                !Here W_VALUES must be ZERO if ALE part of linear matrix
                W_VALUE=0.0_DP
                mhs=0
                DO mh=1,(FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1)
                  MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                  JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                    & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                  DXI_DX=0.0_DP

                  DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                    DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      DXI_DX(mi,ni)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%DXI_DX(mi,ni)
                    ENDDO
                  ENDDO

                  DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    !note mh value derivative 
                    SUM=0.0_DP
                    ! Convective form
                    DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      SUM=SUM+RHO*(PHIMS)*( & 
                        & (U_VALUE(1))*(U_DERIV(mh,ni)*DXI_DX(ni,1))+ &
                        & (U_VALUE(2))*(U_DERIV(mh,ni)*DXI_DX(ni,2))+ &
                        & (U_VALUE(3))*(U_DERIV(mh,ni)*DXI_DX(ni,3)))
                    END DO !ni

                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)+SUM*JGW

                  END DO !ms
                END DO !mh
              END IF
            END IF

            !------------------------------------------------------------------
            ! R e s i d u a l - b a s e d    S t a b i l i s a t i o n
            !------------------------------------------------------------------
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
              CALL NavierStokes_ResidualBasedStabilisation(EQUATIONS_SET,ELEMENT_NUMBER,ng, &
               & MU,RHO,.FALSE.,ERR,ERROR,*999)
            END IF

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!                                        !!!!!
            !!!!!         1 D  T R A N S I E N T         !!!!!
            !!!!!                                        !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Start with matrix calculations
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              Q_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              A_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              Q_DERIV=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,FIRST_PART_DERIV)
              A_DERIV=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,FIRST_PART_DERIV)
              A0=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
              E =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
              H =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
              kp=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(4,NO_PART_DERIV)
              k1=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(5,NO_PART_DERIV)
              k2=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(6,NO_PART_DERIV)
              k3=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(7,NO_PART_DERIV)
              b1=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(8,NO_PART_DERIV)
              A0_DERIV=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(1,FIRST_PART_DERIV)
              E_DERIV =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(2,FIRST_PART_DERIV)
              H_DERIV =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(3,FIRST_PART_DERIV)
              beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

              !Get material constants
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &                
                & MU,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2, &
                & RHO,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3, &
                & alpha,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4, &
                & Pext,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,5, &
                & lengthScale,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,6, &
                & timeScale,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,7, &
                & massScale,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,8, &
                & G0,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,9, &
                & Fr,err,error,*999)

              !If A goes negative during nonlinear iteration, give ZERO_TOLERANCE value to avoid segfault
              IF (A_VALUE<A0*0.001_DP) A_VALUE=A0*0.001_DP

              mhs=0
              !Loop Over Element Rows
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                  & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(mh)%DOMAIN%TOPOLOGY%ELEMENTS
                DXI_DX=0.0_DP
                !Calculate dxi_dx in 3D
                DO xiIdx=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                  DO coordIdx=1,EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE) &
                    & %PTR%NUMBER_OF_X_DIMENSIONS
                    DXI_DX(1,1)=DXI_DX(1,1)+(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)% &
                      & PTR%DXI_DX(xiIdx,coordIdx))**2
                  ENDDO !coordIdx
                ENDDO !xiIdx
                DXI_DX(1,1)=SQRT(DXI_DX(1,1))
                !Loop Over Element rows
                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  DPHIMS_DXI(1)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,FIRST_PART_DERIV,ng)
                  mhs=mhs+1
                  nhs=0
                  IF (UPDATE_STIFFNESS_MATRIX .OR. UPDATE_DAMPING_MATRIX) THEN
                    !Loop Over Element Columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        DPHINS_DXI(1)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,FIRST_PART_DERIV,ng)
                        nhs=nhs+1

                        !!!-- D A M P I N G  M A T R I X --!!!
                        IF (UPDATE_DAMPING_MATRIX) THEN
                          !Momentum Equation, dQ/dt
                          IF (mh==1 .AND. nh==1) THEN
                            SUM=PHINS*PHIMS
                            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)= &
                              & DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                          ENDIF
                          !Mass Equation, dA/dt
                          IF (mh==2 .AND. nh==2) THEN
                            SUM=PHINS*PHIMS
                            DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)= &
                              & DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                          ENDIF
                        ENDIF

                        !!!-- S T I F F N E S S  M A T R I X --!!!
                        IF (UPDATE_STIFFNESS_MATRIX) THEN
                          !Momentum Equation, gravitational force
                          IF (mh==1 .AND. nh==2) THEN
                            SUM=-PHINS*PHIMS*beta/RHO* &
                               & (k1*A0_DERIV/A0+      &    !dA0/dx (linear part)
                               &  k2*E_DERIV/E+        &    !dE/dx  (linear part)
                               &  k3*H_DERIV/H)*DXI_DX(1,1) !dH/dx  (linear part)
                            !SUM=PHINS*PHIMS*G0*slope
                            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)= &
                              & STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                          ENDIF
                          !Mass Equation, dQ/dX, flow derivative
                          IF (mh==2 .AND. nh==1) THEN 
                            SUM=DPHINS_DXI(1)*DXI_DX(1,1)*PHIMS
                            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)= &
                              & STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*JGW
                          END IF
                        END IF

                      END DO !ns    
                    END DO !nh
                  END IF
                  
                  !!!-- N O N L I N E A R  V E C T O R --!!!
                  IF (UPDATE_NONLINEAR_RESIDUAL) THEN
                    !Momentum Equation
                    IF (mh==1) THEN
                      SUM=((2.0_DP*alpha*(Q_VALUE/A_VALUE)*Q_DERIV-               & !Convective
                        & alpha*((Q_VALUE/A_VALUE)**2)*A_DERIV+                   & !Convective
                        & beta/RHO*(A_DERIV*b1*(A_VALUE/A0)**b1+                  & !dA/dx
                        & k1*A0_DERIV*A_VALUE/A0*(1.0_DP-b1/k1)*(A_VALUE/A0)**b1+ & !dA0/dx (nonlinear part)
                        & k2*E_DERIV*A_VALUE/E*(A_VALUE/A0)**b1+                  & !dE/dx  (nonlinear part)
                        & k3*H_DERIV*A_VALUE/H*(A_VALUE/A0)**b1))*                & !dH/dx  (nonlinear part)
                        & DXI_DX(1,1)+Fr*Q_VALUE/A_VALUE)*PHIMS                     !Viscosity
                      NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)= &
                        & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)+SUM*JGW
                    END IF
                  END IF
                END DO !ms
              END DO !mh
            END IF
          END DO !ng

          IF (EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            IF(UPDATE_NONLINEAR_RESIDUAL) THEN
              ELEMENTS_TOPOLOGY=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%components(1)%domain%topology%elements
              numberOfElementNodes=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_NODES
              numberOfParameters=ELEMENTS_TOPOLOGY%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS  
              firstNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(1)
              lastNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(numberOfElementNodes)
              derivIdx=1
              compIdx=1
              !Loop over nodes in a element
              DO nodeIdx=1,numberOfElementNodes
                nodeNumber=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nodeIdx)
                numberOfVersions=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%NODES%NODES(nodeNumber)% &
                  & DERIVATIVES(derivIdx)%numberOfVersions
                versionIdx=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)% &
                  & elementVersions(derivIdx,nodeIdx)

                !Get current Area values
                CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,2,A_VALUE,err,error,*999)
                !If A goes negative during nonlinear iteration, set to positive value to avoid segfault
                IF (A_VALUE<A0*0.001_DP) A_VALUE=A0*0.001_DP

                !Get materials variables for node on this element
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,1,A0,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,2,E,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,3,H,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,4,kp,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,5,k1,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,6,k2,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,7,k3,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,8,b1,err,error,*999)
                beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

                !!!-- P R E S S U R E    C A L C U L A T I O N --!!!
                !Pressure equation in mmHg
                pressure=(Pext+beta*((A_VALUE/A0)**b1-1.0_DP))/(massScale/(lengthScale*timeScale**2))*0.0075_DP
                !Update the dependent field
                IF (ELEMENT_NUMBER<=DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & versionIdx,derivIdx,nodeNumber,1,pressure,err,error,*999)
                ENDIF

                !!!-- B R A N C H   F L U X   U P W I N D I N G --!!!
                !----------------------------------------------------
                ! In order to enforce conservation of mass and momentum across discontinuous
                !  branching topologies, flux is upwinded against the conservative branch values 
                !  established by the characteristic solver.

                !Find the branch node on this element
                IF (numberOfVersions>1) THEN
                  !Find the wave direction - incoming or outgoing
                  CALL Field_ParameterSetGetLocalNode(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, & 
                    & versionIdx,derivIdx,nodeNumber,compIdx,normalWave,err,error,*999)
                  !Get current Q & A values for node on this element
                  CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & versionIdx,derivIdx,nodeNumber,1,Q_VALUE,err,error,*999)         
                  CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & versionIdx,derivIdx,nodeNumber,2,A_VALUE,err,error,*999)
                  !Get upwind Q & A values based on the branch (characteristics) solver
                  CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
                    & versionIdx,derivIdx,nodeNumber,1,Qupwind,err,error,*999)         
                  CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
                    & versionIdx,derivIdx,nodeNumber,2,Aupwind,err,error,*999)

                  !Momentum Equation: F_upwind-F_current
                  momentum=alpha*(Qupwind**2/Aupwind-Q_VALUE**2/A_VALUE)+beta*A0/RHO* &
                    & b1/(b1+1.0_DP)*((Aupwind/A0)**(b1+1.0_DP)-(A_VALUE/A0)**(b1+1.0_DP))*normalWave
                  !Continuity Equation 
                  mass=(Qupwind-Q_VALUE)*normalWave

                  !Add momentum/mass contributions to first/last node accordingly
                  IF (nodeNumber==firstNode) THEN
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(1)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(1)+momentum
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(numberOfParameters+1)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(numberOfParameters+1)+mass
                  ELSE IF(nodeNumber==lastNode) THEN
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(numberOfParameters)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(numberOfParameters)+momentum
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(numberOfParameters*2)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(numberOfParameters*2)+mass
                  ENDIF
                ENDIF !find branch node
              ENDDO !loop over nodes
            ENDIF !update nonlinear vector
          ENDIF !equation set

          ! F a c e   I n t e g r a t i o n
          IF (RHS_VECTOR%UPDATE_VECTOR) THEN
            !If specified, also perform face integration for neumann boundary conditions
            IF (DEPENDENT_FIELD%DECOMPOSITION%CALCULATE_FACES) THEN
              CALL NavierStokes_FiniteElementFaceIntegrate(EQUATIONS_SET,ELEMENT_NUMBER,FIELD_VARIABLE,ERR,ERROR,*999)
            ENDIF
          ENDIF

          !!!--   A S S E M B L E   M A T R I C E S  &  V E C T O R S   --!!!
          mhs_min=mhs
          mhs_max=nhs
          nhs_min=mhs
          nhs_max=nhs
          IF (EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
            IF(STIFFNESS_MATRIX%FIRST_ASSEMBLY) THEN
              IF(UPDATE_STIFFNESS_MATRIX) THEN
                DO mhs=mhs_min+1,mhs_max
                  DO nhs=1,nhs_min
                    !Transpose pressure type entries for mass equation  
                    STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=-STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(nhs,mhs)
                  END DO
                END DO
              END IF
            END IF
          END IF
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Navier-Stokes equation type of a classical field equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokes_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementResidualEvaluate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NavierStokes_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices and RHS for a Navier-Stokes equation finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS,INDEPENDENT_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: STIFFNESS_MATRIX 
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,INDEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    INTEGER(INTG) :: ng,mh,mhs,mi,ms,nh,nhs,ni,ns,x,xiIdx,coordIdx
    INTEGER(INTG) :: nodeIdx,derivIdx,versionIdx,firstNode,lastNode,nodeNumber
    INTEGER(INTG) :: numberOfElementNodes,numberOfParameters,numberOfVersions,compIdx
    INTEGER(INTG) :: FIELD_VAR_TYPE,MESH_COMPONENT_NUMBER,MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: JGW,SUM,DXI_DX(3,3),DPHIMS_DXI(3),DPHINS_DXI(3),PHIMS,PHINS,G0
    REAL(DP) :: U_VALUE(3),W_VALUE(3),U_DERIV(3,3),Q_VALUE,Q_DERIV,A_VALUE,A_DERIV,alpha,beta,normalWave
    REAL(DP) :: MU,RHO,Fr,A0,A0_DERIV,E,E_DERIV,H,H_DERIV,mass,momentum1,momentum2,kp,k1,k2,k3,b1,muScale
    LOGICAL :: UPDATE_JACOBIAN_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NavierStokes_FiniteElementJacobianEvaluate",ERR,ERROR,*999)

    DXI_DX=0.0_DP
    UPDATE_JACOBIAN_MATRIX=.FALSE.

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      NULLIFY(EQUATIONS)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%specification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          !Set some general and case-specific pointers
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          SELECT CASE(EQUATIONS_SET%specification(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            STIFFNESS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            STIFFNESS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            !SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
            STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX=0.0_DP
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
          CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE)
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX=0.0_DP
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX=0.0_DP
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX=0.0_DP
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
          CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            &  EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            DECOMPOSITION => DEPENDENT_FIELD%DECOMPOSITION
            MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX=0.0_DP
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
          CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
            INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD
            INDEPENDENT_BASIS=>INDEPENDENT_FIELD%DECOMPOSITION%DOMAIN(INDEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)% & 
              & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX=0.0_DP
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
            FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
            IF(ASSOCIATED(JACOBIAN_MATRIX)) UPDATE_JACOBIAN_MATRIX=JACOBIAN_MATRIX%UPDATE_JACOBIAN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_MESH_VELOCITY_SET_TYPE,ELEMENT_NUMBER,EQUATIONS% & 
              & INTERPOLATION%INDEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%specification(3),"*",ERR,ERROR))// &
              & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
              & MU,err,error,*999)
            CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2, &
              & RHO,err,error,*999)
            !Loop over all Gauss points 
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
                & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                & MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
                  & INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  W_VALUE(1)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                  W_VALUE(2)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                  IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                    W_VALUE(3)=EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                  END IF 
              ELSE
                W_VALUE=0.0_DP
              END IF

              ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                ! Note the constant from the U_VARIABLE is a scale factor
                muScale = MU
                ! Get the gauss point based value returned from the CellML solver
                CALL Field_ParameterSetGetLocalGaussPoint(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & ng,ELEMENT_NUMBER,1,MU,ERR,ERROR,*999)
                MU=MU*muScale
              END IF

              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN

                U_VALUE(1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                U_VALUE(2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                U_DERIV(1,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S1)
                U_DERIV(1,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S2)
                U_DERIV(2,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,PART_DERIV_S1)
                U_DERIV(2,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,PART_DERIV_S2)
                IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                  U_VALUE(3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                  U_DERIV(3,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,PART_DERIV_S1)
                  U_DERIV(3,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,PART_DERIV_S2)
                  U_DERIV(3,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(3,PART_DERIV_S3)
                  U_DERIV(1,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,PART_DERIV_S3)
                  U_DERIV(2,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,PART_DERIV_S3)
                ELSE
                  U_VALUE(3)=0.0_DP
                  U_DERIV(3,1)=0.0_DP
                  U_DERIV(3,2)=0.0_DP
                  U_DERIV(3,3)=0.0_DP
                  U_DERIV(1,3)=0.0_DP
                  U_DERIV(2,3)=0.0_DP
                END IF
                !Start with calculation of partial matrices
                !Here W_VALUES must be ZERO if ALE part of linear matrix
                W_VALUE=0.0_DP
              END IF

              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
                !Loop over field components
                mhs=0

                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                  MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                  JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                    & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)

                  DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    nhs=0
                    IF(UPDATE_JACOBIAN_MATRIX) THEN
                      !Loop over element columns
                      DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                        MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                          &(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR                       
                        DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                          nhs=nhs+1
                          !Calculate some general values needed below
                          DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                            DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              DXI_DX(mi,ni)=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                                & DXI_DX(mi,ni)
                            END DO
                            DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                            DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          END DO !ni
                          PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                          SUM=0.0_DP
                          IF(UPDATE_JACOBIAN_MATRIX) THEN
                            !Calculate J1 only
                            DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              SUM=SUM+(PHINS*U_DERIV(mh,ni)*DXI_DX(ni,nh)*PHIMS*RHO)
                            END DO 
                            !Calculate MATRIX  
                            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs) &
                             & +SUM*JGW
                            !Calculate J2 only
                            IF(nh==mh) THEN 
                              SUM=0.0_DP
                              !Calculate SUM 
                              DO x=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO mi=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                  SUM=SUM+RHO*(U_VALUE(x)-W_VALUE(x))*DPHINS_DXI(mi)*DXI_DX(mi,x)*PHIMS
                                END DO !mi
                              END DO !x
                              !Calculate MATRIX
                              JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs) &
                               & +SUM*JGW
                            END IF
                          END IF
                        END DO !ns    
                      END DO !nh
                    END IF
                  END DO !ms
                END DO !mh
                ! Stabilisation terms
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                  & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                  & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
                  CALL NavierStokes_ResidualBasedStabilisation(EQUATIONS_SET,ELEMENT_NUMBER,ng,MU,RHO,.TRUE., &
                    & ERR,ERROR,*999)
                END IF
              END IF

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!                                        !!!!!
              !!!!!         1 D  T R A N S I E N T         !!!!!
              !!!!!                                        !!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              !Start with Matrix Calculations
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE     .OR. &
               & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE     .OR. &
               & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
               & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                Q_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                A_VALUE=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                Q_DERIV=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(1,FIRST_PART_DERIV)
                A_DERIV=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR%VALUES(2,FIRST_PART_DERIV)          
                A0=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                E =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                H =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                kp=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(4,NO_PART_DERIV)
                k1=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(5,NO_PART_DERIV)
                k2=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(6,NO_PART_DERIV)
                k3=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(7,NO_PART_DERIV)
                b1=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(8,NO_PART_DERIV)
                A0_DERIV=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(1,FIRST_PART_DERIV)
                E_DERIV =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(2,FIRST_PART_DERIV)
                H_DERIV =EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(3,FIRST_PART_DERIV)
                beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

                !Get material constants
                CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &                
                  & MU,err,error,*999)
                CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2, &
                  & RHO,err,error,*999)
                CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3, &
                  & alpha,err,error,*999)
                CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,8, &
                  & G0,err,error,*999)
                CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,9, &
                  & Fr,err,error,*999)

                !If A goes negative during nonlinear iteration, give ZERO_TOLERANCE value to avoid segfault
                IF (A_VALUE<A0*0.001_DP) A_VALUE=A0*0.001_DP

                mhs=0
                !Loop Over Element Rows
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                  JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                    & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                  ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(mh)%DOMAIN%TOPOLOGY%ELEMENTS
                  DXI_DX(1,1)=0.0_DP
                  !Calculate dxi_dx in 3D
                  DO xiIdx=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                    DO coordIdx=1,EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE) &
                      & %PTR%NUMBER_OF_X_DIMENSIONS
                      DXI_DX(1,1)=DXI_DX(1,1)+(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)% &
                        & PTR%DXI_DX(xiIdx,coordIdx))**2
                    ENDDO !coordIdx
                  ENDDO !xiIdx
                  DXI_DX(1,1)=SQRT(DXI_DX(1,1))
                  !Loop Over Element rows
                  DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                    PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    DPHIMS_DXI(1)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,FIRST_PART_DERIV,ng)
                    mhs=mhs+1
                    nhs=0
                    !Loop Over Element Columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                        & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        DPHINS_DXI(1)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ns,FIRST_PART_DERIV,ng)
                        nhs=nhs+1
                        IF (UPDATE_JACOBIAN_MATRIX) THEN

                          !Momentum Equation (dF/dQ)
                          IF (mh==1 .AND. nh==1) THEN
                            SUM=((2.0_DP*alpha*PHINS/A_VALUE*Q_DERIV+              & !Convective
                              & 2.0_DP*alpha*Q_VALUE/A_VALUE*DPHINS_DXI(1)-        & !Convective
                              & 2.0_DP*alpha*PHINS*Q_VALUE/(A_VALUE**2)*A_DERIV)*  & !Convective
                              & DXI_DX(1,1)+Fr*PHINS/A_VALUE)*PHIMS                  !Viscosity
                            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)= &
                              & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+SUM*JGW
                          ENDIF

                          !Momentum Equation (dF/dA)
                          IF (mh==1 .AND. nh==2) THEN
                            SUM=((2.0_DP*alpha*PHINS*(Q_VALUE**2)/(A_VALUE**3)*A_DERIV-             & !Convective
                              & 2.0_DP*alpha*PHINS*Q_VALUE/(A_VALUE**2)*Q_DERIV-                    & !Convective
                              & alpha*(Q_VALUE/A_VALUE)**2*DPHINS_DXI(1)+                           & !Convective
                              & beta/RHO*(b1*(A_VALUE/A0)**b1*DPHINS_DXI(1)+                        & !dA/dx
                              & b1*PHINS/A0*b1*(A_VALUE/A0)**(b1-1.0_DP)*A_DERIV+                   & !dA/dx
                              & k1*PHINS/A0*((b1+1.0_DP)*(1.0_DP-b1/k1)*(A_VALUE/A0)**b1)*A0_DERIV+ & !dA0/dx
                              & k2*PHINS/E*((b1+1.0_DP)*(A_VALUE/A0)**b1)*E_DERIV+                  & !dE/dx
                              & k3*PHINS/H*((b1+1.0_DP)*(A_VALUE/A0)**b1)*H_DERIV))*                & !dH/dx
                              & DXI_DX(1,1)-Fr*PHINS*Q_VALUE/A_VALUE**2)*PHIMS                        !Viscosity
                            JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)= &
                              & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+SUM*JGW
                          END IF
                        END IF
                      END DO !ns
                    END DO !nh
                  END DO !ms
                END DO !mh
              END IF
            END DO !ng

          IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
           & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
           & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
           & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN

            ELEMENTS_TOPOLOGY=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%components(1)%domain%topology%elements
            numberOfElementNodes=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_NODES
            numberOfParameters=ELEMENTS_TOPOLOGY%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS  
            firstNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(1)
            lastNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(numberOfElementNodes)
            compIdx=1
            derivIdx=1
            DO nodeIdx=1,numberOfElementNodes
              nodeNumber=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nodeIdx)
              numberOfVersions=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%NODES%NODES(nodeNumber)%DERIVATIVES(derivIdx)%numberOfVersions

              !!!-- B R A N C H   F L U X   U P W I N D I N G --!!!
              !----------------------------------------------------
              ! In order to enforce conservation of mass and momentum across discontinuous
              ! branching topologies, flux is upwinded against the conservative branch values 
              ! established by the characteristic solver.

              !Find the branch node on this element
              IF(numberOfVersions>1) THEN                  
                versionIdx=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)% &
                  & elementVersions(derivIdx,nodeIdx)
                !Find the wave direction - incoming or outgoing
                CALL Field_ParameterSetGetLocalNode(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, & 
                  & versionIdx,derivIdx,nodeNumber,compIdx,normalWave,err,error,*999)
                !Get materials variables for node on this element
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,1,A0,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,2,E,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,3,H,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,4,kp,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,5,k1,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,6,k2,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,7,k3,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,8,b1,err,error,*999)
                beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)
                !Get current Q & A values
                CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,1,Q_VALUE,err,error,*999)         
                CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivIdx,nodeNumber,2,A_VALUE,err,error,*999)

                !Momentum Equation, d/dQ
                momentum1=(-2.0_DP*alpha*Q_VALUE/A_VALUE)*normalWave
                !Momentum Equation, d/dA
                momentum2=(alpha*(Q_VALUE/A_VALUE)**2-beta*b1/RHO*(A_VALUE/A0)**b1)*normalWave
                !Continuity Equation, d/dQ
                mass=(-1.0_DP)*normalWave

                !Add momentum/mass contributions to first/last node accordingly
                IF(nodeNumber==firstNode) THEN
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1,1)= &
                    & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1,1)+momentum1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1,numberOfParameters+1)= &
                    & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(1,numberOfParameters+1)+momentum2
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(numberOfParameters+1,1)= &
                    & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(numberOfParameters+1,1)+mass
                ELSE IF (nodeNumber==lastNode) THEN
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(numberOfParameters,numberOfParameters)= &
                    & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(numberOfParameters,numberOfParameters)+momentum1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(numberOfParameters,2*numberOfParameters)= &
                    & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(numberOfParameters,2*numberOfParameters)+momentum2
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(2*numberOfParameters,numberOfParameters)= &
                    & JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(2*numberOfParameters,numberOfParameters)+mass
                END IF
              END IF !find brach node
            END DO !loop over nodes
          END IF !equation set
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokes_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementJacobianEvaluate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NavierStokes_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes problem post solve.
  SUBROUTINE NAVIER_STOKES_POST_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    INTEGER(INTG) :: iteration,timestep,outputIteration,equationsSetNumber
    REAL(DP) :: startTime,stopTime,currentTime,timeIncrement
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NAVIER_STOKES_POST_SOLVE",ERR,ERROR,*999)
    NULLIFY(SOLVER2)
    NULLIFY(SOLVERS)
    NULLIFY(dependentField)
    NULLIFY(fieldVariable)

    IF (ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
            IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
              SELECT CASE(SOLVER%SOLVE_TYPE)
              CASE(SOLVER_NONLINEAR_TYPE)
                ! Characteristic solver- copy branch Q,A values to new parameter set
                dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                END IF
                iteration=CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
                IF(iteration==1) THEN
                  CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
                END IF
              CASE(SOLVER_DYNAMIC_TYPE)
                ! Navier-Stokes solver: do nothing
              CASE DEFAULT
                LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for a 1D Navier-Stokes problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                  CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)

                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  END IF
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                IF(SOLVER%GLOBAL_NUMBER == 1) THEN
                  ! DAE solver- do nothing
                ELSE
                  LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the CellML DAE simple loop of a 1D0D coupled Navier-Stokes problem."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END IF
              ELSE
                LOCAL_ERROR="The control loop type for solver "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                  CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)

                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  END IF
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                IF(SOLVER%GLOBAL_NUMBER == 1) THEN
                  ! DAE solver- do nothing
                ELSE
                  LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the CellML DAE simple loop of a 1D0D coupled Navier-Stokes problem."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END IF
              ELSE
                LOCAL_ERROR="The control loop type for solver "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
              SELECT CASE(SOLVER%GLOBAL_NUMBER)
              CASE(1)
                ! Characteristic solver- copy branch Q,A values to new parameter set
                dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                 & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
              CASE(2)
                ! check characteristic/ N-S convergence at branches
  !                CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              CASE(3)
                ! Advection solver output data if necessary
                IF(CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP .EQV. .FALSE.) THEN
                  ! 1D NSE solver output data if N-S/Chars converged
                  CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
                END IF
              CASE DEFAULT
                LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for a 1D Navier-Stokes and Advection problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                  CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  END IF
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                ! DAE and advection solvers - output data if post advection solve
                IF(SOLVER%SOLVERS%CONTROL_LOOP%SUB_LOOP_INDEX == 3) THEN
                  CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
                END IF
              ELSE
                LOCAL_ERROR="The control loop type for solver "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                  CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SETS_COPY(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  END IF
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                ! DAE and advection solvers - output data if post advection solve
                IF(SOLVER%SOLVERS%CONTROL_LOOP%SUB_LOOP_INDEX == 3) THEN
                  CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
                END IF
              ELSE
                LOCAL_ERROR="The control loop type for solver "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
              CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,startTime,stopTime,currentTime,timeIncrement, &
               & timestep,outputIteration,ERR,ERROR,*999)
              CALL NavierStokes_CalculateBoundaryFlux(SOLVER,ERR,ERROR,*999)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              CALL NavierStokes_CalculateBoundaryFlux(SOLVER,ERR,ERROR,*999)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
              DO equationsSetNumber=1,SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                ! If this is a coupled constitutive (non-Newtonian) viscosity problem, update shear rate values
                !  to be passed to the CellML solver at beginning of next timestep
                IF(SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equationsSetNumber)%PTR% &
                 &  EQUATIONS%EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                 CALL NavierStokes_ShearRateCalculate(SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                  & EQUATIONS_SETS(equationsSetNumber)%PTR%EQUATIONS%EQUATIONS_SET,err,error,*999)
                END IF
              END DO
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Post solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement post solve... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,ERR,ERROR,*999)
                IF(ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                  SOLVER2%DYNAMIC_SOLVER%ALE=.TRUE.
                ELSE  
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",ERR,ERROR,*999)
                END IF
              !Post solve for the dynamic solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"ALE Navier-Stokes post solve... ",ERR,ERROR,*999)
                CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="The third problem specification of  "// &
                & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%specification(3),"*",ERR,ERROR))// &
                & " is not valid for a Navier-Stokes fluid mechanics problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
          END IF
        ELSE
          CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    END IF

    EXITS("NAVIER_STOKES_POST_SOLVE")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NAVIER_STOKES_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Navier-Stokes flow pre solve
  SUBROUTINE NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET,SOLID_EQUATIONS_SET,FLUID_EQUATIONS_SET
    TYPE(EQUATIONS_SET_DEPENDENT_TYPE), POINTER :: SOLID_DEPENDENT
    TYPE(EQUATIONS_SET_GEOMETRY_TYPE), POINTER :: FLUID_GEOMETRIC
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS,SOLID_EQUATIONS,FLUID_EQUATIONS
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD,SOLID_DEPENDENT_FIELD,FLUID_GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentFieldVariable,independentFieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,SOLID_SOLVER_EQUATIONS,FLUID_SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING,SOLID_SOLVER_MAPPING,FLUID_SOLVER_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: Solver2
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    INTEGER(INTG) :: nodeIdx,derivIdx,versionIdx,variableIdx,numberOfSourceTimesteps,timeIdx,compIdx
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE,GLOBAL_DERIV_INDEX,node_idx,variable_type
    INTEGER(INTG) :: variable_idx,local_ny,ANALYTIC_FUNCTION_TYPE,component_idx,deriv_idx,dim_idx,version_idx
    INTEGER(INTG) :: element_idx,en_idx,I,J,K,number_of_nodes_xic(3),search_idx,localDof,globalDof,componentBC,previousNodeNumber
    INTEGER(INTG) :: componentNumberVelocity,numberOfDimensions,numberOfNodes,numberOfGlobalNodes,currentLoopIteration
    INTEGER(INTG) :: dependentVariableType,independentVariableType,dependentDof,independentDof,userNodeNumber,localNodeNumber
    INTEGER(INTG) :: EquationsSetIndex,SolidNodeNumber,FluidNodeNumber
    INTEGER(INTG), ALLOCATABLE :: InletNodes(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,DISPLACEMENT_VALUE,VALUE,XI_COORDINATES(3),timeData,QP,QPP,componentValues(3)
    REAL(DP) :: T_COORDINATES(20,3),MU,RHO,X(3),FluidGFValue,SolidDFValue,NewLaplaceBoundaryValue,lengthScale,timeScale,massScale
    REAL(DP), POINTER :: MESH_VELOCITY_VALUES(:), GEOMETRIC_PARAMETERS(:), BOUNDARY_VALUES(:)
    REAL(DP), POINTER :: TANGENTS(:,:),NORMAL(:),TIME,ANALYTIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    REAL(DP), POINTER :: independentParameters(:),dependentParameters(:) 
    REAL(DP), ALLOCATABLE :: nodeData(:,:),qSpline(:),qValues(:),tValues(:),BoundaryValues(:)
    LOGICAL :: ghostNode,nodeExists,importDataFromFile,ALENavierStokesEquationsSetFound=.FALSE.
    LOGICAL :: SolidEquationsSetFound=.FALSE.,SolidNodeFound=.FALSE.,FluidEquationsSetFound=.FALSE.
    CHARACTER(70) :: inputFile,tempString
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(ANALYTIC_FIELD)
    NULLIFY(DEPENDENT_FIELD)
    NULLIFY(GEOMETRIC_FIELD)
    NULLIFY(MATERIALS_FIELD)
    NULLIFY(INDEPENDENT_FIELD)
    NULLIFY(ANALYTIC_VARIABLE)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(GEOMETRIC_VARIABLE)
    NULLIFY(MATERIALS_VARIABLE)
    NULLIFY(DOMAIN)
    NULLIFY(DOMAIN_NODES)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(MESH_VELOCITY_VALUES)
    NULLIFY(GEOMETRIC_PARAMETERS)
    NULLIFY(BOUNDARY_VALUES)
    NULLIFY(TANGENTS)
    NULLIFY(NORMAL)
    NULLIFY(TIME)
    NULLIFY(ANALYTIC_PARAMETERS)
    NULLIFY(MATERIALS_PARAMETERS)
    NULLIFY(independentParameters)
    NULLIFY(dependentParameters)

    ENTERS("NavierStokes_PreSolveUpdateBoundaryConditions",ERR,ERROR,*999)

    IF (ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)        
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF (ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET

                  !Fitting boundary condition- get values from file
                  !TODO: this should be generalised with input filenames specified from the example file when IO is improved
                  IF (ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN                
                    !Read in field values to independent field
                    NULLIFY(independentFieldVariable)
                    NULLIFY(dependentFieldVariable)
                    INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                    DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    independentVariableType=INDEPENDENT_FIELD%VARIABLES(1)%VARIABLE_TYPE
                    CALL FIELD_VARIABLE_GET(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,independentFieldVariable,ERR,ERROR,*999)
                    dependentVariableType=DEPENDENT_FIELD%VARIABLES(1)%VARIABLE_TYPE
                    CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,dependentFieldVariable,ERR,ERROR,*999)
                    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                      & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                    !Read in field data from file
                    !Loop over nodes and update independent field values - if a fixed fitted boundary, also update dependent
                    IF (ASSOCIATED(INDEPENDENT_FIELD)) THEN
                      componentNumberVelocity=1
                      numberOfDimensions=dependentFieldVariable%NUMBER_OF_COMPONENTS - 1
                      !Get the nodes on this computational domain
                      IF (independentFieldVariable%COMPONENTS(componentNumberVelocity)%INTERPOLATION_TYPE== &
                        & FIELD_NODE_BASED_INTERPOLATION) THEN
                        domain=>independentFieldVariable%COMPONENTS(componentNumberVelocity)%DOMAIN
                        IF (ASSOCIATED(domain)) THEN
                          IF (ASSOCIATED(domain%TOPOLOGY)) THEN
                            DOMAIN_NODES=>domain%TOPOLOGY%NODES
                            IF (ASSOCIATED(DOMAIN_NODES)) THEN
                              numberOfNodes=DOMAIN_NODES%NUMBER_OF_NODES
                              numberOfGlobalNodes=DOMAIN_NODES%NUMBER_OF_GLOBAL_NODES
                            ELSE
                              CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                            END IF
                          ELSE
                            CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                          END IF
                        ELSE
                          CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                      END IF

                      !Construct the filename based on the computational node and time step
                      currentLoopIteration=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                      WRITE(tempString,"(I4.4)") currentLoopIteration
                      inputFile = './../interpolatedData/fitData' // tempString(1:4) // '.dat'

                      INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                      IF (importDataFromFile) THEN
                        !Read fitted data from input file (if exists)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Updating independent field and boundary nodes from "//inputFile, &
                          & ERR,ERROR,*999)
                        OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')                  
                        !Loop over local nodes and update independent field and (and dependent field for any FIXED_FITTED nodes)
                        previousNodeNumber=0
                        DO nodeIdx=1,numberOfNodes
                          userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                          CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(domain%Topology,userNodeNumber,nodeExists,localNodeNumber, &
                            & ghostNode,err,error,*999)
                          IF(nodeExists .AND. .NOT. ghostNode) THEN
                            ! Move to line in file for this node (dummy read)
                            ! NOTE: this takes advantage of the user number increasing ordering of domain nodes 
                            DO search_idx=1,userNodeNumber-previousNodeNumber-1
                              READ(10,*)
                            END DO
                            ! Read in the node data for this timestep file
                            READ(10,*) (componentValues(compIdx), compIdx=1,numberOfDimensions)
                            DO compIdx=1,numberOfDimensions
                              dependentDof=dependentFieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                & NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                              independentDof=independentFieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                              VALUE = componentValues(compIdx)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INDEPENDENT_FIELD,independentVariableType, &
                                & FIELD_VALUES_SET_TYPE,independentDof,VALUE,ERR,ERROR,*999)
                              CALL FIELD_COMPONENT_DOF_GET_USER_NODE(DEPENDENT_FIELD,dependentVariableType,1,1,userNodeNumber, & 
                                & compIdx,localDof,globalDof,err,error,*999)
                              BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(globalDof)
                              IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,dependentVariableType, &
                                  & FIELD_VALUES_SET_TYPE,localDof,VALUE,ERR,ERROR,*999)
                              END IF
                            END DO !compIdx
                            previousNodeNumber=userNodeNumber
                          END IF !ghost/exist check
                        END DO !nodeIdx
                        CLOSE(UNIT=10)
                      END IF !check import file exists
                    ELSE
                      CALL FlagError("Equations set independent field is not associated.",ERR,ERROR,*999)
                    END IF
                  END IF !Equations set independent

                  !Analytic equations
                  IF (ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                    !Standard analytic functions
                    IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                      ! Update analytic time value with current time
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                      !Calculate analytic values
                      BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                      END IF
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                      & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE) THEN
                      IF (ASSOCIATED(EQUATIONS_SET)) THEN
                        IF (ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          IF (ASSOCIATED(DEPENDENT_FIELD)) THEN
                            GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                            IF (ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                              !Geometric parameters
                              CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR, & 
                                & ERROR,*999)
                              NULLIFY(GEOMETRIC_VARIABLE)
                              NULLIFY(GEOMETRIC_PARAMETERS)
                              CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
                              CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                              !Analytic parameters
                              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                              NULLIFY(ANALYTIC_VARIABLE)
                              NULLIFY(ANALYTIC_PARAMETERS)
                              IF (ASSOCIATED(ANALYTIC_FIELD)) THEN
                                CALL FIELD_VARIABLE_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,ANALYTIC_VARIABLE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_DATA_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & ANALYTIC_PARAMETERS,ERR,ERROR,*999)           
                              END IF
                              ! Materials parameters
                              NULLIFY(MATERIALS_FIELD)
                              NULLIFY(MATERIALS_VARIABLE)
                              NULLIFY(MATERIALS_PARAMETERS)
                              IF (ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
                                MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                                CALL FIELD_VARIABLE_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_DATA_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & MATERIALS_PARAMETERS,ERR,ERROR,*999)           
                              END IF
                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                                variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                                FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                                IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                                  CALL Field_ParameterSetEnsureCreated(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                    IF (FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
                                      & FIELD_NODE_BASED_INTERPOLATION) THEN
                                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                      IF (ASSOCIATED(DOMAIN)) THEN
                                        IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                          IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                            !Should be replaced by boundary node flag
                                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                                !Default to version 1 of each node derivative
                                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                              END DO !dim_idx

                                              !Loop over the derivatives
                                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                                  & GLOBAL_DERIVATIVE_INDEX
                                                CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(ANALYTIC_FUNCTION_TYPE,X, &
                                                  & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX,compIdx, &
                                                  & NUMBER_OF_DIMENSIONS,FIELD_VARIABLE%NUMBER_OF_COMPONENTS,ANALYTIC_PARAMETERS, &
                                                  & MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
                                                DO version_idx=1, &
                                                  & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%numberOfVersions
                                                  local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                                    & VERSIONS(version_idx)
                                                  !Set analytic values
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                                  CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                                                    & DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR, &
                                                    & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                                                  IF (EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                                                    & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN) THEN
                                                    !Taylor-Green boundary conditions update
                                                    IF (ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                                      BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                        & CONDITION_TYPES(local_ny)
                                                      IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED .AND. &
                                                        & component_idx<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, &
                                                          & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
                                                          & VALUE,ERR,ERROR,*999)
                                                      END IF !Boundary condition fixed
                                                    END IF !Boundary condition variable
                                                  END IF ! Taylor-Green
                                                END DO !version_idx
                                              END DO !deriv_idx
                                            END DO !node_idx
                                          ELSE
                                            CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                          END IF
                                        ELSE
                                          CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                        END IF
                                      ELSE
                                        CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                    END IF
                                  END DO !component_idx
                                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                ELSE
                                  CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
                                END IF
                              END DO !variable_idx
                              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                                & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                            ELSE
                              CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
                            END IF
                          ELSE
                            CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                          END IF
                        ELSE
                          CALL FlagError("Equations set analytic is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                      END IF
                      ! Unit shape analytic functions
                    ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5 .OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                      IF (ASSOCIATED(EQUATIONS_SET)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF (ASSOCIATED(DEPENDENT_FIELD)) THEN
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF (ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                            ! Geometric parameters
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                             & ERR,ERROR,*999)
                            NULLIFY(GEOMETRIC_VARIABLE)
                            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
                            NULLIFY(GEOMETRIC_PARAMETERS)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, & 
                              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                            ! Analytic parameters
                            ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                            ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                            NULLIFY(ANALYTIC_VARIABLE)
                            NULLIFY(ANALYTIC_PARAMETERS)
                            IF (ASSOCIATED(ANALYTIC_FIELD)) THEN
                              CALL FIELD_VARIABLE_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,ANALYTIC_VARIABLE,ERR,ERROR,*999)
                              CALL FIELD_PARAMETER_SET_DATA_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & ANALYTIC_PARAMETERS,ERR,ERROR,*999)           
                            END IF
                            ! Materials parameters
                            NULLIFY(MATERIALS_FIELD)
                            NULLIFY(MATERIALS_VARIABLE)
                            NULLIFY(MATERIALS_PARAMETERS)
                            IF (ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
                              MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                              CALL FIELD_VARIABLE_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,ERR,ERROR,*999)
                              CALL FIELD_PARAMETER_SET_DATA_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & MATERIALS_PARAMETERS,ERR,ERROR,*999)           
                            END IF
                            TIME=EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME
                            ! Interpolation parameters
                            NULLIFY(INTERPOLATION_PARAMETERS)
                            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(GEOMETRIC_FIELD,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                            NULLIFY(INTERPOLATED_POINT) 
                            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                             & ERR,ERROR,*999)
                            DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                              IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                                DO compIdx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  IF (FIELD_VARIABLE%COMPONENTS(compIdx)%INTERPOLATION_TYPE== & 
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
                                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(compIdx)%DOMAIN
                                    IF (ASSOCIATED(DOMAIN)) THEN
                                      IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                        IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                          !Should be replaced by boundary node flag
                                          DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                            element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surrounding_elements(1)
                                            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                                              & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                                            en_idx=0
                                            XI_COORDINATES=0.0_DP
                                            number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)% &
                                              & basis%number_of_nodes_xic(1)
                                            number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)% & 
                                              & basis%number_of_nodes_xic(2)
                                            IF (NUMBER_OF_DIMENSIONS==3) THEN
                                              number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis% &
                                                & number_of_nodes_xic(3)
                                            ELSE
                                              number_of_nodes_xic(3)=1
                                            END IF
!\todo: change definitions as soon as adjacent elements / boundary elements calculation works for simplex
                                            IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4 .OR. &
                                              & DOMAIN%topology%elements%maximum_number_of_element_parameters==9 .OR. &
                                              & DOMAIN%topology%elements%maximum_number_of_element_parameters==16 .OR. &
                                              & DOMAIN%topology%elements%maximum_number_of_element_parameters==8 .OR. &
                                              & DOMAIN%topology%elements%maximum_number_of_element_parameters==27 .OR. &
                                              & DOMAIN%topology%elements%maximum_number_of_element_parameters==64) THEN
                                              DO K=1,number_of_nodes_xic(3)
                                                DO J=1,number_of_nodes_xic(2)
                                                  DO I=1,number_of_nodes_xic(1)
                                                    en_idx=en_idx+1
                                                    IF (DOMAIN%topology%elements%elements(element_idx)% & 
                                                      & element_nodes(en_idx)==node_idx) EXIT
                                                    XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                                  END DO !I
                                                  IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                    & element_nodes(en_idx)==node_idx) EXIT
                                                  XI_COORDINATES(1)=0.0_DP
                                                  XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                                END DO !J
                                                IF(DOMAIN%topology%elements%elements(element_idx)% & 
                                                  & element_nodes(en_idx)==node_idx) EXIT
                                                XI_COORDINATES(1)=0.0_DP
                                                XI_COORDINATES(2)=0.0_DP
                                                IF (number_of_nodes_xic(3)/=1) THEN
                                                  XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                                END IF
                                              END DO !K
                                              CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES, &
                                                & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                                            ELSE
                                              !\todo: Use boundary flag
                                              IF (DOMAIN%topology%elements%maximum_number_of_element_parameters==3) THEN
                                                T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                              ELSE IF (DOMAIN%topology%elements%maximum_number_of_element_parameters==6) THEN
                                                T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                                T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                                T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                              ELSE IF (DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. & 
                                                & NUMBER_OF_DIMENSIONS==2) THEN
                                                T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                                T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                                T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                              ELSE IF (DOMAIN%topology%elements%maximum_number_of_element_parameters==4) THEN
                                                T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                              ELSE IF (DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. & 
                                                & NUMBER_OF_DIMENSIONS==3) THEN
                                                T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                                T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                                T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                                T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                                T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                              ELSE IF (DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
                                                T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                                T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                                T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                                T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                                T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                              END IF
                                              DO K=1,DOMAIN%topology%elements%maximum_number_of_element_parameters
                                                IF(DOMAIN%topology%elements%elements(element_idx)%element_nodes(K)==node_idx) EXIT
                                              END DO !K
                                              IF(NUMBER_OF_DIMENSIONS==2) THEN
                                                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                                  & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                                              ELSE IF (NUMBER_OF_DIMENSIONS==3) THEN
                                                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                                  & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                                              END IF
                                            END IF
                                            X=0.0_DP
                                            DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                              X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(dim_idx,1)
                                            END DO !dim_idx
                                            !Loop over the derivatives
                                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                              ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                              GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                                & GLOBAL_DERIVATIVE_INDEX
                                              MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                                              !Define MU, density=1
                                              MU=MATERIALS_FIELD%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                                & parameters%cmiss%data_dp(1)
                                              !Define RHO, density=2
                                              RHO=MATERIALS_FIELD%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                                & parameters%cmiss%data_dp(2)
                                              CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(ANALYTIC_FUNCTION_TYPE,X, &
                                                & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX,compIdx,NUMBER_OF_DIMENSIONS,& 
                                                & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,ANALYTIC_PARAMETERS, &
                                                & MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*999)
                                              !Default to version 1 of each node derivative
                                              local_ny=FIELD_VARIABLE%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_FIELD% &
                                                & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR,BOUNDARY_CONDITIONS_VARIABLE, &
                                                & ERR,ERROR,*999)
                                              IF (ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                                BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                  & CONDITION_TYPES(local_ny)
                                                IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, &
                                                   & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
                                                   & VALUE,ERR,ERROR,*999)
                                                END IF
                                              ELSE
                                                CALL FlagError("Boundary conditions U variable is not associated.", &
                                                  & ERR,ERROR,*999)
                                              END IF
                                            END DO !deriv_idx
                                          END DO !node_idx
                                        ELSE
                                          CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                        END IF
                                      ELSE
                                        CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                  END IF
                                END DO !compIdx
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                              ELSE
                                CALL FlagError("Field variable is not associated.",ERR,ERROR,*999)
                              END IF
                            END DO !variable_idx
                            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                             & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Equations set geometric field is not associated.",ERR,ERROR,*999)
                          END IF
                        ELSE
                          CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                      END IF
                    END IF !Standard/unit analytic subtypes

                  END IF ! Analytic boundary conditions

                  !TODO implement non-analytic time-varying boundary conditions (i.e. from file)
                ELSE
                  CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
              END IF
              CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                !If analytic flow waveform, calculate and update
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          ! Calculate analytic values
                          CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                        END IF
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                        ! Perform spline interpolation of values from a file
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                        DO variableIdx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                          dependentVariableType=DEPENDENT_FIELD%VARIABLES(variableIdx)%VARIABLE_TYPE
                          NULLIFY(dependentFieldVariable)
                          CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,dependentVariableType,dependentFieldVariable,ERR,ERROR,*999)
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
                           & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                            IF(ASSOCIATED(dependentFieldVariable)) THEN
                              DO compIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                IF(dependentFieldVariable%COMPONENTS(compIdx)%INTERPOLATION_TYPE== &
                                 & FIELD_NODE_BASED_INTERPOLATION) THEN
                                  domain=>dependentFieldVariable%COMPONENTS(compIdx)%DOMAIN
                                  IF(ASSOCIATED(domain)) THEN
                                    IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                      DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                        !Loop over the local nodes excluding the ghosts.
                                        DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                          userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                          DO derivIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                            DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)%DERIVATIVES(derivIdx)%numberOfVersions
                                              ! Update analytic field if file exists and dependent field if boundary condition set
                                              inputFile = './input/interpolatedData/1D/' 
                                              IF(dependentVariableType == FIELD_U_VARIABLE_TYPE) THEN
                                                inputFile = TRIM(inputFile) // 'U/component' 
                                              END IF
                                              WRITE(tempString,"(I1.1)") compIdx 
                                              inputFile = TRIM(inputFile) // tempString(1:1) // '/derivative'
                                              WRITE(tempString,"(I1.1)") derivIdx 
                                              inputFile = TRIM(inputFile) // tempString(1:1) // '/version'
                                              WRITE(tempString,"(I1.1)") versionIdx 
                                              inputFile = TRIM(inputFile) // tempString(1:1) // '/'
                                              WRITE(tempString,"(I4.4)") userNodeNumber
                                              inputFile = TRIM(inputFile) // tempString(1:4) // '.dat'
                                              inputFile = TRIM(inputFile)
                                              INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                                              IF(importDataFromFile) THEN
                                                ! Create the analytic field values type on the dependent field if it does not exist
                                                IF(.NOT.ASSOCIATED(dependentFieldVariable%PARAMETER_SETS% &
                                                  & SET_TYPE(FIELD_ANALYTIC_VALUES_SET_TYPE)%PTR)) &
                                                  & CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,dependentVariableType, &
                                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                                !Read fitted data from input file (if exists)
                                                OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')                  
                                                ! Header timeData = numberOfTimesteps
                                                READ(10,*) timeData
                                                numberOfSourceTimesteps = INT(timeData)
                                                ALLOCATE(nodeData(numberOfSourceTimesteps,2))
                                                ALLOCATE(qValues(numberOfSourceTimesteps))
                                                ALLOCATE(tValues(numberOfSourceTimesteps))
                                                ALLOCATE(qSpline(numberOfSourceTimesteps))
                                                nodeData = 0.0_DP                                            
                                                ! Read in time and dependent value
                                                DO timeIdx=1,numberOfSourceTimesteps
                                                  READ(10,*) (nodeData(timeIdx,component_idx), component_idx=1,2)
                                                END DO
                                                CLOSE(UNIT=10)
                                                tValues = nodeData(:,1)
                                                qValues = nodeData(:,2)
                                                CALL spline_cubic_set(numberOfSourceTimesteps,tValues,qValues,2,0.0_DP,2,0.0_DP, &
                                                  & qSpline,err,error,*999)
                                                CALL spline_cubic_val(numberOfSourceTimesteps,tValues,qValues,qSpline, & 
                                                  & CURRENT_TIME,VALUE,QP,QPP,err,error,*999)

                                                DEALLOCATE(nodeData)
                                                DEALLOCATE(qSpline)
                                                DEALLOCATE(qValues)
                                                DEALLOCATE(tValues)

                                                dependentDof = dependentFieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                                 & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,dependentVariableType, &
                                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,dependentDof,VALUE,ERR,ERROR,*999)
                                                ! Update dependent field value if this is a splint BC
                                                BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                  & CONDITION_TYPES(dependentDof)
                                                IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_FITTED) THEN 
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,dependentVariableType, &
                                                    & FIELD_VALUES_SET_TYPE,dependentDof,VALUE,ERR,ERROR,*999)
                                                END IF
                                              END IF ! check if import data file exists
                                            END DO !versionIdx
                                          END DO !derivIdx
                                        END DO !nodeIdx
                                      ELSE
                                        CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                  END IF
                                ELSE
                                  CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                END IF
                              END DO !compIdx
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",ERR,ERROR,*999)
                            END IF
                          END IF
                        END DO !variableIdx
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART)
                        ! Using heart lumped parameter model for input
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                        CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,5, &
                          & lengthScale,err,error,*999)
                        CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,6, &
                          & timeScale,err,error,*999)
                        CALL FIELD_PARAMETER_SET_GET_CONSTANT(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,7, &
                          & massScale,err,error,*999)
                        DO variableIdx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                          dependentVariableType=DEPENDENT_FIELD%VARIABLES(variableIdx)%VARIABLE_TYPE
                          NULLIFY(dependentFieldVariable)
                          CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,dependentVariableType,dependentFieldVariable,ERR,ERROR,*999)
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,dependentFieldVariable, &
                            & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                            IF(ASSOCIATED(dependentFieldVariable)) THEN
                              DO compIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                IF(dependentFieldVariable%COMPONENTS(compIdx)%INTERPOLATION_TYPE== &
                                 & FIELD_NODE_BASED_INTERPOLATION) THEN
                                  domain=>dependentFieldVariable%COMPONENTS(compIdx)%DOMAIN
                                  IF(ASSOCIATED(domain)) THEN
                                    IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                      DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                        !Loop over the local nodes excluding the ghosts.
                                        DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                          userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                          DO derivIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                            DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)%DERIVATIVES(derivIdx)%numberOfVersions
                                              dependentDof=dependentFieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
                                              BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                & CONDITION_TYPES(dependentDof)
                                              IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                                CALL Field_ParameterSetGetLocalNode(DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                                                  & FIELD_VALUES_SET_TYPE,versionIdx,derivIdx,userNodeNumber,1,VALUE, &
                                                  & ERR,ERROR,*999)
                                                ! Convert Q from ml/s to non-dimensionalised form.
                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,dependentVariableType, &
                                                  & FIELD_VALUES_SET_TYPE,dependentDof,VALUE*(lengthScale**3.0)/timeScale, &
                                                  & ERR,ERROR,*999)
                                              END IF
                                            END DO !versionIdx
                                          END DO !derivIdx
                                        END DO !nodeIdx
                                      ELSE
                                        CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                  END IF
                                ELSE
                                  CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                END IF
                              END DO !compIdx
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",ERR,ERROR,*999)
                            END IF
                          END IF
                        END DO !variableIdx
                      CASE DEFAULT
                        ! Do nothing (might have another use for analytic equations)
                      END SELECT
                    END IF ! Check for analytic equations
                  ELSE
                    CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                END IF
              END IF ! solver equations associated
              ! Update any multiscale boundary values (coupled 0D or non-reflecting)
              CALL NavierStokes_UpdateMultiscaleBoundary(solver,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                      FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                          & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("Field U variable is not associated",ERR,ERROR,*999)
                      END IF
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                       CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                       NULLIFY(BOUNDARY_VALUES)
                       CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                         & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                       CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_NONLINEAR_TYPE,BOUNDARY_VALUES, &
                         & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                         & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP)
                       DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                         variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                         FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                         IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                           DO compIdx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                             DOMAIN=>FIELD_VARIABLE%COMPONENTS(compIdx)%DOMAIN
                             IF(ASSOCIATED(DOMAIN)) THEN
                               IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                 DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                 IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                   !Loop over the local nodes excluding the ghosts.
                                   DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                     DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                       !Default to version 1 of each node derivative
                                       local_ny=FIELD_VARIABLE%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                         & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                       BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                         & CONDITION_TYPES(local_ny)
                                       IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                         CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                           & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                           & BOUNDARY_VALUES(local_ny),ERR,ERROR,*999)
                                       END IF
                                     END DO !deriv_idx
                                   END DO !node_idx
                                 END IF
                               END IF
                             END IF
                           END DO !compIdx
                         END IF
                       END DO !variable_idx

  !\todo: This part should be read in out of a file eventually
                     ELSE
                       CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                     END IF
                   ELSE
                     CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                   END IF
                 ELSE
                   CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                 END IF
               ELSE
                 CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
               END IF
              ELSE
                CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
              END IF
              CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
             !Pre solve for the dynamic solver
             IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
               CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF (ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF (ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      IF (ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                        IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                            & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Field U variable is not associated",ERR,ERROR,*999)
                        END IF
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                          NULLIFY(MESH_VELOCITY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,ERR,ERROR,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, & 
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO compIdx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(compIdx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          DISPLACEMENT_VALUE=0.0_DP
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% & 
                                            & CONDITION_TYPES(local_ny)
                                          IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & MESH_VELOCITY_VALUES(local_ny),ERR,ERROR,*999)
                                          ELSE IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & BOUNDARY_VALUES(local_ny),ERR,ERROR,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !compIdx
                            END IF
                          END DO !variable_idx
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                  END IF                
                ELSE
                  CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                END IF  
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              IF (SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF (ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF (ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      IF (ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                        IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                            & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Field U variable is not associated",ERR,ERROR,*999)
                        END IF
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, & 
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_MOVED_WALL,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                            IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF (ASSOCIATED(DOMAIN)) THEN
                                  IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% & 
                                            & CONDITION_TYPES(local_ny)
                                          IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & BOUNDARY_VALUES(local_ny),ERR,ERROR,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                          END DO !variable_idx
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                          !\todo: This part should be read in out of a file eventually
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                !Pre solve for the dynamic solver
              ELSE IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF (ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF (ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                      IF (ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                        IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                            & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Field U variable is not associated",ERR,ERROR,*999)
                        END IF
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                          NULLIFY(MESH_VELOCITY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,ERR,ERROR,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, & 
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                            IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF (ASSOCIATED(DOMAIN)) THEN
                                  IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          DISPLACEMENT_VALUE=0.0_DP
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% & 
                                            & CONDITION_TYPES(local_ny)
                                          IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & MESH_VELOCITY_VALUES(local_ny),ERR,ERROR,*999)
                                          ELSE IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & BOUNDARY_VALUES(local_ny),ERR,ERROR,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                          END DO !variable_idx
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,ERR,ERROR,*999)
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              ENDIF
              !do nothing ???
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
                NULLIFY(Solver2)
                !Pre solve for the linear solver
                IF (SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",ERR,ERROR,*999)
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                    EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF (ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                      IF (ASSOCIATED(EQUATIONS_SET)) THEN
                        BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                        IF (ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                          IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                              & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                          ELSE
                            CALL FlagError("Field U variable is not associated",ERR,ERROR,*999)
                          END IF
                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                            !Update moving wall nodes from solid/fluid gap (as we solve for displacements of the mesh
                            !in Laplacian smoothing step).
                            CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,Solver2,ERR,ERROR,*999)
                            IF(.NOT.ASSOCIATED(Solver2)) CALL FlagError("Dynamic solver is not associated.",Err,Error,*999)
                            !Find the FiniteElasticity equations set as there is a NavierStokes equations set too
                            SOLID_SOLVER_EQUATIONS=>Solver2%SOLVER_EQUATIONS
                            IF (ASSOCIATED(SOLID_SOLVER_EQUATIONS)) THEN
                              SOLID_SOLVER_MAPPING=>SOLID_SOLVER_EQUATIONS%SOLVER_MAPPING
                              IF (ASSOCIATED(SOLID_SOLVER_MAPPING)) THEN
                                EquationsSetIndex=1
                                SolidEquationsSetFound=.FALSE.
                                DO WHILE (.NOT.SolidEquationsSetFound &
                                  & .AND.EquationsSetIndex<=SOLID_SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                                  SOLID_EQUATIONS=>SOLID_SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(EquationsSetIndex)%EQUATIONS
                                  IF (ASSOCIATED(SOLID_EQUATIONS)) THEN
                                    SOLID_EQUATIONS_SET=>SOLID_EQUATIONS%EQUATIONS_SET
                                    IF(ASSOCIATED(SOLID_EQUATIONS_SET)) THEN
                                      IF(SOLID_EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_ELASTICITY_CLASS &
                                        & .AND.SOLID_EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE &
                                        & .AND.((SOLID_EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE).OR. &
                                        & (SOLID_EQUATIONS_SET%SPECIFICATION(3)== &
                                        & EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE))) THEN
                                        SolidEquationsSetFound=.TRUE.
                                      ELSE
                                        EquationsSetIndex=EquationsSetIndex+1
                                      END IF
                                    ELSE
                                      CALL FlagError("Solid equations set is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Solid equations not associated.",Err,Error,*999)
                                  END IF
                                END DO
                                IF(SolidEquationsSetFound.EQV..FALSE.) THEN
                                  LOCAL_ERROR="Solid equations set not found when trying to update boundary conditions."
                                  CALL FlagError(LOCAL_ERROR,Err,Error,*999)
                                END IF
                              ELSE
                                CALL FlagError("Solid solver mapping is not associated.",Err,Error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Solver equations for solid equations set not associated.",Err,Error,*999)
                            END IF
                            SOLID_DEPENDENT=>SOLID_EQUATIONS_SET%DEPENDENT
                            IF(.NOT.ASSOCIATED(SOLID_DEPENDENT%DEPENDENT_FIELD)) THEN
                              CALL FlagError("Solid equations set dependent field is not associated.",Err,Error,*999)
                            END IF
                            SOLID_DEPENDENT_FIELD=>SOLID_DEPENDENT%DEPENDENT_FIELD
                            !Find the NavierStokes equations set as there is a FiniteElasticity equations set too
                            FLUID_SOLVER_EQUATIONS=>Solver2%SOLVER_EQUATIONS
                            IF (ASSOCIATED(FLUID_SOLVER_EQUATIONS)) THEN
                              FLUID_SOLVER_MAPPING=>FLUID_SOLVER_EQUATIONS%SOLVER_MAPPING
                              IF (ASSOCIATED(FLUID_SOLVER_MAPPING)) THEN
                                EquationsSetIndex=1
                                FluidEquationsSetFound=.FALSE.
                                DO WHILE (.NOT.FluidEquationsSetFound &
                                  & .AND.EquationsSetIndex<=FLUID_SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                                  FLUID_EQUATIONS=>FLUID_SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(EquationsSetIndex)%EQUATIONS
                                  IF (ASSOCIATED(SOLID_EQUATIONS)) THEN
                                    FLUID_EQUATIONS_SET=>FLUID_EQUATIONS%EQUATIONS_SET
                                    IF(ASSOCIATED(FLUID_EQUATIONS_SET)) THEN
                                      IF(FLUID_EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                                        & .AND.FLUID_EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                                        & .AND.FLUID_EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE) THEN
                                        FluidEquationsSetFound=.TRUE.
                                      ELSE
                                        EquationsSetIndex=EquationsSetIndex+1
                                      END IF
                                    ELSE
                                      CALL FlagError("Fluid equations set is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Fluid equations not associated.",Err,Error,*999)
                                  END IF
                                END DO
                                IF(FluidEquationsSetFound.EQV..FALSE.) THEN
                                  LOCAL_ERROR="Fluid equations set not found when trying to update boundary conditions."
                                  CALL FlagError(LOCAL_ERROR,Err,Error,*999)
                                END IF
                              ELSE
                                CALL FlagError("Fluid solver mapping is not associated.",Err,Error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Fluid equations for fluid equations set not associated.",Err,Error,*999)
                            END IF
                            FLUID_GEOMETRIC=>FLUID_EQUATIONS_SET%GEOMETRY
                            IF(.NOT.ASSOCIATED(FLUID_GEOMETRIC%GEOMETRIC_FIELD)) THEN
                              CALL FlagError("Fluid equations set geometric field is not associated",Err,Error,*999)
                            END IF
                            FLUID_GEOMETRIC_FIELD=>FLUID_GEOMETRIC%GEOMETRIC_FIELD
                            !DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_idx=1
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                            IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF (ASSOCIATED(DOMAIN)) THEN
                                  IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% & 
                                            & CONDITION_TYPES(local_ny)
                                          !Update moved wall nodes only
                                          IF (BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            !NOTE: assuming same mesh and mesh nodes for fluid domain and moving mesh domain
                                            FluidNodeNumber=node_idx
                                            DO search_idx=1,SIZE(Solver2%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                                              & INTERFACE_CONDITIONS(1)%PTR%INTERFACE% &
                                              & NODES%COUPLED_NODES(2,:))
                                              IF (Solver2%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                                                & INTERFACE_CONDITIONS(1)%PTR%INTERFACE% &
                                                & NODES%COUPLED_NODES(2,search_idx)==node_idx) THEN
                                                SolidNodeNumber=Solver2%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                                                  & INTERFACE_CONDITIONS(1)%PTR%INTERFACE% &
                                                  & NODES%COUPLED_NODES(1,search_idx)!might wanna put a break here
                                                SolidNodeFound=.TRUE.
                                              END IF
                                            END DO
                                            IF(.NOT.SolidNodeFound &
                                              & .OR.FluidNodeNumber==0) CALL FlagError("Solid interface node not found.", &
                                              & Err,Error,*999)
                                            !Default to version number 1
                                            IF (variable_idx==1) THEN
                                              CALL FIELD_PARAMETER_SET_GET_NODE(FLUID_GEOMETRIC_FIELD,variable_type, &
                                                & FIELD_VALUES_SET_TYPE,1,deriv_idx, &
                                                & FluidNodeNumber,component_idx,FluidGFValue,Err,Error,*999)
                                            ELSE
                                              FluidGFValue=0.0_DP
                                            END IF
                                            CALL FIELD_PARAMETER_SET_GET_NODE(SOLID_DEPENDENT_FIELD,variable_type, &
                                              & FIELD_VALUES_SET_TYPE,1,deriv_idx, &
                                              & SolidNodeNumber,component_idx,SolidDFValue,Err,Error,*999)
                                            NewLaplaceBoundaryValue=SolidDFValue-FluidGFValue
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & NewLaplaceBoundaryValue,Err,Error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                            !END DO !variable_idx
                          ELSE
                            CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                          END IF
                        ELSE
                          CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  !Pre solve for the dynamic solver
                ELSE IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Velocity field change boundary conditions... ",ERR,ERROR,*999)
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                    !Find the NavierStokes equations set as there is a finite elasticity equations set too 
                    EquationsSetIndex=1
                    ALENavierStokesEquationsSetFound=.FALSE.
                    DO WHILE (.NOT.ALENavierStokesEquationsSetFound &
                      & .AND.EquationsSetIndex<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                      EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(EquationsSetIndex)%EQUATIONS
                      IF (ASSOCIATED(EQUATIONS)) THEN
                        EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                        IF(ASSOCIATED(EQUATIONS_SET)) THEN
                          IF(EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                            & .AND.EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                            & .AND.EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE) THEN
                            ALENavierStokesEquationsSetFound=.TRUE.
                          ELSE
                            EquationsSetIndex=EquationsSetIndex+1
                          END IF
                        ELSE
                          CALL FlagError("ALE Navier-Stokes equations set is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("ALE equations not associated.",Err,Error,*999)
                      END IF
                    END DO
                    IF(ALENavierStokesEquationsSetFound.EQV..FALSE.) THEN
                      LOCAL_ERROR="ALE NavierStokes equations set not found when trying to update boundary conditions."
                      CALL FlagError(LOCAL_ERROR,Err,Error,*999)
                    END IF
                    !Get boundary conditions
                    BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                    IF (ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                      FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                      IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                          & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("Field U variable is not associated",ERR,ERROR,*999)
                      END IF
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                        NULLIFY(MESH_VELOCITY_VALUES)
                        CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                          & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,ERR,ERROR,*999)
                        NULLIFY(BOUNDARY_VALUES)
                        CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                          & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                        !Get update for time-dependent boundary conditions
                        IF (CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER==1) THEN
                          componentBC=1
                          CALL FluidMechanics_IO_UpdateBoundaryConditionUpdateNodes(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                            & SOLVER%SOLVE_TYPE,InletNodes, &
                            & BoundaryValues,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CURRENT_TIME,CONTROL_LOOP%TIME_LOOP%STOP_TIME)
                          DO node_idx=1,SIZE(InletNodes)
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,InletNodes(node_idx),componentBC, &
                              & BoundaryValues(node_idx),ERR,ERROR,*999)
                          END DO
                        ELSE
                          !Figure out which component we're applying BC at
                          IF (CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER==2) THEN
                            componentBC=1
                          ELSE
                            componentBC=2
                          END IF
                          !Get inlet nodes and the corresponding velocities
                          CALL FluidMechanics_IO_UpdateBoundaryConditionUpdateNodes(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                            & SOLVER%SOLVE_TYPE,InletNodes, &
                            & BoundaryValues,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CURRENT_TIME,CONTROL_LOOP%TIME_LOOP%STOP_TIME)
                          DO node_idx=1,SIZE(InletNodes)
                            CALL FIELD_PARAMETER_SET_UPDATE_NODE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,InletNodes(node_idx),componentBC, &
                              & BoundaryValues(node_idx),ERR,ERROR,*999)
                          END DO
                        END IF
                        CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                          & FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,ERR,ERROR,*999)
                        CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                          & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("Boundary condition variable is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                ENDIF
                ! do nothing ???
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for a FiniteElasticity-NavierStokes problem type of a multi physics problem class."
                CALL FlagError(LOCAL_ERROR,Err,Error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
                & " is not valid for NAVIER_STOKES_PRE_SOLVE of a multi physics problem class."
              CALL FlagError(LOCAL_ERROR,Err,Error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The first problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%specification(1),"*",ERR,ERROR))// &
              & " is not valid for NavierStokes_PreSolveUpdateBoundaryConditions."
            CALL FlagError(LOCAL_ERROR,Err,Error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    END IF

    EXITS("NavierStokes_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORS("NavierStokes_PreSolveUpdateBoundaryConditions",ERR,ERROR)
    EXITS("NavierStokes_PreSolveUpdateBoundaryConditions")
    RETURN 1
    
  END SUBROUTINE NavierStokes_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !

  !>Update mesh velocity and move mesh for ALE Navier-Stokes problem
  SUBROUTINE NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_LAPLACE, EQUATIONS_SET_ALE_NAVIER_STOKES
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_LAPLACE, INDEPENDENT_FIELD_ALE_NAVIER_STOKES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_LAPLACE, SOLVER_EQUATIONS_ALE_NAVIER_STOKES
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_LAPLACE, SOLVER_MAPPING_ALE_NAVIER_STOKES
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_NAVIER_STOKES, SOLVER_LAPLACE
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    INTEGER(INTG) :: I,NUMBER_OF_DIMENSIONS_LAPLACE,NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,INPUT_TYPE,INPUT_OPTION,EquationsSetIndex
    INTEGER(INTG) :: component_idx,deriv_idx,local_ny,node_idx,variable_idx,variable_type
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    LOGICAL :: ALENavierStokesEquationsSetFound=.FALSE.
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH",ERR,ERROR,*999)

    IF (ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
        NULLIFY(SOLVER_LAPLACE)
        NULLIFY(SOLVER_ALE_NAVIER_STOKES)
        NULLIFY(INDEPENDENT_FIELD_ALE_NAVIER_STOKES)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE,PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
              !Update mesh within the dynamic solver
              IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                !Get the independent field for the ALE Navier-Stokes problem
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                SOLVER_EQUATIONS_ALE_NAVIER_STOKES=>SOLVER_ALE_NAVIER_STOKES%SOLVER_EQUATIONS
                IF (ASSOCIATED(SOLVER_EQUATIONS_ALE_NAVIER_STOKES)) THEN
                  SOLVER_MAPPING_ALE_NAVIER_STOKES=>SOLVER_EQUATIONS_ALE_NAVIER_STOKES%SOLVER_MAPPING
                  IF (ASSOCIATED(SOLVER_MAPPING_ALE_NAVIER_STOKES)) THEN
                    EQUATIONS_SET_ALE_NAVIER_STOKES=>SOLVER_MAPPING_ALE_NAVIER_STOKES%EQUATIONS_SETS(1)%PTR
                    IF (ASSOCIATED(EQUATIONS_SET_ALE_NAVIER_STOKES)) THEN
                      INDEPENDENT_FIELD_ALE_NAVIER_STOKES=>EQUATIONS_SET_ALE_NAVIER_STOKES%INDEPENDENT%INDEPENDENT_FIELD
                    ELSE
                      CALL FlagError("ALE Navier-Stokes equations set is not associated.",ERR,ERROR,*999)
                    END IF
                    !Get the data
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                    !\todo: Introduce user calls instead of hard-coding 42/1
                    !Copy input to Navier-Stokes' independent field
                    INPUT_TYPE=42
                    INPUT_OPTION=1
                    NULLIFY(MESH_DISPLACEMENT_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, & 
                      & NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES,INPUT_TYPE,INPUT_OPTION, &
                      & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_NAVIER_STOKES%INDEPENDENT%INDEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_NAVIER_STOKES%INDEPENDENT%INDEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("ALE Navier-Stokes solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("ALE Navier-Stokes solver equations are not associated.",ERR,ERROR,*999)
                END IF
                !Use calculated values to update mesh
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                !CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                !  & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                EQUATIONS=>SOLVER_MAPPING_ALE_NAVIER_STOKES%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF (ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                  IF (ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    DO variable_idx=1,EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%NUMBER_OF_VARIABLES
                      variable_type=EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                      FIELD_VARIABLE=>EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                          DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                          IF(ASSOCIATED(DOMAIN)) THEN
                            IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                              IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                !Loop over the local nodes excluding the ghosts.
                                DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                  DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                    !Default to version 1 of each node derivative
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY% &
                                      & GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                      & MESH_DISPLACEMENT_VALUES(local_ny),ERR,ERROR,*999)
                                  END DO !deriv_idx
                                END DO !node_idx
                              END IF
                            END IF
                          END IF
                        END DO !compIdx
                      END IF
                    END DO !variable_idx
                  ELSE
                    CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                !Now use displacement values to calculate velocity values
                TIME_INCREMENT=CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT
                ALPHA=1.0_DP/TIME_INCREMENT
                CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,ERR,ERROR,*999)
              ELSE  
                CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Update mesh within the dynamic solver
              IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                IF (SOLVER%DYNAMIC_SOLVER%ALE) THEN
                  !Get the dependent field for the three component Laplace problem
                  CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_LAPLACE,ERR,ERROR,*999)
                  SOLVER_EQUATIONS_LAPLACE=>SOLVER_LAPLACE%SOLVER_EQUATIONS
                  IF (ASSOCIATED(SOLVER_EQUATIONS_LAPLACE)) THEN
                    SOLVER_MAPPING_LAPLACE=>SOLVER_EQUATIONS_LAPLACE%SOLVER_MAPPING
                    IF (ASSOCIATED(SOLVER_MAPPING_LAPLACE)) THEN
                      EQUATIONS_SET_LAPLACE=>SOLVER_MAPPING_LAPLACE%EQUATIONS_SETS(1)%PTR
                      IF (ASSOCIATED(EQUATIONS_SET_LAPLACE)) THEN
                        DEPENDENT_FIELD_LAPLACE=>EQUATIONS_SET_LAPLACE%DEPENDENT%DEPENDENT_FIELD
                      ELSE
                        CALL FlagError("Laplace equations set is not associated.",ERR,ERROR,*999)
                      END IF
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_LAPLACE%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & NUMBER_OF_DIMENSIONS_LAPLACE,ERR,ERROR,*999)
                    ELSE
                      CALL FlagError("Laplace solver mapping is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Laplace solver equations are not associated.",ERR,ERROR,*999)
                  END IF
                  !Get the independent field for the ALE Navier-Stokes problem
                  CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                  SOLVER_EQUATIONS_ALE_NAVIER_STOKES=>SOLVER_ALE_NAVIER_STOKES%SOLVER_EQUATIONS
                  IF (ASSOCIATED(SOLVER_EQUATIONS_ALE_NAVIER_STOKES)) THEN
                    SOLVER_MAPPING_ALE_NAVIER_STOKES=>SOLVER_EQUATIONS_ALE_NAVIER_STOKES%SOLVER_MAPPING
                    IF (ASSOCIATED(SOLVER_MAPPING_ALE_NAVIER_STOKES)) THEN
                      EQUATIONS_SET_ALE_NAVIER_STOKES=>SOLVER_MAPPING_ALE_NAVIER_STOKES%EQUATIONS_SETS(1)%PTR
                      IF (ASSOCIATED(EQUATIONS_SET_ALE_NAVIER_STOKES)) THEN
                        INDEPENDENT_FIELD_ALE_NAVIER_STOKES=>EQUATIONS_SET_ALE_NAVIER_STOKES%INDEPENDENT%INDEPENDENT_FIELD
                      ELSE
                        CALL FlagError("ALE Navier-Stokes equations set is not associated.",ERR,ERROR,*999)
                      END IF
                      CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                        & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                    ELSE
                      CALL FlagError("ALE Navier-Stokes solver mapping is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("ALE Navier-Stokes solver equations are not associated.",ERR,ERROR,*999)
                  END IF
                  !Copy result from Laplace mesh movement to Navier-Stokes' independent field
                  IF (NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES==NUMBER_OF_DIMENSIONS_LAPLACE) THEN
                    DO I=1,NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES
                      CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_LAPLACE, & 
                        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,INDEPENDENT_FIELD_ALE_NAVIER_STOKES, & 
                        & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,I,ERR,ERROR,*999)
                    ENDDO
                  ELSE
                    CALL FlagError("Dimension of Laplace and ALE Navier-Stokes equations set is not consistent.",ERR,ERROR,*999)
                  END IF
                  !Use calculated values to update mesh
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  NULLIFY(MESH_DISPLACEMENT_VALUES)
                  CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                  EQUATIONS=>SOLVER_MAPPING_LAPLACE%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF (ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                    IF (ASSOCIATED(EQUATIONS_MAPPING)) THEN
                      DO variable_idx=1,EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%NUMBER_OF_VARIABLES
                        variable_type=EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                        FIELD_VARIABLE=>EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD% &
                          & VARIABLE_TYPE_MAP(variable_type)%PTR
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                            DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                            IF(ASSOCIATED(DOMAIN)) THEN
                              IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                  !Loop over the local nodes excluding the ghosts.
                                  DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                    DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                      !Default to version 1 of each node derivative
                                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                      CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY% &
                                        & GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                        & MESH_DISPLACEMENT_VALUES(local_ny),ERR,ERROR,*999)
                                    END DO !deriv_idx
                                  END DO !node_idx
                                END IF
                              END IF
                            END IF
                          END DO !compIdx
                        END IF
                      END DO !variable_idx
                    ELSE
                      CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*999)
                    END IF
                    CALL FIELD_PARAMETER_SET_DATA_RESTORE(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                      & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                  ELSE
                    CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                  END IF
                  CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  !Now use displacement values to calculate velocity values
                  TIME_INCREMENT=CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT
                  ALPHA=1.0_DP/TIME_INCREMENT
                  CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,ERR,ERROR,*999)
                ELSE  
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
                END IF
              ELSE  
                CALL FlagError("Mesh update is not defined for non-dynamic problems.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
                !Update mesh within the dynamic solver
                IF (SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  IF (SOLVER%DYNAMIC_SOLVER%ALE) THEN
                    !Get the dependent field for the Laplace problem
                    !     CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_LAPLACE,ERR,ERROR,*999)
                    CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_LAPLACE,ERR,ERROR,*999)
                    SOLVER_EQUATIONS_LAPLACE=>SOLVER_LAPLACE%SOLVER_EQUATIONS
                    IF (ASSOCIATED(SOLVER_EQUATIONS_LAPLACE)) THEN
                      SOLVER_MAPPING_LAPLACE=>SOLVER_EQUATIONS_LAPLACE%SOLVER_MAPPING
                      IF (ASSOCIATED(SOLVER_MAPPING_LAPLACE)) THEN
                        EQUATIONS_SET_LAPLACE=>SOLVER_MAPPING_LAPLACE%EQUATIONS_SETS(1)%PTR
                        IF (ASSOCIATED(EQUATIONS_SET_LAPLACE)) THEN
                          DEPENDENT_FIELD_LAPLACE=>EQUATIONS_SET_LAPLACE%DEPENDENT%DEPENDENT_FIELD
                        ELSE
                          CALL FlagError("Laplace equations set is not associated.",ERR,ERROR,*999)
                        END IF
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_LAPLACE%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS_LAPLACE,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("Laplace solver mapping is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Laplace solver equations are not associated.",ERR,ERROR,*999)
                    END IF
                    !Get the independent field for the ALE Navier-Stokes problem
                    !    CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                    CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                    SOLVER_EQUATIONS_ALE_NAVIER_STOKES=>SOLVER_ALE_NAVIER_STOKES%SOLVER_EQUATIONS
                    IF (ASSOCIATED(SOLVER_EQUATIONS_ALE_NAVIER_STOKES)) THEN
                      SOLVER_MAPPING_ALE_NAVIER_STOKES=>SOLVER_EQUATIONS_ALE_NAVIER_STOKES%SOLVER_MAPPING
                      IF (ASSOCIATED(SOLVER_MAPPING_ALE_NAVIER_STOKES)) THEN
                        EquationsSetIndex=1
                        ALENavierStokesEquationsSetFound=.FALSE.
                        !Find the NavierStokes equations set as there is a finite elasticity equations set too
                        DO WHILE (.NOT.ALENavierStokesEquationsSetFound &
                          & .AND.EquationsSetIndex<=SOLVER_MAPPING_ALE_NAVIER_STOKES%NUMBER_OF_EQUATIONS_SETS)
                          EQUATIONS_SET_ALE_NAVIER_STOKES=>SOLVER_MAPPING_ALE_NAVIER_STOKES%EQUATIONS_SETS(EquationsSetIndex)%PTR
                          IF(ASSOCIATED(EQUATIONS_SET_ALE_NAVIER_STOKES)) THEN
                            IF(EQUATIONS_SET_ALE_NAVIER_STOKES%SPECIFICATION(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                              & .AND.EQUATIONS_SET_ALE_NAVIER_STOKES%SPECIFICATION(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                              & .AND.EQUATIONS_SET_ALE_NAVIER_STOKES%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE) THEN
                              INDEPENDENT_FIELD_ALE_NAVIER_STOKES=>EQUATIONS_SET_ALE_NAVIER_STOKES%INDEPENDENT%INDEPENDENT_FIELD
                              IF (ASSOCIATED(INDEPENDENT_FIELD_ALE_NAVIER_STOKES)) ALENavierStokesEquationsSetFound=.TRUE.
                            ELSE
                              EquationsSetIndex=EquationsSetIndex+1
                            END IF
                          ELSE
                            CALL FlagError("ALE Navier-Stokes equations set is not associated.",ERR,ERROR,*999)
                          END IF
                        END DO
                        IF(ALENavierStokesEquationsSetFound.EQV..FALSE.) THEN
                          CALL FlagError("ALE NavierStokes equations set not found when trying to update ALE mesh.",Err,Error,*999)
                        END IF
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES,ERR,ERROR,*999)
                      ELSE
                        CALL FlagError("ALE Navier-Stokes solver mapping is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("ALE Navier-Stokes solver equations are not associated.",ERR,ERROR,*999)
                    END IF
                    !Copy result from Laplace mesh movement to Navier-Stokes' independent field
                    IF (NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES==NUMBER_OF_DIMENSIONS_LAPLACE) THEN
                      DO I=1,NUMBER_OF_DIMENSIONS_ALE_NAVIER_STOKES
                        CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_LAPLACE, & 
                          & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,INDEPENDENT_FIELD_ALE_NAVIER_STOKES, & 
                          & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,I,ERR,ERROR,*999)
                      ENDDO
                    ELSE
                      CALL FlagError("Dimension of Laplace and ALE Navier-Stokes equations set is not consistent.",ERR,ERROR,*999)
                    END IF
                    !Use calculated values to update mesh
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    NULLIFY(MESH_DISPLACEMENT_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                      & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                    EQUATIONS=>SOLVER_MAPPING_LAPLACE%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF (ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                      IF (ASSOCIATED(EQUATIONS_MAPPING)) THEN
                        DO variable_idx=1,EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%NUMBER_OF_VARIABLES
                          variable_type=EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(variable_idx)% &
                            & VARIABLE_TYPE
                          FIELD_VARIABLE=>EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD% &
                            & VARIABLE_TYPE_MAP(variable_type)%PTR
                          IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                            DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                              DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                              IF (ASSOCIATED(DOMAIN)) THEN
                                IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                  DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                  IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                    !Loop over the local nodes excluding the ghosts.
                                    DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                      DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                        !Default to version 1 of each node derivative
                                        local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                          & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                        CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY% &
                                          & GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                          & MESH_DISPLACEMENT_VALUES(local_ny),ERR,ERROR,*999)
                                      END DO !deriv_idx
                                    END DO !node_idx
                                  END IF
                                END IF
                              END IF
                            END DO !component_idx
                          END IF
                        END DO !variable_idx
                      ELSE
                        CALL FlagError("Equations mapping is not associated.",ERR,ERROR,*999)
                      END IF
                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                    ELSE
                      CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                    END IF
                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_ALE_NAVIER_STOKES%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                    !Now use displacement values to calculate velocity values
                    TIME_INCREMENT=CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT
                    ALPHA=1.0_DP/TIME_INCREMENT
                    CALL FIELD_PARAMETER_SETS_COPY(INDEPENDENT_FIELD_ALE_NAVIER_STOKES,FIELD_U_VARIABLE_TYPE, & 
                      & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,ERR,ERROR,*999)
                  ELSE  
                    CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
                  END IF
                ELSE  
                  CALL FlagError("Mesh update is not defined for non-dynamic problems.",ERR,ERROR,*999)
                END IF
              CASE DEFAULT
                LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
                & " is not valid for NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH of a multi physics problem class."
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem class "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(1),"*",ERR,ERROR))// &
              & " is not valid for NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NAVIER_STOKES_PRE_SOLVE_ALE_UPDATE_MESH

  !
  !================================================================================================================================
  !
  !>Update mesh parameters for Laplace problem
  SUBROUTINE NavierStokes_PreSolveALEUpdateParameters(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: component_idx,deriv_idx,local_ny,node_idx,variable_idx,variable_type
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP), POINTER :: MESH_STIFF_VALUES(:)

    ENTERS("NavierStokes_PreSolveALEUpdateParameters",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              !do nothing ???
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE,PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              IF (SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                !Get the independent field for the ALE Navier-Stokes problem
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF (ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                    NULLIFY(MESH_STIFF_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,ERR,ERROR,*999)
                    IF (ASSOCIATED(EQUATIONS_SET)) THEN
                      EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                      IF (ASSOCIATED(EQUATIONS)) THEN
                        INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                        IF (ASSOCIATED(INDEPENDENT_FIELD)) THEN
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          !Calculation of K values dependent on current mesh topology
                                          MESH_STIFF_VALUES(local_ny)=1.0_DP
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
                                            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                            & MESH_STIFF_VALUES(local_ny),ERR,ERROR,*999)
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                          END DO !variable_idx
                        ELSE
                          CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                    END IF
                    CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,ERR,ERROR,*999)                     
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                END IF
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE)
                IF (SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                  !Get the independent field for the ALE Navier-Stokes problem
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                    IF (ASSOCIATED(SOLVER_MAPPING)) THEN
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                      NULLIFY(MESH_STIFF_VALUES)
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,ERR,ERROR,*999)
                      IF (ASSOCIATED(EQUATIONS_SET)) THEN
                        EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                        IF (ASSOCIATED(EQUATIONS)) THEN
                          INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                          IF (ASSOCIATED(INDEPENDENT_FIELD)) THEN
                            DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                              FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                              IF (ASSOCIATED(FIELD_VARIABLE)) THEN
                                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                  IF (ASSOCIATED(DOMAIN)) THEN
                                    IF (ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                      DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                      IF (ASSOCIATED(DOMAIN_NODES)) THEN
                                        !Loop over the local nodes excluding the ghosts.
                                        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                          DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                            !Default to version 1 of each node derivative
                                            local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                              & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                            !Calculation of K values dependent on current mesh topology
                                            MESH_STIFF_VALUES(local_ny)=1.0_DP
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT% &
                                              & INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, & 
                                              & MESH_STIFF_VALUES(local_ny),ERR,ERROR,*999)
                                          END DO !deriv_idx
                                        END DO !node_idx
                                      END IF
                                    END IF
                                  END IF
                                END DO !component_idx
                              END IF
                            END DO !variable_idx
                          ELSE
                            CALL FlagError("Independent field is not associated.",ERR,ERROR,*999)
                          END IF
                        ELSE
                          CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                      END IF
                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                        & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,ERR,ERROR,*999)                     
                    ELSE
                      CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
                  END IF
                ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",ERR,ERROR,*999)
                END IF
              CASE DEFAULT
                LOCAL_ERROR="The third problem specification of "// &
                  & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                  & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The second problem specification of "// &
                & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%specification(2),"*",ERR,ERROR))// &
                & " is not valid for NavierStokes_PreSolveALEUpdateParameters of a multi physics problem."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The first problem specification of "// &
              & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%specification(1),"*",ERR,ERROR))// &
              & " is not valid for NavierStokes_PreSolveALEUpdateParameters."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NavierStokes_PreSolveALEUpdateParameters")
    RETURN
999 ERRORSEXITS("NavierStokes_PreSolveALEUpdateParameters",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NavierStokes_PreSolveALEUpdateParameters

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELDS_TYPE), POINTER :: Fields
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR,METHOD,VFileName,FILENAME
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,FileNameLength
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,START_TIME,STOP_TIME
    LOGICAL :: EXPORT_FIELD
    CHARACTER(14) :: FILE,OUTPUT_FILE

    NULLIFY(Fields) 

    ENTERS("NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF (ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))

          CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)

            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    METHOD="FORTRAN"
                    EXPORT_FIELD=.TRUE.
                    IF(EXPORT_FIELD) THEN          
                      OUTPUT_FILE = "StaticSolution"
                      FileNameLength = LEN_TRIM(OUTPUT_FILE)
                      VFileName = OUTPUT_FILE(1:FileNameLength)
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                      Fields=>EQUATIONS_SET%REGION%FIELDS
                      CALL FIELD_IO_NODES_EXPORT(Fields,VFileName,METHOD,ERR,ERROR,*999)
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export elements... ",ERR,ERROR,*999)
                      CALL FIELD_IO_ELEMENTS_EXPORT(Fields,VFileName,METHOD,ERR,ERROR,*999)
                      NULLIFY(Fields)
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                    END IF
                  END DO
                END IF 
              END IF

          CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)

            CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF (ASSOCIATED(SOLVER_MAPPING)) THEN
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                  OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER
                  IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                    IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
                      WRITE(OUTPUT_FILE,'("TimeStep_",I0)') CURRENT_LOOP_ITERATION
                      FILE=OUTPUT_FILE
                      METHOD="FORTRAN"
                      EXPORT_FIELD=.TRUE.
                      IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                        !Use standard field IO routines (also only export nodes after first step as not a moving mesh case)
                        FileNameLength=LEN_TRIM(OUTPUT_FILE)
                        VFileName=OUTPUT_FILE(1:FileNameLength)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                        Fields=>EQUATIONS_SET%REGION%FIELDS
                        CALL FIELD_IO_NODES_EXPORT(Fields,VFileName,METHOD,ERR,ERROR,*999)
                        IF(CURRENT_LOOP_ITERATION==0) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export elements... ",ERR,ERROR,*999)
                          CALL FIELD_IO_ELEMENTS_EXPORT(Fields,VFileName,METHOD,ERR,ERROR,*999)
                        END IF
                        NULLIFY(Fields)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                      END IF                          
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                          CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,ERR,ERROR,*999)
                        END IF
                      END IF
                    END IF 
                  END IF
                END DO
              END IF
            END IF

          CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)

            CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT, &
              & CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,ERR,ERROR,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF (ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF (ASSOCIATED(SOLVER_MAPPING)) THEN
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                    IF(CURRENT_TIME<=STOP_TIME) THEN
                      IF(CURRENT_LOOP_ITERATION<10) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF (CURRENT_LOOP_ITERATION<100) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF (CURRENT_LOOP_ITERATION<1000) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF (CURRENT_LOOP_ITERATION<10000) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                      ENDIF
                      DEPENDENT_REGION=>EQUATIONS_SET%REGION
                      FILE=OUTPUT_FILE
                      FILENAME="./output/"//"MainTime_"//TRIM(NUMBER_TO_VSTRING(CURRENT_LOOP_ITERATION,"*",ERR,ERROR))
                      METHOD="FORTRAN"
                      EXPORT_FIELD=.TRUE.
                      IF (EXPORT_FIELD) THEN          
                        IF (MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
                          CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                          CALL FIELD_IO_ELEMENTS_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILENAME,ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                        END IF
                      END IF 
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                          CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,ERR,ERROR,*999)
                        END IF
                      END IF
                    END IF 
                  END IF
                END DO
              END IF
            END IF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",ERR,ERROR,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Sets up analytic parameters and calls NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE to evaluate solutions to analytic problems
  SUBROUTINE NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoint(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_TYPE), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable,analyticVariable,materialsVariable
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: compIdx,derivIdx,dimensionIdx,local_ny,nodeIdx,numberOfDimensions,variableIdx,variableType,I,J,K
    INTEGER(INTG) :: numberOfNodesXiCoord(3),elementIdx,en_idx,boundaryCount,analyticFunctionType,globalDerivativeIndex,versionIdx
    INTEGER(INTG) :: boundaryConditionsCheckVariable,numberOfXi,nodeNumber,userNodeNumber,localDof,globalDof
    INTEGER(INTG) :: parameterIdx,numberOfParameters
    REAL(DP) :: TIME,VALUE,X(3),xiCoordinates(3),initialValue,T_COORDINATES(20,3),nodeAnalyticParameters(10)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)

    ENTERS("NavierStokes_BoundaryConditionsAnalyticCalculate",err,error,*999)

    boundaryCount=0
    xiCoordinates(3)=0.0_DP

    IF (ASSOCIATED(equationsSet)) THEN
      IF (ASSOCIATED(equationsSet%ANALYTIC)) THEN
        dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
        IF (ASSOCIATED(dependentField)) THEN
          geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
          IF (ASSOCIATED(geometricField)) THEN     
            !Geometric parameters
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            NULLIFY(geometricVariable)
            CALL FIELD_VARIABLE_GET(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
            NULLIFY(geometricParameters)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,geometricParameters, &
              & err,error,*999)
            !Analytic parameters
            analyticFunctionType=equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE
            analyticField=>equationsSet%ANALYTIC%ANALYTIC_FIELD
            NULLIFY(analyticVariable)
            NULLIFY(analyticParameters)
            IF (ASSOCIATED(analyticField)) THEN
              CALL FIELD_VARIABLE_GET(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & analyticParameters,err,error,*999)           
            END IF
            !Materials parameters
            NULLIFY(materialsField)
            NULLIFY(materialsVariable)
            NULLIFY(materialsParameters)
            IF (ASSOCIATED(equationsSet%MATERIALS)) THEN
              materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
              CALL FIELD_VARIABLE_GET(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & materialsParameters,err,error,*999)           
            END IF
            TIME=equationsSet%ANALYTIC%ANALYTIC_TIME
            !Interpolation parameters
            NULLIFY(interpolationParameters)
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,interpolationParameters,err,error,*999)
            NULLIFY(interpolatedPoint) 
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoint,err,error,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    IF (ASSOCIATED(boundaryConditions)) THEN
      DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
        variableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
        fieldVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%PTR
        IF (ASSOCIATED(fieldVariable)) THEN
          IF (.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_ANALYTIC_VALUES_SET_TYPE)%PTR)) &
            & CALL FIELD_PARAMETER_SET_CREATE(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          DO compIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
            boundaryCount=0
            IF (fieldVariable%COMPONENTS(compIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              domain=>fieldVariable%COMPONENTS(compIdx)%DOMAIN
              IF (ASSOCIATED(domain)) THEN
                IF (ASSOCIATED(domain%TOPOLOGY)) THEN
                  domainNodes=>domain%TOPOLOGY%NODES
                  IF (ASSOCIATED(domainNodes)) THEN
                    !Loop over the local nodes excluding the ghosts.
                    DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                      nodeNumber = domainNodes%NODES(nodeIdx)%local_number
                      userNodeNumber = domainNodes%NODES(nodeIdx)%user_number
                      elementIdx=domain%topology%nodes%nodes(nodeNumber)%surrounding_elements(1)
                      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementIdx, &
                        & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                      en_idx=0
                      xiCoordinates=0.0_DP
                      numberOfXi=domain%topology%elements%elements(elementIdx)%basis%number_of_xi
                      numberOfNodesXiCoord(1)=domain%topology%elements%elements(elementIdx)%basis%number_of_nodes_xic(1)
                      IF (numberOfXi>1) THEN
                        numberOfNodesXiCoord(2)=domain%topology%elements%elements(elementIdx)%basis%number_of_nodes_xic(2)
                      ELSE
                        numberOfNodesXiCoord(2)=1
                      END IF
                      IF(numberOfXi>2) THEN
                        numberOfNodesXiCoord(3)=domain%topology%elements%elements(elementIdx)%basis%number_of_nodes_xic(3)
                      ELSE
                        numberOfNodesXiCoord(3)=1
                      END IF

                      SELECT CASE(analyticFunctionType)
                      ! --- Calculate analytic profile for validation ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                        IF (variableIdx<3) THEN
                          ! Get geometric position info for this node
                          DO dimensionIdx=1,numberOfDimensions
                            local_ny=geometricVariable%COMPONENTS(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                              & NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                            X(dimensionIdx)=geometricParameters(local_ny)
                          END DO !dimensionIdx
                          DO derivIdx=1,domainNodes%NODES(nodeNumber)%NUMBER_OF_DERIVATIVES
                            globalDerivativeIndex=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)% &
                              & GLOBAL_DERIVATIVE_INDEX
                            DO versionIdx=1,domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)%numberOfVersions
                              CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                & globalDerivativeIndex,compIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                                & analyticParameters,materialsParameters,VALUE,err,error,*999)
                              local_ny=fieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                            END DO !versionIdx
                          END DO !derivIdx
                        END IF ! variableIdx < 3

                      ! --- Set velocity boundary conditions with analytic value ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                        ! Get geometric position info for this node
                        DO dimensionIdx=1,numberOfDimensions
                          local_ny=geometricVariable%COMPONENTS(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(nodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                          X(dimensionIdx)=geometricParameters(local_ny)
                        END DO !dimensionIdx
                        !Loop over the derivatives
                        DO derivIdx=1,domainNodes%NODES(nodeNumber)%NUMBER_OF_DERIVATIVES
                          globalDerivativeIndex=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)% &
                            & GLOBAL_DERIVATIVE_INDEX
                          IF(compIdx<=numberOfXi .OR. &
                            &  analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                            DO versionIdx=1,domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)%numberOfVersions
                              ! Get global and local dof indices
                              CALL FIELD_COMPONENT_DOF_GET_USER_NODE(dependentField,variableType,versionIdx,derivIdx, &
                               & userNodeNumber,compIdx,localDof,globalDof,err,error,*999)
                              IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                                CALL FIELD_NUMBER_OF_COMPONENTS_GET(analyticField,FIELD_U_VARIABLE_TYPE, &
                                 & numberOfParameters,err,error,*999)
                                DO parameterIdx=1,numberOfParameters
                                  ! populate nodeAnalyticParameters
                                  CALL Field_ParameterSetGetLocalNode(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, & 
                                   & versionIdx,derivIdx,nodeIdx,parameterIdx,nodeAnalyticParameters(parameterIdx), &
                                   & err,error,*999)
                                END DO
                                CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                  & globalDerivativeIndex,compIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                                  & nodeAnalyticParameters,materialsParameters,VALUE,err,error,*999)
                              ELSE
                                CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                  & globalDerivativeIndex,compIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                                  & analyticParameters,materialsParameters,VALUE,err,error,*999)
                              END IF
                              ! update analytic field values
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,localDof,VALUE,err,error,*999)
                              IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                                IF(domainNodes%NODES(nodeNumber)%BOUNDARY_NODE) THEN
                                  CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable, &
                                   & boundaryConditionsVariable,err,error,*999)
                                  IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                                    boundaryConditionsCheckVariable=boundaryConditionsVariable% &
                                     & CONDITION_TYPES(globalDof)
                                    ! update dependent field values if fixed inlet or pressure BC
                                    IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                      ! Set velocity/flowrate values
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                       & FIELD_VALUES_SET_TYPE,localDof,VALUE,err,error,*999)
                                    ELSE IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_PRESSURE) THEN
                                      ! Set neumann boundary pressure value on pressure nodes
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                                       & FIELD_PRESSURE_VALUES_SET_TYPE,1,1,nodeIdx,compIdx,VALUE,err,error,*999)
                                    END IF
                                  END IF
                                END IF
                              END IF
                            END DO !versionIdx
                          END IF
                        END DO !derivIdx

                      ! --- Set Flow rate boundary conditions with analytic value ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                        ! Get geometric position info for this node
                        DO dimensionIdx=1,numberOfDimensions
                          local_ny=geometricVariable%COMPONENTS(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                          X(dimensionIdx)=geometricParameters(local_ny)
                        END DO !dimensionIdx
                        !Loop over the derivatives
                        DO derivIdx=1,domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                          globalDerivativeIndex=domainNodes%NODES(nodeIdx)%DERIVATIVES(derivIdx)% &
                            & GLOBAL_DERIVATIVE_INDEX
                          IF (compIdx==1 .AND. variableType==FIELD_U_VARIABLE_TYPE) THEN
                            DO versionIdx=1,domainNodes%NODES(nodeIdx)%DERIVATIVES(derivIdx)%numberOfVersions
                              local_ny=fieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
                              IF (domainNodes%NODES(nodeIdx)%BOUNDARY_NODE) THEN
                                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable, &
                                 & boundaryConditionsVariable,err,error,*999)
                                IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                                  boundaryConditionsCheckVariable=boundaryConditionsVariable%CONDITION_TYPES(local_ny)
                                  IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                    CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                      & globalDerivativeIndex,compIdx,numberOfXi,fieldVariable%NUMBER_OF_COMPONENTS, &
                                      & analyticParameters,materialsParameters,VALUE,err,error,*999)
                                    !If we are a boundary node then set the analytic value on the boundary
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                      & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                  ELSE
                                    CALL Field_ParameterSetGetLocalNode(dependentField,variableType,FIELD_VALUES_SET_TYPE, & 
                                      & versionIdx,derivIdx,nodeIdx,compIdx,VALUE,err,error,*999)
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                                  END IF
                                END IF
                              END IF
                            END DO !versionIdx
                          END IF
                        END DO !derivIdx

                      ! --- Legacy unit shape testing types ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
                        !Quad/Hex
                        !\todo: Use boundary flag
                        IF (domain%topology%elements%maximum_number_of_element_parameters==4.AND.numberOfDimensions==2.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==9.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==16.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==8.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==27.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==64) THEN
                          DO K=1,numberOfNodesXiCoord(3)
                            DO J=1,numberOfNodesXiCoord(2)
                              DO I=1,numberOfNodesXiCoord(1)
                                en_idx=en_idx+1
                                IF (domain%topology%elements%elements(elementIdx)%element_nodes(en_idx)==nodeIdx) EXIT
                                xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXiCoord(1)-1))
                              END DO
                                IF(domain%topology%elements%elements(elementIdx)%element_nodes(en_idx)==nodeIdx) EXIT
                                xiCoordinates(1)=0.0_DP
                                xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXiCoord(2)-1))
                            END DO
                            IF(domain%topology%elements%elements(elementIdx)%element_nodes(en_idx)==nodeIdx) EXIT
                            xiCoordinates(1)=0.0_DP
                            xiCoordinates(2)=0.0_DP
                            IF (numberOfNodesXiCoord(3)/=1) THEN
                              xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXiCoord(3)-1))
                            END IF
                          END DO
                          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,xiCoordinates, &
                            & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                        !Tri/Tet
                        !\todo: Use boundary flag
                        ELSE
                          IF (domain%topology%elements%maximum_number_of_element_parameters==3) THEN
                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                          ELSE IF (domain%topology%elements%maximum_number_of_element_parameters==6) THEN
                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                            T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                            T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                            T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                          ELSE IF (domain%topology%elements%maximum_number_of_element_parameters==10.AND. &
                            & numberOfDimensions==2) THEN
                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                            T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                          ELSE IF (domain%topology%elements%maximum_number_of_element_parameters==4) THEN
                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                          ELSE IF (domain%topology%elements%maximum_number_of_element_parameters==10.AND. &
                            & numberOfDimensions==3) THEN 
                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                            T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                            T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                            T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                            T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                          ELSE IF (domain%topology%elements%maximum_number_of_element_parameters==20) THEN
                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                          END IF
                          DO K=1,domain%topology%elements%maximum_number_of_element_parameters
                            IF(domain%topology%elements%elements(elementIdx)%element_nodes(K)==nodeIdx) EXIT
                          END DO
                          IF(numberOfDimensions==2) THEN
                            CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                              & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                          ELSE IF (numberOfDimensions==3) THEN
                            CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                              & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                          END IF
                        END IF
                        X=0.0_DP
                        DO dimensionIdx=1,numberOfDimensions
                          X(dimensionIdx)=interpolatedPoint(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(dimensionIdx,1)
                        END DO !dimensionIdx
                        !Loop over the derivatives
                        DO derivIdx=1,domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                          globalDerivativeIndex=domainNodes%NODES(nodeIdx)%DERIVATIVES(derivIdx)% &
                            & GLOBAL_DERIVATIVE_INDEX
                          CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                            & globalDerivativeIndex,compIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                            & analyticParameters,materialsParameters,VALUE,err,error,*999)
                          DO versionIdx=1,domainNodes%NODES(nodeIdx)%DERIVATIVES(derivIdx)%numberOfVersions
                            local_ny=fieldVariable%COMPONENTS(compIdx)%PARAM_TO_DOF_MAP% &
                              & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                              & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                            IF (variableType==FIELD_U_VARIABLE_TYPE) THEN
                              IF (domainNodes%NODES(nodeIdx)%BOUNDARY_NODE) THEN
                                !If we are a boundary node then set the analytic value on the boundary
                                CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField,variableType, &
                                  & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                ! \todo: This is just a workaround for linear pressure fields in simplex element components
                                IF (compIdx>numberOfDimensions) THEN
                                  IF (domain%topology%elements%maximum_number_of_element_parameters==3) THEN
                                    IF (analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                                      IF (-0.001_DP<X(1).AND.X(1)<0.001_DP.AND.-0.001_DP<X(2).AND.X(2)<0.001_DP.OR. &
                                        &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.-0.001_DP<X(2).AND. &
                                        & X(2)<0.001_DP.OR. &
                                        &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<10.0_DP+0.001_DP.OR. &
                                        &  -0.001_DP<X(1).AND.X(1)<0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<10.0_DP+0.001_DP) THEN
                                          CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField, &
                                            & variableType,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                          boundaryCount=boundaryCount+1
                                      END IF
                                    END IF
                                  ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==4.AND. &
                                    & numberOfDimensions==3) THEN
                                    IF (analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                      IF (-5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+ 0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP) THEN
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField, &
                                          & variableType,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                        boundaryCount=boundaryCount+1
                                      END IF
                                    END IF
                                    ! \todo: This is how it should be if adjacent elements would be working
                                  ELSE IF (boundaryCount==0) THEN
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField,variableType,&
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                    boundaryCount=boundaryCount+1
                                  END IF
                                END IF
                              ELSE
                                !Set the initial condition.
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                  & FIELD_VALUES_SET_TYPE,local_ny,initialValue,err,error,*999)
                              END IF
                            END IF
                          END DO !versionIdx
                        END DO !derivIdx

                      CASE DEFAULT
                        localError="Analytic Function Type "//TRIM(NUMBER_TO_VSTRING(analyticFunctionType,"*",ERR,ERROR))// &
                          & " is not yet implemented for a Navier-Stokes problem."
                        CALL FlagError(localError,ERR,ERROR,*999)
                      END SELECT

                    END DO !nodeIdx
                  ELSE
                    CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Domain topology is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Domain is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
            END IF
          END DO !compIdx
          CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
        ELSE
          CALL FlagError("Field variable is not associated.",err,error,*999)
        END IF
      END DO !variableIdx
      CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & geometricParameters,err,error,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoint,err,error,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
    ELSE
      CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("NavierStokes_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("NavierStokes_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE NavierStokes_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !
  !>Calculates the various analytic values for NSE examples with exact solutions
  SUBROUTINE NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(ANALYTIC_FUNCTION_TYPE,X,TIME,VARIABLE_TYPE,GLOBAL_DERIV_INDEX, &
    & componentNumber,NUMBER_OF_DIMENSIONS,NUMBER_OF_COMPONENTS,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE,VARIABLE_TYPE,GLOBAL_DERIV_INDEX
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_COMPONENTS,componentNumber
    REAL(DP), INTENT(IN) :: X(:),TIME,ANALYTIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    REAL(DP), INTENT(OUT) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,j,n,m
    REAL(DP) :: H_PARAM,MU_PARAM,RHO_PARAM,startTime,stopTime,delta(1000),t(1000),q(1000),s
    REAL(DP) :: componentCoeff(4),L_PARAM,U_PARAM,P_PARAM,NU_PARAM,K_PARAM,Qo,tt,tmax
    REAL(DP) :: amplitude,yOffset,period,phaseShift,frequency,INTERNAL_TIME,CURRENT_TIME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR,*999)

    !\todo: Introduce user-defined or default values instead for density and viscosity
    INTERNAL_TIME=TIME
    CURRENT_TIME=TIME

    SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
      !For fully developed 2D laminar flow through a channel, NSE should yield a parabolic profile, 
      !U = Umax(1-y^2/H^2), Umax = (-dP/dx)*(H^2/(2*MU)), Umax = (3/2)*Umean
      !Note: assumes a flat inlet profile (U_PARAM = Umean).
      !Nonlinear terms from NSE will effectively be 0 for Poiseuille flow 
      IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          L_PARAM = ANALYTIC_PARAMETERS(1) ! channel length in x-direction
          H_PARAM = ANALYTIC_PARAMETERS(2) ! channel height in y-direction
          U_PARAM = ANALYTIC_PARAMETERS(3) ! mean (inlet) velocity
          P_PARAM = ANALYTIC_PARAMETERS(4) ! pressure value at outlet 
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=(3.0_DP/2.0_DP)*U_PARAM*(1.0_DP-((X(2)-H_PARAM)**2)/(H_PARAM**2))
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=0.0_DP
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE = (3.0_DP*MU_PARAM*U_PARAM*(X(1)-L_PARAM))/(H_PARAM**2)+P_PARAM
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE( NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
      !Exact solution to 2D laminar, dynamic, nonlinear Taylor-Green vortex decay  
      IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        NU_PARAM = MU_PARAM/RHO_PARAM ! kinematic viscosity
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          U_PARAM = ANALYTIC_PARAMETERS(1) ! characteristic velocity (initial amplitude)
          L_PARAM = ANALYTIC_PARAMETERS(2) ! length scale for square
          K_PARAM = 2.0_DP*PI/L_PARAM   ! scale factor for equations
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=-1.0_DP*U_PARAM*COS(K_PARAM*X(1))*SIN(K_PARAM*X(2))*EXP(-2.0_DP*(K_PARAM**2)*NU_PARAM*CURRENT_TIME)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=U_PARAM*SIN(K_PARAM*X(1))*COS(K_PARAM*X(2))*EXP(-2.0_DP*(K_PARAM**2)*NU_PARAM*CURRENT_TIME)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE =-1.0_DP*(U_PARAM**2)*(RHO_PARAM/4.0_DP)*(COS(2.0_DP*K_PARAM*X(1))+ &
                & COS(2.0_DP*K_PARAM*X(2)))*(EXP(-4.0_DP*(K_PARAM**2)*NU_PARAM*CURRENT_TIME))
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE( NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA)
      SELECT CASE(NUMBER_OF_DIMENSIONS)
      CASE(1)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !Input function
              period = 800
              tt=MOD(TIME,period)
              tmax=150.0_DP
              Qo=100000.0_DP
              VALUE=1.0_DP*(Qo*tt/(tmax**2.0_DP))*EXP(-(tt**2.0_DP)/(2.0_DP*(tmax**2.0_DP)))
            ELSE
              CALL FlagError("Incorrect component specification for Aorta flow rate waveform ",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE) 
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE)
          ! Do nothing
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Aorta flowrate waveform for "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
          & " dimension problem has not yet been implemented."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
      SELECT CASE(NUMBER_OF_DIMENSIONS)
      CASE(1)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !Olufsen Aorta
!              t(1)= 0.0011660 ; q(1)= 17.39051
!              t(2)= 0.0215840 ; q(2)= 10.41978
!              t(3)= 0.0340860 ; q(3)= 18.75892
!              t(4)= 0.0731370 ; q(4)= 266.3842
!              t(5)= 0.0857710 ; q(5)= 346.3755
!              t(6)= 0.1029220 ; q(6)= 413.8419
!              t(7)= 0.1154270 ; q(7)= 424.2680
!              t(8)= 0.1483530 ; q(8)= 429.1147
!              t(9)= 0.1698860 ; q(9)= 411.0127
!              t(10)= 0.220794 ; q(10)= 319.151
!              t(11)= 0.264856 ; q(11)= 207.816
!              t(12)= 0.295415 ; q(12)= 160.490
!              t(13)= 0.325895 ; q(13)= 70.0342
!              t(14)= 0.346215 ; q(14)= 10.1939
!              t(15)= 0.363213 ; q(15)= -5.1222
!              t(16)= 0.383666 ; q(16)= 6.68963
!              t(17)= 0.405265 ; q(17)= 24.0659
!              t(18)= 0.427988 ; q(18)= 35.8762
!              t(19)= 0.455272 ; q(19)= 58.8137
!              t(20)= 0.477990 ; q(20)= 67.8414
!              t(21)= 0.502943 ; q(21)= 57.3893
!              t(22)= 0.535816 ; q(22)= 33.7142
!              t(23)= 0.577789 ; q(23)= 20.4676
!              t(24)= 0.602753 ; q(24)= 16.2763
!              t(25)= 0.639087 ; q(25)= 22.5119
!              t(26)= 0.727616 ; q(26)= 18.9721
!              t(27)= 0.783235 ; q(27)= 18.9334
!              t(28)= 0.800000 ; q(28)= 16.1121
              
              t(1)=0.0026864; q(1)=1.7576
              t(2)=0.018499 ; q(2)=1.73771
              t(3)=0.0220232; q(3)=7.77101
              t(4)=0.0246908; q(4)=26.6356
              t(5)=0.0308954; q(5)=59.0806
              t(6)=0.0371116; q(6)=98.3181
              t(7)=0.0433458; q(7)=148.122
              t(8)=0.0495826; q(8)=199.434
              t(9)=0.0549268; q(9)=242.447
              t(10)=0.059411; q(10)=296.781
              t(11)=0.069157; q(11)=345.07
              t(12)=0.074504; q(12)=389.592
              t(13)=0.081611; q(13)=436.375
              t(14)=0.086072; q(14)=476.37
              t(15)=0.093169; q(15)=517.115
              t(16)=0.099381; q(16)=554.089
              t(17)=0.112609; q(17)=583.506
              t(18)=0.128445; q(18)=597.071
              t(19)=0.140728; q(19)=587.999
              t(20)=0.151236; q(20)=568.363
              t(21)=0.167866; q(21)=532.116
              t(22)=0.186242; q(22)=489.829
              t(23)=0.20112 ; q(23)=457.357
              t(24)=0.214231; q(24)=418.095
              t(25)=0.228209; q(25)=372.795
              t(26)=0.240446; q(26)=336.553
              t(27)=0.260584; q(27)=297.282
              t(28)=0.279854; q(28)=264.05
              t(29)=0.303502; q(29)=222.511
              t(30)=0.319267; q(30)=194.567
              t(31)=0.336785; q(31)=164.356
              t(32)=0.350805; q(32)=143.206
              t(33)=0.358661; q(33)=113.763
              t(34)=0.369132; q(34)=72.24
              t(35)=0.378743; q(35)=42.0393
              t(36)=0.386592; q(36)=8.06715
              t(37)=0.390934; q(37)=-21.3723
              t(38)=0.394398; q(38)=-50.8106
              t(39)=0.400493; q(39)=-82.5164
              t(40)=0.404843; q(40)=-107.427
              t(41)=0.409199; q(41)=-128.565
              t(42)=0.411801; q(42)=-148.191
              t(43)=0.414423; q(43)=-156.496
              t(44)=0.421485; q(44)=-136.128
              t(45)=0.425057; q(45)=-102.17
              t(46)=0.428627; q(46)=-68.9669
              t(47)=0.433945; q(47)=-41.8037
              t(48)=0.436613; q(48)=-22.1844
              t(49)=0.452445; q(49)=-10.8836
              t(50)=0.470906; q(50)=-3.35966
              t(51)=0.488486; q(51)=2.65596
              t(52)=0.510453; q(52)=5.6472
              t(53)=0.530659; q(53)=6.3765
              t(54)=0.557884; q(54)=1.81395
              t(55)=0.601808; q(55)=1.7587
              t(56)=0.668572; q(56)=1.67472
              t(57)=0.725673; q(57)=1.6029
              t(58)=0.768718; q(58)=1.54875
              t(59)=0.828453; q(59)=0.718895
              t(60)=0.884676; q(60)=0.648175
              t(61)=0.937384; q(61)=0.581875
              t(62)=0.994487; q(62)=1.26477

              t(1)=0.000  ;q(1)=0.0
              t(2)=0.001  ;q(2)=0.0
              t(3)=0.002  ;q(3)=0.0
              t(4)=0.003  ;q(4)=0.0
              t(5)=0.004  ;q(5)=0.0
              t(6)=0.005  ;q(6)=0.0
              t(7)=0.006  ;q(7)=0.0
              t(8)=0.007  ;q(8)=0.0
              t(9)=0.008  ;q(9)=0.0
              t(10)=0.009 ;q(10)=0.0
              t(11)=0.010 ;q(11)=0.0
              t(12)=0.011 ;q(12)=0.0
              t(13)=0.012 ;q(13)=0.0
              t(14)=0.013 ;q(14)=0.0
              t(15)=0.014 ;q(15)=0.0
              t(16)=0.015 ;q(16)=0.0
              t(17)=0.016 ;q(17)=0.0
              t(18)=0.017 ;q(18)=0.0
              t(19)=0.018 ;q(19)=0.0
              t(20)=0.019 ;q(20)=0.0
              t(21)=0.020 ;q(21)=0.0
              t(22)=0.021 ;q(22)=0.0
              t(23)=0.022 ;q(23)=0.0
              t(24)=0.023 ;q(24)=0.0
              t(25)=0.024 ;q(25)=0.0
              t(26)=0.025 ;q(26)=0.0
              t(27)=0.026 ;q(27)=0.0
              t(28)=0.027 ;q(28)=0.0
              t(29)=0.028 ;q(29)=0.0
              t(30)=0.029 ;q(30)=0.0
              t(31)=0.030 ;q(31)=0.7248229683
              t(32)=0.031 ;q(32)=6.6225270489
              t(33)=0.032 ;q(33)=15.025980663
              t(34)=0.033 ;q(34)=27.7409438793
              t(35)=0.034 ;q(35)=43.3417402917
              t(36)=0.035 ;q(36)=60.1022796466
              t(37)=0.036 ;q(37)=76.2964716907
              t(38)=0.037 ;q(38)=90.3303569383
              t(39)=0.038 ;q(39)=103.4647748655
              t(40)=0.039 ;q(40)=116.4847535959
              t(41)=0.040 ;q(41)=129.3569670119
              t(42)=0.041 ;q(42)=142.0480889957
              t(43)=0.042 ;q(43)=154.5247934298
              t(44)=0.043 ;q(44)=166.7537541964
              t(45)=0.044 ;q(45)=178.7016451779
              t(46)=0.045 ;q(46)=190.3351402567
              t(47)=0.046 ;q(47)=201.620913315
              t(48)=0.047 ;q(48)=212.5256382351
              t(49)=0.048 ;q(49)=223.0159888995
              t(50)=0.049 ;q(50)=233.0672155411
              t(51)=0.050 ;q(51)=242.8570583364
              t(52)=0.051 ;q(52)=252.463212555
              t(53)=0.052 ;q(53)=261.8794127274
              t(54)=0.053 ;q(54)=271.0993933839
              t(55)=0.054 ;q(55)=280.1168890551
              t(56)=0.055 ;q(56)=288.9256342713
              t(57)=0.056 ;q(57)=297.5193635629
              t(58)=0.057 ;q(58)=305.8918114603
              t(59)=0.058 ;q(59)=314.0367124939
              t(60)=0.059 ;q(60)=321.9478011942
              t(61)=0.060 ;q(61)=329.6188120916
              t(62)=0.061 ;q(62)=337.0434797164
              t(63)=0.062 ;q(63)=344.2155385991
              t(64)=0.063 ;q(64)=351.1287232701
              t(65)=0.064 ;q(65)=357.7817524418
              t(66)=0.065 ;q(66)=364.2182991294
              t(67)=0.066 ;q(67)=370.4576786525
              t(68)=0.067 ;q(68)=376.5095235019
              t(69)=0.068 ;q(69)=382.3834661685
              t(70)=0.069 ;q(70)=388.0891391431
              t(71)=0.070 ;q(71)=393.6361749164
              t(72)=0.071 ;q(72)=399.0342059792
              t(73)=0.072 ;q(73)=404.2928648224
              t(74)=0.073 ;q(74)=409.4217839368
              t(75)=0.074 ;q(75)=414.430595813
              t(76)=0.075 ;q(76)=419.328932942
              t(77)=0.076 ;q(77)=424.1264278145
              t(78)=0.077 ;q(78)=428.8332084822
              t(79)=0.078 ;q(79)=433.4971521949
              t(80)=0.079 ;q(80)=438.1326539139
              t(81)=0.080 ;q(81)=442.7221864366
              t(82)=0.081 ;q(82)=447.2482225606
              t(83)=0.082 ;q(83)=451.693235083
              t(84)=0.083 ;q(84)=456.0396968015
              t(85)=0.084 ;q(85)=460.2700805133
              t(86)=0.085 ;q(86)=464.3668590158
              t(87)=0.086 ;q(87)=468.3125051065
              t(88)=0.087 ;q(88)=472.0894915827
              t(89)=0.088 ;q(89)=475.6802912419
              t(90)=0.089 ;q(90)=479.0673768814
              t(91)=0.090 ;q(91)=482.2332212986
              t(92)=0.091 ;q(92)=485.1692108415
              t(93)=0.092 ;q(93)=488.0305733133
              t(94)=0.093 ;q(94)=490.8733292206
              t(95)=0.094 ;q(95)=493.6893030979
              t(96)=0.095 ;q(96)=496.4703194794
              t(97)=0.096 ;q(97)=499.2082028995
              t(98)=0.097 ;q(98)=501.8947778924
              t(99)=0.098 ;q(99)=504.5218689926
              t(100)=0.099;q(100)=507.0813007343
              t(101)=0.100;q(101)=509.5648976519
              t(102)=0.101;q(102)=511.9644842797
              t(103)=0.102;q(103)=514.2718851519
              t(104)=0.103;q(104)=516.4789248031
              t(105)=0.104;q(105)=518.5774277674
              t(106)=0.105;q(106)=520.5592185792
              t(107)=0.106;q(107)=522.4161217728
              t(108)=0.107;q(108)=524.1399618826
              t(109)=0.108;q(109)=525.7225634428
              t(110)=0.109;q(110)=527.1557509879
              t(111)=0.110;q(111)=528.4313490521
              t(112)=0.111;q(112)=529.5411821697
              t(113)=0.112;q(113)=530.4770748751
              t(114)=0.113;q(114)=531.2706129745
              t(115)=0.114;q(115)=532.0417084103
              t(116)=0.115;q(116)=532.7967212634
              t(117)=0.116;q(117)=533.5328892761
              t(118)=0.117;q(118)=534.2474501905
              t(119)=0.118;q(119)=534.9376417488
              t(120)=0.119;q(120)=535.6007016933
              t(121)=0.120;q(121)=536.2338677661
              t(122)=0.121;q(122)=536.8343777095
              t(123)=0.122;q(123)=537.3994692657
              t(124)=0.123;q(124)=537.9263801768
              t(125)=0.124;q(125)=538.4123481851
              t(126)=0.125;q(126)=538.8546110327
              t(127)=0.126;q(127)=539.250406462
              t(128)=0.127;q(128)=539.596972215
              t(129)=0.128;q(129)=539.891546034
              t(130)=0.129;q(130)=540.1313656612
              t(131)=0.130;q(131)=540.3136688388
              t(132)=0.131;q(132)=540.435693309
              t(133)=0.132;q(133)=540.494676814
              t(134)=0.133;q(134)=540.4909793267
              t(135)=0.134;q(135)=540.4366607062
              t(136)=0.135;q(136)=540.3348775531
              t(137)=0.136;q(137)=540.1875229525
              t(138)=0.137;q(138)=539.9964899894
              t(139)=0.138;q(139)=539.7636717489
              t(140)=0.139;q(140)=539.490961316
              t(141)=0.140;q(141)=539.1802517756
              t(142)=0.141;q(142)=538.833436213
              t(143)=0.142;q(143)=538.4524077131
              t(144)=0.143;q(144)=538.0390593609
              t(145)=0.144;q(145)=537.5952842415
              t(146)=0.145;q(146)=537.1229754399
              t(147)=0.146;q(147)=536.6240260412
              t(148)=0.147;q(148)=536.1003291305
              t(149)=0.148;q(149)=535.5537777926
              t(150)=0.149;q(150)=534.9862651128
              t(151)=0.150;q(151)=534.399684176
              t(152)=0.151;q(152)=533.7959280673
              t(153)=0.152;q(153)=533.1768898717
              t(154)=0.153;q(154)=532.5444626742
              t(155)=0.154;q(155)=531.9005395599
              t(156)=0.155;q(156)=531.2470136139
              t(157)=0.156;q(157)=530.5823272151
              t(158)=0.157;q(158)=529.825055976
              t(159)=0.158;q(159)=528.9444056288
              t(160)=0.159;q(160)=527.9461200374
              t(161)=0.160;q(161)=526.835943066
              t(162)=0.161;q(162)=525.6196185787
              t(163)=0.162;q(163)=524.3028904396
              t(164)=0.163;q(164)=522.8915025127
              t(165)=0.164;q(165)=521.3911986622
              t(166)=0.165;q(166)=519.8077227522
              t(167)=0.166;q(167)=518.1468186468
              t(168)=0.167;q(168)=516.41423021
              t(169)=0.168;q(169)=514.615701306
              t(170)=0.169;q(170)=512.7569757989
              t(171)=0.170;q(171)=510.8437975527
              t(172)=0.171;q(172)=508.8819104315
              t(173)=0.172;q(173)=506.8770582995
              t(174)=0.173;q(174)=504.8349850208
              t(175)=0.174;q(175)=502.7614344594
              t(176)=0.175;q(176)=500.6621504794
              t(177)=0.176;q(177)=498.542876945
              t(178)=0.177;q(178)=496.4093577202
              t(179)=0.178;q(179)=494.2673366691
              t(180)=0.179;q(180)=492.1225576558
              t(181)=0.180;q(181)=489.9807645445
              t(182)=0.181;q(182)=487.8477011992
              t(183)=0.182;q(183)=485.729111484
              t(184)=0.183;q(184)=483.6245705444
              t(185)=0.184;q(185)=481.4709026788
              t(186)=0.185;q(186)=479.2517184466
              t(187)=0.186;q(187)=476.96904946
              t(188)=0.187;q(188)=474.6249273313
              t(189)=0.188;q(189)=472.2213836728
              t(190)=0.189;q(190)=469.7604500965
              t(191)=0.190;q(191)=467.2441582149
              t(192)=0.191;q(192)=464.6745396402
              t(193)=0.192;q(193)=462.0536259846
              t(194)=0.193;q(194)=459.3834488603
              t(195)=0.194;q(195)=456.6660398796
              t(196)=0.195;q(196)=453.9034306548
              t(197)=0.196;q(197)=451.097652798
              t(198)=0.197;q(198)=448.2507379216
              t(199)=0.198;q(199)=445.3647176378
              t(200)=0.199;q(200)=442.4416235588
              t(201)=0.200;q(201)=439.4834872969
              t(202)=0.201;q(202)=436.4923404643
              t(203)=0.202;q(203)=433.4702146733
              t(204)=0.203;q(204)=430.419141536
              t(205)=0.204;q(205)=427.3411526649
              t(206)=0.205;q(206)=424.238279672
              t(207)=0.206;q(207)=421.1125541697
              t(208)=0.207;q(208)=417.9660077702
              t(209)=0.208;q(209)=414.8006720858
              t(210)=0.209;q(210)=411.6185787286
              t(211)=0.210;q(211)=408.421759311
              t(212)=0.211;q(212)=405.2122454452
              t(213)=0.212;q(213)=401.9920687434
              t(214)=0.213;q(214)=398.7632608179
              t(215)=0.214;q(215)=395.5278532808
              t(216)=0.215;q(216)=392.2878777446
              t(217)=0.216;q(217)=389.0453658214
              t(218)=0.217;q(218)=385.8023491234
              t(219)=0.218;q(219)=382.560859263
              t(220)=0.219;q(220)=379.3229278523
              t(221)=0.220;q(221)=376.0905865036
              t(222)=0.221;q(222)=372.8658668292
              t(223)=0.222;q(223)=369.6508004412
              t(224)=0.223;q(224)=366.447418952
              t(225)=0.224;q(225)=363.2577539738
              t(226)=0.225;q(226)=360.0838371188
              t(227)=0.226;q(227)=356.9276999992
              t(228)=0.227;q(228)=353.7913742275
              t(229)=0.228;q(229)=350.6768914156
              t(230)=0.229;q(230)=347.586283176
              t(231)=0.230;q(231)=344.5215811209
              t(232)=0.231;q(232)=341.4848168624
              t(233)=0.232;q(233)=338.478022013
              t(234)=0.233;q(234)=335.5032281847
              t(235)=0.234;q(235)=332.5624669899
              t(236)=0.235;q(236)=329.6577700407
              t(237)=0.236;q(237)=326.7911689495
              t(238)=0.237;q(238)=323.9646953285
              t(239)=0.238;q(239)=321.18038079
              t(240)=0.239;q(240)=318.4402569461
              t(241)=0.240;q(241)=315.7463554091
              t(242)=0.241;q(242)=313.1007077913
              t(243)=0.242;q(243)=310.5053457049
              t(244)=0.243;q(244)=307.9623007622
              t(245)=0.244;q(245)=305.4736045754
              t(246)=0.245;q(246)=303.0389089214
              t(247)=0.246;q(247)=300.6308183643
              t(248)=0.247;q(248)=298.2407695524
              t(249)=0.248;q(249)=295.8684611365
              t(250)=0.249;q(250)=293.5135917679
              t(251)=0.250;q(251)=291.1758600976
              t(252)=0.251;q(252)=288.8549647765
              t(253)=0.252;q(253)=286.5506044558
              t(254)=0.253;q(254)=284.2624777864
              t(255)=0.254;q(255)=281.9902834195
              t(256)=0.255;q(256)=279.733720006
              t(257)=0.256;q(257)=277.4924861971
              t(258)=0.257;q(258)=275.2662806436
              t(259)=0.258;q(259)=273.0548019968
              t(260)=0.259;q(260)=270.8577489077
              t(261)=0.260;q(261)=268.6748200272
              t(262)=0.261;q(262)=266.5057140064
              t(263)=0.262;q(263)=264.3501294965
              t(264)=0.263;q(264)=262.2077651483
              t(265)=0.264;q(265)=260.078319613
              t(266)=0.265;q(266)=257.9614915416
              t(267)=0.266;q(267)=255.8569795852
              t(268)=0.267;q(268)=253.7644823947
              t(269)=0.268;q(269)=251.6836986213
              t(270)=0.269;q(270)=249.614326916
              t(271)=0.270;q(271)=247.5560659298
              t(272)=0.271;q(272)=245.5086143138
              t(273)=0.272;q(273)=243.471670719
              t(274)=0.273;q(274)=241.4449337965
              t(275)=0.274;q(275)=239.4281021972
              t(276)=0.275;q(276)=237.4208745723
              t(277)=0.276;q(277)=235.4229495728
              t(278)=0.277;q(278)=233.4340258498
              t(279)=0.278;q(279)=231.4538020542
              t(280)=0.279;q(280)=229.4819768372
              t(281)=0.280;q(281)=227.5182488497
              t(282)=0.281;q(282)=225.5623167428
              t(283)=0.282;q(283)=223.6138791676
              t(284)=0.283;q(284)=221.6726347751
              t(285)=0.284;q(285)=219.7382822163
              t(286)=0.285;q(286)=217.8105201423
              t(287)=0.286;q(287)=215.8890472042
              t(288)=0.287;q(288)=213.9735620529
              t(289)=0.288;q(289)=212.0637633396
              t(290)=0.289;q(290)=210.1593497152
              t(291)=0.290;q(291)=208.2600198309
              t(292)=0.291;q(292)=206.3654723376
              t(293)=0.292;q(293)=204.4822126615
              t(294)=0.293;q(294)=202.6173895501
              t(295)=0.294;q(295)=200.7702541112
              t(296)=0.295;q(296)=198.9400494428
              t(297)=0.296;q(297)=197.126018643
              t(298)=0.297;q(298)=195.3274048098
              t(299)=0.298;q(299)=193.5434510413
              t(300)=0.299;q(300)=191.7734004356
              t(301)=0.300;q(301)=190.0164960908
              t(302)=0.301;q(302)=188.2719811047
              t(303)=0.302;q(303)=186.5390985757
              t(304)=0.303;q(304)=184.8170916016
              t(305)=0.304;q(305)=183.1052032805
              t(306)=0.305;q(306)=181.4026767106
              t(307)=0.306;q(307)=179.7087549898
              t(308)=0.307;q(308)=178.0226812162
              t(309)=0.308;q(309)=176.3436984879
              t(310)=0.309;q(310)=174.6710499029
              t(311)=0.310;q(311)=173.0041586407
              t(312)=0.311;q(312)=171.3431337921
              t(313)=0.312;q(313)=169.6882406982
              t(314)=0.313;q(314)=168.0397446348
              t(315)=0.314;q(315)=166.3979108778
              t(316)=0.315;q(316)=164.7630047029
              t(317)=0.316;q(317)=163.1352913859
              t(318)=0.317;q(318)=161.5150362028
              t(319)=0.318;q(319)=159.9025044293
              t(320)=0.319;q(320)=158.2979613412
              t(321)=0.320;q(321)=156.7016722143
              t(322)=0.321;q(322)=155.1139023245
              t(323)=0.322;q(323)=153.5349169476
              t(324)=0.323;q(324)=151.9649813593
              t(325)=0.324;q(325)=150.4043608356
              t(326)=0.325;q(326)=148.8533206522
              t(327)=0.326;q(327)=147.3121260849
              t(328)=0.327;q(328)=145.7810424096
              t(329)=0.328;q(329)=144.260334902
              t(330)=0.329;q(330)=142.7502688381
              t(331)=0.330;q(331)=141.2511094935
              t(332)=0.331;q(332)=139.7631221442
              t(333)=0.332;q(333)=138.2865720659
              t(334)=0.333;q(334)=136.8280814549
              t(335)=0.334;q(335)=135.3989811512
              t(336)=0.335;q(336)=133.9965718797
              t(337)=0.336;q(337)=132.6177600653
              t(338)=0.337;q(338)=131.259452133
              t(339)=0.338;q(339)=129.9185545076
              t(340)=0.339;q(340)=128.5919736141
              t(341)=0.340;q(341)=127.2766158775
              t(342)=0.341;q(342)=125.9693877226
              t(343)=0.342;q(343)=124.6671955744
              t(344)=0.343;q(344)=123.3669458579
              t(345)=0.344;q(345)=122.0655449979
              t(346)=0.345;q(346)=120.7598994194
              t(347)=0.346;q(347)=119.4469155473
              t(348)=0.347;q(348)=118.1234998066
              t(349)=0.348;q(349)=116.7865586222
              t(350)=0.349;q(350)=115.4329984189
              t(351)=0.350;q(351)=114.0597256219
              t(352)=0.351;q(352)=112.6636466559
              t(353)=0.352;q(353)=111.2416679459
              t(354)=0.353;q(354)=109.7906959169
              t(355)=0.354;q(355)=108.3076369937
              t(356)=0.355;q(356)=106.7893976013
              t(357)=0.356;q(357)=105.2515034784
              t(358)=0.357;q(358)=103.7458309618
              t(359)=0.358;q(359)=102.2655664064
              t(360)=0.359;q(360)=100.799435745
              t(361)=0.360;q(361)=99.3361649102
              t(362)=0.361;q(362)=97.8644798347
              t(363)=0.362;q(363)=96.3731064511
              t(364)=0.363;q(364)=94.8507706922
              t(365)=0.364;q(365)=93.2861984906
              t(366)=0.365;q(366)=91.6681157789
              t(367)=0.366;q(367)=89.9852484899
              t(368)=0.367;q(368)=88.2263225562
              t(369)=0.368;q(369)=86.3800639104
              t(370)=0.369;q(370)=84.4351984853
              t(371)=0.370;q(371)=82.3804522136
              t(372)=0.371;q(372)=80.2045510278
              t(373)=0.372;q(373)=77.8962208607
              t(374)=0.373;q(374)=75.444187645
              t(375)=0.374;q(375)=72.7688670192
              t(376)=0.375;q(376)=69.6593077452
              t(377)=0.376;q(377)=66.144878867
              t(378)=0.377;q(378)=62.2742322619
              t(379)=0.378;q(379)=58.0960198071
              t(380)=0.379;q(380)=53.6588933796
              t(381)=0.380;q(381)=49.0115048567
              t(382)=0.381;q(382)=44.2025061155
              t(383)=0.382;q(383)=39.2805490333
              t(384)=0.383;q(384)=34.2942854872
              t(385)=0.384;q(385)=29.2631052947
              t(386)=0.385;q(386)=23.7010207386
              t(387)=0.386;q(387)=17.4754843361
              t(388)=0.387;q(388)=10.6594620234
              t(389)=0.388;q(389)=3.3259197367
              t(390)=0.389;q(390)=-4.4521765878
              t(391)=0.390;q(391)=-12.6018610139
              t(392)=0.391;q(392)=-21.0501676054
              t(393)=0.392;q(393)=-29.7241304262
              t(394)=0.393;q(394)=-38.5507835401
              t(395)=0.394;q(395)=-47.4571610108
              t(396)=0.395;q(396)=-56.3702969022
              t(397)=0.396;q(397)=-65.2172252782
              t(398)=0.397;q(398)=-73.9249802024
              t(399)=0.398;q(399)=-82.4205957388
              t(400)=0.399;q(400)=-90.6311059512
              t(401)=0.400;q(401)=-98.4835449033
              t(402)=0.401;q(402)=-105.904946659
              t(403)=0.402;q(403)=-112.8223452821
              t(404)=0.403;q(404)=-119.1627748364
              t(405)=0.404;q(405)=-124.8532693858
              t(406)=0.405;q(406)=-129.820862994
              t(407)=0.406;q(407)=-133.9925897248
              t(408)=0.407;q(408)=-137.2954836422
              t(409)=0.408;q(409)=-139.6565788098
              t(410)=0.409;q(410)=-141.0029092915
              t(411)=0.410;q(411)=-141.2128903389
              t(412)=0.411;q(412)=-139.3941327416
              t(413)=0.412;q(413)=-135.5107556468
              t(414)=0.413;q(414)=-129.9246387691
              t(415)=0.414;q(415)=-122.9976618231
              t(416)=0.415;q(416)=-115.0917045235
              t(417)=0.416;q(417)=-106.5686465848
              t(418)=0.417;q(418)=-97.7903677216
              t(419)=0.418;q(419)=-89.1187476486
              t(420)=0.419;q(420)=-80.9156660803
              t(421)=0.420;q(421)=-73.5178158057
              t(422)=0.421;q(422)=-65.5825856848
              t(423)=0.422;q(423)=-56.5169011857
              t(424)=0.423;q(424)=-46.7935106867
              t(425)=0.424;q(425)=-36.8851625662
              t(426)=0.425;q(426)=-27.2646052024
              t(427)=0.426;q(427)=-18.4045869736
              t(428)=0.427;q(428)=-10.7778562581
              t(429)=0.428;q(429)=-4.8571614342
              t(430)=0.429;q(430)=-1.1152508802
              t(431)=0.430;q(431)=0.0
              t(432)=1.000;q(432)=0.0

t(433)=1.001;q(433)=0
t(434)=1.002;q(434)=0
t(435)=1.003;q(435)=0
t(436)=1.004;q(436)=0
t(437)=1.005;q(437)=0
t(438)=1.006;q(438)=0
t(439)=1.007;q(439)=0
t(440)=1.008;q(440)=0
t(441)=1.009;q(441)=0
t(442)=1.01;q(442)=0
t(443)=1.011;q(443)=0
t(444)=1.012;q(444)=0
t(445)=1.013;q(445)=0
t(446)=1.014;q(446)=0
t(447)=1.015;q(447)=0
t(448)=1.016;q(448)=0
t(449)=1.017;q(449)=0
t(450)=1.018;q(450)=0
t(451)=1.019;q(451)=0
t(452)=1.02;q(452)=0
t(453)=1.021;q(453)=0
t(454)=1.022;q(454)=0
t(455)=1.023;q(455)=0
t(456)=1.024;q(456)=0
t(457)=1.025;q(457)=0
t(458)=1.026;q(458)=0
t(459)=1.027;q(459)=0
t(460)=1.028;q(460)=0
t(461)=1.029;q(461)=0
t(462)=1.03;q(462)=0.7248229683
t(463)=1.031;q(463)=6.6225270489
t(464)=1.032;q(464)=15.025980663
t(465)=1.033;q(465)=27.7409438793
t(466)=1.034;q(466)=43.3417402917
t(467)=1.035;q(467)=60.1022796466
t(468)=1.036;q(468)=76.2964716907
t(469)=1.037;q(469)=90.3303569383
t(470)=1.038;q(470)=103.4647748655
t(471)=1.039;q(471)=116.4847535959
t(472)=1.04;q(472)=129.3569670119
t(473)=1.041;q(473)=142.0480889957
t(474)=1.042;q(474)=154.5247934298
t(475)=1.043;q(475)=166.7537541964
t(476)=1.044;q(476)=178.7016451779
t(477)=1.045;q(477)=190.3351402567
t(478)=1.046;q(478)=201.620913315
t(479)=1.047;q(479)=212.5256382351
t(480)=1.048;q(480)=223.0159888995
t(481)=1.049;q(481)=233.0672155411
t(482)=1.05;q(482)=242.8570583364
t(483)=1.051;q(483)=252.463212555
t(484)=1.052;q(484)=261.8794127274
t(485)=1.053;q(485)=271.0993933839
t(486)=1.054;q(486)=280.1168890551
t(487)=1.055;q(487)=288.9256342713
t(488)=1.056;q(488)=297.5193635629
t(489)=1.057;q(489)=305.8918114603
t(490)=1.058;q(490)=314.0367124939
t(491)=1.059;q(491)=321.9478011942
t(492)=1.06;q(492)=329.6188120916
t(493)=1.061;q(493)=337.0434797164
t(494)=1.062;q(494)=344.2155385991
t(495)=1.063;q(495)=351.1287232701
t(496)=1.064;q(496)=357.7817524418
t(497)=1.065;q(497)=364.2182991294
t(498)=1.066;q(498)=370.4576786525
t(499)=1.067;q(499)=376.5095235019
t(500)=1.068;q(500)=382.3834661685
t(501)=1.069;q(501)=388.0891391431
t(502)=1.07;q(502)=393.6361749164
t(503)=1.071;q(503)=399.0342059792
t(504)=1.072;q(504)=404.2928648224
t(505)=1.073;q(505)=409.4217839368
t(506)=1.074;q(506)=414.430595813
t(507)=1.075;q(507)=419.328932942
t(508)=1.076;q(508)=424.1264278145
t(509)=1.077;q(509)=428.8332084822
t(510)=1.078;q(510)=433.4971521949
t(511)=1.079;q(511)=438.1326539139
t(512)=1.08;q(512)=442.7221864366
t(513)=1.081;q(513)=447.2482225606
t(514)=1.082;q(514)=451.693235083
t(515)=1.083;q(515)=456.0396968015
t(516)=1.084;q(516)=460.2700805133
t(517)=1.085;q(517)=464.3668590158
t(518)=1.086;q(518)=468.3125051065
t(519)=1.087;q(519)=472.0894915827
t(520)=1.088;q(520)=475.6802912419
t(521)=1.089;q(521)=479.0673768814
t(522)=1.09;q(522)=482.2332212986
t(523)=1.091;q(523)=485.1692108415
t(524)=1.092;q(524)=488.0305733133
t(525)=1.093;q(525)=490.8733292206
t(526)=1.094;q(526)=493.6893030979
t(527)=1.095;q(527)=496.4703194794
t(528)=1.096;q(528)=499.2082028995
t(529)=1.097;q(529)=501.8947778924
t(530)=1.098;q(530)=504.5218689926
t(531)=1.099;q(531)=507.0813007343
t(532)=1.1;q(532)=509.5648976519
t(533)=1.101;q(533)=511.9644842797
t(534)=1.102;q(534)=514.2718851519
t(535)=1.103;q(535)=516.4789248031
t(536)=1.104;q(536)=518.5774277674
t(537)=1.105;q(537)=520.5592185792
t(538)=1.106;q(538)=522.4161217728
t(539)=1.107;q(539)=524.1399618826
t(540)=1.108;q(540)=525.7225634428
t(541)=1.109;q(541)=527.1557509879
t(542)=1.11;q(542)=528.4313490521
t(543)=1.111;q(543)=529.5411821697
t(544)=1.112;q(544)=530.4770748751
t(545)=1.113;q(545)=531.2706129745
t(546)=1.114;q(546)=532.0417084103
t(547)=1.115;q(547)=532.7967212634
t(548)=1.116;q(548)=533.5328892761
t(549)=1.117;q(549)=534.2474501905
t(550)=1.118;q(550)=534.9376417488
t(551)=1.119;q(551)=535.6007016933
t(552)=1.12;q(552)=536.2338677661
t(553)=1.121;q(553)=536.8343777095
t(554)=1.122;q(554)=537.3994692657
t(555)=1.123;q(555)=537.9263801768
t(556)=1.124;q(556)=538.4123481851
t(557)=1.125;q(557)=538.8546110327
t(558)=1.126;q(558)=539.250406462
t(559)=1.127;q(559)=539.596972215
t(560)=1.128;q(560)=539.891546034
t(561)=1.129;q(561)=540.1313656612
t(562)=1.13;q(562)=540.3136688388
t(563)=1.131;q(563)=540.435693309
t(564)=1.132;q(564)=540.494676814
t(565)=1.133;q(565)=540.4909793267
t(566)=1.134;q(566)=540.4366607062
t(567)=1.135;q(567)=540.3348775531
t(568)=1.136;q(568)=540.1875229525
t(569)=1.137;q(569)=539.9964899894
t(570)=1.138;q(570)=539.7636717489
t(571)=1.139;q(571)=539.490961316
t(572)=1.14;q(572)=539.1802517756
t(573)=1.141;q(573)=538.833436213
t(574)=1.142;q(574)=538.4524077131
t(575)=1.143;q(575)=538.0390593609
t(576)=1.144;q(576)=537.5952842415
t(577)=1.145;q(577)=537.1229754399
t(578)=1.146;q(578)=536.6240260412
t(579)=1.147;q(579)=536.1003291305
t(580)=1.148;q(580)=535.5537777926
t(581)=1.149;q(581)=534.9862651128
t(582)=1.15;q(582)=534.399684176
t(583)=1.151;q(583)=533.7959280673
t(584)=1.152;q(584)=533.1768898717
t(585)=1.153;q(585)=532.5444626742
t(586)=1.154;q(586)=531.9005395599
t(587)=1.155;q(587)=531.2470136139
t(588)=1.156;q(588)=530.5823272151
t(589)=1.157;q(589)=529.825055976
t(590)=1.158;q(590)=528.9444056288
t(591)=1.159;q(591)=527.9461200374
t(592)=1.16;q(592)=526.835943066
t(593)=1.161;q(593)=525.6196185787
t(594)=1.162;q(594)=524.3028904396
t(595)=1.163;q(595)=522.8915025127
t(596)=1.164;q(596)=521.3911986622
t(597)=1.165;q(597)=519.8077227522
t(598)=1.166;q(598)=518.1468186468
t(599)=1.167;q(599)=516.41423021
t(600)=1.168;q(600)=514.615701306
t(601)=1.169;q(601)=512.7569757989
t(602)=1.17;q(602)=510.8437975527
t(603)=1.171;q(603)=508.8819104315
t(604)=1.172;q(604)=506.8770582995
t(605)=1.173;q(605)=504.8349850208
t(606)=1.174;q(606)=502.7614344594
t(607)=1.175;q(607)=500.6621504794
t(608)=1.176;q(608)=498.542876945
t(609)=1.177;q(609)=496.4093577202
t(610)=1.178;q(610)=494.2673366691
t(611)=1.179;q(611)=492.1225576558
t(612)=1.18;q(612)=489.9807645445
t(613)=1.181;q(613)=487.8477011992
t(614)=1.182;q(614)=485.729111484
t(615)=1.183;q(615)=483.6245705444
t(616)=1.184;q(616)=481.4709026788
t(617)=1.185;q(617)=479.2517184466
t(618)=1.186;q(618)=476.96904946
t(619)=1.187;q(619)=474.6249273313
t(620)=1.188;q(620)=472.2213836728
t(621)=1.189;q(621)=469.7604500965
t(622)=1.19;q(622)=467.2441582149
t(623)=1.191;q(623)=464.6745396402
t(624)=1.192;q(624)=462.0536259846
t(625)=1.193;q(625)=459.3834488603
t(626)=1.194;q(626)=456.6660398796
t(627)=1.195;q(627)=453.9034306548
t(628)=1.196;q(628)=451.097652798
t(629)=1.197;q(629)=448.2507379216
t(630)=1.198;q(630)=445.3647176378
t(631)=1.199;q(631)=442.4416235588
t(632)=1.2;q(632)=439.4834872969
t(633)=1.201;q(633)=436.4923404643
t(634)=1.202;q(634)=433.4702146733
t(635)=1.203;q(635)=430.419141536
t(636)=1.204;q(636)=427.3411526649
t(637)=1.205;q(637)=424.238279672
t(638)=1.206;q(638)=421.1125541697
t(639)=1.207;q(639)=417.9660077702
t(640)=1.208;q(640)=414.8006720858
t(641)=1.209;q(641)=411.6185787286
t(642)=1.21;q(642)=408.421759311
t(643)=1.211;q(643)=405.2122454452
t(644)=1.212;q(644)=401.9920687434
t(645)=1.213;q(645)=398.7632608179
t(646)=1.214;q(646)=395.5278532808
t(647)=1.215;q(647)=392.2878777446
t(648)=1.216;q(648)=389.0453658214
t(649)=1.217;q(649)=385.8023491234
t(650)=1.218;q(650)=382.560859263
t(651)=1.219;q(651)=379.3229278523
t(652)=1.22;q(652)=376.0905865036
t(653)=1.221;q(653)=372.8658668292
t(654)=1.222;q(654)=369.6508004412
t(655)=1.223;q(655)=366.447418952
t(656)=1.224;q(656)=363.2577539738
t(657)=1.225;q(657)=360.0838371188
t(658)=1.226;q(658)=356.9276999992
t(659)=1.227;q(659)=353.7913742275
t(660)=1.228;q(660)=350.6768914156
t(661)=1.229;q(661)=347.586283176
t(662)=1.23;q(662)=344.5215811209
t(663)=1.231;q(663)=341.4848168624
t(664)=1.232;q(664)=338.478022013
t(665)=1.233;q(665)=335.5032281847
t(666)=1.234;q(666)=332.5624669899
t(667)=1.235;q(667)=329.6577700407
t(668)=1.236;q(668)=326.7911689495
t(669)=1.237;q(669)=323.9646953285
t(670)=1.238;q(670)=321.18038079
t(671)=1.239;q(671)=318.4402569461
t(672)=1.24;q(672)=315.7463554091
t(673)=1.241;q(673)=313.1007077913
t(674)=1.242;q(674)=310.5053457049
t(675)=1.243;q(675)=307.9623007622
t(676)=1.244;q(676)=305.4736045754
t(677)=1.245;q(677)=303.0389089214
t(678)=1.246;q(678)=300.6308183643
t(679)=1.247;q(679)=298.2407695524
t(680)=1.248;q(680)=295.8684611365
t(681)=1.249;q(681)=293.5135917679
t(682)=1.25;q(682)=291.1758600976
t(683)=1.251;q(683)=288.8549647765
t(684)=1.252;q(684)=286.5506044558
t(685)=1.253;q(685)=284.2624777864
t(686)=1.254;q(686)=281.9902834195
t(687)=1.255;q(687)=279.733720006
t(688)=1.256;q(688)=277.4924861971
t(689)=1.257;q(689)=275.2662806436
t(690)=1.258;q(690)=273.0548019968
t(691)=1.259;q(691)=270.8577489077
t(692)=1.26;q(692)=268.6748200272
t(693)=1.261;q(693)=266.5057140064
t(694)=1.262;q(694)=264.3501294965
t(695)=1.263;q(695)=262.2077651483
t(696)=1.264;q(696)=260.078319613
t(697)=1.265;q(697)=257.9614915416
t(698)=1.266;q(698)=255.8569795852
t(699)=1.267;q(699)=253.7644823947
t(700)=1.268;q(700)=251.6836986213
t(701)=1.269;q(701)=249.614326916
t(702)=1.27;q(702)=247.5560659298
t(703)=1.271;q(703)=245.5086143138
t(704)=1.272;q(704)=243.471670719
t(705)=1.273;q(705)=241.4449337965
t(706)=1.274;q(706)=239.4281021972
t(707)=1.275;q(707)=237.4208745723
t(708)=1.276;q(708)=235.4229495728
t(709)=1.277;q(709)=233.4340258498
t(710)=1.278;q(710)=231.4538020542
t(711)=1.279;q(711)=229.4819768372
t(712)=1.28;q(712)=227.5182488497
t(713)=1.281;q(713)=225.5623167428
t(714)=1.282;q(714)=223.6138791676
t(715)=1.283;q(715)=221.6726347751
t(716)=1.284;q(716)=219.7382822163
t(717)=1.285;q(717)=217.8105201423
t(718)=1.286;q(718)=215.8890472042
t(719)=1.287;q(719)=213.9735620529
t(720)=1.288;q(720)=212.0637633396
t(721)=1.289;q(721)=210.1593497152
t(722)=1.29;q(722)=208.2600198309
t(723)=1.291;q(723)=206.3654723376
t(724)=1.292;q(724)=204.4822126615
t(725)=1.293;q(725)=202.6173895501
t(726)=1.294;q(726)=200.7702541112
t(727)=1.295;q(727)=198.9400494428
t(728)=1.296;q(728)=197.126018643
t(729)=1.297;q(729)=195.3274048098
t(730)=1.298;q(730)=193.5434510413
t(731)=1.299;q(731)=191.7734004356
t(732)=1.3;q(732)=190.0164960908
t(733)=1.301;q(733)=188.2719811047
t(734)=1.302;q(734)=186.5390985757
t(735)=1.303;q(735)=184.8170916016
t(736)=1.304;q(736)=183.1052032805
t(737)=1.305;q(737)=181.4026767106
t(738)=1.306;q(738)=179.7087549898
t(739)=1.307;q(739)=178.0226812162
t(740)=1.308;q(740)=176.3436984879
t(741)=1.309;q(741)=174.6710499029
t(742)=1.31;q(742)=173.0041586407
t(743)=1.311;q(743)=171.3431337921
t(744)=1.312;q(744)=169.6882406982
t(745)=1.313;q(745)=168.0397446348
t(746)=1.314;q(746)=166.3979108778
t(747)=1.315;q(747)=164.7630047029
t(748)=1.316;q(748)=163.1352913859
t(749)=1.317;q(749)=161.5150362028
t(750)=1.318;q(750)=159.9025044293
t(751)=1.319;q(751)=158.2979613412
t(752)=1.32;q(752)=156.7016722143
t(753)=1.321;q(753)=155.1139023245
t(754)=1.322;q(754)=153.5349169476
t(755)=1.323;q(755)=151.9649813593
t(756)=1.324;q(756)=150.4043608356
t(757)=1.325;q(757)=148.8533206522
t(758)=1.326;q(758)=147.3121260849
t(759)=1.327;q(759)=145.7810424096
t(760)=1.328;q(760)=144.260334902
t(761)=1.329;q(761)=142.7502688381
t(762)=1.33;q(762)=141.2511094935
t(763)=1.331;q(763)=139.7631221442
t(764)=1.332;q(764)=138.2865720659
t(765)=1.333;q(765)=136.8280814549
t(766)=1.334;q(766)=135.3989811512
t(767)=1.335;q(767)=133.9965718797
t(768)=1.336;q(768)=132.6177600653
t(769)=1.337;q(769)=131.259452133
t(770)=1.338;q(770)=129.9185545076
t(771)=1.339;q(771)=128.5919736141
t(772)=1.34;q(772)=127.2766158775
t(773)=1.341;q(773)=125.9693877226
t(774)=1.342;q(774)=124.6671955744
t(775)=1.343;q(775)=123.3669458579
t(776)=1.344;q(776)=122.0655449979
t(777)=1.345;q(777)=120.7598994194
t(778)=1.346;q(778)=119.4469155473
t(779)=1.347;q(779)=118.1234998066
t(780)=1.348;q(780)=116.7865586222
t(781)=1.349;q(781)=115.4329984189
t(782)=1.35;q(782)=114.0597256219
t(783)=1.351;q(783)=112.6636466559
t(784)=1.352;q(784)=111.2416679459
t(785)=1.353;q(785)=109.7906959169
t(786)=1.354;q(786)=108.3076369937
t(787)=1.355;q(787)=106.7893976013
t(788)=1.356;q(788)=105.2515034784
t(789)=1.357;q(789)=103.7458309618
t(790)=1.358;q(790)=102.2655664064
t(791)=1.359;q(791)=100.799435745
t(792)=1.36;q(792)=99.3361649102
t(793)=1.361;q(793)=97.8644798347
t(794)=1.362;q(794)=96.3731064511
t(795)=1.363;q(795)=94.8507706922
t(796)=1.364;q(796)=93.2861984906
t(797)=1.365;q(797)=91.6681157789
t(798)=1.366;q(798)=89.9852484899
t(799)=1.367;q(799)=88.2263225562
t(800)=1.368;q(800)=86.3800639104
t(801)=1.369;q(801)=84.4351984853
t(802)=1.37;q(802)=82.3804522136
t(803)=1.371;q(803)=80.2045510278
t(804)=1.372;q(804)=77.8962208607
t(805)=1.373;q(805)=75.444187645
t(806)=1.374;q(806)=72.7688670192
t(807)=1.375;q(807)=69.6593077452
t(808)=1.376;q(808)=66.144878867
t(809)=1.377;q(809)=62.2742322619
t(810)=1.378;q(810)=58.0960198071
t(811)=1.379;q(811)=53.6588933796
t(812)=1.38;q(812)=49.0115048567
t(813)=1.381;q(813)=44.2025061155
t(814)=1.382;q(814)=39.2805490333
t(815)=1.383;q(815)=34.2942854872
t(816)=1.384;q(816)=29.2631052947
t(817)=1.385;q(817)=23.7010207386
t(818)=1.386;q(818)=17.4754843361
t(819)=1.387;q(819)=10.6594620234
t(820)=1.388;q(820)=3.3259197367
t(821)=1.389;q(821)=-4.4521765878
t(822)=1.39;q(822)=-12.6018610139
t(823)=1.391;q(823)=-21.0501676054
t(824)=1.392;q(824)=-29.7241304262
t(825)=1.393;q(825)=-38.5507835401
t(826)=1.394;q(826)=-47.4571610108
t(827)=1.395;q(827)=-56.3702969022
t(828)=1.396;q(828)=-65.2172252782
t(829)=1.397;q(829)=-73.9249802024
t(830)=1.398;q(830)=-82.4205957388
t(831)=1.399;q(831)=-90.6311059512
t(832)=1.4;q(832)=-98.4835449033
t(833)=1.401;q(833)=-105.904946659
t(834)=1.402;q(834)=-112.8223452821
t(835)=1.403;q(835)=-119.1627748364
t(836)=1.404;q(836)=-124.8532693858
t(837)=1.405;q(837)=-129.820862994
t(838)=1.406;q(838)=-133.9925897248
t(839)=1.407;q(839)=-137.2954836422
t(840)=1.408;q(840)=-139.6565788098
t(841)=1.409;q(841)=-141.0029092915
t(842)=1.41;q(842)=-141.2128903389
t(843)=1.411;q(843)=-139.3941327416
t(844)=1.412;q(844)=-135.5107556468
t(845)=1.413;q(845)=-129.9246387691
t(846)=1.414;q(846)=-122.9976618231
t(847)=1.415;q(847)=-115.0917045235
t(848)=1.416;q(848)=-106.5686465848
t(849)=1.417;q(849)=-97.7903677216
t(850)=1.418;q(850)=-89.1187476486
t(851)=1.419;q(851)=-80.9156660803
t(852)=1.42;q(852)=-73.5178158057
t(853)=1.421;q(853)=-65.5825856848
t(854)=1.422;q(854)=-56.5169011857
t(855)=1.423;q(855)=-46.7935106867
t(856)=1.424;q(856)=-36.8851625662
t(857)=1.425;q(857)=-27.2646052024
t(858)=1.426;q(858)=-18.4045869736
t(859)=1.427;q(859)=-10.7778562581
t(860)=1.428;q(860)=-4.8571614342
t(861)=1.429;q(861)=-1.1152508802
t(862)=1.43;q(862)=0.0
t(863)=2.00;q(863)=0.0


              !Initialize variables
              period = 1000
              m=1
              n=863
              !Compute derivation
              DO i=1,n-1
                delta(i)=(q(i+1)-q(i))/(t(i+1)-t(i))
              END DO
              delta(n)=delta(n-1)+(delta(n-1)-delta(n-2))/(t(n-1)-t(n-2))*(t(n)-t(n-1))
              !Find subinterval
              DO j=1,n-1
                IF(t(j) <= (TIME/period)) THEN
                  m=j
                END IF
              END DO
              !Evaluate interpolant
              s=(TIME/period)-t(m)
              VALUE=(q(m)+s*delta(m))
            ELSE
              CALL FlagError("Incorrect component specification for Olufsen flow rate waveform ",ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE) 
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE)
          ! Do nothing
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Olufsen flowrate waveform for "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
          & " dimension problem has not yet been implemented."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
      ! Returns a sinusoidal value for boundary nodes
      SELECT CASE(NUMBER_OF_DIMENSIONS)
      CASE(2,3)
        componentCoeff(1) = ANALYTIC_PARAMETERS(1) 
        componentCoeff(2) = ANALYTIC_PARAMETERS(2) 
        componentCoeff(3) = ANALYTIC_PARAMETERS(3) 
        componentCoeff(4) = ANALYTIC_PARAMETERS(4) 
        amplitude = ANALYTIC_PARAMETERS(5) 
        yOffset = ANALYTIC_PARAMETERS(6) 
        frequency = ANALYTIC_PARAMETERS(7) 
        phaseShift = ANALYTIC_PARAMETERS(8) 
        startTime = ANALYTIC_PARAMETERS(9) 
        stopTime = ANALYTIC_PARAMETERS(10) 
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(CURRENT_TIME > startTime - ZERO_TOLERANCE .AND. &
              &  CURRENT_TIME < stopTime + ZERO_TOLERANCE) THEN
              VALUE= componentCoeff(componentNumber)*(yOffset + amplitude*SIN(frequency*CURRENT_TIME+phaseShift))
            ELSE
              VALUE= componentCoeff(componentNumber)*(yOffset + amplitude*SIN(frequency*stopTime+phaseShift))
            END IF
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE) 
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Sinusoidal analytic types for "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
          & " dimensional problems have not yet been implemented."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1)
      IF(NUMBER_OF_DIMENSIONS==1.AND.NUMBER_OF_COMPONENTS==3) THEN
        !Polynomial function
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate Q
              VALUE=X(1)**2/10.0_DP**2
            ELSE IF(componentNumber==2) THEN
              !calculate A
              VALUE=X(1)**2/10.0_DP**2
            ELSE IF(componentNumber==3) THEN
              !calculate P
              VALUE=X(1)**2/10.0_DP**2
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE) 
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1)
      IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
        !Polynomial function
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=X(2)**2/10.0_DP**2
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=X(1)**2/10.0_DP**2
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=2.0_DP/3.0_DP*X(1)*(3.0_DP*MU_PARAM*10.0_DP**2-RHO_PARAM*X(1)**2*X(2))/(10.0_DP ** 4)
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2)
      IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
        !Exponential function
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE= EXP((X(1)-X(2))/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE= EXP((X(1)-X(2))/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE= 2.0_DP*MU_PARAM/10.0_DP*EXP((X(1)-X(2))/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE= 0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE= 0.0_DP
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE= 0.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3)
      IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
        !Sine and cosine function
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=COS(2.0_DP*PI*X(1)/10.0_DP)*COS(2.0_DP*PI*X(2)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=4.0_DP*MU_PARAM*PI/10.0_DP*SIN(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)+ &
                & 0.5_DP*RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_index,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=16.0_DP*MU_PARAM*PI**2/10.0_DP**2*COS(2.0_DP*PI*X(2)/ 10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4,EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5)
      IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
        !Taylor-Green vortex solution
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(X(1)/10.0_DP*2.0_DP*PI)*COS(X(2)/10.0_DP*2.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
              VALUE=SIN(X(1)/10.0_DP*PI)*COS(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
              !                      VALUE=SIN(X(1))*COS(X(2))
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=-COS(X(1)/10.0_DP*2.0_DP*PI)*SIN(X(2)/10.0_DP*2.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
              VALUE=-COS(X(1)/10.0_DP*PI)*SIN(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
              !                      VALUE=-COS(X(1))*SIN(X(2))
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1)/10.0_DP*2.0_DP*PI)+COS(2.0_DP*X(2)/10.0_DP*2.0_DP*PI))* &
                & EXP(-4.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)                      
              VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1)/10.0_DP*PI)+COS(2.0_DP*X(2)/10.0_DP*PI))* &
                & EXP(-4.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)                      
              !                      VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1))+COS(2.0_DP*X(2)))
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=0.0_DP         
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1)
      IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
        !Polynomial function
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=X(2)**2/10.0_DP**2+X(3)**2/10.0_DP**2
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=X(1)**2/10.0_DP**2+X(3)**2/10.0_DP** 2
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=X(1)**2/10.0_DP**2+X(2)**2/10.0_DP** 2
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=2.0_DP/3.0_DP*X(1)*(6.0_DP*MU_PARAM*10.0_DP**2-RHO_PARAM*X(2)*X(1)**2-3.0_DP* & 
                & RHO_PARAM*X(2)* &
                & X(3)**2-RHO_PARAM*X(3)*X(1)**2-3.0_DP*RHO_PARAM*X(3)*X(2)**2)/(10.0_DP**4)
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2)
      IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
        !Exponential function
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(3)-X(1))/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=EXP((X(3)-X(1))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=1.0_DP/10.0_DP*(2.0_DP*MU_PARAM*EXP((X(1)-X(2))/10.0_DP)- & 
                & 2.0_DP*MU_PARAM*EXP((X(3)-X(1))/10.0_DP)+RHO_PARAM*10.0_DP*EXP((X(1)-X(3))/10.0_DP)+ & 
                & RHO_PARAM*10.0_DP*EXP((X(2)-X(1))/10.0_DP))
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=-2.0_DP*MU_PARAM*(2.0_DP*EXP(X(1)-X(2))+EXP(X(2)-X(3)))
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=-2.0_DP*MU_PARAM*(2.0_DP*EXP(X(3)-X(1))+EXP(X(2)-X(3)))
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3)
      IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
        !Sine/cosine function
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(3)/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=2.0_DP*COS(2.0_DP*PI*x(1)/10.0_DP)*SIN(2.0_DP*PI*X(3)/10.0_DP)*COS(2.0_DP*PI*X(2)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=-COS(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(3)/10.0_DP)
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=-COS(2.0_DP*PI*X(1)/10.0_DP)*(-12.0_DP*MU_PARAM*PI*SIN(2.0_DP*PI*X(2)/10.0_DP)* & 
                & SIN(2.0_DP*PI*X(3)/10.0_DP)-RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*10.0_DP+ &
                & 2.0_DP*RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*X(3)/10.0_DP)**2- &
                & RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)**2)/10.0_DP/2.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=36*MU_PARAM*PI**2/10.0_DP**2*COS(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(3)/10.0_DP)* & 
                & COS(2.0_DP*PI*X(1)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=0.0_DP
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
      
    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4,EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
      IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
        !Taylor-Green vortex solution
        MU_PARAM = MATERIALS_PARAMETERS(1)
        RHO_PARAM = MATERIALS_PARAMETERS(2)
        SELECT CASE(variable_type)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(X(1)/10.0_DP*PI)*COS(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=-COS(X(1)/10.0_DP*PI)*SIN(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
            ELSE IF(componentNumber==3) THEN
              !calculate v
              VALUE=0.0_DP
              !                      VALUE=-COS(X(1))*SIN(X(2))
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1)/10.0_DP*PI)+COS(2.0_DP*X(2)/10.0_DP*PI))* &
                & EXP(-4.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)                      
              !                      VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1))+COS(2.0_DP*X(2)))
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIV_INDEX)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=0.0_DP         
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",ERR,ERROR,*999)
            END IF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)                                    
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
              & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE 
        LOCAL_ERROR="The number of components does not correspond to the number of dimensions."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END IF
    CASE DEFAULT
      LOCAL_ERROR="The analytic function type of "// &
        & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
        & " is invalid."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
    
    EXITS("NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Update SUPG parameters for Navier-Stokes equation
  SUBROUTINE NavierStokes_ResidualBasedStabilisation(equationsSet,elementNumber,gaussNumber,mu,rho,jacobianFlag,err,error,*)

    !Argument variables                               
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: elementNumber,gaussNumber
    REAL(DP), INTENT(IN) :: mu,rho
    LOGICAL, INTENT(IN) ::  jacobianFlag
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables                        
    TYPE(BASIS_TYPE), POINTER :: basisVelocity,basisPressure
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(FIELD_TYPE), POINTER :: equationsSetField,dependentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureVelocity,quadraturePressure
    INTEGER(INTG) :: fieldVariableType,meshComponent1,meshComponent2,numberOfDimensions,stabilisationType
    INTEGER(INTG) :: i,j,k,l,mhs,nhs,ms,ns,nh,mh,nj,ni,pressureIndex,numberOfElementParameters(4)
    REAL(DP) :: PHIMS,PHINS,DXI_DX(3,3),tauC,tauMp,tauMu
    REAL(DP) :: dPhi_dX_Velocity(27,3),dPhi_dX_Pressure(27,3),DPHINS2_DXI(3,3)
    REAL(DP) :: jacobianMomentum(3),jacobianContinuity,residualMomentum(3),residualContinuity
    REAL(DP) :: velocity(3),velocityPrevious(3),velocityDeriv(3,3),velocity2Deriv(3,3,3),pressure,pressureDeriv(3)
    REAL(DP) :: JGW,SUM,SUM2,SUPG,PSPG,LSIC,crossStress,reynoldsStress,momentumTerm
    REAL(DP) :: uDotGu,doubleDotG,tauSUPS,traceG,nuLSIC,timeIncrement,elementInverse,C1,stabilisationValueDP
    LOGICAL :: linearElement
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_ResidualBasedStabilisation",err,error,*999)

    ! Nullify all local pointers
    NULLIFY(basisVelocity)
    NULLIFY(basisPressure)
    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(nonlinearMapping)
    NULLIFY(equationsEquationsSetField)
    NULLIFY(equationsSetField)
    NULLIFY(quadratureVelocity)
    NULLIFY(quadraturePressure)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(fieldVariable)
    NULLIFY(equationsMatrices)
    NULLIFY(nonlinearMatrices)
    NULLIFY(jacobianMatrix)

    IF(ASSOCIATED(equationsSet))THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        equations=>equationsSet%EQUATIONS
        IF (ASSOCIATED(equations)) THEN
          !Set general and specific pointers
          equationsMapping=>equations%EQUATIONS_MAPPING
          equationsMatrices=>equations%EQUATIONS_MATRICES
          nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
          nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
          jacobianMatrix=>nonlinearMatrices%JACOBIANS(1)%PTR
          fieldVariable=>nonlinearMapping%RESIDUAL_VARIABLES(1)%PTR
          fieldVariableType=fieldVariable%VARIABLE_TYPE
          geometricField=>equations%INTERPOLATION%GEOMETRIC_FIELD
          numberOfDimensions=fieldVariable%NUMBER_OF_COMPONENTS - 1
          equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
          meshComponent1=fieldVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
          meshComponent2=fieldVariable%COMPONENTS(fieldVariable%NUMBER_OF_COMPONENTS)%MESH_COMPONENT_NUMBER
          dependentField=>equations%INTERPOLATION%DEPENDENT_FIELD
          basisVelocity=>dependentField%DECOMPOSITION%DOMAIN(meshComponent1)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
          basisPressure=>dependentField%DECOMPOSITION%DOMAIN(meshComponent2)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS

          IF(basisVelocity%INTERPOLATION_ORDER(1).LE.1) THEN
            linearElement = .TRUE.
          ELSE
            ! higher order element type- can calculate 2nd order terms
            linearElement = .FALSE.
          END IF

          quadratureVelocity=>basisVelocity%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          quadraturePressure=>basisPressure%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
          IF(ASSOCIATED(equationsSetField)) THEN
            ! Stabilisation type (default 1 for RBS, 2 for RBVM, 0 for none)
            CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
             & 4,stabilisationValueDP,err,error,*999)
            stabilisationType=NINT(stabilisationValueDP)
            ! Skip if type 0
            IF(stabilisationType > 0) THEN
              ! Get time step size and calc time derivative
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & 3,timeIncrement,err,error,*999)
              ! TODO: put this somewhere more sensible. This is a workaround since we don't have access to the dynamic solver values
              !       at this level in the element loop
              IF(equationsSet%specification(3)/=EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE &
                & .AND. timeIncrement < ZERO_TOLERANCE) THEN
                CALL FlagError("Please set the equations set field time increment to a value > 0.",ERR,ERROR,*999)                
              END IF
              ! Stabilisation type (default 1 for RBS)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & 4,stabilisationValueDP,err,error,*999)
              stabilisationType=NINT(stabilisationValueDP)
              ! User specified or previously calculated C1
              CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & elementNumber,10,elementInverse,err,error,*999)

              ! Get previous timestep values
              velocityPrevious=0.0_DP
              IF(equationsSet%specification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,elementNumber,equations% &
                 & INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,EQUATIONS%INTERPOLATION% &
                 & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                velocityPrevious=0.0_DP
                DO i=1,numberOfDimensions
                  velocityPrevious(i)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR%VALUES(i,NO_PART_DERIV)
                END DO
              END IF

              ! Interpolate current solution velocity/pressure field values
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,EQUATIONS%INTERPOLATION% &
               & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)              
              IF(linearElement) THEN
                ! Get 1st order derivatives for current timestep value
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,EQUATIONS%INTERPOLATION% &
                 & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              ELSE
                ! Get 2nd order derivatives for current timestep value
                CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,EQUATIONS%INTERPOLATION%&
                 & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              END IF
              velocity=0.0_DP
              velocityDeriv=0.0_DP
              velocity2Deriv=0.0_DP
              pressure=0.0_DP
              pressureDeriv=0.0_DP
              DO i=1,numberOfDimensions
                velocity(i)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR%VALUES(i,NO_PART_DERIV)
                velocityDeriv(i,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR%VALUES(i,PART_DERIV_S1)
                velocityDeriv(i,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR%VALUES(i,PART_DERIV_S2)
                IF(.NOT. linearElement) THEN
                  velocity2Deriv(i,1,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                   & VALUES(i,PART_DERIV_S1_S1)
                  velocity2Deriv(i,1,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                   & VALUES(i,PART_DERIV_S1_S2)
                  velocity2Deriv(i,2,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                   & VALUES(i,PART_DERIV_S1_S2)
                  velocity2Deriv(i,2,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                   & VALUES(i,PART_DERIV_S2_S2)
                END IF
                IF(numberOfDimensions > 2) THEN
                  velocityDeriv(i,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR%VALUES(i,PART_DERIV_S3)
                  IF(.NOT. linearElement) THEN                  
                    velocity2Deriv(i,1,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                     & VALUES(i,PART_DERIV_S1_S3)
                    velocity2Deriv(i,2,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                     & VALUES(i,PART_DERIV_S2_S3)
                    velocity2Deriv(i,3,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                     & VALUES(i,PART_DERIV_S1_S3)
                    velocity2Deriv(i,3,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                     & VALUES(i,PART_DERIV_S2_S3)
                    velocity2Deriv(i,3,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR% &
                     & VALUES(i,PART_DERIV_S3_S3)
                  END IF
                END IF
              END DO
              pressureIndex = numberOfDimensions + 1
              pressure=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)%PTR%VALUES(pressureIndex,NO_PART_DERIV)
              pressureDeriv(1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)% &
               & PTR%VALUES(pressureIndex,PART_DERIV_S1)
              pressureDeriv(2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)% &
               & PTR%VALUES(pressureIndex,PART_DERIV_S2)
              IF(numberOfDimensions > 2) THEN
                pressureDeriv(3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(fieldVariableType)% &
                 & PTR%VALUES(pressureIndex,PART_DERIV_S3)
              END IF
              DXI_DX=0.0_DP
              DO i=1,numberOfDimensions
                DO j=1,numberOfDimensions
                  DXI_DX(j,i)=equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                    & DXI_DX(j,i)
                END DO
              END DO

              ! Get number of element parameters for each dependent component
              numberOfElementParameters=0
              DO i=1,numberOfDimensions
                numberOfElementParameters(i)=basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS              
              END DO
              numberOfElementParameters(numberOfDimensions+1)=basisPressure%NUMBER_OF_ELEMENT_PARAMETERS              
              ! Calculate dPhi/dX
              dPhi_dX_Velocity=0.0_DP
              dPhi_dX_Pressure=0.0_DP
              DO ms=1,numberOfElementParameters(1)
                DO nj=1,numberOfDimensions
                  dPhi_dX_Velocity(ms,nj)=0.0_DP
                  DO ni=1,numberOfDimensions
                    dPhi_dX_Velocity(ms,nj)=dPhi_dX_Velocity(ms,nj) + &
                     & quadratureVelocity%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),gaussNumber)* &
                     & DXI_DX(ni,nj)
                  END DO
                END DO
              END DO
              DO ms=1,numberOfElementParameters(numberOfDimensions+1)
                DO nj=1,numberOfDimensions
                  dPhi_dX_Pressure(ms,nj)=0.0_DP
                  DO ni=1,numberOfDimensions
                    dPhi_dX_Pressure(ms,nj)=dPhi_dX_Pressure(ms,nj) + &
                     & quadraturePressure%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),gaussNumber)* &
                     & DXI_DX(ni,nj)
                  END DO
                END DO
              END DO
              JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
               & quadratureVelocity%GAUSS_WEIGHTS(gaussNumber)

              !----------------------------------------------------------------------------------
              ! C a l c u l a t e   d i s c r e t e   r e s i d u a l s 
              !----------------------------------------------------------------------------------
              SUM = 0.0_DP
              residualMomentum = 0.0_DP
              residualContinuity = 0.0_DP
              ! Calculate momentum residual
              DO i=1,numberOfDimensions    
                SUM = 0.0_DP
                ! velocity time derivative              
                IF(equationsSet%specification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  SUM = rho*(velocity(i)-velocityPrevious(i))/timeIncrement
                END IF
                DO j=1,numberOfDimensions
                  ! pressure gradient
                  SUM = SUM + pressureDeriv(j)*DXI_DX(j,i)
                  DO k=1,numberOfDimensions
                    !Convective term
                    SUM = SUM +rho*((velocity(j))*(velocityDeriv(i,k)*DXI_DX(k,j)))
                    IF(.NOT. linearElement) THEN
                      DO l=1,numberOfDimensions
                        ! viscous stress: only if quadratic or higher basis defined for laplacian
                        SUM = SUM - mu*(velocity2Deriv(i,k,l)*DXI_DX(k,j)*DXI_DX(l,j))
                      END DO
                    END IF
                  END DO
                END DO
                residualMomentum(i) = SUM
              END DO
              ! Calculate continuity residual
              SUM = 0.0_DP
              DO i=1,numberOfDimensions
                DO j=1,numberOfDimensions
                  SUM= SUM + velocityDeriv(i,j)*DXI_DX(j,i)
                END DO
              END DO
              residualContinuity = SUM

              ! Constant of element inverse inequality
              IF(elementInverse > -ZERO_TOLERANCE) THEN
                ! Use user-defined value if specified (default -1)
                C1 = elementInverse
              ELSE IF(linearElement) THEN
                C1=3.0_DP
              ELSE
                IF(numberOfDimensions==2 .AND. basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS==9 &
                  & .AND. basisVelocity%INTERPOLATION_ORDER(1)==2) THEN
                  C1=24.0_DP
                ELSE IF(numberOfDimensions==3 .AND. basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS==27 &
                  & .AND. basisVelocity%INTERPOLATION_ORDER(1)==2) THEN
                  C1=12.0_DP
                !TODO: Expand C1 for more element types
                ELSE
                  CALL FlagError("Element inverse estimate undefined on element " &
                   & //TRIM(NUMBER_TO_VSTRING(elementNumber,"*",err,error)),err,error,*999)                              
                END IF
              END IF
              ! Update element inverse value if calculated
              IF(ABS(C1-elementInverse) > ZERO_TOLERANCE) THEN
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & elementNumber,10,C1,err,error,*999)
              END IF

              !----------------------------------------------------------
              ! S t a b i l i z a t i o n    C o n s t a n t s    (Taus)
              !----------------------------------------------------------
              IF(stabilisationType == 1 .OR. stabilisationType == 2) THEN
                ! Bazilevs method for calculating tau
                pointMetrics=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
                uDotGu = 0.0_DP
                DO i=1,numberOfDimensions
                  DO j=1,numberOfDimensions
                    uDotGu = uDotGu + velocity(i)*pointMetrics%GU(i,j)*velocity(j)
                  END DO
                END DO
                doubleDotG = 0.0_DP
                DO i=1,numberOfDimensions
                  DO j=1,numberOfDimensions
                    doubleDotG = doubleDotG + pointMetrics%GU(i,j)*pointMetrics%GU(i,j)
                  END DO
                END DO
                ! Calculate tauSUPS (used for both PSPG and SUPG weights)
                IF(equationsSet%specification(3) == EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  tauSUPS = (uDotGu + (C1*((mu/rho)**2.0_DP)*doubleDotG))**(-0.5_DP)
                ELSE
                  tauSUPS = ((4.0_DP/(timeIncrement**2.0_DP)) + uDotGu + (C1*((mu/rho)**2.0_DP)*doubleDotG))**(-0.5_DP)
                END IF

                ! Calculate nu_LSIC (Least-squares incompressibility constraint)
                traceG = 0.0_DP
                DO i=1,numberOfDimensions
                  traceG = traceG + pointMetrics%GU(i,i)
                END DO
                nuLSIC = 1.0_DP/(tauSUPS*traceG)            

                tauMp = tauSUPS
                tauMu = tauSUPS
                tauC = nuLSIC

              ELSE 
                CALL FlagError("A tau factor has not been defined for the stabilisation type of " &
                 & //TRIM(NUMBER_TO_VSTRING(stabilisationType,"*",err,error)),err,error,*999)
              END IF

              !-------------------------------------------------------------------------------------------------
              ! A d d   s t a b i l i z a t i o n   f a c t o r s   t o   e l e m e n t   m a t r i c e s
              !-------------------------------------------------------------------------------------------------
              jacobianMomentum = 0.0_DP
              jacobianContinuity = 0.0_DP
              mhs = 0
              DO mh=1,numberOfDimensions+1
                DO ms=1,numberOfElementParameters(mh)
                  mhs = mhs + 1
                  IF(mh <= numberOfDimensions) THEN
                    PHIMS=quadratureVelocity%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussNumber)
                    JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                      & quadratureVelocity%GAUSS_WEIGHTS(gaussNumber)
                  ELSE
                    PHIMS=quadraturePressure%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussNumber)
                    JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                      & quadraturePressure%GAUSS_WEIGHTS(gaussNumber)
                  END IF
                  !------------------
                  ! J A C O B I A N 
                  !------------------
                  IF(jacobianFlag) THEN
                    nhs = 0
                    DO nh=1,numberOfDimensions+1
                      DO ns=1,numberOfElementParameters(nh)
                        nhs=nhs+1
                        ! Note that we still need to assemble the vector momentum jacobian for PSPG in the continuity row
                        IF(nh <= numberOfDimensions) THEN
                          PHINS=quadratureVelocity%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,gaussNumber)
                        ELSE
                          PHINS=quadraturePressure%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,gaussNumber)
                        END IF

                        ! Calculate jacobians of the discrete residual terms
                        jacobianMomentum = 0.0_DP
                        IF(nh == numberOfDimensions+1) THEN
                          ! d(Momentum(mh))/d(Pressure)
                          DO i=1,numberOfDimensions
                            jacobianMomentum(i) = dPhi_dX_Pressure(ns,i)
                          END DO
                          jacobianContinuity=0.0_DP
                        ELSE
                          DPHINS2_DXI=0.0_DP
                          IF(.NOT. linearElement) THEN
                            DPHINS2_DXI(1,1)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,gaussNumber)
                            DPHINS2_DXI(1,2)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,gaussNumber)
                            DPHINS2_DXI(2,1)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,gaussNumber)
                            DPHINS2_DXI(2,2)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,gaussNumber)
                            IF(numberOfDimensions > 2) THEN
                              DPHINS2_DXI(1,3)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,gaussNumber)
                              DPHINS2_DXI(2,3)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,gaussNumber)
                              DPHINS2_DXI(3,1)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,gaussNumber)
                              DPHINS2_DXI(3,2)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,gaussNumber)
                              DPHINS2_DXI(3,3)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,gaussNumber)
                            END IF
                          END IF
                          ! d(Momentum)/d(Velocity(nh))
                          jacobianMomentum = 0.0_DP
                          DO i=1,numberOfDimensions 
                            SUM = 0.0_DP
                            !Note: Convective term split using product rule
                            !Convective term 1: applies to all velocity components
                            DO j=1,numberOfDimensions
                              SUM = SUM + rho*PHINS*velocityDeriv(i,j)*DXI_DX(j,nh)
                            END DO
                            !Diagonal terms
                            IF(i==nh) THEN
                              !Transient
                              IF(equationsSet%specification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
                                SUM = SUM + rho*PHINS/timeIncrement
                              END IF
                              !Convective 2: nh component only
                              DO j=1,numberOfDimensions
                                SUM = SUM + rho*velocity(j)*dPhi_dX_Velocity(ns,j)
                              END DO
                              IF(.NOT. linearElement) THEN
                                !Viscous laplacian term
                                DO j=1,numberOfDimensions
                                  DO k=1,numberOfDimensions
                                    DO l=1,numberOfDimensions
                                      SUM=SUM-mu*DPHINS2_DXI(k,l)*DXI_DX(k,j)*DXI_DX(l,j)
                                    END DO 
                                  END DO 
                                END DO 
                              END IF
                            END IF
                            jacobianMomentum(i)=SUM
                          END DO
                          ! Continuity/velocity
                          jacobianContinuity = dPhi_dX_Velocity(ns,nh)
                        END IF
                        ! Calculate jacobian of discrete residual * RBS factors (apply product rule if neccesary)

                        ! PSPG: Pressure stabilising Petrov-Galerkin
                        IF(mh == numberOfDimensions+1) THEN
                          PSPG = 0.0_DP
                          SUM = 0.0_DP
                          DO i=1,numberOfDimensions
                            SUM = SUM + dPhi_dX_Pressure(ms,i)*jacobianMomentum(i)
                          END DO
                          PSPG = tauMp*SUM/rho*JGW

                          jacobianMatrix%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=jacobianMatrix%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+PSPG

                        ! SUPG: Streamline upwind/Petrov-Galerkin
                        ! LSIC: Least-squares incompressibility constraint
                        ELSE
                          SUPG=0.0_DP
                          LSIC=0.0_DP

                          SUM=0.0_DP
                          IF(nh <= numberOfDimensions) THEN
                            SUPG= SUPG + PHINS*dPhi_dX_Velocity(ms,nh)*residualMomentum(mh)                            
                          END IF
                          DO i=1,numberOfDimensions
                            SUM = SUM + velocity(i)*dPhi_dX_Velocity(ms,i)
                          END DO
                          SUPG = tauMu*(SUPG + SUM*jacobianMomentum(mh))

                          SUM=0.0_DP
                          DO i=1,numberOfDimensions
                            SUM = SUM + dPhi_dX_Velocity(ms,i)
                          END DO
                          LSIC = tauC*rho*dPhi_dX_Velocity(ms,mh)*jacobianContinuity                        

                          momentumTerm = (SUPG + LSIC)*JGW

                          IF(stabilisationType == 2) THEN
                            ! Additional terms for RBVM
                            crossStress=0.0_DP
                            reynoldsStress=0.0_DP
                            crossStress = 0.0_DP
                            IF(nh <= numberOfDimensions) THEN
                              IF(mh == nh) THEN
                                DO i=1,numberOfDimensions
                                  crossStress= crossStress + dPhi_dX_Velocity(ns,i)*residualMomentum(i)                            
                                END DO
                              END IF
                            END IF
                            SUM2=0.0_DP
                            DO i=1,numberOfDimensions
                              SUM=0.0_DP
                              ! dU_mh/dX_i
                              DO j=1,numberOfDimensions                        
                                SUM= SUM + velocityDeriv(mh,j)*DXI_DX(j,i)
                              END DO
                              ! Jm_i*dU_mh/dX_i
                              SUM2 = SUM2 + jacobianMomentum(i)*SUM
                            END DO
                            crossStress = -tauMu*(crossStress + SUM2)                          

                            reynoldsStress = 0.0_DP
                            SUM = 0.0_DP
                            !Rm_mh.Rm_i.dPhi/dX_i
                            DO i=1,numberOfDimensions
                              SUM = SUM + jacobianMomentum(mh)*residualMomentum(i)*dPhi_DX_Velocity(ms,i)
                              SUM = SUM + jacobianMomentum(i)*residualMomentum(mh)*dPhi_DX_Velocity(ms,i)
                            END DO
                            reynoldsStress = -tauMu*tauMu*SUM

                            momentumTerm = momentumTerm + (crossStress + reynoldsStress)*JGW
                          END IF
                          
                          ! Add stabilisation to element jacobian
                          jacobianMatrix%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)= &
                           & jacobianMatrix%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+momentumTerm

                        END IF
                      END DO
                    END DO                    

                  !-----------------
                  ! R E S I D U A L
                  !-----------------
                  ELSE
                    ! PSPG: Pressure stabilising Petrov-Galerkin
                    IF(mh == numberOfDimensions+1) THEN
                      SUM = 0.0_DP
                      DO i=1,numberOfDimensions
                        SUM = SUM + dPhi_dX_Pressure(ms,i)*residualMomentum(i)
                      END DO
                      PSPG = SUM*(tauMp/rho)*JGW
                      nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR(mhs)= &
                       & nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR(mhs) + PSPG

                    ! SUPG: Streamline upwind/Petrov-Galerkin
                    ! LSIC: Least-squares incompressibility constraint
                    ELSE
                      SUPG=0.0_DP
                      LSIC=0.0_DP

                      ! u_i*Rm_mh*dv_mh/dx_i
                      SUM=0.0_DP
                      DO i=1,numberOfDimensions
                        SUM = SUM + velocity(i)*dPhi_dX_Velocity(ms,i)
                      END DO
                      SUPG = tauMu*SUM*residualMomentum(mh)

                      LSIC = tauC*rho*dPhi_dX_Velocity(ms,mh)*residualContinuity                        
                      momentumTerm = (SUPG + LSIC)*JGW

                      IF(stabilisationType ==2) THEN
                        ! Additional terms for RBVM
                        crossStress=0.0_DP
                        reynoldsStress=0.0_DP
                        SUM2 = 0.0_DP
                        DO i=1,numberOfDimensions                        
                          SUM = 0.0_DP
                          ! dU_mh/dX_i
                          DO j=1,numberOfDimensions                        
                            SUM= SUM + velocityDeriv(mh,j)*DXI_DX(j,i)
                          END DO
                          ! Rm_i.dU_mh/dX_i
                          SUM2= SUM2 + residualMomentum(i)*SUM
                        END DO
                        crossStress= -tauMu*PHIMS*SUM2

                        reynoldsStress = 0.0_DP
                        SUM = 0.0_DP
                        !Rm_mh.Rm_i.dPhi/dX_i
                        DO i=1,numberOfDimensions
                          SUM = SUM + dPhi_dX_Velocity(ms,i)*residualMomentum(i)*residualMomentum(mh)
                        END DO
                        reynoldsStress = -SUM*(tauMu*tauMu)/rho
                        momentumTerm = momentumTerm + (crossStress + reynoldsStress)*JGW
                      END IF

                      ! Add stabilisation to element residual
                      nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR(mhs)= &
                       & nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR(mhs) + momentumTerm
                    END IF
                  END IF ! jacobian/residual
                END DO !ms
              END DO !mh

            END IF ! check stabilisation type               
          ELSE
            CALL FlagError("Equations equations set field is not associated.",err,error,*999)
          END IF               
        ELSE
          CALL FlagError("Equations set equations is not associated.",err,error,*999)
        END IF               
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
          & " is not a valid subtype to use SUPG weighting functions."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF
    
    EXITS("NavierStokes_ResidualBasedStabilisation")
    RETURN
999 ERRORSEXITS("NavierStokes_ResidualBasedStabilisation",err,error)
    RETURN 1
    
  END SUBROUTINE NavierStokes_ResidualBasedStabilisation

  !
  !================================================================================================================================
  !

  !>Calculate element-level scale factors: CFL, cell Reynolds number
  SUBROUTINE NavierStokes_CalculateElementMetrics(equationsSet,elementNumber,err,error,*)

    !Argument variables                               
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: elementNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables                        
    TYPE(BASIS_TYPE), POINTER :: basisVelocity
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureVelocity
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    INTEGER(INTG) :: fieldVariableType,meshComponent1
    INTEGER(INTG) :: numberOfDimensions,mh
    INTEGER(INTG) :: gaussNumber
    INTEGER(INTG) :: i,j,ms
    INTEGER(INTG) :: numberOfElementParameters
    INTEGER(INTG) :: LWORK,INFO
    REAL(DP) :: cellReynoldsNumber,cellCourantNumber,timeIncrement
    REAL(DP) :: dPhi_dX_Velocity(27,3)
    REAL(DP) :: DXI_DX(3,3)
    REAL(DP) :: velocity(3),avgVelocity(3),velocityNorm,velocityPrevious(3),velocityDeriv(3,3)
    REAL(DP) :: PHIMS,JGW,SUM,SUM2,mu,rho,normCMatrix,normKMatrix,normMMatrix,muScale
    REAL(DP) :: CMatrix(27,3),KMatrix(27,3),MMatrix(27,3)
    REAL(DP) :: svd(3),U(27,27),VT(3,3)
    REAL(DP), ALLOCATABLE :: WORK(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_CalculateElementMetrics",err,error,*999)

    ! Nullify all local pointers
    NULLIFY(basisVelocity)
    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(nonlinearMapping)
    NULLIFY(equationsEquationsSetField)
    NULLIFY(equationsSetField)
    NULLIFY(quadratureVelocity)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(fieldVariable)

    IF(ASSOCIATED(equationsSet))THEN
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
        equations=>equationsSet%EQUATIONS
        IF(ASSOCIATED(equations)) THEN
          !Set general and specific pointers
          equationsMapping=>equations%EQUATIONS_MAPPING
          nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
          fieldVariable=>nonlinearMapping%RESIDUAL_VARIABLES(1)%PTR
          fieldVariableType=fieldVariable%VARIABLE_TYPE
          geometricField=>equations%INTERPOLATION%GEOMETRIC_FIELD
          numberOfDimensions=fieldVariable%NUMBER_OF_COMPONENTS - 1
          equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
          meshComponent1=fieldVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
          dependentField=>equations%INTERPOLATION%DEPENDENT_FIELD
          basisVelocity=>dependentField%DECOMPOSITION%DOMAIN(meshComponent1)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
          quadratureVelocity=>basisVelocity%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          IF(ASSOCIATED(equationsEquationsSetField)) THEN
            equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
            IF (ASSOCIATED(equationsSetField)) THEN

              ! Get time step size
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & 3,timeIncrement,err,error,*999)

              ! Loop over gauss points
              CMatrix = 0.0_DP
              MMatrix = 0.0_DP
              KMatrix = 0.0_DP
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,equations%INTERPOLATION% &
                & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)

              avgVelocity = 0.0_DP
              DO gaussNumber = 1,quadratureVelocity%NUMBER_OF_GAUSS

                ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
                IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                  ! Note the constant from the U_VARIABLE is a scale factor
                  muScale = mu
                  ! Get the gauss point based value returned from the CellML solver
                  CALL Field_ParameterSetGetLocalGaussPoint(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,gaussNumber,elementNumber,1,mu,err,error,*999)
                  mu=mu*muScale
                END IF
   
                ! Get previous timestep values
                velocityPrevious=0.0_DP
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,elementNumber,equations% &
                 & INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,EQUATIONS%INTERPOLATION%&
                 & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                velocityPrevious=0.0_DP
                DO i=1,numberOfDimensions
                  velocityPrevious(i)=equations%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR% &
                   & VALUES(i,NO_PART_DERIV)
                END DO

                ! Interpolate current solution velocity and first deriv field values
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
                 & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)              
                ! Get 1st order derivatives for current timestep value
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber, &
                 & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                velocity=0.0_DP
                velocityDeriv=0.0_DP
                DO i=1,numberOfDimensions
                  velocity(i)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(i,NO_PART_DERIV)
                  velocityDeriv(i,1)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                   & PTR%VALUES(i,PART_DERIV_S1)
                  velocityDeriv(i,2)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                   & PTR%VALUES(i,PART_DERIV_S2)
                  IF(numberOfDimensions > 2) THEN
                    velocityDeriv(i,3)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                     & PTR%VALUES(i,PART_DERIV_S3)
                  END IF
                END DO

                ! get dXi/dX deriv
                CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%INTERPOLATION%&
                & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                DXI_DX=0.0_DP
                DO i=1,numberOfDimensions
                  DO j=1,numberOfDimensions
                    DXI_DX(j,i)=equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                      & DXI_DX(j,i)
                  END DO
                END DO

                numberOfElementParameters=basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS              
                ! Calculate dPhi/dX
                dPhi_dX_Velocity=0.0_DP
                DO ms=1,numberOfElementParameters
                  DO i=1,numberOfDimensions
                    dPhi_dX_Velocity(ms,i)=0.0_DP
                    DO j=1,numberOfDimensions
                      dPhi_dX_Velocity(ms,i)=dPhi_dX_Velocity(ms,i) + &
                       & quadratureVelocity%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(j),gaussNumber)* &
                       & DXI_DX(j,i)
                    END DO
                  END DO
                END DO

                JGW=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
                 & quadratureVelocity%GAUSS_WEIGHTS(gaussNumber)
                DO mh=1,numberOfDimensions
                  DO ms=1,numberOfElementParameters
                    PHIMS=quadratureVelocity%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussNumber)

                    ! c_(a,i)
                    SUM=0.0_DP
                    DO i=1,numberOfDimensions
                      DO j=1,numberOfDimensions                        
                        SUM = SUM + velocity(i)*velocityDeriv(mh,j)*DXI_DX(j,i)
                      END DO
                    END DO
                    CMatrix(ms,mh)=CMatrix(ms,mh) + rho*PHIMS*SUM*JGW

                    ! ~k_(a,i)
                    SUM=0.0_DP
                    DO i=1,numberOfDimensions
                      SUM = SUM + velocity(i)*dPhi_dX_Velocity(ms,i)
                    END DO
                    SUM2=0.0_DP
                    DO i=1,numberOfDimensions
                      DO j=1,numberOfDimensions                        
                        SUM2 = SUM2 + velocity(i)*velocityDeriv(mh,j)*DXI_DX(j,i)
                      END DO
                    END DO
                    KMatrix(ms,mh)=KMatrix(ms,mh)+rho*SUM*SUM2*JGW

                    ! m_(a,i)
                    MMatrix(ms,mh)=MMatrix(ms,mh)+rho*PHIMS*(velocity(mh)-velocityPrevious(mh))/timeIncrement*JGW

                  END DO !ms
                END DO !mh

                avgVelocity= avgVelocity + velocity/quadratureVelocity%NUMBER_OF_GAUSS
              END DO ! gauss loop

              LWORK=MAX(1,3*MIN(numberOfElementParameters,numberOfDimensions)+ &
               & MAX(numberOfElementParameters,numberOfDimensions),5*MIN(numberOfElementParameters,numberOfDimensions))
              ALLOCATE(WORK(LWORK))

              ! compute the singular value decomposition (SVD) using LAPACK
              CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,CMatrix,numberOfElementParameters,svd, &
               & U,numberOfElementParameters,VT,numberOfDimensions,WORK,LWORK,INFO)
              normCMatrix=svd(1)
              IF(INFO /= 0) THEN
                localError="Error calculating SVD on element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END IF

              CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,KMatrix,numberOfElementParameters,svd, &
               & U,numberOfElementParameters,VT,numberOfDimensions,WORK,LWORK,INFO)
              normKMatrix=svd(1)
              IF(INFO /= 0) THEN
                localError="Error calculating SVD on element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END IF

              CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,MMatrix,numberOfElementParameters,svd, &
               & U,numberOfElementParameters,VT,numberOfDimensions,WORK,LWORK,INFO)
              normMMatrix=svd(1)
              IF(INFO /= 0) THEN
                localError="Error calculating SVD on element "//TRIM(NUMBER_TO_VSTRING(elementNumber,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END IF
              DEALLOCATE(WORK)

              velocityNorm = L2NORM(avgVelocity)
              cellReynoldsNumber = 0.0_DP
              cellCourantNumber = 0.0_DP
              IF(velocityNorm > ZERO_TOLERANCE) THEN
                IF(normKMatrix > ZERO_TOLERANCE) THEN
                  cellReynoldsNumber = velocityNorm**2.0_DP/(mu/rho)*normCMatrix/normKMatrix
                END IF
                IF(normMMatrix > ZERO_TOLERANCE) THEN
                  cellCourantNumber = timeIncrement/2.0_DP*normCMatrix/normMMatrix
                END IF
              END IF
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & elementNumber,2,velocityNorm,err,error,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & elementNumber,3,cellCourantNumber,err,error,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & elementNumber,4,cellReynoldsNumber,err,error,*999)

            ELSE
              CALL FlagError("Equations set field field is not associated.",err,error,*999)
            END IF               
          ELSE
            CALL FlagError("Equations equations set field is not associated.",err,error,*999)
          END IF               
        ELSE
          CALL FlagError("Equations set equations is not associated.",err,error,*999)
        END IF               
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
          & " is not a valid subtype to use SUPG weighting functions."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF
    
    EXITS("NavierStokes_CalculateElementMetrics")
    RETURN
999 ERRORSEXITS("NavierStokes_CalculateElementMetrics",err,error)
    RETURN 1
    
  END SUBROUTINE NavierStokes_CalculateElementMetrics

  !
  !================================================================================================================================
  !

  !>Calculates the face integration term of the finite element formulation for Navier-Stokes equation,
  !>required for pressure and multidomain boundary conditions. 
  !>portions based on DarcyEquation_FiniteElementFaceIntegrate by Adam Reeve.
  SUBROUTINE NavierStokes_FiniteElementFaceIntegrate(equationsSet,elementNumber,dependentVariable,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: elementNumber
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local variables
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: decompElement
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: face
    TYPE(BASIS_TYPE), POINTER :: faceBasis
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: faceQuadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    INTEGER(INTG) :: faceIdx, faceNumber
    INTEGER(INTG) :: compIdx, gaussIdx
    INTEGER(INTG) :: elementBaseDofIdx, faceNodeIdx, elementNodeIdx
    INTEGER(INTG) :: faceNodederivIdx, meshComponentNumber, nodederivIdx,elementParameterIdx
    INTEGER(INTG) :: faceParameterIdx,elementDof,normalcompIdx
    INTEGER(INTG) :: numberOfDimensions,boundaryType
    REAL(DP) :: pressure,density,jacobianGaussWeights,beta,normalFlow
    REAL(DP) :: velocity(3),normalProjection(3),unitNormal(3),stabilisationTerm(3),boundaryNormal(3)
    REAL(DP) :: boundaryValue,normalDifference,normalTolerance,boundaryPressure
    REAL(DP) :: dUDXi(3,3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    LOGICAL :: integratedBoundary

    REAL(DP), POINTER :: geometricParameters(:)

    ENTERS("NavierStokes_FiniteElementFaceIntegrate",err,error,*999)

    NULLIFY(decomposition)
    NULLIFY(decompElement)
    NULLIFY(dependentBasis)
    NULLIFY(geometricVariable)
    NULLIFY(geometricParameters)
    NULLIFY(equations)
    NULLIFY(equationsSetField)
    NULLIFY(equationsEquationsSetField)
    NULLIFY(equationsMatrices)
    NULLIFY(face)
    NULLIFY(faceBasis)
    NULLIFY(faceQuadratureScheme)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(geometricInterpolationParameters)
    NULLIFY(rhsVector)
    NULLIFY(nonlinearMatrices)
    NULLIFY(dependentField)
    NULLIFY(geometricField)

    ! Get pointers and perform sanity checks
    IF(ASSOCIATED(equationsSet)) THEN
      dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
      IF(.NOT.ASSOCIATED(dependentField)) THEN
        CALL FlagError("Dependent field is not associated.",err,error,*999)
      END IF
      equations=>equationsSet%EQUATIONS
      IF (ASSOCIATED(equations)) THEN
        equationsMatrices=>equations%EQUATIONS_MATRICES
        IF (ASSOCIATED(equationsMatrices)) THEN
          rhsVector=>equationsMatrices%RHS_VECTOR
          nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
          IF(.NOT. ASSOCIATED(nonlinearMatrices)) THEN
            CALL FlagError("Nonlinear Matrices not associated.",err,error,*999)
          END IF
        END IF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)

      !Check for the equations set field
      equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
      IF(ASSOCIATED(equationsEquationsSetField)) THEN
        equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
        IF(.NOT.ASSOCIATED(equationsSetField)) THEN
          CALL FlagError("Equations set field (EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Equations set field (EQUATIONS_EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
      END IF

      ! Check whether this element contains an integrated boundary type
      CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementNumber,9,boundaryValue,err,error,*999)
      boundaryType=NINT(boundaryValue)
      integratedBoundary = .FALSE.
      IF(boundaryType == BOUNDARY_CONDITION_PRESSURE) integratedBoundary = .TRUE.

      !Get the mesh decomposition and basis
      decomposition=>dependentVariable%FIELD%DECOMPOSITION
      !Check that face geometric parameters have been calculated
      IF(decomposition%CALCULATE_FACES .AND. integratedBoundary) THEN
        meshComponentNumber=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        dependentBasis=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
          & ELEMENTS(elementNumber)%BASIS

        decompElement=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)
        !Get the dependent interpolation parameters
        dependentInterpolationParameters=>equations%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS( &
          & dependentVariable%VARIABLE_TYPE)%PTR
        dependentInterpolatedPoint=>equations%INTERPOLATION%DEPENDENT_INTERP_POINT( &
          & dependentVariable%VARIABLE_TYPE)%PTR
        !Get the geometric interpolation parameters
        geometricInterpolationParameters=>equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS( &
          & FIELD_U_VARIABLE_TYPE)%PTR
        geometricInterpolatedPoint=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
        CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,ERR,ERROR,*999)
        !Get access to geometric coordinates
        CALL FIELD_VARIABLE_GET(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,ERR,ERROR,*999)
        meshComponentNumber=geometricVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        !Get the geometric distributed vector
        CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & geometricParameters,ERR,ERROR,*999)
        fieldVariable=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
        !Get the density
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & 2,density,err,error,*999)

        ! Get the boundary element parameters
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & 1,beta,err,error,*999)
        boundaryNormal = 0.0_DP
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementNumber,5,boundaryNormal(1),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementNumber,6,boundaryNormal(2),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementNumber,7,boundaryNormal(3),err,error,*999)

        DO faceIdx=1,dependentBasis%NUMBER_OF_LOCAL_FACES
          !Get the face normal and quadrature information
          IF (ALLOCATED(decompElement%ELEMENT_FACES)) THEN
            faceNumber=decompElement%ELEMENT_FACES(faceIdx)
          ELSE
            CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
          END IF
          face=>decomposition%TOPOLOGY%FACES%FACES(faceNumber)
          !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
          !correspond to the other element.
          IF (.NOT.(face%BOUNDARY_FACE)) CYCLE

          SELECT CASE(dependentBasis%TYPE)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            normalcompIdx=ABS(face%XI_DIRECTION)
          CASE DEFAULT
            LOCAL_ERROR="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis%TYPE,"*",ERR,ERROR))// &
              & " is not yet implemented for Navier-Stokes."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

          faceBasis=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
          faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          DO gaussIdx=1,faceQuadratureScheme%NUMBER_OF_GAUSS
            !Get interpolated geometry
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
              & geometricInterpolatedPoint,err,error,*999)
            !Calculate point metrics 
            pointMetrics=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)

            ! TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
            !Get face normal projection
            DO compIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
              normalProjection(compIdx)=DOT_PRODUCT(pointMetrics%GU(normalcompIdx,:),pointMetrics%DX_DXI(compIdx,:))
              IF(face%XI_DIRECTION<0) THEN
                normalProjection(compIdx)=-normalProjection(compIdx)
              END IF
            END DO
            IF(L2NORM(normalProjection)>ZERO_TOLERANCE) THEN
               unitNormal=normalProjection/L2NORM(normalProjection)
            ELSE
               unitNormal=0.0_DP
            END IF

            ! Stabilisation term to correct for possible retrograde flow divergence.
            ! See: Moghadam et al 2011 “A comparison of outlet boundary treatments for prevention of backflow divergence..." and
            !      Ismail et al 2014 "A stable approach for coupling multidimensional cardiovascular and pulmonary networks..."
            ! Note: beta is a relative scaling factor 0 <= beta <= 1; default 1.0
            stabilisationTerm = 0.0_DP
            normalDifference=L2NORM(boundaryNormal-unitNormal)
            normalTolerance=0.1_DP
            IF(normalDifference < normalTolerance) THEN
              normalFlow = DOT_PRODUCT(velocity,normalProjection)
              !normalFlow = DOT_PRODUCT(velocity,boundaryNormal)
              IF(normalFlow < -ZERO_TOLERANCE) THEN
                DO compIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
                  stabilisationTerm(compIdx) = 0.5_DP*beta*density*velocity(compIdx)*(normalFlow - ABS(normalFlow))
                END DO
              ELSE
                stabilisationTerm = 0.0_DP
              END IF
            ELSE
              ! Not the correct boundary face - go to next face
              EXIT
            END IF

            ! Interpolate applied boundary pressure value
            boundaryPressure=0.0_DP
            !Get the pressure value interpolation parameters
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PRESSURE_VALUES_SET_TYPE,elementNumber,equations% &
             & INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
             & dependentInterpolatedPoint,ERR,ERROR,*999)
            boundaryPressure=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(4,NO_PART_DERIV)

            ! Interpolate current solution velocity/pressure field values
            pressure=0.0_DP
            velocity=0.0_DP
            dUDXi=0.0_DP
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,equations%INTERPOLATION% &
             & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)              
            !Get interpolated velocity and pressure 
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
              & dependentInterpolatedPoint,err,error,*999)
            velocity(1)=dependentInterpolatedPoint%values(1,NO_PART_DERIV) 
            velocity(2)=dependentInterpolatedPoint%values(2,NO_PART_DERIV) 
            velocity(3)=dependentInterpolatedPoint%values(3,NO_PART_DERIV) 
            dUDXi(1:3,1)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S1)
            dUDXi(1:3,2)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S2)
            dUDXi(1:3,3)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S3)
            pressure=dependentInterpolatedPoint%values(4,NO_PART_DERIV) 

            ! Keep this here for now: not using for Pressure BC but may want for traction BC
            ! ! Calculate viscous term
            ! dXiDX=0.0_DP
            ! dXiDX=pointMetrics%DXI_DX(:,:)
            ! CALL MATRIX_PRODUCT(dUDXi,dXiDX,gradU,err,error,*999)
            ! DO i=1,numberOfDimensions 
            !   SUM1 = 0.0_DP
            !   SUM2 = 0.0_DP
            !   DO j=1,numberOfDimensions
            !      SUM1 = normalProjection(j)*gradU(i,j)
            !      SUM2 = normalProjection(j)*gradU(j,i)
            !   END DO
            !   normalViscousTerm(i) = viscosity*(SUM1 + SUM2)
            ! END DO

            !Jacobian and Gauss weighting term
            jacobianGaussWeights=pointMetrics%JACOBIAN*faceQuadratureScheme%GAUSS_WEIGHTS(gaussIdx)

            !Loop over field components
            DO compIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
              !Work out the first index of the rhs vector for this element - (i.e. the number of previous)
              elementBaseDofIdx=dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(compIdx-1)
              DO faceNodeIdx=1,faceBasis%NUMBER_OF_NODES
                elementNodeIdx=dependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(faceNodeIdx,faceIdx)
                DO faceNodederivIdx=1,faceBasis%NUMBER_OF_DERIVATIVES(faceNodeIdx)
                  nodederivIdx=dependentBasis%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(faceNodederivIdx,faceNodeIdx,faceIdx)
                  elementParameterIdx=dependentBasis%ELEMENT_PARAMETER_INDEX(nodederivIdx,elementNodeIdx)
                  faceParameterIdx=faceBasis%ELEMENT_PARAMETER_INDEX(faceNodederivIdx,faceNodeIdx)
                  elementDof=elementBaseDofIdx+elementParameterIdx

                  rhsVector%ELEMENT_VECTOR%VECTOR(elementDof) = rhsVector%ELEMENT_VECTOR%VECTOR(elementDof) - &
                    &  (boundaryPressure*normalProjection(compIdx) - stabilisationTerm(compIdx))* &
                    &  faceQuadratureScheme%GAUSS_BASIS_FNS(faceParameterIdx,NO_PART_DERIV,gaussIdx)* &
                    &  jacobianGaussWeights

                END DO !nodederivIdx
              END DO !faceNodeIdx
            END DO !compIdx
          END DO !gaussIdx
        END DO !faceIdx

        !Restore the distributed geometric data used for the normal calculation
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & geometricParameters,ERR,ERROR,*999)
      ENDIF !decomposition%calculate_faces

    CASE DEFAULT
      ! Do nothing for other equation set subtypes
    END SELECT

    EXITS("NavierStokes_FiniteElementFaceIntegrate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementFaceIntegrate",err,error)
    RETURN 1
  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !> Calculate the fluid flux through 3D boundaries for use in problems with coupled solutions (e.g. multidomain)
  SUBROUTINE NavierStokes_CalculateBoundaryFlux(solver,err,error,*)

    !Argument variables

    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: decompElement
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(BASIS_TYPE), POINTER :: dependentBasis2
    TYPE(BASIS_TYPE), POINTER :: geometricFaceBasis
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: face
    TYPE(BASIS_TYPE), POINTER :: faceBasis
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: faceQuadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: rhsVector
    INTEGER(INTG) :: faceIdx, faceNumber,elementIdx,nodeNumber,versionNumber
    INTEGER(INTG) :: compIdx,gaussIdx
    INTEGER(INTG) :: elementBaseDofIdx, faceNodeIdx, elementNodeIdx
    INTEGER(INTG) :: faceNodederivIdx, meshComponentNumber, nodederivIdx, parameterIdx
    INTEGER(INTG) :: faceParameterIdx, elementDofIdx,normalcompIdx
    INTEGER(INTG) :: boundaryID
    INTEGER(INTG) :: MPI_IERROR,numberOfComputationalNodes
    REAL(DP) :: gaussWeight, normalProjection,elementNormal(3)
    REAL(DP) :: normalDifference,normalTolerance,faceFlux
    REAL(DP) :: courant,maxCourant,toleranceCourant
    REAL(DP) :: velocityGauss(3),faceNormal(3),unitNormal(3),boundaryValue,faceArea,faceVelocity
    REAL(DP) :: localBoundaryFlux(10),localBoundaryArea(10),globalBoundaryFlux(10),globalBoundaryArea(10)
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: correctFace

    ENTERS("NavierStokes_CalculateBoundaryFlux",err,error,*999)

    NULLIFY(decomposition)
    NULLIFY(geometricDecomposition)
    NULLIFY(geometricParameters)
    NULLIFY(decompElement)
    NULLIFY(dependentBasis)
    NULLIFY(dependentBasis2)
    NULLIFY(geometricFaceBasis)
    NULLIFY(geometricVariable)
    NULLIFY(equations)
    NULLIFY(equationsMatrices)
    NULLIFY(face)
    NULLIFY(faceBasis)
    NULLIFY(faceQuadratureScheme)
    NULLIFY(fieldVariable)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(geometricInterpolationParameters)
    NULLIFY(rhsVector)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(equationsEquationsSetField)
    NULLIFY(equationsSetField)

    !Some preliminary sanity checks
    IF (ASSOCIATED(SOLVER)) THEN
      solvers=>SOLVER%SOLVERS
      IF (ASSOCIATED(SOLVERS)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          SELECT CASE(controlLoop%PROBLEM%specification(3))
          CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
            &  PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE) 
            solverEquations=>solver%SOLVER_EQUATIONS
            IF (ASSOCIATED(solverEquations)) THEN
              solverMapping=>solverEquations%SOLVER_MAPPING
              IF (ASSOCIATED(solverMapping)) THEN
                equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
                IF (ASSOCIATED(equationsSet)) THEN
                  equations=>equationsSet%EQUATIONS
                    IF(ASSOCIATED(equations)) THEN
                      dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                      IF(.NOT.ASSOCIATED(dependentField)) THEN
                        CALL FlagError("Dependent field is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set equations is not associated.",err,error,*999)
                    END IF
                ELSE
                  CALL FlagError("Equations set is not associated.",err,error,*999)
                END IF
                equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
                IF (ASSOCIATED(equationsEquationsSetField)) THEN
                  equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
                  IF(.NOT.ASSOCIATED(equationsSetField)) THEN
                    CALL FlagError("Equations set field (EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations set field (EQUATIONS_EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver equations is not associated.",err,error,*999)
            END IF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
              & " is not valid for boundary flux calculation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    END IF

    localBoundaryArea=0.0_DP
    localBoundaryFlux=0.0_DP
    faceFlux=0.0_DP
    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)

      dependentVariable=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
      !Get the mesh decomposition and mapping
      decomposition=>dependentVariable%FIELD%DECOMPOSITION
      elementsMapping=>decomposition%DOMAIN(decomposition%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
      ! Get constant max Courant (CFL) number (default 1.0)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & 2,toleranceCourant,err,error,*999)

      ! Loop over elements to locate boundary elements
      maxCourant = 0.0_DP
      DO elementIdx=1,elementsMapping%TOTAL_NUMBER_OF_LOCAL
        meshComponentNumber=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        dependentBasis=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
          & ELEMENTS(elementIdx)%BASIS
        decompElement=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)
        ! Calculate element metrics (courant #, cell Reynolds number)
        CALL NavierStokes_CalculateElementMetrics(equationsSet,elementIdx,err,error,*999)            

        ! C F L  c o n d i t i o n   c h e c k
        ! ------------------------------------
        ! Get element metrics
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,3,courant,err,error,*999)
        IF(courant < -ZERO_TOLERANCE) THEN
          CALL FLAG_WARNING("Negative Courant (CFL) number.",ERR,ERROR,*999)
        END IF
        IF(courant > maxCourant) maxCourant = courant
        ! Check if element CFL number below specified tolerance
        IF(courant > toleranceCourant) THEN
          LOCAL_ERROR="Element "//TRIM(NUMBER_TO_VSTRING(decompElement%user_number, &
            & "*",ERR,ERROR))//" has violated the CFL condition "//TRIM(NUMBER_TO_VSTRING(courant, &
            & "*",ERR,ERROR))//" <= "//TRIM(NUMBER_TO_VSTRING(toleranceCourant,"*",ERR,ERROR))// &
            & ". Decrease timestep or increase CFL tolerance for the 3D Navier-Stokes problem."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END IF

        ! B o u n d a r y   n o r m a l   a n d   I D
        ! ----------------------------------------------
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,5,elementNormal(1),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,6,elementNormal(2),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,7,elementNormal(3),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,8,boundaryValue,err,error,*999)
        !Check if is a non-wall boundary element
        boundaryID=NINT(boundaryValue)
        IF(boundaryID>1) THEN
          faceArea=0.0_DP
          faceVelocity=0.0_DP
          !Get the dependent interpolation parameters
          dependentInterpolationParameters=>equations%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS( &
            & dependentVariable%VARIABLE_TYPE)%PTR
          dependentInterpolatedPoint=>equations%INTERPOLATION%DEPENDENT_INTERP_POINT( &
            & dependentVariable%VARIABLE_TYPE)%PTR
          ! Loop over faces to determine the boundary face contribution
          DO faceIdx=1,dependentBasis%NUMBER_OF_LOCAL_FACES
            !Get the face normal and quadrature information
            IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
              faceNumber=decompElement%ELEMENT_FACES(faceIdx)
            ELSE
              CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
            END IF
            face=>decomposition%TOPOLOGY%FACES%FACES(faceNumber)
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(face%BOUNDARY_FACE)) CYCLE

            SELECT CASE(dependentBasis%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              normalcompIdx=ABS(face%XI_DIRECTION)
            CASE(BASIS_SIMPLEX_TYPE)
              CALL FLAG_WARNING("Boundary flux calculation not yet set up for simplex element types.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis%TYPE,"*",ERR,ERROR))// &
                & " is not yet implemented for Navier-Stokes."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT

            CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceNumber, &
              & dependentInterpolationParameters,err,error,*999)
            faceBasis=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
            faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

            ! Loop over face gauss points
            DO gaussIdx=1,faceQuadratureScheme%NUMBER_OF_GAUSS
              !Use the geometric field to find the face normal and Jacobian for the face integral
              geometricInterpolationParameters=>equations%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS( &
                & FIELD_U_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementIdx, &
                & geometricInterpolationParameters,err,error,*999)
              geometricInterpolatedPoint=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                & geometricInterpolatedPoint,err,error,*999)
              pointMetrics=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)

              gaussWeight=faceQuadratureScheme%GAUSS_WEIGHTS(gaussIdx)
              !Get interpolated velocity
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                & dependentInterpolatedPoint,err,error,*999)
              !Interpolated values at gauss point
              velocityGauss=dependentInterpolatedPoint%values(1:3,NO_PART_DERIV)

              ! TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
              elementBaseDofIdx=0
              SELECT CASE(dependentBasis%TYPE)
              CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                correctFace=.TRUE.
                ! Make sure this is the boundary face that corresponds with boundaryID (could be a wall rather than inlet/outlet)
                DO compIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
                  normalProjection=DOT_PRODUCT(pointMetrics%GU(normalcompIdx,:),pointMetrics%DX_DXI(compIdx,:))
                  IF (face%XI_DIRECTION<0) THEN
                    normalProjection=-normalProjection
                  ENDIF
                  faceNormal(compIdx)=normalProjection
                ENDDO !compIdx
                unitNormal=faceNormal/L2NORM(faceNormal)
                normalDifference=L2NORM(elementNormal-unitNormal)
                normalTolerance=0.1_DP
                IF(normalDifference>normalTolerance) EXIT
              CASE(BASIS_SIMPLEX_TYPE)
                faceNormal=unitNormal
              CASE DEFAULT
                LOCAL_ERROR="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis%TYPE,"*",ERR,ERROR))// &
                  & " is not yet implemented for Navier-Stokes."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT

              ! Integrate face area and velocity
              DO compIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
                normalProjection=faceNormal(compIdx)
                IF(ABS(normalProjection)<ZERO_TOLERANCE) CYCLE
                !Work out the first index of the rhs vector for this element - 1
                elementBaseDofIdx=dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(compIdx-1)
                DO faceNodeIdx=1,faceBasis%NUMBER_OF_NODES
                  elementNodeIdx=dependentBasis%NODE_NUMBERS_IN_LOCAL_FACE(faceNodeIdx,faceIdx)
                  DO faceNodederivIdx=1,faceBasis%NUMBER_OF_DERIVATIVES(faceNodeIdx)
                    ! Integrate
                    nodederivIdx=1
                    parameterIdx=dependentBasis%ELEMENT_PARAMETER_INDEX(nodederivIdx,elementNodeIdx)
                    faceParameterIdx=faceBasis%ELEMENT_PARAMETER_INDEX(faceNodederivIdx,faceNodeIdx)
                    elementDofIdx=elementBaseDofIdx+parameterIdx
                    faceArea=faceArea + normalProjection*gaussWeight*pointMetrics%JACOBIAN* &
                     & faceQuadratureScheme%GAUSS_BASIS_FNS(faceParameterIdx,NO_PART_DERIV,gaussIdx)
                    faceVelocity=faceVelocity+velocityGauss(compIdx)*normalProjection*gaussWeight* &
                     & pointMetrics%JACOBIAN*faceQuadratureScheme%GAUSS_BASIS_FNS(faceParameterIdx,NO_PART_DERIV,gaussIdx)
                  END DO !nodederivIdx
                END DO !faceNodeIdx
              END DO !compIdx
            END DO !gaussIdx
          END DO !faceIdx
          localBoundaryFlux(boundaryID) = localBoundaryFlux(boundaryID) + faceVelocity
          localBoundaryArea(boundaryID) = localBoundaryArea(boundaryID) + faceArea
        END IF !boundaryIdentifier
      END DO !elementIdx                 

      ! Need to add boundary flux for any boundaries split accross computational nodes
      globalBoundaryFlux = 0.0_DP
      globalBoundaryArea = 0.0_DP
      numberOfComputationalNodes=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
      IF(numberOfComputationalNodes>1) THEN !use mpi
        CALL MPI_ALLREDUCE(localBoundaryFlux,globalBoundaryFlux,10,MPI_DOUBLE_PRECISION,MPI_SUM, &
          & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(localBoundaryArea,globalBoundaryArea,10,MPI_DOUBLE_PRECISION,MPI_SUM, &
          & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
      ELSE
        globalBoundaryFlux = localBoundaryFlux
        globalBoundaryArea = localBoundaryArea
      END IF
      
      ! Loop over elements again to allocate flux terms to boundary nodes
      DO elementIdx=1,elementsMapping%TOTAL_NUMBER_OF_LOCAL!elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,5,elementNormal(1),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,6,elementNormal(2),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,7,elementNormal(3),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,8,boundaryValue,err,error,*999)
        boundaryID=NINT(boundaryValue)
        IF(boundaryID>1) THEN
          meshComponentNumber=2
          decompElement=>decomposition%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)
          dependentBasis2=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
            & ELEMENTS(elementIdx)%BASIS
          DO faceIdx=1,dependentBasis2%NUMBER_OF_LOCAL_FACES
            !Get the face normal and quadrature information
            IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
              faceNumber=decompElement%ELEMENT_FACES(faceIdx)
            ELSE
              CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
            END IF
            face=>decomposition%TOPOLOGY%FACES%FACES(faceNumber)
            IF(.NOT.(face%BOUNDARY_FACE)) CYCLE
            faceBasis=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
            faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

            SELECT CASE(dependentBasis2%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              normalcompIdx=ABS(face%XI_DIRECTION)
              pointMetrics=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)
              DO compIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
                normalProjection=DOT_PRODUCT(pointMetrics%GU(normalcompIdx,:),pointMetrics%DX_DXI(compIdx,:))
                IF(face%XI_DIRECTION<0) THEN
                  normalProjection=-normalProjection
                END IF
                faceNormal(compIdx)=normalProjection
              END DO !compIdx
              unitNormal=faceNormal/L2NORM(faceNormal)
            CASE(BASIS_SIMPLEX_TYPE)
              !still have faceNormal/unitNormal
            CASE DEFAULT
              LOCAL_ERROR="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis2%TYPE,"*",ERR,ERROR))// &
                & " is not yet implemented for Navier-Stokes."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            normalDifference=L2NORM(elementNormal-unitNormal)
            normalTolerance=0.1_DP
            IF(normalDifference>normalTolerance) CYCLE

            ! Update local nodes with integrated boundary flow values
            DO faceNodeIdx=1,faceBasis%NUMBER_OF_NODES
              elementNodeIdx=dependentBasis2%NODE_NUMBERS_IN_LOCAL_FACE(faceNodeIdx,faceIdx)
              DO faceNodederivIdx=1,faceBasis%NUMBER_OF_DERIVATIVES(faceNodeIdx)
                nodeNumber=decomposition%DOMAIN(meshComponentNumber)%PTR% &
                 & TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)%ELEMENT_NODES(elementNodeIdx)
                versionNumber=1
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionNumber,faceNodederivIdx,nodeNumber,1,globalBoundaryFlux(boundaryID),err,error,*999) 
              END DO !nodederivIdx
            END DO !faceNodeIdx

          END DO !faceIdx
        END IF !boundaryIdentifier
      END DO !elementIdx                 

    CASE DEFAULT
      LOCAL_ERROR="Boundary flux calcluation for the third equations set specification of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",ERR,ERROR))// &
        & " is not yet implemented for Navier-Stokes."
      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)     
    END SELECT

    EXITS("NavierStokes_CalculateBoundaryFlux")
    RETURN
999 ERRORSEXITS("NavierStokes_CalculateBoundaryFlux",err,error)
    RETURN 1
  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !> Update the solution for the 1D solver with boundary conditions from a lumped parameter model defined by CellML.
  !> For more information please see chapter 11 of: L. Formaggia, A. Quarteroni, and A. Veneziani, Cardiovascular mathematics:
  !> modeling and simulation of the circulatory system. Milan; New York: Springer, 2009.
  SUBROUTINE NavierStokes_Couple1D0D(controlLoop,solver,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: iterativeLoop
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,materialsField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations  
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping 
    TYPE(SOLVER_TYPE), POINTER :: solver1D
    INTEGER(INTG) :: nodeNumber,nodeIdx,derivIdx,versionIdx,compIdx,numberOfLocalNodes,dependentDof
    INTEGER(INTG) :: solver1dNavierStokesNumber,solverNumber,MPI_IERROR,timestep,iteration
    INTEGER(INTG) :: boundaryNumber,numberOfBoundaries,numberOfComputationalNodes,boundaryConditionType
    REAL(DP) :: qPrevious,aPrevious,q1d,a1d,qError,aError,couplingTolerance,pCellML
    LOGICAL :: boundaryConverged(100),localConverged,MPI_LOGICAL
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_Couple1D0D",ERR,ERROR,*999)

    !Get solvers based on the problem type
    SELECT CASE(controlLoop%PROBLEM%specification(3))
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      solverNumber=solver%GLOBAL_NUMBER
      !In the Navier-Stokes/Characteristic subloop, the Navier-Stokes solver should be the second solver
      solver1dNavierStokesNumber=2
      versionIdx=1
      compIdx=1
      derivIdx=1
      IF(solverNumber==solver1dNavierStokesNumber) THEN
        solver1D=>controlLoop%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%PTR
        iterativeLoop=>controlLoop%WHILE_LOOP
        iteration=iterativeLoop%ITERATION_NUMBER
        timestep=controlLoop%PARENT_LOOP%TIME_LOOP%ITERATION_NUMBER
      ELSE 
        localError="The solver number of "//TRIM(NUMBER_TO_VSTRING(solverNumber,"*",err,error))// &
         & " does not correspond with the Navier-Stokes solver number for 1D-0D fluid coupling."
        CALL FlagError(localError,ERR,ERROR,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
        & " is not valid for 1D-0D Navier-Stokes fluid coupling."
      CALL FlagError(localError,ERR,ERROR,*999)
    END SELECT

    couplingTolerance=iterativeLoop%ABSOLUTE_TOLERANCE

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver1D)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          solverEquations=>solver1D%SOLVER_EQUATIONS
          IF(ASSOCIATED(solverEquations)) THEN
            solverMapping=>solverEquations%SOLVER_MAPPING
            IF(ASSOCIATED(solverMapping)) THEN
              equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
              IF(ASSOCIATED(equationsSet)) THEN
                materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
                dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver mapping is not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Solver equations is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control Loop is not associated.",err,error,*999)
    END IF

    !Get the number of local nodes
    domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
      & TOPOLOGY%NODES
    IF(ASSOCIATED(domainNodes)) THEN
      numberOfLocalNodes=domainNodes%NUMBER_OF_NODES
    ELSE
      CALL FlagError("Domain nodes are not associated.",err,error,*999)
    END IF

    boundaryNumber=0
    boundaryConverged=.TRUE.
    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,numberOfLocalNodes
      nodeNumber=domainNodes%NODES(nodeIdx)%local_number
      !Get the boundary condition type for the dependent field primitive variables (Q,A)
      boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
      NULLIFY(fieldVariable)
      CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
      dependentDof=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeNumber)% &
        & DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
        & err,error,*999)
      boundaryConditionType=boundaryConditionsVariable%CONDITION_TYPES(dependentDof)                       

      !!!-- F i n d   B o u n d a r y   N o d e s --!!!
      IF(boundaryConditionType>0) THEN
        boundaryNumber=boundaryNumber+1
        boundaryConverged(boundaryNumber)=.FALSE.  

        ! C u r r e n t   Q 1 D , A 1 D , p C e l l M L   V a l u e s
        ! ------------------------------------------------------------
        !Get Q1D
        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivIdx,nodeNumber,1,q1d,err,error,*999)         
        !Get A1D
        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivIdx,nodeNumber,2,a1d,err,error,*999)
        !Get pCellML
        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivIdx,nodeNumber,1,pCellML,err,error,*999)

        ! C h e c k  1 D / 0 D   C o n v e r g e n c e   f o r   t h i s   n o d e
        ! -------------------------------------------------------------------------
        IF (iteration==1 .AND. timestep==0) THEN
          !Create the previous iteration field values type on the dependent field if it does not exist
          NULLIFY(fieldVariable)
          CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE)%PTR)) THEN
            CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ERR,ERROR,*999)
          ENDIF
          NULLIFY(fieldVariable)
          CALL FIELD_VARIABLE_GET(dependentField,FIELD_U1_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE)%PTR)) THEN
            CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
              & ERR,ERROR,*999)
          ENDIF
        ELSE
          !Get previous Q1D and A1D
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivIdx,nodeNumber,1,qPrevious,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivIdx,nodeNumber,2,aPrevious,err,error,*999)
          !Check if the boundary interface values have converged
          qError=ABS(qPrevious-q1d)
          aError=ABS(aPrevious-a1d)
          IF(qError<couplingTolerance .AND. aError<couplingTolerance) THEN
            boundaryConverged(boundaryNumber)=.TRUE.
          ENDIF
        ENDIF

        !Store current Q and p Boundary values as previous iteration value
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
          & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivIdx,nodeNumber,1,q1d,err,error,*999)          
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
          & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivIdx,nodeNumber,2,a1d,err,error,*999)          
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U1_VARIABLE_TYPE, &
          & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivIdx,nodeNumber,1,pCellML,err,error,*999)          

      ENDIF !find boundary nodes
    ENDDO !loop over nodes 
    numberOfBoundaries=boundaryNumber

    IF(solverNumber==solver1dNavierStokesNumber) THEN
      ! ------------------------------------------------------------------
      ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
      ! ------------------------------------------------------------------
      !Check whether all boundaries on the local process have converged
      IF(numberOfBoundaries==0 .OR. ALL(boundaryConverged(1:numberOfBoundaries))) THEN
        localConverged=.TRUE.
      ELSE
        localConverged=.FALSE.
      ENDIF
      !Need to check that boundaries have converged globally (on all domains) if this is a parallel problem
      numberOfComputationalNodes=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
      IF(numberOfComputationalNodes>1) THEN !use mpi
        !Allocate array for mpi communication
        ALLOCATE(globalConverged(numberOfComputationalNodes),STAT=ERR) 
        IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
        CALL MPI_ALLGATHER(localConverged,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL, &
          & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
        IF(ALL(globalConverged)) THEN
          !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"1D/0D coupling converged; # iterations: ", &
          !  & iteration,err,error,*999)
          iterativeLoop%CONTINUE_LOOP=.FALSE.
        END IF
        DEALLOCATE(globalConverged)
      ELSE
        IF(localConverged) THEN
          !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"1D/0D coupling converged; # iterations: ", &
          !  & iteration,err,error,*999)
          iterativeLoop%CONTINUE_LOOP=.FALSE.
        END IF
      END IF

      !If the solution hasn't converged, need to revert field values to pre-solve state
      ! before continued iteration. This will counteract the field updates that occur
      ! in SOLVER_DYNAMIC_MEAN_PREDICTED_CALCULATE. Ignore for initialisation
      IF(timestep==0) THEN
        iterativeLoop%CONTINUE_LOOP=.FALSE.
      END IF
      IF(iterativeLoop%CONTINUE_LOOP .EQV. .TRUE. ) THEN
        CALL FIELD_PARAMETER_SETS_COPY(dependentField,equationsSet%EQUATIONS%EQUATIONS_MAPPING%DYNAMIC_MAPPING% &
          & DYNAMIC_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
        CALL FIELD_PARAMETER_SETS_COPY(dependentField,equationsSet%EQUATIONS%EQUATIONS_MAPPING%DYNAMIC_MAPPING% &
         & DYNAMIC_VARIABLE_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_RESIDUAL_SET_TYPE,1.0_DP,ERR,ERROR,*999)
      END IF
    END IF

    EXITS("NavierStokes_Couple1D0D")
    RETURN
999 ERRORSEXITS("NavierStokes_Couple1D0D",err,error)
    RETURN 1

  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !> Check convergence of 
  SUBROUTINE NavierStokes_CoupleCharacteristics(controlLoop,solver,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: iterativeLoop
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,independentField,materialsField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations  
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping 
    TYPE(SOLVER_TYPE), POINTER :: solver1DNavierStokes
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeNumber,nodeIdx,derivIdx,i,compIdx,timestep,numberOfLocalNodes
    INTEGER(INTG) :: solver1dNavierStokesNumber,solverNumber,MPI_IERROR,iteration,outputIteration
    INTEGER(INTG) :: branchNumber,numberOfBranches,numberOfComputationalNodes,numberOfVersions
    REAL(DP) :: wNavierStokes(4),wCharacteristic(4),wNext(4)
    REAL(DP) :: l2ErrorQ(100),qCharacteristic(4),qNavierStokes(4)
    REAL(DP) :: l2ErrorA(100),aCharacteristic(4),aNavierStokes(4)
    REAL(DP) :: startTime,stopTime,currentTime,timeIncrement,kp,k1,k2,k3,b1
    REAL(DP) :: RHO,normalWave,A0,E,H,beta,couplingTolerance
    LOGICAL :: branchConverged(100),localConverged,MPI_LOGICAL
    LOGICAL, ALLOCATABLE :: globalConverged(:)

    ENTERS("NavierStokes_CoupleCharacteristics",ERR,ERROR,*999)

    SELECT CASE(controlLoop%PROBLEM%specification(3))
    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
      solver1dNavierStokesNumber=2
      solver1DNavierStokes=>controlLoop%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%PTR
      CALL CONTROL_LOOP_TIMES_GET(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
        & timestep,outputIteration,ERR,ERROR,*999)
      iteration=controlLoop%WHILE_LOOP%ITERATION_NUMBER
      iterativeLoop=>controlLoop%WHILE_LOOP
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      solver1dNavierStokesNumber=2
      solver1DNavierStokes=>controlLoop%PARENT_LOOP%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%PTR
      iterativeLoop=>controlLoop%WHILE_LOOP
      iteration=iterativeLoop%ITERATION_NUMBER
      timestep=controlLoop%PARENT_LOOP%PARENT_LOOP%TIME_LOOP%ITERATION_NUMBER
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
        & " is not valid for 1D-0D Navier-Stokes fluid coupling."
      CALL FlagError(localError,ERR,ERROR,*999)
    END SELECT

    solverNumber=solver%GLOBAL_NUMBER
    couplingTolerance=iterativeLoop%ABSOLUTE_TOLERANCE

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver1DNavierStokes)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          solverEquations=>solver1DNavierStokes%SOLVER_EQUATIONS
          IF(ASSOCIATED(solverEquations)) THEN
            solverMapping=>solverEquations%SOLVER_MAPPING
            IF(ASSOCIATED(solverMapping)) THEN
              equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
              IF(ASSOCIATED(equationsSet)) THEN
                dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
                materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver mapping is not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Solver equations is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control Loop is not associated.",err,error,*999)
    END IF

    !Get the number of local nodes
    domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
      & TOPOLOGY%NODES
    IF(ASSOCIATED(domainNodes)) THEN
      numberOfLocalNodes=domainNodes%NUMBER_OF_NODES
    ELSE
      CALL FlagError("Domain nodes are not associated.",err,error,*999)
    END IF

    compIdx=1
    derivIdx=1
    branchNumber=0
    branchConverged=.TRUE.
    l2ErrorQ=0.0_DP
    l2ErrorA=0.0_DP

    !Get material constants
    CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & 2,RHO,err,error,*999)

    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,numberOfLocalNodes
      nodeNumber=domainNodes%NODES(nodeIdx)%local_number
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)%numberOfVersions             

      !Find branch nodes
      IF(numberOfVersions>1) THEN
        branchNumber=branchNumber+1
        branchConverged(branchNumber)=.FALSE.
        !Loop over versions
        DO i=1,numberOfVersions
          CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,i,derivIdx,nodeNumber,compIdx,normalWave,err,error,*999)
          !Get node-based material parameters
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,A0,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,E,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,3,H,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,4,kp,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,5,k1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,6,k2,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,7,k3,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,8,b1,err,error,*999)
          beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

          !Get current Q,A values based on N-S solve
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,qNavierStokes(i),err,error,*999)         
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,aNavierStokes(i),err,error,*999)
          !Get characteristic (flux conserving) Q,A values
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,qCharacteristic(i),err,error,*999)         
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,aCharacteristic(i),err,error,*999)

          !Calculate the characteristic based on the values converged upon by the N-S solver at this iteration.
          wNavierStokes(i)=qNavierStokes(i)/aNavierStokes(i)+normalWave*SQRT(4.0_DP*beta/RHO/b1)* &
            & ((aNavierStokes(i)/A0)**(b1/2.0_DP)-1.0_DP)
          !Calculate the characteristic based on the upwind values
          wCharacteristic(i)=qCharacteristic(i)/aCharacteristic(i)+normalWave*SQRT(4.0_DP*beta/RHO/b1)* &        
            & ((aCharacteristic(i)/A0)**(b1/2.0_DP)-1.0_DP)
        ENDDO

        !Evaluate error between current and previous Q,A values
        l2ErrorQ(branchNumber)=L2NORM(qNavierStokes-qCharacteristic)
        l2ErrorA(branchNumber)=L2NORM(aNavierStokes-aCharacteristic)

        !Check if the branch values have converged
        IF((ABS(l2ErrorQ(branchNumber))<couplingTolerance) .AND. (ABS(l2ErrorA(branchNumber))<couplingTolerance)) THEN
          branchConverged(branchNumber)=.TRUE.
        ENDIF

        wNext=0.5_DP*(wNavierStokes+wCharacteristic)
        !If N-S/C w values did not converge re-solve with new w.
        DO i=1,numberOfVersions      
          !Update W value
          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,compIdx,wNext(i),err,error,*999)
        ENDDO
      ENDIF !find branch nodes
    ENDDO !loop over nodes 
    numberOfBranches=branchNumber

    ! ------------------------------------------------------------------
    ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
    ! ------------------------------------------------------------------
    !Check whether all branches on the local process have converged
    IF(numberOfBranches==0 .OR. ALL(branchConverged(1:numberOfBranches))) THEN
      localConverged=.TRUE.
    ELSE
      localConverged=.FALSE.
    ENDIF
    !Need to check that boundaries have converged globally (on all domains) if this is a parallel problem
    numberOfComputationalNodes=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
    IF(numberOfComputationalNodes>1) THEN !use mpi
      !Allocate array for mpi communication
      ALLOCATE(globalConverged(numberOfComputationalNodes),STAT=ERR) 
      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
      CALL MPI_ALLGATHER(localConverged,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL, &
        & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
      IF(ALL(globalConverged)) THEN
        !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Navier-Stokes/Characteristic converged; # iterations: ", &
        !  & iteration,err,error,*999)
        controlLoop%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
      ENDIF
      DEALLOCATE(globalConverged)
    ELSE
      IF(localConverged) THEN
        !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Navier-Stokes/Characteristic converged; # iterations: ", &
        !  & iteration,err,error,*999)
        controlLoop%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
      END IF
    END IF

    EXITS("NavierStokes_CoupleCharacteristics")
    RETURN
999 ERRORSEXITS("NavierStokes_CoupleCharacteristics",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_CoupleCharacteristics

  !
  !================================================================================================================================
  !

  !>Calculated the rate of deformation (shear rate) for a navier-stokes finite element equations set.
  SUBROUTINE NavierStokes_ShearRateCalculate(equationsSet,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable,fieldVariable
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint,geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(FIELD_TYPE), POINTER :: dependentField,materialsField
    INTEGER(INTG) :: elementIdx,decompositionLocalElementNumber,gaussIdx
    INTEGER(INTG) :: meshComponentNumber,numberOfDimensions,i,j,userElementNumber
    INTEGER(INTG) :: localElementNumber,startElement,stopElement
    REAL(DP) :: gaussWeight,shearRate,secondInvariant,strainRate,shearRateDefault
    REAL(DP) :: dUdXi(3,3),dXidX(3,3),dUdX(3,3),dUdXTrans(3,3),rateOfDeformation(3,3),velocityGauss(3)
    LOGICAL :: ghostElement,elementExists,defaultUpdate
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_ShearRateCalculate",err,error,*999)

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...Calculating shear rate...",ERR,ERROR,*999)

    NULLIFY(decomposition)
    NULLIFY(dependentBasis)
    NULLIFY(equations)
    NULLIFY(quadratureScheme)
    NULLIFY(fieldVariable)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(dependentField)
    NULLIFY(materialsField)

    !Some preliminary sanity checks
    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
        IF(ASSOCIATED(equations)) THEN
          dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
          IF(.NOT.ASSOCIATED(dependentField)) THEN
            CALL FlagError("Dependent field is not associated.",err,error,*999)
          END IF
          materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
          IF(.NOT.ASSOCIATED(materialsField)) THEN
            CALL FlagError("Materials field is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set equations is not associated.",err,error,*999)
        END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
      dependentVariable=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
      meshComponentNumber=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
      !Get the mesh decomposition and mapping
      decomposition=>dependentVariable%FIELD%DECOMPOSITION
      elementsMapping=>decomposition%DOMAIN(decomposition%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
      fieldVariable=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
      numberOfDimensions=fieldVariable%NUMBER_OF_COMPONENTS - 1
      dependentInterpolatedPoint=>equations%INTERPOLATION%DEPENDENT_INTERP_POINT(dependentVariable%VARIABLE_TYPE)%PTR
      geometricInterpolatedPoint=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
      defaultUpdate=.FALSE.

      ! Loop over internal and boundary elements, skipping ghosts
      startElement = elementsMapping%INTERNAL_START
      stopElement = elementsMapping%BOUNDARY_FINISH
      ! Loop over internal and boundary elements 
      DO elementIdx=startElement,stopElement
        localElementNumber=elementsMapping%DOMAIN_LIST(elementIdx)
        userElementNumber = elementsMapping%LOCAL_TO_GLOBAL_MAP(localElementNumber)
        !Check computational node for elementIdx
        elementExists=.FALSE.
        ghostElement=.TRUE.
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decomposition%TOPOLOGY, &
          & userElementNumber,elementExists,decompositionLocalElementNumber,ghostElement,ERR,ERROR,*999)              
        IF(ghostElement) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Ghost: ",userElementNumber,err,error,*999)            
        END IF

        IF(elementExists) THEN
          dependentBasis=>decomposition%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
          quadratureScheme=>dependentBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber,equations%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber,equations%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)

          ! Loop over gauss points
          DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
            !Get interpolated velocity
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & dependentInterpolatedPoint,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & geometricInterpolatedPoint,err,error,*999)
            pointMetrics=>equations%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)
            gaussWeight=quadratureScheme%GAUSS_WEIGHTS(gaussIdx)
            !Interpolated values at gauss point
            dXidX=0.0_DP
            dUdXi=0.0_DP
            velocityGauss=dependentInterpolatedPoint%values(1:3,NO_PART_DERIV)

            dUdXi(1:3,1)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S1)
            dUdXi(1:3,2)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S2)
            IF(numberOfDimensions == 3) THEN
              dUdXi(1:3,3)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S3)
            ELSE
              dUdXi(1:3,3)=0.0_DP
            END IF
            dXidX=pointMetrics%DXI_DX(:,:)
            dUdX=0.0_DP
            dUdXTrans=0.0_DP
            strainRate=0.0_DP

            CALL MATRIX_PRODUCT(dUdXi,dXidX,dUdX,ERR,ERROR,*999) !dU/dX = dU/dxi * dxi/dX (deformation gradient tensor)
            CALL MATRIX_TRANSPOSE(dUdX,dUdXTrans,ERR,ERROR,*999)
            DO i=1,3
              DO j=1,3
                strainRate = strainRate + (dUdX(i,j)*dUdXTrans(i,j))
                rateOfDeformation(i,j) = (dUdX(i,j) + dUdXTrans(i,j))/2.0_DP
              END DO
            END DO            
            secondInvariant= - rateOfDeformation(1,2)**2.0_DP - &
               & rateOfDeformation(2,3)**2.0_DP - rateOfDeformation(1,3)**2.0_DP

            IF(secondInvariant > -1.0E-30_DP) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE, &
               & "WARNING: positive second invariant of rate of deformation tensor: ",secondInvariant,err,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Element number: ",userElementNumber,err,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"   Gauss point number: ",gaussIdx,err,error,*999)
              defaultUpdate=.TRUE.
              EXIT
            ELSE
              shearRate=SQRT(-4.0_DP*secondInvariant)
            END IF
            CALL Field_ParameterSetUpdateLocalGaussPoint(materialsField,FIELD_V_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,gaussIdx,localElementNumber,2,shearRate,ERR,ERROR,*999)

          END DO !gaussIdx
        END IF ! check for ghost element
        IF(defaultUpdate .EQV. .TRUE.) THEN
          EXIT
        END IF
      END DO !elementIdx

      IF(defaultUpdate .EQV. .TRUE.) THEN
        shearRateDefault=1.0E-10_DP
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Setting default shear field values...", &
         & shearRateDefault,err,error,*999)
        CALL FIELD_COMPONENT_VALUES_INITIALISE(materialsField,FIELD_V_VARIABLE_TYPE, &
         & FIELD_VALUES_SET_TYPE,1,shearRateDefault,err,error,*999)
      END IF

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",ERR,ERROR))// &
        & " is not valid for shear rate calculation in a Navier-Stokes equation type of a classical field equations set class."
      CALL FlagError(localError,ERR,ERROR,*999)
    END SELECT

    EXITS("NavierStokes_ShearRateCalculate")
    RETURN
999 ERRORSEXITS("NavierStokes_ShearRateCalculate",ERR,ERROR)
    RETURN 1

  END SUBROUTINE NavierStokes_ShearRateCalculate

  !
  !================================================================================================================================
  !

  !>Pre-residual evaluation a navier-stokes finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementPreResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NavierStokes_FiniteElementPreResidualEvaluate",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        ! Shear rate should either be calculated here to update at each minor iteration
        ! or during post solve so it is updated once per timestep
        !CALL NavierStokes_ShearRateCalculate(equationsSet,err,error,*999)
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        !Do nothing
      CASE DEFAULT
        LOCAL_ERROR="The third equations set specification of "// &
          & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid mechanics equations set."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_FiniteElementPreResidualEvaluate")
    RETURN
999 ERRORS("NavierStokes_FiniteElementPreResidualEvaluate",err,error)
    EXITS("NavierStokes_FiniteElementPreResidualEvaluate")
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementPreResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE NavierStokes_ControlLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_TYPE), POINTER :: navierStokesSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_ControlLoopPostLoop",err,error,*999)

    NULLIFY(dependentField)
    NULLIFY(fieldVariable)

    IF(ASSOCIATED(controlLoop)) THEN
      SELECT CASE(controlLoop%PROBLEM%specification(3))
      CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, & 
         & PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        !Do nothing
      CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          !Do nothing
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !Global time loop - export data
          navierStokesSolver=>controlLoop%SUB_LOOPS(1)%PTR%SOLVERS%SOLVERS(2)%PTR
          CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(navierStokesSolver,err,error,*999)
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          navierStokesSolver=>controlLoop%SOLVERS%SOLVERS(2)%PTR
          CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,ERR,ERROR,*999)
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NUMBER_TO_VSTRING(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          !CellML simple loop - do nothing 
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !Global time loop - export data
          navierStokesSolver=>controlLoop%SUB_LOOPS(1)%PTR%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(2)%PTR
          CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(navierStokesSolver,err,error,*999)
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          !Couple 1D/0D loop
          IF(controlLoop%CONTROL_LOOP_LEVEL==2) THEN
            navierStokesSolver=>controlLoop%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(2)%PTR
            !Update 1D/0D coupling parameters and check convergence
            CALL NavierStokes_Couple1D0D(controlLoop,navierStokesSolver,err,error,*999)
          !Couple Navier-Stokes/Characteristics loop
          ELSE IF(controlLoop%CONTROL_LOOP_LEVEL==3) THEN
            navierStokesSolver=>controlLoop%SOLVERS%SOLVERS(2)%PTR
            CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,ERR,ERROR,*999)
          ELSE
            localError="The while loop level of "//TRIM(NUMBER_TO_VSTRING(controlLoop%CONTROL_LOOP_LEVEL,"*",err,error))// &
              & " is invalid for a Coupled 1D0D Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          END IF
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NUMBER_TO_VSTRING(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_ControlLoopPostLoop")
    RETURN
999 ERRORSEXITS("NavierStokes_ControlLoopPostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE NavierStokes_ControlLoopPostLoop

  !
  !================================================================================================================================
  !

  !>Updates boundary conditions for multiscale fluid problems
  SUBROUTINE NavierStokes_UpdateMultiscaleBoundary(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,parentLoop,streeLoop
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: dependentDomain
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,streeEquationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations,streeEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,materialsField,streeMaterialsField,independentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations,streeSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,streeSolverMapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: streeSolver
    TYPE(VARYING_STRING) :: localError
    REAL(DP) :: A0,H,E,RHO,beta,Pext,lengthScale,timeScale,massScale,currentTime,timeIncrement
    REAL(DP) :: pCellml,qCellml,ABoundary,W1,W2,ACellML,normalWave,kp,k1,k2,k3,b1
    REAL(DP), POINTER :: Impedance(:),Flow(:)
    INTEGER(INTG) :: nodeIdx,versionIdx,derivIdx,compIdx,numberOfLocalNodes
    INTEGER(INTG) :: dependentDof,boundaryConditionType,k,nodeNumber
    LOGICAL :: boundaryNode

    ENTERS("NavierStokes_UpdateMultiscaleBoundary",err,error,*999)

    NULLIFY(dependentDomain)
    NULLIFY(equationsSet)
    NULLIFY(equations)
    NULLIFY(geometricField)
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(materialsField)
    NULLIFY(fieldVariable)
    NULLIFY(solverEquations)
    NULLIFY(solverMapping)
    normalWave=0.0_DP
    versionIdx=1
    derivIdx=1
    compIdx=1

    !Preliminary checks; get field and domain pointers
    IF(ASSOCIATED(solver)) THEN
      solvers=>solver%SOLVERS
      IF(ASSOCIATED(solvers)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        parentLoop=>controlLoop%PARENT_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop, &
          & currentTime,timeIncrement,ERR,ERROR,*999)
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          SELECT CASE(controlLoop%PROBLEM%specification(3))
          CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            solverEquations=>solver%SOLVER_EQUATIONS
            IF(ASSOCIATED(solverEquations)) THEN
              solverMapping=>solverEquations%SOLVER_MAPPING
              IF(ASSOCIATED(solverMapping)) THEN
                equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(equationsSet)) THEN
                  equations=>equationsSet%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
                    independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
                    dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(dependentField)) THEN
                      dependentDomain=>dependentField%DECOMPOSITION%DOMAIN(dependentField% &
                        & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
                      IF(.NOT.ASSOCIATED(dependentDomain)) THEN
                        CALL FlagError("Dependent domain is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Geometric field is not associated.",err,error,*999)
                    END IF
                    materialsField=>equations%INTERPOLATION%MATERIALS_FIELD
                    IF(.NOT.ASSOCIATED(materialsField)) THEN
                      CALL FlagError("Materials field is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set equations is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations set is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver equations is not associated.",err,error,*999)
            END IF
          CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            streeLoop=>parentLoop%SUB_LOOPS(1)%PTR
            streeSolver=>streeLoop%SOLVERS%SOLVERS(1)%PTR
            solverEquations=>solver%SOLVER_EQUATIONS
            streeSolverEquations=>streeSolver%SOLVER_EQUATIONS
            IF(ASSOCIATED(solverEquations)) THEN
              solverMapping=>solverEquations%SOLVER_MAPPING
              streeSolverMapping=>streeSolverEquations%SOLVER_MAPPING
              IF(ASSOCIATED(solverMapping)) THEN
                equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
                streeEquationsSet=>streeSolverMapping%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(equationsSet)) THEN
                  equations=>equationsSet%EQUATIONS
                  streeEquations=>streeEquationsSet%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
                    independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
                    dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                    streeMaterialsField=>streeEquationsSet%MATERIALS%MATERIALS_FIELD
                    IF(ASSOCIATED(dependentField)) THEN
                      dependentDomain=>dependentField%DECOMPOSITION%DOMAIN(dependentField% &
                        & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
                      IF(.NOT.ASSOCIATED(dependentDomain)) THEN
                        CALL FlagError("Dependent domain is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Geometric field is not associated.",err,error,*999)
                    END IF
                    materialsField=>equations%INTERPOLATION%MATERIALS_FIELD
                    IF(.NOT.ASSOCIATED(materialsField)) THEN
                      CALL FlagError("Materials field is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set equations is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations set is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver equations is not associated.",err,error,*999)
            END IF
          CASE DEFAULT
           localError="The third problem specification of "// &
             & TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
             & " is not valid for boundary flux calculation."
           CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    END IF
              
    SELECT CASE(equationsSet%specification(3))
    !!!-- 1 D    E q u a t i o n s   S e t --!!! 
    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)

      independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
      !Get the number of local nodes
      domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
        & TOPOLOGY%NODES
      IF(ASSOCIATED(domainNodes)) THEN
        numberOfLocalNodes=domainNodes%NUMBER_OF_NODES
      ELSE
        CALL FlagError("Domain nodes are not associated.",err,error,*999)
      END IF

      !Get constant material parameters
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,2,RHO,err,error,*999)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,4,Pext,err,error,*999)
      !Get materials scale factors
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,5,lengthScale,err,error,*999)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,6,timeScale,err,error,*999)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,7,massScale,err,error,*999)

      !!!--  L o o p   o v e r   l o c a l    n o d e s  --!!!
      DO nodeIdx=1,numberOfLocalNodes
        nodeNumber=domainNodes%NODES(nodeIdx)%local_number
        !Check for the boundary node
        boundaryNode=domainNodes%NODES(nodeNumber)%BOUNDARY_NODE

        !!!-- F i n d   b o u n d a r y    n o d e s --!!!
        IF(boundaryNode) THEN
          !Get material parameters
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,1,A0,err,error,*999) 
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,2,E,err,error,*999) 
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,3,H,err,error,*999) 
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,4,kp,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,5,k1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,6,k2,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,7,k3,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
            & derivIdx,nodeNumber,8,b1,err,error,*999)
          beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

          !Get normal wave direction
          CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & versionIdx,derivIdx,nodeNumber,compIdx,normalWave,err,error,*999)

          !Get the boundary condition type for the dependent field primitive variables (Q,A)
          boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
          NULLIFY(fieldVariable)
          CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          dependentDof=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeNumber)% &
            & DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
            & err,error,*999)
          boundaryConditionType=boundaryConditionsVariable%CONDITION_TYPES(dependentDof)                       
          SELECT CASE(boundaryConditionType)
          ! ----------------------------------------------------
          ! N o n - r e f l e c t i n g   B o u n d a r y
          ! ----------------------------------------------------
          CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING)
            ! Outlet - set W2 to 0, get W1 from the extrapolated value
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivIdx,nodeNumber,1,W1,err,error,*999)                
            W2=0.0_DP
            !Calculate new area value based on W1, W2 and update dof
            ABoundary=A0*((ABS(W1-W2)/2.0_DP)*SQRT(RHO*b1/beta/4.0_DP)+1.0_DP)**(2.0_DP/b1)
            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,versionIdx,derivIdx,nodeNumber,2,ABoundary,err,error,*999)
          ! ------------------------------------------------------------
          ! C o u p l e d   C e l l M L  ( 0 D )   B o u n d a r y
          ! ------------------------------------------------------------
          CASE(BOUNDARY_CONDITION_FIXED_CELLML)
            !Get qCellML used in pCellML calculation
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivIdx,nodeNumber,1,QCellML,err,error,*999)
            !Get pCellML if this is a coupled problem
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivIdx,nodeNumber,2,pCellml,err,error,*999)
            !Convert pCellML from SI base units specified in CellML file to scaled units (e.g., kg/(m.s^2) --> g/(mm.ms^2))
            pCellml=0.1_DP*pCellml*massScale/(lengthScale*(timeScale**2.0_DP))
            !Convert pCellML --> A0D
            ACellML=A0*((pCellml-Pext)/beta+1.0_DP)**(1.0_DP/b1)
            !  O u t l e t
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivIdx,nodeNumber,1,W1,err,error,*999)                    
            !Calculate W2 from 0D domain
            W2=QCellml/ACellml-normalWave*SQRT(4.0_DP*beta/RHO/b1)*((ACellml/A0)**(b1/2.0_DP)-1.0_DP)
            !Calculate new area value based on W1,W2 and update dof
            ABoundary=A0*((ABS(W1-W2)/2.0_DP)*SQRT(RHO*b1/beta/4.0_DP)+1.0_DP)**(2.0_DP/b1)
            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,versionIdx,derivIdx,nodeNumber,2,ABoundary,err,error,*999)
          ! ------------------------------------------------------------
          ! S t r u c t u r e d   T r e e   B o u n d a r y
          ! ------------------------------------------------------------
          CASE(BOUNDARY_CONDITION_FIXED_STREE)
            NULLIFY(Impedance)
            NULLIFY(Flow)
            !Get qCellML used in pCellML calculation
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivIdx,nodeNumber,1,QCellML,err,error,*999)
            !Get flow function
            CALL FIELD_PARAMETER_SET_DATA_GET(streeMaterialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & Impedance,err,error,*999)
            !Get impedance function
            CALL FIELD_PARAMETER_SET_DATA_GET(streeMaterialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & Flow,err,error,*999)
            pCellml=0.0_DP
            DO k=1,size(Flow)
              pCellml=pCellml+Flow(k)*Impedance(k)*timeIncrement
            ENDDO
            !Convert pCellML --> A0D 
            ACellML=A0*((pCellml-Pext)/beta+1.0_DP)**(1.0_DP/b1)
            !  O u t l e t
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivIdx,nodeNumber,1,W1,err,error,*999)                    
            !Calculate W2 from 0D domain
            W2=QCellml/ACellml-normalWave*SQRT(4.0_DP*beta/RHO/b1)*((ACellml/A0)**(b1/2.0_DP)-1.0_DP)
            !Calculate new area value based on W1, W2 and update dof
            ABoundary=A0*((ABS(W1-W2)/2.0_DP)*SQRT(RHO*b1/beta/4.0_DP)+1.0_DP)**(2.0_DP/b1)
            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,versionIdx,derivIdx,nodeNumber,2,ABoundary,err,error,*999)
          ! ------------------------------------------------------------
          CASE(BOUNDARY_CONDITION_FREE, &
             & BOUNDARY_CONDITION_FIXED, &
             & BOUNDARY_CONDITION_FIXED_INLET, &
             & BOUNDARY_CONDITION_FIXED_OUTLET, &
             & BOUNDARY_CONDITION_FIXED_FITTED)
            !Do nothing

          CASE DEFAULT
            localError="The boundary conditions type "//TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
              & " is not valid for a coupled characteristic problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        END IF !boundary node
      END DO !Loop over nodes

    !!!-- 3 D    E q u a t i o n s   S e t --!!! 
    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
      !Do nothing

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
       & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_UpdateMultiscaleBoundary")
    RETURN
999 ERRORSEXITS("NavierStokes_UpdateMultiscaleBoundary",err,error)
    RETURN 1
    
  END SUBROUTINE NavierStokes_UpdateMultiscaleBoundary

  !
  !================================================================================================================================
  !

END MODULE NAVIER_STOKES_EQUATIONS_ROUTINES
