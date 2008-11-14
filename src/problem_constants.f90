!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all problem wide constants.
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

!> This module handles all problem wide constants.
MODULE PROBLEM_CONSTANTS
  
  USE KINDS

  IMPLICIT NONE

  !Module parameters

  !Problem Classes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_CLASS=0  
  INTEGER(INTG), PARAMETER :: PROBLEM_ELASTICITY_CLASS=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FLUID_MECHANICS_CLASS=2
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROMAGNETICS_CLASS=3
  INTEGER(INTG), PARAMETER :: PROBLEM_CLASSICAL_FIELD_CLASS=4  
  INTEGER(INTG), PARAMETER :: PROBLEM_MODAL_CLASS=5
  INTEGER(INTG), PARAMETER :: PROBLEM_FITTING_CLASS=6
  INTEGER(INTG), PARAMETER :: PROBLEM_OPTIMISATION_CLASS=7

  !Problem types
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_TYPE=0
  !Elasticity class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTICITY_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_FINITE_ELASTICITY_TYPE=2
  !Fluid mechanics class
  INTEGER(INTG), PARAMETER :: PROBLEM_STOKES_FLUID_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NAVIER_STOKES_FLUID_TYPE=2
  !Electromagnetics class
  INTEGER(INTG), PARAMETER :: PROBLEM_ELECTROSTATIC_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_MAGNETOSTATIC_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_MAXWELLS_EQUATIONS_TYPE=3
  !Classical field class
  INTEGER(INTG), PARAMETER :: PROBLEM_LAPLACE_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_POISSON_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: PROBLEM_HELMHOLTZ_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: PROBLEM_WAVE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: PROBLEM_DIFFUSION_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE=6
  INTEGER(INTG), PARAMETER :: PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE=7
  INTEGER(INTG), PARAMETER :: PROBLEM_BIHARMONIC_EQUATION_TYPE=8
  !Modal class
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_ELASTIC_MODAL_TYPE=1
  
  !Problem subtypes
  INTEGER(INTG), PARAMETER :: PROBLEM_NO_SUBTYPE=0
  !Elasticity class
  !Fluid mechanics class
  !Electromagnetics class
  !Classical field class
  !  Laplace equation
!!TODO: We don't really have two problem types here? Maybe a different type with nonlinear boundary conditions???
  INTEGER(INTG), PARAMETER :: PROBLEM_STANDARD_LAPLACE_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_GENERALISED_LAPLACE_SUBTYPE=2
  !  Poisson equation
  INTEGER(INTG), PARAMETER :: PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE=2
  !  Helmholtz equation
  !  Wave equation
  !  Diffusion equation
  !Modal class

    
  !> \addtogroup PROBLEM_CONSTANTS_SetupTypes PROBLEM_CONSTANTS::SetupTypes
  !> \brief Setup type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_INITIAL_TYPE=1 !<Initial setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_CONTROL_TYPE=2 !<Solver setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_ROU
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_SOLUTION_TYPE=3 !<Solution parameters setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_SOLVER_TYPE=4 !<Solver setup for a problem. \see PROBLEM_CONSTANTS_SetupTypes,PROBLEM_CONSTANTS
  !>@}
  
  !> \addtogroup PROBLEM_CONSTANTS_SetupActionTypes PROBLEM_CONSTANTS::SetupActionTypes
  !> \brief Setup action type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_START_ACTION=1 !<Start setup action. \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_FINISH_ACTION=2 !<Finish setup action. \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
  INTEGER(INTG), PARAMETER :: PROBLEM_SETUP_DO_ACTION=3 !<Do setup action. \see PROBLEM_CONSTANTS_SetupActionTypes,CONSTANTS_ROUTINES
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_SolutionOutputTypes PROBLEM_CONSTANTS::SolutionOutputTypes
  !> \brief Solution output types
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_NO_OUTPUT=0 !<No output \see PROBLEM_CONSTANTS_SolutionOutputTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_TIMING_OUTPUT=1 !<Timing information output \see PROBLEM_CONSTANTS_SolutionOutputTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_MATRIX_OUTPUT=2 !<All below and solution matrices output \see PROBLEM_CONSTANTS_SolutionOutputTypes,PROBLEM_CONSTANTS
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_SolutionLinearityTypes PROBLEM_CONSTANTS::SolutionLinearityTypes
  !> \brief The solution linearity type parameters
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_LINEAR=1 !<The problem solution is linear. \see PROBLEM_CONSTANTS_SolutionLinearityTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLUTION_NONLINEAR=2 !<The problem solution is nonlinear. \see PROBLEM_CONSTANTS_SolutionLinearityTypes,PROBLEM_CONSTANTS
  !>@}


  !> \addtogroup PROBLEM_CONSTANTS_SolverOutputTypes PROBLEM_CONSTANTS::SolverOutputTypes
  !> \brief Solver output types
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_NO_OUTPUT=0 !<No output \see PROBLEM_CONSTANTS_SolverOutputTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_TIMING_OUTPUT=1 !<Timing information output \see PROBLEM_CONSTANTS_SolverOutputTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_SOLVER_OUTPUT=2 !<All below and solver output \see PROBLEM_CONSTANTS_SolverOutputTypes,PROBLEM_CONSTANTS
  !>@}

  !> \addtogroup PROBLEM_CONSTANTS_SolverSparsityTypes PROBLEM_CONSTANTS::SolverSparsityTypes
  !> \brief Solution output types
  !> \see PROBLEM_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices \see PROBLEM_CONSTANTS_SolverSparsityTypes,PROBLEM_CONSTANTS
  INTEGER(INTG), PARAMETER :: PROBLEM_SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices \see PROBLEM_CONSTANTS_SolverSparsityTypes,PROBLEM_CONSTANTS
 !>@}
  
END MODULE PROBLEM_CONSTANTS
