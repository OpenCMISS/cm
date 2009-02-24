!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all constants shared across equations set routines.
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

!> This module defines all constants shared across equations set routines.
MODULE EQUATIONS_SET_CONSTANTS

  USE KINDS

  IMPLICIT NONE

  !Problem Classes
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NO_CLASS=0
  
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_ELASTICITY_CLASS=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FLUID_MECHANICS_CLASS=2
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_ELECTROMAGNETICS_CLASS=3
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_CLASSICAL_FIELD_CLASS=4
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_BIOELECTRICS_CLASS=5
  
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_MODAL_CLASS=6
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FITTING_CLASS=7
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_OPTIMISATION_CLASS=8

  !Problem types
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NO_TYPE=0
  !Elasticity class
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LINEAR_ELASTICITY_TYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FINITE_ELASTICITY_TYPE=2
  !Fluid mechanics class
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_STOKES_FLUID_TYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NAVIER_STOKES_FLUID_TYPE=2
  !Electromagnetics class
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_ELECTROSTATIC_TYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_MAGNETOSTATIC_TYPE=2
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_MAXWELLS_EQUATIONS_TYPE=3
  !Classical field class
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LAPLACE_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_POISSON_EQUATION_TYPE=2
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE=3
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_WAVE_EQUATION_TYPE=4
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_DIFFUSION_EQUATION_TYPE=5
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE=6
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE=7
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_BIHARMONIC_EQUATION_TYPE=8
  !Bioelectrics class
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE=2
  !Modal class
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LINEAR_ELASTIC_MODAL_TYPE=1
  
  !Problem subtypes
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NO_SUBTYPE=0
  !Elasticity class
  !  Linear elasticity
  !  Finite elasticity
  !Fluid mechanics class
  !Electromagnetics class
  !Classical field class
  !  Laplace equation
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE=2
  !  Poisson equation
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_CONSTANT_SOURCE_POISSON_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE=2
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE=3
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE=4
  !  Helmholtz equation
  !  Wave equation
  !  Diffusion equation
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE=1
  !Bioelectrics class
  !  Monodomain equation
  !  Bidomain equation
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE=1
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE=2
  !Modal class


  !Module parameters
  !> \addtogroup EQUATIONS_SET_CONSTANTS_SetupTypes EQUATIONS_SET_CONSTANTS::SetupTypes
  !> \brief Setup type parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_INITIAL_TYPE=1 !<Initial setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_GEOMETRY_TYPE=2 !<Geometry setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_DEPENDENT_TYPE=3 !<Dependent variables. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_MATERIALS_TYPE=4 !<Materials setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_SOURCE_TYPE=5 !<Source setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_SOURCE_MATERIALS_TYPE=6 !<Source materials setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_ANALYTIC_TYPE=7 !<Analytic setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_FIXED_CONDITIONS_TYPE=8 !<Fixed conditions. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_EQUATIONS_TYPE=9 !<Equations setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_FINAL_TYPE=9 !<Final setup. \see EQUATIONS_SET_CONSTANTS_SetupTypes,EQUATIONS_SET_CONSTANTS
  !>@}
  
  !> \addtogroup EQUATIONS_SET_CONSTANTS_SetupActionTypes EQUATIONS_SET_CONSTANTS::SetupActionTypes
  !> \brief Setup action type parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_START_ACTION=1 !<Start setup action. \see EQUATIONS_SET_CONSTANTS_SetupActionTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SETUP_FINISH_ACTION=2 !<Finish setup action. \see EQUATIONS_SET_CONSTANTS_SetupActionTypes,EQUATIONS_SET_CONSTANTS
  !>@}

  !> \addtogroup EQUATIONS_SET_CONSTANTS_FixedConditions EQUATIONS_SET_CONSTANTS::FixedConditions
  !> \brief Fixed conditons parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NOT_FIXED=0 !<The dof is not fixed. \see EQUATIONS_SET_CONSTANTS_FixedConditions,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FIXED_BOUNDARY_CONDITION=1 !<The dof is fixed as a boundary condition. \see EQUATIONS_SET_CONSTANTS_FixedConditions,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_MIXED_BOUNDARY_CONDITION=2 !<The dof is set as a mixed boundary condition. \see EQUATIONS_SET_CONSTANTS_FixedConditions,EQUATIONS_SET_CONSTANTS
  !>@}
  
  !> \addtogroup EQUATIONS_SET_CONSTANTS_LinearityTypes EQUATIONS_SET_CONSTANTS::LinearityTypes
  !> \brief The problem linearity type parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_EQUATIONS_SET_LINEARITY_TYPES=3 !<The number of problem linearity types defined. \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LINEAR=1 !<The problem is linear. \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NONLINEAR=2 !<The problem is non-linear. \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NONLINEAR_BCS=3 !<The problem has non-linear boundary conditions. \see EQUATIONS_SET_CONSTANTS_LinearityTypes,EQUATIONS_SET_CONSTANTS
  !>@}

  !> \addtogroup EQUATIONS_SET_CONSTANTS_JacobianCalculationTypes EQUATIONS_SET_CONSTANTS::JacobianCalculationTypes
  !> \brief The Jacobian types for nonlinear equations sets
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_JACOBIAN_NOT_CALCULATED=1 !<The Jacobian values will not be calculated for the nonlinear equations set \see EQUATIONS_SET_CONSTANTS_JacobianCalculationTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_JACOBIAN_ANALTYIC_CALCULATED=2 !<The Jacobian values will be calculated analytically for the nonlinear equations set \see EQUATIONS_SET_CONSTANTS_JacobianCalculationTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_JACOBIAN_FD_CALCULATED=3 !<The Jacobian values will be calcualted using finite differences for the nonlinear equations set \see EQUATIONS_SET_CONSTANTS_JacobianCalculationTypes,EQUATIONS_SET_CONSTANTS
  !>@}
  
  
  !> \addtogroup EQUATIONS_SET_CONSTANTS_TimeDepedenceTypes EQUATIONS_SET_CONSTANTS::TimeDepedenceTypes
  !> \brief The problem time dependence type parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_EQUATIONS_SET_TIME_TYPES=4 !<The number of equations set time dependence types defined. \see EQUATIONS_SET_CONSTANTS_TimeDependenceTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_STATIC=1 !<The equations set is static and has no time dependence \see EQUATIONS_SET_CONSTANTS_TimeDependenceTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_QUASISTATIC=2 !<The equations set is quasi-static \see EQUATIONS_SET_CONSTANTS_TimeDependenceTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FIRST_ORDER_DYNAMIC=3 !<The equations set is a first order dynamic set \see EQUATIONS_SET_CONSTANTS_TimeDependenceTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SECOND_ORDER_DYNAMIC=4 !<The equations set is a second order dynamic set \see EQUATIONS_SET_CONSTANTS_TimeDependenceTypes,EQUATIONS_SET_CONSTANTS
  !>@}

  !> \addtogroup EQUATIONS_SET_CONSTANTS_SolutionMethods EQUATIONS_SET_CONSTANTS::SolutionMethods
  !> \brief The solution method parameters
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_EQUATIONS_SET_SOLUTION_METHODS=6 !<The number of solution methods defined. \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FEM_SOLUTION_METHOD=1 !<Finite Element Method solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_BEM_SOLUTION_METHOD=2 !<Boundary Element Method solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FD_SOLUTION_METHOD=3 !<Finite Difference solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FV_SOLUTION_METHOD=4 !<Finite Volume solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_GFEM_SOLUTION_METHOD=5 !<Grid-based Finite Element Method solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_GFV_SOLUTION_METHOD=6 !<Grid-based Finite Volume solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
  !>@}

  !> \addtogroup EQUATIONS_SET_CONSTANTS_OutputTypes EQUATIONS_SET_CONSTANTS::OutputTypes
  !> \brief Equations output types
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_NO_OUTPUT=0 !<No output \see EQUATIONS_SET_CONSTANTS_OutputTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_TIMING_OUTPUT=1 !<Timing information output \see EQUATIONS_SET_CONSTANTS_OutputTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_MATRIX_OUTPUT=2 !<All below and equation matrices output \see EQUATIONS_SET_CONSTANTS_OutputTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_ELEMENT_MATRIX_OUTPUT=3 !<All below and Element matrices output \see EQUATIONS_SET_CONSTANTS_OutputTypes,EQUATIONS_SET_CONSTANTS
  !>@}

  !> \addtogroup EQUATIONS_SET_CONSTANTS_EquationsMatrixTypes EQUATIONS_SET_CONSTANTS::EquationsMatrixTypes
  !> \brief Equations matrix types
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LINEAR_MATRIX=1 !<The equations matrix is a linear equations set matrix \see EQUATIONS_SET_CONSTANTS_EquationsMatrixTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_DYNAMIC_MATRIX=2 !<The equations matrix is a dynamic equaitons set matrix \see EQUATIONS_SET_CONSTANTS_EquationsMatrixTypes,EQUATIONS_SET_CONSTANTS
 !>@}
 
  !> \addtogroup EQUATIONS_SET_CONSTANTS_SparsityTypes EQUATIONS_SET_CONSTANTS::SparsityTypes
  !> \brief Equations matrices sparsity types
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_SPARSE_MATRICES=1 !<Use sparse matrices for the equations set \see EQUATIONS_SET_CONSTANTS_SparsityTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_FULL_MATRICES=2 !<Use fully populated matrices for the equations set \see EQUATIONS_SET_CONSTANTS_SparsityTypes,EQUATIONS_SET_CONSTANTS
 !>@}
 
  !> \addtogroup EQUATIONS_SET_CONSTANTS_LumpingTypes EQUATIONS_SET_CONSTANTS::LumpingTypes
  !> \brief Equations matrices lumping types
  !> \see EQUATIONS_SET_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_UNLUMPED_MATRICES=1 !<The matrix is not lumped \see EQUATIONS_SET_CONSTANTS_SparsityTypes,EQUATIONS_SET_CONSTANTS
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LUMPED_MATRICES=2 !<The matrix is "mass" lumped \see EQUATIONS_SET_CONSTANTS_SparsityTypes,EQUATIONS_SET_CONSTANTS
 !>@}
 
  !> \addtogroup EQUATIONS_SET_CONSTANTS_EquationAnalyticFunction EQUATIONS_SET_CONSTANTS::EquationAnalyticFunction
  !> \brief the analytic representive of laplace equation
  !> \see LAPLACE_EQUATIONS_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1=1 !<u=x**2+2*x*y-y**2 
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2=2 !<u=cos(x)cosh(y)
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1=3 !<u=x**2-2*y**2+z**2
  INTEGER(INTG), PARAMETER :: EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2=4 !<u=cos(x)*cosh(y)*z
 !>@}
  
END MODULE EQUATIONS_SET_CONSTANTS
