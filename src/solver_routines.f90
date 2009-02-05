!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all solver routines.
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

!> This module handles all solver routines.
MODULE SOLVER_ROUTINES

  USE BASE_ROUTINES
  USE CMISS_PETSC
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE SOLVER_MAPPING_ROUTINES
  USE SOLVER_MATRICES_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup SOLVER_ROUTINES_SolverTypes SOLVER_ROUTINES::SolverTypes
  !> \brief The types of a problem solver
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_TYPE=1 !<A linear solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_TYPE=2 !<A nonlinear solver  \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_TYPE=3 !<A dynamic solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_INTEGRATION_TYPE=4 !<A integration solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EIGENPROBLEM_TYPE=5 !<A eigenproblem type \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SolverLibraries SOLVER_ROUTINES::SolverLibraries
  !> \brief The types of solver libraries
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_CMISS_LIBRARY=LIBRARY_CMISS_TYPE !<CMISS (internal) solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_PETSC_LIBRARY=LIBRARY_PETSC_TYPE !<PETSc solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_LinearSolverTypes SOLVER_ROUTINES::LinearSolverTypes
  !> \brief The types of linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_DIRECT_SOLVE_TYPE=1 !<Direct linear solver type \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE=2 !<Iterative linear solver type \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DirectLinearSolverTypes SOLVER_ROUTINES::DirectLinearSolverTypes
  !> \brief The types of direct linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_LU=1 !<LU direct linear solver \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_CHOLESKY=2 !<Cholesky direct linear solver \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_SVD=3 !<SVD direct linear solver \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_IterativeLinearSolverTypes SOLVER_ROUTINES::IterativeLinearSolverTypes
  !> \brief The types of iterative linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_RICHARDSON=1 !<Richardson iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CHEBYCHEV=2 !<Chebychev iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CONJUGATE_GRADIENT=3 !<Conjugate gradient iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BICONJUGATE_GRADIENT=4 !<Bi-conjugate gradient iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_GMRES=5 !<Generalised minimum residual iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BiCGSTAB=6 !<Stabalised bi-conjugate gradient iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_CONJGRAD_SQUARED=7 !<Conjugate gradient squared iterative solver type \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_IterativePreconditionerTypes SOLVER_ROUTINES::IterativePreconditionerTypes
  !> \brief The types of iterative preconditioners
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_NO_PRECONDITIONER=0 !<No preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_JACOBI_PRECONDITIONER=1 !<Jacobi preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER=2 !<Iterative block Jacobi preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_SOR_PRECONDITIONER=3 !<Successive over relaxation preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER=4 !<Incomplete Cholesky preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER=5 !<Incomplete LU preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER=6 !<Additive Schwrz preconditioner type \see SOLVER_ROUTINES_IterativePreconditionerTypes,SOLVER_ROUTINES
  !>@}


  !> \addtogroup SOLVER_ROUTINES_NonlinearSolverTypes SOLVER_ROUTINES::NonlinearSolverTypes
  !> \brief The types of nonlinear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_NEWTON=1 !<Newton nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_BFGS_INVERSE=2 !<BFGS inverse nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_SQP=3 !<Sequential Quadratic Program nonlinear solver type \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_NewtonSolverTypes SOLVER_ROUTINES::NewtonSolverTypes
  !> \brief The types of nonlinear Newton solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH=1 !<Newton line search nonlinear solver type \see SOLVER_ROUTINES_NewtonSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_TRUSTREGION=2 !<Newton trust region nonlinear solver type \see SOLVER_ROUTINES_NewtonSolverTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_NewtonLineSearchTypes SOLVER_ROUTINES::NewtonLineSearchTypes
  !> \brief The types line search techniques for Newton line search nonlinear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_NONORMS=1 !<No norms line search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_NONE=2 !<No line search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_QUADRATIC=3 !<Quadratic search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_LINESEARCH_CUBIC=4!<Cubic search for Newton line search nonlinear solves \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_JacobianCalculationTypes SOLVER_ROUTINES::JacobianCalculationTypes
  !> \brief The Jacobian calculation types for a nonlinear solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED=1 !<The Jacobian values will not be calculated for the nonlinear equations set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_ANALTYIC_CALCULATED=2 !<The Jacobian values will be calculated analytically for the nonlinear equations set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NEWTON_JACOBIAN_FD_CALCULATED=3 !<The Jacobian values will be calcualted using finite differences for the nonlinear equations set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
  !>@}  

   !> \addtogroup SOLVER_ROUTINES_DynamicOrderTypes SOLVER_ROUTINES::DynamicOrderTypes
  !> \brief The order types for a dynamic solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_FIRST_ORDER=1 !<Dynamic solver has first order terms \see SOLVER_ROUTINES_DynamicOrderTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_ORDER=2 !<Dynamic solver has second order terms \see SOLVER_ROUTINES_DynamicOrderTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DynamicLinearityTypes SOLVER_ROUTINES::DynamicLinearityTypes
  !> \brief The time linearity types for a dynamic solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_LINEAR=1 !<Dynamic solver has linear terms \see SOLVER_ROUTINES_DynamicLinearityTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NONLINEAR=2 !<Dynamic solver has nonlinear terms \see SOLVER_ROUTINES_DynamicLinearityTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_DynamicDegreeTypes SOLVER_ROUTINES::DynamicDegreeTypes
  !> \brief The time interpolation polynomial degree types for a dynamic solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_FIRST_DEGREE=1 !<Dynamic solver uses a first degree polynomial for time interpolation \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE=2 !<Dynamic solver uses a second degree polynomial for time interpolation \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE=3 !<Dynamic solver uses a third degree polynomial for time interpolation \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
  !>@}    
  
  !> \addtogroup SOLVER_ROUTINES_DynamicSchemeTypes SOLVER_ROUTINES::DynamicSchemeTypes
  !> \brief The types of dynamic solver scheme
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_EULER_SCHEME=1 !<Euler (explicit) dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME=2 !<Backward Euler (implicit) dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_CRANK_NICHOLSON_SCHEME=3 !<Crank-Nicholson dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_GALERKIN_SCHEME=4 !<Galerkin dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_ZLAMAL_SCHEME=5 !<Zlamal dynamic solver \see SOLVER_ROUTINES_DynamicorderTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME=6 !<2nd degree Gear dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME=7 !<1st 2nd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME=8 !<2nd 2nd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK1_SCHEME=9 !<1st Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK2_SCHEME=10 !<2nd Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_NEWMARK3_SCHEME=11 !<3rd Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME=12 !<3rd degree Gear dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME=13 !<1st 3rd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME=14 !<2nd 3rd degree Liniger dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HOUBOLT_SCHEME=15 !<Houbolt dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_WILSON_SCHEME=16 !<Wilson dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME=17 !<1st Bossak-Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME=18 !<2nd Bossak-Newmark dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME=19 !<1st Hilbert-Hughes-Taylor dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME=20 !<1st Hilbert-Hughes-Taylor dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_DYNAMIC_USER_DEFINED_SCHEME=21 !<User specified degree and theta dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_IntegratorTypes SOLVER_ROUTINES::IntegratorTypes
  !> \brief The integration types for a integration solver 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_INTEGRATION_EULER=1 !<Euler integrator \see SOLVER_ROUTINES_IntegratorType,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_INTEGRATION_IMPROVED_EULER=2 !<Improved Euler integrator \see SOLVER_ROUTINES_IntegratorType,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_INTEGRATION_4TH_RUNGE_KUTTA=3 !<4the order Runge-Kutta integrator \see SOLVER_ROUTINES_IntegratorType,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_INTEGRATION_ADAMS_MOULTON=4 !<Adams-Moulton integrator \see SOLVER_ROUTINES_IntegratorType,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_INTEGRATION_LSODA=5 !<LSODA integrator \see SOLVER_ROUTINES_IntegratorType,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_OutputTypes SOLVER_ROUTINES::OutputTypes
  !> \brief The types of output
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NO_OUTPUT=0 !<No output from the solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_PROGRESS_OUTPUT=1 !<Progress output from solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_TIMING_OUTPUT=2 !<Timing output from the solver routines plus below \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SOLVER_OUTPUT=3 !<Solver specific output from the solver routines plus below \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRIX_OUTPUT=4 !<SolVER matrices output from the solver routines plus below\see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SparsityTypes SOLVER_ROUTINES::SparsityTypes
  !> \brief The types of sparse solver matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SelectMatricesTypes SOLVER_ROUTINES::SelectMatricesTypes
  !> \brief The types of selection available for the solver matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_ALL=1 !<Select all the solver matrices and vectors \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_LINEAR_ONLY=2 !<Select only the linear solver matrices and vectors \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_NONLINEAR_ONLY=3 !<Select only the nonlinear solver matrices and vectors \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_JACOBIAN_ONLY=4 !<Select only the Jacobian solver matrix \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RESIDUAL_ONLY=5 !<Select only the residual solver vector \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_ONLY=6 !<Select only the RHS solver vector \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRICES_RHS_RESIDUAL_ONLY=7 !<Select only the residual and RHS solver vectors \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
  !>@}
  
  !> \addtogroup SOLVER_ROUTINES_EquationsLinearityTypes SOLVER_ROUTINES::EquationsLinearityTypes
  !> \brief The solver equations linearity types 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_LINEAR=1 !<Solver equations are linear \see SOLVER_ROUTINES_EquationLinearityTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_NONLINEAR=2 !<Solver equations are nonlinear \see SOLVER_ROUTINES_EquationLinearityTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_EquationsTimeDependenceTypes SOLVER_ROUTINES::EquationsTimeDependenceTypes
  !> \brief The solver equations time dependence types 
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_STATIC=1 !<Solver equations are static \see SOLVER_ROUTINES_EquationTimeDependenceTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_QUASISTATIC=2 !<Solver equations are quasistatic \see SOLVER_ROUTINES_EquationTimeDependenceTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC=3 !<Solver equations are first order dynamic \see SOLVER_ROUTINES_EquationTimeDependenceTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC=4 !<Solver equations are second order dynamic \see SOLVER_ROUTINES_EquationTimeDependenceTypes,SOLVER_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  INTERFACE SOLVER_DYNAMIC_THETA_SET
    MODULE PROCEDURE SOLVER_DYNAMIC_THETA_SET_DP1
    MODULE PROCEDURE SOLVER_DYNAMIC_THETA_SET_DP
  END INTERFACE !SOLVER_DYNAMIC_THETA_SET  

  PUBLIC SOLVER_LINEAR_TYPE,SOLVER_NONLINEAR_TYPE,SOLVER_DYNAMIC_TYPE

  PUBLIC SOLVER_CMISS_LIBRARY,SOLVER_PETSC_LIBRARY

  PUBLIC SOLVER_LINEAR_DIRECT_SOLVE_TYPE,SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
  
  PUBLIC SOLVER_DIRECT_LU,SOLVER_DIRECT_CHOLESKY,SOLVER_DIRECT_SVD

  PUBLIC SOLVER_ITERATIVE_RICHARDSON,SOLVER_ITERATIVE_CHEBYCHEV,SOLVER_ITERATIVE_CONJUGATE_GRADIENT, &
    & SOLVER_ITERATIVE_BICONJUGATE_GRADIENT,SOLVER_ITERATIVE_GMRES,SOLVER_ITERATIVE_BiCGSTAB,SOLVER_ITERATIVE_CONJGRAD_SQUARED

  PUBLIC SOLVER_ITERATIVE_NO_PRECONDITIONER,SOLVER_ITERATIVE_JACOBI_PRECONDITIONER,SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER, &
    & SOLVER_ITERATIVE_SOR_PRECONDITIONER,SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER, &
    & SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER,SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER

  PUBLIC SOLVER_NONLINEAR_NEWTON

  PUBLIC SOLVER_NEWTON_LINESEARCH,SOLVER_NEWTON_TRUSTREGION

  PUBLIC SOLVER_NEWTON_LINESEARCH_NONORMS,SOLVER_NEWTON_LINESEARCH_NONE,SOLVER_NEWTON_LINESEARCH_QUADRATIC, &
    & SOLVER_NEWTON_LINESEARCH_CUBIC

  PUBLIC SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED,SOLVER_NEWTON_JACOBIAN_ANALTYIC_CALCULATED, &
    & SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
  
  PUBLIC SOLVER_DYNAMIC_LINEAR,SOLVER_DYNAMIC_NONLINEAR

  PUBLIC SOLVER_DYNAMIC_FIRST_ORDER,SOLVER_DYNAMIC_SECOND_ORDER

  PUBLIC SOLVER_DYNAMIC_FIRST_DEGREE,SOLVER_DYNAMIC_SECOND_DEGREE,SOLVER_DYNAMIC_THIRD_DEGREE

  PUBLIC SOLVER_DYNAMIC_EULER_SCHEME,SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,SOLVER_DYNAMIC_CRANK_NICHOLSON_SCHEME, &
    & SOLVER_DYNAMIC_GALERKIN_SCHEME,SOLVER_DYNAMIC_ZLAMAL_SCHEME,SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME, &
    & SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME,SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME, &
    & SOLVER_DYNAMIC_NEWMARK1_SCHEME,SOLVER_DYNAMIC_NEWMARK2_SCHEME,SOLVER_DYNAMIC_NEWMARK3_SCHEME, &
    & SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME,SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME, &
    & SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME,SOLVER_DYNAMIC_HOUBOLT_SCHEME,SOLVER_DYNAMIC_WILSON_SCHEME, &
    & SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME,SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME, &
    & SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME,SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME, &
    & SOLVER_DYNAMIC_USER_DEFINED_SCHEME

  PUBLIC SOLVER_INTEGRATION_EULER,SOLVER_INTEGRATION_IMPROVED_EULER,SOLVER_INTEGRATION_4TH_RUNGE_KUTTA, &
    & SOLVER_INTEGRATION_ADAMS_MOULTON,SOLVER_INTEGRATION_LSODA
  
  PUBLIC SOLVER_NO_OUTPUT,SOLVER_PROGRESS_OUTPUT,SOLVER_TIMING_OUTPUT,SOLVER_SOLVER_OUTPUT,SOLVER_MATRIX_OUTPUT
  
  PUBLIC SOLVER_SPARSE_MATRICES,SOLVER_FULL_MATRICES

  PUBLIC SOLVER_EQUATIONS_LINEAR,SOLVER_EQUATIONS_NONLINEAR

  PUBLIC SOLVER_EQUATIONS_STATIC,SOLVER_EQUATIONS_QUASISTATIC,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC, &
    & SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC

  PUBLIC SOLVER_OUTPUT_TYPE_SET
  
  PUBLIC SOLVER_LIBRARY_SET,SOLVER_SOLVE

  PUBLIC SOLVER_DYNAMIC_DEGREE_SET,SOLVER_DYNAMIC_MONITOR,SOLVER_DYNAMIC_ORDER_SET,SOLVER_DYNAMIC_SCHEME_SET, &
    & SOLVER_DYNAMIC_THETA_SET,SOLVER_DYNAMIC_TIMES_SET

  PUBLIC SOLVER_EQUATIONS_CREATE_FINISH,SOLVER_EQUATIONS_CREATE_START,SOLVER_EQUATIONS_DESTROY, &
    & SOLVER_EQUATIONS_EQUATIONS_SET_ADD,SOLVER_EQUATIONS_LINEARITY_TYPE_SET,SOLVER_EQUATIONS_SPARSITY_TYPE_SET, &
    & SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET
  
  PUBLIC SOLVER_LINEAR_TYPE_SET
  
  PUBLIC SOLVER_LINEAR_DIRECT_TYPE_SET

  PUBLIC SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET,SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET, &
    & SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET,SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET, &
    & SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET,SOLVER_LINEAR_ITERATIVE_TYPE_SET

  PUBLIC SOLVER_MATRICES_ALL,SOLVER_MATRICES_LINEAR_ONLY,SOLVER_MATRICES_NONLINEAR_ONLY,SOLVER_MATRICES_JACOBIAN_ONLY, &
    & SOLVER_MATRICES_RESIDUAL_ONLY,SOLVER_MATRICES_RHS_ONLY,SOLVER_MATRICES_RHS_RESIDUAL_ONLY

  PUBLIC SOLVER_MATRICES_DYNAMIC_ASSEMBLE,SOLVER_MATRICES_STATIC_ASSEMBLE

  PUBLIC SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET

  PUBLIC SOLVER_NONLINEAR_MONITOR
  
  PUBLIC SOLVER_SOLVER_EQUATIONS_GET,SOLVER_TYPE_SET,SOLVER_VARIABLES_UPDATE

  PUBLIC SOLVERS_CREATE_FINISH,SOLVERS_CREATE_START,SOLVERS_DESTROY,SOLVERS_NUMBER_SET,SOLVERS_SOLVER_GET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver 
  SUBROUTINE SOLVER_CREATE_FINISH(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        !Set the finished flag. The final solver finish will be done once the solver equations have been finished.
        SOLVER%SOLVER_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a dynamic solver 
  SUBROUTINE SOLVER_DYNAMIC_CREATE_FINISH(DYNAMIC_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER !<A pointer to the dynamic solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_DYNAMIC_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN

    ELSE
      CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_DYNAMIC_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_DYNAMIC_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_DYNAMIC_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the degree of the polynomial used to interpolate time for a dynamic solver.
  SUBROUTINE SOLVER_DYNAMIC_DEGREE_SET(SOLVER,DEGREE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the dynamic solver to set the theta value for
    INTEGER(INTG), INTENT(IN) :: DEGREE !<The degree of the polynomial used for time interpolation in a dynamic solver \see SOLVER_ROUTINES_DynamicDegreeTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: degree_idx
    REAL(DP), ALLOCATABLE :: OLD_THETA(:)
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_DYNAMIC_DEGREE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("The solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
          IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
            IF(DEGREE/=DYNAMIC_SOLVER%DEGREE) THEN
              IF(DEGREE>=DYNAMIC_SOLVER%ORDER) THEN
                SELECT CASE(DEGREE)
                CASE(SOLVER_DYNAMIC_FIRST_DEGREE,SOLVER_DYNAMIC_SECOND_DEGREE,SOLVER_DYNAMIC_THIRD_DEGREE)
                  ALLOCATE(OLD_THETA(DYNAMIC_SOLVER%DEGREE),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old theta.",ERR,ERROR,*999)
                  OLD_THETA=DYNAMIC_SOLVER%THETA
                  IF(ALLOCATED(DYNAMIC_SOLVER%THETA)) DEALLOCATE(DYNAMIC_SOLVER%THETA)
                  ALLOCATE(DYNAMIC_SOLVER%THETA(DEGREE),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate theta.",ERR,ERROR,*999)
                  IF(DEGREE>DYNAMIC_SOLVER%DEGREE) THEN
                    DO degree_idx=1,DYNAMIC_SOLVER%DEGREE
                      DYNAMIC_SOLVER%THETA(degree_idx)=OLD_THETA(degree_idx)
                    ENDDO !degree_idx
                    DO degree_idx=DYNAMIC_SOLVER%DEGREE+1,DEGREE
                      DYNAMIC_SOLVER%THETA(degree_idx)=1.0_DP
                    ENDDO !degree_idx
                  ELSE
                    DO degree_idx=1,DEGREE
                      DYNAMIC_SOLVER%THETA(degree_idx)=OLD_THETA(degree_idx)
                    ENDDO !degree_idx
                  ENDIF
                  IF(ALLOCATED(OLD_THETA)) DEALLOCATE(OLD_THETA)
                  DYNAMIC_SOLVER%DEGREE=DEGREE
                CASE DEFAULT
                  LOCAL_ERROR="The specified degree of "//TRIM(NUMBER_TO_VSTRING(DEGREE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                LOCAL_ERROR="Invalid dynamic solver setup. The specfied degree of "// &
                  & TRIM(NUMBER_TO_VSTRING(DEGREE,"*",ERR,ERROR))//" must be >= the current dynamic order of "// &
                  & TRIM(NUMBER_TO_VSTRING(DYNAMIC_SOLVER%ORDER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified solver is not a dynamic solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_DYNAMIC_DEGREE_SET")
    RETURN
999 IF(ALLOCATED(OLD_THETA)) DEALLOCATE(OLD_THETA)
    CALL ERRORS("SOLVER_DYNAMIC_DEGREE_SET",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_DEGREE_SET")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_DEGREE_SET

  !
  !================================================================================================================================
  !

  !>Finalise a dynamic solver and deallocates all memory
  SUBROUTINE SOLVER_DYNAMIC_FINALISE(DYNAMIC_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER !<A pointer the dynamic solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_DYNAMIC_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
      IF(ALLOCATED(DYNAMIC_SOLVER%THETA)) DEALLOCATE(DYNAMIC_SOLVER%THETA)
      CALL PETSC_TSFINALISE(DYNAMIC_SOLVER%TS,ERR,ERROR,*999)
      DEALLOCATE(DYNAMIC_SOLVER)
    ENDIF
        
    CALL EXITS("SOLVER_DYNAMIC_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_DYNAMIC_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_DYNAMIC_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a dynamic solver for a problem solver
  SUBROUTINE SOLVER_DYNAMIC_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the dynamic solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_DYNAMIC_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%DYNAMIC_SOLVER)) THEN
        CALL FLAG_ERROR("Dynamic solver is already associated for this solver.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(SOLVER%DYNAMIC_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver dynamic solver.",ERR,ERROR,*999)
        SOLVER%DYNAMIC_SOLVER%SOLVER=>SOLVER
        SOLVER%DYNAMIC_SOLVER%SOLVER_LIBRARY=SOLVER_CMISS_LIBRARY
        SOLVER%DYNAMIC_SOLVER%LINEARITY=SOLVER_DYNAMIC_LINEAR
        SOLVER%DYNAMIC_SOLVER%ORDER=SOLVER_DYNAMIC_FIRST_ORDER
        SOLVER%DYNAMIC_SOLVER%DEGREE=SOLVER_DYNAMIC_FIRST_DEGREE
        SOLVER%DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_CRANK_NICHOLSON_SCHEME
        ALLOCATE(SOLVER%DYNAMIC_SOLVER%THETA(1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate theta.",ERR,ERROR,*999)
        SOLVER%DYNAMIC_SOLVER%THETA(1)=1.0_DP/2.0_DP
        SOLVER%DYNAMIC_SOLVER%EXPLICIT=.FALSE.
        SOLVER%DYNAMIC_SOLVER%CURRENT_TIME=0.0_DP
        SOLVER%DYNAMIC_SOLVER%TIME_INCREMENT=0.01_DP
        CALL PETSC_TSINITIALISE(SOLVER%DYNAMIC_SOLVER%TS,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_DYNAMIC_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_DYNAMIC_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_DYNAMIC_INITIALISE

  !
  !================================================================================================================================
  !

  !>Monitors the dynamic solve.
  SUBROUTINE SOLVER_DYNAMIC_MONITOR(DYNAMIC_SOLVER,STEPS,TIME,ERR,ERROR,*)

   !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER !<A pointer to the dynamic solver to monitor
    INTEGER(INTG), INTENT(IN) :: STEPS !<The number of iterations
    REAL(DP), INTENT(IN) :: TIME !<The current time
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_DYNAMIC_MONITOR",ERR,ERROR,*999)

    IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
        
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Dynamic solve monitor: ",ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of steps = ",STEPS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Current time    = ",TIME,ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)      
        
    ELSE
      CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
    ENDIF
     
    CALL EXITS("SOLVER_DYNAMIC_MONITOR")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_MONITOR",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_MONITOR")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_MONITOR

  !
  !================================================================================================================================
  !

  !>Sets/changes the order for a dynamic solver.
  SUBROUTINE SOLVER_DYNAMIC_ORDER_SET(SOLVER,ORDER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the dynamic solver to set the theta value for
    INTEGER(INTG), INTENT(IN) :: ORDER !<The order of the dynamic solver \see SOLVER_ROUTINES_DynamicOrderTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_DYNAMIC_ORDER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("The solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
          IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
            IF(ORDER==SOLVER_DYNAMIC_SECOND_ORDER.AND.DYNAMIC_SOLVER%DEGREE==SOLVER_DYNAMIC_FIRST_DEGREE) THEN
              LOCAL_ERROR="Invalid dynamic solver degree. You must have at least a second degree polynomial "// &
                & "interpolation for a second order dynamic solver."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              SELECT CASE(ORDER)
              CASE(SOLVER_DYNAMIC_FIRST_ORDER)
                DYNAMIC_SOLVER%ORDER=SOLVER_DYNAMIC_FIRST_ORDER
              CASE(SOLVER_DYNAMIC_SECOND_ORDER)
                DYNAMIC_SOLVER%ORDER=SOLVER_DYNAMIC_SECOND_ORDER
              CASE DEFAULT
                LOCAL_ERROR="The specified order of "//TRIM(NUMBER_TO_VSTRING(ORDER,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
          ELSE
            CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified solver is not a dynamic solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_DYNAMIC_ORDER_SET")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_ORDER_SET",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_ORDER_SET")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_ORDER_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the scheme for a dynamic solver.
  SUBROUTINE SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SCHEME,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the dynamic solver to set the scheme for
    INTEGER(INTG), INTENT(IN) :: SCHEME !<The scheme used for a dynamic solver \see SOLVER_ROUTINES_DynamicSchemeTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: ALPHA,BETA,GAMMA,THETA
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_DYNAMIC_SCHEME_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("The solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
          IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
            SELECT CASE(SCHEME)
            CASE(SOLVER_DYNAMIC_EULER_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_EULER_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,0.0_DP,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,1.0_DP,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_CRANK_NICHOLSON_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_CRANK_NICHOLSON_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,1.0_DP/2.0_DP,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_GALERKIN_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_GALERKIN_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,2.0_DP/3.0_DP,ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_ZLAMAL_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_ZLAMAL_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/5.0_DP/6.0_DP,2.0_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/3.0_DP/2.0_DP,2.0_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.0848_DP,1.0_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.2184_DP,1.292_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_NEWMARK1_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_NEWMARK1_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              BETA=0.5_DP
              GAMMA=2.0_DP
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/GAMMA,2.0_DP*BETA/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_NEWMARK2_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_NEWMARK2_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              BETA=0.3025_DP
              GAMMA=0.6_DP
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/GAMMA,2.0_DP*BETA/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_NEWMARK3_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_NEWMARK3_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              BETA=0.25_DP
              GAMMA=0.5_DP
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/GAMMA,2.0_DP*BETA/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/2.0_DP,11.0_DP/3.0_DP,6.0_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.84_DP,3.07_DP,4.5_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/0.80_DP,1.03_DP,1.29_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_HOUBOLT_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_HOUBOLT_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/2.0_DP,11.0_DP/3.0_DP,6.0_DP/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_WILSON_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_WILSON_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              THETA=1.4_DP
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/THETA,THETA**2,THETA**3/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              ALPHA=-0.1_DP
              BETA=0.3025_DP
              GAMMA=0.5_DP-ALPHA
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.0_DP-ALPHA,2.0_DP/3.0_DP-ALPHA+2.0_DP*BETA,6.0_DP*BETA/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              ALPHA=-0.1_DP
              BETA=1.0_DP/6.0_DP-1.0_DP/2.0_DP*ALPHA
              GAMMA=1.0_DP/2.0_DP-ALPHA
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.0_DP-ALPHA,1.0_DP-2.0_DP*ALPHA,1.0_DP-3.0_DP*ALPHA/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              ALPHA=-0.1_DP
              BETA=0.3025_DP
              GAMMA=0.5_DP-ALPHA
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.0_DP,2.0_DP/3.0_DP+2.0_DP*BETA-2.0_DP*ALPHA**2, &
                & 6.0_DP*BETA*(1.0_DP+ALPHA)/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_THIRD_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_SECOND_ORDER,ERR,ERROR,*999)
              ALPHA=-0.3_DP
              BETA=0.3025_DP
              GAMMA=0.5_DP-ALPHA
              CALL SOLVER_DYNAMIC_THETA_SET(SOLVER,(/1.0_DP,2.0_DP/3.0_DP+2.0_DP*BETA-2.0_DP*ALPHA**2, &
                & 6.0_DP*BETA*(1.0_DP+ALPHA)/),ERR,ERROR,*999)
            CASE(SOLVER_DYNAMIC_USER_DEFINED_SCHEME)
              DYNAMIC_SOLVER%SCHEME=SOLVER_DYNAMIC_USER_DEFINED_SCHEME
            CASE DEFAULT
              LOCAL_ERROR="The specified scheme of "//TRIM(NUMBER_TO_VSTRING(SCHEME,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified solver is not a dynamic solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_DYNAMIC_SCHEME_SET")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_SCHEME_SET",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_SCHEME_SET")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_SCHEME_SET

  !
  !================================================================================================================================
  !

  !>Solve a dynamic solver 
  SUBROUTINE SOLVER_DYNAMIC_SOLVE(DYNAMIC_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER !<A pointer to the dynamic solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_DYNAMIC_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
      !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_DYNAMIC_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_DYNAMIC_SOLVE")
    RETURN 1
    
  END SUBROUTINE SOLVER_DYNAMIC_SOLVE
        
  !
  !================================================================================================================================
  !

  !>Sets/changes a single theta value for a dynamic solver.
  SUBROUTINE SOLVER_DYNAMIC_THETA_SET_DP1(SOLVER,THETA,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the dynamic solver to set the theta value for
    REAL(DP), INTENT(IN) :: THETA !<The theta value to set for the first degree polynomial
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   
    CALL ENTERS("SOLVER_DYNAMIC_THETA_SET_DP1",ERR,ERROR,*999)

    CALL SOLVER_DYNAMIC_THETA_SET_DP(SOLVER,(/THETA/),ERR,ERROR,*999)
    
    CALL EXITS("SOLVER_DYNAMIC_THETA_SET_DP1")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_THETA_SET_DP1",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_THETA_SET_DP1")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_THETA_SET_DP1

  !
  !================================================================================================================================
  !

  !>Sets/changes the theta value for a dynamic solver.
  SUBROUTINE SOLVER_DYNAMIC_THETA_SET_DP(SOLVER,THETA,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the dynamic solver to set the theta value for
    REAL(DP), INTENT(IN) :: THETA(:) !<THEATA(degree_idx). The theta value to set for the degree_idx-1'th polynomial
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: degree_idx
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_DYNAMIC_THETA_SET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("The solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
          IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
            IF(SIZE(THETA,1)>=DYNAMIC_SOLVER%DEGREE) THEN
              DO degree_idx=1,DYNAMIC_SOLVER%DEGREE
                IF(THETA(degree_idx)>=0.0_DP) THEN
                  DYNAMIC_SOLVER%THETA(degree_idx)=THETA(degree_idx)
                ELSE
                  LOCAL_ERROR="The specified theta "//TRIM(NUMBER_TO_VSTRING(degree_idx,"*",ERR,ERROR))// &
                    & " value of "//TRIM(NUMBER_TO_VSTRING(THETA(degree_idx),"*",ERR,ERROR))// &
                    & " is invalid. The theta value must be >= 0.0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !degree_idx
            ELSE
              LOCAL_ERROR="Invalid number of the thetas. The supplied number of thetas ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(THETA,1),"*",ERR,ERROR))//") must be equal to the interpolation degree ("// &
                & TRIM(NUMBER_TO_VSTRING(DYNAMIC_SOLVER%DEGREE,"*",ERR,ERROR))//")."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The specified solver is not a dynamic solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_DYNAMIC_THETA_SET_DP")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_THETA_SET_DP",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_THETA_SET_DP")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_THETA_SET_DP
  !
  !================================================================================================================================
  !

  !>Sets/changes the dynamic times for a dynamic solver.
  SUBROUTINE SOLVER_DYNAMIC_TIMES_SET(SOLVER,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the dynamic solver to set the times for
    REAL(DP), INTENT(IN) :: CURRENT_TIME !<The current time to set
    REAL(DP), INTENT(IN) :: TIME_INCREMENT !<The time increment to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_DYNAMIC_TIMES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      !Note: do not check for finished here as we may wish to modify this for multiple solves.
      IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
        DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
        IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
          IF(ABS(TIME_INCREMENT)<=ZERO_TOLERANCE) THEN
            LOCAL_ERROR="The specified time increment of "//TRIM(NUMBER_TO_VSTRING(TIME_INCREMENT,"*",ERR,ERROR))// &
              & " is invalid. The time increment must not be zero."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            DYNAMIC_SOLVER%CURRENT_TIME=CURRENT_TIME
            DYNAMIC_SOLVER%TIME_INCREMENT=TIME_INCREMENT
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dynamic solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The specified solver is not a dynamic solver.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
     
    CALL EXITS("SOLVER_DYNAMIC_TIMES_SET")
    RETURN
999 CALL ERRORS("SOLVER_DYNAMIC_TIMES_SET",ERR,ERROR)
    CALL EXITS("SOLVER_DYNAMIC_TIMES_SET")
    RETURN 1
  END SUBROUTINE SOLVER_DYNAMIC_TIMES_SET

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a eigenproblem solver 
  SUBROUTINE SOLVER_EIGENPROBLEM_CREATE_FINISH(EIGENPROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(EIGENPROBLEM_SOLVER_TYPE), POINTER :: EIGENPROBLEM_SOLVER !<A pointer to the eigenproblem solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_EIGENPROBLEM_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EIGENPROBLEM_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Eigenproblem solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_EIGENPROBLEM_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_EIGENPROBLEM_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_EIGENPROBLEM_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_EIGENPROBLEM_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finalise a eigenproblem solver for a problem solver
  SUBROUTINE SOLVER_EIGENPROBLEM_FINALISE(EIGENPROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(EIGENPROBLEM_SOLVER_TYPE), POINTER :: EIGENPROBLEM_SOLVER !<A pointer the eigenproblem solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_EIGENPROBLEM_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EIGENPROBLEM_SOLVER)) THEN        
      DEALLOCATE(EIGENPROBLEM_SOLVER)
    ENDIF
         
    CALL EXITS("SOLVER_EIGENPROBLEM_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_EIGENPROBLEM_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_EIGENPROBLEM_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_EIGENPROBLEM_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a eigenproblem solver for a problem solver
  SUBROUTINE SOLVER_EIGENPROBLEM_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the eigenproblem solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_EIGENPROBLEM_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%EIGENPROBLEM_SOLVER)) THEN
        CALL FLAG_ERROR("Eigenproblem solver is already associated for this solver.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER%EIGENPROBLEM_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver eigenproblem solver.",ERR,ERROR,*999)
        SOLVER%EIGENPROBLEM_SOLVER%SOLVER=>SOLVER
        SOLVER%EIGENPROBLEM_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_EIGENPROBLEM_INITIALISE")
    RETURN
999 CALL SOLVER_EIGENPROBLEM_FINALISE(SOLVER%EIGENPROBLEM_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_EIGENPROBLEM_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_EIGENPROBLEM_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_EIGENPROBLEM_INITIALISE

  !
  !================================================================================================================================
  !

  !>Solve a eigenproblem solver
  SUBROUTINE SOLVER_EIGENPROBLEM_SOLVE(EIGENPROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(EIGENPROBLEM_SOLVER_TYPE), POINTER :: EIGENPROBLEM_SOLVER !<A pointer the eigenproblem solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_EIGENPROBLEM_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(EIGENPROBLEM_SOLVER)) THEN        
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Eigenproblem solver is not associated.",ERR,ERROR,*999)
    ENDIF
         
    CALL EXITS("SOLVER_EIGENPROBLEM_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_EIGENPROBLEM_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_EIGENPROBLEM_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_EIGENPROBLEM_SOLVE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating solver equations
  SUBROUTINE SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("SOLVER_EQUATIONS_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished.",ERR,ERROR,*998)
      ELSE
        SOLVER=>SOLVER_EQUATIONS%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          !Finish of the solver mapping
          CALL SOLVER_MAPPING_CREATE_FINISH(SOLVER_EQUATIONS%SOLVER_MAPPING,ERR,ERROR,*999)
          !Now finish off with the solver specific actions
          SELECT CASE(SOLVER%SOLVE_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_LINEAR_CREATE_FINISH(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_NONLINEAR_CREATE_FINISH(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_DYNAMIC_TYPE)
            CALL SOLVER_DYNAMIC_CREATE_FINISH(SOLVER%DYNAMIC_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_INTEGRATION_TYPE)
            CALL SOLVER_INTEGRATION_CREATE_FINISH(SOLVER%INTEGRATION_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_EIGENPROBLEM_TYPE)
            CALL SOLVER_EIGENPROBLEM_CREATE_FINISH(SOLVER%EIGENPROBLEM_SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED=.TRUE.
        ELSE
          CALL FLAG_ERROR("Solver equations solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL SOLVER_EQUATIONS_FINALISE(SOLVER_EQUATIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_EQUATIONS_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating solver equations
  SUBROUTINE SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to start the creation of solver equations on
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<On return, A pointer the solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          CALL FLAG_ERROR("Solver equations is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(SOLVER_EQUATIONS)
          CALL SOLVER_EQUATIONS_INITIALISE(SOLVER,ERR,ERROR,*999)
          CALL SOLVER_MAPPING_CREATE_START(SOLVER%SOLVER_EQUATIONS,SOLVER_MAPPING,ERR,ERROR,*999)
          SELECT CASE(SOLVER%SOLVE_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,1,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,1,ERR,ERROR,*999)
          CASE(SOLVER_DYNAMIC_TYPE)
            CALL SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,1,ERR,ERROR,*999)
          CASE(SOLVER_INTEGRATION_TYPE)
            CALL SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,0,ERR,ERROR,*999)
          CASE(SOLVER_EIGENPROBLEM_TYPE)
            CALL SOLVER_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLVER_MAPPING,2,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_CREATE_START",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_CREATE_START")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_CREATE_START
        
  !
  !================================================================================================================================
  !

  !>Destroys the solver equations
  SUBROUTINE SOLVER_EQUATIONS_DESTROY(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      CALL SOLVER_EQUATIONS_FINALISE(SOLVER_EQUATIONS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_DESTROY",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_DESTROY")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_DESTROY
        
  !
  !================================================================================================================================
  !

  !>Adds equations sets to solver equations
  SUBROUTINE SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to add the equations set to.
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_INDEX !<On exit, the index of the equations set that has been added
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_EQUATIONS_EQUATIONS_SET_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN          
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            EQUATIONS=>EQUATIONS_SET%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              IF(EQUATIONS_SET%LINEARITY==SOLVER_EQUATIONS%LINEARITY) THEN
                IF(EQUATIONS_SET%TIME_DEPENDENCE==SOLVER_EQUATIONS%TIME_DEPENDENCE) THEN
                  CALL SOLVER_MAPPING_EQUATIONS_SET_ADD(SOLVER_MAPPING,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="Invalid equations set up. The time dependence of the equations set to add ("// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TIME_DEPENDENCE,"*",ERR,ERROR))// &
                    & ") does not match the solver equations time dependence ("// &
                    & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Invalid equations set up. The linearity of the equations set to add ("// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%LINEARITY,"*",ERR,ERROR))// &
                  & ") does not match the solver equations linearity ("// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%LINEARITY,"*",ERR,ERROR))//")."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_EQUATIONS_EQUATIONS_SET_ADD")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_EQUATIONS_SET_ADD",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_EQUATIONS_SET_ADD")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_EQUATIONS_SET_ADD
        
  !
  !================================================================================================================================
  !

  !>Finalises the solver equations and deallocates all memory.
  SUBROUTINE SOLVER_EQUATIONS_FINALISE(SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_EQUATIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MAPPING)) CALL SOLVER_MAPPING_DESTROY(SOLVER_EQUATIONS%SOLVER_MAPPING,ERR,ERROR,*999)
      IF(ASSOCIATED(SOLVER_EQUATIONS%SOLVER_MATRICES)) CALL SOLVER_MATRICES_DESTROY(SOLVER_EQUATIONS%SOLVER_MATRICES,ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_EQUATIONS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_FINALISE
        
  !
  !================================================================================================================================
  !

  !>Initialises the solver equations for a solver.
  SUBROUTINE SOLVER_EQUATIONS_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_EQUATIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%SOLVER_EQUATIONS)) THEN
        CALL FLAG_ERROR("Solver equations is already associated for this solver.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER%SOLVER_EQUATIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver equations.",ERR,ERROR,*999)
        SOLVER%SOLVER_EQUATIONS%SOLVER=>SOLVER
        SOLVER%SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED=.FALSE.
        SOLVER%SOLVER_EQUATIONS%SPARSITY_TYPE=SOLVER_SPARSE_MATRICES
        NULLIFY(SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING)
        NULLIFY(SOLVER%SOLVER_EQUATIONS%SOLVER_MATRICES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_EQUATIONS_INITIALISE")
    RETURN
999 CALL SOLVER_EQUATIONS_FINALISE(SOLVER%SOLVER_EQUATIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_EQUATIONS_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_INITIALISE
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for solver equations
  SUBROUTINE SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,LINEARITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: LINEARITY_TYPE !<The type of linearity to be set \see SOLVER_ROUTINES_EquationLinearityTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_EQUATIONS_LINEARITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(LINEARITY_TYPE)
        CASE(SOLVER_EQUATIONS_LINEAR)
          SOLVER_EQUATIONS%LINEARITY=SOLVER_EQUATIONS_LINEAR
        CASE(SOLVER_EQUATIONS_NONLINEAR)
          SOLVER_EQUATIONS%LINEARITY=SOLVER_EQUATIONS_NONLINEAR
        CASE DEFAULT
          LOCAL_ERROR="The specified solver equations linearity type of "// &
            & TRIM(NUMBER_TO_VSTRING(LINEARITY_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
   
    CALL EXITS("SOLVER_EQUATIONS_LINEARITY_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_LINEARITY_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_LINEARITY_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_LINEARITY_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for solver equations
  SUBROUTINE SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The type of solver equations sparsity to be set \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished.",ERR,ERROR,*999)
      ELSE
!!TODO: Maybe set the sparsity in the different types of solvers. e.g., a sparse integrator doesn't mean much.
        SELECT CASE(SPARSITY_TYPE)
        CASE(SOLVER_SPARSE_MATRICES)
          SOLVER_EQUATIONS%SPARSITY_TYPE=SOLVER_SPARSE_MATRICES
        CASE(SOLVER_FULL_MATRICES)
          SOLVER_EQUATIONS%SPARSITY_TYPE=SOLVER_FULL_MATRICES
        CASE DEFAULT
          LOCAL_ERROR="The specified solver equations sparsity type of "// &
            & TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_SPARSITY_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for solver equations
  SUBROUTINE SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,TIME_DEPENDENCE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<A pointer the solver equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: TIME_DEPENDENCE_TYPE !<The type of time dependence to be set \see SOLVER_ROUTINES_EquationTimeDependenceTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
      IF(SOLVER_EQUATIONS%SOLVER_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Solver equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TIME_DEPENDENCE_TYPE)
        CASE(SOLVER_EQUATIONS_STATIC)
          SOLVER_EQUATIONS%TIME_DEPENDENCE=SOLVER_EQUATIONS_STATIC
        CASE(SOLVER_EQUATIONS_QUASISTATIC)
          SOLVER_EQUATIONS%TIME_DEPENDENCE=SOLVER_EQUATIONS_QUASISTATIC
        CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC)
          SOLVER_EQUATIONS%TIME_DEPENDENCE=SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC
        CASE(SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
          SOLVER_EQUATIONS%TIME_DEPENDENCE=SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC
        CASE DEFAULT
          LOCAL_ERROR="The specified solver equations time dependence type of "// &
            & TRIM(NUMBER_TO_VSTRING(TIME_DEPENDENCE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
    ENDIF
   
    CALL EXITS("SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Finalises a problem solver and deallocates all memory.
  SUBROUTINE SOLVER_FINALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      CALL SOLVER_LINEAR_FINALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
      CALL SOLVER_NONLINEAR_FINALISE(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
      CALL SOLVER_DYNAMIC_FINALISE(SOLVER%DYNAMIC_SOLVER,ERR,ERROR,*999)        
      CALL SOLVER_INTEGRATION_FINALISE(SOLVER%INTEGRATION_SOLVER,ERR,ERROR,*999)        
      CALL SOLVER_EIGENPROBLEM_FINALISE(SOLVER%EIGENPROBLEM_SOLVER,ERR,ERROR,*999)
      CALL SOLVER_EQUATIONS_FINALISE(SOLVER%SOLVER_EQUATIONS,ERR,ERROR,*999)
      DEALLOCATE(SOLVER)
    ENDIF 
        
    CALL EXITS("SOLVER_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a solver for a control loop
  SUBROUTINE SOLVER_INITIALISE(SOLVERS,SOLVER_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer the solvers to initialise the solver for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The solver index in solvers to initialise the solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVERS)) THEN
      IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
        IF(ALLOCATED(SOLVERS%SOLVERS)) THEN
          IF(ASSOCIATED(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR)) THEN
            CALL FLAG_ERROR("Solver pointer is already associated for this solver index.",ERR,ERROR,*998)
          ELSE
            ALLOCATE(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver.",ERR,ERROR,*999)
            SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%SOLVERS=>SOLVERS
            SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%SOLVER_FINISHED=.FALSE.
            SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%OUTPUT_TYPE=SOLVER_NO_OUTPUT
            NULLIFY(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%LINEAR_SOLVER)
            NULLIFY(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%NONLINEAR_SOLVER)
            NULLIFY(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%DYNAMIC_SOLVER)
            NULLIFY(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%INTEGRATION_SOLVER)
            NULLIFY(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%EIGENPROBLEM_SOLVER)
            NULLIFY(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%SOLVER_EQUATIONS)
            !Default to a linear solver and initialise
            SOLVERS%SOLVERS(SOLVER_INDEX)%PTR%SOLVE_TYPE=SOLVER_LINEAR_TYPE
            CALL SOLVER_LINEAR_INITIALISE(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers solvers is not allocated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="The solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
          & " is invalid. The solver index must be > 0 and <= "// &
          & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_INITIALISE")
    RETURN
999 CALL SOLVER_FINALISE(SOLVERS%SOLVERS(SOLVER_INDEX)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a integration solver 
  SUBROUTINE SOLVER_INTEGRATION_CREATE_FINISH(INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTEGRATION_SOLVER_TYPE), POINTER :: INTEGRATION_SOLVER !<A pointer to the integration solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_INTEGRATION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTEGRATION_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Integration solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_INTEGRATION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_INTEGRATION_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_INTEGRATION_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_INTEGRATION_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finalise a integration solver and deallocate all memory
  SUBROUTINE SOLVER_INTEGRATION_FINALISE(INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTEGRATION_SOLVER_TYPE), POINTER :: INTEGRATION_SOLVER !<A pointer the intergration solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_INTEGRATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTEGRATION_SOLVER)) THEN        
      DEALLOCATE(INTEGRATION_SOLVER)
    ENDIF
         
    CALL EXITS("SOLVER_INTEGRATION_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_INTEGRATION_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_INTEGRATION_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_INTEGRATION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a integration solver for a solver
  SUBROUTINE SOLVER_INTEGRATION_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the integration solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_INTEGRATION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%INTEGRATION_SOLVER)) THEN
        CALL FLAG_ERROR("Integration solver is already associated for this solver.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(SOLVER%INTEGRATION_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver integration solver.",ERR,ERROR,*999)
        SOLVER%INTEGRATION_SOLVER%SOLVER=>SOLVER
        SOLVER%INTEGRATION_SOLVER%SOLVER_LIBRARY=SOLVER_CMISS_LIBRARY
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_INTEGRATION_INITIALISE")
    RETURN
999 CALL SOLVER_INTEGRATION_FINALISE(SOLVER%INTEGRATION_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_INTEGRATION_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_INTEGRATION_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_INTEGRATION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Solve a integration solver
  SUBROUTINE SOLVER_INTEGRATION_SOLVE(INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTEGRATION_SOLVER_TYPE), POINTER :: INTEGRATION_SOLVER !<A pointer the integration solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_INTEGRATION_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTEGRATION_SOLVER)) THEN        
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Integration solver is not associated.",ERR,ERROR,*999)
    ENDIF
         
    CALL EXITS("SOLVER_INTEGRATION_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_INTEGRATION_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_INTEGRATION_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_INTEGRATION_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for the solver
  SUBROUTINE SOLVER_LIBRARY_SET(SOLVER,SOLVER_LIBRARY,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the type of
    INTEGER(INTG), INTENT(IN) :: SOLVER_LIBRARY !<The type of library to use for the solver \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(EIGENPROBLEM_SOLVER_TYPE), POINTER :: EIGENPROBLEM_SOLVER
    TYPE(INTEGRATION_SOLVER_TYPE), POINTER :: INTEGRATION_SOLVER
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: DIRECT_SOLVER
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: ITERATIVE_SOLVER
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LIBRARY_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has alredy been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(SOLVER%SOLVE_TYPE)
        CASE(SOLVER_LINEAR_TYPE)
          LINEAR_SOLVER=>SOLVER%LINEAR_SOLVER
          IF(ASSOCIATED(LINEAR_SOLVER)) THEN
            SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVE_TYPE)
            CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
              DIRECT_SOLVER=>LINEAR_SOLVER%DIRECT_SOLVER
              IF(ASSOCIATED(DIRECT_SOLVER)) THEN
                SELECT CASE(SOLVER_LIBRARY)
                CASE(SOLVER_CMISS_LIBRARY)
                  DIRECT_SOLVER%SOLVER_LIBRARY=SOLVER_CMISS_LIBRARY
                CASE(SOLVER_PETSC_LIBRARY)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT                  
              ELSE
                CALL FLAG_ERROR("Linear solver direct solver is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
              ITERATIVE_SOLVER=>LINEAR_SOLVER%ITERATIVE_SOLVER
              IF(ASSOCIATED(ITERATIVE_SOLVER)) THEN
                SELECT CASE(SOLVER_LIBRARY)
                CASE(SOLVER_CMISS_LIBRARY)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(SOLVER_PETSC_LIBRARY)
                  ITERATIVE_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
                CASE DEFAULT
                  LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVE_TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(SOLVER_NONLINEAR_TYPE)
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            SELECT CASE(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE)
            CASE(SOLVER_NONLINEAR_NEWTON)
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                SELECT CASE(NEWTON_SOLVER%NEWTON_SOLVE_TYPE)
                CASE(SOLVER_NEWTON_LINESEARCH)
                  LINESEARCH_SOLVER=>NEWTON_SOLVER%LINESEARCH_SOLVER
                  IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
                    SELECT CASE(SOLVER_LIBRARY)
                    CASE(SOLVER_CMISS_LIBRARY)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(SOLVER_PETSC_LIBRARY)
                      LINESEARCH_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
                    CASE DEFAULT
                      LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("Newton line search solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE(SOLVER_NEWTON_TRUSTREGION)
                  TRUSTREGION_SOLVER=>NEWTON_SOLVER%TRUSTREGION_SOLVER
                  IF(ASSOCIATED(TRUSTREGION_SOLVER)) THEN
                    SELECT CASE(SOLVER_LIBRARY)
                    CASE(SOLVER_CMISS_LIBRARY)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(SOLVER_PETSC_LIBRARY)
                      TRUSTREGION_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
                    CASE DEFAULT
                      LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("Newton trust region solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                CASE DEFAULT
                  LOCAL_ERROR="The Newton solver type of "// &
                    & TRIM(NUMBER_TO_VSTRING(NEWTON_SOLVER%NEWTON_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(SOLVER_NONLINEAR_SQP)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The nonlinear solver type of "// &
                & TRIM(NUMBER_TO_VSTRING(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(SOLVER_DYNAMIC_TYPE)
          DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
          IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
            SELECT CASE(SOLVER_LIBRARY)
            CASE(SOLVER_CMISS_LIBRARY)
              DYNAMIC_SOLVER%SOLVER_LIBRARY=SOLVER_CMISS_LIBRARY
            CASE(SOLVER_PETSC_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver dynamic solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(SOLVER_INTEGRATION_TYPE)
          INTEGRATION_SOLVER=>SOLVER%INTEGRATION_SOLVER
          IF(ASSOCIATED(INTEGRATION_SOLVER)) THEN
            SELECT CASE(SOLVER_LIBRARY)
            CASE(SOLVER_CMISS_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(SOLVER_PETSC_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver integration solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(SOLVER_EIGENPROBLEM_TYPE)
          EIGENPROBLEM_SOLVER=>SOLVER%EIGENPROBLEM_SOLVER
          IF(ASSOCIATED(EIGENPROBLEM_SOLVER)) THEN
            SELECT CASE(SOLVER_LIBRARY)
            CASE(SOLVER_CMISS_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(SOLVER_PETSC_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solver library type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver eigenproblem solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
     ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LIBRARY_SET")
    RETURN
999 CALL ERRORS("SOLVER_LIBRARY_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LIBRARY_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LIBRARY_SET
  
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a linear solver 
  SUBROUTINE SOLVER_LINEAR_CREATE_FINISH(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linear solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVE_TYPE)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        CALL SOLVER_LINEAR_DIRECT_CREATE_FINISH(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        CALL SOLVER_LINEAR_ITERATIVE_CREATE_FINISH(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVE_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a linear direct solver 
  SUBROUTINE SOLVER_LINEAR_DIRECT_CREATE_FINISH(LINEAR_DIRECT_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: LINEAR_DIRECT_SOLVER !<A pointer to the linear direct solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_DIRECT_CREATE_FINISH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(LINEAR_DIRECT_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_DIRECT_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SELECT CASE(LINEAR_DIRECT_SOLVER%SOLVER_LIBRARY)
            CASE(SOLVER_CMISS_LIBRARY)
              !Create the solver matrices
              CALL SOLVER_MATRICES_CREATE_START(SOLVER_EQUATIONS,SOLVER_MATRICES,ERR,ERROR,*999)
              CALL SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              SELECT CASE(SOLVER_EQUATIONS%SPARSITY_TYPE)
              CASE(SOLVER_SPARSE_MATRICES)
                CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE(SOLVER_FULL_MATRICES)
                CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The specified solver equations sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL SOLVER_MATRICES_CREATE_FINISH(SOLVER_MATRICES,ERR,ERROR,*999)
            CASE(SOLVER_PETSC_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solver library type of "// &
                & TRIM(NUMBER_TO_VSTRING(LINEAR_DIRECT_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Linear solver solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Linear direct solver linear solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear direct solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_DIRECT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_DIRECT_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_DIRECT_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE SOLVER_LINEAR_DIRECT_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finalise a direct linear solver for a linear solver and deallocate all memory.
  SUBROUTINE SOLVER_LINEAR_DIRECT_FINALISE(LINEAR_DIRECT_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: LINEAR_DIRECT_SOLVER !<A pointer to the lienar direct solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_DIRECT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_DIRECT_SOLVER)) THEN
      DEALLOCATE(LINEAR_DIRECT_SOLVER)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_DIRECT_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_DIRECT_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_DIRECT_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_DIRECT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a direct linear solver for a lienar solver
  SUBROUTINE SOLVER_LINEAR_DIRECT_INITIALISE(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer the linear solver to initialise the direct solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_DIRECT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      IF(ASSOCIATED(LINEAR_SOLVER%DIRECT_SOLVER)) THEN
        CALL FLAG_ERROR("Direct solver is already associated for this linear solver.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(LINEAR_SOLVER%DIRECT_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear solver direct solver.",ERR,ERROR,*999)
        LINEAR_SOLVER%DIRECT_SOLVER%LINEAR_SOLVER=>LINEAR_SOLVER
        LINEAR_SOLVER%DIRECT_SOLVER%SOLVER_LIBRARY=SOLVER_CMISS_LIBRARY
        LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_LU
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_DIRECT_INITIALISE")
    RETURN
999 CALL SOLVER_LINEAR_DIRECT_FINALISE(LINEAR_SOLVER%DIRECT_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_LINEAR_DIRECT_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_DIRECT_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_DIRECT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Solve a linear direct solver 
  SUBROUTINE SOLVER_LINEAR_DIRECT_SOLVE(LINEAR_DIRECT_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: LINEAR_DIRECT_SOLVER !<A pointer to the linear direct solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RHS_VECTOR,SOLVER_VECTOR
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_DIRECT_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_DIRECT_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_DIRECT_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
            IF(ASSOCIATED(SOLVER_MATRICES)) THEN
              IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  RHS_VECTOR=>SOLVER_MATRICES%RHS_VECTOR
                  IF(ASSOCIATED(RHS_VECTOR)) THEN
                    SOLVER_VECTOR=>SOLVER_MATRICES%MATRICES(1)%PTR%SOLVER_VECTOR
                    IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                      SELECT CASE(LINEAR_DIRECT_SOLVER%SOLVER_LIBRARY)
                      CASE(SOLVER_CMISS_LIBRARY)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(SOLVER_PETSC_LIBRARY)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The solver library type of "// &
                          & TRIM(NUMBER_TO_VSTRING(LINEAR_DIRECT_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("RHS vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The number of solver matrices of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                  & " is invalid. There should only be one solver matrix for a linear direct solver."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver equations solver matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Linear solver solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Linear direct solver linear solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear direct solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_DIRECT_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_DIRECT_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_DIRECT_SOLVE")
    RETURN 1
    
  END SUBROUTINE SOLVER_LINEAR_DIRECT_SOLVE
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of direct linear solver
  SUBROUTINE SOLVER_LINEAR_DIRECT_TYPE_SET(SOLVER,DIRECT_SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set the direct linear solver type
    INTEGER(INTG), INTENT(IN) :: DIRECT_SOLVER_TYPE !<The type of direct linear solver to set \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_DIRECT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has alredy been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_DIRECT_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)) THEN
                IF(DIRECT_SOLVER_TYPE/=SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE) THEN
                  SELECT CASE(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%SOLVER_LIBRARY)
                  CASE(SOLVER_CMISS_LIBRARY)
                    SELECT CASE(DIRECT_SOLVER_TYPE)
                    CASE(SOLVER_DIRECT_LU)
                      SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_LU
                    CASE(SOLVER_DIRECT_CHOLESKY)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(SOLVER_DIRECT_SVD)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                    CASE DEFAULT
                      LOCAL_ERROR="The direct solver type of "//TRIM(NUMBER_TO_VSTRING(DIRECT_SOLVER_TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT                   
                  CASE(SOLVER_PETSC_LIBRARY)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver direct solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear direct solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_DIRECT_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_DIRECT_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_DIRECT_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_DIRECT_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Finalise a linear solver for a problem solver
  SUBROUTINE SOLVER_LINEAR_FINALISE(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer the linear solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      CALL SOLVER_LINEAR_DIRECT_FINALISE(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CALL SOLVER_LINEAR_ITERATIVE_FINALISE(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      DEALLOCATE(LINEAR_SOLVER)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a linear solver for a problem solver
  SUBROUTINE SOLVER_LINEAR_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the linear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_LINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
        CALL FLAG_ERROR("Linear solver is already associated for this solver.",ERR,ERROR,*998)
      ELSE
        !Allocate and initialise a linear solver
        ALLOCATE(SOLVER%LINEAR_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver linear solver.",ERR,ERROR,*999)
        SOLVER%LINEAR_SOLVER%SOLVER=>SOLVER
        NULLIFY(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)
        NULLIFY(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)
        !Default to a iterative solver
        SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE=SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
        CALL SOLVER_LINEAR_ITERATIVE_INITIALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_INITIALISE")
    RETURN
999 CALL SOLVER_LINEAR_FINALISE(SOLVER%LINEAR_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_LINEAR_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum absolute tolerance for an iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET(SOLVER,ABSOLUTE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set 
    REAL(DP), INTENT(IN) :: ABSOLUTE_TOLERANCE !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(ABSOLUTE_TOLERANCE>ZERO_TOLERANCE) THEN
                  SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE=ABSOLUTE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified absolute tolerance of "//TRIM(NUMBER_TO_VSTRING(ABSOLUTE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The absolute tolerance must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear iterative solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a linear iterative solver 
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_CREATE_FINISH(LINEAR_ITERATIVE_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: LINEAR_ITERATIVE_SOLVER !<A pointer to the linear iterative solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_ITERATIVE_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_ITERATIVE_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SELECT CASE(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY)
            CASE(SOLVER_CMISS_LIBRARY)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(SOLVER_PETSC_LIBRARY)
              !Create the solver matrices and vectors
              NULLIFY(SOLVER_MATRICES)
              CALL SOLVER_MATRICES_CREATE_START(SOLVER_EQUATIONS,SOLVER_MATRICES,ERR,ERROR,*999)
              CALL SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
              SELECT CASE(SOLVER_EQUATIONS%SPARSITY_TYPE)
              CASE(SOLVER_SPARSE_MATRICES)
                CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE(SOLVER_FULL_MATRICES)
                CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The specified solver equations sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL SOLVER_MATRICES_CREATE_FINISH(SOLVER_MATRICES,ERR,ERROR,*999)
              !Create the PETSc KSP solver
              CALL PETSC_KSPCREATE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,LINEAR_ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
              !Set the iterative solver type
              SELECT CASE(LINEAR_ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE)
              CASE(SOLVER_ITERATIVE_RICHARDSON)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPRICHARDSON,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_CHEBYCHEV)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPCHEBYCHEV,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPCG,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_BICONJUGATE_GRADIENT)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPBICG,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_GMRES)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPGMRES,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_BiCGSTAB)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPBCGS,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
                CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,PETSC_KSPCGS,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The iterative solver type of "// &
                  & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              !Get the pre-conditioner
              CALL PETSC_KSPGETPC(LINEAR_ITERATIVE_SOLVER%KSP,LINEAR_ITERATIVE_SOLVER%PC,ERR,ERROR,*999)
              !Set the pre-conditioner type
              SELECT CASE(LINEAR_ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE)
              CASE(SOLVER_ITERATIVE_NO_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCNONE,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_JACOBI_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCJACOBI,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCBJACOBI,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_SOR_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCSOR,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCICC,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCILU,ERR,ERROR,*999)
              CASE(SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER)
                CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PETSC_PCASM,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The iterative preconditioner type of "// &
                  & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              !Set the tolerances for the KSP solver
              CALL PETSC_KSPSETTOLERANCES(LINEAR_ITERATIVE_SOLVER%KSP,LINEAR_ITERATIVE_SOLVER%RELATIVE_TOLERANCE, &
                & LINEAR_ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE,LINEAR_ITERATIVE_SOLVER%DIVERGENCE_TOLERANCE, &
                & LINEAR_ITERATIVE_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
              !Set any further KSP options from the command line options
              CALL PETSC_KSPSETFROMOPTIONS(LINEAR_ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
              !Set the solver matrix to be the KSP matrix
              IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(1)%PTR%MATRIX
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  IF(ASSOCIATED(SOLVER_MATRIX%PETSC)) THEN
                    CALL PETSC_KSPSETOPERATORS(LINEAR_ITERATIVE_SOLVER%KSP,SOLVER_MATRIX%PETSC%MATRIX,SOLVER_MATRIX%PETSC%MATRIX, &
                      & PETSC_DIFFERENT_NONZERO_PATTERN,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Solver matrix PETSc is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrices distributed matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The given number of solver matrices of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                  & " is invalid. There should only be one solver matrix for a linear iterative solver."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The solver library type of "// &
                & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Linear solver solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Linear itreative solver linear solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear itreative solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum divergence tolerance for an iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET(SOLVER,DIVERGENCE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set 
    REAL(DP), INTENT(IN) :: DIVERGENCE_TOLERANCE !<The divergence tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(DIVERGENCE_TOLERANCE>ZERO_TOLERANCE) THEN
                  SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%DIVERGENCE_TOLERANCE=DIVERGENCE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified divergence tolerance of "// &
                    & TRIM(NUMBER_TO_VSTRING(DIVERGENCE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The divergence tolerance must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear iterative solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Finalise an iterative linear solver for a linear solver and deallocate all memory.
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_FINALISE(LINEAR_ITERATIVE_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: LINEAR_ITERATIVE_SOLVER !<A pointer the linear iterative solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_ITERATIVE_SOLVER)) THEN
      CALL PETSC_PCFINALISE(LINEAR_ITERATIVE_SOLVER%PC,ERR,ERROR,*999)
      CALL PETSC_KSPFINALISE(LINEAR_ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
      DEALLOCATE(LINEAR_ITERATIVE_SOLVER)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise an iterative linear solver for a linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_INITIALISE(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer the linear solver to initialise the iterative solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      IF(ASSOCIATED(LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
        CALL FLAG_ERROR("Iterative solver is already associated for this linear solver.",ERR,ERROR,*998)
      ELSE
        !Allocate and initialise a iterative solver
        ALLOCATE(LINEAR_SOLVER%ITERATIVE_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear solver iterative solver.",ERR,ERROR,*999)
        LINEAR_SOLVER%ITERATIVE_SOLVER%LINEAR_SOLVER=>LINEAR_SOLVER
        LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
        LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_GMRES
        LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE=SOLVER_ITERATIVE_JACOBI_PRECONDITIONER
        LINEAR_SOLVER%ITERATIVE_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS=100000
        LINEAR_SOLVER%ITERATIVE_SOLVER%RELATIVE_TOLERANCE=1.0E-05_DP
        LINEAR_SOLVER%ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE=1.0E-10_DP
        LINEAR_SOLVER%ITERATIVE_SOLVER%DIVERGENCE_TOLERANCE=1.0E5_DP
        CALL PETSC_PCINITIALISE(LINEAR_SOLVER%ITERATIVE_SOLVER%PC,ERR,ERROR,*999)
        CALL PETSC_KSPINITIALISE(LINEAR_SOLVER%ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_INITIALISE")
    RETURN
999 CALL SOLVER_LINEAR_ITERATIVE_FINALISE(LINEAR_SOLVER%ITERATIVE_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of iterations for an iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET(SOLVER,MAXIMUM_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set the maximum number of iterations
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_ITERATIONS !<The maximum number of iterations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(MAXIMUM_ITERATIONS>0) THEN
                  SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
                ELSE
                  LOCAL_ERROR="The specified maximum iterations of "//TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                    & " is invalid. The maximum number of iterations must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear iterative solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of preconditioner for an iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET(SOLVER,ITERATIVE_PRECONDITIONER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: ITERATIVE_PRECONDITIONER_TYPE !<The type of iterative preconditioner to set \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(ITERATIVE_PRECONDITIONER_TYPE/=SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE) THEN
                  !Intialise the new preconditioner type
                  SELECT CASE(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY)
                  CASE(SOLVER_PETSC_LIBRARY)
                    SELECT CASE(ITERATIVE_PRECONDITIONER_TYPE)
                    CASE(SOLVER_ITERATIVE_NO_PRECONDITIONER)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE=SOLVER_ITERATIVE_NO_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_JACOBI_PRECONDITIONER)
                      CALL FLAG_ERROR("Iterative Jacobi preconditioning is not implemented for a PETSc library.",ERR,ERROR,*999)
                    CASE(SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_SOR_PRECONDITIONER)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_SOR_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER 
                    CASE(SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER
                   CASE DEFAULT
                      LOCAL_ERROR="The iterative preconditioner type of "// &
                        & TRIM(NUMBER_TO_VSTRING(ITERATIVE_PRECONDITIONER_TYPE,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT                  
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear iterative solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the relative tolerance for an iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET(SOLVER,RELATIVE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set 
    REAL(DP), INTENT(IN) :: RELATIVE_TOLERANCE !<The relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(RELATIVE_TOLERANCE>ZERO_TOLERANCE) THEN
                  SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%RELATIVE_TOLERANCE=RELATIVE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified relative tolerance of "//TRIM(NUMBER_TO_VSTRING(RELATIVE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The relative tolerance must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear iterative solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Solves a linear iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_SOLVE(LINEAR_ITERATIVE_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: LINEAR_ITERATIVE_SOLVER !<A pointer the linear iterative solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: CONVERGED_REASON,NUMBER_ITERATIONS
    REAL(DP) :: RESIDUAL_NORM
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RHS_VECTOR,SOLVER_VECTOR
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_ITERATIVE_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_ITERATIVE_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
            IF(ASSOCIATED(SOLVER_MATRICES)) THEN
              IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
                RHS_VECTOR=>SOLVER_MATRICES%RHS_VECTOR
                IF(ASSOCIATED(RHS_VECTOR)) THEN
                  SOLVER_VECTOR=>SOLVER_MATRICES%MATRICES(1)%PTR%SOLVER_VECTOR
                  IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                    SELECT CASE(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY)
                    CASE(SOLVER_CMISS_LIBRARY)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(SOLVER_PETSC_LIBRARY)
                      IF(ASSOCIATED(RHS_VECTOR%PETSC)) THEN
                        IF(ASSOCIATED(SOLVER_VECTOR%PETSC)) THEN
                          !Solver the linear system
                          CALL PETSC_KSPSOLVE(LINEAR_ITERATIVE_SOLVER%KSP,RHS_VECTOR%PETSC%VECTOR,SOLVER_VECTOR%PETSC%VECTOR, &
                            & ERR,ERROR,*999)
                          !Check for convergence
                          CALL PETSC_KSPGETCONVERGEDREASON(LINEAR_ITERATIVE_SOLVER%KSP,CONVERGED_REASON,ERR,ERROR,*999)
                          SELECT CASE(CONVERGED_REASON)
                          CASE(PETSC_KSP_DIVERGED_NULL)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged null.",ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_ITS)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged its.",ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_DTOL)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged dtol.",ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_BREAKDOWN)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged breakdown.",ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_BREAKDOWN_BICG)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged breakdown BiCG.", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_NONSYMMETRIC)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged nonsymmetric.", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_INDEFINITE_PC)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged indefinite PC.", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_NAN)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged NaN.",ERR,ERROR,*999)
                          CASE(PETSC_KSP_DIVERGED_INDEFINITE_MAT)
                            CALL FLAG_WARNING("Linear iterative solver did not converge. PETSc diverged indefinite mat.", &
                              & ERR,ERROR,*999)
                          END SELECT
                          IF(SOLVER%OUTPUT_TYPE>=SOLVER_SOLVER_OUTPUT) THEN
                            !Output solution characteristics
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Linear iterative solver parameters:",ERR,ERROR,*999)
                            CALL PETSC_KSPGETITERATIONNUMBER(LINEAR_ITERATIVE_SOLVER%KSP,NUMBER_ITERATIONS,ERR,ERROR,*999)
                            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Final number of iterations = ",NUMBER_ITERATIONS, &
                              & ERR,ERROR,*999)
                            CALL PETSC_KSPGETRESIDUALNORM(LINEAR_ITERATIVE_SOLVER%KSP,RESIDUAL_NORM,ERR,ERROR,*999)
                            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Final residual norm = ",RESIDUAL_NORM, &
                              & ERR,ERROR,*999)
                            SELECT CASE(CONVERGED_REASON)
                            CASE(PETSC_KSP_CONVERGED_RTOL)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged RTol",ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_ATOL)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged ATol",ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_ITS)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged its",ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_CG_NEG_CURVE)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged CG neg curve", &
                                & ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_CG_CONSTRAINED)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged CG constrained", &
                                & ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_STEP_LENGTH)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged step length",ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_HAPPY_BREAKDOWN)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged happy breakdown", &
                                & ERR,ERROR,*999)
                            CASE(PETSC_KSP_CONVERGED_ITERATING)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged iterating",ERR,ERROR,*999)
                            END SELECT
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Solver vector PETSc vector is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("RHS vector petsc PETSc is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE DEFAULT
                      LOCAL_ERROR="The solver library type of "// &
                        & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("RHS vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The given number of solver matrices of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                  & " is invalid. There should only be one solver matrix for a linear iterative solver."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Linear solver solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Linear itreative solver linear solver is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Linear iterative solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_SOLVE
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_TYPE_SET(SOLVER,ITERATIVE_SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: ITERATIVE_SOLVER_TYPE !<The type of iterative linear solver to set \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
              IF(ASSOCIATED(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(ITERATIVE_SOLVER_TYPE/=SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE) THEN
                  !Intialise the new solver type
                  SELECT CASE(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY)
                  CASE(SOLVER_PETSC_LIBRARY)
                    SELECT CASE(ITERATIVE_SOLVER_TYPE)
                    CASE(SOLVER_ITERATIVE_RICHARDSON)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_RICHARDSON
                    CASE(SOLVER_ITERATIVE_CHEBYCHEV)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_CHEBYCHEV
                    CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_CONJUGATE_GRADIENT
                    CASE(SOLVER_ITERATIVE_BICONJUGATE_GRADIENT)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_BICONJUGATE_GRADIENT
                    CASE(SOLVER_ITERATIVE_GMRES)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_BiCGSTAB
                    CASE(SOLVER_ITERATIVE_BiCGSTAB)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_BiCGSTAB
                    CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
                      SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_CONJGRAD_SQUARED
                    CASE DEFAULT
                      LOCAL_ERROR="The iterative solver type of "//TRIM(NUMBER_TO_VSTRING(ITERATIVE_SOLVER_TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT                  
                ENDIF
              ELSE
                CALL FLAG_ERROR("The solver linear solver iterative solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The solver is not a linear iterative solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Solve a linear solver 
  SUBROUTINE SOLVER_LINEAR_SOLVE(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer to the linear solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVE_TYPE)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        CALL SOLVER_LINEAR_DIRECT_SOLVE(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        CALL SOLVER_LINEAR_ITERATIVE_SOLVE(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVE_TYPE,"*",ERR,ERROR))// &
          & " is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_SOLVE
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of linear solver
  SUBROUTINE SOLVER_LINEAR_TYPE_SET(SOLVER,LINEAR_SOLVE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the linear solver type
    INTEGER(INTG), INTENT(IN) :: LINEAR_SOLVE_TYPE !<The type of linear solver to set \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*998)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(LINEAR_SOLVE_TYPE/=SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE) THEN
              !Intialise the new solver type
              SELECT CASE(LINEAR_SOLVE_TYPE)
              CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
                CALL SOLVER_LINEAR_DIRECT_INITIALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
                CALL SOLVER_LINEAR_ITERATIVE_INITIALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              !Finalise the old solver type
              SELECT CASE(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE)
              CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
                CALL SOLVER_LINEAR_DIRECT_FINALISE(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
                CALL SOLVER_LINEAR_ITERATIVE_FINALISE(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The linear solver type of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              SOLVER%LINEAR_SOLVER%LINEAR_SOLVE_TYPE=LINEAR_SOLVE_TYPE
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_TYPE_SET")
    RETURN
999 SELECT CASE(LINEAR_SOLVE_TYPE)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
      CALL SOLVER_LINEAR_DIRECT_FINALISE(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
      CALL SOLVER_LINEAR_ITERATIVE_FINALISE(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    END SELECT
998 CALL ERRORS("SOLVER_LINEAR_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Assembles the solver matrices and rhs from the dynamic equations.
  SUBROUTINE SOLVER_MATRICES_DYNAMIC_ASSEMBLE(SOLVER,SELECTION_TYPE,ERR,ERROR,*)

    !Argument variableg
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The type of matrix selection \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,equations_matrix_idx,equations_matrix_number, &
      & equations_row_number,equations_row_number2,equations_set_idx,EQUATIONS_STORAGE_TYPE,field_dof,jacobian_column_idx, &
      & jacobian_column_number,JACOBIAN_STORAGE_TYPE,jacobian_row_number,rhs_boundary_condition, &
      & residual_variable_dof,residual_variable_type,rhs_field_dof,rhs_global_dof,rhs_variable_dof,rhs_variable_type, &
      & variable_boundary_condition,solver_column_idx,solver_column_number,solver_matrix_idx,solver_row_idx, &
      & solver_row_number,variable_dof,variable_field_dof,variable_global_dof,variable_idx,variable_type
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(SP) :: SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),USER_ELAPSED,USER_TIME1(1),USER_TIME2(1)
    REAL(DP) :: column_coupling_coefficient,MATRIX_VALUE,RESIDUAL_VALUE,RHS_VALUE,row_coupling_coefficient,SOURCE_VALUE,VALUE
    REAL(DP), POINTER :: DEPENDENT_PARAMETERS(:),EQUATIONS_MATRIX_DATA(:),JACOBIAN_MATRIX_DATA(:),SOURCE_PARAMETERS(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX,JACOBIAN_DISTRIBUTED_MATRIX, &
      & PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: EQUATIONS_RESIDUAL_VECTOR,EQUATIONS_RHS_VECTOR,EQUATIONS_SOURCE_VECTOR, &
      & SOLVER_RESIDUAL_VECTOR,SOLVER_RHS_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING,RESIDUAL_DOMAIN_MAPPING,RHS_DOMAIN_MAPPING,VARIABLE_DOMAIN_MAPPING
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,SOURCE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,RESIDUAL_VARIABLE,RHS_VARIABLE
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            !Assemble the solver matrices
            NULLIFY(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
              !Assemble solver matrices
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  IF(SOLVER_MATRIX%UPDATE_MATRIX) THEN              
                    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                    IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN                
                      !Initialise matrix to zero
                      CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(SOLVER_DISTRIBUTED_MATRIX,0.0_DP,ERR,ERROR,*999)
                      !Loop over the equations sets
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        !First Loop over the linear equations matrices
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%PTR
                          IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                            EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                              LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
                              IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                                EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                  IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE, &
                                      & ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                    SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                      !Loop over the rows of the equations matrix
                                      DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this equations row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the equations matrix
                                          DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                            !Loop over the solution columns this equations column is mapped to
                                            DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                              & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                                & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                                & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=EQUATIONS_MATRIX_DATA(equations_row_number+(equations_column_number-1)* &
                                                & EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)
                                            ENDDO !solver_column_idx
                                          ENDDO !equations_column_number
                                        ENDDO !solver_row_idx
                                      ENDDO !equations_row_number
                                    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX,ROW_INDICES, &
                                        & COLUMN_INDICES,ERR,ERROR,*999)
                                      !Loop over the rows of the equations matrix
                                      DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this equations row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the equations matrix
                                          DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                            & ROW_INDICES(equations_row_number+1)-1
                                            equations_column_number=COLUMN_INDICES(equations_column_idx)
                                            !Loop over the solution columns this equations column is mapped to
                                            DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                              & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                                & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                                & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=EQUATIONS_MATRIX_DATA(equations_column_idx)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)                                    
                                            ENDDO !solution_column_idx
                                          ENDDO !equations_column_idx
                                        ENDDO !solution_row_idx
                                      ENDDO !equations_row_number
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE DEFAULT
                                      LOCAL_ERROR="The matrix storage type of "// &
                                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                  ELSE
                                    CALL FLAG_ERROR("The equations matrix distributed matrix is not associated",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Linear matrices equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations matrix linear matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("The equations matrix is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("The equations matrix equations to solver map is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !equations_matrix_idx
                        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                          & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
                          !Now set the values from the equations Jacobian
                          JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                          IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                            JACOBIAN_MATRIX=>JACOBIAN_TO_SOLVER_MAP%JACOBIAN_MATRIX
                            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                              NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
                              IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                                EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  JACOBIAN_DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
                                  IF(ASSOCIATED(JACOBIAN_DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_STORAGE_TYPE, &
                                      & ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_DATA_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                    SELECT CASE(JACOBIAN_STORAGE_TYPE)
                                    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                      !Loop over the rows of the Jacobian matrix
                                      DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this Jacobian row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the Jacobian matrix
                                          DO jacobian_column_number=1,JACOBIAN_MATRIX%NUMBER_OF_COLUMNS
                                            !Loop over the solution columns this Jacobian column is mapped to
                                            DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                              & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                                & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                                & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=JACOBIAN_MATRIX_DATA(jacobian_row_number+(jacobian_column_number-1)* &
                                                & EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)
                                            ENDDO !solver_column_idx
                                          ENDDO !jacobian_column_number
                                        ENDDO !solver_row_idx
                                      ENDDO !jacobian_row_number
                                    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(JACOBIAN_DISTRIBUTED_MATRIX,ROW_INDICES, &
                                        & COLUMN_INDICES,ERR,ERROR,*999)
                                      !Loop over the rows of the Jacobian matrix
                                      DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this Jacobian row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the Jacobian matrix
                                          DO jacobian_column_idx=ROW_INDICES(jacobian_row_number), &
                                            & ROW_INDICES(jacobian_row_number+1)-1
                                            jacobian_column_number=COLUMN_INDICES(jacobian_column_idx)
                                            !Loop over the solution columns this equations column is mapped to
                                            DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                              & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                                & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                                & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=JACOBIAN_MATRIX_DATA(jacobian_column_idx)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)                                    
                                            ENDDO !solution_column_idx
                                          ENDDO !jacobian_column_idx
                                        ENDDO !solution_row_idx
                                      ENDDO !jacobian_row_number
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE DEFAULT
                                      LOCAL_ERROR="The matrix storage type of "// &
                                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                  ELSE
                                    CALL FLAG_ERROR("The Jacobian matrix distributed matrix is not associated",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Jacobian matrix nonlinear matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO !equations_set_idx
                      !Update the solver matrix values
                      CALL DISTRIBUTED_MATRIX_UPDATE_START(SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                      IF(ASSOCIATED(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)) THEN
                        CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                      ENDIF
                      PREVIOUS_SOLVER_DISTRIBUTED_MATRIX=>SOLVER_DISTRIBUTED_MATRIX
                    ELSE
                      CALL FLAG_ERROR("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF !Update matrix
                ELSE
                  CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !solver_matrix_idx
              IF(ASSOCIATED(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)) THEN
                CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solver matrices assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solver matrices assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ENDIF
            NULLIFY(SOLVER_RHS_VECTOR)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
              !Assemble rhs vector
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER_MATRICES%UPDATE_RHS_VECTOR) THEN
                SOLVER_RHS_VECTOR=>SOLVER_MATRICES%RHS_VECTOR
                IF(ASSOCIATED(SOLVER_RHS_VECTOR)) THEN
                  !Initialise the RHS to zero
                  CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(SOLVER_RHS_VECTOR,0.0_DP,ERR,ERROR,*999)            
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
                        IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
                          CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                          EQUATIONS=>EQUATIONS_SET%EQUATIONS
                          IF(ASSOCIATED(EQUATIONS)) THEN
                            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                            IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                              RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                              IF(ASSOCIATED(RHS_MAPPING)) THEN
                                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                                SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                                EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
                                  IF(ASSOCIATED(RHS_VECTOR)) THEN
                                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                      LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
                                      IF(.NOT.ASSOCIATED(LINEAR_MATRICES)) &
                                        & CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                    IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                      SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
                                      IF(ASSOCIATED(SOURCE_VECTOR)) THEN
                                        EQUATIONS_SOURCE_VECTOR=>SOURCE_VECTOR%VECTOR
                                        IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
                                          SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                                          IF(ASSOCIATED(SOURCE_FIELD)) THEN
                                            CALL FIELD_PARAMETER_SET_GET(SOURCE_FIELD,FIELD_VALUES_SET_TYPE,SOURCE_PARAMETERS, &
                                              & ERR,ERROR,*999)                                     
                                          ELSE
                                            CALL FLAG_ERROR("Source field is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ELSE
                                          CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Equations matrices source vector is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ENDIF
                                    FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                                    IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
!!TODO: what if the equations set doesn't have a RHS vector???
                                      rhs_variable_type=RHS_MAPPING%RHS_VARIABLE_TYPE
                                      RHS_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
                                      RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
                                      EQUATIONS_RHS_VECTOR=>RHS_VECTOR%VECTOR
                                      !Loop over the rows in the equations set
                                      DO equations_row_number=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
                                        !Add in equations RHS values
                                        CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_RHS_VECTOR,equations_row_number,RHS_VALUE, &
                                          & ERR,ERROR,*999)
                                        !Loop over the solver rows associated with this equations set row
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%COUPLING_COEFFICIENTS( &
                                            & solver_row_idx)
                                          VALUE=RHS_VALUE*row_coupling_coefficient
                                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                            & ERR,ERROR,*999)
                                        ENDDO !solver_row_idx                          
                                        rhs_variable_dof=RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(equations_row_number)
                                        rhs_field_dof=RHS_VARIABLE%DOF_LIST(rhs_variable_dof)
                                        rhs_global_dof=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_field_dof)
                                        rhs_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(rhs_global_dof)
                                        !Apply boundary conditions
                                        SELECT CASE(rhs_boundary_condition)
                                        CASE(EQUATIONS_SET_NOT_FIXED)
                                          !Set Direchlet boundary conditions
                                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                            !Loop over the dependent variables associated with this equations set row
                                            DO variable_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                                              variable_type=LINEAR_MAPPING%MATRIX_VARIABLE_TYPES(variable_idx)
                                              DEPENDENT_VARIABLE=>LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                                & VARIABLE
                                              VARIABLE_DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                                              variable_dof=LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(equations_row_number, &
                                                & variable_idx)
                                              variable_field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                              variable_global_dof=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(variable_field_dof)
                                              variable_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS( &
                                                & variable_global_dof)
                                              IF(variable_boundary_condition==EQUATIONS_SET_FIXED_BOUNDARY_CONDITION) THEN
                                                DO equations_matrix_idx=1,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                                                  & variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                                  equations_matrix_number=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                                                    & variable_type)%EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                                  EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(equations_matrix_number)%PTR
                                                  equations_column_number=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                                                    & variable_type)%DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF( &
                                                    & variable_dof)
                                                  DO equations_row_number2=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
                                                    CALL DISTRIBUTED_MATRIX_VALUES_GET(EQUATIONS_MATRIX%MATRIX, &
                                                      & equations_row_number2,equations_column_number,MATRIX_VALUE,ERR,ERROR,*999)
                                                    DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                      & equations_row_number2)%NUMBER_OF_SOLVER_ROWS
                                                      solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                        & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                        & equations_row_number2)%SOLVER_ROWS(solver_row_idx)
                                                      row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                        & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                        & equations_row_number2)%COUPLING_COEFFICIENTS(solver_row_idx)
                                                      VALUE=-1.0_DP*MATRIX_VALUE*DEPENDENT_PARAMETERS(variable_field_dof)* &
                                                        & row_coupling_coefficient
                                                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number, &
                                                        & VALUE,ERR,ERROR,*999)
                                                    ENDDO !solver_row_idx
                                                  ENDDO !equations_row_number2
                                                ENDDO !matrix_idx
                                              ENDIF
                                            ENDDO !variable_idx
                                          ENDIF
                                        CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                                          !Set Neumann boundary conditions
                                          !Loop over the solver rows associated with this equations set row
                                          DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                            solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                            row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                              & COUPLING_COEFFICIENTS(solver_row_idx)
                                            VALUE=DEPENDENT_PARAMETERS(rhs_field_dof)*row_coupling_coefficient
                                            CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                              & ERR,ERROR,*999)
                                          ENDDO !solver_row_idx
                                        CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                                          !Set Robin or is it Cauchy??? boundary conditions
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global boundary condition of "// &
                                            & TRIM(NUMBER_TO_VSTRING(rhs_boundary_condition,"*",ERR,ERROR))// &
                                            & " for RHS field dof number "//TRIM(NUMBER_TO_VSTRING(rhs_field_dof,"*",ERR,ERROR))// &
                                            & " is invalid."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                          !Add in equations source values
                                          CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_SOURCE_VECTOR,equations_row_number, &
                                            & SOURCE_VALUE,ERR,ERROR,*999)
                                          !Loop over the solver rows associated with this equations set row
                                          DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                            solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                            row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                              & COUPLING_COEFFICIENTS(solver_row_idx)
                                            VALUE=SOURCE_VALUE*row_coupling_coefficient
                                            CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                              & ERR,ERROR,*999)
                                          ENDDO !solver_row_idx                          
                                        ENDIF
                                      ENDDO !equations_row_number
                                    ELSE
                                      CALL FLAG_ERROR("Fixed conditions is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                    IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                      CALL FIELD_PARAMETER_SET_RESTORE(SOURCE_FIELD,FIELD_VALUES_SET_TYPE,SOURCE_PARAMETERS, &
                                        & ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Equations matrice RHS vector is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Equations equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                          ENDIF
                          CALL FIELD_PARAMETER_SET_RESTORE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS, &
                            & ERR,ERROR,*999)
                        ELSE
                          CALL FLAG_ERROR("Dependent field domain mapping is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                  !Start the update the solver RHS vector values
                  CALL DISTRIBUTED_VECTOR_UPDATE_START(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The solver RHS vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solver RHS assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solver RHS assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
              !Assemble residual vector
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVER_RESIDUAL_VECTOR)
              IF(SOLVER_MATRICES%UPDATE_RESIDUAL) THEN            
                SOLVER_RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                  !Initialise the residual to zero              
                  CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(SOLVER_RESIDUAL_VECTOR,0.0_DP,ERR,ERROR,*999)       
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                              EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
                              IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
                                IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                    LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
                                    IF(.NOT.ASSOCIATED(LINEAR_MATRICES)) &
                                      & CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                  FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                                  IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
                                    residual_variable_type=NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                                    RESIDUAL_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLE
                                    RESIDUAL_DOMAIN_MAPPING=>RESIDUAL_VARIABLE%DOMAIN_MAPPING
                                    EQUATIONS_RESIDUAL_VECTOR=>NONLINEAR_MATRICES%RESIDUAL
                                    !Loop over the rows in the equations set
                                    DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                      residual_variable_dof=NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP( &
                                        & equations_row_number)
                                      CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_RESIDUAL_VECTOR,equations_row_number, &
                                        & RESIDUAL_VALUE,ERR,ERROR,*999)
                                      !Loop over the solver rows associated with this equations set residual row
                                      DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                        solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                        row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%COUPLING_COEFFICIENTS( &
                                          & solver_row_idx)
                                        VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                        !Add in nonlinear residual values                                    
                                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                          & ERR,ERROR,*999)
                                      ENDDO !solver_row_idx                          
                                    ENDDO !equations_row_number
                                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                      !Calculate the linear part of the residual
                                      DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                        EQUATIONS_MATRIX=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)% &
                                          & EQUATIONS_MATRIX                                    
                                        IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                                          DEPENDENT_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                            & equations_matrix_idx)%VARIABLE
                                          IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                                            EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                              & EQUATIONS_STORAGE_TYPE,ERR,ERROR,*999)
                                            CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                              & ERR,ERROR,*999)
                                            SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                              !Loop over the rows of the equations matrix
                                              DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                                !Loop over the solution rows this equations row is mapped to
                                                DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                                  solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & SOLVER_ROWS(solver_row_idx)
                                                  row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                    & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & COUPLING_COEFFICIENTS(solver_row_idx)
                                                  VALUE=0.0_DP
                                                  !Loop over the columns of the equations matrix
                                                  DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                                    variable_dof=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                      & equations_matrix_idx)%COLUMN_TO_DOF_MAP(equations_column_number)
                                                    field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                                    VALUE=VALUE+EQUATIONS_MATRIX_DATA(equations_row_number+ &
                                                      & (equations_column_number-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                                                      & row_coupling_coefficient*DEPENDENT_PARAMETERS(field_dof)
                                                  ENDDO !equations_column_number
                                                  VALUE=VALUE*LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                    & equations_matrix_idx)%MATRIX_COEFFICIENT
                                                  !Add in nonlinear residual values                                    
                                                  CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RESIDUAL_VECTOR,solver_row_number, &
                                                    & VALUE,ERR,ERROR,*999)
                                                ENDDO !solver_row_idx
                                              ENDDO !equations_row_number
                                            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                                & ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*999)
                                              !Loop over the rows of the equations matrix
                                              DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                                !Loop over the solution rows this equations row is mapped to
                                                DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                                  solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & SOLVER_ROWS(solver_row_idx)
                                                  row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                    & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & COUPLING_COEFFICIENTS(solver_row_idx)
                                                  VALUE=0.0_DP
                                                  !Loop over the columns of the equations matrix
                                                  DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                                    & ROW_INDICES(equations_row_number+1)-1
                                                    equations_column_number=COLUMN_INDICES(equations_column_idx)
                                                    variable_dof=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                      & equations_matrix_idx)%COLUMN_TO_DOF_MAP(equations_column_number)
                                                    field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                                    !Add in nonlinear residual values
                                                    VALUE=VALUE+EQUATIONS_MATRIX_DATA(equations_column_idx)* &
                                                      & row_coupling_coefficient*DEPENDENT_PARAMETERS(field_dof)
                                                  ENDDO !equations_column_idx
                                                  VALUE=VALUE*LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                    & equations_matrix_idx)%MATRIX_COEFFICIENT
                                                  CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RESIDUAL_VECTOR,solver_row_number, &
                                                    & VALUE,ERR,ERROR,*999)
                                                ENDDO !solution_row_idx
                                              ENDDO !equations_row_number                                          
                                            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                            CASE DEFAULT
                                              LOCAL_ERROR="The matrix storage type of "// &
                                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                            END SELECT
                                            CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX, &
                                              & EQUATIONS_MATRIX_DATA,ERR,ERROR,*999)
                                          ELSE
                                            CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ELSE
                                          CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ENDDO !equations_matrix_idx
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Fixed conditions is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Equations matrics residual vector is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations equations matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                        ENDIF
                        CALL FIELD_PARAMETER_SET_RESTORE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                  !Start the update the solver residual vector values
                  CALL DISTRIBUTED_VECTOR_UPDATE_START(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The solver residual vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solver residual assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solver residual assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(ASSOCIATED(SOLVER_RHS_VECTOR)) THEN
              CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver matrices solver mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_DYNAMIC_ASSEMBLE

  !
  !================================================================================================================================
  !

  !>Assembles the solver matrices and rhs from the static equations.
  SUBROUTINE SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SELECTION_TYPE,ERR,ERROR,*)

    !Argument variableg
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The type of matrix selection \see SOLVER_ROUTINES_SelectMatricesTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,equations_matrix_idx,equations_matrix_number, &
      & equations_row_number,equations_row_number2,equations_set_idx,EQUATIONS_STORAGE_TYPE,field_dof,jacobian_column_idx, &
      & jacobian_column_number,JACOBIAN_STORAGE_TYPE,jacobian_row_number,rhs_boundary_condition, &
      & residual_variable_dof,residual_variable_type,rhs_field_dof,rhs_global_dof,rhs_variable_dof,rhs_variable_type, &
      & variable_boundary_condition,solver_column_idx,solver_column_number,solver_matrix_idx,solver_row_idx, &
      & solver_row_number,variable_dof,variable_field_dof,variable_global_dof,variable_idx,variable_type
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(SP) :: SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),USER_ELAPSED,USER_TIME1(1),USER_TIME2(1)
    REAL(DP) :: column_coupling_coefficient,MATRIX_VALUE,RESIDUAL_VALUE,RHS_VALUE,row_coupling_coefficient,SOURCE_VALUE,VALUE
    REAL(DP), POINTER :: DEPENDENT_PARAMETERS(:),EQUATIONS_MATRIX_DATA(:),JACOBIAN_MATRIX_DATA(:),SOURCE_PARAMETERS(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX,JACOBIAN_DISTRIBUTED_MATRIX, &
      & PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: EQUATIONS_RESIDUAL_VECTOR,EQUATIONS_RHS_VECTOR,EQUATIONS_SOURCE_VECTOR, &
      & SOLVER_RESIDUAL_VECTOR,SOLVER_RHS_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING,RESIDUAL_DOMAIN_MAPPING,RHS_DOMAIN_MAPPING,VARIABLE_DOMAIN_MAPPING
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,SOURCE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,RESIDUAL_VARIABLE,RHS_VARIABLE
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("SOLVER_MATRICES_STATIC_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            !Assemble the solver matrices
            NULLIFY(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
              !Assemble solver matrices
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  IF(SOLVER_MATRIX%UPDATE_MATRIX) THEN              
                    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                    IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN                
                      !Initialise matrix to zero
                      CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(SOLVER_DISTRIBUTED_MATRIX,0.0_DP,ERR,ERROR,*999)
                      !Loop over the equations sets
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        !First Loop over the linear equations matrices
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%PTR
                          IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                            EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                            IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                              LINEAR_MATRICES=>EQUATIONS_MATRIX%LINEAR_MATRICES
                              IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                                EQUATIONS_MATRICES=>LINEAR_MATRICES%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                  IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE, &
                                      & ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                    SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                      !Loop over the rows of the equations matrix
                                      DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this equations row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the equations matrix
                                          DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                            !Loop over the solution columns this equations column is mapped to
                                            DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                              & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                                & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                                & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=EQUATIONS_MATRIX_DATA(equations_row_number+(equations_column_number-1)* &
                                                & EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)
                                            ENDDO !solver_column_idx
                                          ENDDO !equations_column_number
                                        ENDDO !solver_row_idx
                                      ENDDO !equations_row_number
                                    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX,ROW_INDICES, &
                                        & COLUMN_INDICES,ERR,ERROR,*999)
                                      !Loop over the rows of the equations matrix
                                      DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this equations row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the equations matrix
                                          DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                            & ROW_INDICES(equations_row_number+1)-1
                                            equations_column_number=COLUMN_INDICES(equations_column_idx)
                                            !Loop over the solution columns this equations column is mapped to
                                            DO solver_column_idx=1,EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                              & equations_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=EQUATIONS_TO_SOLVER_MAP%EQUATIONS_COL_SOLVER_COLS_MAP( &
                                                & equations_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=EQUATIONS_TO_SOLVER_MAP% &
                                                & EQUATIONS_COL_SOLVER_COLS_MAP(equations_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=EQUATIONS_MATRIX_DATA(equations_column_idx)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)                                    
                                            ENDDO !solution_column_idx
                                          ENDDO !equations_column_idx
                                        ENDDO !solution_row_idx
                                      ENDDO !equations_row_number
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE DEFAULT
                                      LOCAL_ERROR="The matrix storage type of "// &
                                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                  ELSE
                                    CALL FLAG_ERROR("The equations matrix distributed matrix is not associated",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Linear matrices equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations matrix linear matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("The equations matrix is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("The equations matrix equations to solver map is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !equations_matrix_idx
                        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                          & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
                          !Now set the values from the equations Jacobian
                          JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAP
                          IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                            JACOBIAN_MATRIX=>JACOBIAN_TO_SOLVER_MAP%JACOBIAN_MATRIX
                            IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                              NONLINEAR_MATRICES=>JACOBIAN_MATRIX%NONLINEAR_MATRICES
                              IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                                EQUATIONS_MATRICES=>NONLINEAR_MATRICES%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  JACOBIAN_DISTRIBUTED_MATRIX=>JACOBIAN_MATRIX%JACOBIAN
                                  IF(ASSOCIATED(JACOBIAN_DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_STORAGE_TYPE, &
                                      & ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_DATA_GET(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                    SELECT CASE(JACOBIAN_STORAGE_TYPE)
                                    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                      !Loop over the rows of the Jacobian matrix
                                      DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this Jacobian row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the Jacobian matrix
                                          DO jacobian_column_number=1,JACOBIAN_MATRIX%NUMBER_OF_COLUMNS
                                            !Loop over the solution columns this Jacobian column is mapped to
                                            DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                              & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                                & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                                & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=JACOBIAN_MATRIX_DATA(jacobian_row_number+(jacobian_column_number-1)* &
                                                & EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)
                                            ENDDO !solver_column_idx
                                          ENDDO !jacobian_column_number
                                        ENDDO !solver_row_idx
                                      ENDDO !jacobian_row_number
                                    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(JACOBIAN_DISTRIBUTED_MATRIX,ROW_INDICES, &
                                        & COLUMN_INDICES,ERR,ERROR,*999)
                                      !Loop over the rows of the Jacobian matrix
                                      DO jacobian_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                        !Loop over the solution rows this Jacobian row is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(jacobian_row_number)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          !Loop over the columns of the Jacobian matrix
                                          DO jacobian_column_idx=ROW_INDICES(jacobian_row_number), &
                                            & ROW_INDICES(jacobian_row_number+1)-1
                                            jacobian_column_number=COLUMN_INDICES(jacobian_column_idx)
                                            !Loop over the solution columns this equations column is mapped to
                                            DO solver_column_idx=1,JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                              & jacobian_column_number)%NUMBER_OF_SOLVER_COLS
                                              solver_column_number=JACOBIAN_TO_SOLVER_MAP%JACOBIAN_COL_SOLVER_COLS_MAP( &
                                                & jacobian_column_number)%SOLVER_COLS(solver_column_idx)
                                              column_coupling_coefficient=JACOBIAN_TO_SOLVER_MAP% &
                                                & JACOBIAN_COL_SOLVER_COLS_MAP(jacobian_column_number)% &
                                                & COUPLING_COEFFICIENTS(solver_column_idx)
                                              !Add in the solver matrix value
                                              VALUE=JACOBIAN_MATRIX_DATA(jacobian_column_idx)*row_coupling_coefficient* &
                                                & column_coupling_coefficient
                                              CALL DISTRIBUTED_MATRIX_VALUES_ADD(SOLVER_DISTRIBUTED_MATRIX,solver_row_number, &
                                                & solver_column_number,VALUE,ERR,ERROR,*999)                                    
                                            ENDDO !solution_column_idx
                                          ENDDO !jacobian_column_idx
                                        ENDDO !solution_row_idx
                                      ENDDO !jacobian_row_number
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                    CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE DEFAULT
                                      LOCAL_ERROR="The matrix storage type of "// &
                                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                    CALL DISTRIBUTED_MATRIX_DATA_RESTORE(JACOBIAN_DISTRIBUTED_MATRIX,JACOBIAN_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                  ELSE
                                    CALL FLAG_ERROR("The Jacobian matrix distributed matrix is not associated",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Nonlinear matrices equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Jacobian matrix nonlinear matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Jacobian matrix is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDDO !equations_set_idx
                      !Update the solver matrix values
                      CALL DISTRIBUTED_MATRIX_UPDATE_START(SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                      IF(ASSOCIATED(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)) THEN
                        CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                      ENDIF
                      PREVIOUS_SOLVER_DISTRIBUTED_MATRIX=>SOLVER_DISTRIBUTED_MATRIX
                    ELSE
                      CALL FLAG_ERROR("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF !Update matrix
                ELSE
                  CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !solver_matrix_idx
              IF(ASSOCIATED(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)) THEN
                CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solver matrices assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solver matrices assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ENDIF
            NULLIFY(SOLVER_RHS_VECTOR)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
              !Assemble rhs vector
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER_MATRICES%UPDATE_RHS_VECTOR) THEN
                SOLVER_RHS_VECTOR=>SOLVER_MATRICES%RHS_VECTOR
                IF(ASSOCIATED(SOLVER_RHS_VECTOR)) THEN
                  !Initialise the RHS to zero
                  CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(SOLVER_RHS_VECTOR,0.0_DP,ERR,ERROR,*999)            
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
                        IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
                          CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                          EQUATIONS=>EQUATIONS_SET%EQUATIONS
                          IF(ASSOCIATED(EQUATIONS)) THEN
                            EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                            IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                              RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                              IF(ASSOCIATED(RHS_MAPPING)) THEN
                                LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                                SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                                EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
                                IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                  RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
                                  IF(ASSOCIATED(RHS_VECTOR)) THEN
                                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                      LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
                                      IF(.NOT.ASSOCIATED(LINEAR_MATRICES)) &
                                        & CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                    IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                      SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
                                      IF(ASSOCIATED(SOURCE_VECTOR)) THEN
                                        EQUATIONS_SOURCE_VECTOR=>SOURCE_VECTOR%VECTOR
                                        IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
                                          SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
                                          IF(ASSOCIATED(SOURCE_FIELD)) THEN
                                            CALL FIELD_PARAMETER_SET_GET(SOURCE_FIELD,FIELD_VALUES_SET_TYPE,SOURCE_PARAMETERS, &
                                              & ERR,ERROR,*999)                                     
                                          ELSE
                                            CALL FLAG_ERROR("Source field is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ELSE
                                          CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Equations matrices source vector is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ENDIF
                                    FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                                    IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
!!TODO: what if the equations set doesn't have a RHS vector???
                                      rhs_variable_type=RHS_MAPPING%RHS_VARIABLE_TYPE
                                      RHS_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
                                      RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
                                      EQUATIONS_RHS_VECTOR=>RHS_VECTOR%VECTOR
                                      !Loop over the rows in the equations set
                                      DO equations_row_number=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
                                        !Add in equations RHS values
                                        CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_RHS_VECTOR,equations_row_number,RHS_VALUE, &
                                          & ERR,ERROR,*999)
                                        !Loop over the solver rows associated with this equations set row
                                        DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%COUPLING_COEFFICIENTS( &
                                            & solver_row_idx)
                                          VALUE=RHS_VALUE*row_coupling_coefficient
                                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                            & ERR,ERROR,*999)
                                        ENDDO !solver_row_idx                          
                                        rhs_variable_dof=RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(equations_row_number)
                                        rhs_field_dof=RHS_VARIABLE%DOF_LIST(rhs_variable_dof)
                                        rhs_global_dof=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_field_dof)
                                        rhs_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(rhs_global_dof)
                                        !Apply boundary conditions
                                        SELECT CASE(rhs_boundary_condition)
                                        CASE(EQUATIONS_SET_NOT_FIXED)
                                          !Set Direchlet boundary conditions
                                          IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                            !Loop over the dependent variables associated with this equations set row
                                            DO variable_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_MATRIX_VARIABLES
                                              variable_type=LINEAR_MAPPING%MATRIX_VARIABLE_TYPES(variable_idx)
                                              DEPENDENT_VARIABLE=>LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                                & VARIABLE
                                              VARIABLE_DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                                              variable_dof=LINEAR_MAPPING%EQUATIONS_ROW_TO_VARIABLE_DOF_MAPS(equations_row_number, &
                                                & variable_idx)
                                              variable_field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                              variable_global_dof=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(variable_field_dof)
                                              variable_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS( &
                                                & variable_global_dof)
                                              IF(variable_boundary_condition==EQUATIONS_SET_FIXED_BOUNDARY_CONDITION) THEN
                                                DO equations_matrix_idx=1,LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                                                  & variable_type)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                                  equations_matrix_number=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                                                    & variable_type)%EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                                  EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(equations_matrix_number)%PTR
                                                  equations_column_number=LINEAR_MAPPING%VAR_TO_EQUATIONS_MATRICES_MAPS( &
                                                    & variable_type)%DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF( &
                                                    & variable_dof)
                                                  DO equations_row_number2=1,EQUATIONS_MAPPING%TOTAL_NUMBER_OF_ROWS
                                                    CALL DISTRIBUTED_MATRIX_VALUES_GET(EQUATIONS_MATRIX%MATRIX, &
                                                      & equations_row_number2,equations_column_number,MATRIX_VALUE,ERR,ERROR,*999)
                                                    DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                      & equations_row_number2)%NUMBER_OF_SOLVER_ROWS
                                                      solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                        & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                        & equations_row_number2)%SOLVER_ROWS(solver_row_idx)
                                                      row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                        & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                        & equations_row_number2)%COUPLING_COEFFICIENTS(solver_row_idx)
                                                      VALUE=-1.0_DP*MATRIX_VALUE*DEPENDENT_PARAMETERS(variable_field_dof)* &
                                                        & row_coupling_coefficient
                                                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number, &
                                                        & VALUE,ERR,ERROR,*999)
                                                    ENDDO !solver_row_idx
                                                  ENDDO !equations_row_number2
                                                ENDDO !matrix_idx
                                              ENDIF
                                            ENDDO !variable_idx
                                          ENDIF
                                        CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                                          !Set Neumann boundary conditions
                                          !Loop over the solver rows associated with this equations set row
                                          DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                            solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                            row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                              & COUPLING_COEFFICIENTS(solver_row_idx)
                                            VALUE=DEPENDENT_PARAMETERS(rhs_field_dof)*row_coupling_coefficient
                                            CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                              & ERR,ERROR,*999)
                                          ENDDO !solver_row_idx
                                        CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                                          !Set Robin or is it Cauchy??? boundary conditions
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global boundary condition of "// &
                                            & TRIM(NUMBER_TO_VSTRING(rhs_boundary_condition,"*",ERR,ERROR))// &
                                            & " for RHS field dof number "//TRIM(NUMBER_TO_VSTRING(rhs_field_dof,"*",ERR,ERROR))// &
                                            & " is invalid."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                          !Add in equations source values
                                          CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_SOURCE_VECTOR,equations_row_number, &
                                            & SOURCE_VALUE,ERR,ERROR,*999)
                                          !Loop over the solver rows associated with this equations set row
                                          DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                            solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                            row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                              & COUPLING_COEFFICIENTS(solver_row_idx)
                                            VALUE=SOURCE_VALUE*row_coupling_coefficient
                                            CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                              & ERR,ERROR,*999)
                                          ENDDO !solver_row_idx                          
                                        ENDIF
                                      ENDDO !equations_row_number
                                    ELSE
                                      CALL FLAG_ERROR("Fixed conditions is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                    IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                      CALL FIELD_PARAMETER_SET_RESTORE(SOURCE_FIELD,FIELD_VALUES_SET_TYPE,SOURCE_PARAMETERS, &
                                        & ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Equations matrice RHS vector is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Equations equations matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                          ENDIF
                          CALL FIELD_PARAMETER_SET_RESTORE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS, &
                            & ERR,ERROR,*999)
                        ELSE
                          CALL FLAG_ERROR("Dependent field domain mapping is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                  !Start the update the solver RHS vector values
                  CALL DISTRIBUTED_VECTOR_UPDATE_START(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The solver RHS vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solver RHS assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solver RHS assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
              !Assemble residual vector
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              NULLIFY(SOLVER_RESIDUAL_VECTOR)
              IF(SOLVER_MATRICES%UPDATE_RESIDUAL) THEN            
                SOLVER_RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                  !Initialise the residual to zero              
                  CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(SOLVER_RESIDUAL_VECTOR,0.0_DP,ERR,ERROR,*999)       
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                          IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                            NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
                            IF(ASSOCIATED(NONLINEAR_MAPPING)) THEN
                              LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                              EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
                              IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                                NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
                                IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
                                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                    LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
                                    IF(.NOT.ASSOCIATED(LINEAR_MATRICES)) &
                                      & CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                  FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                                  IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
                                    residual_variable_type=NONLINEAR_MAPPING%RESIDUAL_VARIABLE_TYPE
                                    RESIDUAL_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLE
                                    RESIDUAL_DOMAIN_MAPPING=>RESIDUAL_VARIABLE%DOMAIN_MAPPING
                                    EQUATIONS_RESIDUAL_VECTOR=>NONLINEAR_MATRICES%RESIDUAL
                                    !Loop over the rows in the equations set
                                    DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                      residual_variable_dof=NONLINEAR_MAPPING%EQUATIONS_ROW_TO_RESIDUAL_DOF_MAP( &
                                        & equations_row_number)
                                      CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_RESIDUAL_VECTOR,equations_row_number, &
                                        & RESIDUAL_VALUE,ERR,ERROR,*999)
                                      !Loop over the solver rows associated with this equations set residual row
                                      DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                        solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                        row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%COUPLING_COEFFICIENTS( &
                                          & solver_row_idx)
                                        VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                        !Add in nonlinear residual values                                    
                                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                          & ERR,ERROR,*999)
                                      ENDDO !solver_row_idx                          
                                    ENDDO !equations_row_number
                                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                                      !Calculate the linear part of the residual
                                      DO equations_matrix_idx=1,LINEAR_MAPPING%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                                        EQUATIONS_MATRIX=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)% &
                                          & EQUATIONS_MATRIX                                    
                                        IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                                          DEPENDENT_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                            & equations_matrix_idx)%VARIABLE
                                          IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                                            EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                            CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                              & EQUATIONS_STORAGE_TYPE,ERR,ERROR,*999)
                                            CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                              & ERR,ERROR,*999)
                                            SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                              !Loop over the rows of the equations matrix
                                              DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                                !Loop over the solution rows this equations row is mapped to
                                                DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                                  solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & SOLVER_ROWS(solver_row_idx)
                                                  row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                    & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & COUPLING_COEFFICIENTS(solver_row_idx)
                                                  VALUE=0.0_DP
                                                  !Loop over the columns of the equations matrix
                                                  DO equations_column_number=1,EQUATIONS_MATRIX%NUMBER_OF_COLUMNS
                                                    variable_dof=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                      & equations_matrix_idx)%COLUMN_TO_DOF_MAP(equations_column_number)
                                                    field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                                    VALUE=VALUE+EQUATIONS_MATRIX_DATA(equations_row_number+ &
                                                      & (equations_column_number-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)* &
                                                      & row_coupling_coefficient*DEPENDENT_PARAMETERS(field_dof)
                                                  ENDDO !equations_column_number
                                                  VALUE=VALUE*LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                    & equations_matrix_idx)%MATRIX_COEFFICIENT
                                                  !Add in nonlinear residual values                                    
                                                  CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RESIDUAL_VECTOR,solver_row_number, &
                                                    & VALUE,ERR,ERROR,*999)
                                                ENDDO !solver_row_idx
                                              ENDDO !equations_row_number
                                            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                              CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                                & ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*999)
                                              !Loop over the rows of the equations matrix
                                              DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                                !Loop over the solution rows this equations row is mapped to
                                                DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                                  solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & SOLVER_ROWS(solver_row_idx)
                                                  row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                    & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                    & COUPLING_COEFFICIENTS(solver_row_idx)
                                                  VALUE=0.0_DP
                                                  !Loop over the columns of the equations matrix
                                                  DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                                    & ROW_INDICES(equations_row_number+1)-1
                                                    equations_column_number=COLUMN_INDICES(equations_column_idx)
                                                    variable_dof=LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                      & equations_matrix_idx)%COLUMN_TO_DOF_MAP(equations_column_number)
                                                    field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                                    !Add in nonlinear residual values
                                                    VALUE=VALUE+EQUATIONS_MATRIX_DATA(equations_column_idx)* &
                                                      & row_coupling_coefficient*DEPENDENT_PARAMETERS(field_dof)
                                                  ENDDO !equations_column_idx
                                                  VALUE=VALUE*LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS( &
                                                    & equations_matrix_idx)%MATRIX_COEFFICIENT
                                                  CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RESIDUAL_VECTOR,solver_row_number, &
                                                    & VALUE,ERR,ERROR,*999)
                                                ENDDO !solution_row_idx
                                              ENDDO !equations_row_number                                          
                                            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                            CASE DEFAULT
                                              LOCAL_ERROR="The matrix storage type of "// &
                                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                            END SELECT
                                            CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX, &
                                              & EQUATIONS_MATRIX_DATA,ERR,ERROR,*999)
                                          ELSE
                                            CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ELSE
                                          CALL FLAG_ERROR("Equations matrix is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ENDDO !equations_matrix_idx
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Fixed conditions is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Equations matrics residual vector is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Equations equations matrices is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Equations equations mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
                        ENDIF
                        CALL FIELD_PARAMETER_SET_RESTORE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                  !Start the update the solver residual vector values
                  CALL DISTRIBUTED_VECTOR_UPDATE_START(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The solver residual vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solver residual assembly = ",USER_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solver residual assembly = ",SYSTEM_ELAPSED, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(ASSOCIATED(SOLVER_RHS_VECTOR)) THEN
              CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver matrices solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_STATIC_ASSEMBLE")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_STATIC_ASSEMBLE",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_STATIC_ASSEMBLE")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_STATIC_ASSEMBLE

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum absolute tolerance for a nonlinear Newton solver
  SUBROUTINE SOLVER_NEWTON_ABSOLUTE_TOLERANCE_SET(SOLVER,ABSOLUTE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the absolute tolerance for
    REAL(DP), INTENT(IN) :: ABSOLUTE_TOLERANCE !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_ABSOLUTE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(ABSOLUTE_TOLERANCE>ZERO_TOLERANCE) THEN
                  NEWTON_SOLVER%ABSOLUTE_TOLERANCE=ABSOLUTE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified absolute tolerance of "//TRIM(NUMBER_TO_VSTRING(ABSOLUTE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The absolute tolerance must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_ABSOLUTE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_ABSOLUTE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_ABSOLUTE_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_ABSOLUTE_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a Newton solver 
  SUBROUTINE SOLVER_NEWTON_CREATE_FINISH(NEWTON_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer to the Newton solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NEWTON_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(NEWTON_SOLVER)) THEN
      SELECT CASE(NEWTON_SOLVER%NEWTON_SOLVE_TYPE)
      CASE(SOLVER_NEWTON_LINESEARCH)
        CALL SOLVER_NEWTON_LINESEARCH_CREATE_FINISH(NEWTON_SOLVER%LINESEARCH_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_NEWTON_TRUSTREGION)
        CALL SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH(NEWTON_SOLVER%TRUSTREGION_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The Newton solver type of "// &
          & TRIM(NUMBER_TO_VSTRING(NEWTON_SOLVER%NEWTON_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Newton solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finalise a Newton solver and deallocate all memory
  SUBROUTINE SOLVER_NEWTON_FINALISE(NEWTON_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer the Newton solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_NEWTON_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NEWTON_SOLVER)) THEN
      CALL SOLVER_NEWTON_LINESEARCH_FINALISE(NEWTON_SOLVER%LINESEARCH_SOLVER,ERR,ERROR,*999)
      CALL SOLVER_NEWTON_TRUSTREGION_FINALISE(NEWTON_SOLVER%TRUSTREGION_SOLVER,ERR,ERROR,*999)
      DEALLOCATE(NEWTON_SOLVER)
    ENDIF
         
    CALL EXITS("SOLVER_NEWTON_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a Newton solver for a nonlinear solver
  SUBROUTINE SOLVER_NEWTON_INITIALISE(NONLINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer the solver to initialise the Newton solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("SOLVER_NEWTON_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      IF(ASSOCIATED(NONLINEAR_SOLVER%NEWTON_SOLVER)) THEN
        CALL FLAG_ERROR("Newton solver is already associated for this nonlinear solver.",ERR,ERROR,*998)
      ELSE
        !Allocate and initialise a Newton solver
        ALLOCATE(NONLINEAR_SOLVER%NEWTON_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nonlinear solver Newton solver.",ERR,ERROR,*999)
        NONLINEAR_SOLVER%NEWTON_SOLVER%NONLINEAR_SOLVER=>NONLINEAR_SOLVER
        NONLINEAR_SOLVER%NEWTON_SOLVER%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS=0
        NONLINEAR_SOLVER%NEWTON_SOLVER%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS=0
        NONLINEAR_SOLVER%NEWTON_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS=50
        NONLINEAR_SOLVER%NEWTON_SOLVER%MAXIMUM_NUMBER_OF_FUNCTION_EVALUATIONS=1000
        NONLINEAR_SOLVER%NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE=SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
        NONLINEAR_SOLVER%NEWTON_SOLVER%ABSOLUTE_TOLERANCE=1.0E-10_DP
        NONLINEAR_SOLVER%NEWTON_SOLVER%RELATIVE_TOLERANCE=1.0E-05_DP
        NONLINEAR_SOLVER%NEWTON_SOLVER%SOLUTION_TOLERANCE=1.0E-05_DP
        NULLIFY(NONLINEAR_SOLVER%NEWTON_SOLVER%LINESEARCH_SOLVER)
        NULLIFY(NONLINEAR_SOLVER%NEWTON_SOLVER%TRUSTREGION_SOLVER)
        !Default to a Newton linesearch solver
        NONLINEAR_SOLVER%NEWTON_SOLVER%NEWTON_SOLVE_TYPE=SOLVER_NEWTON_LINESEARCH
        CALL SOLVER_NEWTON_LINESEARCH_INITIALISE(NONLINEAR_SOLVER%NEWTON_SOLVER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Nonlinear solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_INITIALISE")
    RETURN
999 CALL SOLVER_NEWTON_FINALISE(NONLINEAR_SOLVER%NEWTON_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_NEWTON_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of Jacobian calculation type for a Newton solver
  SUBROUTINE SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET(SOLVER,JACOBIAN_CALCULATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the Jacobian calculation type
    INTEGER(INTG), INTENT(IN) :: JACOBIAN_CALCULATION_TYPE !<The type of Jacobian calculation type to set \see SOLVER_ROUTINES_JacobianCalculationTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(JACOBIAN_CALCULATION_TYPE/=NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE) THEN
                  SELECT CASE(JACOBIAN_CALCULATION_TYPE)
                  CASE(SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED)
                    NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE=SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED
                  CASE(SOLVER_NEWTON_JACOBIAN_ANALTYIC_CALCULATED)
                    NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE=SOLVER_NEWTON_JACOBIAN_ANALTYIC_CALCULATED
                  CASE(SOLVER_NEWTON_JACOBIAN_FD_CALCULATED)
                    NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE=SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
                  CASE DEFAULT
                    LOCAL_ERROR="The Jacobian calculation type of "// &
                      & TRIM(NUMBER_TO_VSTRING(JACOBIAN_CALCULATION_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT              
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The Solver nonlinear solver is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_JACOBIAN_CALCULATION_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search alpha for a Newton linesearch solver
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_ALPHA_SET(SOLVER,LINESEARCH_ALPHA,ERR,ERROR,*)
    
    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the line search alpha for
    REAL(DP), INTENT(IN) :: LINESEARCH_ALPHA !<The line search alpha to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_ALPHA_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NEWTON_SOLVER%NEWTON_SOLVE_TYPE==SOLVER_NEWTON_LINESEARCH) THEN
                  LINESEARCH_SOLVER=>NEWTON_SOLVER%LINESEARCH_SOLVER
                  IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
                    IF(LINESEARCH_ALPHA>ZERO_TOLERANCE) THEN
                      LINESEARCH_SOLVER%LINESEARCH_ALPHA=LINESEARCH_ALPHA
                    ELSE
                      LOCAL_ERROR="The specified line search alpha of "//TRIM(NUMBER_TO_VSTRING(LINESEARCH_ALPHA,"*",ERR,ERROR))// &
                        & " is invalid. The line search alpha must be > 0."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The Newton solver line search solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The Newton solver is not a line search solver.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_ALPHA_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_ALPHA_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_ALPHA_SET")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_ALPHA_SET
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nonlinear Newton line search solver
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_CREATE_FINISH(LINESEARCH_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER !<A pointer the nonlinear Newton line search solver to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    EXTERNAL :: SNESDefaultComputeJacobianColor
    EXTERNAL :: PROBLEM_SOLVER_JACOBIAN_EVALUATE_PETSC
    EXTERNAL :: PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC
    EXTERNAL :: SOLVER_NONLINEAR_MONITOR_PETSC
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL_VECTOR
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_JACOBIAN
    TYPE(VARYING_STRING) :: LOCAL_ERROR
  
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
      NEWTON_SOLVER=>LINESEARCH_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN
        NONLINEAR_SOLVER=>NEWTON_SOLVER%NONLINEAR_SOLVER
        IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
          SOLVER=>NONLINEAR_SOLVER%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SELECT CASE(LINESEARCH_SOLVER%SOLVER_LIBRARY)
              CASE(SOLVER_CMISS_LIBRARY)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(SOLVER_PETSC_LIBRARY)
                !Create the solver matrices and vectors
                NULLIFY(SOLVER_MATRICES)
                CALL SOLVER_MATRICES_CREATE_START(SOLVER_EQUATIONS,SOLVER_MATRICES,ERR,ERROR,*999)
                CALL SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
                SELECT CASE(SOLVER_EQUATIONS%SPARSITY_TYPE)
                CASE(SOLVER_SPARSE_MATRICES)
                  CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                CASE(SOLVER_FULL_MATRICES)
                  CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified solver equations sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(SOLVER_EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                CALL SOLVER_MATRICES_CREATE_FINISH(SOLVER_MATRICES,ERR,ERROR,*999)
                !Create the PETSc SNES solver
                CALL PETSC_SNESCREATE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,LINESEARCH_SOLVER%SNES,ERR,ERROR,*999)
                !Set the nonlinear solver type to be a Newton line search solver
                CALL PETSC_SNESSETTYPE(LINESEARCH_SOLVER%SNES,PETSC_SNESLS,ERR,ERROR,*999)
                !Set the nonlinear function
                RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                IF(ASSOCIATED(RESIDUAL_VECTOR)) THEN
                  IF(ASSOCIATED(RESIDUAL_VECTOR%PETSC)) THEN                   
                    CALL PETSC_SNESSETFUNCTION(LINESEARCH_SOLVER%SNES,RESIDUAL_VECTOR%PETSC%VECTOR, &
                      & PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC,SOLVER,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("The residual vector PETSc is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrices residual vector is not associated.",ERR,ERROR,*999)
                ENDIF
                !Set the Jacobian
                IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
                  SOLVER_JACOBIAN=>SOLVER_MATRICES%MATRICES(1)%PTR
                  IF(ASSOCIATED(SOLVER_JACOBIAN)) THEN
                    JACOBIAN_MATRIX=>SOLVER_JACOBIAN%MATRIX
                    IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                      IF(ASSOCIATED(JACOBIAN_MATRIX%PETSC)) THEN
                        SELECT CASE(NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE)
                        CASE(SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED)
                          CALL FLAG_ERROR("Cannot have no Jacobian calculation for a PETSc nonlinear linesearch solver.", &
                            & ERR,ERROR,*999)
                        CASE(SOLVER_NEWTON_JACOBIAN_ANALTYIC_CALCULATED)
                          SOLVER_JACOBIAN%UPDATE_MATRIX=.TRUE. !CMISS will fill in the Jacobian values
                          CALL PETSC_SNESSETJACOBIAN(LINESEARCH_SOLVER%SNES,JACOBIAN_MATRIX%PETSC%MATRIX, &
                            & JACOBIAN_MATRIX%PETSC%MATRIX,PROBLEM_SOLVER_JACOBIAN_EVALUATE_PETSC,SOLVER,ERR,ERROR,*999)
                        CASE(SOLVER_NEWTON_JACOBIAN_FD_CALCULATED)
                          SOLVER_JACOBIAN%UPDATE_MATRIX=.FALSE. !Petsc will fill in the Jacobian values
                          CALL DISTRIBUTED_MATRIX_FORM(JACOBIAN_MATRIX,ERR,ERROR,*999)
                          CALL DISTRIBUTED_MATRIX_OUTPUT(GENERAL_OUTPUT_TYPE,JACOBIAN_MATRIX,ERR,ERROR,*999)
                          CALL PETSC_MATGETCOLORING(JACOBIAN_MATRIX%PETSC%MATRIX,PETSC_MATCOLORING_SL,LINESEARCH_SOLVER% &
                            & JACOBIAN_ISCOLORING,ERR,ERROR,*999)
                          CALL PETSC_MATFDCOLORINGCREATE(JACOBIAN_MATRIX%PETSC%MATRIX,LINESEARCH_SOLVER% &
                            & JACOBIAN_ISCOLORING,LINESEARCH_SOLVER%JACOBIAN_FDCOLORING,ERR,ERROR,*999)
                          CALL PETSC_ISCOLORINGDESTROY(LINESEARCH_SOLVER%JACOBIAN_ISCOLORING,ERR,ERROR,*999)
                          CALL PETSC_MATFDCOLORINGSETFUNCTIONSNES(LINESEARCH_SOLVER%JACOBIAN_FDCOLORING, &
                            & PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC,SOLVER,ERR,ERROR,*999)
                          CALL PETSC_MATFDCOLORINGSETFROMOPTIONS(LINESEARCH_SOLVER%JACOBIAN_FDCOLORING,ERR,ERROR,*999)
                          CALL PETSC_SNESSETJACOBIAN(LINESEARCH_SOLVER%SNES,JACOBIAN_MATRIX%PETSC%MATRIX, &
                            & JACOBIAN_MATRIX%PETSC%MATRIX,SNESDefaultComputeJacobianColor,LINESEARCH_SOLVER% &
                            & JACOBIAN_FDCOLORING,ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The Jacobian calculation type of "// &
                            & TRIM(NUMBER_TO_VSTRING(NEWTON_SOLVER%JACOBIAN_CALCULATION_TYPE,"*",ERR,ERROR))// &
                            & " is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        CALL FLAG_ERROR("Jacobian matrix PETSc is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Solver Jacobian matrix is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The solver Jacobian is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Invalid number of solver matrices. The number of solver matrices is "// &
                    & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))//" and it should be 1."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                  !Set the monitor
                  CALL PETSC_SNESMONITORSET(LINESEARCH_SOLVER%SNES,SOLVER_NONLINEAR_MONITOR_PETSC,SOLVER,ERR,ERROR,*999)
                ENDIF
                !Set the line search type
                SELECT CASE(LINESEARCH_SOLVER%LINESEARCH_TYPE)
                CASE(SOLVER_NEWTON_LINESEARCH_NONORMS)
                  CALL PETSC_SNESLINESEARCHSET(LINESEARCH_SOLVER%SNES,PETSC_SNES_LINESEARCH_NONORMS,ERR,ERROR,*999)
                CASE(SOLVER_NEWTON_LINESEARCH_NONE)
                  CALL PETSC_SNESLINESEARCHSET(LINESEARCH_SOLVER%SNES,PETSC_SNES_LINESEARCH_NO,ERR,ERROR,*999)
                CASE(SOLVER_NEWTON_LINESEARCH_QUADRATIC)
                  CALL PETSC_SNESLINESEARCHSET(LINESEARCH_SOLVER%SNES,PETSC_SNES_LINESEARCH_QUADRATIC,ERR,ERROR,*999)
                CASE(SOLVER_NEWTON_LINESEARCH_CUBIC)
                  CALL PETSC_SNESLINESEARCHSET(LINESEARCH_SOLVER%SNES,PETSC_SNES_LINESEARCH_CUBIC,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The nonlinear Newton line search type of "// &
                    & TRIM(NUMBER_TO_VSTRING(LINESEARCH_SOLVER%LINESEARCH_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                !Set line search parameters
                CALL PETSC_SNESLINESEARCHSETPARAMS(LINESEARCH_SOLVER%SNES,LINESEARCH_SOLVER%LINESEARCH_ALPHA, &
                  & LINESEARCH_SOLVER%LINESEARCH_MAXSTEP,LINESEARCH_SOLVER%LINESEARCH_STEPTOLERANCE, &
                  & ERR,ERROR,*999)
                !Set the tolerances for the SNES solver
                CALL PETSC_SNESSETTOLERANCES(LINESEARCH_SOLVER%SNES,NEWTON_SOLVER%ABSOLUTE_TOLERANCE, &
                  & NEWTON_SOLVER%RELATIVE_TOLERANCE,NEWTON_SOLVER%SOLUTION_TOLERANCE, &
                  & NEWTON_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS, &
                  & NEWTON_SOLVER%MAXIMUM_NUMBER_OF_FUNCTION_EVALUATIONS,ERR,ERROR,*999)            
                !Set any further SNES options from the command line options
                CALL PETSC_SNESSETFROMOPTIONS(LINESEARCH_SOLVER%SNES,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solver library type of "// &
                  & TRIM(NUMBER_TO_VSTRING(LINESEARCH_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Nonlinear solver solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Newton solver nonlinear solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Linesearch solver Newton solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Line search solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Finalise a nonlinear Newton line search solver and deallocate all memory
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_FINALISE(LINESEARCH_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER !<A pointer the nonlinear Newton line search solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
      CALL PETSC_ISCOLORINGFINALISE(LINESEARCH_SOLVER%JACOBIAN_ISCOLORING,ERR,ERROR,*999)
      CALL PETSC_MATFDCOLORINGFINALISE(LINESEARCH_SOLVER%JACOBIAN_FDCOLORING,ERR,ERROR,*999)
      CALL PETSC_SNESFINALISE(LINESEARCH_SOLVER%SNES,ERR,ERROR,*999)
      DEALLOCATE(LINESEARCH_SOLVER)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_FINALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a nonlinear Newton line search solver for a Newton solver
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_INITIALISE(NEWTON_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer the nonlinear Newton solver to initialise the Newton line search solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
  
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(NEWTON_SOLVER)) THEN
      IF(ASSOCIATED(NEWTON_SOLVER%LINESEARCH_SOLVER)) THEN
        CALL FLAG_ERROR("Netwon line search solver is already associated for this Newton solver.",ERR,ERROR,*998)
      ELSE
        !Allocate and initialise the Newton linesearch solver
        ALLOCATE(NEWTON_SOLVER%LINESEARCH_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nonlinear solver Newton line search solver.",ERR,ERROR,*999)
        NEWTON_SOLVER%LINESEARCH_SOLVER%NEWTON_SOLVER=>NEWTON_SOLVER
        NEWTON_SOLVER%LINESEARCH_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
        NEWTON_SOLVER%LINESEARCH_SOLVER%LINESEARCH_TYPE=SOLVER_NEWTON_LINESEARCH_CUBIC
        NEWTON_SOLVER%LINESEARCH_SOLVER%LINESEARCH_ALPHA=PETSC_DEFAULT_DOUBLE_PRECISION
        NEWTON_SOLVER%LINESEARCH_SOLVER%LINESEARCH_MAXSTEP=PETSC_DEFAULT_DOUBLE_PRECISION
        NEWTON_SOLVER%LINESEARCH_SOLVER%LINESEARCH_STEPTOLERANCE=PETSC_DEFAULT_DOUBLE_PRECISION
        CALL PETSC_ISCOLORINGINITIALISE(NEWTON_SOLVER%LINESEARCH_SOLVER%JACOBIAN_ISCOLORING,ERR,ERROR,*999)
        CALL PETSC_MATFDCOLORINGINITIALISE(NEWTON_SOLVER%LINESEARCH_SOLVER%JACOBIAN_FDCOLORING,ERR,ERROR,*999)
        CALL PETSC_SNESINITIALISE(NEWTON_SOLVER%LINESEARCH_SOLVER%SNES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Newton solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_INITIALISE")
    RETURN
999 CALL SOLVER_NEWTON_LINESEARCH_FINALISE(NEWTON_SOLVER%LINESEARCH_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the line search maximum step for a nonlinear Newton linesearch solver
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_MAXSTEP_SET(SOLVER,LINESEARCH_MAXSTEP,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the line search maximum step for
    REAL(DP), INTENT(IN) :: LINESEARCH_MAXSTEP !<The line search maximum step to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_MAXSTEP_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NEWTON_SOLVER%NEWTON_SOLVE_TYPE==SOLVER_NEWTON_LINESEARCH) THEN
                  LINESEARCH_SOLVER=>NEWTON_SOLVER%LINESEARCH_SOLVER
                  IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
                    IF(LINESEARCH_MAXSTEP>ZERO_TOLERANCE) THEN
                      LINESEARCH_SOLVER%LINESEARCH_MAXSTEP=LINESEARCH_MAXSTEP
                    ELSE
                      LOCAL_ERROR="The specified line search maximum step of "// &
                        & TRIM(NUMBER_TO_VSTRING(LINESEARCH_MAXSTEP,"*",ERR,ERROR))// &
                        & " is invalid. The line search maximum step must be > 0."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The Newton solver line search solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The Newton solver is not a line search solver.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_MAXSTEP_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_MAXSTEP_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_MAXSTEP_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_MAXSTEP_SET
        
  !
  !================================================================================================================================
  !

  !Solves a nonlinear Newton line search solver 
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_SOLVE(LINESEARCH_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER !<A pointer to the nonlinear Newton line search solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: CONVERGED_REASON,NUMBER_ITERATIONS
    REAL(DP) :: FUNCTION_NORM
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RHS_VECTOR,SOLVER_VECTOR
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
      NEWTON_SOLVER=>LINESEARCH_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN
        NONLINEAR_SOLVER=>NEWTON_SOLVER%NONLINEAR_SOLVER
        IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
          SOLVER=>NONLINEAR_SOLVER%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                IF(SOLVER_MATRICES%NUMBER_OF_MATRICES==1) THEN
                  RHS_VECTOR=>SOLVER_MATRICES%RHS_VECTOR
                  IF(ASSOCIATED(RHS_VECTOR)) THEN
                    SOLVER_VECTOR=>SOLVER_MATRICES%MATRICES(1)%PTR%SOLVER_VECTOR
                    IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                      SELECT CASE(LINESEARCH_SOLVER%SOLVER_LIBRARY)
                      CASE(SOLVER_CMISS_LIBRARY)
                        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                      CASE(SOLVER_PETSC_LIBRARY)
                        !Solve the nonlinear equations
                        CALL PETSC_SNESSOLVE(LINESEARCH_SOLVER%SNES,RHS_VECTOR%PETSC%VECTOR,SOLVER_VECTOR%PETSC%VECTOR, &
                          & ERR,ERROR,*999)
                        !Check for convergence
                        CALL PETSC_SNESGETCONVERGEDREASON(LINESEARCH_SOLVER%SNES,CONVERGED_REASON,ERR,ERROR,*999)
                        SELECT CASE(CONVERGED_REASON)
                        CASE(PETSC_SNES_DIVERGED_FUNCTION_COUNT)
                          CALL FLAG_WARNING("Nonlinear line search solver did not converge. PETSc diverged function count.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_SNES_DIVERGED_LINEAR_SOLVE)
                          CALL FLAG_WARNING("Nonlinear line search solver did not converge. PETSc diverged linear solve.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_SNES_DIVERGED_FNORM_NAN)
                          CALL FLAG_WARNING("Nonlinear line search solver did not converge. PETSc diverged F Norm NaN.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_SNES_DIVERGED_MAX_IT)
                          CALL FLAG_WARNING("Nonlinear line search solver did not converge. PETSc diverged maximum iterations.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_SNES_DIVERGED_LS_FAILURE)
                          CALL FLAG_WARNING("Nonlinear line search solver did not converge. PETSc diverged line search failure.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_SNES_DIVERGED_LOCAL_MIN)
                          CALL FLAG_WARNING("Nonlinear line search solver did not converge. PETSc diverged local minimum.", &
                            & ERR,ERROR,*999)
                        END SELECT
                        IF(SOLVER%OUTPUT_TYPE>=SOLVER_SOLVER_OUTPUT) THEN
                          !Output solution characteristics
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Newton linesearch solver parameters:",ERR,ERROR,*999)
                          CALL PETSC_SNESGETITERATIONNUMBER(LINESEARCH_SOLVER%SNES,NUMBER_ITERATIONS,ERR,ERROR,*999)
                          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Final number of iterations = ",NUMBER_ITERATIONS, &
                            & ERR,ERROR,*999)
                          CALL PETSC_SNESGETFUNCTIONNORM(LINESEARCH_SOLVER%SNES,FUNCTION_NORM,ERR,ERROR,*999)
                          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Final function norm = ",FUNCTION_NORM, &
                            & ERR,ERROR,*999)
                          SELECT CASE(CONVERGED_REASON)
                          CASE(PETSC_SNES_CONVERGED_FNORM_ABS)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged F Norm absolute", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_SNES_CONVERGED_FNORM_RELATIVE)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged F Norm relative", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_SNES_CONVERGED_PNORM_RELATIVE)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged P Norm relative", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_SNES_CONVERGED_ITS)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged its",ERR,ERROR,*999)
                          CASE(PETSC_SNES_CONVERGED_ITERATING)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged iterating",ERR,ERROR,*999)
                          END SELECT
                        ENDIF
                      CASE DEFAULT
                        LOCAL_ERROR="The Newton line search solver library type of "// &
                          & TRIM(NUMBER_TO_VSTRING(LINESEARCH_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Solver RHS vector is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The number of solver matrices of "// &
                    & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                    & " is invalid. There should only be one solver matrix for a Newton linesearch solver."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Nonlinear solver solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Newton solver nonlinear solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Linesearch solver Newton solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linesearch solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_SOLVE
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search step tolerance for a nonlinear Newton linesearch solver
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_STEPTOL_SET(SOLVER,LINESEARCH_STEPTOL,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the line search step tolerance for
    REAL(DP), INTENT(IN) :: LINESEARCH_STEPTOL !<The line search step tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_STEPTOL_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NEWTON_SOLVER%NEWTON_SOLVE_TYPE==SOLVER_NEWTON_LINESEARCH) THEN
                  LINESEARCH_SOLVER=>NEWTON_SOLVER%LINESEARCH_SOLVER
                  IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
                    IF(LINESEARCH_STEPTOL>ZERO_TOLERANCE) THEN
                      LINESEARCH_SOLVER%LINESEARCH_STEPTOLERANCE=LINESEARCH_STEPTOL
                    ELSE
                      LOCAL_ERROR="The specified line search step tolerance of "// &
                        & TRIM(NUMBER_TO_VSTRING(LINESEARCH_STEPTOL,"*",ERR,ERROR))// &
                        & " is invalid. The line search step tolerance must be > 0."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The Newton solver line search solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The Newton solver is not a line search solver.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_STEPTOL_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_STEPTOL_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_STEPTOL_SET")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_STEPTOL_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search type for a nonlinear Newton linesearch solver
  SUBROUTINE SOLVER_NEWTON_LINESEARCH_TYPE_SET(SOLVER,LINESEARCH_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the line search type for
    INTEGER(INTG), INTENT(IN) :: LINESEARCH_TYPE !<The line search type to set \see SOLVER_ROUTINES_NewtonLineSearchTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: LINESEARCH_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_LINESEARCH_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NEWTON_LINESEARCH) THEN
                  LINESEARCH_SOLVER=>NEWTON_SOLVER%LINESEARCH_SOLVER
                  IF(ASSOCIATED(LINESEARCH_SOLVER)) THEN
                    SELECT CASE(LINESEARCH_TYPE)
                    CASE(SOLVER_NEWTON_LINESEARCH_NONORMS)
                      LINESEARCH_SOLVER%LINESEARCH_TYPE=SOLVER_NEWTON_LINESEARCH_NONORMS
                    CASE(SOLVER_NEWTON_LINESEARCH_NONE)
                      LINESEARCH_SOLVER%LINESEARCH_TYPE=SOLVER_NEWTON_LINESEARCH_NONE
                    CASE(SOLVER_NEWTON_LINESEARCH_QUADRATIC)
                      LINESEARCH_SOLVER%LINESEARCH_TYPE=SOLVER_NEWTON_LINESEARCH_QUADRATIC
                    CASE(SOLVER_NEWTON_LINESEARCH_CUBIC)
                      LINESEARCH_SOLVER%LINESEARCH_TYPE=SOLVER_NEWTON_LINESEARCH_CUBIC
                    CASE DEFAULT
                      LOCAL_ERROR="The specified line search type of "//TRIM(NUMBER_TO_VSTRING(LINESEARCH_TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    CALL FLAG_ERROR("The Newton solver line search solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The Newton solver is not a line search solver.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_LINESEARCH_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_LINESEARCH_TYPE_SET")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_LINESEARCH_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of function evaluations for a nonlinear Newton solver
  SUBROUTINE SOLVER_NEWTON_MAXIMUM_FUNCTION_EVALUATIONS_SET(SOLVER,MAXIMUM_FUNCTION_EVALUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the maximum function evaluations for
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_FUNCTION_EVALUATIONS !<The maximum function evaluations to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_MAXIMUM_FUNCTION_EVALUATIONS_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(MAXIMUM_FUNCTION_EVALUATIONS>0) THEN
                  NEWTON_SOLVER%MAXIMUM_NUMBER_OF_FUNCTION_EVALUATIONS=MAXIMUM_FUNCTION_EVALUATIONS
                ELSE
                  LOCAL_ERROR="The specified maximum number of function evaluations of "// &
                    & TRIM(NUMBER_TO_VSTRING(MAXIMUM_FUNCTION_EVALUATIONS,"*",ERR,ERROR))// &
                    & " is invalid. The maximum number of function evaluations must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_MAXIMUM_FUNCTION_EVALUATIONS_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_MAXIMUM_FUNCTION_EVALUATIONS_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_MAXIMUM_FUNCTION_EVALUATIONS_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_MAXIMUM_FUNCTION_EVALUATIONS_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of iterations for a nonlinear Newton solver
  SUBROUTINE SOLVER_NEWTON_MAXIMUM_ITERATIONS_SET(SOLVER,MAXIMUM_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_ITERATIONS !<The maximum iterations to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_MAXIMUM_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(MAXIMUM_ITERATIONS>0) THEN
                  NEWTON_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
                ELSE
                  LOCAL_ERROR="The specified maximum iterations of "//TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                    & " is invalid. The maximum number of iterations must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Nonlinear sovler Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_MAXIMUM_ITERATIONS_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_MAXIMUM_ITERATIONS_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_MAXIMUM_ITEATIONS_SET")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_MAXIMUM_ITERATIONS_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the relative tolerance for a nonlinear Newton solver
  SUBROUTINE SOLVER_NEWTON_RELATIVE_TOLERANCE_SET(SOLVER,RELATIVE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the relative tolerance for
    REAL(DP), INTENT(IN) :: RELATIVE_TOLERANCE !<The relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_RELATIVE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(RELATIVE_TOLERANCE>ZERO_TOLERANCE) THEN
                  NEWTON_SOLVER%RELATIVE_TOLERANCE=RELATIVE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified relative tolerance of "//TRIM(NUMBER_TO_VSTRING(RELATIVE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The relative tolerance must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_RELATIVE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_RELATIVE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_RELATIVE_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_RELATIVE_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution tolerance for a nonlinear Newton solver
  SUBROUTINE SOLVER_NEWTON_SOLUTION_TOLERANCE_SET(SOLVER,SOLUTION_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the solution tolerance for
    REAL(DP), INTENT(IN) :: SOLUTION_TOLERANCE !<The solution tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_SOLUTION_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(SOLUTION_TOLERANCE>ZERO_TOLERANCE) THEN
                  NEWTON_SOLVER%SOLUTION_TOLERANCE=SOLUTION_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified solution tolerance of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The relative tolerance must be > 0."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF            
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_SOLUTION_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_SOLUTION_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_SOLUTION_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_SOLUTION_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !Solves a nonlinear Newton solver 
  SUBROUTINE SOLVER_NEWTON_SOLVE(NEWTON_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer to the nonlinear Newton solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NEWTON_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(NEWTON_SOLVER)) THEN
      SELECT CASE(NEWTON_SOLVER%NEWTON_SOLVE_TYPE)
      CASE(SOLVER_NEWTON_LINESEARCH)
        CALL SOLVER_NEWTON_LINESEARCH_SOLVE(NEWTON_SOLVER%LINESEARCH_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_NEWTON_TRUSTREGION)
        CALL SOLVER_NEWTON_TRUSTREGION_SOLVE(NEWTON_SOLVER%TRUSTREGION_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The nonlinear solver type of "// &
          & TRIM(NUMBER_TO_VSTRING(NEWTON_SOLVER%NEWTON_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Newton solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_SOLVE")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_SOLVE
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nonlinear Newton trust region solver
  SUBROUTINE SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH(TRUSTREGION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER !<A pointer the nonlinear Newton trust region solver to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    EXTERNAL :: PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: RESIDUAL_VECTOR
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
  
    CALL ENTERS("SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(TRUSTREGION_SOLVER)) THEN
      NEWTON_SOLVER=>TRUSTREGION_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN
        NONLINEAR_SOLVER=>NEWTON_SOLVER%NONLINEAR_SOLVER
        IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
          SOLVER=>NONLINEAR_SOLVER%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SELECT CASE(TRUSTREGION_SOLVER%SOLVER_LIBRARY)
              CASE(SOLVER_CMISS_LIBRARY)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(SOLVER_PETSC_LIBRARY)
                !Create the solver matrices and vectors
                CALL SOLVER_MATRICES_CREATE_START(SOLVER_EQUATIONS,SOLVER_MATRICES,ERR,ERROR,*999)
                CALL SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
!!TODO: set up the matrix structure if using an analytic Jacobian            
                CALL SOLVER_MATRICES_CREATE_FINISH(SOLVER_MATRICES,ERR,ERROR,*999)
                !Create the PETSc SNES solver
                CALL PETSC_SNESCREATE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,TRUSTREGION_SOLVER%SNES,ERR,ERROR,*999)
                !Set the nonlinear solver type to be a Newton trust region solver
                CALL PETSC_SNESSETTYPE(TRUSTREGION_SOLVER%SNES,PETSC_SNESTR,ERR,ERROR,*999)
                !Set the nonlinear function
                RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                IF(ASSOCIATED(RESIDUAL_VECTOR)) THEN
                  IF(ASSOCIATED(RESIDUAL_VECTOR%PETSC)) THEN
                    CALL PETSC_SNESSETFUNCTION(TRUSTREGION_SOLVER%SNES,RESIDUAL_VECTOR%PETSC%VECTOR, &
                      & PROBLEM_SOLVER_RESIDUAL_EVALUATE_PETSC,SOLVER,ERR,ERROR,*999)
                    CALL FLAG_ERROR("The residual vector PETSc is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrices residual vector is not associated.",ERR,ERROR,*999)
                ENDIF
                !Set the Jacobian if necessary
                !Set the trust region delta ???
                
                !Set the trust region tolerance
                CALL PETSC_SNESSETTRUSTREGIONTOLERANCE(TRUSTREGION_SOLVER%SNES,TRUSTREGION_SOLVER%TRUSTREGION_TOLERANCE, &
                  & ERR,ERROR,*999)
                !Set the tolerances for the SNES solver
                CALL PETSC_SNESSETTOLERANCES(TRUSTREGION_SOLVER%SNES,NEWTON_SOLVER%ABSOLUTE_TOLERANCE, &
                  & NEWTON_SOLVER%RELATIVE_TOLERANCE,NEWTON_SOLVER%SOLUTION_TOLERANCE, &
                  & NEWTON_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS,NEWTON_SOLVER%MAXIMUM_NUMBER_OF_FUNCTION_EVALUATIONS, &
                  & ERR,ERROR,*999)
                !Set any further SNES options from the command line options
                CALL PETSC_SNESSETFROMOPTIONS(TRUSTREGION_SOLVER%SNES,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solver library type of "// &
                  & TRIM(NUMBER_TO_VSTRING(TRUSTREGION_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Nonlinear solver solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Newton solver nonlinear solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Trust region Newton solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Trust region solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_TRUSTREGION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Sets/changes the trust region delta0 for a nonlinear Newton trust region solver solver
  SUBROUTINE SOLVER_NEWTON_TRUSTREGION_DELTA0_SET(SOLVER,TRUSTREGION_DELTA0,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the trust region delta0 for
    REAL(DP), INTENT(IN) :: TRUSTREGION_DELTA0 !<The trust region delta0 to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_TRUSTREGION_DELTA0_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NEWTON_SOLVER%NEWTON_SOLVE_TYPE==SOLVER_NEWTON_TRUSTREGION) THEN
                  TRUSTREGION_SOLVER=>NEWTON_SOLVER%TRUSTREGION_SOLVER
                  IF(ASSOCIATED(TRUSTREGION_SOLVER)) THEN
                    IF(TRUSTREGION_DELTA0>ZERO_TOLERANCE) THEN
                      TRUSTREGION_SOLVER%TRUSTREGION_DELTA0=TRUSTREGION_DELTA0
                    ELSE
                      LOCAL_ERROR="The specified trust region delta0 of "// &
                        & TRIM(NUMBER_TO_VSTRING(TRUSTREGION_DELTA0,"*",ERR,ERROR))// &
                        & " is invalid. The trust region delta0 must be > 0."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The Newton solver trust region solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The Newton solver is not a trust region solver.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_DELTA0_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_TRUSTREGION_DELTA0_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_DELTA0_SET")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_TRUSTREGION_DELTA0_SET
        
  !
  !================================================================================================================================
  !
  
  !>Finalise a nonlinear Newton trust region solver and deallocate all memory
  SUBROUTINE SOLVER_NEWTON_TRUSTREGION_FINALISE(TRUSTREGION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER !<A pointer the non linear trust region solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    CALL ENTERS("SOLVER_NEWTON_TRUSTREGION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TRUSTREGION_SOLVER)) THEN      
      CALL PETSC_SNESFINALISE(TRUSTREGION_SOLVER%SNES,ERR,ERROR,*999)
      DEALLOCATE(TRUSTREGION_SOLVER)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_TRUSTREGION_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_FINALISE")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_TRUSTREGION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a Newton trust region solver for a nonlinear solver
  SUBROUTINE SOLVER_NEWTON_TRUSTREGION_INITIALISE(NEWTON_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER !<A pointer the Newton solver to initialise the trust region solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
  
    CALL ENTERS("SOLVER_NEWTON_TRUSTREGION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(NEWTON_SOLVER)) THEN
      IF(ASSOCIATED(NEWTON_SOLVER%TRUSTREGION_SOLVER)) THEN
        CALL FLAG_ERROR("Trust region solver is already associated for this nonlinear solver.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(NEWTON_SOLVER%TRUSTREGION_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate Newton solver trust region solver.",ERR,ERROR,*999)
        NEWTON_SOLVER%TRUSTREGION_SOLVER%NEWTON_SOLVER=>NEWTON_SOLVER
        NEWTON_SOLVER%TRUSTREGION_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
!!TODO: set this properly
        NEWTON_SOLVER%TRUSTREGION_SOLVER%TRUSTREGION_DELTA0=0.01_DP
        CALL PETSC_SNESINITIALISE(NEWTON_SOLVER%TRUSTREGION_SOLVER%SNES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Newton solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_INITIALISE")
    RETURN
999 CALL SOLVER_NEWTON_TRUSTREGION_FINALISE(NEWTON_SOLVER%TRUSTREGION_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_NEWTON_TRUSTREGION_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWWTON_TRUSTREGION_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_TRUSTREGION_INITIALISE

  !
  !================================================================================================================================
  !

  !Solves a nonlinear Newton trust region solver 
  SUBROUTINE SOLVER_NEWTON_TRUSTREGION_SOLVE(TRUSTREGION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER !<A pointer to the nonlinear Newton trust region solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NEWTON_TRUSTREGION_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(TRUSTREGION_SOLVER)) THEN
      NEWTON_SOLVER=>TRUSTREGION_SOLVER%NEWTON_SOLVER
      IF(ASSOCIATED(NEWTON_SOLVER)) THEN        
        NONLINEAR_SOLVER=>NEWTON_SOLVER%NONLINEAR_SOLVER
        IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
          SOLVER=>NONLINEAR_SOLVER%SOLVER
          IF(ASSOCIATED(SOLVER)) THEN
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN            
                SELECT CASE(TRUSTREGION_SOLVER%SOLVER_LIBRARY)
                CASE(SOLVER_CMISS_LIBRARY)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(SOLVER_PETSC_LIBRARY)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)              
                CASE DEFAULT
                  LOCAL_ERROR="The nonlinear Newton trust region solver library type of "// &
                    & TRIM(NUMBER_TO_VSTRING(TRUSTREGION_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Solver matrices is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Nonlinear solver solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Newton solver nonlinear solver is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Trust region solver Newton solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Trust region solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_TRUSTREGION_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_SOLVE")
    RETURN 1
    
  END SUBROUTINE SOLVER_NEWTON_TRUSTREGION_SOLVE
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the trust region tolerance for a nonlinear Newton trust region solver
  SUBROUTINE SOLVER_NEWTON_TRUSTREGION_TOLERANCE_SET(SOLVER,TRUSTREGION_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the trust region tolerance for
    REAL(DP), INTENT(IN) :: TRUSTREGION_TOLERANCE !<The trust region tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NEWTON_TRUSTREGION_SOLVER_TYPE), POINTER :: TRUSTREGION_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_TRUSTREGION_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NEWTON_SOLVER%NEWTON_SOLVE_TYPE==SOLVER_NEWTON_TRUSTREGION) THEN
                  TRUSTREGION_SOLVER=>NEWTON_SOLVER%TRUSTREGION_SOLVER
                  IF(ASSOCIATED(TRUSTREGION_SOLVER)) THEN
                    IF(TRUSTREGION_TOLERANCE>ZERO_TOLERANCE) THEN
                      TRUSTREGION_SOLVER%TRUSTREGION_TOLERANCE=TRUSTREGION_TOLERANCE
                    ELSE
                      LOCAL_ERROR="The specified trust region tolerance of "// &
                        & TRIM(NUMBER_TO_VSTRING(TRUSTREGION_TOLERANCE,"*",ERR,ERROR))// &
                        & " is invalid. The trust region tolerance must be > 0."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The Newton solver trust region solver is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The Newton solver is not a trust region solver.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_NEWTON_TRUSTREGION_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_TRUSTREGION_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_TRUSTREGION_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of nonlinear Newton solver
  SUBROUTINE SOLVER_NEWTON_TYPE_SET(SOLVER,NEWTON_SOLVE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the nonlinear Newton solver type
    INTEGER(INTG), INTENT(IN) :: NEWTON_SOLVE_TYPE !<The type of nonlinear solver to set \see SOLVER_ROUTINES_NewtonSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: NEWTON_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NEWTON_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*998)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE==SOLVER_NONLINEAR_NEWTON) THEN
              NEWTON_SOLVER=>NONLINEAR_SOLVER%NEWTON_SOLVER
              IF(ASSOCIATED(NEWTON_SOLVER)) THEN
                IF(NEWTON_SOLVE_TYPE/=NEWTON_SOLVER%NEWTON_SOLVE_TYPE) THEN
                  !Intialise the new solver type
                  SELECT CASE(NEWTON_SOLVE_TYPE)
                  CASE(SOLVER_NEWTON_LINESEARCH)
                    CALL SOLVER_NEWTON_LINESEARCH_INITIALISE(NEWTON_SOLVER,ERR,ERROR,*999)
                  CASE(SOLVER_NEWTON_TRUSTREGION)
                    CALL SOLVER_NEWTON_TRUSTREGION_INITIALISE(NEWTON_SOLVER,ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The Newton solver type of "//TRIM(NUMBER_TO_VSTRING(NEWTON_SOLVE_TYPE,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  !Finalise the old solver type
                  SELECT CASE(NEWTON_SOLVER%NEWTON_SOLVE_TYPE)
                  CASE(SOLVER_NEWTON_LINESEARCH)
                    CALL SOLVER_NEWTON_LINESEARCH_FINALISE(NEWTON_SOLVER%LINESEARCH_SOLVER,ERR,ERROR,*999)
                  CASE(SOLVER_NEWTON_TRUSTREGION)
                    CALL SOLVER_NEWTON_TRUSTREGION_FINALISE(NEWTON_SOLVER%TRUSTREGION_SOLVER,ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The Newton solver type of "// &
                      & TRIM(NUMBER_TO_VSTRING(NEWTON_SOLVER%NEWTON_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  NEWTON_SOLVER%NEWTON_SOLVE_TYPE=NEWTON_SOLVE_TYPE
                ENDIF
              ELSE
                CALL FLAG_ERROR("Nonlinear solver Newton solver is not associated.",ERR,ERROR,*998)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The nonlinear solver is not a Newton solver.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_NEWTON_TYPE_SET")
    RETURN
999 SELECT CASE(NEWTON_SOLVE_TYPE)
    CASE(SOLVER_NEWTON_LINESEARCH)
      CALL SOLVER_NEWTON_LINESEARCH_FINALISE(NEWTON_SOLVER%LINESEARCH_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_NEWTON_TRUSTREGION)
      CALL SOLVER_NEWTON_TRUSTREGION_FINALISE(NEWTON_SOLVER%TRUSTREGION_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    END SELECT
998 CALL ERRORS("SOLVER_NEWTON_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NEWTON_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NEWTON_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a nonlinear solver 
  SUBROUTINE SOLVER_NONLINEAR_CREATE_FINISH(NONLINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer to the nonlinear solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NONLINEAR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      SELECT CASE(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE)
      CASE(SOLVER_NONLINEAR_NEWTON)
        CALL SOLVER_NEWTON_CREATE_FINISH(NONLINEAR_SOLVER%NEWTON_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(SOLVER_NONLINEAR_SQP)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The nonlinear solver type of "// &
          & TRIM(NUMBER_TO_VSTRING(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Nonlinear solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NONLINEAR_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_NONLINEAR_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_NONLINEAR_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_NONLINEAR_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finalise a nonlinear solver for a problem solver
  SUBROUTINE SOLVER_NONLINEAR_FINALISE(NONLINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer the nonlinear solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_NONLINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      CALL SOLVER_NEWTON_FINALISE(NONLINEAR_SOLVER%NEWTON_SOLVER,ERR,ERROR,*999)
      DEALLOCATE(NONLINEAR_SOLVER)
    ENDIF
         
    CALL EXITS("SOLVER_NONLINEAR_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_NONLINEAR_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NONLINEAR_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NONLINEAR_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a nonlinear solver for a solver.
  SUBROUTINE SOLVER_NONLINEAR_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_NONLINEAR_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%NONLINEAR_SOLVER)) THEN
        CALL FLAG_ERROR("Nonlinear solver is already associated for this solver.",ERR,ERROR,*998)
      ELSE
        !Allocate and initialise a Nonlinear solver
        ALLOCATE(SOLVER%NONLINEAR_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver nonlinear solver.",ERR,ERROR,*999)
        SOLVER%NONLINEAR_SOLVER%SOLVER=>SOLVER
        NULLIFY(SOLVER%NONLINEAR_SOLVER%NEWTON_SOLVER)
        !Default to a nonlinear Newton solver
        SOLVER%NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE=SOLVER_NONLINEAR_NEWTON
        CALL SOLVER_NEWTON_INITIALISE(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_NONLINEAR_INITIALISE")
    RETURN
999 CALL SOLVER_NONLINEAR_FINALISE(SOLVER%NONLINEAR_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_NONLINEAR_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NONLINEAR_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NONLINEAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Monitors the nonlinear solve.
  SUBROUTINE SOLVER_NONLINEAR_MONITOR(NONLINEAR_SOLVER,ITS,NORM,ERR,ERROR,*)

   !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer to the nonlinear solver to monitor
    INTEGER(INTG), INTENT(IN) :: ITS !<The number of iterations
    REAL(DP), INTENT(IN) :: NORM !<The residual norm
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("SOLVER_NONLINEAR_MONITOR",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
        
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Nonlinear solve monitor: ",ERR,ERROR,*999)
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Iteration number = ",ITS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Function Norm    = ",NORM,ERR,ERROR,*999)
      !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"    Number of function evaluations = ",NONLINEAR_SOLVER% &
      !  & TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS,ERR,ERROR,*999)
      !CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"    Number of Jacobian evaluations = ",NONLINEAR_SOLVER% &
      !  & TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS,ERR,ERROR,*999)            
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"",ERR,ERROR,*999)      
        
    ELSE
      CALL FLAG_ERROR("Nonlinear solver is not associated.",ERR,ERROR,*999)
    ENDIF
     
    CALL EXITS("SOLVER_NONLINEAR_MONITOR")
    RETURN
999 CALL ERRORS("SOLVER_NONLINEAR_MONITOR",ERR,ERROR)
    CALL EXITS("SOLVER_NONLINEAR_MONITOR")
    RETURN 1
  END SUBROUTINE SOLVER_NONLINEAR_MONITOR

  !
  !================================================================================================================================
  !

  !Solves a nonlinear solver 
  SUBROUTINE SOLVER_NONLINEAR_SOLVE(NONLINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer to the nonlinear solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NONLINEAR_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      SELECT CASE(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE)
      CASE(SOLVER_NONLINEAR_NEWTON)
        CALL SOLVER_NEWTON_SOLVE(NONLINEAR_SOLVER%NEWTON_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(SOLVER_NONLINEAR_SQP)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The nonlinear solver type of "// &
          & TRIM(NUMBER_TO_VSTRING(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Nonlinear solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NONLINEAR_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_NONLINEAR_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_NONLINEAR_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NONLINEAR_SOLVE
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of nonlinear solver
  SUBROUTINE SOLVER_NONLINEAR_TYPE_SET(SOLVER,NONLINEAR_SOLVE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the nonlinear solver type
    INTEGER(INTG), INTENT(IN) :: NONLINEAR_SOLVE_TYPE !<The type of nonlinear solver to set \see SOLVER_ROUTINES_NonlinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_NONLINEAR_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*998)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
            IF(NONLINEAR_SOLVE_TYPE/=NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE) THEN
              !Intialise the new solver type
              SELECT CASE(NONLINEAR_SOLVE_TYPE)
              CASE(SOLVER_NONLINEAR_NEWTON)
                CALL SOLVER_NEWTON_INITIALISE(NONLINEAR_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(SOLVER_NONLINEAR_SQP)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The nonlinear solver type of "//TRIM(NUMBER_TO_VSTRING(NONLINEAR_SOLVE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              !Finalise the old solver type
              SELECT CASE(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE)
              CASE(SOLVER_NONLINEAR_NEWTON)
                CALL SOLVER_NEWTON_FINALISE(NONLINEAR_SOLVER%NEWTON_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                
              CASE(SOLVER_NONLINEAR_SQP)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The nonlinear solver type of "// &
                  & TRIM(NUMBER_TO_VSTRING(NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              NONLINEAR_SOLVER%NONLINEAR_SOLVE_TYPE=NONLINEAR_SOLVE_TYPE
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver nonlinear solver is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a nonlinear solver.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_NONLINEAR_TYPE_SET")
    RETURN
999 SELECT CASE(NONLINEAR_SOLVE_TYPE)
    CASE(SOLVER_NONLINEAR_NEWTON)
      CALL SOLVER_NEWTON_FINALISE(NONLINEAR_SOLVER%NEWTON_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*998)                
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*998)      
    END SELECT
998 CALL ERRORS("SOLVER_NONLINEAR_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_NONLINEAR_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_NONLINEAR_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for a solver
  SUBROUTINE SOLVER_OUTPUT_TYPE_SET(SOLVER,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the output type for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The type of solver output to be set \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*999)
      ELSE        
        SELECT CASE(OUTPUT_TYPE)
        CASE(SOLVER_NO_OUTPUT)
          SOLVER%OUTPUT_TYPE=SOLVER_NO_OUTPUT
        CASE(SOLVER_PROGRESS_OUTPUT)
          SOLVER%OUTPUT_TYPE=SOLVER_PROGRESS_OUTPUT
        CASE(SOLVER_TIMING_OUTPUT)
          SOLVER%OUTPUT_TYPE=SOLVER_TIMING_OUTPUT
        CASE(SOLVER_SOLVER_OUTPUT)
          SOLVER%OUTPUT_TYPE=SOLVER_SOLVER_OUTPUT
        CASE(SOLVER_MATRIX_OUTPUT)
          SOLVER%OUTPUT_TYPE=SOLVER_MATRIX_OUTPUT         
        CASE DEFAULT
          LOCAL_ERROR="The specified solver output type of "// &
            & TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_OUTPUT_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_OUTPUT_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_OUTPUT_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_OUTPUT_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver equations for a solver.
  SUBROUTINE SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLVER_SOLVER_EQUATIONS_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN 
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          CALL FLAG_ERROR("Solver equations is already associated.",ERR,ERROR,*998)
        ELSE
          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
          IF(.NOT.ASSOCIATED(SOLVER_EQUATIONS)) CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solvers has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLVER_SOLVER_EQUATIONS_GET")
    RETURN
999 NULLIFY(SOLVER_EQUATIONS)
998 CALL ERRORS("SOLVER_SOLVER_EQUATIONS_GET",ERR,ERROR)
    CALL EXITS("SOLVER_SOLVER_EQUATIONS_GET")
    RETURN 1
    
  END SUBROUTINE SOLVER_SOLVER_EQUATIONS_GET
  
  !
  !================================================================================================================================
  !

  !>Solve the problem
  SUBROUTINE SOLVER_SOLVE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(SP) :: SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),USER_ELAPSED,USER_TIME1(1),USER_TIME2(1)
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
          CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
          CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)          
        ENDIF
        !Solve the system depending on the solver type
        SELECT CASE(SOLVER%SOLVE_TYPE)
        CASE(SOLVER_LINEAR_TYPE)
          !Assemble the solver matrices
!!TODO: Work out what to assemble
          CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_LINEAR_ONLY,ERR,ERROR,*999)
          !If required output the solver matrices          
          IF(SOLVER%OUTPUT_TYPE>=SOLVER_MATRIX_OUTPUT) THEN
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
              IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRICES,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver equations solver matrices are not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Solve linear system
          CALL SOLVER_LINEAR_SOLVE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_NONLINEAR_TYPE)
          CALL SOLVER_NONLINEAR_SOLVE(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_DYNAMIC_TYPE)
          CALL SOLVER_DYNAMIC_SOLVE(SOLVER%DYNAMIC_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_INTEGRATION_TYPE)
          CALL SOLVER_INTEGRATION_SOLVE(SOLVER%INTEGRATION_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_EIGENPROBLEM_TYPE)
          CALL SOLVER_EIGENPROBLEM_SOLVE(SOLVER%EIGENPROBLEM_SOLVER,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        !If necessary output the timing information
        IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
          CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
          CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
          USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
          SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for solve = ",USER_ELAPSED, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total System time for solve = ",SYSTEM_ELAPSED, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_SOLVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a solver
  SUBROUTINE SOLVER_TYPE_SET(SOLVER,SOLVE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: SOLVE_TYPE !<The type of solver to be set \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*998)
      ELSE
        IF(SOLVE_TYPE/=SOLVER%SOLVE_TYPE) THEN
          !Initialise the new solver type 
          SELECT CASE(SOLVE_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_LINEAR_INITIALISE(SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_NONLINEAR_INITIALISE(SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_DYNAMIC_TYPE)
            CALL SOLVER_DYNAMIC_INITIALISE(SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_INTEGRATION_TYPE)
            CALL SOLVER_INTEGRATION_INITIALISE(SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_EIGENPROBLEM_TYPE)
            CALL SOLVER_EIGENPROBLEM_INITIALISE(SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The specified solve type of "//TRIM(NUMBER_TO_VSTRING(SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Finalise the old solve type
          SELECT CASE(SOLVER%SOLVE_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_LINEAR_FINALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_NONLINEAR_FINALISE(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_DYNAMIC_TYPE)
            CALL SOLVER_DYNAMIC_FINALISE(SOLVER%DYNAMIC_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_INTEGRATION_TYPE)
            CALL SOLVER_INTEGRATION_FINALISE(SOLVER%INTEGRATION_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_EIGENPROBLEM_TYPE)
            CALL SOLVER_EIGENPROBLEM_FINALISE(SOLVER%EIGENPROBLEM_SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver solve type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Set the solve type
          SOLVER%SOLVE_TYPE=SOLVE_TYPE
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_TYPE_SET")
    RETURN
999 SELECT CASE(SOLVE_TYPE)
    CASE(SOLVER_LINEAR_TYPE)
      CALL SOLVER_LINEAR_FINALISE(SOLVER%LINEAR_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_NONLINEAR_TYPE)
      CALL SOLVER_NONLINEAR_FINALISE(SOLVER%NONLINEAR_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_DYNAMIC_TYPE)
      CALL SOLVER_DYNAMIC_FINALISE(SOLVER%DYNAMIC_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_INTEGRATION_TYPE)
      CALL SOLVER_INTEGRATION_FINALISE(SOLVER%INTEGRATION_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      CALL SOLVER_EIGENPROBLEM_FINALISE(SOLVER%EIGENPROBLEM_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    END SELECT
998 CALL ERRORS("SOLVER_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Updates the dependent variables from the solver solution
  SUBROUTINE SOLVER_VARIABLES_UPDATE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to update the variables from
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,equations_set_idx,field_dof,solver_dof_idx,solver_matrix_idx,variable_dof
    REAL(DP) :: additive_constant,VALUE,variable_coefficient
    REAL(DP), POINTER :: SOLVER_DATA(:)
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOLVER_VECTOR
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(SOLVER_DATA)
    
    CALL ENTERS("SOLVER_VARIABLES_UPDATE",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            SOLVER_MAPPING=>SOLVER_MATRICES%SOLVER_MAPPING
            IF(ASSOCIATED(SOLVER_MAPPING)) THEN            
              DO solver_matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
                SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
                  IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                    !Get the solver variables data                  
                    CALL DISTRIBUTED_VECTOR_DATA_GET(SOLVER_VECTOR,SOLVER_DATA,ERR,ERROR,*999)             
                    !Loop over the solver variable dofs
                    DO solver_dof_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)%NUMBER_OF_DOFS
                      !Loop over the equations sets associated with this dof
                      DO equations_set_idx=1,SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                        & SOLVER_DOF_TO_VARIABLE_MAPS(solver_dof_idx)%NUMBER_OF_EQUATIONS_SETS                        
                        DEPENDENT_VARIABLE=>SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                          & SOLVER_DOF_TO_VARIABLE_MAPS(solver_dof_idx)%VARIABLE(equations_set_idx)%PTR
                        IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                          DEPENDENT_FIELD=>DEPENDENT_VARIABLE%FIELD
                          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                            !Get the dependent field dof the solver dof is mapped to
                            variable_dof=SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_DOF_TO_VARIABLE_MAPS(solver_dof_idx)%VARIABLE_DOF(equations_set_idx)
                            field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                            variable_coefficient=SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_DOF_TO_VARIABLE_MAPS(solver_dof_idx)%VARIABLE_COEFFICIENT(equations_set_idx)
                            additive_constant=SOLVER_MAPPING%SOLVER_COL_TO_EQUATIONS_SETS_MAP(solver_matrix_idx)% &
                              & SOLVER_DOF_TO_VARIABLE_MAPS(solver_dof_idx)%ADDITIVE_CONSTANT(equations_set_idx)
                            !Set the dependent field dof
                            VALUE=SOLVER_DATA(solver_dof_idx)*variable_coefficient+additive_constant
                            CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,field_dof,VALUE, &
                              & ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !equations_set_idx
                    ENDDO !solver_dof_idx
                    !Restore the solver dof data
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(SOLVER_VECTOR,SOLVER_DATA,ERR,ERROR,*999)
                    !Start the transfer of the field dofs
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDDO !equations_set_idx
                    !Finish the transfer of the field dofs
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                    ENDDO !equations_set_idx
                  ELSE
                    CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*998)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*998)
                ENDIF
              ENDDO !solver_matrix_idx
            ELSE
              CALL FLAG_ERROR("Solver matrices solution mapping is not associated.",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver equations solver matrices are not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_VARIABLES_UPDATE")
    RETURN
999 IF(ASSOCIATED(SOLVER_DATA)) CALL DISTRIBUTED_VECTOR_DATA_RESTORE(SOLVER_VECTOR,SOLVER_DATA,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_VARIABLES_UPDATE",ERR,ERROR)    
    CALL EXITS("SOLVER_VARIABLES_UPDATE")
    RETURN 1
   
  END SUBROUTINE SOLVER_VARIABLES_UPDATE

  !
  !================================================================================================================================
  !

  !>Finish the creation of solvers.
  SUBROUTINE SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solver_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
   
    CALL ENTERS("SOLVERS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVERS)) THEN
      IF(SOLVERS%SOLVERS_FINISHED) THEN
        CALL FLAG_ERROR("Solvers has already been finished.",ERR,ERROR,*999)
      ELSE        
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN          
          !Finish the solver creation
          IF(ALLOCATED(SOLVERS%SOLVERS)) THEN
            DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
              SOLVER=>SOLVERS%SOLVERS(solver_idx)%PTR
              IF(ASSOCIATED(SOLVER)) THEN
                CALL SOLVER_CREATE_FINISH(SOLVER,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !solver_idx            
            SOLVERS%SOLVERS_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Solvers solvers is not allocated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solvers control loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVERS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVERS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("SOLVERS_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE SOLVERS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a solvers for the control loop. 
  SUBROUTINE SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to create the solvers for
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<On exit, a pointer to the solvers. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVERS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS==0) THEN
          IF(ASSOCIATED(SOLVERS)) THEN
            CALL FLAG_ERROR("Solvers is already associated.",ERR,ERROR,*999)
          ELSE
            NULLIFY(SOLVERS)
            !Initialise the solvers
            CALL SOLVERS_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
            !Return the pointer
            SOLVERS=>CONTROL_LOOP%SOLVERS
          ENDIF
        ELSE
          LOCAL_ERROR="Invalid control loop setup. The specified control loop has "// &
            & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))// &
            & " sub loops. To create solvers the control loop must have 0 sub loops."          
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Control loop has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVERS_CREATE_START")
    RETURN
999 CALL ERRORS("SOLVERS_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLVERS_CREATE_START")
    RETURN 1
  END SUBROUTINE SOLVERS_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroys the solvers
  SUBROUTINE SOLVERS_DESTROY(SOLVERS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVERS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVERS)) THEN
      CALL SOLVERS_FINALISE(SOLVERS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("SOLVERS_DESTROY")
    RETURN
999 CALL ERRORS("SOLVERS_DESTROY",ERR,ERROR)
    CALL EXITS("SOLVERS_DESTROY")
    RETURN 1
    
  END SUBROUTINE SOLVERS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the solvers and deallocates all memory
  SUBROUTINE SOLVERS_FINALISE(SOLVERS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solver_idx
 
    CALL ENTERS("SOLVERS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVERS)) THEN
      IF(ALLOCATED(SOLVERS%SOLVERS)) THEN
        DO solver_idx=1,SIZE(SOLVERS%SOLVERS,1)
          CALL SOLVER_FINALISE(SOLVERS%SOLVERS(solver_idx)%PTR,ERR,ERROR,*999)
        ENDDO !solver_idx
        DEALLOCATE(SOLVERS%SOLVERS)
      ENDIF
      DEALLOCATE(SOLVERS)
    ENDIF
       
    CALL EXITS("SOLVERS_FINALISE")
    RETURN
999 CALL ERRORS("SOLVERS_FINALISE",ERR,ERROR)
    CALL EXITS("SOLVERS_FINALISE")
    RETURN 1
  END SUBROUTINE SOLVERS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the solvers for a control loop.
  SUBROUTINE SOLVERS_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to initialise the solvers for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,solver_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("SOLVERS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%SOLVERS)) THEN
        CALL FLAG_ERROR("Solvers is already allocated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%SOLVERS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate control loop solvers.",ERR,ERROR,*999)
        CONTROL_LOOP%SOLVERS%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%SOLVERS%SOLVERS_FINISHED=.FALSE.
        CONTROL_LOOP%SOLVERS%NUMBER_OF_SOLVERS=1
        ALLOCATE(CONTROL_LOOP%SOLVERS%SOLVERS(CONTROL_LOOP%SOLVERS%NUMBER_OF_SOLVERS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solvers solvers.",ERR,ERROR,*999)
        DO solver_idx=1,CONTROL_LOOP%SOLVERS%NUMBER_OF_SOLVERS
          NULLIFY(CONTROL_LOOP%SOLVERS%SOLVERS(solver_idx)%PTR)
          CALL SOLVER_INITIALISE(CONTROL_LOOP%SOLVERS,solver_idx,ERR,ERROR,*999)
        ENDDO !solver_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLVERS_INITIALISE")
    RETURN
999 CALL SOLVERS_FINALISE(CONTROL_LOOP%SOLVERS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVERS_INITIALISE",ERR,ERROR)
    CALL EXITS("SOLVERS_INITIALISE")
    RETURN 1
    
  END SUBROUTINE SOLVERS_INITIALISE
  
 
  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solvers.
  SUBROUTINE SOLVERS_NUMBER_SET(SOLVERS,NUMBER_OF_SOLVERS,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers to set the number for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SOLVERS !<The number of solvers to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solver_idx
    TYPE(SOLVER_PTR_TYPE), ALLOCATABLE :: OLD_SOLVERS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("SOLVERS_NUMBER_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVERS)) THEN
      IF(SOLVERS%SOLVERS_FINISHED) THEN
        CALL FLAG_ERROR("Solvers have already been finished.",ERR,ERROR,*998)
      ELSE
        IF(NUMBER_OF_SOLVERS>0) THEN
          IF(NUMBER_OF_SOLVERS/=SOLVERS%NUMBER_OF_SOLVERS) THEN
            ALLOCATE(OLD_SOLVERS(SOLVERS%NUMBER_OF_SOLVERS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old solvers.",ERR,ERROR,*999)
            DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
              OLD_SOLVERS(solver_idx)%PTR=>SOLVERS%SOLVERS(solver_idx)%PTR
            ENDDO !solver_idx
            IF(ALLOCATED(SOLVERS%SOLVERS)) DEALLOCATE(SOLVERS%SOLVERS)
            ALLOCATE(SOLVERS%SOLVERS(NUMBER_OF_SOLVERS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solvers.",ERR,ERROR,*999)
            IF(NUMBER_OF_SOLVERS>SOLVERS%NUMBER_OF_SOLVERS) THEN
              DO solver_idx=1,SOLVERS%NUMBER_OF_SOLVERS
                SOLVERS%SOLVERS(solver_idx)%PTR=>OLD_SOLVERS(solver_idx)%PTR
              ENDDO !solver_idx
              DO solver_idx=SOLVERS%NUMBER_OF_SOLVERS+1,NUMBER_OF_SOLVERS
                NULLIFY(SOLVERS%SOLVERS(solver_idx)%PTR)
                CALL SOLVER_INITIALISE(SOLVERS,solver_idx,ERR,ERROR,*999)
              ENDDO !solution_idx
            ELSE
              DO solver_idx=1,NUMBER_OF_SOLVERS
                SOLVERS%SOLVERS(solver_idx)%PTR=>OLD_SOLVERS(solver_idx)%PTR
              ENDDO !solver_idx
              DO solver_idx=NUMBER_OF_SOLVERS+1,SOLVERS%NUMBER_OF_SOLVERS
                CALL SOLVER_FINALISE(OLD_SOLVERS(solver_idx)%PTR,ERR,ERROR,*999)
              ENDDO !solver_idx
            ENDIF
            SOLVERS%NUMBER_OF_SOLVERS=NUMBER_OF_SOLVERS
          ENDIF
        ELSE
          LOCAL_ERROR="The specified number of solvers of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SOLVERS,"*",ERR,ERROR))// &
            & " is invalid. The number of solvers must be > 0."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLVERS_NUMBER_SET")
    RETURN
999 IF(ALLOCATED(OLD_SOLVERS)) DEALLOCATE(OLD_SOLVERS)
998 CALL ERRORS("SOLVERS_NUMBER_SET",ERR,ERROR)
    CALL EXITS("SOLVERS_NUMBER_SET")
    RETURN 1
    
  END SUBROUTINE SOLVERS_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified solver in the list of solvers.
  SUBROUTINE SOLVERS_SOLVER_GET(SOLVERS,SOLVER_INDEX,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS !<A pointer to the solvers to get the solver for
    INTEGER(INTG), INTENT(IN) :: SOLVER_INDEX !<The specified solver to get
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On exit, a pointer to the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("SOLVERS_SOLVER_GET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVERS)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        CALL FLAG_ERROR("Solver is already associated.",ERR,ERROR,*998)
      ELSE
        NULLIFY(SOLVER)
        IF(SOLVER_INDEX>0.AND.SOLVER_INDEX<=SOLVERS%NUMBER_OF_SOLVERS) THEN
          IF(ALLOCATED(SOLVERS%SOLVERS)) THEN
            SOLVER=>SOLVERS%SOLVERS(SOLVER_INDEX)%PTR
            IF(.NOT.ASSOCIATED(SOLVER)) CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Solvers solvers is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified solver index of "//TRIM(NUMBER_TO_VSTRING(SOLVER_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The solver index must be >= 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(SOLVERS%NUMBER_OF_SOLVERS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("SOLVERS_SOLVER_GET")
    RETURN
999 NULLIFY(SOLVER)
998 CALL ERRORS("SOLVERS_SOLVER_GET",ERR,ERROR)
    CALL EXITS("SOLVERS_SOLVER_GET")
    RETURN 1
    
  END SUBROUTINE SOLVERS_SOLVER_GET
  
  !
  !================================================================================================================================
  !
        
END MODULE SOLVER_ROUTINES

!
!================================================================================================================================
!

!>Called from the PETSc TS solvers to monitor the dynamic solver
SUBROUTINE SOLVER_DYNAMIC_MONITOR_PETSC(TS,STEPS,TIME,X,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE ISO_VARYING_STRING
  USE KINDS
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_TS_TYPE), INTENT(INOUT) :: TS !<The PETSc TS type
  INTEGER(INTG), INTENT(INOUT) :: STEPS !<The iteration number
  REAL(DP), INTENT(INOUT) :: TIME !<The current time
  TYPE(PETSC_VEC_TYPE), INTENT(INOUT) :: X !<The current iterate
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
  TYPE(VARYING_STRING) :: ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    IF(CTX%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
      DYNAMIC_SOLVER=>CTX%DYNAMIC_SOLVER

      CALL SOLVER_DYNAMIC_MONITOR(DYNAMIC_SOLVER,STEPS,TIME,ERR,ERROR,*999)

    ELSE
      LOCAL_ERROR="Invalid solve type. The solve type of "//TRIM(NUMBER_TO_VSTRING(CTX%SOLVE_TYPE,"*",ERR,ERROR))// &
        & " does not correspond to a dynamic solver."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF      
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",ERR,ERROR,*999)
  ENDIF
  
  RETURN
999 CALL WRITE_ERROR(ERR,ERROR,*998)
998 CALL FLAG_WARNING("Error monitoring dynamic solve.",ERR,ERROR,*997)
997 RETURN    
END SUBROUTINE SOLVER_DYNAMIC_MONITOR_PETSC

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to monitor the Newton nonlinear solver
SUBROUTINE SOLVER_NONLINEAR_MONITOR_PETSC(SNES,ITS,NORM,CTX,ERR)

  USE BASE_ROUTINES
  USE CMISS_PETSC_TYPES
  USE ISO_VARYING_STRING
  USE KINDS
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PETSC_SNES_TYPE), INTENT(INOUT) :: SNES !<The PETSc SNES type
  INTEGER(INTG), INTENT(INOUT) :: ITS !<The iteration number
  REAL(DP), INTENT(INOUT) :: NORM !<The residual norm
  TYPE(SOLVER_TYPE), POINTER :: CTX !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
  !Local Variables
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
  TYPE(VARYING_STRING) :: ERROR,LOCAL_ERROR

  IF(ASSOCIATED(CTX)) THEN
    IF(CTX%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
      NONLINEAR_SOLVER=>CTX%NONLINEAR_SOLVER

      CALL SOLVER_NONLINEAR_MONITOR(NONLINEAR_SOLVER,ITS,NORM,ERR,ERROR,*999)

    ELSE
      LOCAL_ERROR="Invalid solve type. The solve type of "//TRIM(NUMBER_TO_VSTRING(CTX%SOLVE_TYPE,"*",ERR,ERROR))// &
        & " does not correspond to a nonlinear solver."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF      
  ELSE
    CALL FLAG_ERROR("Solver context is not associated.",ERR,ERROR,*999)
  ENDIF
  
  RETURN
999 CALL WRITE_ERROR(ERR,ERROR,*998)
998 CALL FLAG_WARNING("Error monitoring nonlinear solve.",ERR,ERROR,*997)
997 RETURN    
END SUBROUTINE SOLVER_NONLINEAR_MONITOR_PETSC

