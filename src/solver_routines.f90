!> \file
!> $Id: solver_routines.f90 20 2007-05-28 20:22:52Z cpb $
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
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_TYPE=1 !<Linear solution solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_NONLINEAR_TYPE=2 !<A nonlinear solution solver  \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_TIME_INTEGRATION_TYPE=3 !<A time integration solver \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_EIGENPROBLEM_TYPE=4 !<A eigenproblem type \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
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

  !> \addtogroup SOLVER_ROUTINES_OutputTypes SOLVER_ROUTINES::OutputTypes
  !> \brief The types of output
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_NO_OUTPUT=0 !<No output from the solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_TIMING_OUTPUT=1 !<Timing output from the solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_SOLVER_OUTPUT=2 !<Timing and solver specific output from the solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_MATRIX_OUTPUT=3 !<Timing and solver specific output and solution matrices output from the solver routines \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SparsityTypes SOLVER_ROUTINES::SparsityTypes
  !> \brief The types of sparse solver matrices
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
  !>@}

  !Integration procedures
  INTEGER(INTG), PARAMETER :: SOLVER_EULER_INTEGRATOR=1
  INTEGER(INTG), PARAMETER :: SOLVER_IMPROVED_EULER_INTEGRATOR=2
  INTEGER(INTG), PARAMETER :: SOLVER_4TH_RUNGE_KUTTA_INTEGRATOR=3
  INTEGER(INTG), PARAMETER :: SOLVER_ADAMS_MOULTON_INTEGERATOR=4
  INTEGER(INTG), PARAMETER :: SOLVER_LSODA_INTEGRATOR=5
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC SOLVER_LINEAR_TYPE,SOLVER_NONLINEAR_TYPE,SOLVER_TIME_INTEGRATION_TYPE

  PUBLIC SOLVER_CMISS_LIBRARY,SOLVER_PETSC_LIBRARY

  PUBLIC SOLVER_LINEAR_DIRECT_SOLVE_TYPE,SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
  
  PUBLIC SOLVER_DIRECT_LU,SOLVER_DIRECT_CHOLESKY,SOLVER_DIRECT_SVD

  PUBLIC SOLVER_ITERATIVE_RICHARDSON,SOLVER_ITERATIVE_CHEBYCHEV,SOLVER_ITERATIVE_CONJUGATE_GRADIENT, &
    & SOLVER_ITERATIVE_BICONJUGATE_GRADIENT,SOLVER_ITERATIVE_GMRES,SOLVER_ITERATIVE_BiCGSTAB,SOLVER_ITERATIVE_CONJGRAD_SQUARED

  PUBLIC SOLVER_ITERATIVE_NO_PRECONDITIONER,SOLVER_ITERATIVE_JACOBI_PRECONDITIONER,SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER, &
    & SOLVER_ITERATIVE_SOR_PRECONDITIONER,SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER, &
    & SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER,SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER

  PUBLIC SOLVER_NO_OUTPUT,SOLVER_TIMING_OUTPUT,SOLVER_SOLVER_OUTPUT,SOLVER_MATRIX_OUTPUT

  PUBLIC SOLVER_SPARSE_MATRICES,SOLVER_FULL_MATRICES

  PUBLIC SOLVER_OUTPUT_TYPE_SET,SOLVER_SPARSITY_TYPE_SET

  PUBLIC SOLVER_EULER_INTEGRATOR,SOLVER_IMPROVED_EULER_INTEGRATOR,SOLVER_4TH_RUNGE_KUTTA_INTEGRATOR, &
    & SOLVER_ADAMS_MOULTON_INTEGERATOR
  
  PUBLIC SOLVER_CREATE_FINISH,SOLVER_CREATE_START,SOLVER_DESTROY,SOLVER_LIBRARY_SET,SOLVER_SOLVE

  PUBLIC SOLVER_LINEAR_TYPE_SET
  
  PUBLIC SOLVER_LINEAR_DIRECT_TYPE_SET

  PUBLIC SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET,SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET, &
    & SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET,SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET, &
    & SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET,SOLVER_LINEAR_ITERATIVE_TYPE_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver for a problem solution
  SUBROUTINE SOLVER_CREATE_FINISH(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("SOLVER_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished.",ERR,ERROR,*998)
      ELSE
        SELECT CASE(SOLVER%SOLVE_TYPE)
        CASE(SOLVER_LINEAR_TYPE)
          CALL SOLVER_LINEAR_CREATE_FINISH(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_NONLINEAR_TYPE)
          CALL SOLVER_NONLINEAR_CREATE_FINISH(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_TIME_INTEGRATION_TYPE)
          CALL SOLVER_TIME_INTEGRATION_CREATE_FINISH(SOLVER%TIME_INTEGRATION_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_EIGENPROBLEM_TYPE)
          CALL SOLVER_EIGENPROBLEM_CREATE_FINISH(SOLVER%EIGENPROBLEM_SOLVER,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        SOLVER%SOLVER_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_CREATE_FINISH")
    RETURN
999 CALL SOLVER_FINALISE(SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a solver for a problem solution.
  SUBROUTINE SOLVER_CREATE_START(SOLUTION,SOLVE_TYPE,SOLVER,ERR,ERROR,*)
    
    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer to the solution to create the solver for.
    INTEGER(INTG), INTENT(IN) :: SOLVE_TYPE !<The type of solver to create \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<On exit, a pointer to the problem solver.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLUTION)) THEN             
      IF(ASSOCIATED(SOLVER)) THEN
        CALL FLAG_ERROR("Solver is already associated.",ERR,ERROR,*998)
      ELSE
        NULLIFY(SOLVER)
        !Initialise the solver
        CALL SOLVER_INITIALISE(SOLUTION,SOLVE_TYPE,ERR,ERROR,*999)
        SOLVER=>SOLUTION%SOLVER
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_CREATE_START")
    RETURN
999 CALL SOLVER_FINALISE(SOLUTION%SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLVER_CREATE_START")
    RETURN 1   
  END SUBROUTINE SOLVER_CREATE_START    
  
  !
  !================================================================================================================================
  !

  !>Destroy a problem solver.
  SUBROUTINE SOLVER_DESTROY(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("SOLVER_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      CALL SOLVER_FINALISE(SOLVER,ERR,ERROR,*999)      
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_DESTROY")
    RETURN
999 CALL ERRORS("SOLVER_DESTROY",ERR,ERROR)    
    CALL EXITS("SOLVER_DESTROY")
    RETURN 1
   
  END SUBROUTINE SOLVER_DESTROY

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

    CALL ENTERS("SOLVER_EIGENPROBLEM_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%EIGENPROBLEM_SOLVER)) THEN
        CALL FLAG_ERROR("Eigenproblem solver is already associated for this solver.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(SOLVER%EIGENPROBLEM_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver eigenproblem solver.",ERR,ERROR,*999)
        SOLVER%EIGENPROBLEM_SOLVER%SOLVER=>SOLVER
        SOLVER%EIGENPROBLEM_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_EIGENPROBLEM_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_EIGENPROBLEM_INITIALISE",ERR,ERROR)    
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
      CALL SOLVER_TIME_INTEGRATION_FINALISE(SOLVER%TIME_INTEGRATION_SOLVER,ERR,ERROR,*999)        
      CALL SOLVER_EIGENPROBLEM_FINALISE(SOLVER%EIGENPROBLEM_SOLVER,ERR,ERROR,*999)        
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

  !>Initialise a solver for a problem solution
  SUBROUTINE SOLVER_INITIALISE(SOLUTION,SOLVE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION !<A pointer the solution to initialise the solver for
    INTEGER(INTG), INTENT(IN) :: SOLVE_TYPE !<The type of solver to create \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLUTION)) THEN
      IF(ASSOCIATED(SOLUTION%SOLVER)) THEN
        CALL FLAG_ERROR("Solver is already associated for this solution.",ERR,ERROR,*999)
      ELSE
        SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
        IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
          IF(SOLUTION_MAPPING%SOLUTION_MAPPING_FINISHED) THEN
            ALLOCATE(SOLUTION%SOLVER,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solution solver.",ERR,ERROR,*999)
            SOLUTION%SOLVER%SOLUTION=>SOLUTION
            SOLUTION%SOLVER%SOLVER_FINISHED=.FALSE.
            SOLUTION%SOLVER%SOLUTION_MAPPING=>SOLUTION_MAPPING
            SOLUTION%SOLVER%OUTPUT_TYPE=SOLVER_NO_OUTPUT
            SOLUTION%SOLVER%SPARSITY_TYPE=SOLVER_SPARSE_MATRICES
            NULLIFY(SOLUTION%SOLVER%LINEAR_SOLVER)
            NULLIFY(SOLUTION%SOLVER%NONLINEAR_SOLVER)
            NULLIFY(SOLUTION%SOLVER%TIME_INTEGRATION_SOLVER)
            NULLIFY(SOLUTION%SOLVER%EIGENPROBLEM_SOLVER)
            SELECT CASE(SOLVE_TYPE)
            CASE(SOLVER_LINEAR_TYPE)
              SOLUTION%SOLVER%SOLVE_TYPE=SOLVER_LINEAR_TYPE
              CALL SOLVER_LINEAR_INITIALISE(SOLUTION%SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_NONLINEAR_TYPE)
              SOLUTION%SOLVER%SOLVE_TYPE=SOLVER_NONLINEAR_TYPE
              CALL SOLVER_NONLINEAR_INITIALISE(SOLUTION%SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_TIME_INTEGRATION_TYPE)
              SOLUTION%SOLVER%SOLVE_TYPE=SOLVER_TIME_INTEGRATION_TYPE
              CALL SOLVER_TIME_INTEGRATION_INITIALISE(SOLUTION%SOLVER,ERR,ERROR,*999)
            CASE(SOLVER_EIGENPROBLEM_TYPE)
              SOLUTION%SOLVER%SOLVE_TYPE=SOLVER_EIGENPROBLEM_TYPE
              CALL SOLVER_EIGENPROBLEM_INITIALISE(SOLUTION%SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The specified solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVE_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Solution mapping has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE        
          CALL FLAG_ERROR("Solution solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solution is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_INITIALISE

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
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: DIRECT_SOLVER
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: ITERATIVE_SOLVER
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER
    TYPE(TIME_INTEGRATION_SOLVER_TYPE), POINTER :: TIME_INTEGRATION_SOLVER
    TYPE(EIGENPROBLEM_SOLVER_TYPE), POINTER :: EIGENPROBLEM_SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LIBRARY_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(SOLVER%SOLVE_TYPE)
        CASE(SOLVER_LINEAR_TYPE)
          LINEAR_SOLVER=>SOLVER%LINEAR_SOLVER
          IF(ASSOCIATED(LINEAR_SOLVER)) THEN
            SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVER_TYPE)
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
              LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(SOLVER_NONLINEAR_TYPE)
          NONLINEAR_SOLVER=>SOLVER%NONLINEAR_SOLVER
          IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
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
            CALL FLAG_ERROR("Problem solver nonlinear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(SOLVER_TIME_INTEGRATION_TYPE)
          TIME_INTEGRATION_SOLVER=>SOLVER%TIME_INTEGRATION_SOLVER
          IF(ASSOCIATED(TIME_INTEGRATION_SOLVER)) THEN
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
            CALL FLAG_ERROR("Problem solver time integration solver is not associated.",ERR,ERROR,*999)
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
            CALL FLAG_ERROR("Problem solver eigenproblem solver is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The problem solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
     ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated.",ERR,ERROR,*999)
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
      SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVER_TYPE)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        CALL SOLVER_LINEAR_DIRECT_CREATE_FINISH(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        CALL SOLVER_LINEAR_ITERATIVE_CREATE_FINISH(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))// &
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
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_DIRECT_CREATE_FINISH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(LINEAR_DIRECT_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_DIRECT_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SELECT CASE(LINEAR_DIRECT_SOLVER%SOLVER_LIBRARY)
          CASE(SOLVER_CMISS_LIBRARY)
            !Create the solver matrices
            CALL SOLVER_MATRICES_CREATE_START(SOLVER,SOLVER_MATRICES,ERR,ERROR,*999)
            CALL SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            SELECT CASE(SOLVER%SPARSITY_TYPE)
            CASE(SOLVER_SPARSE_MATRICES)
              CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                & ERR,ERROR,*999)
            CASE(SOLVER_FULL_MATRICES)
              CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/), &
                & ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The specified solver sparsity type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SPARSITY_TYPE,"*",ERR,ERROR))// &
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

    CALL ENTERS("SOLVER_LINEAR_DIRECT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      IF(ASSOCIATED(LINEAR_SOLVER%DIRECT_SOLVER)) THEN
        CALL FLAG_ERROR("Direct solver is already associated for this linear solver.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(LINEAR_SOLVER%DIRECT_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear solver direct solver.",ERR,ERROR,*999)
        LINEAR_SOLVER%DIRECT_SOLVER%LINEAR_SOLVER=>LINEAR_SOLVER
        LINEAR_SOLVER%DIRECT_SOLVER%SOLVER_LIBRARY=SOLVER_CMISS_LIBRARY
        LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_LU
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_DIRECT_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_DIRECT_INITIALISE",ERR,ERROR)    
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
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_DIRECT_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_DIRECT_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_DIRECT_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
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
              LOCAL_ERROR="The given number of solver matrices of "// &
                & TRIM(NUMBER_TO_VSTRING(SOLVER_MATRICES%NUMBER_OF_MATRICES,"*",ERR,ERROR))// &
                & " is invalid. There should only be one solver matrix for a linear direct solver."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_DIRECT_SOLVE_TYPE) THEN
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
                    CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
                      & " is invalid"
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver direct solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear direct solver",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem solver linear solver is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The problem solver is not a linear solver",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
        CALL FLAG_ERROR("Linear solver is already associated for this problems solver.",ERR,ERROR,*999)
      ELSE
        SOLUTION=>SOLVER%SOLUTION
        IF(ASSOCIATED(SOLUTION)) THEN
          SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
          IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
            IF(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES==1) THEN
              ALLOCATE(SOLVER%LINEAR_SOLVER,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solver linear solver.",ERR,ERROR,*999)
              SOLVER%LINEAR_SOLVER%SOLVER=>SOLVER
              NULLIFY(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)
              NULLIFY(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)
              SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE=SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
              CALL SOLVER_LINEAR_ITERATIVE_INITIALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The number of solver matrices in the solution mapping of "// &
                & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
                & " is invalid for a linear solver. There should only be one solver matrix."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solution solution mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum absolute tolerance for an iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET(SOLVER,ABSOLUTE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set 
    REAL(DP), INTENT(IN) :: ABSOLUTE_TOLERANCE !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
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
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_ITERATIVE_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_ITERATIVE_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SELECT CASE(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY)
          CASE(SOLVER_CMISS_LIBRARY)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(SOLVER_PETSC_LIBRARY)
            !Create the solver matrices
            CALL SOLVER_MATRICES_CREATE_START(SOLVER,SOLVER_MATRICES,ERR,ERROR,*999)
            CALL SOLVER_MATRICES_LIBRARY_TYPE_SET(SOLVER_MATRICES,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
            SELECT CASE(SOLVER%SPARSITY_TYPE)
            CASE(SOLVER_SPARSE_MATRICES)
              CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                & ERR,ERROR,*999)
            CASE(SOLVER_FULL_MATRICES)
              CALL SOLVER_MATRICES_STORAGE_TYPE_SET(SOLVER_MATRICES,(/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/), &
                & ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The specified solver sparsity type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SPARSITY_TYPE,"*",ERR,ERROR))// &
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
            !Set any further KSP options
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
              & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
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
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
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

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      IF(ASSOCIATED(LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
        CALL FLAG_ERROR("Iterative solver is already associated for this linear solver.",ERR,ERROR,*999)
      ELSE
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
      CALL FLAG_ERROR("Linear solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_INITIALISE",ERR,ERROR)    
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
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
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
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
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
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
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
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_ITERATIVE_SOLVER)) THEN
      LINEAR_SOLVER=>LINEAR_ITERATIVE_SOLVER%LINEAR_SOLVER
      IF(ASSOCIATED(LINEAR_SOLVER)) THEN
        SOLVER=>LINEAR_SOLVER%SOLVER
        IF(ASSOCIATED(SOLVER)) THEN
          SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
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
                        CASE(PETSC_KSP_DIVERGED_ITS)
                          CALL FLAG_WARNING("Linear iterative solver did not converge. PETsc diverged its.",ERR,ERROR,*999)
                        CASE(PETSC_KSP_DIVERGED_DTOL)
                          CALL FLAG_WARNING("Linear iterative solver did not converge. PETsc diverged dtol.",ERR,ERROR,*999)
                        CASE(PETSC_KSP_DIVERGED_BREAKDOWN)
                          CALL FLAG_WARNING("Linear iterative solver did not converge. PETsc diverged breakdown.",ERR,ERROR,*999)
                        CASE(PETSC_KSP_DIVERGED_BREAKDOWN_BICG)
                          CALL FLAG_WARNING("Linear iterative solver did not converge. PETsc diverged breakdown BiCG.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_KSP_DIVERGED_NONSYMMETRIC)
                          CALL FLAG_WARNING("Linear iterative solver did not converge. PETsc diverged nonsymmetric.", &
                            & ERR,ERROR,*999)
                        CASE(PETSC_KSP_DIVERGED_INDEFINITE_PC)
                          CALL FLAG_WARNING("Linear iterative solver did not converge. PETsc diverged indefinite PC.", &
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
                          CASE(PETSC_KSP_CONVERGED_QCG_NEG_CURVE)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged QCG neg curve", &
                              & ERR,ERROR,*999)
                          CASE(PETSC_KSP_CONVERGED_QCG_CONSTRAINED)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged QCG constrained", &
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
                      CALL SOLVER_VARIABLES_UPDATE(SOLVER,ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) THEN
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
      SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVER_TYPE)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        CALL SOLVER_LINEAR_DIRECT_SOLVE(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        CALL SOLVER_LINEAR_ITERATIVE_SOLVE(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))// &
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
  SUBROUTINE SOLVER_LINEAR_TYPE_SET(SOLVER,LINEAR_SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to set the linear solver type
    INTEGER(INTG), INTENT(IN) :: LINEAR_SOLVER_TYPE !<The type of linear solver to set \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*998)
      ELSE
        IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
            IF(LINEAR_SOLVER_TYPE/=SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE) THEN
              !Intialise the new solver type
              SELECT CASE(LINEAR_SOLVER_TYPE)
              CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
                CALL SOLVER_LINEAR_DIRECT_INITIALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
                CALL SOLVER_LINEAR_ITERATIVE_INITIALISE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              !Finalise the old solver type
              SELECT CASE(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE)
              CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
                CALL SOLVER_LINEAR_DIRECT_FINALISE(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
                CALL SOLVER_LINEAR_ITERATIVE_FINALISE(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The linear solver type of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE=LINEAR_SOLVER_TYPE
            ENDIF
          ELSE
            CALL FLAG_ERROR("The solver linear solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The solver is not a linear solver.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_TYPE_SET")
    RETURN
999 IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%LINEAR_SOLVER)) THEN
        SELECT CASE(LINEAR_SOLVER_TYPE)
        CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
          CALL SOLVER_LINEAR_DIRECT_FINALISE(SOLVER%LINEAR_SOLVER%DIRECT_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
        CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
          CALL SOLVER_LINEAR_ITERATIVE_FINALISE(SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
        END SELECT
      ENDIF
    ENDIF
998 CALL ERRORS("SOLVER_LINEAR_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Assembles the solver matrices and rhs from the equations.
  SUBROUTINE SOLVER_MATRICES_ASSEMBLE(SOLVER,ERR,ERROR,*)

    !Argument variableg
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,equations_matrix_idx,equations_matrix_number, &
      & equations_row_number,equations_row_number2,equations_set_idx,EQUATIONS_STORAGE_TYPE,rhs_boundary_condition, &
      & rhs_field_dof,rhs_global_dof,rhs_variable_dof,rhs_variable_type,variable_boundary_condition,solver_column_idx, &
      & solver_column_number,solver_matrix_idx,solver_row_idx,solver_row_number,variable_dof,variable_field_dof, &
      & variable_global_dof,variable_idx,variable_type
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(SP) :: SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),USER_ELAPSED,USER_TIME1(1),USER_TIME2(1)
    REAL(DP) :: column_coupling_coefficient,MATRIX_VALUE,RHS_VALUE,row_coupling_coefficient,VALUE
    REAL(DP), POINTER :: DEPENDENT_PARAMETERS(:),EQUATIONS_MATRIX_DATA(:)
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: EQUATIONS_RHS_VECTOR,SOLVER_RHS_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: RHS_DOMAIN_MAPPING,VARIABLE_DOMAIN_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,RHS_VARIABLE
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(SOURCE_EQUATIONS_MATRICES_MAP_TYPE), POINTER :: SOURCE_MAPPINGS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("SOLVER_MATRICES_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOlUTION_MAPPING=>SOLVER%SOLUTION_MAPPING
      IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
        SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
        IF(ASSOCIATED(SOLVER_MATRICES)) THEN
          !Assemble solver matrices
          IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
            CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
            CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
          ENDIF
          !Assemble the solver matrices
          DO solver_matrix_idx=1,SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES
            SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
            IF(ASSOCIATED(SOLVER_MATRIX)) THEN
              IF(SOLVER_MATRIX%UPDATE_MATRIX) THEN              
                SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN                
                  !Initialise matrix to zero
                  CALL DISTRIBUTED_MATRIX_ALL_VALUES_SET(SOLVER_DISTRIBUTED_MATRIX,0.0_DP,ERR,ERROR,*999)
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    !Loop over the equations matrices
                    DO equations_matrix_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                      & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_EQUATIONS_MATRICES
                      EQUATIONS_TO_SOLVER_MAP=>SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                        & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                        & equations_matrix_idx)%PTR
                      IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                        EQUATIONS_MATRIX=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                        IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                          EQUATIONS_MATRICES=>EQUATIONS_MATRIX%EQUATIONS_MATRICES
                          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                            EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                            IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                              CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE, &
                                & ERR,ERROR,*999)
                              CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA,ERR,ERROR,*999)
                              SELECT CASE(EQUATIONS_STORAGE_TYPE)
                              CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                !Loop over the rows of the equations matrix
                                DO equations_row_number=1,EQUATIONS_MATRICES%NUMBER_OF_ROWS
                                  !Loop over the solution rows this equations row is mapped to
                                  DO solver_row_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
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
                                          & EQUATIONS_MATRICES%NUMBER_OF_ROWS)*row_coupling_coefficient*column_coupling_coefficient
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
                                  DO solver_row_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                    solver_row_number=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & SOLVER_ROWS(solver_row_idx)
                                    row_coupling_coefficient=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                      & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & COUPLING_COEFFICIENTS(solver_row_idx)
                                    !Loop over the columns of the equations matrix
                                    DO equations_column_idx=ROW_INDICES(equations_row_number),ROW_INDICES(equations_row_number+1)-1
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
                            CALL FLAG_ERROR("The equations matrix equations matrices is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("The equations matrix is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("The equations matrix equations to solver map is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDDO !equations_matrix_idx
                  ENDDO !equations_set_idx
                  !Update the solver matrix values
                  CALL DISTRIBUTED_MATRIX_UPDATE_START(SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                  CALL DISTRIBUTED_MATRIX_UPDATE_FINISH(SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF !Update matrix
            ELSE
              CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDDO !solver_matrix_idx
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
          !Assemble right hand side vector
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
              DO equations_set_idx=1,SOLUTION_MAPPING%NUMBER_OF_EQUATIONS_SETS
                EQUATIONS_SET=>SOLUTION_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                    EQUATIONS=>EQUATIONS_SET%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                      IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
                        IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
                          FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                          IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
                            SOURCE_MAPPINGS=>EQUATIONS_MAPPING%SOURCE_MAPPINGS
                            IF(ASSOCIATED(SOURCE_MAPPINGS)) THEN
!!TODO: source term
                            ENDIF
                            rhs_variable_type=EQUATIONS_MAPPING%RHS_VARIABLE_TYPE
                            RHS_VARIABLE=>EQUATIONS_MAPPING%RHS_VARIABLE
                            RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
                            EQUATIONS_RHS_VECTOR=>EQUATIONS_MATRICES%VECTOR
                            !Loop over the rows in the equations set
                            DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
!!TODO: what if the equations set doesn't have a RHS vector???
                              rhs_variable_dof=EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(equations_row_number)% &
                                & ROW_TO_RHS_DOF
                              rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                              rhs_field_dof=RHS_VARIABLE%GLOBAL_DOF_LIST(rhs_global_dof)
                              rhs_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(rhs_field_dof)
                              !Apply boundary conditions
                              SELECT CASE(rhs_boundary_condition)
                              CASE(EQUATIONS_SET_NOT_FIXED)
                                !Set Direchlet boundary conditions
                                !Loop over the dependent variables associated with this equations set row
                                DO variable_idx=1,EQUATIONS_MAPPING%NUMBER_OF_MATRIX_VARIABLES
                                  variable_type=EQUATIONS_MAPPING%MATRIX_VARIABLE_TYPES(variable_idx)
                                  DEPENDENT_VARIABLE=>EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS(variable_type)% &
                                    & VARIABLE
                                  VARIABLE_DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                                  variable_dof=EQUATIONS_MAPPING%EQUATIONS_ROW_TO_VARIABLES_MAPS(equations_row_number)% &
                                    & ROW_TO_DOFS_MAP(variable_idx)
                                  variable_global_dof=VARIABLE_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(variable_dof)
                                  variable_field_dof=DEPENDENT_VARIABLE%GLOBAL_DOF_LIST(variable_global_dof)
                                  variable_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(variable_field_dof)
                                  IF(variable_boundary_condition==EQUATIONS_SET_FIXED_BOUNDARY_CONDITION) THEN
                                    DO equations_matrix_idx=1,EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS( &
                                      & variable_type)%NUMBER_OF_EQUATIONS_MATRICES
                                      equations_matrix_number=EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS( &
                                        & variable_type)%EQUATIONS_MATRIX_NUMBERS(equations_matrix_idx)
                                      EQUATIONS_MATRIX=>EQUATIONS_MATRICES%MATRICES(equations_matrix_number)%PTR
                                      equations_column_number=EQUATIONS_MAPPING%VARIABLE_TO_EQUATIONS_MATRICES_MAPS( &
                                        & variable_type)%DOF_TO_COLUMNS_MAPS(equations_matrix_idx)%COLUMN_DOF(variable_dof)
                                      DO equations_row_number2=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                        CALL DISTRIBUTED_MATRIX_VALUES_GET(EQUATIONS_MATRIX%MATRIX,equations_row_number2, &
                                          & equations_column_number,MATRIX_VALUE,ERR,ERROR,*999)
                                        DO solver_row_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number2)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number2)%SOLVER_ROWS(solver_row_idx)
                                          row_coupling_coefficient=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                            & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number2)% &
                                            & COUPLING_COEFFICIENTS(solver_row_idx)
                                          VALUE=-1.0_DP*MATRIX_VALUE*DEPENDENT_PARAMETERS(variable_field_dof)* &
                                            & row_coupling_coefficient
                                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                            & ERR,ERROR,*999)
                                        ENDDO !solver_row_idx
                                      ENDDO !equations_row_number2
                                    ENDDO !matrix_idx
                                  ENDIF
                                ENDDO !variable_idx
                              CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                                !Set Neumann boundary conditions
                                !Loop over the solver rows associated with this equations set row
                                DO solver_row_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                  solver_row_number=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                  row_coupling_coefficient=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%COUPLING_COEFFICIENTS(solver_row_idx)
                                  VALUE=DEPENDENT_PARAMETERS(rhs_field_dof)*row_coupling_coefficient
                                  CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE,ERR,ERROR,*999)
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
                              !Add in equations RHS and source values
                              CALL DISTRIBUTED_VECTOR_VALUES_GET(EQUATIONS_RHS_VECTOR,equations_row_number,RHS_VALUE,ERR,ERROR,*999)
                              !Loop over the solver rows associated with this equations set row
                              DO solver_row_idx=1,SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                solver_row_number=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS(solver_row_idx)
                                row_coupling_coefficient=SOLUTION_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%COUPLING_COEFFICIENTS(solver_row_idx)
                                VALUE=RHS_VALUE*row_coupling_coefficient
                                CALL DISTRIBUTED_VECTOR_VALUES_ADD(SOLVER_RHS_VECTOR,solver_row_number,VALUE,ERR,ERROR,*999)
                              ENDDO !solver_row_idx                          
                              IF(ASSOCIATED(SOURCE_MAPPINGS)) THEN
!!TODO: source term
                              ENDIF
                            ENDDO !equations_row_number
                          ELSE
                            CALL FLAG_ERROR("Fixed conditions are not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations equations matrices is not associated.",ERR,ERROR,*999)
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
              !Update the solver RHS vector values
              CALL DISTRIBUTED_VECTOR_UPDATE_START(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
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
        ELSE
          CALL FLAG_ERROR("Solver solver matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver matrices solution mapping is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_MATRICES_ASSEMBLE")
    RETURN
999 CALL ERRORS("SOLVER_MATRICES_ASSEMBLE",ERR,ERROR)
    CALL EXITS("SOLVER_MATRICES_ASSEMBLE")
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_ASSEMBLE

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

    CALL ENTERS("SOLVER_NONLINEAR_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
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

  !>Initialise a nonlinear solver for a problem solver
  SUBROUTINE SOLVER_NONLINEAR_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_NONLINEAR_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%NONLINEAR_SOLVER)) THEN
        CALL FLAG_ERROR("Nonlinear solver is already associated for this solver.",ERR,ERROR,*999)
      ELSE
        SOLUTION=>SOLVER%SOLUTION
        IF(ASSOCIATED(SOLUTION)) THEN
          SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
          IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
            IF(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES==1) THEN
              ALLOCATE(SOLVER%NONLINEAR_SOLVER,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver nonlinear solver.",ERR,ERROR,*999)
              SOLVER%NONLINEAR_SOLVER%SOLVER=>SOLVER
              SOLVER%NONLINEAR_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
            ELSE
              LOCAL_ERROR="The number of solver matrices in the solution mapping of "// &
                & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
                & " is invalid for a nonlinear solver. There should only be one solver matrix."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solution solution mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_NONLINEAR_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_NONLINEAR_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_NONLINEAR_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_NONLINEAR_INITIALISE

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

    CALL ENTERS("SOLVER_NONLINEAR_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
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

  !>Sets/changes the output type for a solver
  SUBROUTINE SOLVER_OUTPUT_TYPE_SET(SOLVER,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The type of solver output to be set \see SOLVER_ROUTINES_OutputTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE        
        SELECT CASE(OUTPUT_TYPE)
        CASE(SOLVER_NO_OUTPUT)
          SOLVER%OUTPUT_TYPE=SOLVER_NO_OUTPUT
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

  !>Sets/changes the sparsity type for a solver
  SUBROUTINE SOLVER_SPARSITY_TYPE_SET(SOLVER,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the problem solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The type of solver sparsity to be set \see SOLVER_ROUTINES_SparsityTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has been finished.",ERR,ERROR,*999)
      ELSE
!!TODO: Maybe set the sparsity in the different types of solvers. e.g., a sparse integrator doesn't mean much.
        SELECT CASE(SPARSITY_TYPE)
        CASE(SOLVER_SPARSE_MATRICES)
          SOLVER%SPARSITY_TYPE=SOLVER_SPARSE_MATRICES
        CASE(SOLVER_FULL_MATRICES)
          SOLVER%SPARSITY_TYPE=SOLVER_FULL_MATRICES
        CASE DEFAULT
          LOCAL_ERROR="The specified solver sparsity type of "// &
            & TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("SOLVER_SPARSITY_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_SPARSITY_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_SPARSITY_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_SPARSITY_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a time integration solver 
  SUBROUTINE SOLVER_TIME_INTEGRATION_CREATE_FINISH(TIME_INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(TIME_INTEGRATION_SOLVER_TYPE), POINTER :: TIME_INTEGRATION_SOLVER !<A pointer to the time integration solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_TIME_INTEGRATION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_INTEGRATION_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Time integration solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_TIME_INTEGRATION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("SOLVER_TIME_INTEGRATION_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_TIME_INTEGRATION_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_TIME_INTEGRATION_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Finalise a time integration solver for a problem solver
  SUBROUTINE SOLVER_TIME_INTEGRATION_FINALISE(TIME_INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(TIME_INTEGRATION_SOLVER_TYPE), POINTER :: TIME_INTEGRATION_SOLVER !<A pointer the time integration solver to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_TIME_INTEGRATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_INTEGRATION_SOLVER)) THEN        
      DEALLOCATE(TIME_INTEGRATION_SOLVER)
    ENDIF
        
    CALL EXITS("SOLVER_TIME_INTEGRATION_FINALISE")
    RETURN
999 CALL ERRORS("SOLVER_TIME_INTEGRATION_FINALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_TIME_INTEGRATION_FINALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_TIME_INTEGRATION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialise a time integration solver for a problem solver
  SUBROUTINE SOLVER_TIME_INTEGRATION_INITIALISE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to initialise the time integration solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_TIME_INTEGRATION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(ASSOCIATED(SOLVER%TIME_INTEGRATION_SOLVER)) THEN
        CALL FLAG_ERROR("Time integration solver is already associated for this solver.",ERR,ERROR,*999)
      ELSE
        SOLUTION=>SOLVER%SOLUTION
        IF(ASSOCIATED(SOLUTION)) THEN
          SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
          IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
            IF(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES==1) THEN
              ALLOCATE(SOLVER%TIME_INTEGRATION_SOLVER,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate solver time integration solver.",ERR,ERROR,*999)
              SOLVER%TIME_INTEGRATION_SOLVER%SOLVER=>SOLVER
              SOLVER%TIME_INTEGRATION_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
            ELSE
              LOCAL_ERROR="The number of solver matrices in the solution mapping of "// &
                & TRIM(NUMBER_TO_VSTRING(SOLUTION_MAPPING%NUMBER_OF_SOLVER_MATRICES,"*",ERR,ERROR))// &
                & " is invalid for a time integration solver. There should only be one solver matrix."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Problem solution solution mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver problem solution is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_TIME_INTEGRATION_INITIALISE")
    RETURN
999 CALL ERRORS("SOLVER_TIME_INTEGRATION_INITIALISE",ERR,ERROR)    
    CALL EXITS("SOLVER_TIME_INTEGRATION_INITIALISE")
    RETURN 1
   
  END SUBROUTINE SOLVER_TIME_INTEGRATION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Solve a time integration solver 
  SUBROUTINE SOLVER_TIME_INTEGRATION_SOLVE(TIME_INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(TIME_INTEGRATION_SOLVER_TYPE), POINTER :: TIME_INTEGRATION_SOLVER !<A pointer to the time integration solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_TIME_INTEGRATION_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_INTEGRATION_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Time integration solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_TIME_INTEGRATION_SOLVE")
    RETURN
999 CALL ERRORS("SOLVER_TIME_INTEGRATION_SOLVE",ERR,ERROR)    
    CALL EXITS("SOLVER_TIME_INTEGRATION_SOLVE")
    RETURN 1
    
  END SUBROUTINE SOLVER_TIME_INTEGRATION_SOLVE
        
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
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        !Assemble the solver matrices
        CALL SOLVER_MATRICES_ASSEMBLE(SOLVER,ERR,ERROR,*999)
        !If required output the solver matrices
        IF(SOLVER%OUTPUT_TYPE>=SOLVER_TIMING_OUTPUT) THEN
          CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
          CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
          IF(SOLVER%OUTPUT_TYPE>=SOLVER_MATRIX_OUTPUT) THEN
            SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
            IF(ASSOCIATED(SOLVER_MATRICES)) THEN
              CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRICES,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Solver solver matrices are not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ENDIF
        !Solve the system depending on the solver type
        SELECT CASE(SOLVER%SOLVE_TYPE)
        CASE(SOLVER_LINEAR_TYPE)
          CALL SOLVER_LINEAR_SOLVE(SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_NONLINEAR_TYPE)
          CALL SOLVER_NONLINEAR_SOLVE(SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
        CASE(SOLVER_TIME_INTEGRATION_TYPE)
          CALL SOLVER_TIME_INTEGRATION_SOLVE(SOLVER%TIME_INTEGRATION_SOLVER,ERR,ERROR,*999)
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

  !>Updates the dependent variables from the solver solution
  SUBROUTINE SOLVER_VARIABLES_UPDATE(SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer the solver to update the variables from
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: solver_dof_idx,solver_matrix_idx
    REAL(DP), POINTER :: SOLVER_DATA(:)
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOLVER_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX

    CALL ENTERS("SOLVER_VARIABLES_UPDATE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      IF(SOLVER%SOLVER_FINISHED) THEN
        SOLVER_MATRICES=>SOLVER%SOLVER_MATRICES
        IF(ASSOCIATED(SOLVER_MATRICES)) THEN
          SOLUTION_MAPPING=>SOLVER_MATRICES%SOLUTION_MAPPING
          IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
            DO solver_matrix_idx=1,SOLVER_MATRICES%NUMBER_OF_MATRICES
              SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
              IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                SOLVER_VECTOR=>SOLVER_MATRIX%SOLVER_VECTOR
                IF(ASSOCIATED(SOLVER_VECTOR)) THEN
                  DOMAIN_MAPPING=>SOLVER_VECTOR%DOMAIN_MAPPING
                  IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_GET(SOLVER_VECTOR,SOLVER_DATA,ERR,ERROR,*999)
                    DO solver_dof_idx=1,DOMAIN_MAPPING%NUMBER_OF_LOCAL
                    ENDDO !solver_dof_idx
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(SOLVER_VECTOR,SOLVER_DATA,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Solver vector domain mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver matrix is not associated.",ERR,ERROR,*999)
              ENDIF
            ENDDO !solver_matrix_idx
          ELSE
            CALL FLAG_ERROR("Solver matrices solution mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver solver matrices are not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_VARIABLES_UPDATE")
    RETURN
999 CALL ERRORS("SOLVER_VARIABLES_UPDATE",ERR,ERROR)    
    CALL EXITS("SOLVER_VARIABLES_UPDATE")
    RETURN 1
   
  END SUBROUTINE SOLVER_VARIABLES_UPDATE

  !
  !================================================================================================================================
  !
        
END MODULE SOLVER_ROUTINES
