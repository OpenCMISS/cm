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
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE STRINGS
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
  !>@}

  !> \addtogroup SOLVER_ROUTINES_SolverLibraries SOLVER_ROUTINES::SolverLibraries
  !> \brief The types of solver libraries
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_CMISS_LIBRARY=1 !<CMISS (internal) solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_PETSC_LIBRARY=2 !<PETSc solver library \see SOLVER_ROUTINES_SolverLibraries,SOLVER_ROUTINES
  !>@}

  !> \addtogroup SOLVER_ROUTINES_LinearSolverTypes SOLVER_ROUTINES::LinearSolverTypes
  !> \brief The types of linear solvers
  !> \see SOLVER_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_DIRECT_SOLVE=1 !<Direct linear solver type \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
  INTEGER(INTG), PARAMETER :: SOLVER_LINEAR_ITERATIVE_SOLVE=2 !<Iterative linear solver type \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
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

  PUBLIC SOLVER_LINEAR_DIRECT_SOLVE,SOLVER_LINEAR_ITERATIVE_SOLVE

  PUBLIC SOLVER_DIRECT_LU,SOLVER_DIRECT_CHOLESKY,SOLVER_DIRECT_SVD

  PUBLIC SOLVER_ITERATIVE_RICHARDSON,SOLVER_ITERATIVE_CHEBYCHEV,SOLVER_ITERATIVE_CONJUGATE_GRADIENT, &
    & SOLVER_ITERATIVE_BICONJUGATE_GRADIENT,SOLVER_ITERATIVE_GMRES,SOLVER_ITERATIVE_BiCGSTAB,SOLVER_ITERATIVE_CONJGRAD_SQUARED

  PUBLIC SOLVER_ITERATIVE_NO_PRECONDITIONER,SOLVER_ITERATIVE_JACOBI_PRECONDITIONER,SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER, &
    & SOLVER_ITERATIVE_SOR_PRECONDITIONER,SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER, &
    & SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER,SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER

  PUBLIC SOLVER_NO_OUTPUT,SOLVER_TIMING_OUTPUT,SOLVER_SOLVER_OUTPUT

  PUBLIC SOLVER_EULER_INTEGRATOR,SOLVER_IMPROVED_EULER_INTEGRATOR,SOLVER_4TH_RUNGE_KUTTA_INTEGRATOR, &
    & SOLVER_ADAMS_MOULTON_INTEGERATOR
  
  PUBLIC SOLVER_CREATE_FINISH,SOLVER_CREATE_START,SOLVER_DESTROY,SOLVER_TYPE_SET,SOLVER_SOLVE

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
  SUBROUTINE SOLVER_CREATE_FINISH(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("SOLVER_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Solver has already been finished",ERR,ERROR,*998)
      ELSE
        PROBLEM_SOLUTION=>PROBLEM_SOLVER%PROBLEM_SOLUTION
        IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
          SELECT CASE(PROBLEM_SOLVER%SOLVER_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_LINEAR_CREATE_FINISH(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_NONLINEAR_CREATE_FINISH(PROBLEM_SOLVER%NONLINEAR_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_TIME_INTEGRATION_TYPE)
            CALL SOLVER_TIME_INTEGRATION_CREATE_FINISH(PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SOLVER%SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          PROBLEM_SOLVER%SOLVER_FINISHED=.TRUE.
        ELSE
          CALL FLAG_ERROR("Problem solver problem solution is not associated",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_CREATE_FINISH")
    RETURN
999 CALL SOLVER_FINALISE(PROBLEM_SOLUTION,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("SOLVER_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE SOLVER_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a solver for a problem solution.
  SUBROUTINE SOLVER_CREATE_START(PROBLEM_SOLUTION,PROBLEM_SOLVER,ERR,ERROR,*)
    
    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer to the problem solution to create the solver for.
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<On exit, a pointer to the problem solver.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("SOLVER_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(PROBLEM_SOLUTION%SOLUTION_FINISHED) THEN
        CALL FLAG_ERROR("Problem solution has already been finished",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
          CALL FLAG_ERROR("Problem solver is already associated",ERR,ERROR,*998)
        ELSE
          NULLIFY(PROBLEM_SOLVER)
          IF(ASSOCIATED(PROBLEM_SOLUTION%PROBLEM)) THEN
            !Initialise the problem solver
            CALL SOLVER_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*999)
            PROBLEM_SOLVER=>PROBLEM_SOLUTION%SOLVER
          ELSE
            CALL FLAG_ERROR("Problem solution problem is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_CREATE_START")
    RETURN
999 CALL SOLVER_FINALISE(PROBLEM_SOLUTION,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("SOLVER_CREATE_START",ERR,ERROR)
    CALL EXITS("SOLVER_CREATE_START")
    RETURN 1   
  END SUBROUTINE SOLVER_CREATE_START    
  
  !
  !================================================================================================================================
  !

  !>Destroy a problem solver.
  SUBROUTINE SOLVER_DESTROY(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION

    CALL ENTERS("SOLVER_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      PROBLEM_SOLUTION=>PROBLEM_SOLVER%PROBLEM_SOLUTION
      IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
        CALL SOLVER_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Problem solver problem solution is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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

  !>Finalises a problem solver for a problem solution and deallocates all memory.
  SUBROUTINE SOLVER_FINALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer the problem solution to finalise the solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER)) THEN
        IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER%LINEAR_SOLVER)) &
          & CALL SOLVER_LINEAR_FINALISE(PROBLEM_SOLUTION%SOLVER,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER%NONLINEAR_SOLVER)) &
          & CALL SOLVER_NONLINEAR_FINALISE(PROBLEM_SOLUTION%SOLVER,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER%TIME_INTEGRATION_SOLVER)) &
          & CALL SOLVER_TIME_INTEGRATION_FINALISE(PROBLEM_SOLUTION%SOLVER,ERR,ERROR,*999)        
        DEALLOCATE(PROBLEM_SOLUTION%SOLVER)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solution is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_INITIALISE(PROBLEM_SOLUTION,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLUTION_TYPE), POINTER :: PROBLEM_SOLUTION !<A pointer the solution to initialise the solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLUTION)) THEN
      IF(ASSOCIATED(PROBLEM_SOLUTION%SOLVER)) THEN
        CALL FLAG_ERROR("Solver is already associated for this problems solution",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLUTION%SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solution solver",ERR,ERROR,*999)
        PROBLEM_SOLUTION%SOLVER%PROBLEM_SOLUTION=>PROBLEM_SOLUTION
        PROBLEM_SOLUTION%SOLVER%SOLVER_FINISHED=.FALSE.
        NULLIFY(PROBLEM_SOLUTION%SOLVER%LINEAR_SOLVER)
        NULLIFY(PROBLEM_SOLUTION%SOLVER%NONLINEAR_SOLVER)
        NULLIFY(PROBLEM_SOLUTION%SOLVER%TIME_INTEGRATION_SOLVER)
        PROBLEM_SOLUTION%SOLVER%SOLVER_TYPE=SOLVER_LINEAR_TYPE
        CALL SOLVER_LINEAR_INITIALISE(PROBLEM_SOLUTION%SOLVER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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
      CASE(SOLVER_LINEAR_DIRECT_SOLVE)
        CALL SOLVER_LINEAR_DIRECT_CREATE_FINISH(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE)
        CALL SOLVER_LINEAR_ITERATIVE_CREATE_FINISH(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))// &
          & " is invalid"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated",ERR,ERROR,*999)
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

    CALL ENTERS("SOLVER_LINEAR_DIRECT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_DIRECT_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Linear direct solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_DIRECT_FINALISE(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer the linear solver to finalise the direct solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_DIRECT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      IF(ASSOCIATED(LINEAR_SOLVER%DIRECT_SOLVER)) THEN        
        DEALLOCATE(LINEAR_SOLVER%DIRECT_SOLVER)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Direct solver is already associated for this linear solver",ERR,ERROR,*999)
      ELSE
        ALLOCATE(LINEAR_SOLVER%DIRECT_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear solver direct solver",ERR,ERROR,*999)
        LINEAR_SOLVER%DIRECT_SOLVER%LINEAR_SOLVER=>LINEAR_SOLVER
        LINEAR_SOLVER%DIRECT_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
        LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_LU
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated",ERR,ERROR,*999)
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

  !>Sets/changes the type of direct linear solver
  SUBROUTINE SOLVER_LINEAR_DIRECT_TYPE_SET(PROBLEM_SOLVER,DIRECT_SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set the direct linear solver type
    INTEGER(INTG), INTENT(IN) :: DIRECT_SOLVER_TYPE !<The type of direct linear solver to set \see SOLVER_ROUTINES_DirectLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_DIRECT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_DIRECT_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)) THEN
                IF(DIRECT_SOLVER_TYPE/=PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE) THEN
                  SELECT CASE(PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%SOLVER_LIBRARY)
                  CASE(SOLVER_PETSC_LIBRARY)
                    CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                    !Intialise the new solver type
                    !SELECT CASE(DIRECT_SOLVER_TYPE)
                    !CASE(SOLVER_DIRECT_LU)
                    !  PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_LU
                    !CASE(SOLVER_DIRECT_CHOLESKY)
                    !  PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_CHOLESKY
                    !CASE(SOLVER_DIRECT_SVD)
                    !  PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER%DIRECT_SOLVER_TYPE=SOLVER_DIRECT_SVD
                    !CASE DEFAULT
                    !  LOCAL_ERROR="The direct solver type of "//TRIM(NUMBER_TO_VSTRING(DIRECT_SOLVER_TYPE,"*",ERR,ERROR))// &
                    !    & " is invalid"
                    !  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    !END SELECT
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
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
  SUBROUTINE SOLVER_LINEAR_FINALISE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to finalise the linear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
        IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)) &
          & CALL SOLVER_LINEAR_DIRECT_FINALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
        IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) &
          & CALL SOLVER_LINEAR_ITERATIVE_FINALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
        DEALLOCATE(PROBLEM_SOLVER%LINEAR_SOLVER)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_INITIALISE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to initialise the linear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
        CALL FLAG_ERROR("Linear solver is already associated for this problems solver",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLVER%LINEAR_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solver linear solver",ERR,ERROR,*999)
        PROBLEM_SOLVER%LINEAR_SOLVER%PROBLEM_SOLVER=>PROBLEM_SOLVER
        NULLIFY(PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)
        NULLIFY(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)
        PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE=SOLVER_LINEAR_ITERATIVE_SOLVE
        CALL SOLVER_LINEAR_ITERATIVE_INITIALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET(PROBLEM_SOLVER,ABSOLUTE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set 
    REAL(DP), INTENT(IN) :: ABSOLUTE_TOLERANCE !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(ABSOLUTE_TOLERANCE>ZERO_TOLERANCE) THEN
                  PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE=ABSOLUTE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified absolute tolerance of "//TRIM(NUMBER_TO_VSTRING(ABSOLUTE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The absolute tolerance must be > 0"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver iterative solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear iterative solver",ERR,ERROR,*999)
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
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_ITERATIVE_SOLVER)) THEN
      SELECT CASE(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL PETSC_KSPCREATE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,LINEAR_ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
        CALL PETSC_KSPGETPC(LINEAR_ITERATIVE_SOLVER%KSP,LINEAR_ITERATIVE_SOLVER%PC,ERR,ERROR,*999)
        SELECT CASE(LINEAR_ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE)
        CASE(SOLVER_ITERATIVE_RICHARDSON)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPRICHARDSON_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_CHEBYCHEV)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPCHEBYCHEV_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPCG_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_BICONJUGATE_GRADIENT)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPBICG_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_GMRES)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPGMRES_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_BiCGSTAB)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPBCGS_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
          CALL PETSC_KSPSETTYPE(LINEAR_ITERATIVE_SOLVER%KSP,KSPCGS_,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The iterative solver type of "// &
            & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        SELECT CASE(LINEAR_ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE)
        CASE(SOLVER_ITERATIVE_NO_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCNONE_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_JACOBI_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCJACOBI_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCBJACOBI_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_SOR_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCSOR_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCICC_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCILU_,ERR,ERROR,*999)
        CASE(SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER)
          CALL PETSC_PCSETTYPE(LINEAR_ITERATIVE_SOLVER%PC,PCASM_,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The iterative preconditioner type of "// &
            & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        CALL PETSC_KSPSETTOLERANCES(LINEAR_ITERATIVE_SOLVER%KSP,LINEAR_ITERATIVE_SOLVER%RELATIVE_TOLERANCE, &
          & LINEAR_ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE,LINEAR_ITERATIVE_SOLVER%DIVERGENCE_TOLERANCE, &
          & LINEAR_ITERATIVE_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*999)
        CALL PETSC_KSPSETFROMOPTIONS(LINEAR_ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The solver library type of "// &
          & TRIM(NUMBER_TO_VSTRING(LINEAR_ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))//" is invalid"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT      
    ELSE
      CALL FLAG_ERROR("Linear itreative solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET(PROBLEM_SOLVER,DIVERGENCE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set 
    REAL(DP), INTENT(IN) :: DIVERGENCE_TOLERANCE !<The divergence tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_DIVERGENCE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(DIVERGENCE_TOLERANCE>ZERO_TOLERANCE) THEN
                  PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%DIVERGENCE_TOLERANCE=DIVERGENCE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified divergence tolerance of "// &
                    & TRIM(NUMBER_TO_VSTRING(DIVERGENCE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The divergence tolerance must be > 0"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver iterative solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear iterative solver",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_FINALISE(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer the linear solver to finalise the iterative solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      IF(ASSOCIATED(LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
        CALL PETSC_KSPDESTROY(LINEAR_SOLVER%ITERATIVE_SOLVER%KSP,ERR,ERROR,*999)
        DEALLOCATE(LINEAR_SOLVER%ITERATIVE_SOLVER)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Iterative solver is already associated for this linear solver",ERR,ERROR,*999)
      ELSE
        ALLOCATE(LINEAR_SOLVER%ITERATIVE_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate linear solver iterative solver",ERR,ERROR,*999)
        LINEAR_SOLVER%ITERATIVE_SOLVER%LINEAR_SOLVER=>LINEAR_SOLVER
        LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
        LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_GMRES
        LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE=SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER
        LINEAR_SOLVER%ITERATIVE_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS=100000
        LINEAR_SOLVER%ITERATIVE_SOLVER%RELATIVE_TOLERANCE=1.0E-05_DP
        LINEAR_SOLVER%ITERATIVE_SOLVER%ABSOLUTE_TOLERANCE=1.0E-10_DP
        LINEAR_SOLVER%ITERATIVE_SOLVER%DIVERGENCE_TOLERANCE=1.0E5_DP
      ENDIF
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET(PROBLEM_SOLVER,MAXIMUM_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set the maximum number of iterations
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_ITERATIONS !<The maximum number of iterations
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(MAXIMUM_ITERATIONS>0) THEN
                  PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
                ELSE
                  LOCAL_ERROR="The specified maximum iterations of "//TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                    & " is invalid. The maximum number of iterations must be > 0"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver iterative solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear iterative solver",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET(PROBLEM_SOLVER,ITERATIVE_PRECONDITIONER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: ITERATIVE_PRECONDITIONER_TYPE !<The type of iterative preconditioner to set \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_PRECONDITIONER_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(ITERATIVE_PRECONDITIONER_TYPE/=PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE) THEN
                  !Intialise the new preconditioner type
                  SELECT CASE(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY)
                  CASE(SOLVER_PETSC_LIBRARY)
                    SELECT CASE(ITERATIVE_PRECONDITIONER_TYPE)
                    CASE(SOLVER_ITERATIVE_NO_PRECONDITIONER)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_NO_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_JACOBI_PRECONDITIONER)
                      CALL FLAG_ERROR("Iterative Jacobi preconditioning is not implemented for a PETSc library",ERR,ERROR,*999)
                    CASE(SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_SOR_PRECONDITIONER)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_SOR_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER 
                    CASE(SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER
                    CASE(SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_PRECONDITIONER_TYPE= &
                        & SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER
                   CASE DEFAULT
                      LOCAL_ERROR="The iterative preconditioner type of "// &
                        & TRIM(NUMBER_TO_VSTRING(ITERATIVE_PRECONDITIONER_TYPE,"*",ERR,ERROR))//" is invalid"
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
                      & " is invalid"
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT                  
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver iterative solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear iterative solver",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET(PROBLEM_SOLVER,RELATIVE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set 
    REAL(DP), INTENT(IN) :: RELATIVE_TOLERANCE !<The relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(RELATIVE_TOLERANCE>ZERO_TOLERANCE) THEN
                  PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%RELATIVE_TOLERANCE=RELATIVE_TOLERANCE
                ELSE
                  LOCAL_ERROR="The specified relative tolerance of "//TRIM(NUMBER_TO_VSTRING(RELATIVE_TOLERANCE,"*",ERR,ERROR))// &
                    & " is invalid. The relative tolerance must be > 0"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver iterative solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear iterative solver",ERR,ERROR,*999)
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
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_RELATIVE_TOLERANCE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of iterative linear solver
  SUBROUTINE SOLVER_LINEAR_ITERATIVE_TYPE_SET(PROBLEM_SOLVER,ITERATIVE_SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: ITERATIVE_SOLVER_TYPE !<The type of iterative linear solver to set \see SOLVER_ROUTINES_IterativeLinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_ITERATIVE_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*999)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE==SOLVER_LINEAR_ITERATIVE_SOLVE) THEN
              IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) THEN
                IF(ITERATIVE_SOLVER_TYPE/=PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE) THEN
                  !Intialise the new solver type
                  SELECT CASE(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY)
                  CASE(SOLVER_PETSC_LIBRARY)
                    SELECT CASE(ITERATIVE_SOLVER_TYPE)
                    CASE(SOLVER_ITERATIVE_RICHARDSON)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_RICHARDSON
                    CASE(SOLVER_ITERATIVE_CHEBYCHEV)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_CHEBYCHEV
                    CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_CONJUGATE_GRADIENT
                    CASE(SOLVER_ITERATIVE_BICONJUGATE_GRADIENT)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_BICONJUGATE_GRADIENT
                    CASE(SOLVER_ITERATIVE_GMRES)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_BiCGSTAB
                    CASE(SOLVER_ITERATIVE_BiCGSTAB)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_BiCGSTAB
                    CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
                      PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%ITERATIVE_SOLVER_TYPE=SOLVER_ITERATIVE_CONJGRAD_SQUARED
                    CASE DEFAULT
                      LOCAL_ERROR="The iterative solver type of "//TRIM(NUMBER_TO_VSTRING(ITERATIVE_SOLVER_TYPE,"*",ERR,ERROR))// &
                        & " is invalid"
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  CASE DEFAULT
                    LOCAL_ERROR="The solver library type of "// &
                      & TRIM(NUMBER_TO_VSTRING(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER%SOLVER_LIBRARY,"*",ERR,ERROR))// &
                      & " is invalid"
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT                  
                ENDIF
              ELSE
                CALL FLAG_ERROR("The problem solver linear solver iterative solver is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The problem solver is not a linear iterative solver",ERR,ERROR,*999)
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
    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_TYPE_SET")
    RETURN
999 CALL ERRORS("SOLVER_LINEAR_ITERATIVE_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_ITERATIVE_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_ITERATIVE_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of linear solver
  SUBROUTINE SOLVER_LINEAR_TYPE_SET(PROBLEM_SOLVER,LINEAR_SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the problem solver to set the linear solver type
    INTEGER(INTG), INTENT(IN) :: LINEAR_SOLVER_TYPE !<The type of linear solver to set \see SOLVER_ROUTINES_LinearSolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("SOLVER_LINEAR_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*998)
      ELSE
        IF(PROBLEM_SOLVER%SOLVER_TYPE==SOLVER_LINEAR_TYPE) THEN
          IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) THEN
            IF(LINEAR_SOLVER_TYPE/=PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE) THEN
              !Intialise the new solver type
              SELECT CASE(LINEAR_SOLVER_TYPE)
              CASE(SOLVER_LINEAR_DIRECT_SOLVE)
                CALL SOLVER_LINEAR_DIRECT_INITIALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_LINEAR_ITERATIVE_SOLVE)
                CALL SOLVER_LINEAR_ITERATIVE_INITIALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The linear solver type of "//TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              !Finalise the old solver type
              SELECT CASE(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE)
              CASE(SOLVER_LINEAR_DIRECT_SOLVE)
                CALL SOLVER_LINEAR_DIRECT_FINALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE(SOLVER_LINEAR_ITERATIVE_SOLVE)
                CALL SOLVER_LINEAR_ITERATIVE_FINALISE(PROBLEM_SOLVER%LINEAR_SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The linear solver type of "// &
                  & TRIM(NUMBER_TO_VSTRING(PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              PROBLEM_SOLVER%LINEAR_SOLVER%LINEAR_SOLVER_TYPE=LINEAR_SOLVER_TYPE             
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem solver linear solver is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The problem solver is not a linear solver",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("SOLVER_LINEAR_TYPE_SET")
    RETURN
999 SELECT CASE(LINEAR_SOLVER_TYPE)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE)
      IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%DIRECT_SOLVER)) &
        & CALL SOLVER_LINEAR_DIRECT_FINALISE(PROBLEM_SOLVER%LINEAR_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE)
      IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER%ITERATIVE_SOLVER)) &
        & CALL SOLVER_LINEAR_ITERATIVE_FINALISE(PROBLEM_SOLVER%LINEAR_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    END SELECT
998 CALL ERRORS("SOLVER_LINEAR_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_LINEAR_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_LINEAR_TYPE_SET
        
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
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Nonlinear solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_NONLINEAR_FINALISE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to finalise the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_NONLINEAR_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(ASSOCIATED(PROBLEM_SOLVER%NONLINEAR_SOLVER)) THEN        
        DEALLOCATE(PROBLEM_SOLVER%NONLINEAR_SOLVER)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_NONLINEAR_INITIALISE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to initialise the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_NONLINEAR_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(ASSOCIATED(PROBLEM_SOLVER%NONLINEAR_SOLVER)) THEN
        CALL FLAG_ERROR("Nonlinear solver is already associated for this problems solver",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLVER%NONLINEAR_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solver nonlinear solver",ERR,ERROR,*999)
        PROBLEM_SOLVER%NONLINEAR_SOLVER%PROBLEM_SOLVER=>PROBLEM_SOLVER
        PROBLEM_SOLVER%NONLINEAR_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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

  !>Finishes the process of creating a time integration solver 
  SUBROUTINE SOLVER_TIME_INTEGRATION_CREATE_FINISH(TIME_INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(TIME_INTEGRATION_SOLVER_TYPE), POINTER :: TIME_INTEGRATION_SOLVER !<A pointer to the time integration solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_TIME_INTEGRATION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_INTEGRATION_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Time integration solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_TIME_INTEGRATION_FINALISE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to finalise the time integration solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_TIME_INTEGRATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(ASSOCIATED(PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER)) THEN        
        DEALLOCATE(PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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
  SUBROUTINE SOLVER_TIME_INTEGRATION_INITIALISE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to initialise the time integration solver for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_TIME_INTEGRATION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(ASSOCIATED(PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER)) THEN
        CALL FLAG_ERROR("Time integration solver is already associated for this problems solver",ERR,ERROR,*999)
      ELSE
        ALLOCATE(PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem solver time integration solver",ERR,ERROR,*999)
        PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER%PROBLEM_SOLVER=>PROBLEM_SOLVER
        PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER%SOLVER_LIBRARY=SOLVER_PETSC_LIBRARY
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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

  !>Sets/changes the type of solver
  SUBROUTINE SOLVER_TYPE_SET(PROBLEM_SOLVER,SOLVER_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to set the type of
    INTEGER(INTG), INTENT(IN) :: SOLVER_TYPE !<The type of solver to set \see SOLVER_ROUTINES_SolverTypes,SOLVER_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("SOLVER_TYPE_SET",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        CALL FLAG_ERROR("Problem solver has been finished",ERR,ERROR,*998)
      ELSE
        IF(SOLVER_TYPE/=PROBLEM_SOLVER%SOLVER_TYPE) THEN
          !Intialise the new solver type
          SELECT CASE(SOLVER_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_LINEAR_INITIALISE(PROBLEM_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_NONLINEAR_INITIALISE(PROBLEM_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_TIME_INTEGRATION_TYPE)
            CALL SOLVER_TIME_INTEGRATION_INITIALISE(PROBLEM_SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Finalise the old solver type
          SELECT CASE(PROBLEM_SOLVER%SOLVER_TYPE)
          CASE(SOLVER_LINEAR_TYPE)
            CALL SOLVER_LINEAR_FINALISE(PROBLEM_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL SOLVER_NONLINEAR_FINALISE(PROBLEM_SOLVER,ERR,ERROR,*999)
          CASE(SOLVER_TIME_INTEGRATION_TYPE)
            CALL SOLVER_TIME_INTEGRATION_FINALISE(PROBLEM_SOLVER,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The solver type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SOLVER%SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT          
          PROBLEM_SOLVER%SOLVER_TYPE=SOLVER_TYPE
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*998)
    ENDIF
        
    CALL EXITS("SOLVER_TYPE_SET")
    RETURN
999 SELECT CASE(SOLVER_TYPE)
    CASE(SOLVER_LINEAR_TYPE)
      IF(ASSOCIATED(PROBLEM_SOLVER%LINEAR_SOLVER)) CALL SOLVER_LINEAR_FINALISE(PROBLEM_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_NONLINEAR_TYPE)
      IF(ASSOCIATED(PROBLEM_SOLVER%NONLINEAR_SOLVER)) CALL SOLVER_NONLINEAR_FINALISE(PROBLEM_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    CASE(SOLVER_TIME_INTEGRATION_TYPE)
      IF(ASSOCIATED(PROBLEM_SOLVER%TIME_INTEGRATION_SOLVER)) &
        & CALL SOLVER_TIME_INTEGRATION_FINALISE(PROBLEM_SOLVER,DUMMY_ERR,DUMMY_ERROR,*998)
    END SELECT
998 CALL ERRORS("SOLVER_TYPE_SET",ERR,ERROR)    
    CALL EXITS("SOLVER_TYPE_SET")
    RETURN 1
   
  END SUBROUTINE SOLVER_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Solve the problem
  SUBROUTINE SOLVER_SOLVE(PROBLEM_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_SOLVER_TYPE), POINTER :: PROBLEM_SOLVER !<A pointer the solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM_SOLVER)) THEN
      IF(PROBLEM_SOLVER%SOLVER_FINISHED) THEN
        
      ELSE
        CALL FLAG_ERROR("Problem solver has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem solver is not associated",ERR,ERROR,*999)
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

  !>Solve the linear problem
  SUBROUTINE SOLVER_SOLVE_LINEAR(LINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_SOLVER_TYPE), POINTER :: LINEAR_SOLVER !<A pointer the linear solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("SOLVER_SOLVE_LINEAR",ERR,ERROR,*999)

    IF(ASSOCIATED(LINEAR_SOLVER)) THEN
      SELECT CASE(LINEAR_SOLVER%LINEAR_SOLVER_TYPE)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE)
        CALL SOLVER_SOLVE_LINEAR_DIRECT(LINEAR_SOLVER%DIRECT_SOLVER,ERR,ERROR,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE)
        CALL SOLVER_SOLVE_LINEAR_ITERATIVE(LINEAR_SOLVER%ITERATIVE_SOLVER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The linear solver type of "// &
          & TRIM(NUMBER_TO_VSTRING(LINEAR_SOLVER%LINEAR_SOLVER_TYPE,"*",ERR,ERROR))//" is invalid"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Linear solver is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_SOLVE_LINEAR")
    RETURN
999 CALL ERRORS("SOLVER_SOLVE_LINEAR",ERR,ERROR)    
    CALL EXITS("SOLVER_SOLVE_LINEAR")
    RETURN 1
   
  END SUBROUTINE SOLVER_SOLVE_LINEAR
        
  !
  !================================================================================================================================
  !

  !>Solve the direct linear problem
  SUBROUTINE SOLVER_SOLVE_LINEAR_DIRECT(DIRECT_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_DIRECT_SOLVER_TYPE), POINTER :: DIRECT_SOLVER !<A pointer the direct linear solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_SOLVE_LINEAR_DIRECT",ERR,ERROR,*999)

    IF(ASSOCIATED(DIRECT_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Direct linear solver is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_SOLVE_LINEAR_DIRECT")
    RETURN
999 CALL ERRORS("SOLVER_SOLVE_LINEAR_DIRECT",ERR,ERROR)    
    CALL EXITS("SOLVER_SOLVE_LINEAR_DIRECT")
    RETURN 1
   
  END SUBROUTINE SOLVER_SOLVE_LINEAR_DIRECT
        
  !
  !================================================================================================================================
  !

  !>Solve the iterative linear problem
  SUBROUTINE SOLVER_SOLVE_LINEAR_ITERATIVE(ITERATIVE_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(LINEAR_ITERATIVE_SOLVER_TYPE), POINTER :: ITERATIVE_SOLVER !<A pointer the direct iterative solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_SOLVE_LINEAR_ITERATIVE",ERR,ERROR,*999)

    IF(ASSOCIATED(ITERATIVE_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Direct iterative solver is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_SOLVE_LINEAR_ITERATIVE")
    RETURN
999 CALL ERRORS("SOLVER_SOLVE_LINEAR_ITERATIVE",ERR,ERROR)    
    CALL EXITS("SOLVER_SOLVE_LINEAR_ITERATIVE")
    RETURN 1
   
  END SUBROUTINE SOLVER_SOLVE_LINEAR_ITERATIVE
        
  !
  !================================================================================================================================
  !

  !>Solve the nonlinear problem
  SUBROUTINE SOLVER_SOLVE_NONLINEAR(NONLINEAR_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: NONLINEAR_SOLVER !<A pointer the nonlinear solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_SOLVE_NONLINEAR",ERR,ERROR,*999)

    IF(ASSOCIATED(NONLINEAR_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)      
    ELSE
      CALL FLAG_ERROR("Nonlinear solver is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_SOLVE_NONLINEAR")
    RETURN
999 CALL ERRORS("SOLVER_SOLVE_NONLINEAR",ERR,ERROR)    
    CALL EXITS("SOLVER_SOLVE_NONLINEAR")
    RETURN 1
   
  END SUBROUTINE SOLVER_SOLVE_NONLINEAR
        
  !
  !================================================================================================================================
  !

  !>Solve the time integration problem
  SUBROUTINE SOLVER_SOLVE_TIME_INTEGRATION(TIME_INTEGRATION_SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(TIME_INTEGRATION_SOLVER_TYPE), POINTER :: TIME_INTEGRATION_SOLVER !<A pointer the time integration solver to solve
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("SOLVER_SOLVE_TIME_INTEGRATION",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_INTEGRATION_SOLVER)) THEN
      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Time integration solver is not associated",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("SOLVER_SOLVE_TIME_INTEGRATION")
    RETURN
999 CALL ERRORS("SOLVER_SOLVE_TIME_INTEGRATION",ERR,ERROR)    
    CALL EXITS("SOLVER_SOLVE_TIME_INTEGRATION")
    RETURN 1
   
  END SUBROUTINE SOLVER_SOLVE_TIME_INTEGRATION
        
  !
  !================================================================================================================================
  !
        
END MODULE SOLVER_ROUTINES
