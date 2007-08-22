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
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Solver libraries
  INTEGER(INTG), PARAMETER :: SOLVER_PETSC_LIBRARY=1

  !Linear solution methods
  INTEGER(INTG), PARAMETER :: SOLVER_DIRECT_SOLVE=1
  INTEGER(INTG), PARAMETER :: SOLVER_ITERATIVE_SOLVE=2

  !Direct solver types
  INTEGER(INTG), PARAMETER :: SOLVER_LU_DECOMPOSITION=1
  INTEGER(INTG), PARAMETER :: SOLVER_CHOLESKY_DECOMPOSITION=2
  INTEGER(INTG), PARAMETER :: SOLVER_SVD_DECOMPOSITION=3
  
  !Iterative method types
  INTEGER(INTG), PARAMETER :: SOLVER_RICHARDSON_ITERATIVE=1
  INTEGER(INTG), PARAMETER :: SOLVER_CHEBYCHEV_ITERATIVE=2
  INTEGER(INTG), PARAMETER :: SOLVER_CONJUGATE_GRADIENT_ITERATIVE=3
  INTEGER(INTG), PARAMETER :: SOLVER_BICONJUGATE_GRADIENT_ITERATIVE=4
  INTEGER(INTG), PARAMETER :: SOLVER_GMRES_ITERATIVE=5
  INTEGER(INTG), PARAMETER :: SOLVER_BiCGSTAB_ITERATIVE=6
  INTEGER(INTG), PARAMETER :: SOLVER_CONJGRAD_SQUARED_ITERATIVE=7

  !Preconditioner types
  INTEGER(INTG), PARAMETER :: SOLVER_NO_PRECONDITIONER=0
  INTEGER(INTG), PARAMETER :: SOLVER_JACOBI_PRECONDITIONER=1
  INTEGER(INTG), PARAMETER :: SOLVER_BLOCK_JACOBI_PRECONDITIONER=2
  INTEGER(INTG), PARAMETER :: SOLVER_SOR_PRECONDITIONER=3
  INTEGER(INTG), PARAMETER :: SOLVER_INCOMPLETE_CHOLESKY_PRECONDITIONER=4
  INTEGER(INTG), PARAMETER :: SOLVER_INCOMPLETE_LU_PRECONDITIONER=5

  !Time output parameters
  INTEGER(INTG), PARAMETER :: SOLVER_NO_OUTPUT=0
  INTEGER(INTG), PARAMETER :: SOLVER_TIMING_OUTPUT=1
  INTEGER(INTG), PARAMETER :: SOLVER_SOLVER_OUTPUT=2
  INTEGER(INTG), PARAMETER :: SOLVER_GLOBAL_SOLUTION_MATRIX=3
  INTEGER(INTG), PARAMETER :: SOLVER_GLOBAL_STIFFNESS_MATRIX=4
  INTEGER(INTG), PARAMETER :: SOLVER_ELEMENT_MATRICES=5

  !Integration procedures
  INTEGER(INTG), PARAMETER :: SOLVER_EULER_INTEGRATOR=1
  INTEGER(INTG), PARAMETER :: SOLVER_IMPROVED_EULER_INTEGRATOR=2
  INTEGER(INTG), PARAMETER :: SOLVER_4TH_RUNGE_KUTTA_INTEGRATOR=3
  INTEGER(INTG), PARAMETER :: SOLVER_ADAMS_MOULTON_INTEGERATOR=4
  INTEGER(INTG), PARAMETER :: SOLVER_LSODA_INTEGRATOR=5
  
  !Module types

  !Module variables

  !Interfaces

!CONTAINS

  !
  !================================================================================================================================
  !

  !
  !================================================================================================================================
  !

END MODULE SOLVER_ROUTINES
