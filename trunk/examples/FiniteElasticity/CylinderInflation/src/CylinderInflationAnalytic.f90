!> \file
!> $Id: CylinderInflationExample.f90  $
!> \author Jack Lee
!> \brief This is an example program to solve a finite elasticity equation using openCMISS calls.
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
!> Contributor(s): Jack Lee
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


MODULE CYLINDERINFLATIONANALYTIC
  USE CONSTANTS
  USE KINDS  
  USE TYPES

  IMPLICIT NONE

  !  MODULE PROCEDURES
  PUBLIC :: CYLINDER_INFLATION_SOLVE
!   PRIVATE :: 

  CONTAINS
  
  ! =============================================================================
  SUBROUTINE CYLINDER_INFLATION_SOLVE(R_OUTER,R_INNER,C1,C2,P_INNER,P_OUTER,LAMBDA,MU1,MU2)
    REAL(DP),INTENT(IN) :: R_OUTER   !<initial outer radius
    REAL(DP),INTENT(IN) :: R_INNER   !<initial inner radius
    REAL(DP),INTENT(IN) :: C1        !<Mooney-Rivlin parameter for c_1*(I_1-3) term
    REAL(DP),INTENT(IN) :: C2        !<Mooney-Rivlin parameter for c_2*(I_2-3) term
    REAL(DP),INTENT(IN) :: P_INNER   !<pressure on the inner surface
    REAL(DP),INTENT(IN) :: P_OUTER   !<pressure on the outer surface
    REAL(DP),INTENT(IN) :: LAMBDA    !<stretch ratio in axial direction
    REAL(DP),INTENT(OUT) :: MU1      !<R_outer=MU1*R1
    REAL(DP),INTENT(OUT) :: MU2      !<R_inner=MU2*R2
    ! local variables
    REAL(DP) :: B,BP,DB,DMU,DMU1,DIFF
    REAL(DP) :: PHI ! MAKE THIS AN INPUT
    REAL(DP) :: TOL=1E-11_DP

    ! We're gonna solve for this by using the equations in Rivlin (1950) Phil. Trans A242, 173-95.
    ! The equations are nonlinear so will use Newton-Raphson method (finite difference jacobian,
    ! since I can't be arsed)
    write(*,*) "*** calculating analytic solutions ***"

    DMU = 1.0e-3_DP
    
PHI=0.0_DP  ! HACK

    ! calculate MU1 first
    MU1 = 1.0_DP ! initial estimate
    DO
      ! calculate residual
      CALL CYLINDER_INFLATION_CALC_RESIDUAL(R_OUTER,R_INNER,C1,C2,P_INNER,P_OUTER,LAMBDA,PHI,MU1,B)

      ! calculate perturbed residual
      CALL CYLINDER_INFLATION_CALC_RESIDUAL(R_OUTER,R_INNER,C1,C2,P_INNER,P_OUTER,LAMBDA,PHI,MU1+DMU,BP)

      ! calculate the new estimate for MU1
      DB = (BP-B)/DMU  ! gradient of function
      DMU1=-B/DB
      MU1 = MU1 + DMU1
      ! when do we exit?
      DIFF = ABS(DMU1)
!       WRITE(*,*) "relative norm = ",DIFF
      IF(DIFF<TOL) EXIT
      IF(ISNAN(DIFF)) EXIT      ! fuck
    ENDDO

    ! calculate MU2 now
    MU2 = SQRT((R_OUTER**2*(LAMBDA*MU1**2-1.0_DP)/R_INNER**2+1.0_DP)/LAMBDA)

    write(*,*) "*** finished calculating analytic solutions ***"

  END SUBROUTINE CYLINDER_INFLATION_SOLVE

  ! =============================================================================
  SUBROUTINE CYLINDER_INFLATION_CALC_RESIDUAL(R1,R2,C1,C2,P_INNER,P_OUTER,LAMBDA,PHI,MU1,B)
    REAL(DP),INTENT(IN) :: R1   !<initial outer radius
    REAL(DP),INTENT(IN) :: R2   !<initial inner radius
    REAL(DP),INTENT(IN) :: C1   !<Mooney-Rivlin parameter for c_1*(I_1-3) term
    REAL(DP),INTENT(IN) :: C2   !<Mooney-Rivlin parameter for c_2*(I_2-3) terms
    REAL(DP),INTENT(IN) :: P_INNER   !<pressure on the inner surface
    REAL(DP),INTENT(IN) :: P_OUTER   !<pressure on the outer surface
    REAL(DP),INTENT(IN) :: LAMBDA    !<stretch ratio in axial direction
    REAL(DP),INTENT(IN) :: PHI       !<angle of torsion
    REAL(DP),INTENT(IN) :: MU1       !<current estimate of MU1
    REAL(DP),INTENT(OUT) :: B        !<on exit, contains the residual
    ! local variables
    REAL(DP) :: R,K,MU,P,SIGMARR   ! SOME DUMMIES

    ! equation: Sigma(rr) + p_inner = 0 (=residual)
    R=R2  ! undeformed R for this problem = R2
    K=R1**2*(LAMBDA*MU1**2-1.0_DP)
    MU=SQRT(1.0_DP/LAMBDA*(1.0_DP+K/R**2))
    P=P_OUTER-(C1/LAMBDA+C2*LAMBDA)*(1.0_DP/LAMBDA/MU1**2-R**2/(R**2+K) &
        & +2.0_DP*LOG(MU/MU1))+C1*PHI**2*LAMBDA*(R**2-R1**2) &
        & -2.0_DP*(C1/LAMBDA**2/MU**2+C2*(1.0_DP/LAMBDA**2+1.0_DP/MU**2+PHI**2*R**2))
    SIGMARR=2.0_dp*(C1/LAMBDA**2/MU**2+C2*(1.0_DP/LAMBDA**2+1.0_DP/MU**2+PHI**2*R**2))+P

    B=SIGMARR+P_INNER


  END SUBROUTINE CYLINDER_INFLATION_CALC_RESIDUAL

END MODULE CYLINDERINFLATIONANALYTIC