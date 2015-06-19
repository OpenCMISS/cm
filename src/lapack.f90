!> \file
!> \author Chris Bradley
!> \brief This module contains the interface descriptions to the LAPACK routines.
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

!>This module contains the interface descriptions to the LAPACK routines.
MODULE LAPACK
  
  USE KINDS
  
  IMPLICIT NONE

  INTERFACE

    SUBROUTINE DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      INTEGER(INTG), INTENT(IN) :: NRHS
      INTEGER(INTG), INTENT(IN) :: LDA
      REAL(DP), INTENT(INOUT) :: A(LDA,*)
      INTEGER(INTG), INTENT(OUT) :: IPIV(*)
      INTEGER(INTG), INTENT(IN) :: LDB
      REAL(DP), INTENT(INOUT) :: B(LDB,*)
      INTEGER(INTG), INTENT(OUT) :: INFO
    END SUBROUTINE DGESV

    ! DGESVD - compute the singular value decomposition (SVD) of a
    !  real M-by-N matrix A, optionally computing the left and/or
    !  right singular vectors
    SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
      USE KINDS
      CHARACTER(1) :: JOBU ! Specifies options for computing all or part of the matrix U (options: A,S,O,N)
      CHARACTER(1) :: JOBVT ! Specifies options for computing all or part of the matrix V**T (options: A,S,O,N)
      INTEGER(INTG), INTENT(IN) :: M ! Number of rows in A
      INTEGER(INTG), INTENT(IN) :: N ! Number of columns in A
      REAL(DP), INTENT(INOUT) :: A(LDA,*) ! The matrix to perform the SVD on
      INTEGER(INTG), INTENT(IN) :: LDA ! Leading dimension of A
      REAL(DP), INTENT(OUT) :: S(MIN(M,N)) ! Singular values of A, sorted S(i) >= S(i+1)
      REAL(DP), INTENT(OUT) :: U(LDU,*) ! If JOBU = 'A', U contains the M-by-M orthogonal matrix U
      INTEGER(INTG), INTENT(IN) :: LDU ! Leading dimension of U
      REAL(DP), INTENT(OUT) :: VT(LDVT,N) ! If JOBVT = 'A', VT contains the N-by-N orthogonal matrix V**T
      INTEGER(INTG), INTENT(IN) :: LDVT ! The leading dimension of the array VT
      REAL(DP), INTENT(INOUT) :: WORK(*) ! On exit, if INFO = 0, WORK(1) returns the optimal LWORK
      INTEGER(INTG), INTENT(IN) :: LWORK ! The dimension of the array WORK
      INTEGER(INTG), INTENT(OUT) :: INFO ! 0 if successful exit; < 0 if INFO = -i (the i-th argument had an illegal value); > 0 if DBDSQR did not converge
    END SUBROUTINE DGESVD
    
  END INTERFACE

END MODULE LAPACK
