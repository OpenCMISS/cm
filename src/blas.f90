!> \file
!> $Id: blas.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module contains the interface descriptions to the BLAS routines.
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

!>This module contains the interface descriptions to the BLAS routines.
MODULE BLAS
  
  USE KINDS
  
  IMPLICIT NONE

  INTERFACE

    !Level 1 BLAS

    SUBROUTINE SASUM(N, X, INCX)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
    END SUBROUTINE SASUM
    
    SUBROUTINE DASUM(N, X, INCX)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
    END SUBROUTINE DASUM
    
    SUBROUTINE SAXPY(N, A, X, INCX, Y, INCY)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(IN) :: A
      REAL(SP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(SP), INTENT(OUT) :: Y(*)
      INTEGER(INTG), INTENT(IN) :: INCY
    END SUBROUTINE SAXPY

    SUBROUTINE DAXPY(N, A, X, INCX, Y, INCY)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(DP), INTENT(OUT) :: Y(*)
      INTEGER(INTG), INTENT(IN) :: INCY
    END SUBROUTINE DAXPY

    SUBROUTINE SCOPY(N, DX, INCX, DY, INCY)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(IN) :: DX(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(SP), INTENT(OUT) :: DY(*)
      INTEGER(INTG), INTENT(IN) :: INCY
    END SUBROUTINE SCOPY

    SUBROUTINE DCOPY(N, DX, INCX, DY, INCY)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: DX(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(DP), INTENT(OUT) :: DY(*)
      INTEGER(INTG), INTENT(IN) :: INCY
    END SUBROUTINE DCOPY

    FUNCTION SDOT(N, X, INCX, Y, INCY)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(SP), INTENT(IN) :: Y(*)
      INTEGER(INTG), INTENT(IN) :: INCY
      REAL(SP) :: SDOT
    END FUNCTION SDOT

    FUNCTION DDOT(N, X, INCX, Y, INCY)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(DP), INTENT(IN) :: Y(*)
      INTEGER(INTG), INTENT(IN) :: INCY
      REAL(DP) :: DDOT
    END FUNCTION DDOT

    FUNCTION SNRM2(N, X, INCX)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(SP) :: SNRM2
    END FUNCTION SNRM2

    FUNCTION DNRM2(N, X, INCX)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(DP) :: DNRM2
    END FUNCTION DNRM2

    SUBROUTINE SROT(N, DX, INCX, DY, INCY, C, S)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(OUT) :: DX(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(SP), INTENT(OUT) :: DY(*)
      INTEGER(INTG), INTENT(IN) :: INCY
      REAL(SP), INTENT(IN) :: C
      REAL(SP), INTENT(IN) :: S
    END SUBROUTINE SROT

    SUBROUTINE DROT(N, DX, INCX, DY, INCY, C, S)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(OUT) :: DX(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(DP), INTENT(OUT) :: DY(*)
      INTEGER(INTG), INTENT(IN) :: INCY
      REAL(DP), INTENT(IN) :: C
      REAL(DP), INTENT(IN) :: S
    END SUBROUTINE DROT

    SUBROUTINE SROTG(DA, DB, C, S)
      USE KINDS
      REAL(SP), INTENT(IN) :: DA
      REAL(SP), INTENT(IN) :: DB
      REAL(SP), INTENT(IN) :: C
      REAL(SP), INTENT(IN) :: S
    END SUBROUTINE SROTG

    SUBROUTINE DROTG(DA, DB, C, S)
      USE KINDS
      REAL(DP), INTENT(IN) :: DA
      REAL(DP), INTENT(IN) :: DB
      REAL(DP), INTENT(IN) :: C
      REAL(DP), INTENT(IN) :: S
    END SUBROUTINE DROTG

    SUBROUTINE SSCAL(N, A, X, INCX)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(SP), INTENT(IN) :: A
      REAL(SP), INTENT(INOUT) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
    END SUBROUTINE SSCAL

    SUBROUTINE DSCAL(N, A, X, INCX)
      USE KINDS
      INTEGER(INTG), INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: A
      REAL(DP), INTENT(INOUT) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
    END SUBROUTINE DSCAL

    !Level 2 BLAS
    
    SUBROUTINE SGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y ,&
      & INCY)
      USE KINDS
      CHARACTER(LEN=1), INTENT(IN) :: TRANS
      INTEGER(INTG), INTENT(IN) :: M, N
      REAL(SP), INTENT(IN) :: ALPHA
      INTEGER(INTG), INTENT(IN) :: LDA
      REAL(SP), INTENT(IN) :: A(LDA,*)
      REAL(SP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(SP), INTENT(IN) :: BETA
      REAL(SP), INTENT(INOUT) :: Y(*)
      INTEGER(INTG), INTENT(IN) :: INCY      
    END SUBROUTINE SGEMV

    SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y ,&
      & INCY)
      USE KINDS
      CHARACTER(LEN=1), INTENT(IN) :: TRANS
      INTEGER(INTG), INTENT(IN) :: M, N
      REAL(DP), INTENT(IN) :: ALPHA
      INTEGER(INTG), INTENT(IN) :: LDA
      REAL(DP), INTENT(IN) :: A(LDA,*)
      REAL(DP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
      REAL(DP), INTENT(IN) :: BETA
      REAL(DP), INTENT(INOUT) :: Y(*)
      INTEGER(INTG), INTENT(IN) :: INCY      
    END SUBROUTINE DGEMV

    SUBROUTINE STRSV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
      USE KINDS
      CHARACTER(LEN=1), INTENT(IN) :: UPLO, TRANS, DIAG
      INTEGER(INTG), INTENT(IN) :: N, LDA
      REAL(SP), INTENT(IN) :: A(LDA, *)
      REAL(SP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
    END SUBROUTINE STRSV

    SUBROUTINE DTRSV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
      USE KINDS
      CHARACTER(LEN=1), INTENT(IN) :: UPLO, TRANS, DIAG
      INTEGER(INTG), INTENT(IN) :: N, LDA
      REAL(DP), INTENT(IN) :: A(LDA, *)
      REAL(DP), INTENT(IN) :: X(*)
      INTEGER(INTG), INTENT(IN) :: INCX
    END SUBROUTINE DTRSV

    !Level 3 BLAS
    
  END INTERFACE

END MODULE BLAS
