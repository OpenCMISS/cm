!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all mathematics support routines.
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
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s): Kumar Mithraratne
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

!> This module contains all mathematics support routines.
MODULE MATHS

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  INTERFACE CROSS_PRODUCT
    MODULE PROCEDURE CROSS_PRODUCT_INTG
    MODULE PROCEDURE CROSS_PRODUCT_SP
    MODULE PROCEDURE CROSS_PRODUCT_DP
  END INTERFACE !CROSS_PRODUCT

  INTERFACE D_CROSS_PRODUCT
    MODULE PROCEDURE D_CROSS_PRODUCT_INTG
    MODULE PROCEDURE D_CROSS_PRODUCT_SP
    MODULE PROCEDURE D_CROSS_PRODUCT_DP
  END INTERFACE !D_CROSS_PRODUCT

  INTERFACE DETERMINANT
    MODULE PROCEDURE DETERMINANT_FULL_INTG
    MODULE PROCEDURE DETERMINANT_FULL_SP
    MODULE PROCEDURE DETERMINANT_FULL_DP
  END INTERFACE !DETERMINANT
    
  INTERFACE EDP
    MODULE PROCEDURE EDP_DP
    MODULE PROCEDURE EDP_SP
  END INTERFACE !EDP

  INTERFACE EIGENVALUE
    MODULE PROCEDURE EIGENVALUE_FULL_SP
    MODULE PROCEDURE EIGENVALUE_FULL_DP
  END INTERFACE !EIGENVALUE

  INTERFACE EIGENVECTOR
    MODULE PROCEDURE EIGENVECTOR_FULL_SP
    MODULE PROCEDURE EIGENVECTOR_FULL_DP
  END INTERFACE !EIGENVECTOR

  INTERFACE I0
    MODULE PROCEDURE I0_DP
    MODULE PROCEDURE I0_SP
  END INTERFACE !I0

  INTERFACE I1
    MODULE PROCEDURE I1_DP
    MODULE PROCEDURE I1_SP
  END INTERFACE !I1

  INTERFACE INVERT
    MODULE PROCEDURE INVERT_FULL_SP
    MODULE PROCEDURE INVERT_FULL_DP
  END INTERFACE !INVERT

  INTERFACE K0
    MODULE PROCEDURE K0_DP
    MODULE PROCEDURE K0_SP
  END INTERFACE !K0

  INTERFACE K1
    MODULE PROCEDURE K1_DP
    MODULE PROCEDURE K1_SP
  END INTERFACE !K1

  INTERFACE L2NORM
    MODULE PROCEDURE L2NORM_SP
    MODULE PROCEDURE L2NORM_DP
  END INTERFACE !L2NORM

  INTERFACE MATRIX_PRODUCT
    MODULE PROCEDURE MATRIX_PRODUCT_SP
    MODULE PROCEDURE MATRIX_PRODUCT_DP
  END INTERFACE !MATRIX_PRODUCT
  
  INTERFACE MATRIX_TRANSPOSE
    MODULE PROCEDURE MATRIX_TRANSPOSE_SP
    MODULE PROCEDURE MATRIX_TRANSPOSE_DP
  END INTERFACE !MATRIX_TRANSPOSE 

  INTERFACE NORMALISE
    MODULE PROCEDURE NORMALISE_SP
    MODULE PROCEDURE NORMALISE_DP
  END INTERFACE !NORMALISE

  INTERFACE SOLVE_SMALL_LINEAR_SYSTEM
    MODULE PROCEDURE SOLVE_SMALL_LINEAR_SYSTEM_SP
    MODULE PROCEDURE SOLVE_SMALL_LINEAR_SYSTEM_DP
  END INTERFACE !SOLVE_SMALL_LINEAR_SYSTEM

  PUBLIC CROSS_PRODUCT,D_CROSS_PRODUCT,DETERMINANT,EIGENVALUE,EIGENVECTOR,INVERT,L2NORM,MATRIX_PRODUCT,MATRIX_TRANSPOSE, &
    & NORMALISE,SOLVE_SMALL_LINEAR_SYSTEM
  
  
CONTAINS

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: CROSS_PRODUCT
  !###  Description:
  !###    Calculates the vector cross product of two vectors
  !###  Child-subroutines: CROSS_PRODUCT_INTG,CROSS_PRODUCT_SP,CROSS_PRODUCT_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE CROSS_PRODUCT_INTG(A,B,C,ERR,ERROR,*)
  
    !#### Subroutine: CROSS_PRODUCT_INTG
    !###  Description:
    !###    Calculates and returns the vector cross-prouct of the integer vectors A*B in C.
    !###  Parent-function: CROSS_PRODUCT
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:),B(:)
    INTEGER(INTG), INTENT(OUT) :: C(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    
    CALL ENTERS("CROSS_PRODUCT_INTG",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(B,1)) THEN
      IF(SIZE(C,1)==3) THEN
        SELECT CASE(SIZE(A,1))
        CASE(3)
          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid vector size",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("The vector C is not the correct size",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("The vectors A and B are not the same size",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CROSS_PRODUCT_INTG")
    RETURN
999 CALL ERRORS("CROSS_PRODUCT_INTG",ERR,ERROR)
    CALL EXITS("CROSS_PRODUCT_INTG")
    RETURN 1
  END SUBROUTINE CROSS_PRODUCT_INTG
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE CROSS_PRODUCT_SP(A,B,C,ERR,ERROR,*)
  
    !#### Subroutine: CROSS_PRODUCT_SP
    !###  Description:
    !###    Calculates and returns the vector cross-prouct of the single precision vectors A*B in C.
    !###  Parent-function: CROSS_PRODUCT
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:),B(:)
    REAL(SP), INTENT(OUT) :: C(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    
    CALL ENTERS("CROSS_PRODUCT_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(B,1)) THEN
      IF(SIZE(C,1)==3) THEN
        SELECT CASE(SIZE(A,1))
        CASE(3)
          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid vector size",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("The vector C is not the correct size",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("The vectors A and B are not the same size",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CROSS_PRODUCT_SP")
    RETURN
999 CALL ERRORS("CROSS_PRODUCT_SP",ERR,ERROR)
    CALL EXITS("CROSS_PRODUCT_SP")
    RETURN 1
  END SUBROUTINE CROSS_PRODUCT_SP
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE CROSS_PRODUCT_DP(A,B,C,ERR,ERROR,*)
  
    !#### Subroutine: CROSS_PRODUCT_DP
    !###  Description:
    !###    Calculates and returns the vector cross-prouct of the double precision vectors A*B in C.
    !###  Parent-function: CROSS_PRODUCT
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:),B(:)
    REAL(DP), INTENT(OUT) :: C(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    
    CALL ENTERS("CROSS_PRODUCT_DP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(B,1)) THEN
      IF(SIZE(C,1)==3) THEN
        SELECT CASE(SIZE(A,1))
        CASE(3)
          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        CASE DEFAULT
          CALL FLAG_ERROR("Invalid vector size",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("The vector C is not the correct size",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("The vectors A and B are not the same size",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("CROSS_PRODUCT_DP")
    RETURN
999 CALL ERRORS("CROSS_PRODUCT_DP",ERR,ERROR)
    CALL EXITS("CROSS_PRODUCT_DP")
    RETURN 1
  END SUBROUTINE CROSS_PRODUCT_DP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: D_CROSS_PRODUCT
  !###  Description:
  !###    Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
  !###    derivatives D_A and D_B of A and B
  !###  Child-subroutines: D_CROSS_PRODUCT_INTG,D_CROSS_PRODUCT_SP,D_CROSS_PRODUCT_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE D_CROSS_PRODUCT_INTG(N,A,B,C,D_A,D_B,D_C,ERR,ERROR,*)
  
    !#### Subroutine: D_CROSS_PRODUCT_INTG
    !###  Description:
    !###    Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
    !###    derivatives D_A and D_B of A and B for integer vectors.
    !###  Parent-function: D_CROSS_PRODUCT
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N,A(:),B(:)
    INTEGER(INTG), INTENT(OUT) :: C(:)
    INTEGER(INTG), INTENT(IN) :: D_A(:,:),D_B(:,:)
    INTEGER(INTG), INTENT(OUT) :: D_C(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: ni
    
    CALL ENTERS("D_CROSS_PRODUCT_INTG",ERR,ERROR,*999)

    CALL CROSS_PRODUCT(A,B,C,ERR,ERROR,*999)
    IF(SIZE(D_A,1)==SIZE(D_B,1).AND.SIZE(A,1)==SIZE(D_A,1).AND.SIZE(B,1)==SIZE(D_B,1)) THEN
      IF(SIZE(D_A,2)==N.AND.SIZE(D_B,2)==N) THEN
        IF(SIZE(C,1)==3) THEN
          SELECT CASE(SIZE(A,1))
          CASE(3)
            DO ni=1,N
              D_C(1,ni)=D_A(2,ni)*B(3)-D_A(3,ni)*B(2)+A(2)*D_B(3,ni)-A(3)*D_B(2,ni)
              D_C(2,ni)=D_A(3,ni)*B(1)-D_A(1,ni)*B(3)+A(3)*D_B(1,ni)-A(1)*D_B(3,ni)
              D_C(3,ni)=D_A(1,ni)*B(2)-D_A(2,ni)*B(1)+A(1)*D_B(2,ni)-A(2)*D_B(1,ni)
            ENDDO !ni
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid vector size",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("The vector C is not the correct size",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The number of derivative vectors is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("The vectors for D_A and D_B are not the same size",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("D_CROSS_PRODUCT_INTG")
    RETURN
999 CALL ERRORS("D_CROSS_PRODUCT_INTG",ERR,ERROR)
    CALL EXITS("D_CROSS_PRODUCT_INTG")
    RETURN 1
  END SUBROUTINE D_CROSS_PRODUCT_INTG
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE D_CROSS_PRODUCT_SP(N,A,B,C,D_A,D_B,D_C,ERR,ERROR,*)
  
    !#### Subroutine: D_CROSS_PRODUCT_SP
    !###  Description:
    !###    Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
    !###    derivatives D_A and D_B of A and B for single precision vectors.
    !###  Parent-function: D_CROSS_PRODUCT
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N
    REAL(SP), INTENT(IN) :: A(:),B(:)
    REAL(SP), INTENT(OUT) :: C(:)
    REAL(SP), INTENT(IN) :: D_A(:,:),D_B(:,:)
    REAL(SP), INTENT(OUT) :: D_C(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: ni
    
    CALL ENTERS("D_CROSS_PRODUCT_SP",ERR,ERROR,*999)

    CALL CROSS_PRODUCT(A,B,C,ERR,ERROR,*999)
    IF(SIZE(D_A,1)==SIZE(D_B,1).AND.SIZE(A,1)==SIZE(D_A,1).AND.SIZE(B,1)==SIZE(D_B,1)) THEN
      IF(SIZE(D_A,2)==N.AND.SIZE(D_B,2)==N) THEN
        IF(SIZE(C,1)==3) THEN
          SELECT CASE(SIZE(A,1))
          CASE(3)
            DO ni=1,N
              D_C(1,ni)=D_A(2,ni)*B(3)-D_A(3,ni)*B(2)+A(2)*D_B(3,ni)-A(3)*D_B(2,ni)
              D_C(2,ni)=D_A(3,ni)*B(1)-D_A(1,ni)*B(3)+A(3)*D_B(1,ni)-A(1)*D_B(3,ni)
              D_C(3,ni)=D_A(1,ni)*B(2)-D_A(2,ni)*B(1)+A(1)*D_B(2,ni)-A(2)*D_B(1,ni)
            ENDDO !ni
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid vector size",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("The vector C is not the correct size",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The number of derivative vectors is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("The vectors for D_A and D_B are not the same size",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("D_CROSS_PRODUCT_SP")
    RETURN
999 CALL ERRORS("D_CROSS_PRODUCT_SP",ERR,ERROR)
    CALL EXITS("D_CROSS_PRODUCT_SP")
    RETURN 1
  END SUBROUTINE D_CROSS_PRODUCT_SP
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE D_CROSS_PRODUCT_DP(N,A,B,C,D_A,D_B,D_C,ERR,ERROR,*)
  
    !#### Subroutine: D_CROSS_PRODUCT_DP
    !###  Description:
    !###    Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
    !###    derivatives D_A and D_B of A and B for double precision vectors.
    !###  Parent-function: D_CROSS_PRODUCT
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N
    REAL(DP), INTENT(IN) :: A(:),B(:)
    REAL(DP), INTENT(OUT) :: C(:)
    REAL(DP), INTENT(IN) :: D_A(:,:),D_B(:,:)
    REAL(DP), INTENT(OUT) :: D_C(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: ni
    
    CALL ENTERS("D_CROSS_PRODUCT_DP",ERR,ERROR,*999)

    CALL CROSS_PRODUCT(A,B,C,ERR,ERROR,*999)
    IF(SIZE(D_A,1)==SIZE(D_B,1).AND.SIZE(A,1)==SIZE(D_A,1).AND.SIZE(B,1)==SIZE(D_B,1)) THEN
      IF(SIZE(D_A,2)==N.AND.SIZE(D_B,2)==N) THEN
        IF(SIZE(C,1)==3) THEN
          SELECT CASE(SIZE(A,1))
          CASE(3)
            DO ni=1,N
              D_C(1,ni)=D_A(2,ni)*B(3)-D_A(3,ni)*B(2)+A(2)*D_B(3,ni)-A(3)*D_B(2,ni)
              D_C(2,ni)=D_A(3,ni)*B(1)-D_A(1,ni)*B(3)+A(3)*D_B(1,ni)-A(1)*D_B(3,ni)
              D_C(3,ni)=D_A(1,ni)*B(2)-D_A(2,ni)*B(1)+A(1)*D_B(2,ni)-A(2)*D_B(1,ni)
            ENDDO !ni
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid vector size",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("The vector C is not the correct size",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The number of derivative vectors is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("The vectors for D_A and D_B are not the same size",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("D_CROSS_PRODUCT_DP")
    RETURN
999 CALL ERRORS("D_CROSS_PRODUCT_DP",ERR,ERROR)
    CALL EXITS("D_CROSS_PRODUCT_DP")
    RETURN 1
  END SUBROUTINE D_CROSS_PRODUCT_DP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: DETERMINANT
  !###  Description:
  !###    Returns the determinant of a matrix
  !###  Child-functions: DETERMINANT_FULL_INTG,DETERMINANT_FULL_SP,DETERMINANT_FULL_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION DETERMINANT_FULL_INTG(A,ERR,ERROR)
  
    !#### Function: DETERMINANT_FULL_INTG
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Returns the determinant of a full single precision matrix A.
    !###  Parent-function: DETERMINANT
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    INTEGER(INTG) :: DETERMINANT_FULL_INTG
    
    CALL ENTERS("DETERMINANT_FULL_INTG",ERR,ERROR,*999)

    DETERMINANT_FULL_INTG=0
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        DETERMINANT_FULL_INTG=A(1,1)
      CASE(2)
        DETERMINANT_FULL_INTG=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      CASE(3)
        DETERMINANT_FULL_INTG=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)- &
          A(2,1)*A(1,2)*A(3,3)-A(3,1)*A(2,2)*A(1,3)
      CASE DEFAULT
        CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DETERMINANT_FULL_INTG")
    RETURN
999 CALL ERRORS("DETERMINANT_FULL_INTG",ERR,ERROR)
    CALL EXITS("DETERMINANT_FULL_INTG")
    RETURN 
  END FUNCTION DETERMINANT_FULL_INTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION DETERMINANT_FULL_SP(A,ERR,ERROR)
  
    !#### Function: DETERMINANT_FULL_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Returns the determinant of a full single precision matrix A.
    !###  Parent-function: DETERMINANT
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(SP) :: DETERMINANT_FULL_SP
    
    CALL ENTERS("DETERMINANT_FULL_SP",ERR,ERROR,*999)

    DETERMINANT_FULL_SP=0.0_SP
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        DETERMINANT_FULL_SP=A(1,1)
      CASE(2)
        DETERMINANT_FULL_SP=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      CASE(3)
        DETERMINANT_FULL_SP=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)- &
          A(2,1)*A(1,2)*A(3,3)-A(3,1)*A(2,2)*A(1,3)
      CASE DEFAULT
        CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DETERMINANT_FULL_SP")
    RETURN
999 CALL ERRORS("DETERMINANT_FULL_SP",ERR,ERROR)
    CALL EXITS("DETERMINANT_FULL_SP")
    RETURN 
  END FUNCTION DETERMINANT_FULL_SP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION DETERMINANT_FULL_DP(A,ERR,ERROR)
  
    !#### Function: DETERMINANT_FULL_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Returns the determinant of a full double precision matrix A
    !###  Parent-function: DETERMINANT
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: DETERMINANT_FULL_DP
    
    CALL ENTERS("DETERMINANT_FULL_DP",ERR,ERROR,*999)

    DETERMINANT_FULL_DP=0.0_DP

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        DETERMINANT_FULL_DP=A(1,1)
      CASE(2)
        DETERMINANT_FULL_DP=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      CASE(3)
        DETERMINANT_FULL_DP=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)- &
          A(2,1)*A(1,2)*A(3,3)-A(3,1)*A(2,2)*A(1,3)
      CASE DEFAULT
        CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DETERMINANT_FULL_DP")
    RETURN
999 CALL ERRORS("DETERMINANT_FULL_DP",ERR,ERROR)
    CALL EXITS("DETERMINANT_FULL_DP")
    RETURN    
  END FUNCTION DETERMINANT_FULL_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: EDP
  !###  Description:
  !###    Calculates the elliptic integral of the second kind - E(m).
  !###  Child-functions: EDP_DP,EDP_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION EDP_DP(X)
  
    !#### Function: EDP_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Calculates the elliptic integral of the second kind - E(m), for a double precision argument.
    !###  Parent-function: EDP
    
    !Argument variables
    REAL(DP), INTENT(IN) :: X
    !Function variable
    REAL(DP) :: EDP_DP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: A1=0.44325141463_DP
    REAL(DP), PARAMETER :: A2=0.06260601220_DP
    REAL(DP), PARAMETER :: A3=0.04757383546_DP
    REAL(DP), PARAMETER :: A4=0.01736506451_DP
    REAL(DP), PARAMETER :: B1=0.24998368310_DP
    REAL(DP), PARAMETER :: B2=0.09200180037_DP
    REAL(DP), PARAMETER :: B3=0.04069697526_DP
    REAL(DP), PARAMETER :: B4=0.00526449639_DP
    REAL(DP) :: TERM1,TERM2,X1
    
    X1=1.0_DP-X
    TERM1=1.0_DP+(A1+(A2+(A3+A4*X1)*X1)*X1)*X1
    TERM2=(B1+(B2+(B3+B4*X1)*X1)*X1)*X1
    EDP_DP=TERM1+TERM2*LOG(1.0_DP/X1)
    
    RETURN 
  END FUNCTION EDP_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION EDP_SP(X)
  
    !#### Function: EDP_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Calculates the elliptic integral of the second kind - E(m), for a single precision argument.
    !###  Parent-function: EDP
    
    !Argument variables
    REAL(SP), INTENT(IN) :: X
    !Function variable
    REAL(SP) :: EDP_SP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: A1=0.44325141463_SP
    REAL(SP), PARAMETER :: A2=0.06260601220_SP
    REAL(SP), PARAMETER :: A3=0.04757383546_SP
    REAL(SP), PARAMETER :: A4=0.01736506451_SP
    REAL(SP), PARAMETER :: B1=0.24998368310_SP
    REAL(SP), PARAMETER :: B2=0.09200180037_SP
    REAL(SP), PARAMETER :: B3=0.04069697526_SP
    REAL(SP), PARAMETER :: B4=0.00526449639_SP
    REAL(SP) :: TERM1,TERM2,X1
    
    X1=1.0_SP-X
    TERM1=1.0_SP+(A1+(A2+(A3+A4*X1)*X1)*X1)*X1
    TERM2=(B1+(B2+(B3+B4*X1)*X1)*X1)*X1
    EDP_SP=TERM1+TERM2*LOG(1.0_SP/X1)
        
    RETURN 
  END FUNCTION EDP_SP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: EIGENVALUE
  !###  Description:
  !###    Returns the eigenvalues of a matrix.
  !###  Child-subroutines: EIGENVALUE_FULL_SP,EIGNEVECTOR_FULL_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE EIGENVALUE_FULL_SP(A,EVALUES,ERR,ERROR,*)
  
    !#### Subroutine: EIGENVALUE_FULL_SP
    !###  Description:
    !###    Returns the eigenvalues of a full single precision matrix A.
    !###  Parent-subroutine: EIGENVALUE
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:)
    REAL(SP), INTENT(OUT) :: EVALUES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i
    REAL(SP) :: ANGLE,B2,B3,C1,C2,D,Q,Q3,R,RI1,RI2,RI3,RI4,RQ,TEMP,THETA
    
    CALL ENTERS("EIGENVALUE_FULL_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)<=SIZE(EVALUES,1)) THEN
        SELECT CASE(SIZE(A,1))
        CASE(1)
          EVALUES(1)=A(1,1)
        CASE(2)
          IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
            RI1=A(1,1)+A(2,2)
            RI2=A(1,1)*A(2,2)-A(1,2)**2
            B2=RI1/2.0_SP
            C1=RI1*RI1
            C2=4.0_SP*RI2
            IF(C2>C1) CALL FLAG_ERROR("Complex roots found in quadratic equation",ERR,ERROR,*999)
            B3=SQRT(C1-C2)/2.0_SP
            EVALUES(1)=B2+B3
            EVALUES(2)=B2-B3
          ELSE
            EVALUES(1)=A(1,1)
            EVALUES(2)=A(2,2)
          ENDIF
          IF(ABS(EVALUES(2))>ABS(EVALUES(1))) THEN
            TEMP=EVALUES(1)
            EVALUES(1)=EVALUES(2)
            EVALUES(2)=TEMP
          ENDIF
        CASE(3)
          RI1=A(1,1)+A(2,2)+A(3,3)
          RI2=A(1,1)*A(2,2)+A(2,2)*A(3,3)+A(3,3)*A(1,1)-(A(1,2)**2+A(2,3)**2+A(3,1)**2)
          RI3=DETERMINANT(A,ERR,ERROR)
          IF(ERR /=0) GOTO 999
          RI4=RI1/3.0_SP
          Q=RI4*RI4-RI2/3.0_SP   
          R=RI4*(RI4*RI4-RI2/2.0_SP)+RI3/2.0_SP
          Q3=Q*Q*Q
          D=R*R-Q3
          IF(ABS(D)>ZERO_TOLERANCE_SP) CALL FLAG_ERROR("Complex roots found in solution of cubic equation",ERR,ERROR,*999)
          RQ=SQRT(ABS(Q))
          IF(ABS(Q)<ZERO_TOLERANCE_SP) THEN
            THETA=0.0_SP
          ELSE
            THETA=ACOS(R/SQRT(ABS(Q3)))/3.0_SP
          ENDIF
          ANGLE=2.0_SP*REAL(PI,SP)/3.0_SP
          EVALUES(1)=2.0_SP*RQ*COS(THETA)+RI4
          EVALUES(2)=2.0_SP*RQ*COS(THETA+ANGLE)+RI4
          EVALUES(3)=2.0_SP*RQ*COS(THETA+2.0_SP*ANGLE)+RI4
          DO i=1,2
            IF(ABS(EVALUES(3))>ABS(EVALUES(i))) THEN
              TEMP=EVALUES(i)
              EVALUES(i)=EVALUES(3)
              EVALUES(3)=TEMP
            ENDIF
          ENDDO !i
        CASE DEFAULT
          CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Evalues is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EIGENVALUE_FULL_SP")
    RETURN
999 CALL ERRORS("EIGENVALUE_FULL_SP",ERR,ERROR)
    CALL EXITS("EIGENVALUE_FULL_SP")
    RETURN 1
  END SUBROUTINE EIGENVALUE_FULL_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE EIGENVALUE_FULL_DP(A,EVALUES,ERR,ERROR,*)
  
    !#### Subroutine: EIGENVALUE_FULL_DP
    !###  Description:
    !###    Returns the eigenvalues of a full double precision matrix A.
    !###  Parent-subroutine: EIGENVALUE
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:)
    REAL(DP), INTENT(OUT) :: EVALUES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i
    REAL(DP) :: ANGLE,B2,B3,C1,C2,D,Q,Q3,R,RI1,RI2,RI3,RI4,RQ,TEMP,THETA
    
    CALL ENTERS("EIGENVALUE_FULL_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)<= SIZE(EVALUES,1)) THEN
        SELECT CASE(SIZE(A,1))
        CASE(1)
          EVALUES(1)=A(1,1)
        CASE(2)
          IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
            RI1=A(1,1)+A(2,2)
            RI2=A(1,1)*A(2,2)-A(1,2)**2
            B2=RI1/2.0_DP
            C1=RI1*RI1
            C2=4.0_DP*RI2
            IF(C2>C1) CALL FLAG_ERROR("Complex roots found in quadratic equation",ERR,ERROR,*999)
            B3=SQRT(C1-C2)/2.0_DP
            EVALUES(1)=B2+B3
            EVALUES(2)=B2-B3
          ELSE
            EVALUES(1)=A(1,1)
            EVALUES(2)=A(2,2)
          ENDIF
          IF(ABS(EVALUES(2))>ABS(EVALUES(1))) THEN
            TEMP=EVALUES(1)
            EVALUES(1)=EVALUES(2)
            EVALUES(2)=TEMP
          ENDIF
        CASE(3)
          RI1=A(1,1)+A(2,2)+A(3,3)
          RI2=A(1,1)*A(2,2)+A(2,2)*A(3,3)+A(3,3)*A(1,1)-(A(1,2)**2+A(2,3)**2+A(3,1)**2)
          RI3=DETERMINANT(A,ERR,ERROR)
          IF(ERR /=0) GOTO 999
          RI4=RI1/3.0_DP
          Q=RI4*RI4-RI2/3.0_DP   
          R=RI4*(RI4*RI4-RI2/2.0_DP)+RI3/2.0_DP
          Q3=Q*Q*Q
          D=R*R-Q3
          IF(ABS(D)>ZERO_TOLERANCE) CALL FLAG_ERROR("Complex roots found in solution of cubic equation",ERR,ERROR,*999)
          RQ=SQRT(ABS(Q))
          IF(ABS(Q)<ZERO_TOLERANCE) THEN
            THETA=0.0_DP
          ELSE
            THETA=ACOS(R/SQRT(ABS(Q3)))/3.0_DP
          ENDIF
          ANGLE=2.0_DP*REAL(PI,SP)/3.0_DP
          EVALUES(1)=2.0_DP*RQ*COS(THETA)+RI4
          EVALUES(2)=2.0_DP*RQ*COS(THETA+ANGLE)+RI4
          EVALUES(3)=2.0_DP*RQ*COS(THETA+2.0_DP*ANGLE)+RI4
          DO i=1,2
            IF(ABS(EVALUES(3))>ABS(EVALUES(i))) THEN
              TEMP=EVALUES(i)
              EVALUES(i)=EVALUES(3)
              EVALUES(3)=TEMP
            ENDIF
          ENDDO !i
        CASE DEFAULT
          CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Evalues is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EIGENVALUE_FULL_DP")
    RETURN
999 CALL ERRORS("EIGENVALUE_FULL_DP",ERR,ERROR)
    CALL EXITS("EIGENVALUE_FULL_DP")
    RETURN 1
  END SUBROUTINE EIGENVALUE_FULL_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: EIGENVECTOR
  !###  Description:
  !###    Returns the eigenvectors of a matrix.
  !###  Child-subroutines: EIGENVECTOR_FULL_SP,EIGENVECTOR_FULL_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE EIGENVECTOR_FULL_SP(A,EVALUE,EVECTOR,ERR,ERROR,*)
  
    !#### Subroutine: EIGENVECTOR_FULL_SP
    !###  Description:
    !###    Returns the normalised eigenvector of a full single precision symmetric matrix A that corresponds to the eigenvalue
    !###    EVALUE. 
    !###  Parent-subroutine: EIGENVECTOR
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:), EVALUE
    REAL(SP), INTENT(OUT) :: EVECTOR(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,i1,i2,i3,ICYCLE(3,3)
    REAL(SP) :: AL,b(SIZE(A,1)),SUM,U(SIZE(A,1),SIZE(A,2)),x(SIZE(A,1))

    DATA ICYCLE /2,1,1,3,1,2,1,2,1/
    
    CALL ENTERS("EIGENVECTOR_FULL_SP",ERR,ERROR,*999)

!!THIS NEEDS TO BE CHECKED
    
    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)<=SIZE(EVECTOR,1)) THEN
        SELECT CASE(SIZE(A,1))
        CASE(1)
          EVECTOR(1)=1.0_SP
        CASE(2)
          IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
            IF(ABS(A(1,1)-EVALUE)>ABS(A(2,2)-EVALUE)) THEN
              AL=SQRT(A(1,2)**2+(A(1,1)-EVALUE)**2)
              EVECTOR(1)=A(1,2)/AL
              EVECTOR(2)=(EVALUE-A(1,1))/AL
            ELSE
              AL=SQRT(A(1,2)**2+(A(2,2)-EVALUE)**2)
              EVECTOR(1)=(EVALUE-A(2,2))/AL
              EVECTOR(2)=A(1,2)/AL
            ENDIF
          ELSE IF(EVALUE==A(1,1)) THEN
            EVECTOR(1)=1.0_SP
            EVECTOR(2)=0.0_SP
          ELSE IF(EVALUE==A(2,2)) THEN
            EVECTOR(1)=0.0_SP
            EVECTOR(2)=1.0_DP
          ENDIF
        CASE(3)
          IF(ABS(A(1,2))<ZERO_TOLERANCE_SP.AND.ABS(A(1,3))<ZERO_TOLERANCE_SP.AND.ABS(A(2,3))<ZERO_TOLERANCE_SP) THEN
            EVECTOR=0.0_SP
            CALL FLAG_ERROR("Zero matrix?? Eigenvectors undetermined",ERR,ERROR,*999)
          ELSE
            DO i=1,3
              U(i,:)=A(i,:)
              U(i,i)=U(i,i)-EVALUE
            ENDDO !i
            DO i=1,3
              x(i)=1.0_SP
              i1=ICYCLE(i,1)
              i2=ICYCLE(i,2)
              i3=ICYCLE(i,3)
              b(1)=-1.0_SP*U(i1,i)
              b(2)=-1.0_SP*U(i2,i)
              CALL SOLVE_SMALL_LINEAR_SYSTEM(U(i1:i2:i3,i1:i2:i3),x,b,ERR,ERROR,*999)
              SUM=DOT_PRODUCT(U(i,:),X)
              IF(ABS(SUM)<ZERO_TOLERANCE) THEN
                EVECTOR=NORMALISE(X,ERR,ERROR)
                IF(ERR /= 0) GOTO 999
              ENDIF
            ENDDO !i
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Evector is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EIGENVECTOR_FULL_SP")
    RETURN
999 CALL ERRORS("EIGENVECTOR_FULL_SP",ERR,ERROR)
    CALL EXITS("EIGENVECTOR_FULL_SP")
    RETURN 1
  END SUBROUTINE EIGENVECTOR_FULL_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE EIGENVECTOR_FULL_DP(A,EVALUE,EVECTOR,ERR,ERROR,*)
  
    !#### Subroutine: EIGENVECTOR_FULL_DP
    !###  Description:
    !###    Returns the normalised eigenvector of a full double precision symmetric matrix A that corresponds to the eigenvalue
    !###    EVALUE.
    !###  Parent-subroutine: EIGENVECTOR
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:), EVALUE
    REAL(DP), INTENT(OUT) :: EVECTOR(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,i1,i2,i3,ICYCLE(3,3)
    REAL(DP) :: AL,b(SIZE(A,1)),SUM,U(SIZE(A,1),SIZE(A,2)),x(SIZE(A,1))

    DATA ICYCLE /2,1,1,3,1,2,1,2,1/
    
    CALL ENTERS("EIGENVECTOR_FULL_SP",ERR,ERROR,*999)

!!THIS NEEDS TO BE CHECKED

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)<=SIZE(EVECTOR,1)) THEN
        SELECT CASE(SIZE(A,1))
        CASE(1)
          EVECTOR(1)=1.0_DP
        CASE(2)
          IF(ABS(A(1,2))>ZERO_TOLERANCE_SP) THEN
            IF(ABS(A(1,1)-EVALUE)>ABS(A(2,2)-EVALUE)) THEN
              AL=SQRT(A(1,2)**2+(A(1,1)-EVALUE)**2)
              EVECTOR(1)=A(1,2)/AL
              EVECTOR(2)=(EVALUE-A(1,1))/AL
            ELSE
              AL=SQRT(A(1,2)**2+(A(2,2)-EVALUE)**2)
              EVECTOR(1)=(EVALUE-A(2,2))/AL
              EVECTOR(2)=A(1,2)/AL
            ENDIF
          ELSE IF(EVALUE==A(1,1)) THEN
            EVECTOR(1)=1.0_DP
            EVECTOR(2)=0.0_DP
          ELSE IF(EVALUE==A(2,2)) THEN
            EVECTOR(1)=0.0_DP
            EVECTOR(2)=1.0_DP
          ENDIF
        CASE(3)
          IF(ABS(A(1,2))<ZERO_TOLERANCE.AND.ABS(A(1,3))<ZERO_TOLERANCE.AND.ABS(A(2,3))<ZERO_TOLERANCE) THEN
            EVECTOR=0.0_DP
            CALL FLAG_ERROR("Zero matrix?? Eigenvectors undetermined",ERR,ERROR,*999)
          ELSE
            DO i=1,3
              U(i,:)=A(i,:)
              U(i,i)=U(i,i)-EVALUE
            ENDDO !i
            DO i=1,3
              x(i)=1.0_DP
              i1=ICYCLE(i,1)
              i2=ICYCLE(i,2)
              i3=ICYCLE(i,3)
              b(1)=-1.0_DP*U(i1,i)
              b(2)=-1.0_DP*U(i2,i)
              CALL SOLVE_SMALL_LINEAR_SYSTEM(U(i1:i2:i3,i1:i2:i3),x,b,ERR,ERROR,*999)
              SUM=DOT_PRODUCT(U(i,:),X)
              IF(ABS(SUM)<ZERO_TOLERANCE) THEN
                EVECTOR=NORMALISE(X,ERR,ERROR)
                IF(ERR /= 0) GOTO 999
              ENDIF
            ENDDO !i
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Evector is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EIGENVECTOR_FULL_DP")
    RETURN
999 CALL ERRORS("EIGENVECTOR_FULL_DP",ERR,ERROR)
    CALL EXITS("EIGENVECTOR_FULL_DP")
    RETURN 1
  END SUBROUTINE EIGENVECTOR_FULL_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: I0
  !###  Description:
  !###    Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun.
  !###  Child-functions: I0_DP,I0_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION I0_DP(X)
  
    !#### Function: I0_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun,
    !###    for a double precision argument.
    !###  Parent-function: I0
    
    !Argument variables
    REAL(DP), INTENT(IN) :: X
    !Function variable
    REAL(DP) :: I0_DP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: A1=1.0_DP
    REAL(DP), PARAMETER :: A2=3.5156229_DP
    REAL(DP), PARAMETER :: A3=3.0899424_DP
    REAL(DP), PARAMETER :: A4=1.2067492_DP
    REAL(DP), PARAMETER :: A5=0.2659732_DP
    REAL(DP), PARAMETER :: A6=0.0360768_DP
    REAL(DP), PARAMETER :: A7=0.0045813_DP
    REAL(DP) :: T

    !Calculate I0(x) for x < 3.75
    T=X*X/(3.75_DP*3.75_DP)
    I0_DP=A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T
    
    RETURN 
  END FUNCTION I0_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION I0_SP(X)
  
    !#### Function: I0_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun,
    !###    for a single precision argument.
    !###  Parent-function: I0
    
    !Argument variables
    REAL(SP), INTENT(IN) :: X
    !Function variable
    REAL(SP) :: I0_SP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: A1=1.0_SP
    REAL(SP), PARAMETER :: A2=3.5156229_SP
    REAL(SP), PARAMETER :: A3=3.0899424_SP
    REAL(SP), PARAMETER :: A4=1.2067492_SP
    REAL(SP), PARAMETER :: A5=0.2659732_SP
    REAL(SP), PARAMETER :: A6=0.0360768_SP
    REAL(SP), PARAMETER :: A7=0.0045813_SP
    REAL(SP) :: T

    !Calculate I0(x) for x < 3.75
    T=X*X/(3.75_SP*3.75_SP)
    I0_SP=A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T
    
    RETURN 
  END FUNCTION I0_SP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: I1
  !###  Description:
  !###    Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun.
  !###  Child-functions: I1_DP,I1_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION I1_DP(X)
  
    !#### Function: I1_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun,
    !###    for a double precision argument.
    !###  Parent-function: I1
    
    !Argument variables
    REAL(DP), INTENT(IN) :: X
    !Function variable
    REAL(DP) :: I1_DP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: A1=0.5_DP
    REAL(DP), PARAMETER :: A2=0.87890594_DP
    REAL(DP), PARAMETER :: A3=0.51498869_DP
    REAL(DP), PARAMETER :: A4=0.15084934_DP
    REAL(DP), PARAMETER :: A5=0.02658733_DP
    REAL(DP), PARAMETER :: A6=0.00301532_DP
    REAL(DP), PARAMETER :: A7=0.00032411_DP
    REAL(DP) :: T

    !Calculate I1(x)
    T=(X/3.75_DP)*(X/3.75_DP)
    I1_DP=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*X
    
    RETURN 
  END FUNCTION I1_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION I1_SP(X)
  
    !#### Function: I1_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun,
    !###    for a single precision argument.
    !###  Parent-function: I1
    
    !Argument variables
    REAL(SP), INTENT(IN) :: X
    !Function variable
    REAL(SP) :: I1_SP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: A1=0.5_SP
    REAL(SP), PARAMETER :: A2=0.87890594_SP
    REAL(SP), PARAMETER :: A3=0.51498869_SP
    REAL(SP), PARAMETER :: A4=0.15084934_SP
    REAL(SP), PARAMETER :: A5=0.02658733_SP
    REAL(SP), PARAMETER :: A6=0.00301532_SP
    REAL(SP), PARAMETER :: A7=0.00032411_SP
    REAL(SP) :: T

    !Calculate I1(x)
    T=(X/3.75_SP)*(X/3.75_SP)
    I1_SP=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*X
    
    RETURN 
  END FUNCTION I1_SP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-Subroutine: INVERT
  !###  Description:
  !###    Returns the inverse of a matrix.
  !###  Child-Subroutines: INVERT_FULL_SP,INVERT_FULL_DP

  !
  !================================================================================================================================
  !

  SUBROUTINE INVERT_FULL_SP(A,B,DET,ERR,ERROR,*)
    
    !#### Subroutine: INVERT_FULL_SP
    !###  Description:
    !###    INVERT_FULL_SP inverts a full single precision matrix A to give matrix B and returns the determinant of A in DET.
    !###  Parent-Subroutine: INVERT

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:)
    REAL(SP), INTENT(OUT) :: B(:,:)
    REAL(SP), INTENT(OUT) :: DET
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables

    CALL ENTERS("INVERT_FULL_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(B,1)==SIZE(A,1).AND.SIZE(B,2)==SIZE(A,2)) THEN
        SELECT CASE(SIZE(A,1)) 
        CASE(1)
          DET=A(1,1)
          IF(ABS(DET)>ZERO_TOLERANCE_SP) THEN
            B(1,1)=1.0_SP/A(1,1)
          ELSE
            CALL FLAG_WARNING("Matrix A is zero and cannot be inverted",ERR,ERROR,*999)
            B(1,1)=0.0_DP
          ENDIF
        CASE(2)
          DET=A(1,1)*A(2,2)-A(1,2)*A(2,1)
          IF(ABS(DET)>ZERO_TOLERANCE_SP) THEN
            B(1,1)=A(2,2)/DET
            B(1,2)=-A(1,2)/DET
            B(2,1)=-A(2,1)/DET
            B(2,2)=A(1,1)/DET
          ELSE
            CALL FLAG_WARNING("Zero Determinant for matrix A",ERR,ERROR,*999)
            B=0.0_DP
          ENDIF
        CASE(3)
          DET=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)-A(2,1)*A(1,2)*A(3,3)- &
            & A(3,1)*A(2,2)*A(1,3)
          IF(ABS(DET)>ZERO_TOLERANCE_SP) THEN
            B(1,1)=(A(2,2)*A(3,3)-A(3,2)*A(2,3))/DET
            B(2,1)=(A(2,3)*A(3,1)-A(3,3)*A(2,1))/DET
            B(3,1)=(A(2,1)*A(3,2)-A(3,1)*A(2,2))/DET
            B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))/DET
            B(2,2)=(A(3,3)*A(1,1)-A(1,3)*A(3,1))/DET
            B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/DET
            B(1,3)=(A(1,2)*A(2,3)-A(2,2)*A(1,3))/DET
            B(2,3)=(A(1,3)*A(2,1)-A(2,3)*A(1,1))/DET
            B(3,3)=(A(1,1)*A(2,2)-A(2,1)*A(1,2))/DET
          ELSE
            CALL FLAG_WARNING("Zero Determinant for matrix A",ERR,ERROR,*999)
            B=0.0_DP
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Matrix size is not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Matrix B is not the same size as matrix A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix A is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("INVERT_FULL_SP")
    RETURN
999 CALL ERRORS("INVERT_FULL_SP",ERR,ERROR)
    CALL EXITS("INVERT_FULL_SP")
    RETURN 1
  END SUBROUTINE INVERT_FULL_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE INVERT_FULL_DP(A,B,DET,ERR,ERROR,*)
    
    !#### Subroutine: INVERT_FULL_DP
    !###  Description:
    !###    INVERT_FULL_DP inverts a full double precision matrix A to give matrix B and returns the determinant of A in DET.
    !###  Parent-Subroutine: INVERT

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:)
    REAL(DP), INTENT(OUT) :: B(:,:)
    REAL(DP), INTENT(OUT) :: DET
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables

    CALL ENTERS("INVERT_FULL_DP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(B,1)==SIZE(A,1).AND.SIZE(B,2)==SIZE(A,2)) THEN
        SELECT CASE(SIZE(A,1)) 
        CASE(1)
          DET=A(1,1)
          IF(ABS(DET)>ZERO_TOLERANCE) THEN
            B(1,1)=1.0_DP/A(1,1)
          ELSE
            CALL FLAG_WARNING("Matrix A is zero and cannot be inverted",ERR,ERROR,*999)
            B(1,1)=0.0_DP
          ENDIF
        CASE(2)
          DET=A(1,1)*A(2,2)-A(2,1)*A(1,2)
          IF(ABS(DET)>ZERO_TOLERANCE) THEN
            B(1,1)=A(2,2)/DET
            B(1,2)=-A(1,2)/DET
            B(2,1)=-A(2,1)/DET
            B(2,2)=A(1,1)/DET
          ELSE
            CALL FLAG_WARNING("Zero Determinant for matrix A",ERR,ERROR,*999)
            B=0.0_DP
          ENDIF
        CASE(3)
          DET=A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,1)*A(3,2)*A(2,3)-A(2,1)*A(1,2)*A(3,3)- &
            & A(3,1)*A(2,2)*A(1,3)
          IF(ABS(DET)>ZERO_TOLERANCE) THEN
            B(1,1)=(A(2,2)*A(3,3)-A(3,2)*A(2,3))/DET
            B(2,1)=(A(2,3)*A(3,1)-A(3,3)*A(2,1))/DET
            B(3,1)=(A(2,1)*A(3,2)-A(3,1)*A(2,2))/DET
            B(1,2)=(A(3,2)*A(1,3)-A(1,2)*A(3,3))/DET
            B(2,2)=(A(3,3)*A(1,1)-A(1,3)*A(3,1))/DET
            B(3,2)=(A(3,1)*A(1,2)-A(1,1)*A(3,2))/DET
            B(1,3)=(A(1,2)*A(2,3)-A(2,2)*A(1,3))/DET
            B(2,3)=(A(1,3)*A(2,1)-A(2,3)*A(1,1))/DET
            B(3,3)=(A(1,1)*A(2,2)-A(2,1)*A(1,2))/DET
          ELSE
            CALL FLAG_WARNING("Zero Determinant for matrix A",ERR,ERROR,*999)
            B=0.0_DP
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Matrix size is not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Matrix B is not the same size as matrix A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix A is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("INVERT_FULL_DP")
    RETURN
999 CALL ERRORS("INVERT_FULL_DP",ERR,ERROR)
    CALL EXITS("INVERT_FULL_DP")
    RETURN 1
  END SUBROUTINE INVERT_FULL_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: K0
  !###  Description:
  !###    Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun.
  !###  Child-functions: K0_DP,K0_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION K0_DP(X)
  
    !#### Function: K0_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
    !###    for a double precision argument.
    !###  Parent-function: K0
    
    !Argument variables
    REAL(DP), INTENT(IN) :: X
    !Function variable
    REAL(DP) :: K0_DP
    !Local variables
    REAL(DP) :: A1,A2,A3,A4,A5,A6,A7,T

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K0(x)
    IF(X<=2.0_DP) THEN
      T=X*X/4.0_DP
      A1=-0.57721566_DP
      A2=0.42278420_DP
      A3=0.23069756_DP
      A4=0.03488590_DP
      A5=0.00262698_DP
      A6=0.00010750_DP
      A7=0.00000740_DP
      K0_DP=-LOG(X/2.0_DP)*I0(X)+(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)
    ELSE
      IF(X>174.0_DP) THEN
        K0_DP=0.0_DP
      ELSE
        T=2.0_DP/X
        A1=1.25331414_DP
        A2=-0.07832358_DP
        A3=0.02189568_DP
        A4=-0.01062446_DP
        A5=0.00587872_DP
        A6=-0.00251540_DP
        A7=0.00053208_DP
        K0_DP=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*EXP(-X)/SQRT(X)
      ENDIF
    ENDIF
    
    RETURN 
  END FUNCTION K0_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION K0_SP(X)
  
    !#### Function: K0_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Calculates the modified Bessel function of the second kind of order 0 using the approximation of Abromowitz and Stegun,
    !###    for a single precision argument.
    !###  Parent-function: K0
    
    !Argument variables
    REAL(SP), INTENT(IN) :: X
    !Function variable
    REAL(SP) :: K0_SP
    !Local variables
    REAL(SP) :: A1,A2,A3,A4,A5,A6,A7,T

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K0(x)
    IF(X<=2.0_SP) THEN
      T=X*X/4.0_SP
      A1=-0.57721566_SP
      A2=0.42278420_SP
      A3=0.23069756_SP
      A4=0.03488590_SP
      A5=0.00262698_SP
      A6=0.00010750_SP
      A7=0.00000740_SP
      K0_SP=-LOG(X/2.0_SP)*I0(X)+(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)
    ELSE
      IF(X>174.0_SP) THEN
        K0_SP=0.0_SP
      ELSE
        T=2.0_SP/X
        A1=1.25331414_SP
        A2=-0.07832358_SP
        A3=0.02189568_SP
        A4=-0.01062446_SP
        A5=0.00587872_SP
        A6=-0.00251540_SP
        A7=0.00053208_SP
        K0_SP=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*EXP(-X)/SQRT(X)
      ENDIF
    ENDIF
    
    RETURN 
  END FUNCTION K0_SP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: K1
  !###  Description:
  !###    Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun.
  !###  Child-functions: K1_DP,K1_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION K1_DP(X)
  
    !#### Function: K1_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
    !###    for a double precision argument.
    !###  Parent-function: K1
    
    !Argument variables
    REAL(DP), INTENT(IN) :: X
    !Function variable
    REAL(DP) :: K1_DP
    !Local variables
    REAL(DP) :: A1,A2,A3,A4,A5,A6,A7,A8,T

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K1(x)
    IF(X<=2.0_DP) THEN
      T=(X/2.0_DP)*(X/2.0_DP)
      A1=LOG(X/2.0_DP)*I1(X)
      A2=1.0_DP
      A3=0.15443144_DP
      A4=-0.67278579_DP
      A5=-0.18156897_DP
      A6=-0.01919402_DP
      A7=-0.00110404_DP
      A8=-0.00004686_DP
      K1_DP=A1+A2/X+(A3+(A4+(A5+(A6+(A7+A8*T)*T)*T)*T)*T)*X/4
    ELSE
      IF (X>174.0_DP) THEN
        K1_DP=0.0_DP
      ELSE
        T=2.0_DP/X
        A1=1.25331414_DP
        A2=0.23498619_DP
        A3=-0.03655620_DP
        A4=0.01504268_DP
        A5=-0.00780355_DP
        A6=0.00325614_DP
        A7=-0.00068245_DP
        K1_DP=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*EXP(-X)/SQRT(X)
      ENDIF
    ENDIF
    
    RETURN 
  END FUNCTION K1_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION K1_SP(X)
  
    !#### Function: K1_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
    !###    for a single precision argument.
    !###  Parent-function: K1
    
    !Argument variables
    REAL(SP), INTENT(IN) :: X
    !Function variable
    REAL(SP) :: K1_SP
    !Local variables
    REAL(SP) :: A1,A2,A3,A4,A5,A6,A7,A8,T

    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    !Calculate K1(x)
    IF(X<=2.0_SP) THEN
      T=(X/2.0_SP)*(X/2.0_SP)
      A1=LOG(X/2.0_SP)*I1(X)
      A2=1.0_SP
      A3=0.15443144_SP
      A4=-0.67278579_SP
      A5=-0.18156897_SP
      A6=-0.01919402_SP
      A7=-0.00110404_SP
      A8=-0.00004686_SP
      K1_SP=A1+A2/X+(A3+(A4+(A5+(A6+(A7+A8*T)*T)*T)*T)*T)*X/4
    ELSE
      IF (X>174.0_SP) THEN
        K1_SP=0.0_SP
      ELSE
        T=2.0_SP/X
        A1=1.25331414_SP
        A2=0.23498619_SP
        A3=-0.03655620_SP
        A4=0.01504268_SP
        A5=-0.00780355_SP
        A6=0.00325614_SP
        A7=-0.00068245_SP
        K1_SP=(A1+(A2+(A3+(A4+(A5+(A6+A7*T)*T)*T)*T)*T)*T)*EXP(-X)/SQRT(X)
      ENDIF
    ENDIF
    
    RETURN 
  END FUNCTION K1_SP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: KDP
  !###  Description:
  !###    Calculates the elliptic integral of the first kind - K(m).
  !###  Child-functions: KDP_DP,KDP_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION KDP_DP(X)
  
    !#### Function: KDP_DP
    !###  Type: REAL(DP)
    !###  Description:
    !###    Calculates the elliptic integral of the first kind - K(m), for a double precision argument.
    !###  Parent-function: KDP
    
    !Argument variables
    REAL(DP), INTENT(IN) :: X
    !Function variable
    REAL(DP) :: KDP_DP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(DP), PARAMETER :: A0=1.38629436112_DP
    REAL(DP), PARAMETER :: A1=0.09666344259_DP
    REAL(DP), PARAMETER :: A2=0.03590092383_DP
    REAL(DP), PARAMETER :: A3=0.03742563713_DP
    REAL(DP), PARAMETER :: A4=0.01451196212_DP
    REAL(DP), PARAMETER :: B0=0.5_DP
    REAL(DP), PARAMETER :: B1=0.12498593597_DP
    REAL(DP), PARAMETER :: B2=0.06880248576_DP
    REAL(DP), PARAMETER :: B3=0.03328355346_DP
    REAL(DP), PARAMETER :: B4=0.00441787012_DP
    REAL(DP) :: TERM1,TERM2,X1

    X1=1.0_DP-X
    TERM1=A0+(A1+(A2+(A3+A4*X)*X)*X)*X
    TERM2=B0+(B1+(B2+(B3+B4*X)*X)*X)*X
    KDP_DP=TERM1+TERM2*LOG(1.0_DP/X)    
    
    RETURN 
  END FUNCTION KDP_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION KDP_SP(X)
  
    !#### Function: KDP_SP
    !###  Type: REAL(SP)
    !###  Description:
    !###    Calculates the elliptic integral of the first kind - K(m), for a single precision argument.
    !###  Parent-function: KDP
    
    !Argument variables
    REAL(SP), INTENT(IN) :: X
    !Function variable
    REAL(SP) :: KDP_SP
    !Local variables
    !!TODO:: CHECK PRECISION ON THESE CONSTANTS
    REAL(SP), PARAMETER :: A0=1.38629436112_SP
    REAL(SP), PARAMETER :: A1=0.09666344259_SP
    REAL(SP), PARAMETER :: A2=0.03590092383_SP
    REAL(SP), PARAMETER :: A3=0.03742563713_SP
    REAL(SP), PARAMETER :: A4=0.01451196212_SP
    REAL(SP), PARAMETER :: B0=0.5_SP
    REAL(SP), PARAMETER :: B1=0.12498593597_SP
    REAL(SP), PARAMETER :: B2=0.06880248576_SP
    REAL(SP), PARAMETER :: B3=0.03328355346_SP
    REAL(SP), PARAMETER :: B4=0.00441787012_SP
    REAL(SP) :: TERM1,TERM2,X1

    X1=1.0_SP-X
    TERM1=A0+(A1+(A2+(A3+A4*X)*X)*X)*X
    TERM2=B0+(B1+(B2+(B3+B4*X)*X)*X)*X
    KDP_SP=TERM1+TERM2*LOG(1.0_SP/X)    
    
    RETURN 
  END FUNCTION KDP_SP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: L2NORM
  !###  Description:
  !###    Returns the L2 norm of a vector.
  !###  Child-functions: L2NORM_SP,L2NORM_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION L2NORM_SP(A)

    !#### Function: L2NORM_SP
    !###  Description:
    !###    L2NORM_SP returns the L2-norm of the single precision vector A.
    !###  Parent-function: L2NORM

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:)
    !Function variable
    REAL(SP) :: L2NORM_SP
    !Local variables
    REAL(SP) :: ASUM
    
    ASUM=SUM(A*A,1)
    L2NORM_SP=SQRT(ASUM)

    RETURN
  END FUNCTION L2NORM_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION L2NORM_DP(A)

    !#### Function: L2NORM_DP
    !###  Description:
    !###    L2NORM_DP returns the L2-norm of the double precision vector A.
    !###  Parent-function: L2NORM

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:)
    !Function variable
    REAL(DP) :: L2NORM_DP
    !Local variables
    REAL(DP) :: ASUM

    ASUM=SUM(A*A,1)
    L2NORM_DP=SQRT(ASUM)

    RETURN
  END FUNCTION L2NORM_DP

  !
  !================================================================================================================================
  !

  SUBROUTINE MATRIX_PRODUCT_SP(A,B,C,ERR,ERROR,*)
  
    !#### Subroutine: MATRIX_PRODUCT_SP
    !###  Description:
    !###    Calculates and returns the matrix-prouct of the single precision matrix A*B in C.
    !###  Parent-function: MATRIX_PRODUCT
    !Author: Kumar Mithraratne
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(3,3),B(3,3)
    REAL(SP), INTENT(OUT) :: C(3,3)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,j,k 
    
    CALL ENTERS("MATRIX_PRODUCT_SP",ERR,ERROR,*999)

    DO i=1,3,1
      DO j=1,3,1
        C(i,j)=0.0_SP
	DO k=1,3,1	      
	  C(i,j)=C(i,j)+A(i,k)*B(k,j) 
	ENDDO  	    
      ENDDO
    ENDDO

    CALL EXITS("MATRIX_PRODUCT_SP")
    RETURN
999 CALL ERRORS("MATRIX_PRODUCT_SP",ERR,ERROR)
    CALL EXITS("MATRIX_PRODUCT_SP")
    RETURN 1
  END SUBROUTINE MATRIX_PRODUCT_SP

  !
  !================================================================================================================================
  !

  SUBROUTINE MATRIX_PRODUCT_DP(A,B,C,ERR,ERROR,*)
  
    !#### Subroutine: MATRIX_PRODUCT_DP
    !###  Description:
    !###    Calculates and returns the matrix-prouct of the double precision matrix A*B in C.
    !###  Parent-function: MATRIX_PRODUCT
    !Author: Kumar Mithraratne
        
    !Argument variables
    REAL(DP), INTENT(IN) :: A(3,3),B(3,3)
    REAL(DP), INTENT(OUT) :: C(3,3)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,j,k 
        
    CALL ENTERS("MATRIX_PRODUCT_DP",ERR,ERROR,*999)

    DO i=1,3,1
      DO j=1,3,1
        C(i,j)=0.0_DP
	DO k=1,3,1	      
	  C(i,j)=C(i,j)+A(i,k)*B(k,j) 
	ENDDO  	    
      ENDDO
    ENDDO

    CALL EXITS("MATRIX_PRODUCT_DP")
    RETURN
999 CALL ERRORS("MATRIX_PRODUCT_DP",ERR,ERROR)
    CALL EXITS("MATRIX_PRODUCT_DP")
    RETURN 1
  END SUBROUTINE MATRIX_PRODUCT_DP

  !
  !================================================================================================================================
  !

  SUBROUTINE MATRIX_TRANSPOSE_SP(A,AT,ERR,ERROR,*)
  
    !#### Subroutine: MATRIX_TRANSPOSE_SP
    !###  Description:
    !###    Returns the transpose of a single precision matrix Ain AT.
    !###  Parent-function: MATRIX_TRANSPOSE
    !Author: Kumar Mithraratne
        
    !Argument variables
    REAL(SP), INTENT(IN) :: A(3,3)
    REAL(SP), INTENT(OUT) :: AT(3,3)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,j
    
    CALL ENTERS("MATRIX_TRANSPOSE_SP",ERR,ERROR,*999)

    DO i=1,3,1
      DO j=1,3,1
	AT(i,j)=A(j,i)
      ENDDO
    ENDDO

    CALL EXITS("MATRIX_TRANSPOSE_SP")
    RETURN
999 CALL ERRORS("MATRIX_TRANSPOSE_SP",ERR,ERROR)
    CALL EXITS("MATRIX_TRANSPOSE_SP")
    RETURN 1
  END SUBROUTINE MATRIX_TRANSPOSE_SP

  !
  !================================================================================================================================
  !

  SUBROUTINE MATRIX_TRANSPOSE_DP(A,AT,ERR,ERROR,*)
  
    !#### Subroutine: MATRIX_TRANSPOSE_DP
    !###  Description:
    !###    Returns the transpose of a double precision matrix Ain AT.
    !###  Parent-function: MATRIX_TRANSPOSE
    !Author: Kumar Mithraratne
        
    !Argument variables
    REAL(DP), INTENT(IN) :: A(3,3)
    REAL(DP), INTENT(OUT) :: AT(3,3)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: i,j
        
    CALL ENTERS("MATRIX_TRANSPOSE_DP",ERR,ERROR,*999)

    DO i=1,3,1
      DO j=1,3,1
	AT(i,j)=A(j,i)
      ENDDO
    ENDDO

    CALL EXITS("MATRIX_TRANSPOSE_DP")
    RETURN
999 CALL ERRORS("MATRIX_TRANSPOSE_DP",ERR,ERROR)
    CALL EXITS("MATRIX_TRANSPOSE_DP")
    RETURN 1
  END SUBROUTINE MATRIX_TRANSPOSE_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-Function: NORMALISE
  !###  Description:
  !###    Normalises a vector
  !###  Child-functions: NORMALISE_SP, NORMALISE_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION NORMALISE_SP(A,ERR,ERROR)

    !#### Function: NORMALISE_SP
    !###  Description:
    !###    NORMALISE_SP normalises a real single precision vector A.
    !###  Parent-function: NORMALISE

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(SP) :: NORMALISE_SP(SIZE(A,1))
    !Local variables
    REAL(SP) :: LENGTH
    
    CALL ENTERS("NORMALISE_SP",ERR,ERROR,*999)

    LENGTH=L2NORM(A)
    IF(ABS(LENGTH)<ZERO_TOLERANCE_SP) THEN
      NORMALISE_SP=A/LENGTH
    ELSE
      NORMALISE_SP=A
      CALL FLAG_ERROR("Length of vector to normalise is zero",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("NORMALISE_SP")
    RETURN
999 CALL ERRORS("NORMALISE_SP",ERR,ERROR)
    CALL EXITS("NORMALISE_SP")
    RETURN    
  END FUNCTION NORMALISE_SP

  !
  !================================================================================================================================
  !
  
  FUNCTION NORMALISE_DP(A,ERR,ERROR)

    !#### Function: NORMALISE_DP
    !###  Description:
    !###    NORMALISE_DP normalises a real double precision vector A.
    !###  Parent-function: NORMALISE

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: NORMALISE_DP(SIZE(A,1))
    !Local variables
    REAL(DP) :: LENGTH
    
    CALL ENTERS("NORMALISE_DP",ERR,ERROR,*999)

    LENGTH=L2NORM(A)
    IF(ABS(LENGTH)<ZERO_TOLERANCE) THEN
      NORMALISE_DP=A/LENGTH
    ELSE
      NORMALISE_DP=A
      CALL FLAG_ERROR("Length of vector to normalise is zero",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("NORMALISE_DP")
    RETURN
999 CALL ERRORS("NORMALISE_DP",ERR,ERROR)
    CALL EXITS("NORMALISE_DP")
    RETURN    
  END FUNCTION NORMALISE_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: SOLVE_SMALL_LINEAR_SYSTEM
  !###  Description:
  !###    Solves a small linear system Ax=b.
  !###  Child-subroutines: SOLVE_SMALL_LINEAR_SYSTEM_SP,SOLVE_SMALL_LINEAR_SYSTEM_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_SP(A,x,b,ERR,ERROR,*)
  
    !#### Subroutine: SOLVE_SMALL_LINEAR_SYSTEM_SP
    !###  Description:
    !###    Finds the solution to a small single precision linear system Ax=b.
    !###  Parent-subroutine: SOLVE_SMALL_LINEAR_SYSTEM
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:)
    REAL(SP), INTENT(OUT) :: x(:)
    REAL(SP), INTENT(IN) :: b(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    REAL(SP) :: AINV(SIZE(A,1),SIZE(A,2)),ADET
    
    CALL ENTERS("SOLVE_SMALL_LINEAR_SYSTEM_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)==SIZE(b,1)) THEN
        IF(SIZE(A,1)<=SIZE(x,1)) THEN
          SELECT CASE(SIZE(A,1))
          CASE(1:3)
            CALL INVERT(A,AINV,ADET,ERR,ERROR,*999)
            x=MATMUL(AINV,b)
          CASE DEFAULT
            CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("x is too small",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Size of b is not the same as the number of rows in A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("SOLVE_SMALL_LINEAR_SYSTEM_SP")
    RETURN
999 CALL ERRORS("SOLVE_SMALL_LINEAR_SYSTEM_SP",ERR,ERROR)
    CALL EXITS("SOLVE_SMALL_LINEAR_SYSTEM_SP")
    RETURN 1
  END SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_SP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_DP(A,x,b,ERR,ERROR,*)
  
    !#### Subroutine: SOLVE_SMALL_LINEAR_SYSTEM_DP
    !###  Description:
    !###    Finds the solution to a small double precision linear system Ax=b.
    !###  Parent-subroutine: SOLVE_SMALL_LINEAR_SYSTEM
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:)
    REAL(DP), INTENT(OUT) :: x(:)
    REAL(DP), INTENT(IN) :: b(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    REAL(DP) :: AINV(SIZE(A,1),SIZE(A,2)),ADET
    
    CALL ENTERS("SOLVE_SMALL_LINEAR_SYSTEM_DP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)==SIZE(b,1)) THEN
        IF(SIZE(A,1)<=SIZE(x,1)) THEN
          SELECT CASE(SIZE(A,1))
          CASE(1:3)
            CALL INVERT(A,AINV,ADET,ERR,ERROR,*999)
            x=MATMUL(AINV,b)
          CASE DEFAULT
            CALL FLAG_ERROR("Matrix size not implemented",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("x is too small",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Size of b is not the same as the number of rows in A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("SOLVE_SMALL_LINEAR_SYSTEM_DP")
    RETURN
999 CALL ERRORS("SOLVE_SMALL_LINEAR_SYSTEM_DP",ERR,ERROR)
    CALL EXITS("SOLVE_SMALL_LINEAR_SYSTEM_DP")
    RETURN 1
  END SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_DP

  !
  !================================================================================================================================
  !
  
END MODULE MATHS
