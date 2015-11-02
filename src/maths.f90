!> \file
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
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
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
MODULE Maths

  USE BASE_ROUTINES
  USE Constants
  USE Kinds
  USE ISO_VARYING_STRING
  USE Strings

#include "macros.h"  
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  !>Calculates the vector cross product of two vectors
  INTERFACE CROSS_PRODUCT
    MODULE PROCEDURE CROSS_PRODUCT_INTG
    MODULE PROCEDURE CROSS_PRODUCT_SP
    MODULE PROCEDURE CROSS_PRODUCT_DP
  END INTERFACE CROSS_PRODUCT

  !>Calculates the vector cross product of two vectors
  INTERFACE CrossProduct
    MODULE PROCEDURE CROSS_PRODUCT_INTG
    MODULE PROCEDURE CROSS_PRODUCT_SP
    MODULE PROCEDURE CROSS_PRODUCT_DP
  END INTERFACE CrossProduct

  !>Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the derivatives D_A and D_B of A and B
  INTERFACE D_CROSS_PRODUCT
    MODULE PROCEDURE D_CROSS_PRODUCT_INTG
    MODULE PROCEDURE D_CROSS_PRODUCT_SP
    MODULE PROCEDURE D_CROSS_PRODUCT_DP
  END INTERFACE D_CROSS_PRODUCT

  !>Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the derivatives D_A and D_B of A and B
  INTERFACE dCrossProduct
    MODULE PROCEDURE D_CROSS_PRODUCT_INTG
    MODULE PROCEDURE D_CROSS_PRODUCT_SP
    MODULE PROCEDURE D_CROSS_PRODUCT_DP
  END INTERFACE dCrossProduct

  !>Calculates and returns the MATRIX-VECTOR-prouct of the double precision VECTOR A*B in C.
  INTERFACE MATRIX_VECTOR_PRODUCT
    MODULE PROCEDURE MATRIX_VECTOR_PRODUCT_SP
    MODULE PROCEDURE MATRIX_VECTOR_PRODUCT_DP
  END INTERFACE !MATRIX_VECTOR_PRODUCT

  !>Returns the determinant of a matrix
  INTERFACE Determinant
    MODULE PROCEDURE DETERMINANT_FULL_INTG
    MODULE PROCEDURE DETERMINANT_FULL_SP
    MODULE PROCEDURE DETERMINANT_FULL_DP
  END INTERFACE Determinant
        
  !>Calculates the elliptic integral of the second kind - E(m).
  INTERFACE Edp
    MODULE PROCEDURE EDP_DP
    MODULE PROCEDURE EDP_SP
  END INTERFACE Edp

  !>Returns the eigenvalues of a matrix.
  INTERFACE Eigenvalue
    MODULE PROCEDURE EIGENVALUE_FULL_SP
    MODULE PROCEDURE EIGENVALUE_FULL_DP
  END INTERFACE Eigenvalue

  !>Returns the eigenvectors of a matrix.
  INTERFACE Eigenvector
    MODULE PROCEDURE EIGENVECTOR_FULL_SP
    MODULE PROCEDURE EIGENVECTOR_FULL_DP
  END INTERFACE Eigenvector

  !>Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun.
  INTERFACE I0
    MODULE PROCEDURE I0_DP
    MODULE PROCEDURE I0_SP
  END INTERFACE !I0

  !>Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun.
  INTERFACE I1
    MODULE PROCEDURE I1_DP
    MODULE PROCEDURE I1_SP
  END INTERFACE !I1

  !>Returns the inverse of a matrix.
  INTERFACE Invert
    MODULE PROCEDURE INVERT_FULL_SP
    MODULE PROCEDURE INVERT_FULL_DP
  END INTERFACE Invert

  !>Calculates the modified Bessel FUNCTION of the second kind of order 0 using the approximation of Abromowitz and Stegun.
  INTERFACE K0
    MODULE PROCEDURE K0_DP
    MODULE PROCEDURE K0_SP
  END INTERFACE !K0

  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun.
  INTERFACE K1
    MODULE PROCEDURE K1_DP
    MODULE PROCEDURE K1_SP
  END INTERFACE !K1

  !>Calculates the elliptic integral of the first kind - K(m).
  INTERFACE KDP
    MODULE PROCEDURE KDP_DP
    MODULE PROCEDURE K1_SP
  END INTERFACE !KDP

  !>Returns the identity matrix.
  INTERFACE IdentityMatrix
    MODULE PROCEDURE IdentityMatrixSP
    MODULE PROCEDURE IdentityMatrixDP
  END INTERFACE IdentityMatrix

  !>Returns the L2 norm of a vector.
  INTERFACE L2Norm
    MODULE PROCEDURE L2NORM_SP
    MODULE PROCEDURE L2NORM_DP
  END INTERFACE L2Norm

  !>Calculates and returns the matrix-prouct A*B in the matrix C.
  INTERFACE MATRIX_PRODUCT
    MODULE PROCEDURE MATRIX_PRODUCT_SP
    MODULE PROCEDURE MATRIX_PRODUCT_DP
  END INTERFACE MATRIX_PRODUCT
  
  !>Calculates and returns the matrix-prouct A*B in the matrix C.
  INTERFACE MatrixProduct
    MODULE PROCEDURE MATRIX_PRODUCT_SP
    MODULE PROCEDURE MATRIX_PRODUCT_DP
  END INTERFACE MatrixProduct
  
  !>Returns the transpose of a matrix A in A^T.
  INTERFACE MATRIX_TRANSPOSE
    MODULE PROCEDURE MATRIX_TRANSPOSE_SP
    MODULE PROCEDURE MATRIX_TRANSPOSE_DP
  END INTERFACE MATRIX_TRANSPOSE

  !>Returns the transpose of a matrix A in A^T.
  INTERFACE MatrixTranspose
    MODULE PROCEDURE MATRIX_TRANSPOSE_SP
    MODULE PROCEDURE MATRIX_TRANSPOSE_DP
  END INTERFACE MatrixTranspose

  !>Calculates and returns the matrix-vector-prouct A*b in the vector c.
  INTERFACE MATRIX_VECTOR_PRODUCT
    MODULE PROCEDURE MATRIX_VECTOR_PRODUCT_SP
    MODULE PROCEDURE MATRIX_VECTOR_PRODUCT_DP
  END INTERFACE MATRIX_VECTOR_PRODUCT

  !>Calculates and returns the matrix-vector-prouct A*b in the vector c.
  INTERFACE MatrixVectorProduct
    MODULE PROCEDURE MATRIX_VECTOR_PRODUCT_SP
    MODULE PROCEDURE MATRIX_VECTOR_PRODUCT_DP
  END INTERFACE MatrixVectorProduct
  
  !>Normalises a vector
  INTERFACE Normalise
    MODULE PROCEDURE NORMALISE_SP
    MODULE PROCEDURE NORMALISE_DP
  END INTERFACE Normalise

  !>Calculates the normalised vector cross product of two vectors
  INTERFACE NORM_CROSS_PRODUCT
    MODULE PROCEDURE NORM_CROSS_PRODUCT_SP
    MODULE PROCEDURE NORM_CROSS_PRODUCT_DP
  END INTERFACE NORM_CROSS_PRODUCT

  !>Calculates the normalised vector cross product of two vectors
  INTERFACE NormCrossProduct
    MODULE PROCEDURE NORM_CROSS_PRODUCT_SP
    MODULE PROCEDURE NORM_CROSS_PRODUCT_DP
  END INTERFACE NormCrossProduct

  !>Solves a small linear system Ax=b.
  INTERFACE SOLVE_SMALL_LINEAR_SYSTEM
    MODULE PROCEDURE SOLVE_SMALL_LINEAR_SYSTEM_SP
    MODULE PROCEDURE SOLVE_SMALL_LINEAR_SYSTEM_DP
  END INTERFACE SOLVE_SMALL_LINEAR_SYSTEM

  !>Solves a small linear system Ax=b.
  INTERFACE SolveSmallLinearSystem
    MODULE PROCEDURE SOLVE_SMALL_LINEAR_SYSTEM_SP
    MODULE PROCEDURE SOLVE_SMALL_LINEAR_SYSTEM_DP
  END INTERFACE SolveSmallLinearSystem

  !>Returns hyperbolic cotangent of argument
  INTERFACE COTH
    MODULE PROCEDURE COTH_SP
    MODULE PROCEDURE COTH_DP
  END INTERFACE COTH 

  PUBLIC CROSS_PRODUCT,CrossProduct,D_CROSS_PRODUCT,dCrossProduct,Determinant,Eigenvalue,Eigenvector,IdentityMatrix,Invert, &
    & L2Norm,MATRIX_PRODUCT,MatrixProduct,MATRIX_TRANSPOSE,MatrixTranspose,Normalise,NORM_CROSS_PRODUCT,NormCrossProduct, &
    & SOLVE_SMALL_LINEAR_SYSTEM,SolveSmallLinearSystem,Coth,spline_cubic_set,s3_fs,spline_cubic_val,MATRIX_VECTOR_PRODUCT, &
    & MatrixVectorProduct
  
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector cross-prouct of the integer vectors A*B in C.
  SUBROUTINE CROSS_PRODUCT_INTG(A,B,C,ERR,ERROR,*)
      
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:) !<The first vector in the cross product
    INTEGER(INTG), INTENT(IN) :: B(:) !<The second vector in the cross product
    INTEGER(INTG), INTENT(OUT) :: C(:) !<On exit, the cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("CROSS_PRODUCT_INTG",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(B,1)) THEN
      IF(SIZE(C,1)==3) THEN
        SELECT CASE(SIZE(A,1))
        CASE(3)
          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        CASE DEFAULT
          CALL FlagError("Invalid vector size",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("The vector C is not the correct size",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("The vectors A and B are not the same size",ERR,ERROR,*999)
    ENDIF

    EXITS("CROSS_PRODUCT_INTG")
    RETURN
999 ERRORSEXITS("CROSS_PRODUCT_INTG",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CROSS_PRODUCT_INTG
  
  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector cross-prouct of the single precision vectors A*B in C.
  SUBROUTINE CROSS_PRODUCT_SP(A,B,C,ERR,ERROR,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:) !<The first vector in the cross product
    REAL(SP), INTENT(IN) :: B(:) !<The second vector in the cross product
    REAL(SP), INTENT(OUT) :: C(:) !<On exit, the cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("CROSS_PRODUCT_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(B,1)) THEN
      IF(SIZE(C,1)==3) THEN
        SELECT CASE(SIZE(A,1))
        CASE(3)
          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        CASE DEFAULT
          CALL FlagError("Invalid vector size",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("The vector C is not the correct size",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("The vectors A and B are not the same size",ERR,ERROR,*999)
    ENDIF

    EXITS("CROSS_PRODUCT_SP")
    RETURN
999 ERRORSEXITS("CROSS_PRODUCT_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CROSS_PRODUCT_SP
  
  !
  !================================================================================================================================
  !

  !>Calculates and returns the vector cross-prouct of the double precision vectors A*B in C.
  SUBROUTINE CROSS_PRODUCT_DP(A,B,C,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:) !<The first vector in the cross product
    REAL(DP), INTENT(IN) :: B(:) !<The second vector in the cross product
    REAL(DP), INTENT(OUT) :: C(:) !<On exit, the cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("CROSS_PRODUCT_DP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(B,1)) THEN
      IF(SIZE(C,1)==3) THEN
        SELECT CASE(SIZE(A,1))
        CASE(3)
          C(1)=A(2)*B(3)-A(3)*B(2)
          C(2)=A(3)*B(1)-A(1)*B(3)
          C(3)=A(1)*B(2)-A(2)*B(1)
        CASE DEFAULT
          CALL FlagError("Invalid vector size",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("The vector C is not the correct size",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("The vectors A and B are not the same size",ERR,ERROR,*999)
    ENDIF

    EXITS("CROSS_PRODUCT_DP")
    RETURN
999 ERRORSEXITS("CROSS_PRODUCT_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CROSS_PRODUCT_DP
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
  !>derivatives D_A and D_B of A and B for integer vectors.
  SUBROUTINE D_CROSS_PRODUCT_INTG(N,A,B,C,D_A,D_B,D_C,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N !<The number of derivatives
    INTEGER(INTG), INTENT(IN) :: A(:) !<The A vector
    INTEGER(INTG), INTENT(IN) :: B(:) !<The B vector
    INTEGER(INTG), INTENT(OUT) :: C(:) !<On exit, the cross product of A*B
    INTEGER(INTG), INTENT(IN) :: D_A(:,:) !<The N derivatives of A
    INTEGER(INTG), INTENT(IN) :: D_B(:,:) !<The N derivatives of B
    INTEGER(INTG), INTENT(OUT) :: D_C(:,:) !<On exit, the derivatives of C
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: ni
    
    ENTERS("D_CROSS_PRODUCT_INTG",ERR,ERROR,*999)

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
            CALL FlagError("Invalid vector size",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("The vector C is not the correct size",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("The number of derivative vectors is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("The vectors for D_A and D_B are not the same size",ERR,ERROR,*999)
    ENDIF

    EXITS("D_CROSS_PRODUCT_INTG")
    RETURN
999 ERRORSEXITS("D_CROSS_PRODUCT_INTG",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE D_CROSS_PRODUCT_INTG
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
  !>derivatives D_A and D_B of A and B for single precision vectors.
  SUBROUTINE D_CROSS_PRODUCT_SP(N,A,B,C,D_A,D_B,D_C,ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N !<The number of derivatives
    REAL(SP), INTENT(IN) :: A(:) !<The A vector
    REAL(SP), INTENT(IN) :: B(:) !<The B vector
    REAL(SP), INTENT(OUT) :: C(:) !<On exit, the cross product of A*B
    REAL(SP), INTENT(IN) :: D_A(:,:) !<The N derivatives of A
    REAL(SP), INTENT(IN) :: D_B(:,:) !<The N derivatives of B
    REAL(SP), INTENT(OUT) :: D_C(:,:) !<On exit, the derivatives of C
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: ni
    
    ENTERS("D_CROSS_PRODUCT_SP",ERR,ERROR,*999)

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
            CALL FlagError("Invalid vector size",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("The vector C is not the correct size",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("The number of derivative vectors is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("The vectors for D_A and D_B are not the same size",ERR,ERROR,*999)
    ENDIF

    EXITS("D_CROSS_PRODUCT_SP")
    RETURN
999 ERRORSEXITS("D_CROSS_PRODUCT_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE D_CROSS_PRODUCT_SP
  
  !
  !================================================================================================================================
  !

  !>Calculates the the vector cross product of A*B in C and the N derivatives, D_C, of the vector cross product given the 
  !>derivatives D_A and D_B of A and B for double precision vectors.
  SUBROUTINE D_CROSS_PRODUCT_DP(N,A,B,C,D_A,D_B,D_C,ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: N !<The number of derivatives
    REAL(DP), INTENT(IN) :: A(:) !<The A vector
    REAL(DP), INTENT(IN) :: B(:) !<The B vector
    REAL(DP), INTENT(OUT) :: C(:) !<On exit, the cross product of A*B
    REAL(DP), INTENT(IN) :: D_A(:,:) !<The N derivatives of A
    REAL(DP), INTENT(IN) :: D_B(:,:) !<The N derivatives of B
    REAL(DP), INTENT(OUT) :: D_C(:,:) !<On exit, the derivatives of C
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: ni
    
    ENTERS("D_CROSS_PRODUCT_DP",ERR,ERROR,*999)

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
            CALL FlagError("Invalid vector size",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("The vector C is not the correct size",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("The number of derivative vectors is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("The vectors for D_A and D_B are not the same size",ERR,ERROR,*999)
    ENDIF

    EXITS("D_CROSS_PRODUCT_DP")
    RETURN
999 ERRORSEXITS("D_CROSS_PRODUCT_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE D_CROSS_PRODUCT_DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the MATRIX-VECTOR-prouct of the single precision VECTOR A*B in C.
  SUBROUTINE MATRIX_VECTOR_PRODUCT_SP(A,B,C,err,error,*)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:)  !<The A MATRIX
    REAL(SP), INTENT(IN) :: B(:)    !<The B VECTOR
    REAL(SP), INTENT(OUT) :: C(:)   !<On exit, the product VECTOR C=A*B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("MATRIX_VECTOR_PRODUCT_SP",err,error,*999)

    IF(SIZE(A,2)==SIZE(B,1).AND.SIZE(A,1)==SIZE(C,1)) THEN
       SELECT CASE(SIZE(A,1))
       CASE(1)
         C(1)=A(1,1)*B(1)
       CASE(2)
         C(1)=A(1,1)*B(1)+A(1,2)*B(2)
         C(2)=A(2,1)*B(1)+A(2,2)*B(2)
       CASE(3)
         C(1)=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
         C(2)=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
         C(3)=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
       CASE DEFAULT
         CALL FlagError("Invalid matrix and/or vector size.",err,error,*999)
       END SELECT
     ELSE
       CALL FlagError("Invalid matrix sizes.",err,error,*999)
     ENDIF

     EXITS("MATRIX_VECTOR_PRODUCT_SP")
     RETURN
 999 ERRORSEXITS("MATRIX_VECTOR_PRODUCT_SP",err,error)
     RETURN 1
     
  END SUBROUTINE MATRIX_VECTOR_PRODUCT_SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the MATRIX-VECTOR-prouct of the double precision VECTOR A*B in C.
  SUBROUTINE MATRIX_VECTOR_PRODUCT_DP(A,B,C,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:)  !<The A MATRIX
    REAL(DP), INTENT(IN) :: B(:)    !<The B VECTOR
    REAL(DP), INTENT(OUT) :: C(:)   !<On exit, the product VECTOR C=A*B
    INTEGER(INTG) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("MATRIX_VECTOR_PRODUCT_DP",err,error,*999)

    IF(SIZE(A,2)==SIZE(B,1).AND.SIZE(A,1)==SIZE(C,1)) THEN
       SELECT CASE(SIZE(A,1))
       CASE(1)
         C(1)=A(1,1)*B(1)
       CASE(2)
         C(1)=A(1,1)*B(1)+A(1,2)*B(2)
         C(2)=A(2,1)*B(1)+A(2,2)*B(2)
       CASE(3)
         C(1)=A(1,1)*B(1)+A(1,2)*B(2)+A(1,3)*B(3)
         C(2)=A(2,1)*B(1)+A(2,2)*B(2)+A(2,3)*B(3)
         C(3)=A(3,1)*B(1)+A(3,2)*B(2)+A(3,3)*B(3)
       CASE DEFAULT
         CALL FlagError("Invalid matrix and/or vector size.",err,error,*999)
       END SELECT
     ELSE
       CALL FlagError("Invalid matrix sizes.",err,error,*999)
     ENDIF

     EXITS("MATRIX_VECTOR_PRODUCT_DP")
     RETURN
 999 ERRORSEXITS("MATRIX_VECTOR_PRODUCT_DP",err,error)
     RETURN 1
     
  END SUBROUTINE MATRIX_VECTOR_PRODUCT_DP
  
  !
  !================================================================================================================================
  !

  !>Returns the determinant of a full integer matrix A.
  FUNCTION DETERMINANT_FULL_INTG(A,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:,:) !<The matrix to find the determinant of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: DETERMINANT_FULL_INTG
    
    ENTERS("DETERMINANT_FULL_INTG",ERR,ERROR,*999)

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
        CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("DETERMINANT_FULL_INTG")
    RETURN
999 ERRORSEXITS("DETERMINANT_FULL_INTG",ERR,ERROR)
    RETURN
    
  END FUNCTION DETERMINANT_FULL_INTG
  
  !
  !================================================================================================================================
  !

  !>Returns the determinant of a full single precision matrix A.
  FUNCTION DETERMINANT_FULL_SP(A,ERR,ERROR)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the determinant of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(SP) :: DETERMINANT_FULL_SP
    
    ENTERS("DETERMINANT_FULL_SP",ERR,ERROR,*999)

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
        CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("DETERMINANT_FULL_SP")
    RETURN
999 ERRORSEXITS("DETERMINANT_FULL_SP",ERR,ERROR)
    RETURN
    
  END FUNCTION DETERMINANT_FULL_SP
  
  !
  !================================================================================================================================
  !

  !>Returns the determinant of a full double precision matrix A
  FUNCTION DETERMINANT_FULL_DP(A,ERR,ERROR)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the determinant of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: DETERMINANT_FULL_DP
    
    ENTERS("DETERMINANT_FULL_DP",ERR,ERROR,*999)

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
        CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("DETERMINANT_FULL_DP")
    RETURN
999 ERRORSEXITS("DETERMINANT_FULL_DP",ERR,ERROR)
    RETURN
    
  END FUNCTION DETERMINANT_FULL_DP

  !
  !================================================================================================================================
  !

  !>Calculates the elliptic integral of the second kind - E(m), for a double precision argument.
  PURE FUNCTION EDP_DP(X)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: X !<The value to evaluate the function at
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

  !>Calculates the elliptic integral of the second kind - E(m), for a single precision argument.
  PURE FUNCTION EDP_SP(X)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: X !<The value to evaluate the function at
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

  !>Returns the eigenvalues of a full single precision matrix A.
  SUBROUTINE EIGENVALUE_FULL_SP(A,EVALUES,ERR,ERROR,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the eignenvalues for
    REAL(SP), INTENT(OUT) :: EVALUES(:) !<On exit, the eignevalues
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(SP) :: ANGLE,B2,B3,C1,C2,D,Q,Q3,R,RI1,RI2,RI3,RI4,RQ,TEMP,THETA
    
    ENTERS("EIGENVALUE_FULL_SP",ERR,ERROR,*999)

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
            IF(C2>C1) CALL FlagError("Complex roots found in quadratic equation",ERR,ERROR,*999)
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
          IF(ABS(D)>ZERO_TOLERANCE_SP) CALL FlagError("Complex roots found in solution of cubic equation",ERR,ERROR,*999)
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
          CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Evalues is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("EIGENVALUE_FULL_SP")
    RETURN
999 ERRORSEXITS("EIGENVALUE_FULL_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EIGENVALUE_FULL_SP

  !
  !================================================================================================================================
  !

  !>Returns the eigenvalues of a full double precision matrix A.
  SUBROUTINE EIGENVALUE_FULL_DP(A,EVALUES,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the eigenvalues of
    REAL(DP), INTENT(OUT) :: EVALUES(:) !<On exit, the eigenvalues of the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(DP) :: ANGLE,B2,B3,C1,C2,D,Q,Q3,R,RI1,RI2,RI3,RI4,RQ,TEMP,THETA
    
    ENTERS("EIGENVALUE_FULL_SP",ERR,ERROR,*999)

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
            IF(C2>C1) CALL FlagError("Complex roots found in quadratic equation",ERR,ERROR,*999)
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
          IF(ABS(D)>ZERO_TOLERANCE) CALL FlagError("Complex roots found in solution of cubic equation",ERR,ERROR,*999)
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
          CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Evalues is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("EIGENVALUE_FULL_DP")
    RETURN
999 ERRORSEXITS("EIGENVALUE_FULL_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EIGENVALUE_FULL_DP

  !
  !================================================================================================================================
  !

  !>Returns the normalised eigenvector of a full single precision symmetric matrix A that corresponds to the eigenvalue EVALUE. 
  SUBROUTINE EIGENVECTOR_FULL_SP(A,EVALUE,EVECTOR,ERR,ERROR,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to find the eignevectors for
    REAL(SP), INTENT(IN) :: EVALUE !<The eigenvalue to find the eignevector for
    REAL(SP), INTENT(OUT) :: EVECTOR(:) !<On exit, the eigenvector corresponding the the eigenvalue
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,i1,i2,i3,ICYCLE(3,3)
    REAL(SP) :: AL,b(SIZE(A,1)),SUM,U(SIZE(A,1),SIZE(A,2)),x(SIZE(A,1))

    DATA ICYCLE /2,1,1,3,1,2,1,2,1/
    
    ENTERS("EIGENVECTOR_FULL_SP",ERR,ERROR,*999)

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
          ELSE IF(ABS(EVALUE-A(1,1))<ZERO_TOLERANCE_SP) THEN
            EVECTOR(1)=1.0_SP
            EVECTOR(2)=0.0_SP
          ELSE IF(ABS(EVALUE-A(2,2))<ZERO_TOLERANCE_SP) THEN
            EVECTOR(1)=0.0_SP
            EVECTOR(2)=1.0_DP
          ENDIF
        CASE(3)
          IF(ABS(A(1,2))<ZERO_TOLERANCE_SP.AND.ABS(A(1,3))<ZERO_TOLERANCE_SP.AND.ABS(A(2,3))<ZERO_TOLERANCE_SP) THEN
            EVECTOR=0.0_SP
            CALL FlagError("Zero matrix?? Eigenvectors undetermined",ERR,ERROR,*999)
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
          CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Evector is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("EIGENVECTOR_FULL_SP")
    RETURN
999 ERRORSEXITS("EIGENVECTOR_FULL_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EIGENVECTOR_FULL_SP

  !
  !================================================================================================================================
  !

  !>Returns the normalised eigenvector of a full double precision symmetric matrix A that corresponds to the eigenvalue EVALUE.
  SUBROUTINE EIGENVECTOR_FULL_DP(A,EVALUE,EVECTOR,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to find the eignevectors for
    REAL(DP), INTENT(IN) :: EVALUE !<The eigenvalue to find the eignevector for
    REAL(DP), INTENT(OUT) :: EVECTOR(:) !<On exit, the eigenvector corresponding the the eigenvalue
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,i1,i2,i3,ICYCLE(3,3)
    REAL(DP) :: AL,b(SIZE(A,1)),SUM,U(SIZE(A,1),SIZE(A,2)),x(SIZE(A,1))

    DATA ICYCLE /2,1,1,3,1,2,1,2,1/
    
    ENTERS("EIGENVECTOR_FULL_SP",ERR,ERROR,*999)

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
          ELSE IF(ABS(EVALUE-A(1,1))<ZERO_TOLERANCE_SP) THEN
            EVECTOR(1)=1.0_DP
            EVECTOR(2)=0.0_DP
          ELSE IF(ABS(EVALUE-A(2,2))<ZERO_TOLERANCE_SP) THEN
            EVECTOR(1)=0.0_DP
            EVECTOR(2)=1.0_DP
          ENDIF
        CASE(3)
          IF(ABS(A(1,2))<ZERO_TOLERANCE.AND.ABS(A(1,3))<ZERO_TOLERANCE.AND.ABS(A(2,3))<ZERO_TOLERANCE) THEN
            EVECTOR=0.0_DP
            CALL FlagError("Zero matrix?? Eigenvectors undetermined",ERR,ERROR,*999)
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
          CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Evector is too small",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("EIGENVECTOR_FULL_DP")
    RETURN
999 ERRORSEXITS("EIGENVECTOR_FULL_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE EIGENVECTOR_FULL_DP

  !
  !================================================================================================================================
  !
  
  !>Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION I0_DP(X)
      
    !Argument variables
    REAL(DP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the modified Bessel function of the first kind of order 0 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION I0_SP(X)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION I1_DP(X)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the modified Bessel function of the first kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION I1_SP(X)
  
   !Argument variables
    REAL(SP), INTENT(IN) :: X !<The value to evaluate the function at
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

  !>Returns an identity matrix
  SUBROUTINE IdentityMatrixSP(A,err,error,*)
    
    !Argument variables
    REAL(SP), INTENT(OUT) :: A(:,:) !<On exit, the identity matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("IdentityMatrixSP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      SELECT CASE(SIZE(A,1)) 
      CASE(1)
        A(1,1)=1.0_SP
      CASE(2)
        A(1,1)=1.0_SP
        A(2,1)=0.0_SP
        A(1,2)=0.0_SP
        A(2,2)=1.0_SP
      CASE(3)
        A(1,1)=1.0_SP
        A(2,1)=0.0_SP
        A(3,1)=0.0_SP
        A(1,2)=0.0_SP
        A(2,2)=1.0_SP
        A(3,2)=0.0_SP
        A(1,3)=0.0_SP
        A(2,3)=0.0_SP
        A(3,3)=1.0_SP
      CASE DEFAULT
        A=0.0_DP
        DO i=1,SIZE(A,1)
          A(i,i)=1.0_SP
        ENDDO !i
      END SELECT
    ELSE
      CALL FlagError("Matrix A is not square",err,error,*999)
    ENDIF

    EXITS("IdentityMatrixSP")
    RETURN
999 ERRORSEXITS("IdentityMatrixSP",err,error)
    RETURN 1
    
  END SUBROUTINE IdentityMatrixSP

  !
  !================================================================================================================================
  !

  !>Returns an identity matrix
  SUBROUTINE IdentityMatrixDP(A,err,error,*)
    
    !Argument variables
    REAL(DP), INTENT(OUT) :: A(:,:) !<On exit, the identity matrix
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    
    ENTERS("IdentityMatrixDP",err,error,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      SELECT CASE(SIZE(A,1)) 
      CASE(1)
        A(1,1)=1.0_DP
      CASE(2)
        A(1,1)=1.0_DP
        A(2,1)=0.0_DP
        A(1,2)=0.0_DP
        A(2,2)=1.0_DP
      CASE(3)
        A(1,1)=1.0_DP
        A(2,1)=0.0_DP
        A(3,1)=0.0_DP
        A(1,2)=0.0_DP
        A(2,2)=1.0_DP
        A(3,2)=0.0_DP
        A(1,3)=0.0_DP
        A(2,3)=0.0_DP
        A(3,3)=1.0_DP
      CASE DEFAULT
        A=0.0_DP
        DO i=1,SIZE(A,1)
          A(i,i)=1.0_DP
        ENDDO !i
      END SELECT
    ELSE
      CALL FlagError("Matrix A is not square",err,error,*999)
    ENDIF

    EXITS("IdentityMatrixDP")
    RETURN
999 ERRORSEXITS("IdentityMatrixDP",err,error)
    RETURN 1
    
  END SUBROUTINE IdentityMatrixDP

  !
  !================================================================================================================================
  !

  !>Inverts a full single precision matrix A to give matrix B and returns the determinant of A in DET.
  SUBROUTINE INVERT_FULL_SP(A,B,DET,ERR,ERROR,*)
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix to invert
    REAL(SP), INTENT(OUT) :: B(:,:) !<On exit, the inverse of A
    REAL(SP), INTENT(OUT) :: DET !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    ENTERS("INVERT_FULL_SP",ERR,ERROR,*999)

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
          CALL FlagError("Matrix size is not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Matrix B is not the same size as matrix A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix A is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("INVERT_FULL_SP")
    RETURN
999 ERRORSEXITS("INVERT_FULL_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INVERT_FULL_SP

  !
  !================================================================================================================================
  !

  !>Inverts a full double precision matrix A to give matrix B and returns the determinant of A in DET.
  SUBROUTINE INVERT_FULL_DP(A,B,DET,ERR,ERROR,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix A to invert
    REAL(DP), INTENT(OUT) :: B(:,:) !<On exit, the inverse of A
    REAL(DP), INTENT(OUT) :: DET !<On exit, the determinant of A
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

    ENTERS("INVERT_FULL_DP",ERR,ERROR,*999)

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
          CALL FlagError("Matrix size is not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Matrix B is not the same size as matrix A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix A is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("INVERT_FULL_DP")
    RETURN
999 ERRORSEXITS("INVERT_FULL_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE INVERT_FULL_DP
  
  !
  !================================================================================================================================
  !

  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION K0_DP(X)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the modified Bessel function of the second kind of order 0 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION K0_SP(X)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a double precision argument.
  PURE FUNCTION K1_DP(X)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the modified Bessel function of the second kind of order 1 using the approximation of Abromowitz and Stegun,
  !>for a single precision argument.
  PURE FUNCTION K1_SP(X)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the elliptic integral of the first kind - K(m), for a double precision argument.
  PURE FUNCTION KDP_DP(X)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Calculates the elliptic integral of the first kind - K(m), for a single precision argument.
  PURE FUNCTION KDP_SP(X)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: X !<The value to evaluate the function at
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
  
  !>Returns the L2-norm of the single precision vector A.
  PURE FUNCTION L2NORM_SP(A)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:) !<The vector to calculate the L2 norm of
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

  !>Returns the L2-norm of the double precision vector A.
  FUNCTION L2NORM_DP(A)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:) !<The vector to calculate the L2 norm of
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

  !>Calculates and returns the matrix-prouct of the single precision matrix A*B in C for single precision arguments.
  SUBROUTINE MATRIX_PRODUCT_SP(A,B,C,ERR,ERROR,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The first matrix A
    REAL(SP), INTENT(IN) :: B(:,:) !<The second matrix B
    REAL(SP), INTENT(OUT) :: C(:,:) !<On exit, the product matrix C=A*B
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("MATRIX_PRODUCT_SP",ERR,ERROR,*999)

    IF(SIZE(A,2)==SIZE(B,1).AND.SIZE(A,1)==SIZE(C,1).AND.SIZE(B,2)==SIZE(C,2)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)        
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
        C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      CASE DEFAULT
        CALL FlagError("Invalid matrix size.",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Invalid matrix sizes.",ERR,ERROR,*999)
    ENDIF

    EXITS("MATRIX_PRODUCT_SP")
    RETURN
999 ERRORSEXITS("MATRIX_PRODUCT_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MATRIX_PRODUCT_SP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the matrix-prouct of the double precision matrix A*B in C.
  SUBROUTINE MATRIX_PRODUCT_DP(A,B,C,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(IN) :: B(:,:) !<The B matrix
    REAL(DP), INTENT(OUT) :: C(:,:) !<On exit, the product matrix C=A*B
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
        
    ENTERS("MATRIX_PRODUCT_DP",ERR,ERROR,*999)
    
   IF(SIZE(A,2)==SIZE(B,1).AND.SIZE(A,1)==SIZE(C,1).AND.SIZE(B,2)==SIZE(C,2)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        C(1,1)=A(1,1)*B(1,1)
      CASE(2)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)
      CASE(3)
        C(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)
        C(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)
        C(1,3)=A(1,1)*B(1,3)+A(1,2)*B(2,3)+A(1,3)*B(3,3)        
        C(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)
        C(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)
        C(2,3)=A(2,1)*B(1,3)+A(2,2)*B(2,3)+A(2,3)*B(3,3)
        C(3,1)=A(3,1)*B(1,1)+A(3,2)*B(2,1)+A(3,3)*B(3,1)
        C(3,2)=A(3,1)*B(1,2)+A(3,2)*B(2,2)+A(3,3)*B(3,2)
        C(3,3)=A(3,1)*B(1,3)+A(3,2)*B(2,3)+A(3,3)*B(3,3)
      CASE DEFAULT
        CALL FlagError("Invalid matrix size.",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Invalid matrix sizes.",ERR,ERROR,*999)
    ENDIF

    EXITS("MATRIX_PRODUCT_DP")
    RETURN
999 ERRORSEXITS("MATRIX_PRODUCT_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MATRIX_PRODUCT_DP

  !
  !================================================================================================================================
  !

  !>Returns the transpose of a single precision matrix A in AT.
  SUBROUTINE MATRIX_TRANSPOSE_SP(A,AT,ERR,ERROR,*)
          
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The matrix to take the transpose of
    REAL(SP), INTENT(OUT) :: AT(:,:) !<On exit, the transpose of the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("MATRIX_TRANSPOSE_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(AT,2).AND.SIZE(A,2)==SIZE(AT,1)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        AT(1,1)=A(1,1)
      CASE(2)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
      CASE(3)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(1,3)=A(3,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
        AT(2,3)=A(3,2)
        AT(3,1)=A(1,3)
        AT(3,2)=A(2,3)
        AT(3,3)=A(3,3)
      CASE DEFAULT
        CALL FlagError("Invalid matrix size.",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Invalid matrix size.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("MATRIX_TRANSPOSE_SP")
    RETURN
999 ERRORSEXITS("MATRIX_TRANSPOSE_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MATRIX_TRANSPOSE_SP

  !
  !================================================================================================================================
  !

  !>Returns the transpose of a double precision matrix A in AT.
  SUBROUTINE MATRIX_TRANSPOSE_DP(A,AT,ERR,ERROR,*)
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The matrix to take the transpose of
    REAL(DP), INTENT(OUT) :: AT(:,:) !<On exit, the transpose of the matrix
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
        
    ENTERS("MATRIX_TRANSPOSE_DP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(AT,2).AND.SIZE(A,2)==SIZE(AT,1)) THEN
      SELECT CASE(SIZE(A,1))
      CASE(1)
        AT(1,1)=A(1,1)
      CASE(2)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
      CASE(3)
        AT(1,1)=A(1,1)
        AT(1,2)=A(2,1)
        AT(1,3)=A(3,1)
        AT(2,1)=A(1,2)
        AT(2,2)=A(2,2)
        AT(2,3)=A(3,2)
        AT(3,1)=A(1,3)
        AT(3,2)=A(2,3)
        AT(3,3)=A(3,3)
      CASE DEFAULT
        CALL FlagError("Invalid matrix size.",ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Invalid matrix size.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MATRIX_TRANSPOSE_DP")
    RETURN
999 ERRORSEXITS("MATRIX_TRANSPOSE_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MATRIX_TRANSPOSE_DP

  !
  !================================================================================================================================
  !
  
  !>Normalises a real single precision vector A.
  FUNCTION NORMALISE_SP(A,ERR,ERROR)

    !Argument variables
    REAL(SP), INTENT(IN) :: A(:) !<The vector to normalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(SP) :: NORMALISE_SP(SIZE(A,1))
    !Local variables
    REAL(SP) :: LENGTH
    
    ENTERS("NORMALISE_SP",ERR,ERROR,*999)

    LENGTH=L2NORM(A)
    IF(ABS(LENGTH)<ZERO_TOLERANCE_SP) THEN
        NORMALISE_SP=A
        CALL FlagError("Length of vector to normalise is zero",ERR,ERROR,*999)
    ELSE
        NORMALISE_SP=A/LENGTH
    ENDIF

    EXITS("NORMALISE_SP")
    RETURN
999 ERRORSEXITS("NORMALISE_SP",ERR,ERROR)
    RETURN
    
  END FUNCTION NORMALISE_SP

  !
  !================================================================================================================================
  !

  !>Normalises a real double precision vector A.
  FUNCTION NORMALISE_DP(A,ERR,ERROR)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(:) !<The vector to normalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: NORMALISE_DP(SIZE(A,1))
    !Local variables
    REAL(DP) :: LENGTH
    
    ENTERS("NORMALISE_DP",ERR,ERROR,*999)

    LENGTH=L2NORM(A)
    IF(ABS(LENGTH)<ZERO_TOLERANCE) THEN
      NORMALISE_DP=A
      CALL FlagError("Length of vector to normalise is zero",ERR,ERROR,*999)
    ELSE
      NORMALISE_DP=A/LENGTH
    ENDIF

    EXITS("NORMALISE_DP")
    RETURN
999 ERRORSEXITS("NORMALISE_DP",ERR,ERROR)
    RETURN
    
  END FUNCTION NORMALISE_DP

  !
  !================================================================================================================================
  !

  !>Calculates and returns the normalised vector cross-prouct of the single precision vectors A*B in C.
  SUBROUTINE NORM_CROSS_PRODUCT_SP(A,B,C,ERR,ERROR,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:) !<The first vector in the cross product
    REAL(SP), INTENT(IN) :: B(:) !<The second vector in the cross product
    REAL(SP), INTENT(OUT) :: C(:) !<On exit, the normalised cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("NORM_CROSS_PRODUCT_SP",ERR,ERROR,*999)

    CALL CROSS_PRODUCT(A,B,C,ERR,ERROR,*999)
    C=NORMALISE(C,ERR,ERROR)
    IF(ERR/=0) GOTO 999

    EXITS("NORM_CROSS_PRODUCT_SP")
    RETURN
999 ERRORSEXITS("NORM_CROSS_PRODUCT_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NORM_CROSS_PRODUCT_SP
  
  !
  !================================================================================================================================
  !

  !>Calculates and returns the normalised vector cross-prouct of the double precision vectors A*B in C.
  SUBROUTINE NORM_CROSS_PRODUCT_DP(A,B,C,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:) !<The first vector in the cross product
    REAL(DP), INTENT(IN) :: B(:) !<The second vector in the cross product
    REAL(DP), INTENT(OUT) :: C(:) !<On exit, the normalised cross product of the first and second vectors
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    ENTERS("NORM_CROSS_PRODUCT_DP",ERR,ERROR,*999)

    CALL CROSS_PRODUCT(A,B,C,ERR,ERROR,*999)
    C=NORMALISE(C,ERR,ERROR)
    IF(ERR/=0) GOTO 999

    EXITS("NORM_CROSS_PRODUCT_DP")
    RETURN
999 ERRORSEXITS("NORM_CROSS_PRODUCT_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NORM_CROSS_PRODUCT_DP
  
  !
  !================================================================================================================================
  !

  !>Finds the solution to a small single precision linear system Ax=b.
  SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_SP(A,x,b,ERR,ERROR,*)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(SP), INTENT(OUT) :: x(:) !<On exit, the solution vector x
    REAL(SP), INTENT(IN) :: b(:) !<The RHS vector b
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(SP) :: AINV(SIZE(A,1),SIZE(A,2)),ADET
    
    ENTERS("SOLVE_SMALL_LINEAR_SYSTEM_SP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)==SIZE(b,1)) THEN
        IF(SIZE(A,1)<=SIZE(x,1)) THEN
          SELECT CASE(SIZE(A,1))
          CASE(1:3)
            CALL INVERT(A,AINV,ADET,ERR,ERROR,*999)
            x=MATMUL(AINV,b)
          CASE DEFAULT
            CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("x is too small",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Size of b is not the same as the number of rows in A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("SOLVE_SMALL_LINEAR_SYSTEM_SP")
    RETURN
999 ERRORSEXITS("SOLVE_SMALL_LINEAR_SYSTEM_SP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_SP

  !
  !================================================================================================================================
  !

  !>Finds the solution to a small double precision linear system Ax=b.
  SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_DP(A,x,b,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:,:) !<The A matrix
    REAL(DP), INTENT(OUT) :: x(:) !<On exit, the solution vector x
    REAL(DP), INTENT(IN) :: b(:) !<The RHS vector b
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: AINV(SIZE(A,1),SIZE(A,2)),ADET
    
    ENTERS("SOLVE_SMALL_LINEAR_SYSTEM_DP",ERR,ERROR,*999)

    IF(SIZE(A,1)==SIZE(A,2)) THEN
      IF(SIZE(A,1)==SIZE(b,1)) THEN
        IF(SIZE(A,1)<=SIZE(x,1)) THEN
          SELECT CASE(SIZE(A,1))
          CASE(1:3)
            CALL INVERT(A,AINV,ADET,ERR,ERROR,*999)
            x=MATMUL(AINV,b)
          CASE DEFAULT
            CALL FlagError("Matrix size not implemented",ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("x is too small",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Size of b is not the same as the number of rows in A",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Matrix is not square",ERR,ERROR,*999)
    ENDIF

    EXITS("SOLVE_SMALL_LINEAR_SYSTEM_DP")
    RETURN
999 ERRORSEXITS("SOLVE_SMALL_LINEAR_SYSTEM_DP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE SOLVE_SMALL_LINEAR_SYSTEM_DP

  !
  !================================================================================================================================
  !
  
  !> Calculates single precision hyperbolic cotangent function
  FUNCTION COTH_SP(A)

    !Argument variables
    REAL(SP), INTENT(IN) :: A !<argument to perform coth() on
    !Function variable
    REAL(SP) :: COTH_SP
    
    COTH_SP=(EXP(A)+EXP(-1.0_SP*A))/(EXP(A)-EXP(-1.0_SP*A))

    RETURN
    
  END FUNCTION COTH_SP

  !
  !================================================================================================================================
  !

  !> Calculates double precision hyperbolic cotangent function
  FUNCTION COTH_DP(A)

    !Argument variables
    REAL(DP), INTENT(IN) :: A !<argument to perform coth() on
    !Function variable
    REAL(DP) :: COTH_DP

    COTH_DP=(EXP(A)+EXP(-1.0_DP*A))/(EXP(A)-EXP(-1.0_DP*A))

    RETURN
    
  END FUNCTION COTH_DP

  !
  !================================================================================================================================
  !

  !> Calculates second derivatives of a cubic spline function for a tabulated function y(x). Call  spline_cubic_val to evaluate at t values.
  !> algorithm adapted from John Burkhardt's spline_cubic_set routine from the SPLINE package (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)
  SUBROUTINE spline_cubic_set (n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp, err, error, *)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !< size of x,y arrays to interpolate values from
    REAL(DP), INTENT(IN) :: t(n) !< t array: known values
    REAL(DP), INTENT(IN) :: y(n) !< y array: values to interpolate
    INTEGER(INTG), INTENT(IN) :: ibcbeg !< left boundary condition flag
    REAL(DP), INTENT(IN) :: ybcbeg !< 1st derivative interpolating function at point 1 (left boundary)
    INTEGER(INTG), INTENT(IN) :: ibcend !< right boundary condition flag
    REAL(DP), INTENT(IN) :: ybcend !< 1st derivative interpolating function at point n (right boundary)
    REAL(DP), INTENT(OUT) :: ypp(n) !< 2nd derivatives of interpolating function at x values
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: diag(n)
    REAL(DP) :: sub(2:n)
    REAL(DP) :: sup(1:n-1)
    INTEGER(INTG) :: i
    TYPE(VARYING_STRING) :: localError

    ENTERS("spline_cubic_set",ERR,ERROR,*999)

    ! Sanity checks
    IF ( n <= 1 ) then
      localError="spline interpolation requires at least 2 knots- user supplied "//TRIM(NumberToVString(n,"*",ERR,ERROR))
      CALL FlagError(localError,ERR,ERROR,*999)
    ENDIF
    DO i = 1, n-1
      IF ( t(i) >= t(i+1) ) then
        localError="Non-increasing knots supplied for cubic spline interpolation."
        CALL FlagError(localError,ERR,ERROR,*999)
      ENDIF
    ENDDO

    !  Set the first equation.
    IF ( ibcbeg == 0 ) then
      ypp(1) = 0.0E+00_DP
      diag(1) = 1.0E+00_DP
      sup(1) = -1.0E+00_DP
    ELSE IF ( ibcbeg == 1 ) then
      ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
      diag(1) = ( t(2) - t(1) ) / 3.0E+00_DP
      sup(1) = ( t(2) - t(1) ) / 6.0E+00_DP
    ELSE IF ( ibcbeg == 2 ) then
      ypp(1) = ybcbeg
      diag(1) = 1.0E+00_DP
      sup(1) = 0.0E+00_DP
    ELSE
      localError="The boundary flag IBCBEG must be 0, 1 or 2."
      CALL FlagError(localError,ERR,ERROR,*999)
    ENDIF

    !  Set the intermediate equations.
    DO i = 2, n-1
      ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
       &     - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
      sub(i) = ( t(i) - t(i-1) ) / 6.0E+00_DP
      diag(i) = ( t(i+1) - t(i-1) ) / 3.0E+00_DP
      sup(i) = ( t(i+1) - t(i) ) / 6.0E+00_DP
    ENDDO

    !  Set the last equation.
    IF ( ibcend == 0 ) then
      ypp(n) = 0.0E+00_DP
      sub(n) = -1.0E+00_DP
      diag(n) = 1.0E+00_DP
    ELSE IF ( ibcend == 1 ) then
      ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
      sub(n) = ( t(n) - t(n-1) ) / 6.0E+00_DP
      diag(n) = ( t(n) - t(n-1) ) / 3.0E+00_DP
    ELSE IF ( ibcend == 2 ) then
      ypp(n) = ybcend
      sub(n) = 0.0E+00_DP
      diag(n) = 1.0E+00_DP
    ELSE
      localError="The boundary flag IBCEND must be 0, 1 or 2."
      CALL FlagError(localError,ERR,ERROR,*999)
    ENDIF

    !  Special case:
    !    N = 2, IBCBEG = IBCEND = 0.
    IF ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then
      ypp(1) = 0.0E+00_DP
      ypp(2) = 0.0E+00_DP

    !  Solve the linear system.
    ELSE
      CALL s3_fs ( sub, diag, sup, n, ypp, ypp, err, error, *999 )
    ENDIF

    EXITS("spline_cubic_set")
    RETURN
999 ERRORSEXITS("spline_cubic_set",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE spline_cubic_set

  !
  !================================================================================================================================
  !

  !> S3_FS factors and solves a tridiagonal linear system.
  !> algorithm adapted from John Burkhardt's s3_fs routine from the SPLINE package (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)
  SUBROUTINE s3_fs ( a1, a2, a3, n, b, x, err, error, *)

    !Argument variables
    REAL(DP), INTENT(INOUT) :: a1(2:n) !< IN: nonzero diagonal of linear system OUT: factorization info
    REAL(DP), INTENT(INOUT) :: a2(1:n) !< IN: nonzero diagonal of linear system OUT: factorization info
    REAL(DP), INTENT(INOUT) :: a3(1:n-1) !< IN: nonzero diagonal of linear system OUT: factorization info
    INTEGER(INTG), INTENT(IN) :: n !< size of x,y arrays to interpolate values from
    REAL(DP), INTENT(INOUT) :: b(n) !< IN: RHS of linear system OUT: factorization info
    REAL(DP), INTENT(OUT) :: x(n) !< solution of linear system
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    REAL(DP) :: xmult
    TYPE(VARYING_STRING) :: localError

    ENTERS("s3_fs",ERR,ERROR,*999)

    !  The diagonal entries can't be zero.
    DO i = 1, n
      IF ( ABS(a2(i)) < ZERO_TOLERANCE ) then
        localError="Zero diagonal entry in tridiagonal linear system."
        CALL FlagError(localError,ERR,ERROR,*999)
      ENDIF
    ENDDO

    DO i = 2, n-1
      xmult = a1(i) / a2(i-1)
      a2(i) = a2(i) - xmult * a3(i-1)
      b(i) = b(i) - xmult * b(i-1)
    ENDDO

    xmult = a1(n) / a2(n-1)
    a2(n) = a2(n) - xmult * a3(n-1)
    x(n) = ( b(n) - xmult * b(n-1) ) / a2(n)
    DO i = n-1, 1, -1
      x(i) = ( b(i) - a3(i) * x(i+1) ) / a2(i)
    ENDDO

    EXITS("s3_fs")
    RETURN
999 ERRORSEXITS("s3_fs",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE s3_fs

  !
  !================================================================================================================================
  !
  
  !> Evaluates a cubic spline at a specified point. First call spline_cubic_set to calculate derivatives
  !> algorithm adapted from John Burkhardt's spline_cubic_val routine from the SPLINE package (http://people.sc.fsu.edu/~jburkardt/f_src/spline/spline.html)
  SUBROUTINE spline_cubic_val (n, t, y, ypp, tval, yval, ypval, yppval, err, error, *)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: n !< size of t,y arrays to interpolate values from
    REAL(DP), INTENT(IN) :: t(n) !< t array: known knot values
    REAL(DP), INTENT(IN) :: y(n) !< y array: data values to interpolate at the knots
    REAL(DP), INTENT(IN) :: ypp(n) !< 2nd derivatives of interpolating function at t values
    REAL(DP), INTENT(IN) :: tval !< point in t at which spline is to be evaluated
    REAL(DP), INTENT(OUT) :: yval !< spline interpolated y value at tval
    REAL(DP), INTENT(OUT) :: ypval !< first derivative of spline interpolated y value at tval
    REAL(DP), INTENT(OUT) :: yppval !< second derivative of spline interpolated y value at tval
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: dt
    REAL(DP) :: h
    INTEGER(INTG) :: i
    INTEGER(INTG) :: left
    INTEGER(INTG) :: right
    LOGICAL :: foundInterval

    ENTERS("spline_cubic_val",ERR,ERROR,*999)

    !  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
    !  Values below T(1) or above T(N) use extrapolation.
    foundInterval = .FALSE.
    DO i = 2, n - 1
      IF ( tval < t(i) ) THEN
        foundInterval=.TRUE.
        left = i - 1
        right = i
        EXIT
      ENDIF
    ENDDO
    IF (foundInterval .EQV. .FALSE.) THEN
      left = n - 1
      right = n
    ENDIF

    !  Evaluate the polynomial.
    dt = tval - t(left)
    h = t(right) - t(left)

    yval = y(left) &
     &   + dt * ( ( y(right) - y(left) ) / h &
     &          - ( ypp(right) / 6.0E+00_DP + ypp(left) / 3.0E+00_DP ) * h &
     &   + dt * ( 0.5E+00 * ypp(left) &
     &   + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0E+00_DP * h ) ) ) )

    ypval = ( y(right) - y(left) ) / h &
     &   - ( ypp(right) / 6.0E+00_DP + ypp(left) / 3.0E+00_DP ) * h &
     &   + dt * ( ypp(left) &
     &   + dt * ( 0.5E+00_DP * ( ypp(right) - ypp(left) ) / h ) )

    yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h 

    EXITS("spline_cubic_val")
    RETURN
999 ERRORSEXITS("spline_cubic_val",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE spline_cubic_val

  !
  !================================================================================================================================
  !

END MODULE MATHS

