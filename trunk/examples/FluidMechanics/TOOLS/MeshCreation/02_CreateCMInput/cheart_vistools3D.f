      SUBROUTINE GET_SPATIAL_FIELD_VALUE(MM,AA,ID,ELEM,XI,Pt)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) MM
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I, J, K, ELEM, ORD1(3), ID, NID
        DOUBLE PRECISION XI(3), Pt(AA%B(ID)%DM), Ax
        Pt(1:AA%B(ID)%DM) = 0_DPR; ORD1(1:3) = 0; NID = 0
        IF (AA%SPL .NE. ID) NID = 3
        DO K = 1,AA%B(ID)%n
           CALL POLY(AA,ID,K,XI,ORD1,Ax)
           Pt(1:AA%B(ID)%DM) = Pt(1:AA%B(ID)%DM) + MM%X(MM%T(ELEM,K),NID+1:NID+AA%B(ID)%DM)*Ax
        END DO
      END SUBROUTINE GET_SPATIAL_FIELD_VALUE


      SUBROUTINE GET_GREEN_STRAIN(DOM,MM,AA,ELEM,XI,GRN)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) DOM, MM
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I, J, K, ELEM, ORD1(3), ID, n
        DOUBLE PRECISION XI(3), VM(3,3), GRN(6), Ax
        CALL GET_FIELD_GRADIENT_XO(DOM,MM,AA,AA%VEL,ELEM,XI,VM)
!         PRINT *, VM(1,1:3)
!         PRINT *, VM(2,1:3)
!         PRINT *, VM(3,1:3)
        n = 0
        DO I = 1,3
          DO J = I,3
            n = n + 1
            GRN(n) = VM(1,I)*VM(1,J) + VM(2,I)*VM(2,J) + VM(3,I)*VM(3,J) 
          END DO
        END DO
!         PRINT *, GRN(1:6)
!         READ *, I
      END SUBROUTINE GET_GREEN_STRAIN


      SUBROUTINE GET_FIELD_GRADIENT_XI(MM,AA,ID,ELEM,XI,M)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) MM
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I, J, K, ELEM, ORD1(3), ID, NID
        DOUBLE PRECISION XI(3), M(AA%B(ID)%DM,3), Ax
        M(1:AA%B(ID)%DM,1:3) = 0_DPR; NID = 0
        IF (AA%SPL .NE. ID) NID = 3
!         PRINT *, 'ELEM', ELEM
        DO J = 1,AA%B(ID)%n
!            PRINT *, J, MM%T(ELEM,J)
           DO K = 1,3
              ORD1(1:3) = 0; ORD1(K) = 1
              CALL POLY(AA,ID,J,XI,ORD1,Ax)
              M(1:AA%B(ID)%DM,K) = M(1:AA%B(ID)%DM,K) + MM%X(MM%T(ELEM,J),NID+1:NID+AA%B(ID)%DM)*Ax
           END DO
        END DO
      END SUBROUTINE GET_FIELD_GRADIENT_XI


      SUBROUTINE GET_FIELD_GRADIENT_XO(DOM,MM,AA,ID,ELEM,XI,VM)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) DOM, MM
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I, J, K, ELEM, ORD1(3), ID
        DOUBLE PRECISION XI(3), VM(AA%B(ID)%DM,3), VXI(AA%B(ID)%DM,3), M(3,3), invM(3,3), Ax
        PRINT *, '1her', ELEM, ID
        CALL GET_FIELD_GRADIENT_XI(MM, AA, ID, ELEM, XI, VXI)
!         PRINT *, '2her', VXI(1,1:3)
!         PRINT *, '2her', VXI(2,1:3)
!         PRINT *, '2her', VXI(3,1:3)
!         PRINT *, '2her', ELEM
        CALL GET_FIELD_GRADIENT_XI(DOM,AA,AA%SPL,ELEM,XI,M)
!         PRINT *, '3her', M(1,1:3)
!         PRINT *, '3her', M(2,1:3)
!         PRINT *, '3her', M(3,1:3)
!         PRINT *, '3her', ELEM
        CALL TENSOR_INVERSE(M,invM,K)
!         PRINT *, '4her', invM(1,1:3), K
!         PRINT *, '4her', invM(2,1:3)
!         PRINT *, '4her', invM(3,1:3)
        IF ( K == 1 ) STOP
!         PRINT *, '4her', AA%B(ID)%DM
        CALL MATRIX_MATRIX_PRODUCT(VXI,invM,AA%B(ID)%DM,3,3,VM)
!         READ *, I
      END SUBROUTINE GET_FIELD_GRADIENT_XO


      SUBROUTINE POLY(AA,ID,BAS,X,ORD1,VAL)
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER n,ORD1(3),I,J,P,BAS,K,ID
        DOUBLE PRECISION VAL, X(3)
        IF (AA%HEXA == 1) THEN
           CALL POLY_HEX(AA,ID,BAS,X,ORD1,VAL)
        ELSE
           CALL POLY_TET(AA,ID,BAS,X,ORD1,VAL)
        END IF
      END SUBROUTINE POLY


      SUBROUTINE POLY_HEX(AA,ID,BAS,X,ORD1,VAL)
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER n,ORD1(3),I,J,P,BAS,K,ID
        DOUBLE PRECISION VAL,Ax,Bx,Cx,Bo,X(3),Co(3)
        Co(1:3) = 0_DPR
        DO I = 1,AA%B(ID)%nl
           DO J = 1,3
              IF (ORD1(J) <= AA%B(ID)%P(1,I)) THEN
                 Ax = 1_DPR
                 DO K = 0,ORD1(J)-1
                    Ax = Ax*(REAL(AA%B(ID)%P(1,I))-REAL(K))
                 END DO
                 Ax = Ax*(X(J)**(AA%B(ID)%P(1,I)-ORD1(J)))
                 Co(J) = Co(J) + AA%B(ID)%Q(AA%B(ID)%B_ID(BAS,J),I)*Ax
              END IF
           END DO
        END DO
        VAL = Co(1)*Co(2)*Co(3)
      END SUBROUTINE POLY_HEX


      SUBROUTINE POLY_TET(AA,ID,BAS,X,ORD1,VAL)
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER n,ORD1(3),I,J,P,BAS,K,ID
        DOUBLE PRECISION VAL,Ax,Bx,Cx,Bo,X(3),Co(3)
        VAL = 0_DPR
        DO I = 1,AA%B(ID)%n
           Co(1:3) = 0_DPR
           DO J = 1,3
             IF (ORD1(J) <= AA%B(ID)%P(J,I)) THEN
              Ax = 1_DPR
              DO K = 0,ORD1(J)-1
                 Ax = Ax*(REAL(AA%B(ID)%P(J,I))-REAL(K))
              END DO
              Ax = Ax*(X(J)**(AA%B(ID)%P(J,I)-ORD1(J)))
              Co(J) = Co(J) + Ax
             END IF
           END DO
           VAL = VAL + AA%B(ID)%Q(BAS,I)*Co(1)*Co(2)*Co(3)
        END DO
      END SUBROUTINE POLY_TET


      SUBROUTINE CRSS(P1,P2,P3)
        DOUBLE PRECISION P1(3), P2(3), P3(3)
        P3(1) = P1(2)*P2(3)-P1(3)*P2(2)
        P3(2) = P1(3)*P2(1)-P1(1)*P2(3)
        P3(3) = P1(1)*P2(2)-P1(2)*P2(1)
      END SUBROUTINE CRSS


      SUBROUTINE DOT(P1,P2,n,VAL)
        INTEGER I,n
        INTEGER, PARAMETER :: SPR=SELECTED_REAL_KIND(6,15)
        INTEGER, PARAMETER :: DPR=SELECTED_REAL_KIND(15,307)
        DOUBLE PRECISION P1(n), P2(n), VAL
        VAL = 0_DPR
        DO I = 1,n
           VAL = VAL + P1(I)*P2(I)
        END DO
      END SUBROUTINE DOT


      SUBROUTINE MATRIX_MATRIX_PRODUCT(A,B,n,m,p,C)
        INTEGER I,J,K,n,m,p
        INTEGER, PARAMETER :: SPR=SELECTED_REAL_KIND(6,15)
        INTEGER, PARAMETER :: DPR=SELECTED_REAL_KIND(15,307)
        DOUBLE PRECISION Ax, A(n,m), B(m,p), C(n,p)
        DO I = 1,n
           DO J = 1,p
              Ax = 0_DPR
              DO K = 1,m
                 Ax = Ax + A(I,K)*B(K,J)
              END DO
              C(I,J) = Ax
           END DO
        END DO
      END SUBROUTINE MATRIX_MATRIX_PRODUCT


      SUBROUTINE MATRIX_VECTOR_PRODUCT(A,V,n,m,B)
        INTEGER I,J,K,n,m
        INTEGER, PARAMETER :: SPR=SELECTED_REAL_KIND(6,15)
        INTEGER, PARAMETER :: DPR=SELECTED_REAL_KIND(15,307)
        DOUBLE PRECISION Ax, A(n,m), V(m), B(n)
        DO I = 1,n
           Ax = 0_DPR
           DO K = 1,m
              Ax = Ax + A(I,K)*V(K)
           END DO
           B(I) = Ax
        END DO
      END SUBROUTINE MATRIX_VECTOR_PRODUCT


      SUBROUTINE HADAMARD_PRODUCT(A,B,n,VAL)
        INTEGER I,J,K,n,m
        INTEGER, PARAMETER :: SPR=SELECTED_REAL_KIND(6,15)
        INTEGER, PARAMETER :: DPR=SELECTED_REAL_KIND(15,307)
        DOUBLE PRECISION Ax, A(n,n), B(n,n), VAL
        VAL = 0_DPR
        DO I = 1,n
           Ax = 0_DPR
           DO K = 1,n
              Ax = Ax + A(I,K)*B(I,K)
           END DO
           VAL = VAL + Ax
        END DO
      END SUBROUTINE HADAMARD_PRODUCT


      SUBROUTINE P_NORM(A,n,p,dL)
        INTEGER n,I,p
        INTEGER, PARAMETER :: SPR=SELECTED_REAL_KIND(6,15)
        INTEGER, PARAMETER :: DPR=SELECTED_REAL_KIND(15,307)
        DOUBLE PRECISION A(n), dL
        dL = 0
        IF (p.NE.0) THEN
           DO I = 1,n
           dL = dL + (abs(A(I)))**REAL(p)
           END DO
           dL = dL**(1_DPR / REAL(p))
        ELSE
           DO I = 1,n
              dL = MAX(dL,abs(A(I)))
           END DO
        END IF
      END SUBROUTINE P_NORM


      SUBROUTINE P_NORMALIZE(A,n,p,dL)
        INTEGER n,I
        DOUBLE PRECISION A(n), p, dL
        CALL P_NORM(A,n,p,dL)
        A(1:n) = A(1:n)/dL
      END SUBROUTINE P_NORMALIZE


      SUBROUTINE TENSOR_INVERSE(M,invM,ID)
!       P_NORM - Written by David A. Nordsletten (C) 2008
!       This subroutine calculates the 'little' l_p norms for a vector A(n).
!
!       VARIABLES ::  I, J, K, n (INTEGER), dummy variables
        USE CMGUI_VARS
        INTEGER I,J,ID
        DOUBLE PRECISION M(3,3), invM(3,3), Ax
        ID = 0
        invM(1,1) = M(2,2)*M(3,3) - M(2,3)*M(3,2)
        invM(2,1) = M(2,3)*M(3,1) - M(2,1)*M(3,3)
        invM(3,1) = M(2,1)*M(3,2) - M(2,2)*M(3,1)
        invM(1,2) = M(3,2)*M(1,3) - M(3,3)*M(1,2)
        invM(2,2) = M(3,3)*M(1,1) - M(3,1)*M(1,3)
        invM(3,2) = M(3,1)*M(1,2) - M(3,2)*M(1,1)
        invM(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
        invM(2,3) = M(1,3)*M(2,1) - M(1,1)*M(2,3)
        invM(3,3) = M(1,1)*M(2,2) - M(1,2)*M(2,1)
        Ax = M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1) &
             + M(1,3)*M(2,1)*M(3,2)- (M(1,3)*M(2,2)*M(3,1) &
             + M(1,2)*M(2,1)*M(3,3)+M(1,1)*M(2,3)*M(3,2))
        IF (abs(Ax) > GLOBAL_TOL) THEN
          invM = invM/Ax
        ELSE
          ID = 1
        END IF
      END SUBROUTINE TENSOR_INVERSE




