      SUBROUTINE LOAD_BASE_SET(AA)
!       Subprogram LOAD_BASE_SET - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program reads in all input basis for base AA.
!
        USE CMGUI_VARS
        INTEGER I,STP,D,VFLD
        CHARACTER*60 IN_CHAR
        TYPE(ARRAY_PROBLEM_BASE) AA
        NIMZ = TRIM(NIMZ); AA%n_B = 0; AA%HEXA = 0; AA%DM = 3
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')	! Read base file for initial parameters
          DO WHILE (0 < 1)
              READ(1,*,END=50) IN_CHAR
              IF (INDEX(IN_CHAR,'no_fields!') == 1)         READ(1,*) AA%n_B
              IF (INDEX(IN_CHAR,'no_gauss!') == 1)          READ(1,*) AA%n_pts
              IF (INDEX(IN_CHAR,'volume!') == 1)            READ(1,*) AA%VL
              IF (INDEX(IN_CHAR,'no_gauss_f!') == 1)        READ(1,*) AA%n_pts_f
              IF (INDEX(IN_CHAR,'no_ele_faces!') == 1)      READ(1,*) AA%FACES
              IF (INDEX(IN_CHAR,'no_ele_nodes_f!') == 1)    READ(1,*) AA%FNODES
              IF (INDEX(IN_CHAR,'hexa_basis!') == 1)        AA%HEXA = 1
              IF (INDEX(IN_CHAR,'domain_dimension!') == 1)  READ(1,*) AA%DM
              IF (INDEX(IN_CHAR,'TRI_BASIS!') == 1)  AA%TRI_BASIS = 1
              IF (INDEX(IN_CHAR,'TET_BASIS!') == 1)  AA%TET_BASIS = 1
              IF (INDEX(IN_CHAR,'QUAD_BASIS!') == 1) AA%QUAD_BASIS = 1
              IF (INDEX(IN_CHAR,'HEX_BASIS!') == 1)  AA%HEX_BASIS = 1
            END DO
 50     CLOSE(1)
        IF (AA%n_B == 0) CALL ERRORZ(3); ALLOCATE(AA%B(AA%n_B))		! Allocate base size
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')		! Read call letters
          DO WHILE (0 < 1)
              READ(1,*,END=55) IN_CHAR
              IF (INDEX(IN_CHAR,'call_letters!') == 1)       READ(1,*) AA%B(1:AA%n_B)%CL
          END DO
 55     CLOSE(1)
        AA%B(1:AA%n_B)%DISCONT = 0
        AA%n_ptsl = AA%n_pts
        AA%n_pts_fl = AA%n_pts_f
        IF (AA%HEXA == 0) THEN			! If tensor input is indicated
           D = 3
           ALLOCATE(AA%gpt_f(AA%n_pts_f,D,AA%FACES),AA%gw_f(AA%n_pts_f))
           ALLOCATE(AA%gpt(AA%n_pts,D),AA%gw(AA%n_pts))
        ELSE					! If full input is indicated
           D = 1
           AA%n_pts_f = AA%n_pts**REAL(AA%DM-1); AA%n_pts = AA%n_pts**REAL(AA%DM)
           ALLOCATE(AA%gpt(AA%n_pts,3),AA%gw(AA%n_pts))
           ALLOCATE(AA%gpt_f(AA%n_pts_f,3,AA%FACES),AA%gw_f(AA%n_pts_f))
        END IF
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')	! Read quadrature rule
          DO WHILE (0 < 1)
              READ(1,*,END=65) IN_CHAR
              IF (INDEX(IN_CHAR,'surface_area!') == 1)  READ(1,*) AA%VL_f(1:AA%FACES)
              IF (INDEX(IN_CHAR,'gauss_points!') == 1) THEN
                 DO I = 1,AA%n_ptsl
                   READ(1,*) AA%gpt(I,1:D)
                 END DO
              ELSE IF (INDEX(IN_CHAR,'gauss_weights!') == 1) THEN
                 DO I = 1,AA%n_ptsl
                   READ(1,*) AA%gw(I)
                 END DO
              ELSE IF (INDEX(IN_CHAR,'gauss_points_f!') == 1) THEN
                DO I = 1,AA%n_pts_fl
                  READ(1,*) AA%gpt_f(I,1:3,1)
                END DO
              ELSE IF (INDEX(IN_CHAR,'gauss_weights_f!') == 1) THEN
                DO I = 1,AA%n_pts_fl
                  READ(1,*) AA%gw_f(I)
                END DO
              END IF
            END DO
 65     CLOSE(1)
        IF (AA%HEXA == 1) CALL BUILD_HEXA_QUADRATURE(AA)	! Build Tensor quadrature rule
        CALL BUILD_FACE_GPT_ARRAY(AA)				! Build facet quadrature rule
        DO VFLD = 1,AA%n_B		! Loop over number of basis
        IF (INDEX(AA%B(VFLD)%CL,'M') >= 1) AA%SPL = VFLD
        IF (INDEX(AA%B(VFLD)%CL,'V') >= 1) AA%VEL = VFLD
        IF (INDEX(AA%B(VFLD)%CL,'P') >= 1) AA%PRS = VFLD
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')	! Read dimensions of basis CL
          DO WHILE (0 < 1)
              READ(1,*,END=70) IN_CHAR
              IF (INDEX(IN_CHAR,'no_basis_'//TRIM(AA%B(VFLD)%CL)//'!') == 1) READ(1,*) AA%B(VFLD)%n
              IF (INDEX(IN_CHAR,'dim_field_'//TRIM(AA%B(VFLD)%CL)//'!') == 1) READ(1,*) AA%B(VFLD)%DM
            END DO
 70     CLOSE(1)
        AA%B(VFLD)%nl = AA%B(VFLD)%n
        IF (AA%HEXA == 0) THEN		! If tensor input is indicated
           D = 3
           ALLOCATE(AA%B(VFLD)%Q(AA%B(VFLD)%n,AA%B(VFLD)%n),AA%B(VFLD)%P(3,AA%B(VFLD)%n))
           ALLOCATE( AA%B(VFLD)%XI(AA%B(VFLD)%n,3) )
        ELSE				! If full input is indicated
           D = 1
           AA%B(VFLD)%n  = AA%B(VFLD)%n**REAL(AA%DM)
           ALLOCATE(AA%B(VFLD)%Q(AA%B(VFLD)%nl,AA%B(VFLD)%nl),AA%B(VFLD)%P(D,AA%B(VFLD)%nl))
           ALLOCATE( AA%B(VFLD)%B_ID(AA%B(VFLD)%n,3), AA%B(VFLD)%XI(AA%B(VFLD)%n,3) )
        END IF
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')		! Read in Basis
          DO WHILE (0 < 1)
              READ(1,*,END=80) IN_CHAR
              IF (INDEX(IN_CHAR,TRIM(AA%B(VFLD)%CL)//'_basis!') == 1) THEN
                 DO I = 1,AA%B(VFLD)%nl
                 READ(1,*) AA%B(VFLD)%Q(I,1:AA%B(VFLD)%nl)
                 END DO
              ELSE IF (INDEX(IN_CHAR,TRIM(AA%B(VFLD)%CL)//'_basis_discontinuous!') == 1) THEN
                 AA%B(VFLD)%DISCONT = 1
              ELSE IF (INDEX(IN_CHAR,TRIM(AA%B(VFLD)%CL)//'_xi_coordinates!') == 1) THEN
                 DO I = 1,AA%B(VFLD)%nl
                 READ(1,*) AA%B(VFLD)%XI(I,1:D)
                 END DO
              ELSE IF (INDEX(IN_CHAR,TRIM(AA%B(VFLD)%CL)//'_P_basis!') == 1) THEN
                 DO I = 1,D
                 READ(1,*) AA%B(VFLD)%P(I,1:AA%B(VFLD)%nl)
                 END DO
              ELSE IF (INDEX(IN_CHAR,'basis_ordering_'//TRIM(AA%B(VFLD)%CL)//'!') == 1) THEN
                DO I = 1,AA%B(VFLD)%n
                  READ(1,*) AA%B(VFLD)%B_ID(I,1:3)
                END DO
              END IF
            END DO
 80     CLOSE(1)
        IF (AA%HEXA == 1) CALL BUILD_HEXA_XI(AA,VFLD)		! Build tensor xi points
        CALL LOAD_INTEGRALS(AA,VFLD)				! Calculate basis integrals
        PRINT *, ' ==== >>>   BASIS FIELD  <<< ==== '
        PRINT *, '   CALL LETTER     :: ', AA%B(VFLD)%CL
        PRINT *, '   FIELD DIMENSION :: ', AA%B(VFLD)%DM
        PRINT *, '   DIM  (S ( Tm )) :: ', AA%B(VFLD)%n

        END DO
        RETURN
 90     CALL ERRORZ(5)
      END SUBROUTINE LOAD_BASE_SET


      SUBROUTINE LOAD_INTEGRALS(AA,VFLD)
!       Subprogram LOAD_INTEGRALS - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program loads integrals for base AA%B(VFLD).  The basis is evaluated a given quadrature points.
!
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K,x,y,z,ID,ORD1(3),VFLD,S
        DOUBLE PRECISION Ax(3),Bx,Cx(3),Dx
        ALLOCATE( AA%B(VFLD)%I%Y(AA%B(VFLD)%n,AA%n_pts) )
        ALLOCATE( AA%B(VFLD)%I%Y_f(AA%B(VFLD)%n,AA%n_pts_f,AA%FACES) )
        ALLOCATE( AA%B(VFLD)%I%dY(3,AA%B(VFLD)%n,AA%n_pts) )
        ALLOCATE( AA%B(VFLD)%I%TdY(3,AA%B(VFLD)%n,AA%n_pts) )
        ALLOCATE( AA%B(VFLD)%I%dY_f(3,AA%B(VFLD)%n,AA%n_pts_f,AA%FACES) )
        DO I = 1,AA%B(VFLD)%n		! Loop over basis
           ORD1(1:3) = 0
           DO J = 1,AA%n_pts		! Loop over quadrature points
              CALL POLY(AA,VFLD,I,AA%gpt(J,1:3),ORD1,AA%B(VFLD)%I%Y(I,J))	! Calculate basis
           END DO
           DO K = 1,3
             ORD1(1:3) = 0; ORD1(K) = 1
             DO J = 1,AA%n_pts		! Loop over quadrature points
                CALL POLY(AA,VFLD,I,AA%gpt(J,1:3),ORD1,AA%B(VFLD)%I%dY(K,I,J))	! Calculate basis
             END DO
           END DO
        END DO
        DO S = 1,AA%FACES	! Loop over facets
          DO I = 1,AA%B(VFLD)%n		! Loop over basis
            ORD1(1:3) = 0
            DO J = 1,AA%n_pts_f			! Loop over facet quadrature points
              CALL POLY(AA,VFLD,I,AA%gpt_f(J,1:3,S),ORD1,AA%B(VFLD)%I%Y_f(I,J,S))	! Calculate basis
            END DO
            DO K = 1,3
              ORD1(1:3) = 0; ORD1(K) = 1
              DO J = 1,AA%n_pts_f		! Loop over quadrature points
                CALL POLY(AA,VFLD,I,AA%gpt_f(J,1:3,S),ORD1,AA%B(VFLD)%I%dY_f(K,I,J,S))	! Calculate basis
              END DO
            END DO
          END DO
        END DO
      END SUBROUTINE LOAD_INTEGRALS


      SUBROUTINE BUILD_HEXA_XI(AA,VFLD)
!       Subprogram BUILD_HEXA_XI - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine builds the XI coordinate for the tensor basis AA%B(VFLD)
!
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K,S,m,n,G,VFLD
        DOUBLE PRECISION Bx(AA%B(VFLD)%nl)
        Bx(1:AA%B(VFLD)%nl) = AA%B(VFLD)%XI(1:AA%B(VFLD)%nl,1)
        DO I = 1,AA%B(VFLD)%n		! Loop over basis
           AA%B(VFLD)%XI(I,1:3) = Bx(AA%B(VFLD)%B_ID(I,1:3))	! calculate xi point
        END DO
      END SUBROUTINE BUILD_HEXA_XI


      SUBROUTINE BUILD_HEXA_QUADRATURE(AA)
!       Subprogram BUILD_HEXA_QUADRATURE - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine builds the volume/facet quadrature rule for the tensor basis AA
!
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K,S,m,n,G
        DOUBLE PRECISION Ax(AA%n_ptsl),Bx(AA%n_ptsl)
        Ax(1:AA%n_ptsl) = AA%gpt(1:AA%n_ptsl,1)
        Bx(1:AA%n_ptsl) = AA%gw(1:AA%n_ptsl)
        n = 0
        DO I = 1,AA%n_ptsl	! Loop over basis DIM 3
           IF ((AA%DM == 3).OR.(I == 1)) THEN	! If the dimension is 3 or I = 1
           DO J = 1,AA%n_ptsl	! Loop over basis DIM 2
              IF ((AA%DM >= 2).OR.(J == 1)) THEN	! If the dimension is > 2 or J = 1
              DO K = 1,AA%n_ptsl	! Loop over basis DIM 1
                 IF (AA%DM >= 1) THEN		! If the dimension is > 1 or J = 1
                   n = n + 1
                   AA%gpt(n,1:3) = 0; AA%gw(n) = 0	! Initialize volume quadrature rule, and build
                   IF (AA%DM >= 1) AA%gpt(n,1) = Ax(K)
                   IF (AA%DM >= 2) AA%gpt(n,2) = Ax(J)
                   IF (AA%DM >= 3) AA%gpt(n,3) = Ax(I)
                   IF (AA%DM == 1) AA%gw(n) = Bx(K)
                   IF (AA%DM == 2) AA%gw(n) = Bx(K)*Bx(J)
                   IF (AA%DM == 3) AA%gw(n) = Bx(K)*Bx(J)*Bx(I)
                 END IF
              END DO
              END IF
           END DO
           END IF
        END DO
        m = 0
        DO I = 1,AA%n_ptsl	! Loop over basis DIM 2
           IF ((AA%DM-1 == 2).OR.(I == 1)) THEN		! If the dimension is 2 or I = 1
           DO J = 1,AA%n_ptsl		! Loop over basis DIM 1
              IF ((AA%DM-1 >= 1).OR.(J == 1)) THEN	! If the dimension is > 1 or J = 1
                m = m + 1
                AA%gpt_f(m,1:3,1) = 0; AA%gw_f(m) = 0	! Initialize facet quadrature rule, and build
                IF (AA%DM-1 >= 1) AA%gpt_f(m,1,1) = Ax(J)
                IF (AA%DM-1 >= 2) AA%gpt_f(m,2,1) = Ax(I)
                IF (AA%DM-1 == 1) AA%gw_f(m) = Bx(J)
                IF (AA%DM-1 == 2) AA%gw_f(m) = Bx(J)*Bx(I)
              END IF
           END DO
           END IF
        END DO
        AA%VL_f(1:AA%FACES) = AA%VL**REAL(AA%DM-1)	! Initialize facet volume
        AA%VL = AA%VL**REAL(AA%DM)	! Initialize volume
      END SUBROUTINE BUILD_HEXA_QUADRATURE


      SUBROUTINE BUILD_FACE_GPT_ARRAY(AA)
!       Subprogram BUILD_FACE_GPT_ARRAY - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine builds facet quadrature (must be preceeded by a call to BUILD_HEXA_QUADRATURE,
!       if necessary).
!
        USE CMGUI_VARS
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K
        DOUBLE PRECISION Ax(3),Bx,Cx(3),Dx
        AA%gpt_f(1:AA%n_pts_f,AA%DM:3,1) = 0_DPR	! Initialize first quadrature rule
        DO I = 1,AA%FACES
           AA%nrm(1:3,I) = 0_DPR			! Initialize normal vectors
        END DO
        IF (AA%DM == 3) THEN		! If in R^3
         IF (AA%FACES == 4) THEN	! master element is a tetrahedron, build appropriate facet quadrature
          AA%nrm(3,1) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,2) = 0_DPR
          AA%gpt_f(1:AA%n_pts_f,3,2) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(2,2) = -1_DPR


          AA%gpt_f(1:AA%n_pts_f,1,3) = 0_DPR
          AA%gpt_f(1:AA%n_pts_f,2,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,3) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(1,3) = -1_DPR

          DO I = 1,AA%n_pts_f
          AA%gpt_f(I,1,4) = 1_DPR-AA%gpt_f(I,1,1)-AA%gpt_f(I,2,1)
          END DO
          AA%gpt_f(1:AA%n_pts_f,2,4) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,4) = AA%gpt_f(1:AA%n_pts_f,2,1)
          Ax(1) = 1_DPR
          Ax(2) = 3_DPR
          Ax(2) = Ax(2)**0.5
          AA%nrm(1:3,4) = Ax(1)/Ax(2)
         ELSE IF (AA%FACES == 6) THEN	! master element is a cube, build appropriate facet quadrature
          AA%nrm(3,1) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,2) = 0_DPR
          AA%gpt_f(1:AA%n_pts_f,3,2) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(2,2) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,3) = 0_DPR
          AA%gpt_f(1:AA%n_pts_f,2,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,3) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(1,3) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,4) = 1_DPR
          AA%gpt_f(1:AA%n_pts_f,2,4) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,4) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(1,4) = 1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,5) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,5) = 1_DPR
          AA%gpt_f(1:AA%n_pts_f,3,5) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(2,5) = 1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,6) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,6) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%gpt_f(1:AA%n_pts_f,3,6) = 1_DPR
          AA%nrm(3,6) = 1_DPR
         END IF
        ELSE IF (AA%DM == 2) THEN	! If in R^2
         IF (AA%FACES == 3) THEN	! master element is a triangle, build appropriate facet quadrature
          AA%nrm(2,1) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,2) = 0_DPR
          AA%gpt_f(1:AA%n_pts_f,2,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,2) = 0_DPR
          AA%nrm(1,2) = -1_DPR

          DO I = 1,AA%n_pts_f
          AA%gpt_f(I,1,3) = 1_DPR-AA%gpt_f(I,1,1)
          END DO
          AA%gpt_f(1:AA%n_pts_f,2,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,3) = 0_DPR
          Ax(1) = 1_DPR
          Ax(2) = 2_DPR
          Ax(2) = Ax(2)**0.5
          AA%nrm(1:2,3) = Ax(1)/Ax(2); AA%nrm(3,3) = 0_DPR
         ELSE IF (AA%FACES == 4) THEN	! master element is a square, build appropriate facet quadrature
          AA%nrm(2,1) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,2) = 0_DPR
          AA%gpt_f(1:AA%n_pts_f,3,2) = 0_DPR
          AA%nrm(1,2) = -1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,3) = 1_DPR
          AA%gpt_f(1:AA%n_pts_f,3,3) = 0_DPR
          AA%nrm(2,3) = 1_DPR

          AA%gpt_f(1:AA%n_pts_f,1,4) = 1_DPR
          AA%gpt_f(1:AA%n_pts_f,2,4) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,4) = 0_DPR
          AA%nrm(1,4) = 1_DPR
         END IF
        END IF
      END SUBROUTINE BUILD_FACE_GPT_ARRAY






