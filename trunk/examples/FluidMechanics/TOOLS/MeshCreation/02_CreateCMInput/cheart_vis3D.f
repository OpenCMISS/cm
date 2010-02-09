      MODULE VARSKINDS
!       Module VARSKINDS - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This module defines the flags for single precision (SPR)
!       and double precision (DPR)
!
!          i.e.  1_DPR is a double precision variable 15 digits in length
!                1_SPR is a single precision variable 6 digits in length 
!
        INTEGER, PARAMETER :: SPR=SELECTED_REAL_KIND(6,15)
        INTEGER, PARAMETER :: DPR=SELECTED_REAL_KIND(15,307)
      END MODULE VARSKINDS


      MODULE CMGUI_VARS
        USE VARSKINDS
        INCLUDE 'cheart_cmgui3D.h'
        TYPE (ARRAY_MESH), pointer :: MH(:)
        TYPE (ARRAY_MESH) Ml, Vl, Pr, Cl, Clx, Ll, Ml2, Vl2, Pr2, Cl2
        TYPE (ARRAY_PROBLEM_BASE) E, E2
        INTEGER FLD,subFLD(20),DIMEN
        INTEGER INIT_ELEM, INIT_NODE
        DOUBLE PRECISION GLOBAL_TOL,detT,VOL
        DOUBLE PRECISION, pointer :: F(:,:), F_T(:,:), jac(:), GrU(:,:,:), Uxyz(:,:)
        CHARACTER*90 NIMZ,EXTN
        CHARACTER*90 NAMz,OUT_NODE,OUT_ELEM
        CHARACTER*16 OUT_SETNM,FLD_NAMz(20,20)
        CHARACTER*16 BASE_FILE
        CHARACTER*2 NMs(99),KNOT
      END MODULE CMGUI_VARS


       PROGRAM PRO
       USE CMGUI_VARS
       DOUBLE PRECISION Ax, Bx
       INTEGER I,J,K, OPT
       GLOBAL_TOL = 0.00000000001
       EXTN = ''
       CALL BUILD_NMS
       CALL VISUALIZATION
       END


       SUBROUTINE VISUALIZATION
       USE CMGUI_VARS
       DOUBLE PRECISION THRES,Ax,Pt(3),DCON(1000),dL,dC,DC2(3,1000),P1(3),P2(3)
       INTEGER FL,LENG3,IK,ORD(6,4),m,n,S,L_bnd,N1(4),N2(4),L_SC, INTERP,CHOICE
       INTEGER I,J,K,SKIP,Lo,IJ,OPT,lM,Lo1,CNL,ST,A,B,P,L1,L2,TYP, INDX(100)
       DOUBLE PRECISION, pointer :: U(:,:), DATA2(:,:)
       INTEGER, pointer :: C2(:,:),LS(:),BND(:,:)
       DOUBLE PRECISION, pointer :: CN(:,:)
       INTEGER, pointer :: SEGcon(:,:)
       CHARACTER*30 FIL(250),LIN
       CHARACTER*60 CHAR1,CHAR2
       CHARACTER*16 inJOY, outJOY

       INTEGER CNT
       CHARACTER*60 NAMZ2

       I = 0
       DO WHILE (I < 10)
       I = I + 1
       CALL GET_COMMAND_ARGUMENT(I,NAMZ2,J,K)
           IF (INDEX(NAMZ2,'-iF') == 1) THEN		! If the initial map flag
               I = I + 1
               CALL GET_COMMAND_ARGUMENT(I,NAMZ,J,K)
               NAMZ = TRIM(NAMZ)
	  END IF
      END DO

       WRITE(*,*)NAMZ

! ! !         PRINT *, 'Please Edit Cubit file by removing all header information'
! ! !         PRINT *, 'And header info between nodes and elements'
! ! !         PRINT *, 'Then put two numbers preceding the nodes which are the'
! ! !         PRINT *, 'number of nodes, then number of elements'
! ! !         PRINT *, 'PLEASE chose kind of modified CUBIT file'
! ! ! 	WRITE(*,*)
! ! ! 	WRITE(*,*)'(1) linear hex element...'
! ! ! 	WRITE(*,*)'(2) linear tet element...'
! ! ! 	READ(*,*) CHOICE
! ! ! 	IF(CHOICE==1) THEN
! ! ! 	WRITE(*,*)'Reading input file >> hex_mesh.inp << ...'
! ! !         NAMZ='hex_mesh.inp'
! ! ! 	FLD=8
! ! ! 	ELSE IF(CHOICE==2) THEN
! ! ! 	WRITE(*,*)'Reading input file >> tet_mesh.inp << ...'
! ! !         NAMZ='tet_mesh.inp'
! ! ! 	FLD=4
! ! ! 	ELSE
! ! ! 	STOP
! ! ! 	END IF

	IF(NAMZ=='hex_mesh.inp') THEN
	  CHOICE=1
	ELSE IF(NAMZ=='tet_mesh.inp') THEN
	  CHOICE=2
	ELSE IF(NAMZ=='quad_mesh.inp') THEN
	  CHOICE=3
	ELSE IF(NAMZ=='tri_mesh.inp') THEN
	  CHOICE=4
	ELSE
	WRITE(*,*)'Wrong input file format!'
	STOP
	END IF


	IF(CHOICE==1) THEN
	WRITE(*,*)'Reading input file >> hex_mesh.inp << ...'
	FLD=8
	ELSE IF(CHOICE==2) THEN
	WRITE(*,*)'Reading input file >> tet_mesh.inp << ...'
	FLD=4
	ELSE IF(CHOICE==3) THEN
	WRITE(*,*)'Reading input file >> quad_mesh.inp << ...'
	FLD=4
	ELSE IF(CHOICE==4) THEN
	WRITE(*,*)'Reading input file >> tri_mesh.inp << ...'
	FLD=3
	ELSE
	WRITE(*,*)'Wrong choice!'
	STOP
	END IF


        OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
        READ(1,*) Ml%Lx, Ml%Lt

	WRITE(*,*)'Number of nodes: ',Ml%Lx,' Number of elements: ', Ml%Lt


        ALLOCATE(Ml%X(Ml%Lx,3),Ml%T(Ml%Lt,FLD))
        DO I = 1, Ml%Lx
          READ(1,*) J,Ml%X(I,1:3)
        END DO


        DO I = 1,Ml%Lt
           READ(1,*) J,Ml%T(I,1:FLD)
        END DO



        IF (FLD == 8) THEN
           DO I = 1,Ml%Lt
              INDX(1:FLD) = Ml%T(I,1:FLD)
              INDX(3)  = Ml%T(I,4); INDX(4)  = Ml%T(I,3)
              INDX(7)  = Ml%T(I,8); INDX(8)  = Ml%T(I,7)
              Ml%T(I,1:FLD) = INDX(1:FLD)
           END DO
        END IF
       CLOSE(1)
       WRITE(*,*)
       IF(CHOICE==1) THEN
       WRITE(*,*)'Exporting nodes to file >> hex_linear.X << ...'
       WRITE(*,*)
       NAMZ='hex_linear.C'
       ELSE IF(CHOICE==2) THEN
       WRITE(*,*)'Exporting nodes to file >> hex_linear.X << ...'
       WRITE(*,*)
       NAMZ='tet_linear.C'
       ELSE IF(CHOICE==3) THEN
       WRITE(*,*)'Exporting nodes to file >> quad_linear.X << ...'
       WRITE(*,*)
       NAMZ='quad_linear.C'
       ELSE IF(CHOICE==4) THEN
       WRITE(*,*)'Exporting nodes to file >> hex_linear.X << ...'
       WRITE(*,*)
       NAMZ='tri_linear.C'
       ELSE
	  STOP
       ENDIF
       OPEN(UNIT = 1, FILE=NAMz,STATUS='unknown')
        WRITE(1,*) Ml%Lx
        DO I = 1, Ml%Lx
           WRITE(1,*) Ml%X(I,1:3)
        END DO
       CLOSE(1)
       WRITE(*,*)
       IF(CHOICE==1) THEN      
       WRITE(*,*)'Exporting elements to file >> tet_linear.M << ...'
       WRITE(*,*)
       NAMZ='hex_linear.M'
       ELSE IF(CHOICE==2) THEN
       WRITE(*,*)'Exporting elements to file >> tet_linear.M << ...'
       WRITE(*,*)
       NAMZ='tet_linear.M'
       ELSE IF(CHOICE==3) THEN
       WRITE(*,*)'Exporting elements to file >> quad_linear.M << ...'
       WRITE(*,*)
       NAMZ='quad_linear.M'
       ELSE IF(CHOICE==4) THEN
       WRITE(*,*)'Exporting elements to file >> tri_linear.M << ...'
       WRITE(*,*)
       NAMZ='tri_linear.M'
       ELSE
	  WRITE(*,*)'Exporting data failed...'
	  STOP
       ENDIF

       OPEN(UNIT = 1, FILE=NAMz,STATUS='unknown')
        WRITE(1,*) Ml%Lt
        DO I = 1, Ml%Lt
           WRITE(1,*) Ml%T(I,1:FLD)
        END DO
       CLOSE(1)
       WRITE(*,*)
       WRITE(*,*)'Press ENTER to exit'
       READ(*,*)
       WRITE(*,*)
       WRITE(*,*)'...finished successfully.'
       WRITE(*,*)
       END SUBROUTINE VISUALIZATION


      SUBROUTINE READ_MESH_COREI(A,n,m)
        INTEGER I,n,m,A(n,m)
        DOUBLE PRECISION Bx(m)
        DO I = 1,n
          READ(1,*,END=30) Bx(1:m); A(I,1:m) = INT(Bx(1:m))
        END DO
        RETURN
 30     PRINT *, 'END OF FILE REACHED!'
        STOP
      END SUBROUTINE READ_MESH_COREI


      SUBROUTINE READ_MESH_CORED(A,n,m)
        INTEGER I,n,m
        DOUBLE PRECISION A(n,m)
        DO I = 1,n
          READ(1,*,END=35) A(I,1:m)
        END DO
        RETURN
 35     PRINT *, 'END OF FILE REACHED!'
        STOP
      END SUBROUTINE READ_MESH_CORED


      SUBROUTINE STREAM_POINTS
        USE CMGUI_VARS
        TYPE (ARRAY_MESH) s_Ml(2), s_Vl(2)
        INTEGER I,J,K,n,FL,n_STEPS,NOP(10000),n_OUT,my_INTERVAL,my_OUT
        INTEGER SKIP,Lo,IJ,OPT,lM,Lo1,CNL,ST,A,B,P,L1,L2,TYP
        DOUBLE PRECISION TSTEPS(10000),DT,my_DT,Xi(3),P1(3),Ax,Vp(3),TEMP(10),alp
        CHARACTER * 16 F_NAMz(1000), COMFILE
        FL = 6
        FLD = 2
        DIMEN = 3
        

      END SUBROUTINE STREAM_POINTS


      SUBROUTINE WRITE_NODES(AA)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) AA
        INTEGER I,J,K,n,Fl
        PRINT *, 'DESIRED NODE FILE NAME (including extension):'
        READ *, OUT_NODE
        PRINT *, 'CMGUI data set name:'
        READ *, OUT_SETNM
        PRINT *, 'DESIRED STARTING NODE:'
        READ *, INIT_NODE
        INIT_NODE = INIT_NODE - 1

       OPEN(UNIT=2,FILE=OUT_NODE,STATUS='unknown')
       WRITE(2,*) 'Group name: '//TRIM(OUT_SETNM)
       WRITE(2,*) '#Fields='//TRIM(NMs(FLD))
       n = 0
       DO Fl = 1,FLD
         WRITE(2,*) TRIM(NMs(Fl))//') '//TRIM(FLD_NAMz(Fl,1))// &
           ', coordinate, rectangular cartesian, #Components='//TRIM(NMs(subFLD(Fl)))
         DO I = 1,subFLD(Fl)
            n = n + 1
            WRITE(2,*) ' '//TRIM(FLD_NAMz(Fl,1+I))//'.  Value index='//TRIM(NMs(n))//', #Derivatives=0'
         END DO
       END DO
       DO I = 1,AA%Lx
         WRITE(2,*) 'Node:',I+INIT_NODE
         n = 0
         DO J = 1,FLD
            DO K = 1,subFLD(J)
               n = n + 1
               WRITE(2,*) AA%X(I,n)
            END DO
         END DO
       END DO
       WRITE(2,*) ' '
       CLOSE(2)
      END SUBROUTINE WRITE_NODES


      SUBROUTINE WRITE_ELEMENTS(AA,BB,ID)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) AA
        TYPE(ARRAY_BASE) BB
        INTEGER I,J,K,m,n,Fl,A(3),INITN(BB%n),ID,L_INITN
        REAL SF(BB%n)
        CHARACTER*60 ELEM_TYPE,TYPE2,PWR(3)
        INITN(1:BB%n) = INIT_NODE
        SF(1:BB%n) = 1_DPR
        PWR(1) = 'l'; PWR(2) = 'q'; PWR(3) = 'c'
        PRINT *, 'DESIRED ELEM FILE NAME (including extension):'
        READ *, OUT_ELEM
        PRINT *, 'DESIRED STARTING ELEMENT'
        READ *, INIT_ELEM


        IF (E%TRI_BASIS == 1) THEN
          ELEM_TYPE = 'l.simplex(2)*l.simplex'
          TYPE2 = 'simplex(2)*simplex'
          A(2) = E%FNODES
          A(3) = BB%n
        ELSE IF (E%TET_BASIS == 1) THEN
          IF (BB%n == 4) THEN
            ELEM_TYPE = 'l.simplex(2;3)*l.simplex*l.simplex'
            TYPE2 = 'simplex(2;3)*simplex*simplex'
          ELSE IF ( BB%n == 10 ) THEN
            ELEM_TYPE = 'q.simplex(2;3)*q.simplex*q.simplex'
            TYPE2 = 'simplex(2;3)*simplex*simplex'
          ELSE
            PRINT *, ' >> ELEMENT NOT SUPPORTED IN CMGUI << '; STOP
          END IF



        ELSE IF (E%QUAD_BASIS == 1) THEN
          n = INT(REAL(BB%n)**(0.5))-1; PRINT *, n
          ELEM_TYPE = TRIM(PWR(n))//'.Lagrange*'//TRIM(PWR(n))//'.Lagrange'
        ELSE IF (E%HEX_BASIS == 1) THEN
          n = INT(REAL(BB%n)**(1.0/3.0))-1; PRINT *, n
          ELEM_TYPE = TRIM(PWR(n))//'.Lagrange*'//TRIM(PWR(n))//'.Lagrange*'//TRIM(PWR(n))//'.Lagrange'
        ELSE 
          CALL ERRORZ(1)
        END IF


        PRINT *, ELEM_TYPE
        A(2) = E%FNODES
        A(3) = BB%n
       OPEN(UNIT=2,FILE=OUT_ELEM,STATUS='unknown')
       WRITE(2,*) 'Group name: '//TRIM(OUT_SETNM)
       IF ((E%TRI_BASIS == 1).OR.(E%TET_BASIS == 1)) THEN
         WRITE(2,*) 'Shape.  Dimension='//TRIM(NMs(DIMEN))//', '//TRIM(TYPE2)
         IF ((BB%n == 4).OR.(BB%n == 3)) THEN
           WRITE(2,*) '#Scale factor sets= 1'
           WRITE(2,*) ' '//TRIM(ELEM_TYPE)//', #Scale factors='//TRIM(NMs( A(DIMEN) ))
         ELSE
           WRITE(2,*) '#Scale factor sets= 0'
         END IF
       ELSE
         WRITE(2,*) 'Shape.  Dimension='//TRIM(NMs(DIMEN))
         WRITE(2,*) '#Scale factor sets= 1'
         WRITE(2,*) ' '//TRIM(ELEM_TYPE)//', #Scale factors='//TRIM(NMs( A(DIMEN) ))
       END IF
       WRITE(2,*) '#Nodes= '//TRIM(NMs( A(DIMEN) ))
       WRITE(2,*) '#Fields='//TRIM(NMs(FLD))
       DO Fl = 1,FLD
         WRITE(2,*) TRIM(NMs(Fl))//') '//TRIM(FLD_NAMz(Fl,1))// &
           ', coordinate, rectangular cartesian, #Components='//TRIM(NMs(subFLD(Fl)))
         DO I = 1,subFLD(Fl)
            WRITE(2,*) ' '//TRIM(FLD_NAMz(Fl,1+I))//'.  '//TRIM(ELEM_TYPE)//', no modify, standard node based.'
            WRITE(2,*) '   #Nodes= '//TRIM(NMs( A(DIMEN) ))
            DO J = 1,A(DIMEN)
               WRITE(2,*) '    '//TRIM(NMs(J))//'.  #Values=1'
               WRITE(2,*) '     Value indices:     1'
               IF ((E%FACES == 4).AND.(BB%n == 10)) THEN
                 WRITE(2,*) '     Scale factor indices:    0'
               ELSE
                 WRITE(2,*) '     Scale factor indices:   '//TRIM(NMs(J))
               END IF
            END DO
         END DO
       END DO
       n = INIT_ELEM-1
        CALL BUILD_ELE_CORE(AA,BB,INITN,L_INITN,ID)
        PRINT *, INITN(1:L_INITN)
        DO I = 1,AA%Lt
           n = n + 1
           WRITE(2,*) 'Element:', n,' 0  0'
           WRITE(2,*) '   Nodes:'
           WRITE(2,*) '   ', AA%T(I,INITN(1:L_INITN)) + INIT_NODE
           IF ((E%FACES == 4).AND.(BB%n == 10)) THEN
           ELSE
             WRITE(2,*) '   Scale factors:'
             WRITE(2,*) '   ',SF(1:L_INITN)
           END IF
        END DO
       WRITE(2,*) ' '
       CLOSE(2)
      END SUBROUTINE WRITE_ELEMENTS


      SUBROUTINE BUILD_ELE_CORE(AA,BB,INITN,L_INITN,ID)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) AA
        TYPE(ARRAY_BASE) BB
        INTEGER I,J,K,m,n,Fl,A(3),X1,X2,X3,INITN(BB%n),ID,L_INITN, INDX
        DOUBLE PRECISION Ax, Bx, Cx
        IF (ID == 2) THEN
          L_INITN = E%FNODES
          DO J = 1,L_INITN
             INITN(J) = AA%T(I,J)+INIT_NODE
          END DO
        ELSE
          IF (E%HEXA == 1) THEN
           L_INITN = BB%n; INDX = BB%nl
           IF (E%DM == 2) INDX = 1
           K = 0
           DO X3 = 1,INDX
              DO X2 = 1,BB%nl
                 DO X1 = 1,BB%nl
                    DO J = 1,BB%n
                       IF (BB%B_ID(J,1) == X1) THEN
                          IF (BB%B_ID(J,2) == X2) THEN
                             IF (BB%B_ID(J,3) == X3) THEN
                                K = K + 1
                                INITN(K) = J
                             END IF
                          END IF
                       END IF
                    END DO
                 END DO
              END DO
           END DO
          ELSE
             L_INITN = BB%n
             IF ((BB%n == 3) .OR. (BB%n == 4)) THEN
                DO I = 1,BB%n
                  INITN(I) = I
                END DO
             ELSE IF (BB%n == 10) THEN
                INITN(1) = 1; INITN(2) = 5; INITN(3) = 2;
                INITN(4) = 6; INITN(5) = 8; INITN(6) = 3;
                INITN(7) = 7; INITN(8) = 9; INITN(9) = 10;
                INITN(10) = 4;
             END IF
          END IF
        END IF
      END SUBROUTINE BUILD_ELE_CORE


      SUBROUTINE ELEMENT_TYPE(BB)
        USE CMGUI_VARS
        INTEGER I,J,K
        TYPE (ARRAY_PROBLEM_BASE) BB
        PRINT *, ' >> PLEASE INPUT :: BASE SET <<'
        READ *, NIMZ
        NIMZ = TRIM(EXTN)//TRIM(NIMZ)
        CALL LOAD_BASE_SET(BB)
      END SUBROUTINE ELEMENT_TYPE


      SUBROUTINE BUILD_COORDINATES
        USE CMGUI_VARS
        INTEGER I,J,K
        subFLD(1) = 3
        FLD_NAMz(1,1) = 'coordinates'
        DO I = 1,3
           FLD_NAMz(1,1+I) = 'x'//TRIM(NMs(I))
        END DO
      END SUBROUTINE BUILD_COORDINATES


      SUBROUTINE BUILD_NMS
       USE CMGUI_VARS
       INTEGER I,J,K
       KNOT = '0'
       NMs(1) = '1'
       NMs(2) = '2'
       NMs(3) = '3'
       NMs(4) = '4'
       NMs(5) = '5'
       NMs(6) = '6'
       NMs(7) = '7'
       NMs(8) = '8'
       NMs(9) = '9'
       K = 9
       DO I = 1,9
          K = K + 1
          NMs(K) = TRIM(NMs(I))//TRIM(KNOT)
          DO J = 1,9
             K = K + 1
             NMs(K) = TRIM(NMs(I))//TRIM(NMs(J))
          END DO
       END DO
     END SUBROUTINE BUILD_NMS


      SUBROUTINE BUILD_FACE_BNDRY
        USE CMGUI_VARS
        INTEGER I,J,K,S,ELEM,FACE,m,n,G
        DOUBLE PRECISION Ax
        DO I = 1,Ml%Lb
           ELEM = Ml%B(I,1)
           DO J = 1,E%FACES
              Ax = 0
              DO K = 1,E%B(E%SPL)%n
                 S = Ml%T(ELEM,K)
                 n = 0
                 DO m = 2,E%FNODES+1
                    IF (Ml%B(I,m) == S) THEN
                        n = 1
                        EXIT
                    END IF
                 END DO
                 IF (n == 1) THEN
                    DO G = 1,E%n_pts_f
                       Ax = Ax + E%gw_f(G)*E%B(E%SPL)%I%Y_f(K,G,J)
                    END DO
                 END IF
              END DO
              IF (abs(1_DPR-Ax).LT.0.00000001) THEN
                  FACE = J
                  EXIT
              END IF
           END DO
           Ml%B(I,E%FNODES+3) = FACE
        END DO
      END SUBROUTINE BUILD_FACE_BNDRY


      SUBROUTINE ADD_EXTENSION(NAMZ2)
        USE CMGUI_VARS
        CHARACTER*90 NAMZ2
        IF (TRIM(EXTN).NE.'') THEN
            NAMZ2 = TRIM(EXTN)//TRIM(NAMZ2)
        END IF
      END SUBROUTINE ADD_EXTENSION



      SUBROUTINE MAKE_BOUNDARY_FILE
       USE CMGUI_VARS
       DOUBLE PRECISION THRES,Ax,Pt(3),DCON(1000),dL,dC,DC2(3,1000),P1(3),P2(3)
       INTEGER FL,LENG3,IK,ORD(6,4),m,n,S,L_bnd,N1(4),N2(4),OPT, L_set, ID
       INTEGER I,J,K,SKIP,Lo,IJ,lM,Lo1,CNL,ST,A,B,P,L1,L2,TYP, INDX(100)
       INTEGER LNUM(2000000), L_lnum,DUMM(100)
       DOUBLE PRECISION, pointer :: U(:,:), DATA2(:,:)
       INTEGER, pointer :: C2(:,:),LS(:),BND(:,:)
       CHARACTER*16 iJY, oJY, FIL(5), NAMl(50)

       PRINT *, 'NOTE, this works to identify boundaries by'
       PRINT *, 'neighboring tet.'
       CALL ELEMENT_TYPE(E)

       PRINT *, 'Input 0th facet file name:'
       READ *, NAMz

       OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
       READ(1,*) Ml%Lx
       ALLOCATE(Ml%X(Ml%Lx,5))
       CALL READ_MESH_CORED(Ml%X(1:Ml%Lx,1:3),Ml%Lx,3)
       CLOSE(1)

       PRINT *, 'Input TETRA/HEXA-HEDERALIZATION file name:'
       READ *, NAMz

       OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
       READ(1,*) Ml%Lt
       ALLOCATE(Ml%T(Ml%Lt,E%B(E%SPL)%n+1))
       CALL READ_MESH_COREI(Ml%T(1:Ml%Lt,1:E%B(E%SPL)%n),Ml%Lt,E%B(E%SPL)%n)
       CLOSE(1)

       PRINT *, 'Input TETRA/HEXA-HEDERALIZATION NEIGHBORING file name:'
       READ *, NAMz

       OPEN(UNIT = 1, FILE=NAMz,STATUS='old')
       READ(1,*) I
       ALLOCATE(Ml%NBT(Ml%Lt,E%FACES))
       CALL READ_MESH_COREI(Ml%NBT(1:Ml%Lt,1:E%FACES),Ml%Lt,E%FACES)
       CLOSE(1)

        Ml%Lb = 0
        DO I = 1,Ml%Lt
           DO J = 1,E%FACES
              IF (Ml%NBT(I,J) == 0) Ml%Lb = Ml%Lb + 1
           END DO
        END DO
        ALLOCATE(Ml%B(Ml%Lb,2+E%FNODES)); Ml%B(1:Ml%Lb,2+E%FNODES) = 0; ID = 0
        DO I = 1,Ml%Lt
           DO J = 1,E%FACES
              IF (Ml%NBT(I,J) == 0) THEN
                 m = 1; ID = ID + 1; Ml%B(ID,1) = I
                 DO K = 1,E%B(E%SPL)%n
                    Ax = 0
                    DO n = 1,E%n_pts_f
                      Ax = Ax + E%B(E%SPL)%I%Y_f(K,n,J)
                    END DO
                    IF (abs(Ax) > GLOBAL_TOL) THEN
                       m = m + 1; Ml%B(ID,m) = Ml%T(I,K)
                    END IF
                 END DO
              END IF
           END DO
        END DO
        Ml%Lt = Ml%Lb; DEALLOCATE(Ml%T); ALLOCATE(Ml%T(Ml%Lt,1+E%FNODES))
        DO I = 1,Ml%Lt
          Ml%T(I,1:1+E%FNODES) = Ml%B(I,2:2+E%FNODES)
        END DO
        PRINT *, ' ==>> Boundary elements :: ', Ml%Lb
        DIMEN = 2; FLD = 2; CALL BUILD_COORDINATES
        subFLD(2) = 1; FLD_NAMz(2,1) = 'boundary_num'; FLD_NAMz(2,2) = 'b_no'
        CALL VISUALIZE_BOUNDARY(Ml,E,4)

        PRINT *,' ** You may now look at your boundary in CMGUI'
        PRINT *,' ** Within CMGUI, create new egroups or ngroups labelling boundaries'
        PRINT *,' ** Remember to export either the nodes or elements'

        PRINT *, ' 1] READING IN AN ELEMENT FILE'
        PRINT *, ' 2] READING IN A NODE FILE'
        READ *, OPT
        PRINT *, 'Input CMGUI file name:'
        READ *, NAMZ

        L_set = 0; I = 0; INDX(1) = 1
        OPEN(UNIT=1,FILE=NAMZ,status='old',action='read')
         DO WHILE (I >= 0)
         I = I + 1; FIL(1:3) = ''
         READ(1,*,ERR=10,END=20) FIL(1)
 10      IF (INDEX(FIL(1),'Group') >= 1) THEN
            L_set = L_set + 1; INDX(L_set+1) = I
         END IF
         END DO
 20     CLOSE(1); INDX(L_set+2) = I-1
        OPEN(UNIT=1,FILE=NAMZ,status='old',action='read')
         m = 0
         DO I = 1,L_set
           DO J = 1,INDX(I+1)-INDX(I) - m
             READ(1,*)
           END DO
           m = 1; READ(1,*) FIL(1), FIL(2), NAMl(I)
         END DO
        CLOSE(1)
        PRINT *, ' ==>> Boundary sets :: ', L_set
        DO I = 1,L_set
           PRINT *,  NAMl(I), INDX(I+1)
        END DO

        PRINT *, ' ** You will now be prompted for different set names/boundary IDs'
        iJY = 'y'
        DO WHILE (iJY == 'y')
          PRINT *, 'PLEASE INPUT set name'
          READ *, oJY
          PRINT *, 'What is the associated boundary number?'
          READ *, ID
          Ml%X(1:Ml%Lx,4) = 0
          DO I = 1,L_set
            IF (INDEX(NAMl(I),oJY) >= 1) K = I+1
          END DO
          PRINT *, ' ==>> READING SET '//TRIM(NAMl(I-1)) 

          OPEN(UNIT=1,FILE=NAMZ,status='old',action='read')
           DO I = 1,INDX(K)
              READ(1,*)
           END DO
           I = INDX(K); L_lnum = 0; LNUM(1) = 1
           DO WHILE (I < INDX(K+1))
             I = I + 1; READ(1,*,ERR=30) FIL(1)
 30          IF (OPT == 1) THEN
               IF (INDEX(FIL(1),'Nodes:') >= 1) THEN
                 L_lnum = L_lnum + 1; LNUM(L_lnum+1) = I - INDX(K)
               END IF
             ELSE IF ( OPT == 2) THEN
               IF (INDEX(FIL(1),'Node:') >= 1) THEN
                 L_lnum = L_lnum + 1; LNUM(L_lnum+1) = I - INDX(K)
               END IF
             END IF
           END DO
          CLOSE(1)
          OPEN(UNIT=1,FILE=NAMZ,status='old',action='read')
           DO I = 1,INDX(K)
              READ(1,*)
           END DO
           m = 0
           DO I = 1,L_lnum
             DO J = 1,LNUM(I+1)-LNUM(I)-m
               READ(1,*)
             END DO
             IF (OPT == 1) THEN
                m = 2; READ(1,*)
                READ(1,*) DUMM(1:E%FNODES)
                Ml%X(DUMM(1:E%FNODES),4) = 1
             ELSE IF (OPT == 2) THEN
                m = 1; READ(1,*) FIL(1), K
                IF (K > 0) Ml%X(K,4) = 1
             END IF
           END DO
          CLOSE(1)
          n = 0
          DO I = 1,Ml%Lb
            K = 0
            DO J = 2,E%FNODES+1
              K = K + INT(Ml%X( Ml%B(I,J) ,4))
            END DO
            IF (K == E%FNODES) THEN
              Ml%B(I,2+E%FNODES) = ID; n = n + 1
            END IF
          END DO
          PRINT *,  n,' boundary elements added to group', ID
          PRINT *, 'INPUT another boundary set (y/n)?'
          READ *, iJY
        END DO

        PRINT *, 'PLEASE INPUT BOUNDARY file name'
        READ *, NAMZ
        OPEN(UNIT = 1, FILE=NAMz,STATUS='unknown')
         WRITE(1,*) Ml%Lb
         DO I = 1, Ml%Lb
            WRITE(1,*) Ml%B(I,1:E%FNODES+2)
         END DO
        CLOSE(1)
      END SUBROUTINE MAKE_BOUNDARY_FILE


      SUBROUTINE VISUALIZE_BOUNDARY(AA,EE,FL)
        USE CMGUI_VARS
        TYPE(ARRAY_MESH) AA
        TYPE(ARRAY_PROBLEM_BASE) EE
        INTEGER I,J,K,m,n,FL,INDX(FL)
        IF (FL < 4) STOP
        DO I = 1, AA%Lt
           DO J = 1,EE%FNODES
             AA%X(AA%T(I,J),4:FL) = 0
           END DO
        END DO
        DO I = 1,AA%Lx
           IF (AA%X(I,4) > GLOBAL_TOL) THEN
             AA%X(I,1:3) = -100000; AA%X(I,4:FL) = 0
           END IF
        END DO
        INDX(1:FL) = 0; m = 0
        DO I = 1, AA%Lt
           K = 0
           DO J = 1,m
             IF (INDX(J) == AA%T(I,EE%FNODES+1)) K = J + 4
           END DO
           IF (K == 0) THEN
             m = m + 1; INDX(m) = AA%T(I,EE%FNODES+1); K = m + 4
           END IF
        END DO
        CALL ORDER_SORT_I(m,INDX(1:m))
        DO I = 1, AA%Lt
           K = 0
           DO J = 1,m
             IF (INDX(J) == AA%T(I,EE%FNODES+1)) THEN
                K = J + 4; EXIT
             END IF
           END DO
           IF ( K == 0 ) THEN
              PRINT *, 'ERROR IN IDENTIFYING BOUNDARY'; STOP
           END IF
           DO J = 1,EE%FNODES
             AA%X(AA%T(I,J),4) = AA%T(I,EE%FNODES+1)
             AA%X(AA%T(I,J),K) = 1
           END DO
        END DO
        CALL WRITE_NODES(AA)
!        CALL WRITE_ELEMENTS(AA,EE%B(EE%SPL),2)
      END SUBROUTINE VISUALIZE_BOUNDARY


      SUBROUTINE ORDER_SORT_I(n,ID)
!       ORDER_SORT - Written by David A. Nordsletten (C) 2008
!       This subroutine sorts a list of integers.
!
!       VARIABLES ::  I, J, K, n (INTEGER), dummy variables
!
        USE CMGUI_VARS
        INTEGER n, I, J, ID(n), ID2(n), m, K
        ID2 = ID
        DO I = n,1,-1
           m = -100000
           DO J = 1,n
              IF (m < ID2(J)) THEN
                 m = ID2(J); K = J
              END IF
           END DO
           ID2(K) = -100000
           ID(I) = m
        END DO
      END SUBROUTINE ORDER_SORT_I
