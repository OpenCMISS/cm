!    CMHeart - Preparatory routines
!   =======================================
!     For use in generating general function spaces for nodal lagrange interpolations
!
!     Written by David A. Nordsletten (C) 2008
!     DPhil Student, Comlab, University of Oxford 2006-2008
!     PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!

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


      MODULE VARS_PREP
!       Module VARS_PREP - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This module defines the common variables for the preparatory routines
!
        USE VARSKINDS
        INCLUDE 'cheart_prepvars.h'
        TYPE (ARRAY_MESH) Fl					! Input Spatial Tessellation
        TYPE (ARRAY_MESH), pointer :: Ml(:)			! Output Function Spaces
        TYPE (ARRAY_PROBLEM_BASE) E, E2				! Initial and final Basis sets
        DOUBLE PRECISION GLOBAL_TOL				! Global tolerance parameter
        INTEGER, pointer :: nL(:,:)				! Work array, (T)^-1 map
        INTEGER LC						! Common index for nL
        CHARACTER*60 MAPFILE, COFILE, iBASE, fBASE, FROOT, BNDRYFILE	! Input file names
        CHARACTER*90 NIMZ,EXTN					! Character variales
      END MODULE VARS_PREP


      PROGRAM CMHeart_preparatory
!       Program CMHeart_preparatory - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program executes the preparatory routines for CMHeart
!
       USE VARS_PREP
       IMPLICIT NONE
       INTEGER I,J,K
       GLOBAL_TOL = 0.00000001
       EXTN = ''	! Set extension to basis root file
       CALL COMM_LINE_IO	! Obtain Command line inputs
       CALL INPUT_CORE		! Read command line input files and build basis
       CALL MAKE_NEIGH		! Construct neighboring array for spatial tessellation
       DO I = 1,E2%n_B
         Ml(I)%Lb = 0
         IF (INDEX(E2%B(I)%CL,'M') == 1) K = I
         CALL MAKE_FUNCTION_SPACE( Ml(I), E2%B(I) )	! Construct function space
         CALL PRINTZ(3,I)	! Output completion of function space
       END DO
       IF (Fl%Lb > 0) THEN
          Ml(K)%ID = K; CALL RECALCULATE_BOUNDARY( Ml( K ), E2, Fl, E )
       END IF
       CALL OUTPUT_FILE		! Output updated function spaces
       CALL PRINTZ(100,0)	! Output completion of routines
      END


      SUBROUTINE RECALCULATE_BOUNDARY( NW, Enw, OD, Eod )
!       Subprogram RECALCULATE_BOUNDARY - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
        USE VARS_PREP
        TYPE (ARRAY_MESH) OD, NW
        TYPE (ARRAY_PROBLEM_BASE) Enw, Eod
        INTEGER I,J,K,m,n,n_l,P(3),S,n_o,LS(1000)
        DOUBLE PRECISION XInw(3), XIod(3), Ax, Bx
        NW%Lb = OD%Lb; ALLOCATE( NW%B(NW%Lb, 2 + Enw%FNODES) )
        NW%B(1:NW%Lb,1) = OD%B(1:OD%Lb,1); NW%B(1:NW%Lb,2 + Enw%FNODES) = OD%B(1:OD%Lb,2 + Eod%FNODES)
        DO I = 1,OD%Lb
          Ax = 0; XIod(1:3) = 0
          DO J = 1,Eod%B(OD%ID)%n
            DO K = 2, 1 + Eod%FNODES
             IF (OD%T(OD%B(I,1),J) == OD%B(I,K)) THEN
               XIod(1:3) = XIod(1:3) + Eod%B(OD%ID)%XI(J,1:3); Ax = Ax + 1; EXIT
             END IF
            END DO
          END DO
          XIod(1:3) = XIod(1:3) / Ax
          PRINT *, XIod(1:3)
          DO J = 1,Enw%FACES
            IF (OD%NBT(OD%B(I,1),J) == 0) THEN
              Ax = 0; XInw(1:3) = 0; n = 1
              DO K = 1,Enw%B(NW%ID)%n
                Bx = 0
                DO m = 1,Enw%n_pts_f
                   Bx = Bx + Enw%B(NW%ID)%I%Y_f(K,m,J)	! Calculate approximate integral on facet
                END DO
                IF (abs(Bx) > GLOBAL_TOL) THEN
                   n = n + 1; NW%B(I,n) = NW%T(NW%B(I,1),K)
                   XInw(1:3) = XInw(1:3) + Enw%B(NW%ID)%XI(K,1:3); Ax = Ax + 1;
                END IF
              END DO
              XInw(1:3) = XInw(1:3)/Ax - XIod(1:3)
              PRINT *, XInw(1:3)
              CALL P_NORM(XInw,3,2,Bx)
              IF (Bx < GLOBAL_TOL) EXIT
            END IF
          END DO
        END DO
      END SUBROUTINE RECALCULATE_BOUNDARY


      SUBROUTINE MAKE_NEIGH
!       Subprogram MAKE_NEIGH - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program builds a map of neighboring elements for each element of the
!       spatial tessellation.  This is accomplished by first constructing the inverse
!       map of Fl%T, identifying all elements for each node (nL).  The set of facet 
!       nodes is constructed for each facet (RW).  The neighbor of the facet J, is 
!       a unique element at the intersection of all I, nL(RW(I,J),:).
!
        USE VARS_PREP
        INTEGER I,J,K,m,n,n_l,P(3),S,n_o,LS(1000)
        INTEGER IRW(E%B(Fl%ID)%n,E%FACES)
        INTEGER RW(E%B(Fl%ID)%n,E%FACES),L_RW(E%FACES)
        LC = 100
        ALLOCATE(Fl%NBT(Fl%Lt,E%FACES),nL(Fl%Lx,LC))
        Fl%NBT(1:Fl%Lt,1:E%FACES) = 0		! Initialize neighboring array as zeros
        nL(1:Fl%Lx,LC) = 0			! nL(:,LC), no. of elements, initialized as zero
        DO I = 1,Fl%Lt
           DO J = 1,E%B(Fl%ID)%n
              nL(Fl%T(I,J),LC) = nL(Fl%T(I,J),LC) + 1	! increment no. elements for basis Fl%T(I,J)
              nL(Fl%T(I,J),nL(Fl%T(I,J),LC)) = I	! add element I to list for basis Fl%T(I,J)
           END DO
        END DO
        CALL BUILD_IRW(IRW,L_RW)	! Generate IRW
        DO I = 1,Fl%Lt
          CALL BUILD_RW(I,RW,IRW,L_RW)	! Generate RW and L_RW
          DO J = 1,E%FACES		! Loop over facets
             n_l = nL(RW(1,J),LC)	! set initial list length
             LS(1:n_l) = nL(RW(1,J),1:n_l)	! set initial list
             DO K = 1,n_l
               IF (LS(K) == I) LS(K) = 0	! delete current element from list
             END DO
             DO K = 2,L_RW(J)	! Loop over all basis on facet J
               DO m = 1,n_l	! Loop over all components from initial element list
                 IF (LS(m) > 0) THEN
                   S = 0
                   DO n = 1,nL(RW(K,J),LC)	! Loop over all elements associated with basis RW(K,J)
                     IF (nL(RW(K,J),n) == LS(m)) THEN ! If LS(m) is shared, keep it
                       S = 1; EXIT
                     END IF
                   END DO
                   IF (S == 0) LS(m) = 0	! LS(m) is not shared, delete it
                 END IF
               END DO
             END DO
             n = 0
             DO K = 1,n_l	! Loop over initial list
               IF (LS(K) > 0) THEN
                 n = n + 1; Fl%NBT(I,J) = LS(K)		! If component is at the intersection, set as neighbor
               END IF
             END DO
             IF (n > 1) CALL ERRORZ(6)		! If the number of interesection components is non-unique, quit
          END DO
        END DO
      END SUBROUTINE MAKE_NEIGH


      SUBROUTINE BUILD_IRW(IRW,L_RW)
!       Subprogram BUILD_IRW - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program builds arrays IRW and L_RW for spatial basis.  IRW (I, J) stores the
!       the Ith basis index on the Jth facet.  L_RW(J) is the number of basis on the Jth facet.
!       These arrays are computed by identifying facets based on those with possible values
!       on a facet.  In nodal lagrange interpolants, non-zeros indicate the basis is on
!       a given facet.
!
        USE VARS_PREP
        INTEGER I,J,K,ELEM,IRW(E%B(Fl%ID)%n,E%FACES),L_RW(E%FACES),FACE
        DOUBLE PRECISION Ax
        DO FACE = 1,E%FACES	! Loop over facets
           IRW(1:E%B(Fl%ID)%n,FACE) = 0		! Initialize RW as zeros
           L_RW(FACE) = 0			! Initialize L_RW as zeros
           DO I = 1,E%B(Fl%ID)%n		! Loop over all basis components of ELEM
              Ax = 0
              DO K = 1,E%n_pts_f
                 Ax = Ax + E%B(Fl%ID)%I%Y_f(I,K,FACE)	! Calculate approximate integral on facet
              END DO
              IF (abs(Ax) > GLOBAL_TOL) THEN
                  L_RW(FACE) = L_RW(FACE) + 1
                  IRW(L_RW(FACE),FACE) = I
              END IF
           END DO
        END DO
      END SUBROUTINE BUILD_IRW


      SUBROUTINE BUILD_RW(ELEM,RW,IRW,L_RW)
!       Subprogram BUILD_RW - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program builds arrays RW and L_RW for element ELEM.  RW (I, J) stores the
!       the Ith basis on the Jth facet.  L_RW(J) is the number of basis on the Jth facet.
!       These arrays are computed by identifying facets based on those with possible values
!       on a facet.  In nodal lagrange interpolants, non-zeros indicate the basis is on
!       a given facet.
!
        USE VARS_PREP
        INTEGER I,J,K,ELEM,RW(E%B(Fl%ID)%n,E%FACES),L_RW(E%FACES),FACE
        INTEGER IRW(E%B(Fl%ID)%n,E%FACES),ne
        DOUBLE PRECISION Ax
        DO FACE = 1,E%FACES	! Loop over facets
           DO I = 1,L_RW(FACE)	! Loop over number of basis for facet FACE
              RW(I,FACE) = Fl%T(ELEM,IRW(I,FACE))	! Set value of RW
           END DO
           K = 1; ne = 100000		! Reorder so that smallest element set is first
           DO I = 1,L_RW(FACE)
              IF (nL(RW(I,FACE),LC) < ne) THEN
                 ne = nL(RW(I,FACE),LC); K = I
              END IF
           END DO
           I = RW(1,FACE); J = RW(K,FACE)
           RW(1,FACE) = J; RW(K,FACE) = I
        END DO
      END SUBROUTINE BUILD_RW


      SUBROUTINE MAKE_FUNCTION_SPACE(AA,EE)
!       Subprogram MAKE_FUNCTION_SPACE - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program generates the function space (AA) on the spatial map (Fl) based
!       on the basis (EE).  This is process is straightforward for discontinuous fields.
!       In the case of continuous fields, facet basis will be shared.  Thus, as the 
!       space is computed for an element ELEM, a list of all previously defined elements 
!       ( I < ELEM) which share a corner, edge, or facet with ELEM is composed.  It is then
!       verified that a basis is unique before added to the coefficient list.
!
        USE VARS_PREP
        TYPE(ARRAY_MESH) AA
        TYPE(ARRAY_BASE) EE
        INTEGER I,J,K,m,n,P,S,ELEM,ELE_ARRAY(1000),L_EA
        DOUBLE PRECISION P1(3),P2(3),dL,WGT(8),TOL,ASH1,ASH2
        ALLOCATE(AA%T(Fl%Lt,EE%n),AA%X(EE%n*Fl%Lt,3))
        AA%Lt = Fl%Lt; AA%Lx = 0		! Initialize map size and coefficient size
        DO ELEM = 1,Fl%Lt		! Loop over elements of Fl
           IF (EE%DISCONT == 0) THEN	! If field is continuous
           L_EA = 0		! Set list length of previous elements to zero
           DO I = 1, E%B(Fl%ID)%n	! Loop over basis of ELEM
              DO J = 1,nL(Fl%T(ELEM,I),LC)	! Loop over elements sharing basis Fl%T(ELEM,I)
                 IF (nL(Fl%T(ELEM,I),J) < ELEM) THEN 	! If the element is previously encountered
                   m = 0
                   DO K = 1,L_EA
                    IF (nL(Fl%T(ELEM,I),J) == ELE_ARRAY(K)) THEN	! See if it exists on list ELE_ARRAY
                      m = 1; EXIT
                    END IF
                   END DO
                   IF (m == 0) THEN		! If not, add to list and update counter L_EA
                     L_EA = L_EA + 1;  ELE_ARRAY(L_EA) = nL(Fl%T(ELEM,I),J)
                   END IF
                 END IF
              END DO
           END DO
           END IF
           DO I = 1,EE%n	! Loop over the basis for Function Space
              CALL GET_SPATIAL_COORDINATE(E,ELEM,EE%XI(I,1:3),P1)	! Get a spatial coordinate
              m = 0
              IF (EE%DISCONT == 0) THEN		! If the field is continuous
                DO J = 1,L_EA		! Loop over previous element list
                   DO K = 1,EE%n	! Loop over the basis for Function Space
                     CALL P_NORM(P1(1:3) - AA%X(AA%T(ELE_ARRAY(J),K),1:3),3,2,dL)
                     IF (dL < GLOBAL_TOL) THEN		! Is the basis unique?
                       m = AA%T(ELE_ARRAY(J),K); EXIT		! If no, set m to its redundancy
                     END IF
                   END DO
                   IF (m > 0) EXIT	! If a redundancy has been found, exit
                END DO
              END IF
              DO J = 1,I-1	! Loop over previously defined basis within ELEM
                 CALL P_NORM(P1(1:3) - AA%X(AA%T(ELEM,J),1:3),3,2,dL)
                 IF (dL < GLOBAL_TOL) THEN 	! Is the basis unique?
                   m = AA%T(ELEM,J); EXIT		! If no, set m to its redundancy
                 END IF
              END DO
              IF (m == 0) THEN		! If no redundancy exists, add basis to coefficient map
                AA%Lx = AA%Lx + 1; AA%X(AA%Lx,1:3) = P1(1:3)
                m = AA%Lx
              END IF
              AA%T(ELEM,I) = m		! Add basis to function space
           END DO
        END DO
      END SUBROUTINE MAKE_FUNCTION_SPACE


      SUBROUTINE INPUT_CORE
!       Subprogram INPUT_CORE - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program reads in all input files.
!
        USE VARS_PREP
        INTEGER I,J,K,m,n
        DOUBLE PRECISION Ax,Bx(1000),Cx
        CHARACTER*60 NAMZ
        CHARACTER*1 ID
        NIMZ = TRIM(iBASE); CALL ADD_EXTENSION(NIMZ)	! Add root pathway to initial basis file name
        CALL LOAD_BASE_SET(E)		! Read initial basis
        NAMZ = TRIM(COFILE)
        OPEN(UNIT=1,FILE=NAMZ,STATUS='old',action='read')	! Read initial coefficient file
          READ(1,*) Ax
          Fl%Lx = INT(Ax); ALLOCATE(Fl%X(Fl%Lx,3))
          CALL READ_MESH_CORED(Fl%X,Fl%Lx,3)
        CLOSE(1)
        NAMZ = TRIM(MAPFILE); Fl%ID = 0
        OPEN(UNIT=1,FILE=NAMZ,STATUS='old',action='read')	! Find basis call letter
          READ(1,*,ERR=45) Ax, ID
 45       DO I = 1,E%n_B
             IF (INDEX(ID,E%B(I)%CL) == 1) Fl%ID = I	! If call letter found, set ID = I
          END DO
          IF (Fl%ID == 0) THEN		! If the call letter is not found, search basis for M
             DO I = 1,E%n_B
                IF (INDEX(E%B(I)%CL,'M') == 1) Fl%ID = I
             END DO
             IF (Fl%ID == 0) CALL ERRORZ(2)	! If no basis is found, exit
          END IF
        CLOSE(1)
        OPEN(UNIT=1,FILE=NAMZ,STATUS='old',action='read')	! Read initial spatial map
          READ(1,*) Ax
          Fl%Lt = INT(Ax); ALLOCATE(Fl%T(Fl%Lt,E%B(Fl%ID)%n))
          CALL READ_MESH_COREI(Fl%T,Fl%Lt,E%B(Fl%ID)%n)
        CLOSE(1)
        NIMZ = TRIM(fBASE); CALL ADD_EXTENSION(NIMZ)	! Add root pathway to final basis file name
        Fl%Lb = 0
        IF (BNDRYFILE .NE. '') THEN
           NAMZ = TRIM(BNDRYFILE)
           OPEN(UNIT=1,FILE=NAMZ,STATUS='old',action='read')	! Find basis call letter
             READ(1,*) Ax
             Fl%Lb = INT(Ax); ALLOCATE(Fl%B(Fl%Lb, E%FNODES + 2  ))
             CALL READ_MESH_COREI(Fl%B,Fl%Lb, E%FNODES + 2  )
           CLOSE(1)
        END IF
        CALL LOAD_BASE_SET(E2)		! Read initial basis
        ALLOCATE(Ml(E2%n_B))		! Allocate no. of Function spaces
        CALL PRINTZ(2,0)
      END SUBROUTINE INPUT_CORE


      SUBROUTINE READ_MESH_COREI(A,n,m)
!       Subprogram READ_MESH_COREI - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program reads in integer arrays for file 1.
!
        INTEGER I,n,m,A(n,m)
        DOUBLE PRECISION Bx(m)
        DO I = 1,n
          READ(1,*,END=30) Bx(1:m); A(I,1:m) = INT(Bx(1:m))
        END DO
        RETURN
 30     CALL ERRORZ(4)
      END SUBROUTINE READ_MESH_COREI


      SUBROUTINE READ_MESH_CORED(A,n,m)
!       Subprogram READ_MESH_CORED - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program reads in double arrays for file 1.
!
        INTEGER I,n,m
        DOUBLE PRECISION A(n,m)
        DO I = 1,n
          READ(1,*,END=35) A(I,1:m)
        END DO
        RETURN
 35     CALL ERRORZ(4)
      END SUBROUTINE READ_MESH_CORED


      SUBROUTINE OUTPUT_FILE
!       Subprogram OUTPUT_FILE - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program reads in all output files.
!
        USE VARS_PREP
        INTEGER I,J,K,m,n
        DOUBLE PRECISION Ax,Bx(100),Cx
        CHARACTER*30 NAMZ
        CALL PRINTZ(4,0)
        NAMZ = TRIM(FROOT)//'_FE.C'
        OPEN(UNIT=1,FILE=NAMZ,STATUS='unknown')		! Output updated coefficient file
          WRITE(1,*) Ml(1:E2%n_B)%Lx
          DO J = 1,E2%n_B
            DO I = 1,Ml(J)%Lx
              WRITE(1,*) Ml(J)%X(I,1:3)
            END DO
          END DO
        CLOSE(1)
        NAMZ = TRIM(FROOT)//'_FE.M'
        OPEN(UNIT=1,FILE=NAMZ,STATUS='unknown')		! Output updated map file
          WRITE(1,*) Ml(1:E2%n_B)%Lt, E2%B(1:E2%n_B)%CL
          DO J = 1,E2%n_B
            DO I = 1,Ml(J)%Lt
              WRITE(1,*) Ml(J)%T(I,1:E2%B(J)%n)
            END DO
          END DO
        CLOSE(1)
        NAMZ = TRIM(FROOT)//'_FE.NBT'
        OPEN(UNIT=1,FILE=NAMZ,STATUS='unknown')		! Output updated neighboring file
          WRITE(1,*) Fl%Lt
          DO I = 1,Fl%Lt
            WRITE(1,*) Fl%NBT(I,1:E%FACES)
          END DO
        CLOSE(1)
        DO K = 1,E2%n_B
          IF (Ml(K)%Lb > 0) THEN
            NAMZ = TRIM(FROOT)//'_FE.B'
            OPEN(UNIT=1,FILE=NAMZ,STATUS='unknown')		! Output updated neighboring file
              WRITE(1,*) Ml(K)%Lb
              DO I = 1,Ml(K)%Lb
                WRITE(1,*) Ml(K)%B(I,1:E2%FNODES+2)
              END DO
            CLOSE(1)
          END IF
        END DO
      END SUBROUTINE OUTPUT_FILE


      SUBROUTINE LOAD_BASE_SET(AA)
!       Subprogram LOAD_BASE_SET - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program reads in all input basis for base AA.
!
        USE VARS_PREP
        INTEGER I,STP,D,FLD
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
        DO FLD = 1,AA%n_B		! Loop over number of basis
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')	! Read dimensions of basis CL
          DO WHILE (0 < 1)
              READ(1,*,END=70) IN_CHAR
              IF (INDEX(IN_CHAR,'no_basis_'//TRIM(AA%B(FLD)%CL)//'!') == 1) READ(1,*) AA%B(FLD)%n
              IF (INDEX(IN_CHAR,'dim_field_'//TRIM(AA%B(FLD)%CL)//'!') == 1) READ(1,*) AA%B(FLD)%DM
            END DO
 70     CLOSE(1)
        AA%B(FLD)%nl = AA%B(FLD)%n
        IF (AA%HEXA == 0) THEN		! If tensor input is indicated
           D = 3
           ALLOCATE(AA%B(FLD)%Q(AA%B(FLD)%n,AA%B(FLD)%n),AA%B(FLD)%P(3,AA%B(FLD)%n))
           ALLOCATE( AA%B(FLD)%XI(AA%B(FLD)%n,3) )
        ELSE				! If full input is indicated
           D = 1
           AA%B(FLD)%n  = AA%B(FLD)%n**REAL(AA%DM)
           ALLOCATE(AA%B(FLD)%Q(AA%B(FLD)%nl,AA%B(FLD)%nl),AA%B(FLD)%P(D,AA%B(FLD)%nl))
           ALLOCATE( AA%B(FLD)%B_ID(AA%B(FLD)%n,3), AA%B(FLD)%XI(AA%B(FLD)%n,3) )
        END IF
        OPEN(UNIT=1,FILE=NIMZ,STATUS='old',action='read')		! Read in Basis
          DO WHILE (0 < 1)
              READ(1,*,END=80) IN_CHAR
              IF (INDEX(IN_CHAR,TRIM(AA%B(FLD)%CL)//'_basis!') == 1) THEN
                 DO I = 1,AA%B(FLD)%nl
                 READ(1,*) AA%B(FLD)%Q(I,1:AA%B(FLD)%nl)
                 END DO
              ELSE IF (INDEX(IN_CHAR,TRIM(AA%B(FLD)%CL)//'_basis_discontinuous!') == 1) THEN
                 AA%B(FLD)%DISCONT = 1
              ELSE IF (INDEX(IN_CHAR,TRIM(AA%B(FLD)%CL)//'_xi_coordinates!') == 1) THEN
                 DO I = 1,AA%B(FLD)%nl
                 READ(1,*) AA%B(FLD)%XI(I,1:D)
                 END DO
              ELSE IF (INDEX(IN_CHAR,TRIM(AA%B(FLD)%CL)//'_P_basis!') == 1) THEN
                 DO I = 1,D
                 READ(1,*) AA%B(FLD)%P(I,1:AA%B(FLD)%nl)
                 END DO
              ELSE IF (INDEX(IN_CHAR,'basis_ordering_'//TRIM(AA%B(FLD)%CL)//'!') == 1) THEN
                DO I = 1,AA%B(FLD)%n
                  READ(1,*) AA%B(FLD)%B_ID(I,1:3)
                END DO
              END IF
            END DO
 80     CLOSE(1)
        IF (AA%HEXA == 1) CALL BUILD_HEXA_XI(AA,FLD)		! Build tensor xi points
        CALL LOAD_INTEGRALS(AA,FLD)				! Calculate basis integrals
        END DO
        RETURN
 90     CALL ERRORZ(5)
      END SUBROUTINE LOAD_BASE_SET


      SUBROUTINE LOAD_INTEGRALS(AA,FLD)
!       Subprogram LOAD_BASE_SET - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This program loads integrals for base AA%B(FLD).  The basis is evaluated a given quadrature points.
!
        USE VARS_PREP
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K,x,y,z,ID,ORD1(3),FLD,S
        DOUBLE PRECISION Ax(3),Bx,Cx(3),Dx
        ALLOCATE( AA%B(FLD)%I%Y(AA%B(FLD)%n,AA%n_pts) )
        ALLOCATE( AA%B(FLD)%I%Y_f(AA%B(FLD)%n,AA%n_pts_f,AA%FACES) )
        ORD1(1:3) = 0
        DO I = 1,AA%B(FLD)%n		! Loop over basis
           DO J = 1,AA%n_pts		! Loop over quadrature points
              CALL POLY(AA%B(FLD),AA%HEXA,I,AA%gpt(J,1:3),ORD1,AA%B(FLD)%I%Y(I,J))	! Calculate basis
           END DO
        END DO
        ORD1(1:3) = 0
        DO S = 1,AA%FACES	! Loop over facets
          DO I = 1,AA%B(FLD)%n		! Loop over basis
            DO J = 1,AA%n_pts_f			! Loop over facet quadrature points
              CALL POLY(AA%B(FLD),AA%HEXA,I,AA%gpt_f(J,1:3,S),ORD1,AA%B(FLD)%I%Y_f(I,J,S))	! Calculate basis
            END DO
          END DO
        END DO
      END SUBROUTINE LOAD_INTEGRALS


      SUBROUTINE BUILD_HEXA_XI(AA,FLD)
!       Subprogram BUILD_HEXA_XI - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine builds the XI coordinate for the tensor basis AA%B(FLD)
!
        USE VARS_PREP
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K,S,m,n,G,FLD
        DOUBLE PRECISION Bx(AA%B(FLD)%nl)
        Bx(1:AA%B(FLD)%nl) = AA%B(FLD)%XI(1:AA%B(FLD)%nl,1)
        DO I = 1,AA%B(FLD)%n		! Loop over basis
           AA%B(FLD)%XI(I,1:3) = Bx(AA%B(FLD)%B_ID(I,1:3))	! calculate xi point
        END DO
      END SUBROUTINE BUILD_HEXA_XI


      SUBROUTINE BUILD_HEXA_QUADRATURE(AA)
!       Subprogram BUILD_HEXA_QUADRATURE - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine builds the volume/facet quadrature rule for the tensor basis AA
!
        USE VARS_PREP
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
        USE VARS_PREP
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I,J,K
        DOUBLE PRECISION Ax(3),Bx,Cx(3),Dx,ONN
        ONN = 1.0
        AA%gpt_f(1:AA%n_pts_f,AA%DM:3,1) = 0		! Initialize first quadrature rule
        DO I = 1,AA%FACES
           AA%nrm(1:3,I) = 0		! Initialize normal vectors
        END DO
        IF (AA%DM == 3) THEN		! If in R^3
         IF (AA%FACES == 4) THEN	! master element is a tetrahedron, build appropriate facet quadrature
          AA%nrm(3,1) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,2) = 0
          AA%gpt_f(1:AA%n_pts_f,3,2) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(2,2) = -ONN


          AA%gpt_f(1:AA%n_pts_f,1,3) = 0
          AA%gpt_f(1:AA%n_pts_f,2,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,3) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(1,3) = -ONN

          DO I = 1,AA%n_pts_f
          AA%gpt_f(I,1,4) = ONN-AA%gpt_f(I,1,1)-AA%gpt_f(I,2,1)
          END DO
          AA%gpt_f(1:AA%n_pts_f,2,4) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,4) = AA%gpt_f(1:AA%n_pts_f,2,1)
          Ax(1) = ONN
          Ax(2) = 3.0
          Ax(2) = Ax(2)**0.5
          AA%nrm(1:3,4) = Ax(1)/Ax(2)
         ELSE IF (AA%FACES == 6) THEN	! master element is a cube, build appropriate facet quadrature
          AA%nrm(3,1) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,2) = 0
          AA%gpt_f(1:AA%n_pts_f,3,2) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(2,2) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,3) = 0
          AA%gpt_f(1:AA%n_pts_f,2,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,3) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(1,3) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,4) = ONN
          AA%gpt_f(1:AA%n_pts_f,2,4) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,4) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(1,4) = ONN

          AA%gpt_f(1:AA%n_pts_f,1,5) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,5) = ONN
          AA%gpt_f(1:AA%n_pts_f,3,5) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%nrm(2,5) = ONN

          AA%gpt_f(1:AA%n_pts_f,1,6) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,6) = AA%gpt_f(1:AA%n_pts_f,2,1)
          AA%gpt_f(1:AA%n_pts_f,3,6) = ONN
          AA%nrm(3,6) = ONN
         END IF
        ELSE IF (AA%DM == 2) THEN	! If in R^2
         IF (AA%FACES == 3) THEN	! master element is a triangle, build appropriate facet quadrature
          AA%nrm(2,1) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,2) = 0
          AA%gpt_f(1:AA%n_pts_f,2,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,2) = 0
          AA%nrm(1,2) = -ONN

          DO I = 1,AA%n_pts_f
          AA%gpt_f(I,1,3) = ONN-AA%gpt_f(I,1,1)
          END DO
          AA%gpt_f(1:AA%n_pts_f,2,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,3) = 0
          Ax(1) = ONN
          Ax(2) = 2.0
          Ax(2) = Ax(2)**0.5
          AA%nrm(1:2,3) = Ax(1)/Ax(2); AA%nrm(3,3) = 0
         ELSE IF (AA%FACES == 4) THEN	! master element is a square, build appropriate facet quadrature
          AA%nrm(2,1) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,2) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,2) = 0
          AA%gpt_f(1:AA%n_pts_f,3,2) = 0
          AA%nrm(1,2) = -ONN

          AA%gpt_f(1:AA%n_pts_f,1,3) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,2,3) = ONN
          AA%gpt_f(1:AA%n_pts_f,3,3) = 0
          AA%nrm(2,3) = ONN

          AA%gpt_f(1:AA%n_pts_f,1,4) = ONN
          AA%gpt_f(1:AA%n_pts_f,2,4) = AA%gpt_f(1:AA%n_pts_f,1,1)
          AA%gpt_f(1:AA%n_pts_f,3,4) = 0
          AA%nrm(1,4) = ONN
         END IF
        END IF
      END SUBROUTINE BUILD_FACE_GPT_ARRAY


      SUBROUTINE GET_SPATIAL_COORDINATE(AA,ELEM,XI,Pt)
!       Subprogram GET_SPATIAL_COORDINATE - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine evaluates spatial tessellation (Fl) with basis (AA), at xi point (XI)
!       in element Pt.
!
        USE VARS_PREP
        TYPE(ARRAY_PROBLEM_BASE) AA
        INTEGER I, J, K, ELEM, ORD1(3), ID
        DOUBLE PRECISION XI(3), Pt(3), P1(3), Ax
        Pt(1:3) = 0
        ORD1(1:3) = 0
        DO K = 1,AA%B(Fl%ID)%n	! Loop over basis
           CALL POLY(AA%B(Fl%ID),AA%HEXA,K,XI,ORD1,Ax)	! calculate basis K
           Pt(1:3) = Pt(1:3) + Fl%X(Fl%T(ELEM,K),1:3)*Ax	! Sum and add coefficient weighting
        END DO
      END SUBROUTINE GET_SPATIAL_COORDINATE


      SUBROUTINE POLY(AA,HEXA,BAS,X,ORD1,VAL)
!       Subprogram POLY - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine evaluates basis (BAS) of (AA) with derivatives of order (ORD1) in each spatial
!       dimension.  The answer is provided in VAL
!
        USE VARS_PREP
        TYPE(ARRAY_BASE) AA
        INTEGER n,ORD1(3),I,J,P,BAS,K,HEXA
        DOUBLE PRECISION VAL,X(3)
        IF (HEXA == 1) THEN	! If tensor basis
           CALL POLY_HEX(AA,BAS,X,ORD1,VAL)
        ELSE			! If general basis
           CALL POLY_TET(AA,BAS,X,ORD1,VAL)
        END IF
      END SUBROUTINE POLY


      SUBROUTINE POLY_HEX(AA,BAS,X,ORD1,VAL)
!       Subprogram POLY_HEX - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine evaluates basis (BAS) of (AA) with derivatives of order (ORD1) in each spatial
!       dimension.  The answer is provided in VAL
!
        USE VARS_PREP
        TYPE(ARRAY_BASE) AA
        INTEGER n,ORD1(3),I,J,P,BAS,K
        DOUBLE PRECISION VAL,Ax,Bx,Cx,Bo,X(3),Co(3)
        Co(1:3) = 0
        DO I = 1,AA%nl		! Loop over basis in a given coordinate dimension
           DO J = 1,AA%DM		! Loop over dimension
              IF (ORD1(J) <= AA%P(1,I)) THEN	! If power is great enough, that term is non-zero after diff
                 Ax = 1
                 DO K = 0,ORD1(J)-1	! Calculate the derivative of the term
                    Ax = Ax*(REAL(AA%P(1,I))-REAL(K))
                 END DO
                 Ax = Ax*(X(J)**(AA%P(1,I)-ORD1(J)))
                 Co(J) = Co(J) + AA%Q(AA%B_ID(BAS,J),I)*Ax	! add to current direction weighted by c.tensor
              END IF
           END DO
        END DO
        VAL = Co(1)
        DO I = 2,AA%DM		! Do tensor product
          VAL = VAL*Co(I)
        END DO
      END SUBROUTINE POLY_HEX


      SUBROUTINE POLY_TET(AA,BAS,X,ORD1,VAL)
!       Subprogram POLY_TET - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine evaluates basis (BAS) of (AA) with derivatives of order (ORD1) in each spatial
!       dimension.  The answer is provided in VAL
!
        USE VARS_PREP
        TYPE(ARRAY_BASE) AA
        INTEGER n,ORD1(3),I,J,P,BAS,K
        DOUBLE PRECISION VAL,Ax,Bx,Cx,Bo,X(3),Co(3)
        VAL = 0
        DO I = 1,AA%n			! Loop over basis 
           Co(1:3) = 0; Bx = 1_DPR
           DO J = 1,AA%DM		! Loop over dimension
             IF (ORD1(J) <= AA%P(J,I)) THEN	! If power is great enough, that term is non-zero after diff
              Ax = 1
              DO K = 0,ORD1(J)-1	! Calculate the derivative of the term
                 Ax = Ax*(REAL(AA%P(J,I))-REAL(K))
              END DO
              Ax = Ax*(X(J)**(AA%P(J,I)-ORD1(J)))
              Co(J) = Ax
             END IF
             Bx = Bx*Co(J)	! calculate product of term
           END DO
           VAL = VAL + AA%Q(BAS,I)*Bx	! add to current value weighted by c.tensor
        END DO
      END SUBROUTINE POLY_TET


      SUBROUTINE P_NORM(A,n,p,dL)
!       Subprogram P_NORM - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine evaluates dL = l_p (A(1:n)).  If 0, its the l_infty norm.
!
        USE VARSKINDS
        INTEGER n,I,p
        DOUBLE PRECISION A(n), dL
        dL = 0
        IF (p.NE.0) THEN	! Not zero, calculate norm
           DO I = 1,n
           dL = dL + (abs(A(I)))**REAL(p)
           END DO
           dL = dL**(1_DPR / REAL(p))
        ELSE			! Calculate l_infty norm
           DO I = 1,n
              dL = MAX(dL,abs(A(I)))
           END DO
        END IF
      END SUBROUTINE P_NORM


      SUBROUTINE ADD_EXTENSION(NAMZ)
!       Subprogram ADD_EXTENSION - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine adds pathway to basis files to NAMZ
!
        USE VARS_PREP
        CHARACTER*90 NAMZ
        IF (TRIM(EXTN).NE.'') THEN
            NAMZ = TRIM(EXTN)//TRIM(NAMZ)
        END IF
      END SUBROUTINE ADD_EXTENSION


      SUBROUTINE COMM_LINE_IO
!       Subprogram COMM_LINE_IO - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine gets the command line arguements for the preparatory routines
!
        USE VARS_PREP
        INTEGER I,J,K,CNT, n
        CHARACTER*60 NAMZ2
        CNT = COMMAND_ARGUMENT_COUNT(); I = 0; n = 0; BNDRYFILE = ''
        DO WHILE (I < CNT)
           I = I + 1
           CALL GET_COMMAND_ARGUMENT(I,NAMZ2,J,K)
           IF (INDEX(NAMZ2,'-iM') == 1) THEN		! If the initial map flag
               I = I + 1; n = n + 1
               CALL GET_COMMAND_ARGUMENT(I,MAPFILE,J,K)
               MAPFILE = TRIM(MAPFILE)
           ELSE IF (INDEX(NAMZ2,'-iC') == 1) THEN	! If the initial coefficient flag
               I = I + 1; n = n + 1
               CALL GET_COMMAND_ARGUMENT(I,COFILE,J,K)
               COFILE = TRIM(COFILE)
           ELSE IF (INDEX(NAMZ2,'-iGAMMA') == 1) THEN
               I = I + 1
               CALL GET_COMMAND_ARGUMENT(I,BNDRYFILE,J,K)
               BNDRYFILE = TRIM(BNDRYFILE)
           ELSE IF (INDEX(NAMZ2,'-iB') == 1) THEN	! If the initial basis flag
                I = I + 1; n = n + 1
               CALL GET_COMMAND_ARGUMENT(I,iBASE,J,K)
               iBASE = TRIM(iBASE)
           ELSE IF (INDEX(NAMZ2,'-fB') == 1) THEN	! If the final basis flag
                I = I + 1; n = n + 1
               CALL GET_COMMAND_ARGUMENT(I,fBASE,J,K)
               fBASE = TRIM(fBASE)
           ELSE IF (INDEX(NAMZ2,'-fO') == 1) THEN	! If the output flag
                I = I + 1; n = n + 1
               CALL GET_COMMAND_ARGUMENT(I,FROOT,J,K)
               FROOT = TRIM(FROOT)
           END IF
        END DO
        CALL PRINTZ(1,0)
        IF (n < 5) CALL ERRORZ(1)
      END SUBROUTINE COMM_LINE_IO


      SUBROUTINE PRINTZ(I,K)
!       Subprogram PRINTZ - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine prints messages and status for the preparatory routines
!
        USE VARS_PREP
        INTEGER I,K
        IF (I == 1) THEN
          PRINT *, ' ==>>  COMMAND LINE INPUTS COLLECTED  <<== '
          PRINT *, ' ==>>  MAP FILE  <<== '
          PRINT *, ' ----  '//TRIM(MAPFILE)
          PRINT *, ' ==>>  COEFFICIENT FILE  <<== '
          PRINT *, ' ----  '//TRIM(COFILE)
          PRINT *, ' ==>>  INITIAL BASIS  <<== '
          PRINT *, ' ----  '//TRIM(iBASE)
          PRINT *, ' ==>>  FINAL BASIS  <<== '
          PRINT *, ' ----  '//TRIM(fBASE)
          PRINT *, ' ==>>  ROOT OUTPUT FILE NAME  <<== '
          PRINT *, ' ----  '//TRIM(FROOT)
        END IF
        IF (I == 2) PRINT *, ' ==>>  STATUS :: INPUT FILES READ  <<=='
        IF (I == 3) THEN
          PRINT *, ' ==>>  STATUS :: SPACE '//TRIM(E2%B(K)%CL)//' COMPLETE  <<== '
        END IF
        IF (I == 4) THEN
          PRINT *, ' ==>>  OUTPUT FILES  <<== '
          PRINT *, ' ==>>  UPDATED MAP FILE  <<== '
          PRINT *, ' ----  '//TRIM(FROOT)//'_FE.T'
          PRINT *, ' ==>>  UPDATED COEFFICIENT FILE  <<== '
          PRINT *, ' ----  '//TRIM(FROOT)//'_FE.X'
          PRINT *, ' ==>>  NEIGHBORING FILE  <<== '
          PRINT *, ' ----  '//TRIM(FROOT)//'_FE.NBT'
        END IF

        IF (I == 100) PRINT *, ' ==>>  STATUS :: DONE  <<=='
      END SUBROUTINE PRINTZ


      SUBROUTINE ERRORZ(I)
!       Subprogram ERRORZ - Written by David A. Nordsletten (C) 2008
!       DPhil Student, Comlab, University of Oxford 2006-2008
!       PhD Student, Bioengineering Institute, University of Auckland 2005-2006
!
!       This routine prints error messages and stops execution for the preparatory routines
!
        USE VARS_PREP
        INTEGER I,K
        IF (I == 1) PRINT *, ' ==>>  FORCE QUIT :: INCOMPLETE COMMAND LINE INPUT  <<== '
        IF (I == 2) PRINT *, ' ==>>  FORCE QUIT :: NO CALL LETTER or M BASIS  <<== '
        IF (I == 3) PRINT *, ' ==>>  FORCE QUIT :: NO FIELDS IN BASIS FILE <<== '
        IF (I == 4) PRINT *, ' ==>>  FORCE QUIT :: INCOMPLETE MAP or COEFFICIENT FILE(S) <<== '
        IF (I == 5) PRINT *, ' ==>>  FORCE QUIT :: PROBLEM ENCOUNTERED READING BASIS FILE <<== '
        IF (I == 6) PRINT *, ' ==>>  FORCE QUIT :: MULTIPLE ELEMENTS POINTING TO A SINGLE FACET <<== '
        PRINT *, ' ==>> QUITING <<== '
        STOP
      END SUBROUTINE ERRORZ
