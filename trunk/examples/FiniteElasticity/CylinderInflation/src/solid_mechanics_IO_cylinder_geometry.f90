!>creates automatically a cylinder geometry according to given level of refinements
MODULE SOLID_MECHANICS_IO_CYLINDER_GEOMETRY
  USE CONSTANTS
  USE COMP_ENVIRONMENT
  USE KINDS  
  USE SOLID_MECHANICS_IO_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  ! REMEMBER, THERE ARE THREE COORDINATE DIRECTIONS: (R, THETA, Z)

  !>cylinder type parameters
  TYPE SOLID_MECHANICS_IO_CYLINDER
    REAL(DP) :: R1         !<inner radius
    REAL(DP) :: R2         !<outer radius
    REAL(DP) :: Z0         !<length of cylinder
    INTEGER(INTG) :: NER   !<number of elements in r direction
    INTEGER(INTG) :: NET   !<number of elements in circumferential (theta) direction
    INTEGER(INTG) :: NEZ   !<number of elements in z direction
    !
    REAL(DP) :: DR,DT,DZ   !<extent of a single element in each dimension (use Radians)
    INTEGER(INTG),ALLOCATABLE :: INNER_FACE(:)   !<nodes on the inner face
    INTEGER(INTG),ALLOCATABLE :: OUTER_FACE(:)   !<nodes on the outer face
    INTEGER(INTG),ALLOCATABLE :: TOP_INNER(:)    !<nodes on the top inner ring
    INTEGER(INTG),ALLOCATABLE :: TOP_OUTER(:)    !<nodes on the top outer ring
    INTEGER(INTG),ALLOCATABLE :: BOTTOM_INNER(:) !<nodes on the bottom inner ring
    INTEGER(INTG),ALLOCATABLE :: BOTTOM_OUTER(:) !<nodes on the bottom outer ring
    INTEGER(INTG),ALLOCATABLE :: TOP_FACE(:)     !<nodes on the top face of the cylinder
    INTEGER(INTG),ALLOCATABLE :: BOTTOM_FACE(:)  !<nodes on the bottom face of the cylinder
    INTEGER(INTG) :: NODE_ON_X_AXIS  !<a node that sits on x-axis initially
    INTEGER(INTG) :: NODE_ON_Y_AXIS  !<a node that sits on y-axis initially
  END TYPE

  ! MODULE PROCEDURES
  PUBLIC :: SOLID_MECHANICS_IO_CREATE_CYLINDER
  PUBLIC :: SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS, SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH
  PUBLIC :: SOLID_MECHANICS_IO_CREATE_BC

  CONTAINS

  ! =============================================================================
  !>main routine for creating a cylinder mesh
  SUBROUTINE SOLID_MECHANICS_IO_CREATE_CYLINDER(CYLINDER,MESH,ERR,ERROR,*,ARGS)
    TYPE(SOLID_MECHANICS_IO_CYLINDER),INTENT(OUT) :: CYLINDER     !<optional cylinder struct
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(OUT) :: MESH   !<cylinder mesh container
    INTEGER(INTG),INTENT(OUT) :: ERR                              !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                     !<The error string
    REAL(DP),INTENT(IN),OPTIONAL :: ARGS(:)                       !<Contains command-line arguments
    ! local variables
    CHARACTER*30 :: WORD
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
 
    CALL ENTERS("SOLID_MECHANICS_IO_CREATE_CYLINDER",ERR,ERROR,*999)

    ! not the best place to do this but..
    MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    IF(ERR/=0) GOTO 999

    ! CHECK IF THE ARGUMENTS WERE SPECIFIED FROM COMMAND LINE
    IF(PRESENT(ARGS)) THEN
      CYLINDER%R1=ARGS(1)
      CYLINDER%R2=ARGS(2)
      CYLINDER%Z0=ARGS(3)
      CYLINDER%NER=ARGS(4)
      CYLINDER%NET=ARGS(5)
      CYLINDER%NEZ=ARGS(6)
      ! element sizes
      CYLINDER%DR = (CYLINDER%R2-CYLINDER%R1)/REAL(CYLINDER%NER,DP)
      CYLINDER%DT = TWOPI/REAL(CYLINDER%NET,DP)
      CYLINDER%DZ = CYLINDER%Z0/REAL(CYLINDER%NEZ,DP)
    ELSE
      CALL SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS(CYLINDER,ERR,ERROR,*999)
    ENDIF
    ! TODO: GET MESH%ND(:) FROM BASIS PROPERLY (CURRENTLY, ONLY TRILINEAR ALLOWED)
    ALLOCATE(MESH%ND(3))
    MESH%ND=1   ! HARDCODED
    CALL SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH(CYLINDER,MESH,ERR,ERROR,*999)

    CALL EXITS("SOLID_MECHANICS_IO_CREATE_CYLINDER")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_CREATE_CYLINDER",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_CREATE_CYLINDER")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_CREATE_CYLINDER

  ! =============================================================================
  !>user interface routine to read in parameters for a cylinder
  SUBROUTINE SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS(CYLINDER,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_CYLINDER),INTENT(OUT) :: CYLINDER
    INTEGER(INTG),INTENT(OUT) :: ERR                            !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                   !<The error string    
    ! local variables
    CHARACTER*100 :: WORD

    CALL ENTERS("SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS",ERR,ERROR,*999)

    WRITE(*,*) "*** Enter parameters for cylinder ***"

    ! read in some parameters for the cylinder
    CYLINDER%R1 = 0.5_DP
    DO
      WRITE(*,*) "Enter the inner radius (R1) [0.5]: "
      READ(*,'(A)') WORD
      IF(LEN_TRIM(WORD)==0) EXIT
      CYLINDER%R1 = STRING_TO_DOUBLE(WORD,ERR,ERROR)
      EXIT
    ENDDO

    CYLINDER%R2 = 1.0_DP
    DO
      WRITE(*,*) "Enter the outer radius (R2) [1.0]: "
      READ(*,'(A)') WORD
      IF(LEN_TRIM(WORD)==0) EXIT
      CYLINDER%R2 = STRING_TO_DOUBLE(WORD,ERR,ERROR)
      EXIT
    ENDDO

    CYLINDER%Z0 = 0.1_DP
    DO
      WRITE(*,*) "Enter the height of cylinder (Z0) [0.1]: "
      READ(*,'(A)') WORD
      IF(LEN_TRIM(WORD)==0) EXIT
      CYLINDER%Z0 = STRING_TO_DOUBLE(WORD,ERR,ERROR)
      EXIT
    ENDDO

    CYLINDER%NER = 1
    DO
      WRITE(*,*) "Enter the number of elements in radial direction [1]: "
      READ(*,'(A)') WORD
      IF(LEN_TRIM(WORD)==0) EXIT
      CYLINDER%NER = STRING_TO_INTEGER(WORD,ERR,ERROR)
      EXIT
    ENDDO

    CYLINDER%NET = 4
    DO
      WRITE(*,*) "Enter the number of elements in circumferential direction [4]: "
      READ(*,'(A)') WORD
      IF(LEN_TRIM(WORD)==0) EXIT
      CYLINDER%NET = STRING_TO_INTEGER(WORD,ERR,ERROR)
      EXIT
    ENDDO

    CYLINDER%NEZ = 1
    DO
      WRITE(*,*) "Enter the number of elements in axial direction [1]: "
      READ(*,'(A)') WORD
      IF(LEN_TRIM(WORD)==0) EXIT
      CYLINDER%NEZ = STRING_TO_INTEGER(WORD,ERR,ERROR)
      EXIT
    ENDDO

    ! element sizes
    CYLINDER%DR = (CYLINDER%R2-CYLINDER%R1)/REAL(CYLINDER%NER,DP)
    CYLINDER%DT = TWOPI/REAL(CYLINDER%NET,DP)
    CYLINDER%DZ = CYLINDER%Z0/REAL(CYLINDER%NEZ,DP)

    CALL EXITS("SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS")
    RETURN

999 CALL ERRORS("SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_GET_CYLINDER_PARAMS

  ! =============================================================================
  !>create a cylinder mesh from user-specified parameters and find external faces
  !>(returns a SOLID_MECHANICS_IO_MESH_CONTAINER type)
  SUBROUTINE SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH(CYLINDER,MESH,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_CYLINDER),INTENT(INOUT) :: CYLINDER
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(INOUT) :: MESH !<need to have ND(:) assigned
    INTEGER(INTG),INTENT(OUT) :: ERR                            !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                   !<The error string
    ! local variables
    INTEGER(INTG),ALLOCATABLE :: NIDX(:,:,:)    ! holds the node numbers for given (r,th,z)
    INTEGER(INTG) :: NR,NT,NZ   ! dimensions of NIDX
    INTEGER(INTG) :: I,J,K,R,T,Z  ! loop indices
    INTEGER(INTG) :: NN         ! number of nodes
    INTEGER(INTG) :: NE         ! number of elements
    REAL(DP) :: RAD,THETA       ! actual values of coordinates
    REAL(DP) :: DUM

    CALL ENTERS("SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH",ERR,ERROR,*999)
    ! populate the NIDX array
    NR=CYLINDER%NER+1
    NT=CYLINDER%NET      ! NOTE SPECIAL CASE
    NZ=CYLINDER%NEZ+1
    ALLOCATE(NIDX(NR,NT,NZ))

    ! allocate NODES
    NN=NR*NT*NZ ! TOTAL NUMBER OF NODES
    ALLOCATE(MESH%NODES(NN))
    DO I=1,NN
      IF(MESH%ND(1)/=1) CALL FLAG_ERROR("only implemented for trilinear element type at current time",ERR,ERROR,*999)
      ALLOCATE(MESH%NODES(I)%X(MESH%ND(1)),MESH%NODES(I)%Y(MESH%ND(2)),MESH%NODES(I)%Z(MESH%ND(3)))
    ENDDO

    ! assign NIDX & NODES
    NN=0
    DO Z=1,NZ
      DO T=1,NT
        DO R=1,NR
          ! NIDX
          NN=NN+1
          NIDX(R,T,Z)=NN
          ! NODES
          IF(MESH%ND(1)==1.AND.MESH%ND(2)==1.AND.MESH%ND(3)==1) THEN
            MESH%NODES(NN)%ID=NN
            RAD = CYLINDER%R1 + REAL(R-1,DP)*CYLINDER%DR
            THETA = REAL(T-1,DP)*CYLINDER%DT
            MESH%NODES(NN)%X(1) = RAD*COS(THETA)
            MESH%NODES(NN)%Y(1) = RAD*SIN(THETA)
            MESH%NODES(NN)%Z(1) = REAL(Z-1,DP)*CYLINDER%DZ
          ELSE
            CALL FLAG_ERROR("only implemented for trilinear element type at current time",ERR,ERROR,*999)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! cool. now allocate ELS
    NE=CYLINDER%NER*CYLINDER%NET*CYLINDER%NEZ
    ALLOCATE(MESH%ELS(NE))
    DO I=1,NE
      ALLOCATE(MESH%ELS(I)%NODES(8)) ! again, 3D trilinear or tricubic assumed
    ENDDO

    ! assign ELS
    NE=0
    DO Z=1,CYLINDER%NEZ
      DO T=1,CYLINDER%NET
        DO R=1,CYLINDER%NER
          NE=NE+1
          IF(T<CYLINDER%NET) THEN
            MESH%ELS(NE)%NODES = (/NIDX(R,T,Z),NIDX(R+1,T,Z),NIDX(R,T+1,Z),NIDX(R+1,T+1,Z), &
            & NIDX(R,T,Z+1),NIDX(R+1,T,Z+1),NIDX(R,T+1,Z+1),NIDX(R+1,T+1,Z+1)/)
          ELSE
            ! the last element wraps around, so some adjustment required
            MESH%ELS(NE)%NODES = (/NIDX(R,T,Z),NIDX(R+1,T,Z),NIDX(R,1,Z),NIDX(R+1,1,Z), &
            & NIDX(R,T,Z+1),NIDX(R+1,T,Z+1),NIDX(R,1,Z+1),NIDX(R+1,1,Z+1)/)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! update MESH to return
    MESH%NN=NN
    MESH%NE=NE
    MESH%NC=3   ! TODO: HARDCODED 3 DIMENSIONS
    MESH%FINISHED=.TRUE.

    !-----------------------
    ! let's cram in external/internal face list buildling here
    ALLOCATE(CYLINDER%INNER_FACE(NT*(NZ)))
    ALLOCATE(CYLINDER%OUTER_FACE(NT*(NZ)))
    ALLOCATE(CYLINDER%TOP_INNER(NT))
    ALLOCATE(CYLINDER%TOP_OUTER(NT))
    ALLOCATE(CYLINDER%BOTTOM_INNER(NT))
    ALLOCATE(CYLINDER%BOTTOM_OUTER(NT))
    ALLOCATE(CYLINDER%TOP_FACE(NR*NT))
    ALLOCATE(CYLINDER%BOTTOM_FACE(NR*NT))
    I=0
    J=0
    K=0
    DO T=1,NT
      I=I+1
      CYLINDER%TOP_INNER(I)=NIDX(1,T,NZ)
      CYLINDER%TOP_OUTER(I)=NIDX(NR,T,NZ)
      CYLINDER%BOTTOM_INNER(I)=NIDX(1,T,1)
      CYLINDER%BOTTOM_OUTER(I)=NIDX(NR,T,1)
      DO Z=1,NZ
        ! these should end up in vertical groups of equal theta
        J=J+1
        CYLINDER%INNER_FACE(J)=NIDX(1,T,Z)
        CYLINDER%OUTER_FACE(J)=NIDX(NR,T,Z)
      ENDDO
      DO R=1,NR
        K=K+1
        CYLINDER%TOP_FACE(K)=NIDX(R,T,NZ)
        CYLINDER%BOTTOM_FACE(K)=NIDX(R,T,1)
      ENDDO
    ENDDO

    ! node on x axis
    CYLINDER%NODE_ON_X_AXIS = 1
    ! node on y axis
    DUM=(PI/2)/CYLINDER%DT
    IF(ABS(DUM-ANINT(DUM))<1E-8_DP) THEN
      CYLINDER%NODE_ON_Y_AXIS = NIDX(1,NINT(DUM)+1,1)
    ELSE
      ! there's no node on the y axis!!
      CALL FLAG_ERROR("Choose an even number of elements in circumferential direction.",ERR,ERROR,*999)
    ENDIF

    DEALLOCATE(NIDX)

    CALL EXITS("SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH")
    RETURN
  END SUBROUTINE SOLID_MECHANICS_IO_CREATE_CYLINDER_MESH

  ! =============================================================================
  !>creates a bc object from given inner/outer pressures and extension ratio
  SUBROUTINE SOLID_MECHANICS_IO_CREATE_BC(CYLINDER,BC,ERR,ERROR,*,PROB,ARGS)
    TYPE(SOLID_MECHANICS_IO_CYLINDER),INTENT(IN) :: CYLINDER
    TYPE(SOLID_MECHANICS_IO_BOUNDARY_CONDITION),ALLOCATABLE,INTENT(OUT) :: BC(:)
    INTEGER(INTG),INTENT(OUT) :: ERR            !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR   !<The error string
    REAL(DP),INTENT(OUT),OPTIONAL :: PROB(3)    !<p_inner, p_outer and lambda
    REAL(DP),INTENT(IN),OPTIONAL :: ARGS(:)     !<contains command-line arguments
    ! local variables
    REAL(DP) :: P_INNER,P_OUTER,LAMBDA
    REAL(DP) :: FACE_AREA_INNER, FACE_AREA_OUTER
    REAL(DP) :: FORCE_INNER, FORCE_OUTER
    REAL(DP) :: FORCE_X, FORCE_Y
    INTEGER(INTG) :: NBC,BC_IDX,NEZ,I,J
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
    CHARACTER*30 :: WORD

    CALL ENTERS("SOLID_MECHANICS_IO_CREATE_BC",ERR,ERROR,*999)

    ! CHECK IF THE ARGUMENTS WERE SPECIFIED FROM COMMAND LINE
    IF(PRESENT(ARGS)) THEN
      P_INNER=ARGS(7)
      P_OUTER=ARGS(8)
      LAMBDA=ARGS(9)
    ELSE
      ! ask for inner and outer pressures
      WRITE(*,*) "Enter the value of inner pressure [1.0]: "
      P_INNER=1.0_DP
      DO
        READ(*,'(A)') WORD
        IF(LEN_TRIM(WORD)==0) EXIT
        P_INNER=STRING_TO_DOUBLE(WORD,ERR,ERROR)
        EXIT
      ENDDO

      WRITE(*,*) "Enter the value of outer pressure [0.0]: "
      P_OUTER=0.0_DP
      DO
        READ(*,'(A)') WORD
        IF(LEN_TRIM(WORD)==0) EXIT
        P_OUTER=STRING_TO_DOUBLE(WORD,ERR,ERROR)
        EXIT
      ENDDO

      ! read the axial stretch ratio
      WRITE(*,*) "Enter the value of stretch ratio [1.1]: "
      LAMBDA=1.1_DP
      DO
        READ(*,'(A)') WORD
        IF(LEN_TRIM(WORD)==0) EXIT
        LAMBDA=STRING_TO_DOUBLE(WORD,ERR,ERROR)
        EXIT
      ENDDO
    ENDIF

    PROB(1)=P_INNER
    PROB(2)=P_OUTER
    PROB(3)=LAMBDA

    ! The boundary condtions are:
    ! X-Y: inner and outer surface pressure boundary condition
    ! Z: top and bottom dirichlet condition
    ! rigid body motion: need two points to restrict X and a Y
    NEZ=CYLINDER%NEZ
    ! inner face X&Y (NET), outer face X&Y (NET), top face (1), bottom face (1) , rigid body motion (2)
    NBC=4*CYLINDER%NET+4
    ALLOCATE(BC(NBC))
    ! let's do the faces first
    FACE_AREA_INNER=CYLINDER%R1*CYLINDER%DT*CYLINDER%DZ       ! true geometry
    FACE_AREA_OUTER=CYLINDER%R2*CYLINDER%DT*CYLINDER%DZ       ! true geometry
!     FACE_AREA_INNER=2.0_DP*CYLINDER%R1*SIN(CYLINDER%DT/2.0_DP)*CYLINDER%DZ  ! linearised geom
!     FACE_AREA_OUTER=2.0_DP*CYLINDER%R2*SIN(CYLINDER%DT/2.0_DP)*CYLINDER%DZ  ! linearised geom
    FORCE_INNER=P_INNER*FACE_AREA_INNER
    FORCE_OUTER=P_OUTER*FACE_AREA_OUTER
    BC_IDX=0
    ! INNER FACE
    DO I=1,CYLINDER%NET
      ! do x component first
      BC_IDX=BC_IDX+1
      ALLOCATE(BC(BC_IDX)%NODES(NEZ+1))
      ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
      ALLOCATE(BC(BC_IDX)%DERIVATIVES(1)) ! TODO: TRILINEAR BASIS ONLY AT THE MO
      BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_NEUMANN
      J=(I-1)*(CYLINDER%NEZ+1)+1 ! starting index of nodes in current z-column
      BC(BC_IDX)%NODES = CYLINDER%INNER_FACE(J:J+NEZ)
      BC(BC_IDX)%COMPONENTS = 1  ! X component first
      BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
      ! calculate the X component force and assign
      FORCE_X = FORCE_INNER*COS(CYLINDER%DT*(I-1))
      BC(BC_IDX)%INCREMENT = FORCE_X
      ! do the y component now
      BC_IDX=BC_IDX+1
      ALLOCATE(BC(BC_IDX)%NODES(NEZ+1))
      ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
      ALLOCATE(BC(BC_IDX)%DERIVATIVES(1))
      BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_NEUMANN      
      BC(BC_IDX)%NODES = CYLINDER%INNER_FACE(J:J+NEZ)
      BC(BC_IDX)%COMPONENTS = 2  ! Y component now
      BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
      FORCE_Y = FORCE_INNER*SIN(CYLINDER%DT*(I-1))
      BC(BC_IDX)%INCREMENT = FORCE_Y
    ENDDO

    ! ----------------------------
    ! REPEAT ALL FOR OUTER FACE
    DO I=1,CYLINDER%NET
      ! do x component first
      BC_IDX=BC_IDX+1
      ALLOCATE(BC(BC_IDX)%NODES(NEZ+1))
      ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
      ALLOCATE(BC(BC_IDX)%DERIVATIVES(1)) ! TODO: TRILINEAR BASIS ONLY AT THE MO
      BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_NEUMANN
      J=(I-1)*(CYLINDER%NEZ+1)+1 ! starting index of nodes in current z-column
      BC(BC_IDX)%NODES = CYLINDER%OUTER_FACE(J:J+NEZ)
      BC(BC_IDX)%COMPONENTS = 1  ! X component first
      BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
      ! calculate the X component force and assign
      FORCE_X = FORCE_OUTER*COS(CYLINDER%DT*(I-1))
      BC(BC_IDX)%INCREMENT = FORCE_X
      ! do the y component now
      BC_IDX=BC_IDX+1
      ALLOCATE(BC(BC_IDX)%NODES(NEZ+1))
      ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
      ALLOCATE(BC(BC_IDX)%DERIVATIVES(1))
      BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_NEUMANN      
      BC(BC_IDX)%NODES = CYLINDER%OUTER_FACE(J:J+NEZ)
      BC(BC_IDX)%COMPONENTS = 2  ! Y component now
      BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
      FORCE_Y = FORCE_OUTER*SIN(CYLINDER%DT*(I-1))
      BC(BC_IDX)%INCREMENT = FORCE_Y
    ENDDO

    ! ------------------------
    ! Z: fix bottom nodes
    BC_IDX=BC_IDX+1
    ALLOCATE(BC(BC_IDX)%NODES(CYLINDER%NET*(CYLINDER%NER+1)))
    ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
    ALLOCATE(BC(BC_IDX)%DERIVATIVES(1)) ! TODO: TRILINEAR BASIS ONLY AT THE MO
    BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_DIRICHLET
    BC(BC_IDX)%NODES = CYLINDER%BOTTOM_FACE
    BC(BC_IDX)%COMPONENTS = 3  ! Z
    BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
    BC(BC_IDX)%INCREMENT = 0.0_DP

    ! Z: stretch upper nodes
    BC_IDX=BC_IDX+1
    ALLOCATE(BC(BC_IDX)%NODES(CYLINDER%NET*(CYLINDER%NER+1)))
    ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
    ALLOCATE(BC(BC_IDX)%DERIVATIVES(1)) ! TODO: TRILINEAR BASIS ONLY AT THE MO
    BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_DIRICHLET
    BC(BC_IDX)%NODES = CYLINDER%TOP_FACE
    BC(BC_IDX)%COMPONENTS = 3  ! Z
    BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
    BC(BC_IDX)%INCREMENT = CYLINDER%Z0*(LAMBDA-1.0_DP)

    ! ------------------------
    ! prevent rigid body motion
    ! fix Y component of node 1
    BC_IDX=BC_IDX+1
    ALLOCATE(BC(BC_IDX)%NODES(1))
    ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
    ALLOCATE(BC(BC_IDX)%DERIVATIVES(1))  ! TODO: TRILINEAR BASIS ONLY AT THE MO
    BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_DIRICHLET
    BC(BC_IDX)%NODES = CYLINDER%NODE_ON_X_AXIS  ! nodes 1 & 2 are fixed on x axis
    BC(BC_IDX)%COMPONENTS = 2                   ! Y
    BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
    BC(BC_IDX)%INCREMENT = 0.0_DP
    ! fix X component of node that's on y axis initially* (NEED A NODE THERE!!)
    BC_IDX=BC_IDX+1
    ALLOCATE(BC(BC_IDX)%NODES(1))
    ALLOCATE(BC(BC_IDX)%COMPONENTS(1))
    ALLOCATE(BC(BC_IDX)%DERIVATIVES(1))  ! TODO: TRILINEAR BASIS ONLY AT THE MO
    BC(BC_IDX)%BC_TYPE = SOLID_MECHANICS_IO_BC_DIRICHLET
    BC(BC_IDX)%NODES = CYLINDER%NODE_ON_Y_AXIS  ! nodes 1 & 2 are fixed on x axis
    BC(BC_IDX)%COMPONENTS = 1                   ! X
    BC(BC_IDX)%DERIVATIVES = 1 ! TODO: TRILINEAR ONLY AT THE MOMENT
    BC(BC_IDX)%INCREMENT = 0.0_DP


    CALL EXITS("SOLID_MECHANICS_IO_CREATE_BC")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_CREATE_BC",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_CREATE_BC")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_CREATE_BC


END MODULE SOLID_MECHANICS_IO_CYLINDER_GEOMETRY