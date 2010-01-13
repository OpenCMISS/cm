!> \file
!> $Id: finite_elasticity_routines.f90 644 2009-09-02 22:56:28Z kmith $
!> \author Chris Bradley
!> \brief This module handles all finite elasticity routines.
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

!>This module supplies routines to handle IO for solid mechanics problems
!
! JUST SOME COMMENTS:
! as this is a temporary stop-gap before FIELD_IO routines actually work
! there's a lot this doesn't handle:
! i.e. it's probably best to re-write from cm before running the files
! through these routines. They not super-general. but the are commented
! throughout the code, to identify where the problems are.
MODULE SOLID_MECHANICS_IO_ROUTINES
  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE DOMAIN_MAPPINGS
  USE FIELD_ROUTINES
  USE MESH_ROUTINES
  USE TYPES
  USE KINDS
  USE STRINGS

  IMPLICIT NONE

  ! BARE MINIMAL CUSTOM TYPES
  TYPE NODE
    INTEGER(INTG) :: ID
    REAL(DP),ALLOCATABLE :: X(:) ! YES, this is fixed for RC coordinates
    REAL(DP),ALLOCATABLE :: Y(:)
    REAL(DP),ALLOCATABLE :: Z(:)
  END TYPE

  TYPE ELEMENT
    INTEGER(INTG),ALLOCATABLE :: NODES(:) ! list of nodes in this element
  END TYPE

  TYPE SOLID_MECHANICS_IO_MESH_CONTAINER
    INTEGER(INTG) :: NN                 !>total number of nodes
    INTEGER(INTG) :: NE                 !>total number of elements
    INTEGER(INTG) :: NC                 !>number of coordinates
    INTEGER(INTG),ALLOCATABLE :: ND(:)  !>number of derivatives in each coordinate
    TYPE(NODE),ALLOCATABLE :: NODES(:)  !>subtype for nodal dofs
    TYPE(ELEMENT),ALLOCATABLE :: ELS(:) !>subtype for element topology
    LOGICAL :: FINISHED                 !>true if mesh import is completed
  END TYPE

  TYPE SOLID_MECHANICS_IO_BOUNDARY_CONDITION
    INTEGER(INTG),ALLOCATABLE :: NODES(:)
    INTEGER(INTG) :: BC_TYPE
    INTEGER(INTG),ALLOCATABLE :: COMPONENTS(:)
    INTEGER(INTG),ALLOCATABLE :: DERIVATIVES(:)
    REAL(DP) :: INCREMENT
  END TYPE

  ! MODULE CONSTANTS
  INTEGER(INTG),PARAMETER :: SOLID_MECHANICS_IO_BC_DIRICHLET=1
  INTEGER(INTG),PARAMETER :: SOLID_MECHANICS_IO_BC_NEUMANN=2

  ! MODULE PROCEDURES
  PUBLIC :: SOLID_MECHANICS_IO_READ_NODES, SOLID_MECHANICS_IO_READ_ELEMENTS
  PUBLIC :: SOLID_MECHANICS_IO_WRITE_MESH, SOLID_MECHANICS_IO_READ_MESH
  PUBLIC :: SOLID_MECHANICS_IO_ASSIGN_NODES, SOLID_MECHANICS_IO_ASSIGN_ELEMENTS
  PUBLIC :: SOLID_MECHANICS_IO_READ_BC, SOLID_MECHANICS_IO_ASSIGN_BC
  PUBLIC :: SOLID_MECHANICS_IO_CLEAR_MESH, SOLID_MECHANICS_IO_CLEAR_BC
!   PUBLIC :: SOLID_MECHANICS_IO_WRITE_CMGUI_NODES
!   PUBLIC :: SOLID_MECHANICS_IO_WRITE_CMGUI
  PRIVATE :: COUNT_COMMAS, PRINT_FIELD_INFO!, CREATE_FACE_INFO

  CONTAINS

  ! =============================================================================
  !>main calling routine to import legacy cm ip-files
  SUBROUTINE SOLID_MECHANICS_IO_READ_IPFILES(NAME,MESH,ERR,ERROR,*)
    !Argument variables
    TYPE(VARYING_STRING),INTENT(IN) :: NAME                     !<path of ip files without extension
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(OUT) :: MESH !<on exit, contains the imported mesh
    INTEGER(INTG),INTENT(OUT) :: ERR                            !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                   !<The error string

    CALL ENTERS("SOLID_MECHANICS_IO_READ_IPFILES",ERR,ERROR,*999)

    ! IT'S ASSUMED THAT .ipnode AND .ipelem HAVE THE SAME NAME FOR NOW
    MESH%FINISHED = .FALSE.
    CALL SOLID_MECHANICS_IO_READ_NODES(NAME//".ipnode",MESH,ERR,ERROR,*999)
    CALL SOLID_MECHANICS_IO_READ_ELEMENTS(NAME//".ipelem",MESH,ERR,ERROR,*999)
    MESH%FINISHED = .TRUE.

    CALL EXITS("SOLID_MECHANICS_IO_READ_IPFILES")

    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_READ_NODES",ERR,ERROR)
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_READ_IPFILES

  ! =============================================================================
  !>reads in a legacy ipnode file
  !>imported node information is stored into an intermediate mesh container type
  !>defined in this module
  SUBROUTINE SOLID_MECHANICS_IO_READ_NODES(FILENAME,MESH,ERR,ERROR,*)
    TYPE(VARYING_STRING),INTENT(IN) :: FILENAME
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(INOUT) :: MESH
    INTEGER(INTG),INTENT(OUT) :: ERR                            !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                   !<The error string
    ! LOCAL VARIABLES
    INTEGER(INTG) :: I,J,N,NNI,LENG,ENDL,NDI  ! GENERALLY USED AS DUMMIES
    CHARACTER*128 :: LINE
    CHARACTER*30 :: VAL
    LOGICAL :: OPEN

    CALL ENTERS("SOLID_MECHANICS_IO_READ_NODES",ERR,ERROR,*999)
    ! OPEN FILE
    OPEN(UNIT=5,FILE=CHAR(FILENAME),STATUS='OLD',IOSTAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("error opening file "//CHAR(FILENAME),ERR,ERROR,*999)

    ! READ THE HEADER INFO
    NDI=0
    DO
      READ(5,'(128A)') LINE
      LENG=LEN_TRIM(LINE)
      ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
      VAL=LINE(ENDL+1:LENG)
      IF(INDEX(LINE,'CMISS Version')>0) THEN
        ! FIRST CHECK IF IT'S THE RIGHT FILE
        IF(INDEX(LINE,'ipnode')==0) CALL FLAG_ERROR("ipnode file is invalid",ERR,ERROR,*999)
      ELSEIF(INDEX(LINE,'The number of nodes is')>0) THEN
        MESH%NN = STRING_TO_INTEGER(VAL,ERR,ERROR)
      ELSEIF(INDEX(LINE,'Number of coordinates')>0) THEN
        MESH%NC = STRING_TO_INTEGER(VAL,ERR,ERROR)
        IF(MESH%NC/=3) CALL FLAG_ERROR("currently only 3-dimensional meshes are handled.",ERR,ERROR,*999)
        ALLOCATE(MESH%ND(MESH%NC))
      ELSEIF(INDEX(LINE,'Do you want prompting for different versions of nj')>0) THEN
        ! VERSIONS NOT IMPLEMENTED AT THE MOMENT
      ELSEIF(INDEX(LINE,'The number of derivatives for coordinate')>0) THEN
        NDI=NDI+1
        MESH%ND(NDI)=STRING_TO_INTEGER(VAL,ERR,ERROR)+1
        IF(NDI==MESH%NC) EXIT
      ELSEIF(INDEX(LINE,'Node number [')>0) THEN
        CALL FLAG_ERROR("header block is incomplete",ERR,ERROR,*999)
      ENDIF
    ENDDO

    ! ALLOCATE NODES CONTAINER
    ALLOCATE(MESH%NODES(MESH%NN))
    DO I=1,MESH%NN
      ALLOCATE(MESH%NODES(I)%X(MESH%ND(1)))
      ALLOCATE(MESH%NODES(I)%Y(MESH%ND(2)))
      ALLOCATE(MESH%NODES(I)%Z(MESH%ND(3)))
    ENDDO

    ! READ THE NODE INFO
    NNI=0
    DO
      READ(5,'(128A)',END=50) LINE
      IF(INDEX(LINE,'Node number [')>0) THEN
        LENG=LEN_TRIM(LINE)
        ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
        N=STRING_TO_INTEGER(LINE(ENDL+1:LENG),ERR,ERROR) ! CURRENT NODE NUMBER 
        MESH%NODES(NNI+1)%ID=N                           ! CAN BE OUT OF SEQUENCE/NOT START FROM 1/ETC
        ! NB. A 3 DIMENSIONAL MESH IS ASSUMED HERE
        DO J=1,MESH%ND(1)
          READ(5,'(128A)') LINE
          LENG=LEN_TRIM(LINE)
          ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
          MESH%NODES(NNI+1)%X(J)=STRING_TO_DOUBLE(LINE(ENDL+1:LENG),ERR,ERROR)
        ENDDO
        DO J=1,MESH%ND(2)
          READ(5,'(128A)') LINE
          LENG=LEN_TRIM(LINE)
          ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
          MESH%NODES(NNI+1)%Y(J)=STRING_TO_DOUBLE(LINE(ENDL+1:LENG),ERR,ERROR)
        ENDDO
        DO J=1,MESH%ND(3)
          READ(5,'(128A)') LINE
          LENG=LEN_TRIM(LINE)
          ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
          MESH%NODES(NNI+1)%Z(J)=STRING_TO_DOUBLE(LINE(ENDL+1:LENG),ERR,ERROR)
        ENDDO
        NNI=NNI+1   ! IT'S DONE. SIGN IT OFF
      ENDIF
    ENDDO

50 CLOSE(5)

    ! SIMPLE CHECK TO SEE IF IMPORT WAS COMPLETE (NN==N)
    IF(MESH%NN/=NNI) CALL FLAG_ERROR("nodes import failed before completion",ERR,ERROR,*999)
    CLOSE(5)
    CALL EXITS("SOLID_MECHANICS_IO_READ_NODES")
    RETURN

999 CALL ERRORS("SOLID_MECHANICS_IO_READ_NODES",ERR,ERROR)
    INQUIRE(5,OPENED=OPEN)
    IF(OPEN) CLOSE(5)
    CALL EXITS("SOLID_MECHANICS_IO_READ_NODES")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_READ_NODES

  ! =============================================================================
  SUBROUTINE SOLID_MECHANICS_IO_READ_ELEMENTS(FILENAME,MESH,ERR,ERROR,*)
    TYPE(VARYING_STRING),INTENT(IN) :: FILENAME
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(INOUT) :: MESH
    INTEGER(INTG),INTENT(OUT) :: ERR                            !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                   !<The error string
    ! local variables
    INTEGER(INTG) :: I,E,NEI,LENG,ENDL  ! USED AS DUMMIES
    CHARACTER*128 :: LINE
    CHARACTER*30 :: VAL
    LOGICAL :: OPEN

    CALL ENTERS("SOLID_MECHANICS_IO_READ_ELEMENTS",ERR,ERROR,*999)
    ! OPEN FILE
    OPEN(UNIT=5,FILE=CHAR(FILENAME),STATUS='OLD',IOSTAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("error opening file "//CHAR(FILENAME),ERR,ERROR,*999)

    ! READ THE HEADER INFO
    DO
      READ(5,'(128A)') LINE
      LENG=LEN_TRIM(LINE)
      ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
      VAL=LINE(ENDL+1:LENG)
      IF(INDEX(LINE,'CMISS Version')>0) THEN
        ! FIRST CHECK IF IT'S THE RIGHT FILE
        IF(INDEX(LINE,'ipelem')==0) CALL FLAG_ERROR("ipelem file is invalid",ERR,ERROR,*999)
      ELSEIF(INDEX(LINE,'The number of elements is')>0) THEN
        MESH%NE = STRING_TO_INTEGER(VAL,ERR,ERROR)
        EXIT
      ENDIF
    ENDDO

    ! ALLOCATE THE ELEMENTS
    ! TODO: ONLY LINEAR OR CUBIC HERMITE FOR NOW (8 NODES PER EL)
    IF(.NOT.ALLOCATED(MESH%ELS)) ALLOCATE(MESH%ELS(MESH%NE))
    DO I=1,MESH%NE
      ALLOCATE(MESH%ELS(I)%NODES(8))
    ENDDO

    ! READ IN THE ELEMENTS
    NEI=0
    DO
      READ(5,'(128A)',END=50) LINE
      IF(INDEX(LINE,'Element number [')>0) THEN
        LENG=LEN_TRIM(LINE)
        ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
        E=STRING_TO_INTEGER(LINE(ENDL+1:LENG),ERR,ERROR) ! CURRENT EL NUMBER (CAN BE OUT OF SEQUENCE)
        ! NB. A 3 DIMENSIONAL MESH IS ASSUMED HERE
        ! TODO: typically 2 components are defined one each for U and p, we'll assume they're same for now
        DO
          READ(5,'(128A)',END=50) LINE
          IF(INDEX(LINE,'Enter the 8 global numbers for basis')>0) THEN
            ! ASSUMING THE 8 NUMBERS ARE IN SAME LINE WITH 1+ WHITESPACES IN BETWEEN
            DO I=8,1,-1
              LENG=LEN_TRIM(LINE)
              ENDL=INDEX(TRIM(LINE),' ',.TRUE.)
              MESH%ELS(E)%NODES(I)=STRING_TO_INTEGER(LINE(ENDL+1:LENG),ERR,ERROR)
              ! TRIM AWAY THE LAST NUMBER
              LINE=LINE(1:ENDL)
            ENDDO
            EXIT        ! GO BACK TO READING 'element number ['
          ENDIF
        ENDDO
        NEI=NEI+1   ! IT'S DONE. SIGN IT OFF
      ENDIF
    ENDDO
    
50 CLOSE(5)

    ! SIMPLE CHECK TO SEE IF IMPORT WAS COMPLETE (NN==N)
    IF(MESH%NE/=NEI) CALL FLAG_ERROR("element import failed before completion",ERR,ERROR,*999)
    CLOSE(5)
    CALL EXITS("SOLID_MECHANICS_IO_READ_ELEMENTS")
    RETURN

999 CALL ERRORS("SOLID_MECHANICS_IO_READ_ELEMENTS",ERR,ERROR)
    INQUIRE(5,OPENED=OPEN)
    IF(OPEN) CLOSE(5)
    CALL EXITS("SOLID_MECHANICS_IO_READ_ELEMENTS")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_READ_ELEMENTS
  
  ! =============================================================================
  !>prints the mesh to the given output stream (file or screen)
  SUBROUTINE SOLID_MECHANICS_IO_WRITE_MESH(UNIT,MESH,*)
    INTEGER(INTG),INTENT(IN) :: UNIT                    !<unit 6 is reserved by fortran for screen output
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(IN) :: MESH      !<mesh to be printed out
    ! local variables
    INTEGER(INTG) :: NC_IDX,NN_IDX,NE_IDX
    INTEGER(INTG) :: ERR
    TYPE(VARYING_STRING):: ERROR

    CALL ENTERS("SOLID_MECHANICS_IO_WRITE_MESH",ERR,ERROR,*999)

    IF(.NOT.MESH%FINISHED) THEN
      WRITE(*,*) "SOLID_MECHANICS_IO_WRITE_MESH: MESH IMPORT IS NOT COMPLETE"   ! KEEP GOING?
    ENDIF

    ! HEADER INFO
    WRITE(UNIT,*) "MESH%NN=",MESH%NN
    WRITE(UNIT,*) "MESH%NE=",MESH%NE
    WRITE(UNIT,*) "MESH%NC=",MESH%NC
    DO NC_IDX=1,MESH%NC
      WRITE(UNIT,'(" MESH%ND(",I1,")=",I1)') NC_IDX,MESH%ND(NC_IDX)
    ENDDO

    ! NODES
    DO NN_IDX=1,MESH%NN
      WRITE(UNIT,*) "NODE ",NN_IDX
      WRITE(UNIT,'(100g25.17:," ")') MESH%NODES(NN_IDX)%X
      WRITE(UNIT,'(100g25.17:," ")') MESH%NODES(NN_IDX)%Y
      WRITE(UNIT,'(100g25.17:," ")') MESH%NODES(NN_IDX)%Z
    ENDDO

    ! ELEMENTS
    DO NE_IDX=1,MESH%NE
      WRITE(UNIT,*) "ELEMENT ",NE_IDX
      WRITE(UNIT,'(100I10:," ")') MESH%ELS(NE_IDX)%NODES
    ENDDO

    CALL EXITS("SOLID_MECHANICS_IO_WRITE_MESH")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_WRITE_MESH",ERR,ERROR)
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_WRITE_MESH

  ! =============================================================================
  !>reads back in the text mesh format
  SUBROUTINE SOLID_MECHANICS_IO_READ_MESH(FILENAME,MESH,ERR,ERROR,*)
    TYPE(VARYING_STRING),INTENT(IN) :: FILENAME                 !<full path of the file to be read in
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(OUT) :: MESH !<imported mesh
    INTEGER(INTG), INTENT(OUT) :: ERR                           !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR                   !<The error string
    ! local variables
    INTEGER(INTG) :: IOS,NC_IDX,NN_IDX,NE_IDX
    CHARACTER*128 :: WORD
    
    CALL ENTERS("SOLID_MECHANICS_IO_READ_MESH",ERR,ERROR,*999)

    ! OPEN AND CHECK THE FILE
    OPEN(UNIT=55,FILE=CHAR(FILENAME),STATUS='OLD',IOSTAT=IOS)
    IF(IOS/=0) CALL FLAG_ERROR("error opening file "//CHAR(FILENAME),ERR,ERROR,*999)

    ! READ HEADER INFO
    READ(55,'(A)',ERR=999) WORD
    MESH%NN=STRING_TO_INTEGER(WORD(INDEX(WORD,'=')+1:LEN(WORD)),ERR,ERROR)
    READ(55,'(A)',ERR=999) WORD
    MESH%NE=STRING_TO_INTEGER(WORD(INDEX(WORD,'=')+1:LEN(WORD)),ERR,ERROR)
    READ(55,'(A)',ERR=999) WORD
    MESH%NC=STRING_TO_INTEGER(WORD(INDEX(WORD,'=')+1:LEN(WORD)),ERR,ERROR)
    ALLOCATE(MESH%ND(MESH%NC))
    DO NC_IDX=1,MESH%NC
      READ(55,'(A)',ERR=999) WORD
      MESH%ND(NC_IDX)=STRING_TO_INTEGER(WORD(INDEX(WORD,'=')+1:LEN(WORD)),ERR,ERROR)
    ENDDO

    ! NODES
    ALLOCATE(MESH%NODES(MESH%NN))
    DO NN_IDX=1,MESH%NN
      ALLOCATE(MESH%NODES(NN_IDX)%X(MESH%ND(1)))
      ALLOCATE(MESH%NODES(NN_IDX)%Y(MESH%ND(2)))
      ALLOCATE(MESH%NODES(NN_IDX)%Z(MESH%ND(3)))
      READ(55,'(128A)') WORD
      READ(55,*,ERR=999) MESH%NODES(NN_IDX)%X
      READ(55,*,ERR=999) MESH%NODES(NN_IDX)%Y
      READ(55,*,ERR=999) MESH%NODES(NN_IDX)%Z
    ENDDO

    ! ELEMENTS
    ALLOCATE(MESH%ELS(MESH%NE))
    DO NE_IDX=1,MESH%NE
      ALLOCATE(MESH%ELS(NE_IDX)%NODES(8))
      READ(55,'(128A)') WORD
      READ(55,*,ERR=999) MESH%ELS(NE_IDX)%NODES
    ENDDO

    CALL EXITS("SOLID_MECHANICS_IO_READ_MESH")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_READ_MESH",ERR,ERROR)
    RETURN
  END SUBROUTINE SOLID_MECHANICS_IO_READ_MESH

  ! =============================================================================
  !>assigns the imported nodes to OpenCMISS-native mesh data structure
  SUBROUTINE SOLID_MECHANICS_IO_ASSIGN_NODES(MESH,FIELD,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(IN) :: MESH      !<imported mesh
    TYPE(FIELD_TYPE),POINTER :: FIELD                   !<field to which the mesh is to be assigned
    INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
    TYPE(VARYING_STRING),INTENT(OUT) :: ERROR           !<The error string
    ! local variables
    INTEGER(INTG) :: DERIVATIVE_IDX,NODE_IDX
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT
    TYPE(DECOMPOSITION_TYPE),POINTER :: decomposition
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    LOGICAL :: LOCAL_EXISTS

    CALL ENTERS("SOLID_MECHANICS_IO_ASSIGN_NODES",ERR,ERROR,*999)

    ! IMPORTED MESH SHOULD BE COMPLETE
    IF(.NOT.MESH%FINISHED) CALL FLAG_ERROR("mesh definition is not finished",ERR,ERROR,*999)
    
    IF(ASSOCIATED(field%decomposition)) THEN
      decomposition=>field%decomposition
      ! ah fuck it let's skip the tests
      DO NODE_IDX=1,MESH%NN
        FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(1)%PTR
        MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
        DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
        NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
        CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,NODE_IDX,LOCAL_EXISTS, &
        & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
        IF(LOCAL_EXISTS) THEN
          DO DERIVATIVE_IDX=1,MESH%ND(1) ! FIRST (X) COMPONENT
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & DERIVATIVE_IDX,NODE_IDX,1,MESH%NODES(NODE_IDX)%X(DERIVATIVE_IDX),ERR,ERROR,*999)
          ENDDO
          DO DERIVATIVE_IDX=1,MESH%ND(2) ! SECOND (Y) COMPONENT
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & DERIVATIVE_IDX,NODE_IDX,2,MESH%NODES(NODE_IDX)%Y(DERIVATIVE_IDX),ERR,ERROR,*999)
          ENDDO
          DO DERIVATIVE_IDX=1,MESH%ND(3) ! THIRD (Z) COMPONENT
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & DERIVATIVE_IDX,NODE_IDX,3,MESH%NODES(NODE_IDX)%Z(DERIVATIVE_IDX),ERR,ERROR,*999)
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    CALL EXITS("SOLID_MECHANICS_IO_ASSIGN_NODES")
    RETURN

999 CALL ERRORS("SOLID_MECHANICS_IO_ASSIGN_NODES",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_ASSIGN_NODES")
    RETURN
  END SUBROUTINE SOLID_MECHANICS_IO_ASSIGN_NODES

  ! =============================================================================
  !>assigns the imported elements to OpenCMISS-native mesh data structure
  SUBROUTINE SOLID_MECHANICS_IO_ASSIGN_ELEMENTS(MESH,ELEMENTS,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(IN) :: MESH      !<imported mesh
    TYPE(MESH_ELEMENTS_TYPE),pointer :: ELEMENTS        !<element type to which the mesh is to be assigned
    INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
    ! local variables
    INTEGER(INTG) :: ELEMENT_IDX

    CALL ENTERS("SOLID_MECHANICS_IO_ASSIGN_ELEMENTS",ERR,ERROR,*999)

    DO ELEMENT_IDX=1,MESH%NE
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ELEMENT_IDX,ELEMENTS, &
      & MESH%ELS(ELEMENT_IDX)%NODES,ERR,ERROR,*999)
    ENDDO

    CALL EXITS("SOLID_MECHANICS_IO_ASSIGN_ELEMENTS")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_ASSIGN_ELEMENTS",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_ASSIGN_NODES")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_ASSIGN_ELEMENTS

  ! =============================================================================
  !>clears all data contained in the mesh container and deallocates memory
  SUBROUTINE SOLID_MECHANICS_IO_CLEAR_MESH(MESH,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER),INTENT(INOUT) :: MESH !<imported mesh
    INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
    ! local variables
    INTEGER(INTG) :: NN_IDX,NE_IDX
    
    CALL ENTERS("SOLID_MECHANICS_IO_ASSIGN_ELEMENTS",ERR,ERROR,*999)

    IF(ALLOCATED(MESH%ND)) DEALLOCATE(MESH%ND)
    IF(ALLOCATED(MESH%NODES)) THEN
      DO NN_IDX=1,MESH%NN
        IF(ALLOCATED(MESH%NODES(NN_IDX)%X)) DEALLOCATE(MESH%NODES(NN_IDX)%X)
        IF(ALLOCATED(MESH%NODES(NN_IDX)%Y)) DEALLOCATE(MESH%NODES(NN_IDX)%Y)
        IF(ALLOCATED(MESH%NODES(NN_IDX)%Z)) DEALLOCATE(MESH%NODES(NN_IDX)%Z)
      ENDDO
      DEALLOCATE(MESH%NODES)
    ENDIF
    IF(ALLOCATED(MESH%ELS)) THEN
      DO NE_IDX=1,MESH%NE
        IF(ALLOCATED(MESH%ELS(NE_IDX)%NODES)) DEALLOCATE(MESH%ELS(NE_IDX)%NODES)
      ENDDO
      DEALLOCATE(MESH%ELS)
    ENDIF
    MESH%NN=0
    MESH%NE=0
    MESH%NC=0
    MESH%FINISHED=.FALSE.

    CALL EXITS("SOLID_MECHANICS_IO_CLEAR_MESH")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_CLEAR_MESH",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_CLEAR_MESH")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_CLEAR_MESH

  ! =============================================================================
  !>imports boundary conditions from text file format
  !>TODO: only implemented for nodal bc currently
  SUBROUTINE SOLID_MECHANICS_IO_READ_BC(FILENAME,BC_OUT,ERR,ERROR,*)
    TYPE(VARYING_STRING),INTENT(IN) :: FILENAME
    TYPE(SOLID_MECHANICS_IO_BOUNDARY_CONDITION),ALLOCATABLE,INTENT(OUT) :: BC_OUT(:)
    INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
    ! local variables
    TYPE(SOLID_MECHANICS_IO_BOUNDARY_CONDITION) :: BC(1000) ! HARDCODING ARRAY IS A HACK
    CHARACTER*132 :: LINE,LINE2                             ! BUT SO IS THIS ENTIRE MODULE
    INTEGER(INTG) :: NODES(1000)
    INTEGER(INTG) :: COMPONENTS(3)
    INTEGER(INTG) :: DERIVATIVES(8)
    INTEGER(INTG) :: IDX,NBC,NN,NC,ND     ! COUNTERS
    INTEGER(INTG) :: LENG,POS             ! LINE POSITIONS
    LOGICAL :: READING  ! TRUE IF A BC IS BEING READ
    LOGICAL :: OPEN

    ! OPEN THE FILE
    CALL ENTERS("SOLID_MECHANICS_IO_READ_BC",ERR,ERROR,*999)
    ! OPEN FILE
    OPEN(UNIT=9,FILE=CHAR(FILENAME),STATUS='OLD',IOSTAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("error opening file "//CHAR(FILENAME),ERR,ERROR,*999)

    NBC=0
    READING=.FALSE.
    DO
      READ(9,'(A)',END=50) LINE
      IF(LEN_TRIM(LINE)==0) CYCLE       ! BLANK LINE
      ! SIGH. READ ENTIRE LINE AND DISSECT IT. HATE.
      ! TODO: IF THE NUMBERS ARE SPREAD OVER >1 LINE, WE'LL HAVE TROUBLE
      LINE=ADJUSTL(LINE)

      IF(INDEX(LINE,'!')==1 .OR. INDEX(LINE,'%')==1 .OR. INDEX(LINE,'/')==1 &
         & .OR. INDEX(LINE,'#')==1) CYCLE ! THIS IS A COMMENT. DO NOTHING.

      LINE=CHARACTER_TO_LOWERCASE(LINE)
      LENG=LEN_TRIM(LINE)
      POS=INDEX(TRIM(LINE),' ')
      LINE2=LINE(POS+1:LENG)    ! WITHOUT THE PRECEDING LABEL
      IF(INDEX(LINE,'nodes:')==1) THEN
        ! ARE WE IN THE MIDDLE OF READING ANOTHER BC?
        IF(READING) CALL FLAG_ERROR("Check error in boundary condition file",ERR,ERROR,*999)
        READING=.TRUE.  ! NO. START ONE NOW THEN.
        NBC=NBC+1
        NODES=0 ! RESET BEFORE READING
        ! COUNT THE NUMBER OF COMMAS TO FIGURE OUT HOW MANY ENTRIES THERE ARE
        NN=COUNT_COMMAS(LINE2,LENG-POS)+1
        READ(LINE2,*,ERR=999) NODES(1:NN)
        ALLOCATE(BC(NBC)%NODES(NN))
        BC(NBC)%NODES=NODES(1:NN)
      ELSEIF(INDEX(LINE,'type:')==1) THEN
        IF(INDEX(LINE2,'dirichlet')==1) THEN
          BC(NBC)%BC_TYPE = SOLID_MECHANICS_IO_BC_DIRICHLET
        ELSEIF(INDEX(LINE2,'neumann')==1) THEN
          BC(NBC)%BC_TYPE = SOLID_MECHANICS_IO_BC_NEUMANN
        ELSE
          CALL FLAG_ERROR("Unknown boundary condition type",ERR,ERROR,*999)
        ENDIF
      ELSEIF(INDEX(LINE,'components:')==1) THEN
        COMPONENTS=0 ! RESET BEFORE READING
        NC=COUNT_COMMAS(LINE2,LENG-POS)+1
        READ(LINE2,*,ERR=999) COMPONENTS(1:NC)
        ALLOCATE(BC(NBC)%COMPONENTS(NC))
        BC(NBC)%COMPONENTS = COMPONENTS(1:NC)
      ELSEIF(INDEX(LINE,'derivatives:')==1) THEN
        DERIVATIVES=0 ! RESET BEFORE READING
        ND=COUNT_COMMAS(LINE2,LENG-POS)+1
        READ(LINE2,*,ERR=999) DERIVATIVES(1:ND)
        ALLOCATE(BC(NBC)%DERIVATIVES(ND))
        BC(NBC)%DERIVATIVES = DERIVATIVES(1:ND)
      ELSEIF(INDEX(LINE,'increment:')==1) THEN
        READ(LINE2,*,ERR=999,END=50) BC(NBC)%INCREMENT
        READING=.FALSE.
      ENDIF
    ENDDO

50 CLOSE(9)

    ! LET'S COPY BC TO SIZE
    ALLOCATE(BC_OUT(NBC))
    DO IDX=1,NBC
      ALLOCATE(BC_OUT(IDX)%NODES(SIZE(BC(IDX)%NODES)))
      ALLOCATE(BC_OUT(IDX)%COMPONENTS(SIZE(BC(IDX)%COMPONENTS)))
      ALLOCATE(BC_OUT(IDX)%DERIVATIVES(SIZE(BC(IDX)%DERIVATIVES)))
      BC_OUT(IDX)%NODES=BC(IDX)%NODES
      BC_OUT(IDX)%BC_TYPE=BC(IDX)%BC_TYPE
      BC_OUT(IDX)%COMPONENTS=BC(IDX)%COMPONENTS
      BC_OUT(IDX)%DERIVATIVES=BC(IDX)%DERIVATIVES
      BC_OUT(IDX)%INCREMENT=BC(IDX)%INCREMENT
    ENDDO

    CALL EXITS("SOLID_MECHANICS_IO_READ_BC")
    RETURN

999 CALL ERRORS("SOLID_MECHANICS_IO_READ_BC",ERR,ERROR)
    INQUIRE(9,OPENED=OPEN)
    IF(OPEN) CLOSE(9)
    CALL EXITS("SOLID_MECHANICS_IO_READ_BC")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_READ_BC

  ! =============================================================================
  !>assigns the imported boundary conditions into OpenCMISS data structure
  SUBROUTINE SOLID_MECHANICS_IO_ASSIGN_BC(BC,BCP,MESH,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_BOUNDARY_CONDITION),INTENT(IN) :: BC(:)      !<imported boundary conditions
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BCP      !<OpenCMISS boundary conditions pointer
    TYPE(SOLID_MECHANICS_IO_MESH_CONTAINER) :: MESH     !<imported mesh container
    INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
    ! local variables
    INTEGER(INTG) :: NBC,NN,NC,ND               ! COUNTERS
    INTEGER(INTG) :: BC_IDX,N_IDX,C_IDX,D_IDX   ! LOOP INDICES
    INTEGER(INTG) :: N,C,D                      ! INTERMEDIATE HOLDERS
    REAL(DP) :: VAL                             ! VALUE OF BC
    INTEGER(INTG) :: BC_TYPE
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(FIELD_TYPE),POINTER :: FIELD

    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT
    TYPE(DECOMPOSITION_TYPE),POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    LOGICAL :: LOCAL_EXISTS

    CALL ENTERS("SOLID_MECHANICS_IO_ASSIGN_BC",ERR,ERROR,*999)

    ! basically, it is easier to prescribe increments rather than the actual
    ! values.. for now the increment is added to the initial conditions as listed
    ! in the imported mesh object (which may not be the best)
    NBC=SIZE(BC)        ! NUMBER OF BOUNDARY CONDITION GROUPS
    FIELD=>BCP%equations_set%dependent%dependent_field
! write(*,*) "   N    C    D        val"
    DO BC_IDX=1,NBC
      NN=SIZE(BC(BC_IDX)%NODES)
      DO N_IDX=1,NN
        N=BC(BC_IDX)%NODES(N_IDX)
        NC=SIZE(BC(BC_IDX)%COMPONENTS)
        DO C_IDX=1,NC
          C=BC(BC_IDX)%COMPONENTS(C_IDX)
          ND=SIZE(BC(BC_IDX)%DERIVATIVES)
          DO D_IDX=1,ND
            D=BC(BC_IDX)%DERIVATIVES(D_IDX)
            ! DON'T FORGET THAT INCREMENTS, RATHER THAN THE VALUES ARE GIVEN
            VAL=BC(BC_IDX)%INCREMENT    ! INCREMENT
            ! GRAB INITIAL CONDITION
            IF(BC(BC_IDX)%BC_TYPE==SOLID_MECHANICS_IO_BC_DIRICHLET) THEN
              ! ADD VALUES TO INCREMENT
              IF(C==1) THEN ! X,Y,Z HOW ANNOYING
                VAL=VAL+MESH%NODES(N)%X(D)
              ELSEIF(C==2) THEN
                VAL=VAL+MESH%NODES(N)%Y(D)
              ELSEIF(C==3) THEN
                VAL=VAL+MESH%NODES(N)%Z(D)
              ENDIF
              BC_TYPE=FIELD_U_VARIABLE_TYPE
            ELSEIF(BC(BC_IDX)%BC_TYPE==SOLID_MECHANICS_IO_BC_NEUMANN) THEN
              ! WHAT IS THE INITIAL VALUE FOR THE NEUMANN CONDITION?
              BC_TYPE=FIELD_DELUDELN_VARIABLE_TYPE
              ! TODO: CHECK IF THIS IS CORRECT ??
            ENDIF
            ! test to make it parallel-safe: set only if dof is non-ghost in current domain
DECOMPOSITION=>FIELD%DECOMPOSITION
FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(BC_TYPE)%PTR
MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,N,LOCAL_EXISTS,DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)

if(local_exists) then
CALL FIELD_COMPONENT_DOF_GET_USER_NODE(FIELD,BC_TYPE,D,N,C,local_ny,global_ny,ERR,ERROR,*999)
IF(FIELD%VARIABLE_type_map(BC_TYPE)%ptr%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
            ! ASSIGN
            CALL BOUNDARY_CONDITIONS_SET_NODE(BCP,BC_TYPE,D,N,C, &
              & BOUNDARY_CONDITION_FIXED,VAL,ERR,ERROR,*999)
            ! TODO: ONLY BOUNDARY_CONDITION_FIXED ALLOWED FOR NOW
ENDIF
endif
              ENDDO !D_IDX
            ENDDO !C_IDX
      ENDDO !N_IDX
    ENDDO !BC_IDX

    CALL EXITS("SOLID_MECHANICS_IO_ASSIGN_BC")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_ASSIGN_BC",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_ASSIGN_BC")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_ASSIGN_BC

  ! =============================================================================
  !>clears all data contained in the mesh container and deallocates memory
  SUBROUTINE SOLID_MECHANICS_IO_CLEAR_BC(BC,ERR,ERROR,*)
    TYPE(SOLID_MECHANICS_IO_BOUNDARY_CONDITION),ALLOCATABLE,INTENT(INOUT) :: BC(:)      !<BC to be cleared
    INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
    ! local variables
    INTEGER(INTG) :: BC_IDX

    CALL ENTERS("SOLID_MECHANICS_IO_CLEAR_BC",ERR,ERROR,*999)

    IF(ALLOCATED(BC)) THEN
      DO BC_IDX=1,SIZE(BC)
        IF(ALLOCATED(BC(BC_IDX)%NODES)) DEALLOCATE(BC(BC_IDX)%NODES)
        IF(ALLOCATED(BC(BC_IDX)%COMPONENTS)) DEALLOCATE(BC(BC_IDX)%COMPONENTS)
        IF(ALLOCATED(BC(BC_IDX)%DERIVATIVES)) DEALLOCATE(BC(BC_IDX)%DERIVATIVES)
      ENDDO
      DEALLOCATE(BC)
    ENDIF

    CALL EXITS("SOLID_MECHANICS_IO_CLEAR_BC")
    RETURN
999 CALL ERRORS("SOLID_MECHANICS_IO_CLEAR_BC",ERR,ERROR)
    CALL EXITS("SOLID_MECHANICS_IO_CLEAR_BC")
    RETURN

  END SUBROUTINE SOLID_MECHANICS_IO_CLEAR_BC

  ! =============================================================================
  !>counts the number of occurrences of commas in a string
  FUNCTION COUNT_COMMAS(LINE,LAST)
    CHARACTER(len=*),INTENT(INOUT) :: LINE
    INTEGER(INTG),OPTIONAL :: LAST
    INTEGER(INTG) :: COUNT_COMMAS
    ! local variables
    INTEGER(INTG) :: IDX

    IF(.NOT.PRESENT(LAST)) LAST=LEN(LINE)

    COUNT_COMMAS=0
    DO IDX=1,LAST
      IF(LINE(IDX:IDX).eq.',') COUNT_COMMAS=COUNT_COMMAS+1
    ENDDO

    RETURN
  END FUNCTION COUNT_COMMAS

!   ! =============================================================================
!   !> export existing nodes to a cmgui file
!   SUBROUTINE SOLID_MECHANICS_IO_WRITE_CMGUI_NODES(FILENAME,REGION,ERR,ERROR,*)
!     TYPE(VARYING_STRING), INTENT(IN) :: FILENAME !<file path without extension
!     TYPE(REGION_TYPE), POINTER :: REGION !<pointer to the region to get the coordinate system for
!     INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
!     ! local variables
!     CHARACTER*2 :: WORD
!     INTEGER(INTG) :: VALUE_IDX,I
!     CHARACTER*1 :: XYZ(3) = (/'x','y','z'/)
!     ! grab below from Region
!     INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_FIBRE_COMPONENTS
!     INTEGER(INTG) :: NUMBER_OF_MESH_COMPONENTS
!     INTEGER(INTG),ALLOCATABLE :: NODES_PER_MESH_COMPONENT(:)
!     TYPE(FIELD_PARAMETER_SET_TYPE),POINTER :: COORDINATES,DEFORMED,FIBRES
!     
!     ! TODO: PRESSURE FIELD IS KIND OF DROPPED AT THE MOMENT!? RECTIFY!
!     ! TODO: cubic Hermite IS NOT HANDLED CURRENTLY
! 
!     CALL ENTERS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS",ERR,ERROR,*999)
! 
! 22 FORMAT(I2) ! for use with WORD and WORD2
! 
!     ! scan REGION for parameters
!     ! TODO: CURRENTLY ASSUMES 1 MESH COMPONENT...!
!     NUMBER_OF_DIMENSIONS = REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
!     NUMBER_OF_FIBRE_COMPONENTS = REGION%EQUATIONS_SETS%EQUATIONS_SETS(1)%PTR%GEOMETRY% &
!       & FIBRE_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
!     NUMBER_OF_MESH_COMPONENTS = REGION%MESHES%MESHES(1)%PTR%NUMBER_OF_COMPONENTS
!     ALLOCATE(NODES_PER_MESH_COMPONENT(NUMBER_OF_MESH_COMPONENTS))
!     DO I=1,NUMBER_OF_MESH_COMPONENTS
!       NODES_PER_MESH_COMPONENT(I)=REGION%MESHES%MESHES(1)%PTR%TOPOLOGY(I)%PTR%NODES%NUMBER_OF_NODES
!     END DO
!  
!     ! POINT POINTERS TO DATA
!     NULLIFY(COORDINATES,DEFORMED,FIBRES)
!     COORDINATES => REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field% &
!         & variables(1)%parameter_sets%parameter_sets(1)%PTR
!     DEFORMED => REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field% &
!         & variables(1)%parameter_sets%parameter_sets(1)%PTR
!     FIBRES => REGION%equations_sets%equations_sets(1)%ptr%geometry%fibre_field% &
!         & variables(1)%parameter_sets%parameter_sets(1)%PTR
! 
!     ! BEGIN WRITING TO FILE BY OPENING IT
!     OPEN(UNIT=14, FILE=CHAR(FILENAME)//'.exnode',STATUS='unknown')
!     WRITE(*,*) "Writing Nodes..."
! 
!     ! WRITE HEADER INFORMATION
!     WRITE(14,*) 'Group name: FiniteElasticity'
!     WRITE(14,*) '#Fields=3' ! TODO: FIX THIS AT 3 FOR NOW
! 
!     ! 1) Coordinate field
!     VALUE_IDX=1
!     WRITE(WORD,22) NUMBER_OF_DIMENSIONS
!     WRITE(14,*) ' 1) Coordinates,  coordinate, rectangular cartesian, #Components='//WORD
!     DO I=1,NUMBER_OF_DIMENSIONS
!       ! TODO: ENABLE EXPORTS OF CUBIC HERMITE ELEMENT NODES
!       WRITE(WORD,22) VALUE_IDX
!       WRITE(14,*) '   '//XYZ(I)//'.  Value index= '//TRIM(WORD)//',     #Derivatives= 0'
!       VALUE_IDX = VALUE_IDX + 1
!     END DO
! 
!     ! 2) Deformed coordinate field
!     WRITE(WORD,22) NUMBER_OF_DIMENSIONS
!     WRITE(14,*) ' 2) Deformed,  field,  rectangular cartesian, #Components='//WORD
!     DO I=1,NUMBER_OF_DIMENSIONS
!       WRITE(WORD,22) VALUE_IDX
!       WRITE(14,*)  '   '//XYZ(I)//'.  Value index= '//TRIM(WORD)//',     #Derivatives= 0' 
!       VALUE_IDX = VALUE_IDX + 1
!     END DO
! 
!     ! 3) Fibre field
!     WRITE(14,*) ' 3) Fibre, anatomical, fibre, #Components=3' ! TODO: grab no. of components from Region
!     WRITE(WORD,22) VALUE_IDX
!     WRITE(14,*) '      fibre angle.  Value index='//TRIM(WORD)//',     #Derivatives= 0'
!     VALUE_IDX = VALUE_IDX + 1
!     WRITE(WORD,22) VALUE_IDX
!     WRITE(14,*) '      imbrication angle.  Value index='//TRIM(WORD)//',     #Derivatives= 0'
!     VALUE_IDX = VALUE_IDX + 1
!     WRITE(WORD,22) VALUE_IDX
!     WRITE(14,*) '      sheet angle.  Value index='//TRIM(WORD)//',     #Derivatives= 0'
!     VALUE_IDX = VALUE_IDX + 1
! 
!     ! NOW WRITE NODE INFORMATION
!     DO I = 1,NODES_PER_MESH_COMPONENT(1)
!       WRITE(14,*) ' Node: ',I
!       ! Coordinate field
!       WRITE(14,'("    ", es25.16 )') COORDINATES%PARAMETERS%CMISS%DATA_DP(I)
!       WRITE(14,'("    ", es25.16 )') COORDINATES%PARAMETERS%CMISS%DATA_DP(I+NODES_PER_MESH_COMPONENT(1))
!       IF(NUMBER_OF_DIMENSIONS==3) &
!     & WRITE(14,'("    ", es25.16 )') COORDINATES%PARAMETERS%CMISS%DATA_DP(I+2*NODES_PER_MESH_COMPONENT(1))
!       ! Deformed coordinate field
!       WRITE(14,'("    ", es25.16 )') DEFORMED%PARAMETERS%CMISS%DATA_DP(I)
!       WRITE(14,'("    ", es25.16 )') DEFORMED%PARAMETERS%CMISS%DATA_DP(I+NODES_PER_MESH_COMPONENT(1))
!       IF(NUMBER_OF_DIMENSIONS==3) &
!     & WRITE(14,'("    ", es25.16 )') DEFORMED%PARAMETERS%CMISS%DATA_DP(I+2*NODES_PER_MESH_COMPONENT(1))
!       ! Fibre field
!       WRITE(14,'("    ", es25.16 )') FIBRES%PARAMETERS%CMISS%DATA_DP(I)
!       WRITE(14,'("    ", es25.16 )') FIBRES%PARAMETERS%CMISS%DATA_DP(I+NODES_PER_MESH_COMPONENT(1))
!       IF(NUMBER_OF_DIMENSIONS==3) &
!     & WRITE(14,'("    ", es25.16 )') FIBRES%PARAMETERS%CMISS%DATA_DP(I+2*NODES_PER_MESH_COMPONENT(1))
!     END DO
!  
!     WRITE(14,*)
!     CLOSE(14)
! 
!     DEALLOCATE(NODES_PER_MESH_COMPONENT)
! 
!     CALL EXITS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS")
!     RETURN
! 999 CALL ERRORS("SOLID_MECHANICS_IO_WRITE_CMGUI_NODES",ERR,ERROR)    
!     CALL EXITS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS")
!     RETURN
! 
!   END SUBROUTINE SOLID_MECHANICS_IO_WRITE_CMGUI_NODES

!   ! =============================================================================
!   !> export existing nodes to a cmgui file
!   SUBROUTINE SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS(FILENAME,REGION,ERR,ERROR,*)
!     TYPE(VARYING_STRING), INTENT(IN) :: FILENAME !<file path without extension
!     TYPE(REGION_TYPE), POINTER :: REGION !<pointer to the region to get the coordinate system for
!     INTEGER(INTG), INTENT(OUT) :: ERR                   !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR          !<The error string
!     ! local variables
!     CHARACTER*5 :: WORD
!     INTEGER(INTG) :: VALUE_IDX,I,J
!     CHARACTER*1 :: XYZ(3) = (/'x','y','z'/)
!     ! grab below from Region
!     INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_FIBRE_COMPONENTS
!     INTEGER(INTG) :: NUMBER_OF_MESH_COMPONENTS
!     INTEGER(INTG) :: SCALE_FACTORS_OFFSET
!     INTEGER(INTG),ALLOCATABLE :: NODES_PER_MESH_COMPONENT(:)
!     INTEGER(INTG),ALLOCATABLE :: NODES_PER_ELEMENT(:)
!     TYPE(FIELD_VARIABLE_TYPE),POINTER :: COORDINATES,DEFORMED,FIBRES
!     TYPE(BASIS_TYPE),POINTER :: BASIS
!     INTEGER(INTG),ALLOCATABLE :: BASIS_TYPES(:)
!     TYPE(MESH_TOPOLOGY_TYPE),POINTER :: TOPOLOGY
! 
!     CALL ENTERS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS",ERR,ERROR,*999)
! 
! 22 FORMAT(I5) ! for use with WORD and WORD2
! 
!     ! scan REGION for parameters
!     ! TODO: CURRENTLY ASSUMES 1 MESH COMPONENT...!
!     NUMBER_OF_DIMENSIONS = REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
!     NUMBER_OF_FIBRE_COMPONENTS = REGION%EQUATIONS_SETS%EQUATIONS_SETS(1)%PTR%GEOMETRY% &
!       & FIBRE_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
!     NUMBER_OF_MESH_COMPONENTS = REGION%MESHES%MESHES(1)%PTR%NUMBER_OF_COMPONENTS
! !     ALLOCATE(NODES_PER_MESH_COMPONENT(NUMBER_OF_MESH_COMPONENTS))
! !     ALLOCATE(NODES_PER_ELEMENT(NUMBER_OF_MESH_COMPONENTS))
! !     DO I=1,NUMBER_OF_MESH_COMPONENTS
! !       NODES_PER_MESH_COMPONENT(I)=REGION%MESHES%MESHES(1)%PTR%TOPOLOGY(I)%PTR%NODES%NUMBER_OF_NODES
! !       NODES_PER_ELEMENT(I)=REGION%FIELDS%FIELDS(1)%PTR%GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(1) &
! !         & %PTR%TOPOLOGY%ELEMENTS%ELEMENTS(1)%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
! !     END DO
! 
!     ! POINT POINTERS TO DATA
!     ! these three guys point to FIELD-related data
!     NULLIFY(COORDINATES,DEFORMED,FIBRES)
!     COORDINATES => REGION%equations_sets%equations_sets(1)%ptr%geometry%geometric_field% &
!         & variables(1) !%parameter_sets%parameter_sets(1)%PTR
!     DEFORMED => REGION%equations_sets%equations_sets(1)%ptr%dependent%dependent_field% &
!         & variables(1) !%parameter_sets%parameter_sets(1)%PTR
!     FIBRES => REGION%equations_sets%equations_sets(1)%ptr%geometry%fibre_field% &
!         & variables(1) !%parameter_sets%parameter_sets(1)%PTR
!     ! BASIS pointers, and test them now for unimplemented stuff
!     ALLOCATE(BASIS_TYPES(NUMBER_OF_MESH_COMPONENTS))
!     DO I=1,NUMBER_OF_MESH_COMPONENTS
!       NULLIFY(BASIS)
!       BASIS=>REGION%MESHES%MESHES(1)%PTR%TOPOLOGY(I)%PTR%ELEMENTS%ELEMENTS(1)%BASIS ! look at TYPE and HERMITE
!       BASIS_TYPES(I)=BASIS%INTERPOLATION_TYPE(1)
!       ! test them
!       IF(BASIS%TYPE/=BASIS_LAGRANGE_HERMITE_TP_TYPE) THEN
!         CALL FLAG_ERROR("Only Lagrange Hermite type basis is implemented",ERR,ERROR,*999)
!       ENDIF
!       ! do nothing for implemented bases
!       IF(BASIS_TYPES(I)==BASIS_LAGRANGE_HERMITE_INTERPOLATION) THEN
!       ELSE
!         CALL FLAG_ERROR("this type of basis is not implemented currently",ERR,ERROR,*999)
!       ENDIF
!     ENDDO
!     ! TOPOLOGY pointer
!     TOPOLOGY => region%meshes%meshes(1)%ptr%topology(1)%ptr
! !     TOPOLOGY => REGION%MESHES%MESHES(1)%PTR%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR%TOPOLOGY
! 
!     ! OPEN FILE
!     OPEN(UNIT=14, FILE=CHAR(FILENAME)//'.exelem',STATUS='unknown')
!     WRITE(*,*) "Writing Elements..."
!     WRITE(14,*) " Group name: FiniteElasticity"
! 
!     ! Actually, cmgui can generate lines and faces itself! (no need to export those)
!     ! (do nothing)
! 
!     ! scale factors are basis dependent
!     ! currently, 1) coordinates 2) deformed 3) fibres 4) pressure (not yet?)
!     WRITE(14,*) ' Shape.  Dimension=3'
!     WRITE(WORD,22) NUMBER_OF_MESH_COMPONENTS
!     WRITE(14,*) ' #Scale factor sets= '//TRIM(WORD)
!     DO I=1,NUMBER_OF_MESH_COMPONENTS
!       IF(BASIS_TYPES(I)==BASIS_HERMITE_INTERPOLATION) THEN
!         WRITE(14,*) '   c.Hermite*c.Hermite*c.Hermite, #Scale factors= 64'
!       ELSEIF(BASIS_TYPES(I)==BASIS_LAGRANGE_INTERPOLATION) THEN
!         WRITE(14,*) '   l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors= 8'
!       ENDIF
!     ENDDO
!     ! TODO: BELOW TWO LINES ARE HARD-CODED
!     WRITE(14,*) ' #Nodes=           8'
!     WRITE(14,*) ' #Fields=3'
! 
!     ! ----------------------
!     ! fields are the more complicated
!     WRITE(14,*) ' 1) Coordinates, coordinate, rectangular cartesian, #Components=3'
!     ! grab the field and have a look at its basis
!     I=COORDINATES%COMPONENTS(1)%MESH_COMPONENT_NUMBER
!     IF(BASIS_TYPES(I)==BASIS_HERMITE_INTERPOLATION) THEN
!       DO J=1,3
!         WRITE(14,*) "   "//XYZ(J)//".  c.Hermite*c.Hermite*c.Hermite, no modify, standard node based."
!         CALL PRINT_FIELD_INFO(14,1,0)
!       ENDDO
!     ELSEIF(BASIS_TYPES(I)==BASIS_LAGRANGE_INTERPOLATION) THEN
!       DO J=1,3
!         WRITE(14,*) "   "//XYZ(J)//".  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based."
!         CALL PRINT_FIELD_INFO(14,2,0)
!       ENDDO
!     ENDIF
!     WRITE(14,*) ' 2) Deformed, field, rectangular cartesian, #Components=3'
!     IF(BASIS_TYPES(i)==BASIS_HERMITE_INTERPOLATION) THEN
!       DO J=1,3
!         WRITE(14,*) "   "//XYZ(J)//".  c.Hermite*c.Hermite*c.Hermite, no modify, standard node based."
!         CALL PRINT_FIELD_INFO(14,1,0)
!       ENDDO
!     ELSEIF(BASIS_TYPES(I)==BASIS_LAGRANGE_INTERPOLATION) THEN
!       DO J=1,3
!         WRITE(14,*) "   "//XYZ(J)//".  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based."
!         CALL PRINT_FIELD_INFO(14,2,0)
!       ENDDO
!     ENDIF
!     WRITE(14,*) ' 3) Fibre, anatomical, fibre, #Components=3'
!     I=FIBRES%COMPONENTS(1)%MESH_COMPONENT_NUMBER
!     IF(BASIS_TYPES(I)/=BASIS_TYPES(1).AND.BASIS_TYPES(1)==BASIS_HERMITE_INTERPOLATION) THEN
!       SCALE_FACTORS_OFFSET=64
!     ELSEIF(BASIS_TYPES(I)/=BASIS_TYPES(1).AND.BASIS_TYPES(1)==BASIS_LAGRANGE_INTERPOLATION) THEN
!       SCALE_FACTORS_OFFSET=8
!     ENDIF
!     IF(BASIS_TYPES(I)==BASIS_HERMITE_INTERPOLATION) THEN
!       WRITE(14,*) "   fibre angle.  c.Hermite*c.Hermite*c.Hermite, no modify, standard node based."
!       CALL PRINT_FIELD_INFO(14,1,SCALE_FACTORS_OFFSET)
!       WRITE(14,*) "   imbrication angle.  c.Hermite*c.Hermite*c.Hermite, no modify, standard node based."
!       CALL PRINT_FIELD_INFO(14,1,SCALE_FACTORS_OFFSET)
!       WRITE(14,*) "   sheet angle.  c.Hermite*c.Hermite*c.Hermite, no modify, standard node based."
!       CALL PRINT_FIELD_INFO(14,1,SCALE_FACTORS_OFFSET)
!     ELSEIF(BASIS_TYPES(I)==BASIS_LAGRANGE_INTERPOLATION) THEN
!       WRITE(14,*) "   fibre angle.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based."
!       CALL PRINT_FIELD_INFO(14,2,SCALE_FACTORS_OFFSET)
!       WRITE(14,*) "   imbrication angle.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based."
!       CALL PRINT_FIELD_INFO(14,2,SCALE_FACTORS_OFFSET)
!       WRITE(14,*) "   sheet angle.  l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based."
!       CALL PRINT_FIELD_INFO(14,2,SCALE_FACTORS_OFFSET)
!     ENDIF
! 
!     ! ----------------------
!     ! phew. finally get onto printing element information
!     ! Element: 1 0 0
!     ! Nodes: 1 2 3 4 5 6 7 8
!     ! Scale Factors: etc etc
!     DO I=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
!       WRITE(WORD,22) I
!       WRITE(14,*) " Element:            "//TRIM(WORD)//" 0 0"
!       WRITE(14,*) "   Nodes:"
!       WRITE(14,*) TOPOLOGY%ELEMENTS%ELEMENTS(I)%MESH_ELEMENT_NODES
!       WRITE(14,*) "   Scale factors:"
!       IF(BASIS_TYPES(1)==BASIS_HERMITE_INTERPOLATION) THEN
! !         WRITE(14,*) 
!         !TODO: WRITE OUT THE CORRECT SCALINGS
!       ELSEIF(BASIS_TYPES(1)==BASIS_LAGRANGE_INTERPOLATION) THEN
!         WRITE(14,*) "     1 1 1 1 1 1 1 1"
!       ENDIF
!     ENDDO
! 
!     CALL EXITS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS")
!     RETURN
! 999 CALL ERRORS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS",ERR,ERROR)    
!     CALL EXITS("SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS")
!     RETURN
! 
!   END SUBROUTINE SOLID_MECHANICS_IO_WRITE_CMGUI_ELEMENTS

  ! =============================================================================
  SUBROUTINE PRINT_FIELD_INFO(UNIT,TYPE,OFFSET)
    INTEGER(INTG),INTENT(IN) :: UNIT    ! FILE UNIT
    INTEGER(INTG),INTENT(IN) :: TYPE    ! 1=CUBIC HERMITE 2=TRILINEAR
    INTEGER(INTG),INTENT(IN) :: OFFSET  ! OFFSET FOR SCALE FACTOR INDICES
    ! local variables
    INTEGER(INTG) :: I, N_STEP
    INTEGER(INTG),ALLOCATABLE :: N(:),VALUE_INDICES(:)
    CHARACTER*2 :: WORD,LENG

    IF(TYPE==1) THEN     ! CUBIC HERMITE
      ALLOCATE(N(8),VALUE_INDICES(8))
      N=(/1,2,3,4,5,6,7,8/)
      VALUE_INDICES=N
      N_STEP=8
      WRITE(LENG,'(I2)') 8
    ELSEIF(TYPE==2) THEN ! TRILINEAR
      ALLOCATE(N(1),VALUE_INDICES(1))
      N=OFFSET+1
      VALUE_INDICES=N
      N_STEP=1
      WRITE(LENG,'(I2)') 1
    ELSE
      WRITE(*,*) "SOLID_MECHANICS_IO: ERROR IN PRINT_FIELD_INFO"
    ENDIF

    WRITE(UNIT,*) "     #Nodes= 8"
    DO I=1,8
      WRITE(WORD,'(I2)') I
      WRITE(UNIT,*) "      "//TRIM(WORD)//".  #Values="//LENG
      WRITE(UNIT,'("       Value indices:     ",8I4:)') VALUE_INDICES
      WRITE(UNIT,'("       Scale factor indices: ",8I4:)') N
      N=N+N_STEP
    ENDDO
 
  END SUBROUTINE PRINT_FIELD_INFO

  ! =============================================================================


END MODULE SOLID_MECHANICS_IO_ROUTINES



