!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all Hamilton-Jacobi equations routines.
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

!>This module handles all Hamilton-Jacobi equations routines.
MODULE HAMILTON_JACOBI_EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE COMP_ENVIRONMENT
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE MATHS
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE REGION_ROUTINES
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

!!MERGE: move

  !PUBLIC HJ_EQUATION_FINITE_ELEMENT_CALCULATE,HJ_EQUATION_EQUATIONS_SET_SETUP, &
  !  & HJ_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET,HJ_EQUATION_EQUATIONS_SET_SUBTYPE_SET, &
  !  & HJ_EQUATION_PROBLEM_SUBTYPE_SET,HJ_EQUATION_PROBLEM_SETUP

  PUBLIC HJ_EQUATION_PROBLEM_SUBTYPE_SET

  PUBLIC NUMBER_OF_INPUT_NODES,PRE_PROCESS_INFORMATION,SOLVE_PROBLEM_FMM,SOLVE_PROBLEM_GEODESIC
  PUBLIC SOLVE_PROBLEM_GEODESIC_CONNECTIVITY,SOLVE_PROBLEM_FMM_CONNECTIVITY
  PUBLIC FIND_MINIMAX,POST_PROCESS_DATA
  
  PUBLIC READ_TETGEN_MESH,WRITE_VTK_MESH

CONTAINS

  !
  !================================================================================================================================
  !


  !>Sets/changes the problem subtype for a Laplace equation type .
  SUBROUTINE HJ_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("HJ_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_HJ_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_FMM_CLASS
        PROBLEM%TYPE=PROBLEM_HJ_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_HJ_SUBTYPE     
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a HJ equation type."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("HJ_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("HJ_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("HJ_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE HJ_EQUATION_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Calculates to give back the number of nodes from input file.
  SUBROUTINE NUMBER_OF_INPUT_NODES(INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,&
  &TOTAL_NUMBER_OF_CONNECTIVITY,Err)

    !subroutine variables
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_CONNECTIVITY
    CHARACTER (LEN=300) :: INPUT_FILE_NAME
    CHARACTER (LEN=10) :: INPUT_FILE_FORMAT
    INTEGER(INTG) :: Err

    !Argument variables
!    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string

    !Local variables
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:):: CONNECTIVITY_LIST
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:):: ELEMENT_LIST
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: CONNECTIVITY_NUMBER
    INTEGER(INTG) :: I,J,K,N,NUMBER_OF_NODES_PER_ELEMENT,THERE_IS_IN_CONNECTIVITY_LIST
    CHARACTER (LEN=300) :: STRING
        
!    CALL ENTERS("GENERATE_STATUS_MASK",Err)
    
! """""""""""""""""""""""""""""""""""INPUT OF TABC FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TABC") THEN 

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".tabc"
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
!      PRINT *, STRING
      TOTAL_NUMBER_OF_NODES=-1
      DO WHILE (STRING .ne. "Connectivity") 
        READ(11,*) STRING
        TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES+1
      ENDDO
      
      TOTAL_NUMBER_OF_CONNECTIVITY=0
      DO I=1,TOTAL_NUMBER_OF_NODES
        READ(11,*) STRING,J
        TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+J
      ENDDO
      
      CLOSE (11)
      TOTAL_NUMBER_OF_ELEMENTS = TOTAL_NUMBER_OF_CONNECTIVITY ! SHOULD BE DEFINED
    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".vtk"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_NODES
!      PRINT *, STRING,TOTAL_NUMBER_OF_NODES
      
      DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)
      
        READ(11,*) STRING

      ENDDO
!      READ(11,*) STRING
!      print*,I,INT(TOTAL_NUMBER_OF_NODES/3.0+0.5),STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_ELEMENTS
    
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET1NPL") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".vtk"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_NODES,STRING
!      PRINT *, STRING,TOTAL_NUMBER_OF_NODES,STRING
      
      DO I=1,TOTAL_NUMBER_OF_NODES
      
        READ(11,*) STRING

      ENDDO
!      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_ELEMENTS

      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF CARP FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "CARP") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".pts"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_NODES
!      PRINT *, TOTAL_NUMBER_OF_NODES
      CLOSE (11)
      
      STRING = INPUT_FILE_NAME(1:I)//".elem"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF TETGEN FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TETGEN") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".node"
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_NODES
      CLOSE (11)
      
      STRING = INPUT_FILE_NAME(1:I)//".ele"
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(TOTAL_NUMBER_OF_NODES,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(TOTAL_NUMBER_OF_NODES),STAT=ERR)
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE (11)
      
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

    ENDIF


  END SUBROUTINE NUMBER_OF_INPUT_NODES

  !
  !================================================================================================================================
  !


  !>to READ input file.
  SUBROUTINE PRE_PROCESS_INFORMATION(MATERIAL_BEHAVIOUR,INPUT_FILE_NAME,INPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,&
&INPUT_TYPE_FOR_SEED_VALUE,INPUT_TYPE_FOR_SPEED_FUNCTION,SPEED_FUNCTION_ALONG_EIGEN_VECTOR,INPUT_TYPE_FOR_CONDUCTIVITY,&
&STATUS_MASK,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,&
&SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

    !subroutine variables
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LIST(:,:)

    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: COLUMN_INDEX
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: RAW_INDEX
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    REAL(DP) :: SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS,TOTAL_NUMBER_OF_CONNECTIVITY
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NODES_PER_ELEMENT
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SEED_VALUE
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SPEED_FUNCTION
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_CONDUCTIVITY
    CHARACTER (LEN=300) :: INPUT_FILE_NAME
    CHARACTER (LEN=10)  :: INPUT_FILE_FORMAT
    CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    
    !Local variables
    CHARACTER (LEN=300) :: STRING
    INTEGER(INTG) :: I,J,K,N,FIRST_NODE_NUMBER
    INTEGER(INTG) :: TEXT_LENGTH, THERE_IS_IN_CONNECTIVITY_LIST
    REAL(DP) :: A(3),B(3),C(3)
    REAL(DP) :: DOT_PRODUCT_VALUE
        
!INITIALIZE PARAMETERS:
    DO I=1,TOTAL_NUMBER_OF_NODES
      CONNECTIVITY_NUMBER(I)=0
      DO J=1,3
        NODE_LIST(I,J) = 0.0
        SPEED_FUNCTION_TABLE(I,J) = 0.0
      ENDDO
      DO J=1,9
        CONDUCTIVITY_TENSOR(I,J) = 0.0
      ENDDO
      RAW_INDEX(I)=0
    ENDDO
    RAW_INDEX(TOTAL_NUMBER_OF_NODES+1)=0

    DO I=1,TOTAL_NUMBER_OF_CONNECTIVITY
      COLUMN_INDEX(I)=0
      DO J=1,3
        SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(I,J) = 0.0
      ENDDO
      DO J=1,9
        CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(I,J) = 0.0
      ENDDO
    ENDDO
     
! """""""""""""""""""""""""""""""""""INPUT OF TABC FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TABC") THEN 

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".tabc"
      OPEN (11,FILE=STRING)
      READ(11,*) STRING

! SOSIOISOISOISOISOIS      load data for the case material behaves = ISOTROPIC 
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!      input type for *velocity function* = FILE and *seed points* = FILE
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),SEED_VALUE(I)
            
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
           
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

          ENDDO
        ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SEED_VALUE(I)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

          ENDDO
        ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

            IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
              SEED_VALUE(I) = 1000.0_DP
            ENDIF

          ENDDO
        ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

            IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
              SEED_VALUE(I) = 1000.0_DP
            ENDIF

          ENDDO
        ENDIF

      ENDIF

! ANANANAIANSOANAIANI      load data for the case material behaves = ANISOTROPIC 
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 


!      if conductivity format is TENSOR type i.e. three EigenVectors
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN

!      input type for *velocity function* = FILE and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3),SEED_VALUE(I)

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SEED_VALUE(I)
 
            ENDDO
          ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

        ENDIF


!      if conductivity format is VECTOR type i.e. first EigenVectors
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN

!      input type for *velocity function* = FILE and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3),SEED_VALUE(I)

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SEED_VALUE(I)
 
            ENDDO
          ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,TOTAL_NUMBER_OF_NODES 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

          DO I=1,TOTAL_NUMBER_OF_NODES
!            CALL CALCULATE_SECOND_EIGENVECTOR()
            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            CALL CALCULATE_SECOND_EIGENVECTOR()
          ENDDO

        ENDIF

      ENDIF

! CONNSDONCOCNCNOSKCN      load data for the CONNECTIVITY list 
      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_NODES

        READ(11,*) STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
        
        DO J=1,3
          SPEED_FUNCTION_TABLE(I,J)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
        ENDDO
          
!        PRINT *,STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))

      ENDDO

      CLOSE(11)

    ENDIF


! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"

      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING

      DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-1

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3),(NODE_LIST(3*(I-1)+3,J),J=1,3)

      ENDDO

      I=INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)

      IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 0) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3),(NODE_LIST(3*(I-1)+3,J),J=1,3)

      ENDIF
      IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 1) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3)


      ENDIF
      IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 2) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3)

      ENDIF
      

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 

      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO
        
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO
        
      ENDDO
      
      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        DO I=1,TOTAL_NUMBER_OF_NODES
          STATUS_MASK(I) = ""
        ENDDO

        OPEN (11,FILE=STRING)
        READ (11,*) N

        DO I=1,N 

          READ(11,*) J,SEED_VALUE(J+1)
!          PRINT*,(NODE_LIST(J+1,K),K=1,3),SEED_VALUE(J+1)      
          STATUS_MASK(J+1) = "SEED POINT"

        ENDDO

        CLOSE(11)

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          IF (STRING .EQ. '#') THEN! begin if

            DO WHILE (STRING .NE. 'fiber')
      
              READ(11,*) STRING
        
            ENDDO
      
!      PRINT *,STRING
      
            DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-1

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDDO

            I=INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)

            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 0) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 1) THEN
    
              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 2) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3)

            ENDIF
      
          
          ELSE
            DO I=1,TOTAL_NUMBER_OF_NODES 
            
              READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3) 
              
            ENDDO               
          ENDIF ! end if
          
          DO I=1,TOTAL_NUMBER_OF_NODES 
          
!            READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)

            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 
        
        DO I=1,TOTAL_NUMBER_OF_NODES
          DO J=1,3
            SPEED_FUNCTION_TABLE(I,J) = 0.0_DP
          ENDDO
        ENDDO
        
          DO I=1,TOTAL_NUMBER_OF_NODES
!            DO J=1,CONNECTIVITY_NUMBER(I)
!              IF (CONNECTIVITY_LIST(I,J) > TOTAL_NUMBER_OF_NODES) THEN
!                PRINT*, I,J,CONNECTIVITY_NUMBER(I),CONNECTIVITY_LIST(I,J),TOTAL_NUMBER_OF_NODES
!              ENDIF
              SPEED_FUNCTION_TABLE(I,1) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
              SPEED_FUNCTION_TABLE(I,2) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)
              SPEED_FUNCTION_TABLE(I,3) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),1) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),2) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),3) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
!            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET1NPL") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"

      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_NODES

        READ(11,*) (NODE_LIST(I,J),J=1,3)

      ENDDO

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 

      READ(11,*) STRING

      DO N=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) N

        DO I=1,N 

          READ(11,*) J,SEED_VALUE(J+1)
       
          STATUS_MASK(J+1) = "SEED POINT"

        ENDDO

        CLOSE(11)

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING
!          print *,STRING
        
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          IF (STRING .EQ. '#') THEN! begin if

            DO WHILE (STRING .NE. 'fiber')
      
              READ(11,*) STRING
        
            ENDDO
      
!      PRINT *,STRING
      
            DO I=1,INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-1

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDDO

            I=INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)

            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 0) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 1) THEN
    
              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3)

            ENDIF
            IF (3*INT((TOTAL_NUMBER_OF_NODES/3.0_DP)+0.5_DP)-TOTAL_NUMBER_OF_NODES .EQ. 2) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3)

            ENDIF
      
          
          ELSE
            DO I=1,TOTAL_NUMBER_OF_NODES 
            
              READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3) 
              
            ENDDO               
          ENDIF ! end if
          
          DO I=1,TOTAL_NUMBER_OF_NODES 
          
!            READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF CARP FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "CARP") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".pts"

      OPEN (11,FILE=STRING)
      READ (11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_NODES 

        READ(11,*) (NODE_LIST(I,J),J=1,3)

      ENDDO

      CLOSE(11)

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".elem"

      OPEN (11,FILE=STRING)
      READ (11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      DO N=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        DO I=1,TOTAL_NUMBER_OF_NODES 
 
          READ(11,*) STRING,SEED_VALUE(I)

        ENDDO

        CLOSE(11)

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".lon"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR" .OR. INPUT_TYPE_FOR_CONDUCTIVITY .EQ. " ") THEN
          DO I=1,TOTAL_NUMBER_OF_ELEMENTS
            READ(11,*) (CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),1),CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),2), &
             & CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE
                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),4),CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),5), &
             & CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF
            
            DO J=2,NUMBER_OF_NODES_PER_ELEMENT
              DO K=1,9
                CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,J),K)=CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),K)
              ENDDO
            ENDDO

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF
    
! """""""""""""""""""""""""""""""""""INPUT OF TETGEN FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TETGEN") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".node"

      OPEN (11,FILE=STRING)
      READ (11,*) STRING

      READ(11,*) FIRST_NODE_NUMBER,(NODE_LIST(1,J),J=1,3)

      DO I=2,TOTAL_NUMBER_OF_NODES 


        READ(11,*) STRING,(NODE_LIST(I,J),J=1,3)


      ENDDO

      CLOSE(11)

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".ele"

      OPEN (11,FILE=STRING)
      READ (11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      DO N=1,TOTAL_NUMBER_OF_NODES
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        IF (FIRST_NODE_NUMBER .EQ. 0) THEN        
          DO J=1,NUMBER_OF_NODES_PER_ELEMENT
            ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
          ENDDO
        ENDIF

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,TOTAL_NUMBER_OF_NODES 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        DO I=1,TOTAL_NUMBER_OF_NODES 
 
          READ(11,*) STRING,SEED_VALUE(I)

        ENDDO

        CLOSE(11)

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          DO I=1,TOTAL_NUMBER_OF_NODES 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=(/CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)/)
            B=(/0.0_DP,0.0_DP,1.0_DP/)
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
            ELSE
              B=(/0.0_DP,1.0_DP,0.0_DP/)
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ELSE

                B=(/1.0_DP,0.0_DP,0.0_DP/)
                CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=(/CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)/)
            CALL CROSS_PRODUCT(A,B,C,Err,Error,*999)

            IF (SQRT(C(1)**2+C(2)**2+C(3)**2) .NE. 0.0_DP) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,TOTAL_NUMBER_OF_NODES
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF

    ! to set RAW_INDEX and COLUMN_INDEX
    RAW_INDEX(1) = 0
    DO I=1,TOTAL_NUMBER_OF_NODES
        
      RAW_INDEX(I+1) = RAW_INDEX(I) + CONNECTIVITY_NUMBER(I)
      DO J = 1,CONNECTIVITY_NUMBER(I)
        COLUMN_INDEX(RAW_INDEX(I)+J) = CONNECTIVITY_LIST(I,J)
      ENDDO          

    ENDDO
    
    DO I=1,TOTAL_NUMBER_OF_NODES
        
      DO J = RAW_INDEX(I)+1,RAW_INDEX(I+1)
        DO K = 1,3
          SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,K) = &
          & (SPEED_FUNCTION_TABLE(I,K)+SPEED_FUNCTION_TABLE(COLUMN_INDEX(J),K))/2.0_DP
        ENDDO
        DO K = 1,9
          CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,K) = &
          & (CONDUCTIVITY_TENSOR(I,K)+CONDUCTIVITY_TENSOR(COLUMN_INDEX(J),K))/2.0_DP
        ENDDO
      ENDDO

    ENDDO

!    CALL EXITS("GENERATE_STATUS_MASK")
!    RETURN
999 CALL ERRORS("GENERATE_STATUS_MASK",ERR,ERROR)
    CALL EXITS("GENERATE_STATUS_MASK")
    RETURN

  END SUBROUTINE PRE_PROCESS_INFORMATION


  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_FMM(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,&
  &SEED_VALUE,CONNECTIVITY_NUMBER,CONNECTIVITY_LIST,STATUS_MASK)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER, CHANGED_STATUS, MIN_TRIAL_NODE, TRIAL_STATUS
    REAL(DP), DIMENSION(3)   :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,DET,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW,CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    !Start Program

    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

!    MAIN LOOP 
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=1,CONNECTIVITY_NUMBER(MIN_TRIAL_NODE)
        TIME_NEW=1000

        DO J=1,CONNECTIVITY_NUMBER(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))
          IF (STATUS_MASK(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=(&
                             &/NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),1)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),2)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),3)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3)/)

!            CONDUCTION_RATIO=SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)/SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)
!            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE = (/1.0_DP,0.0_DP,0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO/), SHAPE = (/3,3/))
            CONDUCTION_RATIO= &
            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)/&
            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)
!            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO/), SHAPE = (/3,3/))


!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),4),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),5),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),6),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),7),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),8),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),9)/&
                               &),SHAPE = (/3,3/))

            CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)
!            CALL INVERT(FMFT,INV_FMFT,DET,Err,Error,*999)

!	    PRINT *,F(1,1),F(1,2),F(1,3),F(2,1),F(2,2),F(2,3),F(3,1),F(3,2),F(3,3)

            CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J))+SQRT(ABS(VMV))*&
          &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)
!          &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)


            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .lt. SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))) THEN
          SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = TIME_NEW
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "KNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0.0_DP

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0.AND.STATUS_MASK(I) .EQ. "SEED POINT".AND.SEED_VALUE(I).LT.MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


999 CALL ERRORS("SOLVE_PROBLEM_FMM",ERR,ERROR)
    CALL EXITS("SOLVE_PROBLEM_FMM")
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_FMM

  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_FMM_CONNECTIVITY(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDEX(:)  
    INTEGER(INTG), ALLOCATABLE :: RAW_INDEX(:)    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)

    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER, CHANGED_STATUS, MIN_TRIAL_NODE, TRIAL_STATUS
    REAL(DP), DIMENSION(3)   :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(2)   :: CONDUCTION_RATIO
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,DET,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    !Start Program

    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

! ::::::: MAIN LOOP :::::::	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=RAW_INDEX(MIN_TRIAL_NODE)+1,RAW_INDEX(MIN_TRIAL_NODE+1)
        TIME_NEW=1000

        DO J=RAW_INDEX(COLUMN_INDEX(I))+1,RAW_INDEX(COLUMN_INDEX(I)+1)
          IF (STATUS_MASK(COLUMN_INDEX(J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=(/NODE_LIST(COLUMN_INDEX(I),1)-NODE_LIST(COLUMN_INDEX(J),1)&
                            &,NODE_LIST(COLUMN_INDEX(I),2)-NODE_LIST(COLUMN_INDEX(J),2)&
                            &,NODE_LIST(COLUMN_INDEX(I),3)-NODE_LIST(COLUMN_INDEX(J),3)/)

            CONDUCTION_RATIO(1)= SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,2)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
            CONDUCTION_RATIO(2)= SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,3)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO(1)*CONDUCTION_RATIO(1),0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO(2)*CONDUCTION_RATIO(2)/), SHAPE = (/3,3/))


!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,1),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,2),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,3),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,4),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,5),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,6),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,7),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,8),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,9)/),SHAPE = (/3,3/))

            CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)
!            CALL INVERT(FMFT,INV_FMFT,DET,Err,Error,*999)

            CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(COLUMN_INDEX(J))+SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .LT. SEED_VALUE(COLUMN_INDEX(I))) THEN
          SEED_VALUE(COLUMN_INDEX(I)) = TIME_NEW
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "KNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0 .AND. STATUS_MASK(I) .EQ. "SEED POINT" .AND. SEED_VALUE(I) .LT. MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


999 CALL ERRORS("SOLVE_PROBLEM_FMM_CONNECTIVITY",ERR,ERROR)
    CALL EXITS("SOLVE_PROBLEM_FMM_CONNECTIVITY")
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_FMM_CONNECTIVITY


  !
  !================================================================================================================================
  !

  SUBROUTINE SOLVE_PROBLEM_GEODESIC_CONNECTIVITY(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK,TRACE_NODE)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDEX(:)
    INTEGER(INTG), ALLOCATABLE :: RAW_INDEX(:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: TRACE_NODE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    REAL(DP), ALLOCATABLE :: CONNECTIVITY_WEIGHT(:)
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER,CHANGED_STATUS,MIN_TRIAL_NODE,TRIAL_STATUS,TRACE_NODE_NUMBER
    REAL(DP), DIMENSION(3) :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,DET,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW
    REAL(DP), DIMENSION(2) :: CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    
    !initialize data
    DO I=1,TOTAL_NUMBER_OF_NODES
      TRACE_NODE(I)=0
    ENDDO

    !Start Program
    ALLOCATE(CONNECTIVITY_WEIGHT(TOTAL_NUMBER_OF_CONNECTIVITY),STAT=ERR)
    DO I=1,TOTAL_NUMBER_OF_NODES
      DO J=RAW_INDEX(I)+1,RAW_INDEX(I+1)
        DISTANCE_VECTOR=(/NODE_LIST(I,1)-NODE_LIST(COLUMN_INDEX(J),1)&
                        &,NODE_LIST(I,2)-NODE_LIST(COLUMN_INDEX(J),2)&
                        &,NODE_LIST(I,3)-NODE_LIST(COLUMN_INDEX(J),3)/)

        CONDUCTION_RATIO(1) = SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,2)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
        CONDUCTION_RATIO(2) = SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,3)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

        CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                            &0.0_DP,CONDUCTION_RATIO(1)*CONDUCTION_RATIO(1),0.0_DP,&
                                            &0.0_DP,0.0_DP,CONDUCTION_RATIO(2)*CONDUCTION_RATIO(2)/), SHAPE = (/3,3/))

        F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,1),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,2),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,3),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,4),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,5),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,6),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,7),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,8),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,9)/),SHAPE = (/3,3/))

        CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
        CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
        CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)

        CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err)
        CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

        CONNECTIVITY_WEIGHT(J)=SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
      ENDDO
    ENDDO
          
    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

! ::::::: MAIN LOOP :::::::	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
!    LOOP_NUMBER = 0
    PRINT *,"Running GEODESIC Solver ...."

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=RAW_INDEX(MIN_TRIAL_NODE)+1,RAW_INDEX(MIN_TRIAL_NODE+1)
        TIME_NEW=1000

        DO J=RAW_INDEX(COLUMN_INDEX(I))+1,RAW_INDEX(COLUMN_INDEX(I)+1)
          IF (STATUS_MASK(COLUMN_INDEX(J)) == "KNOWN") THEN 

            TIME_ITER=SEED_VALUE(COLUMN_INDEX(J))+CONNECTIVITY_WEIGHT(J)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
              TRACE_NODE_NUMBER=COLUMN_INDEX(J)
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .LT. SEED_VALUE(COLUMN_INDEX(I))) THEN
          SEED_VALUE(COLUMN_INDEX(I)) = TIME_NEW
          TRACE_NODE(COLUMN_INDEX(I)) = TRACE_NODE_NUMBER
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "KNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE = I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0 .AND. STATUS_MASK(I) .EQ. "SEED POINT" .AND. SEED_VALUE(I) .LT. MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO

999 CALL ERRORS("SOLVE_PROBLEM_GEODESIC_CONNECTIVITY",ERR,ERROR)
    CALL EXITS("SOLVE_PROBLEM_GEODESIC_CONNECTIVITY")
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_GEODESIC_CONNECTIVITY

  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_GEODESIC(TOTAL_NUMBER_OF_NODES,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,&
  & SEED_VALUE,CONNECTIVITY_NUMBER,CONNECTIVITY_LIST,STATUS_MASK,TRACE_NODE)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    INTEGER(INTG), ALLOCATABLE :: TRACE_NODE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER,CHANGED_STATUS,MIN_TRIAL_NODE,TRIAL_STATUS,TRACE_NODE_NUMBER
    REAL(DP), DIMENSION(3) :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,DET,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW,CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0

    !Start Program

    CALL GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

!   4444444444444444 MAIN LOOP 55555555555555555	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=1,CONNECTIVITY_NUMBER(MIN_TRIAL_NODE)
        TIME_NEW=1000

        DO J=1,CONNECTIVITY_NUMBER(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))
          IF (STATUS_MASK(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=(&
                             &/NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),1)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),2)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),3)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3)/)

            CONDUCTION_RATIO=SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)/&
                            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=(/1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO/), SHAPE = (/3,3/))

!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =(/CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),4),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),5),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),6),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),7),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),8),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),9)/&
                               &),SHAPE = (/3,3/))

            CALL MATRIX_TRANSPOSE(F,FT,Err,Error,*999)
            CALL MATRIX_PRODUCT(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MATRIX_PRODUCT(F,MFT,FMFT,Err,Error,*999)

            CALL MATRIX_VECTOR_PRODUCT(FMFT,DISTANCE_VECTOR,MV,Err)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J))+&
                      &SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
              TRACE_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .lt. SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))) THEN
          SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = TIME_NEW
          TRACE_NODE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))=TRACE_NODE_NUMBER
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "KNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000
      CHANGED_STATUS = 0

      DO I=1,TOTAL_NUMBER_OF_NODES

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0.AND.STATUS_MASK(I) .EQ. "SEED POINT".AND.SEED_VALUE(I).LT.MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


999 CALL ERRORS("SOLVE_PROBLEM_GEODESIC",ERR,ERROR)
    CALL EXITS("SOLVE_PROBLEM_GEODESIC")
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_GEODESIC

  
  !
  !================================================================================================================================
  !

  !>Calculates status mask for the local nodes.
  SUBROUTINE GENERATE_STATUS_MASK(TOTAL_NUMBER_OF_NODES,SEED_VALUE,STATUS_MASK,Err)

    !subroutine variables
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)

    !Argument variables
    INTEGER(INTG) :: Err !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string

    !Local variables
    INTEGER(INTG) :: I
        
!    CALL ENTERS("GENERATE_STATUS_MASK",Err)
    
    DO I=1,TOTAL_NUMBER_OF_NODES

      IF (SEED_VALUE(I) .LT. 100.0_DP) THEN
        STATUS_MASK(I) = "SEED POINT"
      ELSE
        STATUS_MASK(I) = "UNKNOWN"
      ENDIF

    ENDDO

!    CALL EXITS("GENERATE_STATUS_MASK")
!    RETURN
!999 CALL ERRORS("GENERATE_STATUS_MASK",ERR,ERROR)
!    CALL EXITS("GENERATE_STATUS_MASK")
!    RETURN 1

  END SUBROUTINE GENERATE_STATUS_MASK


  !
  !================================================================================================================================
  !

  !>Calculates and returns the MATRIX-VECTOR-prouct of the double precision VECTOR A*B in C.
  SUBROUTINE MATRIX_VECTOR_PRODUCT(A,B,C,Err)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(3,3) !<The A MATRIX
    REAL(DP), INTENT(IN) :: B(3) !<The B VECTOR
    REAL(DP), INTENT(OUT) :: C(3) !<On exit, the product VECTOR C=A*B
    INTEGER(INTG) :: ERR !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string

    !Local variables
        
!   CALL ENTERS("MATRIX_VECTOR_PRODUCT",Err)
    
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
!      CASE DEFAULT
!        CALL FLAG_ERROR("Invalid matrix size.",Err)
      END SELECT
!    ELSE
!      CALL FLAG_ERROR("Invalid matrix sizes.",Err)
    ENDIF

!    CALL EXITS("MATRIX_VECTOR_PRODUCT")
!    RETURN
!999 CALL ERRORS("MATRIX_VECTOR_PRODUCT",ERR,ERROR)
!    CALL EXITS("MATRIX_VECTOR_PRODUCT")
!    RETURN 1

  END SUBROUTINE MATRIX_VECTOR_PRODUCT

  !
  !================================================================================================================================
  !

  !>Calculates minimum and maximum value at array A.
  SUBROUTINE FIND_MINIMAX(A,N,MIN_VALUE,MAX_VALUE,Err)

    !Argument variables
    REAL(DP), ALLOCATABLE :: A(:)

    REAL(DP), INTENT(OUT) :: MIN_VALUE
    REAL(DP), INTENT(OUT) :: MAX_VALUE
    INTEGER(INTG), INTENT(IN)  :: N
    INTEGER(INTG) :: ERR !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I
        
!    CALL ENTERS("FIND_MINIMAX",Err)
    
    IF(SIZE(A,1).GT.2) THEN
      MIN_VALUE = A(1)
      MAX_VALUE = A(1)
      DO I=2,N
        IF (MIN_VALUE .GT. A(I)) THEN
          MIN_VALUE = A(I)
        ENDIF
        IF (MAX_VALUE .LT. A(I)) THEN
          MAX_VALUE = A(I)
        ENDIF
      ENDDO
!    ELSE
!      CALL FLAG_ERROR("Invalid matrix sizes.",Err)
    ENDIF

!    CALL EXITS("FIND_MINIMAX")
!    RETURN
!999 CALL ERRORS("FIND_MINIMAX",ERR,ERROR)
!    CALL EXITS("FIND_MINIMAX")
!    RETURN 1

  END SUBROUTINE FIND_MINIMAX

  !
  !================================================================================================================================
  !


  !>Calculates and returns the VECTOR-VECTOR-prouct of the double precision VECTOR A*B in C.
  SUBROUTINE VECTOR_VECTOR_PRODUCT(A,B,C,Err)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(3) !<The A VECTOR
    REAL(DP), INTENT(IN) :: B(3) !<The B VECTOR
    REAL(DP), INTENT(OUT) :: C !<On exit, the product SCALAR C=A*B
    INTEGER(INTG) :: ERR !<The error code
    !    TYPE(VARYING_STRING) :: LOCAL_ERROR !<The error string
    !Local variables
        
!    CALL ENTERS("VECTOR_VECTOR_PRODUCT",Err)
    
    IF(SIZE(A,1)==SIZE(B,1)) THEN
      SELECT CASE(SIZE(A,1))
        CASE(1)
          C=A(1)*B(1)
        CASE(2)
          C=A(1)*B(1)+A(2)*B(2)
        CASE(3)
          C=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
!        CASE DEFAULT
!          CALL FLAG_ERROR("Invalid matrix size.",Err)
      END SELECT
!    ELSE
!      CALL FLAG_ERROR("Invalid matrix sizes.",Err)
    ENDIF

!    CALL EXITS("VECTOR_VECTOR_PRODUCT")
!    RETURN
!999 CALL ERRORS("VECTOR_VECTOR_PRODUCT",ERR,ERROR)
!    CALL EXITS("VECTOR_VECTOR_PRODUCT")
!    RETURN 1

  END SUBROUTINE VECTOR_VECTOR_PRODUCT

  !
  !================================================================================================================================
  !


  !>to EXPORT output.
  SUBROUTINE POST_PROCESS_DATA(MATERIAL_BEHAVIOUR,OUTPUT_FILE_NAME,OUTPUT_FILE_FORMAT,TOTAL_NUMBER_OF_NODES,NODE_LIST,&
&CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,OUTPUT_FILE_FIELD_TITLE,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

    !subroutine variables
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)

    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES_PER_ELEMENT
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_NODES
    CHARACTER (LEN=300) :: OUTPUT_FILE_NAME
    CHARACTER (LEN=300) :: OUTPUT_FILE_FIELD_TITLE
    CHARACTER (LEN=10)  :: OUTPUT_FILE_FORMAT
    CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
    INTEGER(INTG) :: Err
    
    !Local variables
    CHARACTER (LEN=300) :: STRING
    INTEGER(INTG) :: TEXT_LENGTH
    LOGICAL :: EXIST
    INTEGER(INTG) :: I, J

   

! in the case the OUTPUT is in TABC format (2)
    IF (OUTPUT_FILE_FORMAT .EQ. "TABC") THEN 

!      EXPORT NODE TABLE list
      TEXT_LENGTH = INDEX(OUTPUT_FILE_NAME,' ') - 1
      STRING = OUTPUT_FILE_NAME(1:TEXT_LENGTH)//".tabc"
!  INQUIRE(FILE=STRING, EXIST=ex)
!      PRINT *, STRING
      OPEN (12,FILE=STRING)
!      OPEN (12,FILE=OUTPUT_FILE_NAME)

      WRITE(12,*)"VARIABLES=""NODE"",""X"",""Y"",""Z"",""U"",""V"",""W"",""Speed function along fibers"", &
 &                ""Speed function in transverse direction"",""Time"""
      WRITE(12,*)"zone i=",TOTAL_NUMBER_OF_NODES," , DATAPACKING=POINT"

      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,*) I,NODE_LIST(I,1),NODE_LIST(I,2),NODE_LIST(I,3),CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),&
            &CONDUCTIVITY_TENSOR(I,3),SPEED_FUNCTION_TABLE(I,1),&
            &SPEED_FUNCTION_TABLE(I,2),SPEED_FUNCTION_TABLE(I,3),&
            &SEED_VALUE(I)
      ENDDO
!      EXPORT NODE CONNECTIVITY list
      WRITE(12,*) "Connectivity"
      DO I=1,TOTAL_NUMBER_OF_NODES
        WRITE(12,*) I,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
!        PRINT *,STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
      ENDDO
      CLOSE(12)

    ENDIF



!   FILE="cmgui"
!   METHOD="FORTRAN"

!   EXPORT_FIELD=.TRUE.
!   IF(EXPORT_FIELD) THEN
!     WRITE(*,*)'Now export fields...'
!     CALL FLUID_MECHANICS_IO_WRITE_CMGUI(REGION,FILE,Err)
!     CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, Err)  
!     CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, Err)
!     WRITE(*,*)'All fields exported...'
!   ENDIF

!    CALL EXITS("GENERATE_STATUS_MASK")
!    RETURN
!999 CALL ERRORS("GENERATE_STATUS_MASK",ERR,ERROR)
!    CALL EXITS("GENERATE_STATUS_MASK")
!    RETURN 1

  END SUBROUTINE POST_PROCESS_DATA


  SUBROUTINE READ_TETGEN_MESH(INPUT_FILE_NAME,WORLD_REGION,MESH,GEOMETRIC_FIELD,ERR,ERROR,*)
    !subroutine parameters
    TYPE(VARYING_STRING), INTENT(IN) :: INPUT_FILE_NAME
    TYPE(REGION_TYPE), INTENT(IN), POINTER :: WORLD_REGION
    TYPE(MESH_TYPE), INTENT(INOUT), POINTER :: MESH
    TYPE(FIELD_TYPE), INTENT(INOUT), POINTER :: GEOMETRIC_FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    !CHARACTER (LEN=300) :: TEMP_STRING
    INTEGER(INTG) :: TEXT_LENGTH, NUMBER_OF_NODES, NUMBER_OF_ELEMENTS, NUMBER_OF_NODE_COMPONENTS, &
        & NUMBER_OF_ELEMENT_COMPONENTS, ELEMENT_ID, NODE_ID, NUMBER_OF_PROCESSORS, TEMP_INT, i, j
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    REAL(DP), ALLOCATABLE :: NODE_COORDINATES(:)
    TYPE(VARYING_STRING) :: INPUT_FILE_NAME_NODES
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
   
    CALL ENTERS("READ_TETGEN_MESH",ERR,ERROR,*999)
    
    OPEN (11,FILE=CHAR(INPUT_FILE_NAME//".node"))
    READ (11,*) NUMBER_OF_NODES, NUMBER_OF_NODE_COMPONENTS, TEMP_INT, TEMP_INT
    PRINT *,"Read Tetgen Mesh - #Nodes: ",NUMBER_OF_NODES
    PRINT *,"Read Tetgen Mesh - Dimension: ",NUMBER_OF_NODE_COMPONENTS
    
    OPEN (12,FILE=CHAR(INPUT_FILE_NAME//".ele"))
    READ (12,*) NUMBER_OF_ELEMENTS, NUMBER_OF_ELEMENT_COMPONENTS, TEMP_INT
    PRINT *,"Read Tetgen Mesh - #Elements: ",NUMBER_OF_ELEMENTS
    PRINT *,"Read Tetgen Mesh - Element components: ",NUMBER_OF_ELEMENT_COMPONENTS
    
    ! Create coordinate system
    NULLIFY(COORDINATE_SYSTEM)
    CALL COORDINATE_SYSTEM_CREATE_START(77000,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)
    
    ! Create region and assign coordinate system to it
    NULLIFY(REGION)
    CALL REGION_CREATE_START(77000,WORLD_REGION,REGION,ERR,ERROR,*999)
    CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
    CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)
    
    ! Create linear basis
    NULLIFY(BASIS)
    CALL BASIS_CREATE_START(77000,BASIS,ERR,ERROR,*999)
    CALL BASIS_TYPE_SET(BASIS,BASIS_SIMPLEX_TYPE,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_SET(BASIS,(/BASIS_LINEAR_SIMPLEX_INTERPOLATION,BASIS_LINEAR_SIMPLEX_INTERPOLATION, &
        & BASIS_LINEAR_SIMPLEX_INTERPOLATION/),ERR,ERROR,*999)
    CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
    
    !Create mesh
    NULLIFY(MESH)
    CALL MESH_CREATE_START(77000,REGION,NUMBER_OF_NODE_COMPONENTS,MESH,ERR,ERROR,*999)
    
    ! Create nodes
    NULLIFY(NODES)
    CALL NODES_CREATE_START(REGION,NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
    CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
    
    ! Create elements
    CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,1,ERR,ERROR,*999)
    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
    NULLIFY(ELEMENTS)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,1,BASIS,ELEMENTS,ERR,ERROR,*999)
    ALLOCATE(ELEMENT_NODES(NUMBER_OF_ELEMENT_COMPONENTS),STAT=ERR)
    DO i=1,NUMBER_OF_ELEMENTS
      READ (12,*) ELEMENT_ID,ELEMENT_NODES
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ELEMENT_ID,ELEMENTS,ELEMENT_NODES,Err,ERROR,*999)
    ENDDO
    DEALLOCATE(ELEMENT_NODES)
    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*999)
    
    CALL MESH_CREATE_FINISH(MESH,ERR,ERROR,*999)
    
    !Calculate decomposition
    NULLIFY(DECOMPOSITION)
    NUMBER_OF_PROCESSORS = COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
    CALL DECOMPOSITION_CREATE_START(77000,MESH,DECOMPOSITION,ERR,ERROR,*999)
    CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
    CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_PROCESSORS,ERR,ERROR,*999)
    CALL DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*999)

    !Create a field to put the geometry
    NULLIFY(GEOMETRIC_FIELD)
    CALL FIELD_CREATE_START(77000,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
    CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
    CALL FIELD_TYPE_SET(GEOMETRIC_FIELD,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_VARIABLES_SET(GEOMETRIC_FIELD,1,ERR,ERROR,*999)
    CALL FIELD_NUMBER_OF_COMPONENTS_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_NODE_COMPONENTS,Err,ERROR,*999)
    DO i=1,NUMBER_OF_NODE_COMPONENTS
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,i,1,ERR,ERROR,*999)
    ENDDO
    CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)

    !Set node positions
    ALLOCATE(NODE_COORDINATES(NUMBER_OF_NODE_COMPONENTS),STAT=ERR)
    DO i=1,NUMBER_OF_NODES
      READ (11,*) NODE_ID,NODE_COORDINATES
      DO j=1,NUMBER_OF_NODE_COMPONENTS
        CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,i,j, &
          & NODE_COORDINATES(j),ERR,ERROR,*999)
      ENDDO
    ENDDO
    DEALLOCATE(NODE_COORDINATES)
    
    CLOSE(11)
    CLOSE(12)
    
    CALL EXITS("READ_TETGEN_MESH")
    RETURN
999 CALL ERRORS("READ_TETGEN_MESH",ERR,ERROR)
    CALL EXITS("READ_TETGEN_MESH")
    RETURN 1
    
  END SUBROUTINE READ_TETGEN_MESH


  SUBROUTINE WRITE_VTK_MESH(OUTPUT_FILE_NAME,MESH,GEOMETRIC_FIELD,ERR,ERROR,*)
    !subroutine parameters
    TYPE(VARYING_STRING), INTENT(IN) :: OUTPUT_FILE_NAME
    TYPE(MESH_TYPE), INTENT(IN), POINTER :: MESH
    TYPE(FIELD_TYPE), INTENT(IN), POINTER :: GEOMETRIC_FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !local variables
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS, NUMBER_OF_NODES, node_idx, dim_idx, local_ny, ne, &
            & NUMBER_OF_NODES_PER_ELEMENT, i
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: GEOMETRIC_VARIABLE
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS
    
    CALL ENTERS("WRITE_VTK_MESH",ERR,ERROR,*999)
    
    OPEN (12,FILE=CHAR(OUTPUT_FILE_NAME // ".vtk"))
    
    WRITE(12,'(A)')"# vtk DataFile Version 3.0"
    WRITE(12,'(A)')"vtk output"
    WRITE(12,'(A)')"ASCII"
    WRITE(12,'(A)')"DATASET UNSTRUCTURED_GRID"
    
    CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
    NULLIFY(GEOMETRIC_VARIABLE)
    CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
            & ERR,ERROR,*999)
    
    NUMBER_OF_NODES=GEOMETRIC_VARIABLE%COMPONENTS(1)%DOMAIN%TOPOLOGY%NODES%NUMBER_OF_NODES
    WRITE(12,'(A,I8,A6)')"POINTS",NUMBER_OF_NODES,"float"
    
    DO node_idx=1,NUMBER_OF_NODES
      DO dim_idx=1,NUMBER_OF_DIMENSIONS
        local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
        WRITE(12,*) GEOMETRIC_PARAMETERS(local_ny)
      ENDDO
    ENDDO
    
    ELEMENTS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS
    NUMBER_OF_NODES_PER_ELEMENT=size(ELEMENTS%ELEMENTS(1)%ELEMENT_NODES,1)
    WRITE(12,'(A,I8,I8)')"CELLS ",ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,ELEMENTS% &
            & TOTAL_NUMBER_OF_ELEMENTS*(NUMBER_OF_NODES_PER_ELEMENT+1)
    DO ne=1,ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
      WRITE(12,*) NUMBER_OF_NODES_PER_ELEMENT, &
            & ((ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(i)-1),i=1,size(ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES,1))
    ENDDO
    
    WRITE(12,'(A,I8)')"CELL_TYPES",ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
    DO ne=1,ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
      WRITE(12,'(A)') "10"
    ENDDO
    
    WRITE(12,'(A,I8)')"CELL_DATA",ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
    WRITE(12,'(A,I8)')"POINT_DATA",NUMBER_OF_NODES
    
    !      export FIELD information
    !WRITE(12,'(A,A)')"FIELD number"," 1"
    !WRITE(12,'(A,I3,I8,A6)')OUTPUT_FILE_FIELD_TITLE,1,TOTAL_NUMBER_OF_NODES,"float"
    !DO I=1,TOTAL_NUMBER_OF_NODES
    !WRITE(12,'(F15.10)') SEED_VALUE(I)
    !ENDDO
    
    !      export VECTORS information
    !WRITE(12,'(A,A,A6)') "VECTORS ","fiber_vector","float"
    !DO I=1,TOTAL_NUMBER_OF_NODES
    !WRITE(12,'(3F8.5)') (CONDUCTIVITY_TENSOR(I,J),J=1,3)
    !ENDDO
    
    CLOSE(12)
    
    CALL EXITS("WRITE_VTK_MESH")
    RETURN
999 CALL ERRORS("WRITE_VTK_MESH",ERR,ERROR)
    CALL EXITS("WRITE_VTK_MESH")
    RETURN 1
    
  END SUBROUTINE WRITE_VTK_MESH
  
  
  
END MODULE HAMILTON_JACOBI_EQUATIONS_ROUTINES

