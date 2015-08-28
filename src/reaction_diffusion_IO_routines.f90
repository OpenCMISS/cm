!> \file
!> \author Vijay Rajagopal
!> \brief This module handles some mesh/parameter input routines and cmgui output routines for reaction diffusion
!> routines and should be eventually replaces by field_IO_routines.f90 
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

!> Temporary IO routines for fluid mechanics

MODULE REACTION_DIFFUSION_IO_ROUTINES

 USE BASE_ROUTINES
 USE COMP_ENVIRONMENT
 USE EQUATIONS_SET_CONSTANTS
 USE FIELD_ROUTINES
 USE TYPES
 USE INPUT_OUTPUT 
 USE KINDS   
 USE MESH_ROUTINES
 USE MPI

#include "macros.h"  

  IMPLICIT NONE

  PUBLIC REACTION_DIFFUSION_IO_WRITE_CMGUI

CONTAINS

  ! OK
  !================================================================================================================================
  !

  !> Writes solution into cmgui formats exelem and exnode.
  SUBROUTINE REACTION_DIFFUSION_IO_WRITE_CMGUI(REGION, EQUATIONS_SET_GLOBAL_NUMBER, NAME, ERR, ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), INTENT(IN), POINTER :: REGION !<A pointer to the region to get the coordinate system for
    CHARACTER(28), INTENT(IN) :: NAME !<the prefix name of file.
    INTEGER(INTG) :: ERR !<The error code
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_GLOBAL_NUMBER
    TYPE(VARYING_STRING):: ERROR !<The error string

    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(DOMAIN_TYPE), POINTER :: COMPUTATIONAL_DOMAIN
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD
    REAL(DP) :: NodeXValue,NodeYValue,NodeZValue,NodeUValue
    INTEGER(INTG):: myComputationalNodeNumber,NumberOfOutputFields,NumberOfDimensions,NumberOfElements,NumberOfNodes
    INTEGER(INTG):: NumberOfVariableComponents,NumberOfSourceComponents,I,J,K,ValueIndex,NODE_GLOBAL_NUMBER
    INTEGER(INTG) :: NodesInMeshComponent,BasisType,MaxNodesPerElement,NumberOfFieldComponents(3),ELEMENT_GLOBAL_NUMBER
    INTEGER(INTG) :: NODE_LOCAL_NUMBER
    INTEGER(INTG),ALLOCATABLE :: ElementNodes(:,:),SimplexOutputHelp(:)
    REAL(DP), ALLOCATABLE :: ElementNodesScales(:,:)
    LOGICAL :: OUTPUT_SOURCE
    TYPE(VARYING_STRING) :: FILENAME !<the prefix name of file.
    CHARACTER(50) :: INTG_STRING,INTG_STRING2


    ENTERS("REACTION_DIFFUSION_IO_WRITE_CMGUI",ERR,ERROR,*999)

    myComputationalNodeNumber = COMPUTATIONAL_NODE_NUMBER_GET(err,error)

    EQUATIONS_SET => REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr
    NULLIFY(SOURCE_FIELD)
    COMPUTATIONAL_DOMAIN=>REGION%MESHES%MESHES(1) & 
      & %ptr%DECOMPOSITIONS%DECOMPOSITIONS(1)%ptr%DOMAIN(1)%ptr

    myComputationalNodeNumber = COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
    NumberOfDimensions = COMPUTATIONAL_DOMAIN%NUMBER_OF_DIMENSIONS
    NumberOfNodes = COMPUTATIONAL_DOMAIN%TOPOLOGY%NODES%NUMBER_OF_NODES
    NodesInMeshComponent = REGION%meshes%meshes(1)%ptr%topology(1)%ptr%nodes%numberOfNodes
    NumberOfElements = COMPUTATIONAL_DOMAIN%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
    NumberOfVariableComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
      & variables(1)%number_of_components
    NumberOfOutputFields=2
    !determine if there is a source field
    OUTPUT_SOURCE = .FALSE.
    IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) & 
      & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
        & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
          SOURCE_FIELD=>REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field
          IF( ASSOCIATED(SOURCE_FIELD) ) OUTPUT_SOURCE = .FALSE. !currently set to false to rethink how source is accessed for output
     END IF

    IF( OUTPUT_SOURCE ) THEN
      NumberOfSourceComponents=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%source%source_field% &
        & variables(1)%number_of_components
      NumberOfOutputFields = NumberOfOutputFields + 1
    !  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(SOURCE_FIELD,SOURCE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
    !  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(SOURCE_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATED_POINT,ERR,ERROR,*999)

    END IF

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Nodes...",ERR,ERROR,*999)

    FILENAME="./output/"//NAME//".exnode"

    OPEN(UNIT=myComputationalNodeNumber, FILE=CHAR(FILENAME),STATUS='unknown')
    ! WRITING HEADER INFORMATION
    WRITE(myComputationalNodeNumber,*) 'Group name: Cell'
    WRITE(INTG_STRING,'(I0)'),NumberOfOutputFields 
    WRITE(myComputationalNodeNumber,*) '#Fields=',TRIM(INTG_STRING)

    ValueIndex=1
    WRITE(INTG_STRING,'(I0)'),NumberOfDimensions 
    WRITE(myComputationalNodeNumber,*) &
      & ' 1) coordinates,  coordinate, rectangular cartesian, #Components=',TRIM(INTG_STRING)
    DO I=1,NumberOfDimensions
      IF(I==1) THEN
        WRITE(INTG_STRING,'(I0)'),ValueIndex 
        WRITE(myComputationalNodeNumber,*) '   x.  Value index= ',TRIM(INTG_STRING),', #Derivatives= 0'
      ELSE IF(I==2) THEN
        WRITE(INTG_STRING,'(I0)'),ValueIndex 
        WRITE(myComputationalNodeNumber,*) '   y.  Value index= ',TRIM(INTG_STRING),', #Derivatives= 0'
      ELSE
        WRITE(INTG_STRING,'(I0)'),ValueIndex 
        WRITE(myComputationalNodeNumber,*) '   z.  Value index= ',TRIM(INTG_STRING),', #Derivatives= 0'
      END IF
      ValueIndex=ValueIndex+1
    END DO

    WRITE(INTG_STRING,'(I0)'),NumberOfVariableComponents 
    WRITE(myComputationalNodeNumber,*) ' 2) dependent, field, rectangular cartesian, #Components=', & 
      & TRIM(INTG_STRING)

    DO I=1,NumberOfVariableComponents
      WRITE(INTG_STRING,'(I0)'),ValueIndex
      WRITE(INTG_STRING2,'(I0)'),I 
      WRITE(myComputationalNodeNumber,*)  '  ',TRIM(INTG_STRING2),'. Value index= ',TRIM(INTG_STRING), & 
        & ', #Derivatives= 0' 
      ValueIndex=ValueIndex+1
    END DO

    IF( OUTPUT_SOURCE ) THEN !Watch out that no numbering conflict occurs with Analytic: 4.)
      WRITE(INTG_STRING,'(I0)'),NumberOfSourceComponents 
      WRITE(myComputationalNodeNumber,*) ' 3) source, field, rectangular cartesian, #Components=', & 
        & TRIM(INTG_STRING)
      DO I=1,NumberOfSourceComponents
        WRITE(INTG_STRING,'(I0)'),ValueIndex
        WRITE(INTG_STRING2,'(I0)'),I 
        WRITE(myComputationalNodeNumber,*)  '   ',TRIM(INTG_STRING2),'.  Value index= ', & 
          & TRIM(INTG_STRING),', #Derivatives= 0' 
        ValueIndex=ValueIndex+1
      END DO
    END IF

    !WRITE OUT NODE VALUES
    DO I = 1,NumberOfNodes
      NODE_GLOBAL_NUMBER = COMPUTATIONAL_DOMAIN%TOPOLOGY%NODES%NODES(I)%GLOBAL_NUMBER
      NodeXValue = REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
        & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(I)
      IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
        NodeYValue = REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(I+NumberOfNodes)
      ENDIF
      IF(NumberOfDimensions==3) THEN
        NodeZValue = REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%geometry%geometric_field%variables(1) &
          & %parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(I+(2*NumberOfNodes))
      ENDIF
      NodeUValue=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%dependent%dependent_field% &
        & variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(I)

      WRITE(myComputationalNodeNumber,*) ' Node: ',NODE_GLOBAL_NUMBER
      WRITE(myComputationalNodeNumber,'("    ", es25.16 )')NodeXValue

      IF(NumberOfDimensions==2 .OR. NumberOfDimensions==3) THEN
        WRITE(myComputationalNodeNumber,'("    ", es25.16 )')NodeYValue
      END IF

      IF(NumberOfDimensions==3) THEN
        WRITE(myComputationalNodeNumber,'("    ", es25.16 )')NodeZValue
      END IF
      WRITE(myComputationalNodeNumber,'("    ", es25.16 )')NodeUValue

      IF( (EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS) & 
        & .AND.(EQUATIONS_SET%TYPE==EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE) &
          & .AND.(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_REAC_DIFF_SUBTYPE) )THEN
          !source field
          IF( OUTPUT_SOURCE ) THEN
            !NodeSourceValue = SOURCE_INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
            !WRITE(myComputationalNodeNumber,'("    ", es25.16 )')NodeSourceValue
          END IF
      END IF
    END DO !nodes I
    CLOSE(myComputationalNodeNumber)

    !OUTPUT ELEMENTS IN CURRENT DOMAIN
    MaxNodesPerElement=COMPUTATIONAL_DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(1)%basis%number_of_element_parameters
    BasisType = 1
    IF(NumberOfDimensions==2) THEN
      IF(MaxNodesPerElement==4.OR.MaxNodesPerElement==9.OR.MaxNodesPerElement==16) THEN
        BasisType=1
      ELSEIF(MaxNodesPerElement==3.OR.MaxNodesPerElement==6.OR.MaxNodesPerElement==10) THEN
        BasisType=2
      ENDIF
    ELSEIF(NumberOfDimensions==3) THEN
      BasisType=REGION%equations_sets%equations_sets(EQUATIONS_SET_GLOBAL_NUMBER)%ptr%equations% &
        & interpolation%geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr%bases(1)%ptr%type
    ENDIF
    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Writing Elements...",ERR,ERROR,*999)
    FILENAME="./output/"//NAME//".exelem"
    OPEN(UNIT=myComputationalNodeNumber, FILE=CHAR(FILENAME),STATUS='unknown')
    WRITE(myComputationalNodeNumber,*) 'Group name: Cell'
    IF (BasisType==1) THEN !lagrange basis in 1 and 2D
      WRITE(INTG_STRING,'(I0)'),NumberOfDimensions
      WRITE(myComputationalNodeNumber,*) 'Shape.  Dimension= ',TRIM(INTG_STRING)
      WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 1'
      IF(NumberOfDimensions==1) THEN
        WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
        WRITE(myComputationalNodeNumber,*) 'q.Lagrange, #Scale factors=',TRIM(INTG_STRING)
      ELSE IF (NumberOfDimensions==2) THEN
        IF(MaxNodesPerElement==4) THEN

          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & 'l.Lagrange*l.Lagrange, #Scale factors=',TRIM(INTG_STRING) !linear lagrange
        ELSE IF(MaxNodesPerElement==9) THEN
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & 'q.Lagrange*q.Lagrange, #Scale factors=',TRIM(INTG_STRING) !quadratic lagrange
        ELSE IF(MaxNodesPerElement==16) THEN
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & 'c.Lagrange*c.Lagrange, #Scale factors=',TRIM(INTG_STRING) !cubic lagrange
        END IF
      ELSE !three dimensions
        IF(MaxNodesPerElement==8) THEN
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & 'l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=',TRIM(INTG_STRING)
        ELSE IF(MaxNodesPerElement==27) THEN
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & 'q.Lagrange*q.Lagrange*q.Lagrange, #Scale factors=',TRIM(INTG_STRING)
        ELSE IF(MaxNodesPerElement==64) THEN
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & 'c.Lagrange*c.Lagrange*c.Lagrange, #Scale factors=',TRIM(INTG_STRING)
        END IF
      END IF
    ELSEIF(BasisType==2) THEN
      IF(NumberOfDimensions==2) THEN
        WRITE(myComputationalNodeNumber,*) 'Shape.  Dimension=', &
          & NumberOfDimensions,', simplex(2)*simplex'
        IF(MaxNodesPerElement==3) THEN
          WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 1'
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*)  & 
            & ' l.simplex(2)*l.simplex, #Scale factors= ', TRIM(INTG_STRING)
        ELSE IF(MaxNodesPerElement==6) THEN
          WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 1'
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) & 
            & ' l.simplex(2)*l.simplex, #Scale factors= ', TRIM(INTG_STRING)
        ELSE IF (MaxNodesPerElement== 10 ) THEN
          WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 1'
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & ' q.simplex(2)*q.simplex, #Scale factors= ', TRIM(INTG_STRING)
        ENDIF
      ELSE IF(NumberOfDimensions==3) THEN
        WRITE(INTG_STRING2,'(I0)'),NumberOfDimensions
        WRITE(myComputationalNodeNumber,*) &
          & 'Shape.  Dimension=',TRIM(INTG_STRING2),', simplex(2;3)*simplex*simplex'
        IF(MaxNodesPerElement==4) THEN
          WRITE(myComputationalNodeNumber,*) &
            & '#Scale factor sets= 1'
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & ' l.simplex(2;3)*l.simplex*l.simplex, #Scale factors= ', TRIM(INTG_STRING)
        ELSE IF (MaxNodesPerElement== 10 ) THEN
          WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 1'
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', TRIM(INTG_STRING)
        ELSE IF(MaxNodesPerElement==20) THEN
          WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 1'
          WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
          WRITE(myComputationalNodeNumber,*) &
            & ' q.simplex(2;3)*q.simplex*q.simplex, #Scale factors= ', TRIM(INTG_STRING)
        ENDIF      
      ELSE
        WRITE(myComputationalNodeNumber,*) '#Scale factor sets= 0'
      END IF

    END IF
    WRITE(INTG_STRING,'(I0)'),MaxNodesPerElement
    WRITE(myComputationalNodeNumber,*) '#Nodes= ',TRIM(INTG_STRING)
    WRITE(INTG_STRING,'(I0)'),NumberOfOutputFields
    WRITE(myComputationalNodeNumber,*) '#Fields= ',TRIM(INTG_STRING)
    NumberOfFieldComponents(1) = NumberOfDimensions
    NumberOfFieldComponents(2) = NumberOfVariableComponents
    NumberOfFieldComponents(3) = NumberOfSourceComponents
    DO I=1,NumberOfOutputFields
      IF(I==1)THEN
        WRITE(INTG_STRING,'(I0)'),NumberOfDimensions
        WRITE(myComputationalNodeNumber,*) & 
          & ' 1) coordinates,  coordinate, rectangular cartesian, #Components= ',TRIM(INTG_STRING)
      ELSE IF(I==2) THEN
        WRITE(INTG_STRING,'(I0)'),NumberOfVariableComponents
        WRITE(myComputationalNodeNumber,*) &
        & ' 2) dependent,  field,  rectangular cartesian, #Components= ',TRIM(INTG_STRING)
      ELSE IF(I==3) THEN
        WRITE(INTG_STRING,'(I0)'),NumberOfSourceComponents
        WRITE(myComputationalNodeNumber,*) &
          & ' 3) source,  field,  rectangular cartesian, #Components= ',TRIM(INTG_STRING)
      END IF

      DO J=1,NumberOfFieldComponents(I)
        IF(NumberOfDimensions==1) THEN
          IF(I==1)THEN
            IF(J==1) THEN
                WRITE(myComputationalNodeNumber,*)'   x.   l.Lagrange, no modify, standard node based.'
            ELSE IF(J==2) THEN
                WRITE(myComputationalNodeNumber,*)'   y.   l.Lagrange, no modify, standard node based.'
            ELSE IF(J==3) THEN
                WRITE(myComputationalNodeNumber,*)'   z.   l.Lagrange, no modify, standard node based.'
            END IF
          ELSE
            WRITE(myComputationalNodeNumber,*) &
              & '   ',J,'.   l.Lagrange, no modify, standard node based.'
          END IF
        ELSE IF(NumberOfDimensions==2) THEN
          IF(I==1)THEN
            IF(J==1) THEN
              IF(MaxNodesPerElement==4)THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9) THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.   c.Lagrange*c.Lagrange, no modify, standard node based.'

              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   x.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF 
            ELSE IF(J==2) THEN
              IF(MaxNodesPerElement==4) THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   y.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   y.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   y.   c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   y.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   y.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                & '   y.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF
            ELSE IF(J==3) THEN
              IF(MaxNodesPerElement==4) THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   z.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   z.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   z.   c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   z.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   z.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   z.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF
            END IF
          ELSE
              IF(MaxNodesPerElement==4) THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.   l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==9)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   ',J,'.   q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==16)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   ',J,'.   c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==3)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   ',J,'.  l.simplex(2)*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==6)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   ',J,'.  q.simplex(2)*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) & 
                  & '   ',J,'.  c.simplex(2)*c.simplex, no modify, standard node based.'
              END IF
          END IF
        ELSE IF(NumberOfDimensions==3) THEN
          IF(I==1)THEN
            IF(J==1) THEN
              IF(MaxNodesPerElement==8) THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   x.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF 
            ELSE IF(J==2) THEN
              IF(MaxNodesPerElement==8) THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   y.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   y.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   y.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   y.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   y.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF
            ELSE IF(J==3) THEN
              IF(MaxNodesPerElement==8) THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   z.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   z.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   z.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   z.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   z.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF
            END IF
          ELSE
              IF(MaxNodesPerElement==8) THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.   l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==27)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.   q.Lagrange*q.Lagrange*q.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==64)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.   c.Lagrange*c.Lagrange*c.Lagrange, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==4)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==10)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.'
              ELSE IF(MaxNodesPerElement==20)  THEN
                WRITE(myComputationalNodeNumber,*) &
                  & '   ',J,'.  c.simplex(2;3)*c.simplex*c.simplex, no modify, standard node based.'
              END IF
          END IF
        END IF
        WRITE(INTG_STRING,'(I0)') MaxNodesPerElement
        WRITE(myComputationalNodeNumber,*) '   #Nodes= ',TRIM(INTG_STRING)
 
        DO K = 1,MaxNodesPerElement
          WRITE(INTG_STRING,'(I0)'),K
          WRITE(myComputationalNodeNumber,*) '    ',TRIM(INTG_STRING),'.  #Values=1'
          WRITE(myComputationalNodeNumber,*) '     Value indices:     1'
          WRITE(myComputationalNodeNumber,*) '     Scale factor indices:   ',TRIM(INTG_STRING)
        END DO
      END DO !J loop
    END DO !I loop
    IF(.NOT.ALLOCATED(ElementNodes)) ALLOCATE(ElementNodes(NumberOfElements,MaxNodesPerElement))
    IF(.NOT.ALLOCATED(ElementNodesScales)) ALLOCATE(ElementNodesScales(NumberOfElements,MaxNodesPerElement))
    DO I=1,NumberOfElements
      ELEMENT_GLOBAL_NUMBER=COMPUTATIONAL_DOMAIN%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(K)%GLOBAL_NUMBER 
      DO J=1,MaxNodesPerElement
        NODE_LOCAL_NUMBER=COMPUTATIONAL_DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(I)%ELEMENT_NODES(J)
        NODE_GLOBAL_NUMBER=COMPUTATIONAL_DOMAIN%MAPPINGS%NODES%LOCAL_TO_GLOBAL_MAP(NODE_LOCAL_NUMBER)
        ElementNodes(I,J)=NODE_GLOBAL_NUMBER
        ElementNodesScales(I,J)=1.0000000000000000E+00
      END DO
    END DO

    
    DO K=1,NumberOfElements
      ELEMENT_GLOBAL_NUMBER=COMPUTATIONAL_DOMAIN%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(K)%GLOBAL_NUMBER 
      IF (BasisType==1) THEN
        WRITE(INTG_STRING,'(I0)'),ELEMENT_GLOBAL_NUMBER
        WRITE(myComputationalNodeNumber,*) 'Element:     ', TRIM(INTG_STRING),' 0  0'
        WRITE(myComputationalNodeNumber,*) '   Nodes:'
        WRITE(myComputationalNodeNumber,*) '   ', ElementNodes(K,1:MaxNodesPerElement)
        WRITE(myComputationalNodeNumber,*) '   Scale factors:'
        WRITE(myComputationalNodeNumber,*) '   ',ElementNodesScales(K,1:MaxNodesPerElement)

      ELSEIF(BasisType==2) THEN
        IF(.NOT.ALLOCATED(SimplexOutputHelp)) ALLOCATE(SimplexOutputHelp(MaxNodesPerElement))
        IF(NumberOfDimensions==2)THEN
          SimplexOutputHelp(1)=ElementNodes(K,1)
          SimplexOutputHelp(2)=ElementNodes(K,2)
          SimplexOutputHelp(3)=ElementNodes(K,3)
        ELSE IF(NumberOfDimensions==3) THEN
          SimplexOutputHelp(1)=ElementNodes(K,1)
          SimplexOutputHelp(2)=ElementNodes(K,4)
          SimplexOutputHelp(3)=ElementNodes(K,2)
          SimplexOutputHelp(4)=ElementNodes(K,3)
        END IF
        WRITE(INTG_STRING,'(I0)') ELEMENT_GLOBAL_NUMBER
        WRITE(myComputationalNodeNumber,*) 'Element:     ', TRIM(INTG_STRING),' 0  0'
        WRITE(myComputationalNodeNumber,*) '   Nodes:'
        WRITE(myComputationalNodeNumber,*) '   ', SimplexOutputHelp
        WRITE(myComputationalNodeNumber,*) '   Scale factors:'
        WRITE(myComputationalNodeNumber,*) '   ',ElementNodesScales(K,1:MaxNodesPerElement)
      END IF
    ENDDO
    CLOSE(myComputationalNodeNumber)

    EXITS("REACTION_DIFFUSION_IO_WRITE_CMGUI")
    RETURN     
999 ERRORSEXITS("REACTION_DIFFUSION_IO_WRITE_CMGUI",ERR,ERROR)    
    RETURN 1

  END SUBROUTINE REACTION_DIFFUSION_IO_WRITE_CMGUI

  !================================================================================================================================
  !


END MODULE REACTION_DIFFUSION_IO_ROUTINES
