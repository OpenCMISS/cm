!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all node routines.
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
!> The Original Code is openCMISS
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

!> This module handles all node routines.
MODULE NODE_ROUTINES

  USE BASE_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TREES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC NODE_CHECK_EXISTS,NODE_DESTROY,NODE_INITIAL_POSITION_SET,NODE_NUMBER_SET,NODES_INITIALISE,NODES_FINALISE, &
    & NODES_CREATE_START,NODES_CREATE_FINISH

CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE NODE_CHECK_EXISTS(USER_NUMBER,REGION,NODE_EXISTS,GLOBAL_NUMBER,ERR,ERROR,*)

    !#### Subroutine:_NODE_CHECK_EXISTS
    !###  Description:
    !###    Checks that a user node number is defined on the specified region. If the node is defined then NODE_EXISTS is set to
    !###    .TRUE. and GLOBAL_NUMBER is set to the global node number in that region. If the node is not defined then NODE_EXISTS
    !###    is set to .FALSE. and GLOBAL_NUMBER is set to 0.

    !Argument variables
    INTEGER(INTG) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    LOGICAL, INTENT(OUT) :: NODE_EXISTS
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    GLOBAL_NUMBER=0
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(REGION%NODES%NODE_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(REGION%NODES%NODE_TREE,TREE_NODE,GLOBAL_NUMBER,ERR,ERROR,*999)
          NODE_EXISTS=.TRUE.
        ENDIF
      ELSE
        LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " does not have any associated nodes"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("NODE_CHECK_EXISTS")
    RETURN
999 CALL ERRORS("NODE_CHECK_EXISTS",ERR,ERROR)
    CALL EXITS("NODE_CHECK_EXISTS")
    RETURN 1   
  END SUBROUTINE NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  SUBROUTINE NODE_DESTROY(NODE,ERR,ERROR,*)

    !#### Subroutine: NODE_DESTROY
    !###  Description:
    !###    Destroys the a nodes and deallocates all memory

    !Argument variables
    TYPE(NODE_TYPE) :: NODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("NODE_DESTROY",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%INITIAL_POSITION)) DEALLOCATE(NODE%INITIAL_POSITION)

    CALL EXITS("NODE_DESTROY")
    RETURN
999 CALL ERRORS("NODE_DESTROY",ERR,ERROR)
    CALL EXITS("NODE_DESTROY")
    RETURN 1   
  END SUBROUTINE NODE_DESTROY
  
  !
  !================================================================================================================================
  !

  SUBROUTINE NODE_INITIAL_POSITION_SET(GLOBAL_NUMBER,INITIAL_POSITION,NODES,ERR,ERROR,*)

    !#### Subroutine: NODE_INITIAL_POSITION_SET
    !###  Description:
    !###    Changes/sets the initial position of a node identified by a given global number.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER
    REAL(DP), INTENT(IN) :: INITIAL_POSITION(:)
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("NODE_INITIAL_POSITION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        CALL FLAG_ERROR("Nodes have been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(NODES%REGION)) THEN
          IF(ASSOCIATED(NODES%REGION%COORDINATE_SYSTEM)) THEN
            IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
              NUMBER_OF_DIMENSIONS=NODES%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              IF(SIZE(INITIAL_POSITION,1)==NUMBER_OF_DIMENSIONS) THEN
                NODES%NODES(GLOBAL_NUMBER)%INITIAL_POSITION=INITIAL_POSITION
              ELSE
                LOCAL_ERROR="The number of dimensions for the position ("// &
                  & TRIM(NUMBER_TO_VSTRING(SIZE(INITIAL_POSITION,1),"*",ERR,ERROR))// &
                  & ") do not match the number of dimensions in the region ("// &
                  & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999) 
              ENDIF
            ELSE
              LOCAL_ERROR="Global node number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Nodes region coordinate system is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Nodes region is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Nodes is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("NODE_INITIAL_POSITION_SET")
    RETURN
999 CALL ERRORS("NODE_INITIAL_POSITION_SET",ERR,ERROR)    
    CALL EXITS("NODE_INITIAL_POSITION_SET")
    RETURN 1
   
  END SUBROUTINE NODE_INITIAL_POSITION_SET
        
  !
  !================================================================================================================================
  !

  SUBROUTINE NODE_NUMBER_SET(GLOBAL_NUMBER,USER_NUMBER,NODES,ERR,ERROR,*)

    !#### Subroutine: NODE_NUMBER_SET
    !###  Description:
    !###    Changes/sets the user number for a node identified by a given global number.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER,USER_NUMBER
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: INSERT_STATUS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("NODE_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        CALL FLAG_ERROR("Nodes have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          CALL TREE_ITEM_DELETE(NODES%NODE_TREE,NODES%NODES(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
          CALL TREE_ITEM_INSERT(NODES%NODE_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
          IF(INSERT_STATUS/=TREE_NODE_INSERT_SUCESSFUL) CALL FLAG_ERROR("Unsucessful tree insert",ERR,ERROR,*999)
          CALL TREE_ITEM_DELETE(NODES%NODE_TREE,NODES%NODES(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
          CALL TREE_ITEM_INSERT(NODES%NODE_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
          IF(INSERT_STATUS/=TREE_NODE_INSERT_SUCESSFUL) CALL FLAG_ERROR("Unsucessful tree insert",ERR,ERROR,*999)
          NODES%NODES(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
        ELSE
          LOCAL_ERROR="Global node number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Nodes is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("NODE_NUMBER_SET")
    RETURN
999 CALL ERRORS("NODE_NUMBER_SET",ERR,ERROR)    
    CALL EXITS("NODE_NUMBER_SET")
    RETURN 1
   
  END SUBROUTINE NODE_NUMBER_SET
        
  !
  !================================================================================================================================
  !

  SUBROUTINE NODES_CREATE_FINISH(REGION,ERR,ERROR,*)

    !#### Subroutine: NODES_CREATE_FINISH
    !###  Description:
    !###    Finishes the process of creating nodes in the region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: np,NUMBER_OF_DIMENSIONS
    
    CALL ENTERS("NODES_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        REGION%NODES%NODES_FINISHED=.TRUE.
      ELSE
        CALL FLAG_ERROR("Region nodes is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",REGION%NODES%NUMBER_OF_NODES,ERR,ERROR,*999)
      NUMBER_OF_DIMENSIONS=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
      DO np=1,REGION%NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number    = ",REGION%NODES%NODES(np)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number      = ",REGION%NODES%NODES(np)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Label            = ",REGION%NODES%NODES(np)%LABEL, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NUMBER_OF_DIMENSIONS,3,3, &
          & REGION%NODES%NODES(np)%INITIAL_POSITION,'("    Initial Position =",3(X,E13.6))','(18X,3(X,E13.6))',ERR,ERROR,*999)
      ENDDO !np
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"User->Global number tree",ERR,ERROR,*999)
      CALL TREE_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,REGION%NODES%NODE_TREE,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("NODES_CREATE_FINISH")
    RETURN
999 CALL ERRORS("NODES_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("NODES_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE NODES_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  SUBROUTINE NODES_CREATE_START(NUMBER_OF_NODES,REGION,NODES,ERR,ERROR,*)

    !#### Subroutine: NODES_CREATE_START
    !###  Description:
    !###    Starts the process of creating nodes in the region identified by REGION. NUMBER_OF_NODES is the number of nodes
    !###    that will be created and NODES is the pointer to the NODES data structure.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(NODES_TYPE), POINTER :: NODES
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: INSERT_STATUS,np,NUMBER_OF_DIMENSIONS
    TYPE(NODES_TYPE), POINTER :: NEW_NODES

    NULLIFY(NEW_NODES)
    
    CALL ENTERS("NODES_CREATE_START",ERR,ERROR,*999)

    NULLIFY(NODES)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        CALL FLAG_ERROR("Region already has nodes associated",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
          NUMBER_OF_DIMENSIONS=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
          ALLOCATE(NEW_NODES,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh nodes",ERR,ERROR,*999)
          ALLOCATE(NEW_NODES%NODES(NUMBER_OF_NODES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate individual nodes",ERR,ERROR,*999)
          NEW_NODES%NUMBER_OF_NODES=NUMBER_OF_NODES
          NEW_NODES%NODES_FINISHED=.FALSE.
          NULLIFY(NEW_NODES%NODE_TREE)
          CALL TREE_CREATE_START(NEW_NODES%NODE_TREE,ERR,ERROR,*999)
          CALL TREE_INSERT_TYPE_SET(NEW_NODES%NODE_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
          CALL TREE_CREATE_FINISH(NEW_NODES%NODE_TREE,ERR,ERROR,*999)
          !Set default node numbers
          DO np=1,NEW_NODES%NUMBER_OF_NODES
            NEW_NODES%NODES(np)%GLOBAL_NUMBER=np
            NEW_NODES%NODES(np)%USER_NUMBER=np
            NEW_NODES%NODES(np)%LABEL=""
            CALL TREE_ITEM_INSERT(NEW_NODES%NODE_TREE,np,np,INSERT_STATUS,ERR,ERROR,*999)
            !!TODO:: Make initial position optional and not allocate etc. ???
            ALLOCATE(NEW_NODES%NODES(np)%INITIAL_POSITION(NUMBER_OF_DIMENSIONS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate initial position",ERR,ERROR,*999)
            NEW_NODES%NODES(np)%INITIAL_POSITION=0.0_DP
          ENDDO !np
          REGION%NODES=>NEW_NODES
          NEW_NODES%REGION=>REGION
          NODES=>NEW_NODES
        ELSE
          CALL FLAG_ERROR("Region coordinate system is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("NODES_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_NODES)) THEN
      IF(ASSOCIATED(NEW_NODES%NODES)) THEN
        DO np=1,NEW_NODES%NUMBER_OF_NODES
          IF(ALLOCATED(NEW_NODES%NODES(np)%INITIAL_POSITION)) DEALLOCATE(NEW_NODES%NODES(np)%INITIAL_POSITION)
        ENDDO !np
        DEALLOCATE(NEW_NODES%NODES)
      ENDIF
      DEALLOCATE(NEW_NODES)
    ENDIF
    NULLIFY(NODES)
    CALL ERRORS("NODES_CREATE_START",ERR,ERROR)
    CALL EXITS("NODES_CREATE_START")
    RETURN 1
   
  END SUBROUTINE NODES_CREATE_START

  !
  !===============================================================================================================================
  !

  SUBROUTINE NODES_FINALISE(REGION,ERR,ERROR,*)

    !#### Subroutine: NODES_FINALISE
    !###  Description:
    !###    Finalises the nodes data structures for a region and deallocates any memory. 

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: np

    CALL ENTERS("NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        DO np=1,REGION%NODES%NUMBER_OF_NODES
          CALL NODE_DESTROY(REGION%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        DEALLOCATE(REGION%NODES%NODES)
        IF(ASSOCIATED(REGION%NODES%NODE_TREE)) CALL TREE_DESTROY(REGION%NODES%NODE_TREE,ERR,ERROR,*999)
        DEALLOCATE(REGION%NODES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("NODES_FINALISE")
    RETURN
999 CALL ERRORS("NODES_FINALISE",ERR,ERROR)
    CALL EXITS("NODES_FINALISE")
    RETURN 1
  END SUBROUTINE NODES_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE NODES_INITIALISE(REGION,ERR,ERROR,*)

    !#### Subroutine: NODES_INITIALISE
    !###  Description:
    !###    Initialises the nodes in a given region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        CALL FLAG_ERROR("Region has associated nodes",ERR,ERROR,*999)
      ELSE
        NULLIFY(REGION%NODES)
        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("NODES_INITIALISE")
    RETURN
999 CALL ERRORS("NODES_INITIALISE",ERR,ERROR)
    CALL EXITS("NODES_INITIALISE")
    RETURN 1
  END SUBROUTINE NODES_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE NODE_ROUTINES
