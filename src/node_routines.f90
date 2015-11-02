!> \file
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

!> This module handles all node routines.
MODULE NODE_ROUTINES

  USE BASE_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TREES
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  !>Starts the process of creating nodes for an interface or region
  INTERFACE NODES_CREATE_START
    MODULE PROCEDURE NODES_CREATE_START_REGION
    MODULE PROCEDURE NODES_CREATE_START_INTERFACE
  END INTERFACE !NODES_CREATE_START

  !>Initialises nodes for an interface or region
  INTERFACE NODES_INITIALISE
    MODULE PROCEDURE NODES_INITIALISE_REGION
    MODULE PROCEDURE NODES_INITIALISE_INTERFACE
  END INTERFACE !NODES_INITIALIES

  !>Gets the label for a node identified by a given global number.
  INTERFACE NODES_LABEL_GET
    MODULE PROCEDURE NODES_LABEL_GET_C
    MODULE PROCEDURE NODES_LABEL_GET_VS
  END INTERFACE !NODES_LABEL_SET

  !>Changes/sets the label for a node identified by a given global number.
  INTERFACE NODES_LABEL_SET
    MODULE PROCEDURE NODES_LABEL_SET_C
    MODULE PROCEDURE NODES_LABEL_SET_VS
  END INTERFACE !NODES_LABEL_SET

  PUBLIC NODE_CHECK_EXISTS

  PUBLIC NODES_CREATE_FINISH,NODES_CREATE_START,NODES_DESTROY

  PUBLIC NODES_LABEL_GET,NODES_LABEL_SET

  PUBLIC NODES_NUMBER_OF_NODES_GET
  
  PUBLIC NODES_USER_NUMBER_GET,NODES_USER_NUMBER_SET
  
  PUBLIC NodesUserNumbersAllSet

  !PUBLIC NODES_NUMBER_OF_VERSIONS_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks that a user node number is defined on the specified region.
  SUBROUTINE NODE_CHECK_EXISTS(NODES,USER_NUMBER,NODE_EXISTS,GLOBAL_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to check
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_NUMBER !<On exit, if the node exists the global number corresponding to the user node number. If the node does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
   
    ENTERS("NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    GLOBAL_NUMBER=0
    IF(ASSOCIATED(NODES)) THEN
      NULLIFY(TREE_NODE)
      CALL TREE_SEARCH(NODES%NODES_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
      IF(ASSOCIATED(TREE_NODE)) THEN
        CALL TREE_NODE_VALUE_GET(NODES%NODES_TREE,TREE_NODE,GLOBAL_NUMBER,ERR,ERROR,*999)
        NODE_EXISTS=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NODE_CHECK_EXISTS")
    RETURN
999 ERRORSEXITS("NODE_CHECK_EXISTS",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises a node and deallocates all memory
  SUBROUTINE NODE_FINALISE(NODE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(NODE_TYPE),INTENT(OUT) :: NODE !<The node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("NODE_FINALISE",ERR,ERROR,*999)

    NODE%GLOBAL_NUMBER=0
    NODE%USER_NUMBER=0
    NODE%LABEL=""
    
    EXITS("NODE_FINALISE")
    RETURN
999 ERRORSEXITS("NODE_FINALISE",ERR,ERROR)
    RETURN 1   
  END SUBROUTINE NODE_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nodes in the region. \see OPENCMISS::CMISSNodesCreateFinish
  SUBROUTINE NODES_CREATE_FINISH(NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to be finished
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np
    
    ENTERS("NODES_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        CALL FlagError("Nodes have already been finished.",ERR,ERROR,*999)
      ELSE
        NODES%NODES_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",NODES%NUMBER_OF_NODES,ERR,ERROR,*999)
      DO np=1,NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number    = ",NODES%NODES(np)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number      = ",NODES%NODES(np)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Label            = ",NODES%NODES(np)%LABEL, &
          & ERR,ERROR,*999)
      ENDDO !np
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"User->Global number tree",ERR,ERROR,*999)
      CALL TREE_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,NODES%NODES_TREE,ERR,ERROR,*999)
    ENDIF

    EXITS("NODES_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("NODES_CREATE_FINISH",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE NODES_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating generic nodes.
  SUBROUTINE NODES_CREATE_START_GENERIC(NODES,NUMBER_OF_NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<The nodes pointer
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES !<The number of nodes to create
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: INSERT_STATUS,np
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NODES_CREATE_START_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NUMBER_OF_NODES>0) THEN
        ALLOCATE(NODES%NODES(NUMBER_OF_NODES),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate nodes nodes.",ERR,ERROR,*999)
        NODES%NUMBER_OF_NODES=NUMBER_OF_NODES
        CALL TREE_CREATE_START(NODES%NODES_TREE,ERR,ERROR,*999)
        CALL TREE_INSERT_TYPE_SET(NODES%NODES_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
        CALL TREE_CREATE_FINISH(NODES%NODES_TREE,ERR,ERROR,*999)
        !Set default node numbers
        DO np=1,NODES%NUMBER_OF_NODES
          NODES%NODES(np)%GLOBAL_NUMBER=np
          NODES%NODES(np)%USER_NUMBER=np
          NODES%NODES(np)%LABEL=""
          CALL TREE_ITEM_INSERT(NODES%NODES_TREE,np,np,INSERT_STATUS,ERR,ERROR,*999)
        ENDDO !np
      ELSE
        LOCAL_ERROR="The specified number of nodes of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_NODES,"*",ERR,ERROR))// &
          & " is invalid. The number of nodes must be > 0."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_CREATE_GENERIC")
    RETURN
999 ERRORSEXITS("NODES_CREATE_START_GENERIC",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NODES_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating nodes in an interface.
  SUBROUTINE NODES_CREATE_START_INTERFACE(INTERFACE,NUMBER_OF_NODES,NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface in which to create the nodes
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES !<The number of nodes to create
    TYPE(NODES_TYPE), POINTER :: NODES !<On exit, a pointer to the created nodes. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("NODES_CREATE_START_INTERFACE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(NODES)) THEN
        CALL FlagError("Nodes is already associated.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE%NODES)) THEN
          CALL FlagError("Interface already has nodes associated.",ERR,ERROR,*998)
        ELSE
          !Initialise the nodes for the interface
          CALL NODES_INITIALISE(INTERFACE,ERR,ERROR,*999)
          !Create the nodes 
          CALL NODES_CREATE_START_GENERIC(INTERFACE%NODES,NUMBER_OF_NODES,ERR,ERROR,*999)
          !Return the pointer        
          NODES=>INTERFACE%NODES
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("NODES_CREATE_START_INTERFACE")
    RETURN
999 CALL NODES_FINALISE(INTERFACE%NODES,DUMMY_ERR,DUMMY_ERROR,*998)    
998 ERRORSEXITS("NODES_CREATE_START_INTERFACE",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE NODES_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating nodes in the region.  \see OPENCMISS::CMISSNodesCreateStart
  SUBROUTINE NODES_CREATE_START_REGION(REGION,NUMBER_OF_NODES,NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region in which to create the nodes
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES !<The number of nodes to create
    TYPE(NODES_TYPE), POINTER :: NODES !<On exit, a pointer to the created nodes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    ENTERS("NODES_CREATE_START_REGION",ERR,ERROR,*998)

    NULLIFY(NODES)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        CALL FlagError("Region already has nodes associated.",ERR,ERROR,*998)
      ELSE
        IF(ASSOCIATED(NODES)) THEN
          CALL FlagError("Nodes is already associated.",ERR,ERROR,*998)
        ELSE
          !Initialise the nodes for the region
          CALL NODES_INITIALISE(REGION,ERR,ERROR,*999)
          !Create the generic nodes
          CALL NODES_CREATE_START_GENERIC(REGION%NODES,NUMBER_OF_NODES,ERR,ERROR,*999)
          !Return the pointer        
          NODES=>REGION%NODES
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*998)
    ENDIF
    
    EXITS("NODES_CREATE_START_REGION")
    RETURN
999 CALL NODES_FINALISE(REGION%NODES,DUMMY_ERR,DUMMY_ERROR,*998)    
998 ERRORSEXITS("NODES_CREATE_START_REGION",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE NODES_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys nodes. \see OPENCMISS::CMISSNodesDestroy
  SUBROUTINE NODES_DESTROY(NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("NODES_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(ASSOCIATED(NODES%REGION)) THEN
        NULLIFY(NODES%REGION%NODES)
      ELSE
        IF(ASSOCIATED(NODES%INTERFACE)) THEN
          NULLIFY(NODES%INTERFACE%NODES)
        ELSE
          CALL FlagError("Nodes region and interface are not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
      CALL NODES_FINALISE(NODES,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_DESTROY")
    RETURN
999 ERRORSEXITS("NODES_DESTROY",ERR,ERROR)
    RETURN 1
   
  END SUBROUTINE NODES_DESTROY
    
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes and deallocates any memory. 
  SUBROUTINE NODES_FINALISE(NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    ENTERS("NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(ALLOCATED(NODES%NODES)) THEN
        DO np=1,SIZE(NODES%NODES,1)
          CALL NODE_FINALISE(NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        DEALLOCATE(NODES%NODES)
      ENDIF
      IF(ASSOCIATED(NODES%NODES_TREE)) CALL TREE_DESTROY(NODES%NODES_TREE,ERR,ERROR,*999)
      DEALLOCATE(NODES)
    ENDIF
    
    EXITS("NODES_FINALISE")
    RETURN
999 ERRORSEXITS("NODES_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes.
  SUBROUTINE NODES_INITIALISE_GENERIC(NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("NODES_INITIALISE_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      NULLIFY(NODES%REGION)
      NULLIFY(NODES%INTERFACE)
      NODES%NODES_FINISHED=.FALSE.
      NODES%NUMBER_OF_NODES=0
      NULLIFY(NODES%NODES_TREE)
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_INITIALISE_GENERIC")
    RETURN
999 ERRORSEXITS("NODES_INITIALISE_GENERIC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NODES_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given interface.
  SUBROUTINE NODES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("NODES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%NODES)) THEN
        CALL FlagError("Interface already has associated nodes.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate interface nodes.",ERR,ERROR,*999)
        CALL NODES_INITIALISE_GENERIC(INTERFACE%NODES,ERR,ERROR,*999)
        INTERFACE%NODES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FlagError("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_INITIALISE_INTERFACE")
    RETURN
999 ERRORSEXITS("NODES_INITIALISE_INTERFACE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE NODES_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given region.
  SUBROUTINE NODES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("NODES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%NODES)) THEN
        CALL FlagError("Region has associated nodes.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(REGION%NODES,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate region nodes.",ERR,ERROR,*999)
        CALL NODES_INITIALISE_GENERIC(REGION%NODES,ERR,ERROR,*999)
        REGION%NODES%REGION=>REGION
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_INITIALISE_REGION")
    RETURN
999 ERRORSEXITS("NODES_INITIALISE_REGION",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NODES_INITIALISE_REGION

  !
  !================================================================================================================================
  !

  !>Gets the character label for a node identified by a given global number. \see OPENCMISS::CMISSNodesLabelGet
  SUBROUTINE NODES_LABEL_GET_C(NODES,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to get the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On exit, the label of the specified global node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER :: C_LENGTH,VS_LENGTH
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("NODES_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          C_LENGTH=LEN(LABEL)
          VS_LENGTH=LEN_TRIM(NODES%NODES(GLOBAL_NUMBER)%LABEL)
          IF(C_LENGTH>VS_LENGTH) THEN
            LABEL=CHAR(LEN_TRIM(NODES%NODES(GLOBAL_NUMBER)%LABEL))
          ELSE
            LABEL=CHAR(NODES%NODES(GLOBAL_NUMBER)%LABEL,C_LENGTH)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global node number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global node number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Nodes have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_LABEL_GET_C")
    RETURN
999 ERRORSEXITS("NODES_LABEL_GET_C",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_LABEL_GET_C
        
  !
  !================================================================================================================================
  !

  !>Gets the varying string label for a node identified by a given global number. \see OPENCMISS::CMISSNodesLabelGet
  SUBROUTINE NODES_LABEL_GET_VS(NODES,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to get the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On exit, the label of the specified global node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("NODES_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          LABEL=NODES%NODES(GLOBAL_NUMBER)%LABEL
        ELSE
          LOCAL_ERROR="The specified global node number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global node number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Nodes have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_LABEL_GET_VS")
    RETURN
999 ERRORSEXITS("NODES_LABEL_GET_VS",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_LABEL_GET_VS
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the character label for a node identified by a given global number. \see OPENCMISS::CMISSNodesLabelSet
  SUBROUTINE NODES_LABEL_SET_C(NODES,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to set the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("NODES_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        CALL FlagError("Nodes have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          NODES%NODES(GLOBAL_NUMBER)%LABEL=LABEL
        ELSE
          LOCAL_ERROR="The specified global node number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global node number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_LABEL_SET_C")
    RETURN
999 ERRORSEXITS("NODES_LABEL_SET_C",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_LABEL_SET_C
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the varying string label for a node identified by a given global number. \see OPENCMISS::CMISSNodesLabelSet
  SUBROUTINE NODES_LABEL_SET_VS(NODES,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to set the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("NODES_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        CALL FlagError("Nodes have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          NODES%NODES(GLOBAL_NUMBER)%LABEL=LABEL
        ELSE
          LOCAL_ERROR="The specified global node number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global node number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_LABEL_SET_VS")
    RETURN
999 ERRORSEXITS("NODES_LABEL_SET_VS",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_LABEL_SET_VS
        
  !
  !================================================================================================================================
  !

  !>Returns the number of nodes. \see OPENCMISS::CMISSNodesNumberOfNodesGet
  SUBROUTINE NODES_NUMBER_OF_NODES_GET(NODES,NUMBER_OF_NODES,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to get the number of nodes for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NODES !<On return, the number of nodes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    ENTERS("NODES_NUMBER_OF_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        NUMBER_OF_NODES=NODES%NUMBER_OF_NODES
      ELSE
        CALL FlagError("Nodes have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_NUMBER_OF_NODES_GET")
    RETURN
999 ERRORSEXITS("NODES_NUMBER_OF_NODES_GET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_NUMBER_OF_NODES_GET
        
  !
  !================================================================================================================================
  !

  !>Gets the user number for a node identified by a given global number. \see OPENCMISS::CMISSNodesUserNumberGet
  SUBROUTINE NODES_USER_NUMBER_GET(NODES,GLOBAL_NUMBER,USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to get the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the user number for
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<On exit, the user number of the specified global node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("NODES_USER_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          USER_NUMBER=NODES%NODES(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="The specified global node number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global node number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Nodes have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODES_USER_NUMBER_GET")
    RETURN
999 ERRORSEXITS("NODES_USER_NUMBER_GET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_USER_NUMBER_GET
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a node identified by a given global number. \see OPENCMISS::CMISSNodesUserNumberSet
  SUBROUTINE NODES_USER_NUMBER_SET(NODES,GLOBAL_NUMBER,USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: NODES !<A pointer to the nodes to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the user number for
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: INSERT_STATUS,OLD_GLOBAL_NUMBER
    LOGICAL :: NODE_EXISTS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("NODES_USER_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(NODES)) THEN
      IF(NODES%NODES_FINISHED) THEN
        CALL FlagError("Nodes have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=NODES%NUMBER_OF_NODES) THEN
          !Check the node user number is not already used
          CALL NODE_CHECK_EXISTS(NODES,USER_NUMBER,NODE_EXISTS,OLD_GLOBAL_NUMBER,ERR,ERROR,*999)
          IF(NODE_EXISTS) THEN
            IF(OLD_GLOBAL_NUMBER/=GLOBAL_NUMBER) THEN
              LOCAL_ERROR="The specified node user number of "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                & " is already used by global node number "//TRIM(NUMBER_TO_VSTRING(OLD_GLOBAL_NUMBER,"*",ERR,ERROR))// &
                & ". User node numbers must be unique."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL TREE_ITEM_DELETE(NODES%NODES_TREE,NODES%NODES(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
            CALL TREE_ITEM_INSERT(NODES%NODES_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
            IF(INSERT_STATUS/=TREE_NODE_INSERT_SUCESSFUL) CALL FlagError("Unsucessful nodes tree insert.",ERR,ERROR,*999)
            NODES%NODES(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global node number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global node number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(NODES%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("NODE_USER_NUMBER_SET")
    RETURN
999 ERRORSEXITS("NODE_USER_NUMBER_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE NODES_USER_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets the user numbers for all nodes. \see OPENCMISS::CMISSNodes_UserNumbersAllSet
  SUBROUTINE NodesUserNumbersAllSet(nodes,userNumbers,err,error,*)

    !Argument variables
    TYPE(NODES_TYPE), POINTER :: nodes !<A pointer to the nodes to set the numbers for
    INTEGER(INTG), INTENT(IN) :: userNumbers(:) !<The user numbers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nodeIdx,insertStatus
    TYPE(TREE_TYPE), POINTER :: newNodesTree
    TYPE(VARYING_STRING) :: localError

    NULLIFY(newNodesTree)
    
    ENTERS("NodesUserNumbersAllSet",err,error,*999)

    IF(ASSOCIATED(nodes)) THEN
      IF(nodes%NODES_FINISHED) THEN
        CALL FlagError("Nodes have been finished.",err,error,*999)
      ELSE
        IF(nodes%NUMBER_OF_NODES==SIZE(userNumbers,1)) THEN
          !Check the users numbers to ensure that there are no duplicates                    
          CALL TREE_CREATE_START(newNodesTree,err,error,*999)
          CALL TREE_INSERT_TYPE_SET(newNodesTree,TREE_NO_DUPLICATES_ALLOWED,err,error,*999)
          CALL TREE_CREATE_FINISH(newNodesTree,err,error,*999)
          DO nodeIdx=1,nodes%NUMBER_OF_NODES
            CALL TREE_ITEM_INSERT(newNodesTree,userNumbers(nodeIdx),nodeIdx,insertStatus,err,error,*999)
            IF(insertStatus/=TREE_NODE_INSERT_SUCESSFUL) THEN
              localError="The specified user number of "//TRIM(NumberToVstring(userNumbers(nodeIdx),"*",err,error))// &
                & " for global node number "//TRIM(NumberToVstring(nodeIdx,"*",err,error))// &
                & " is a duplicate. The user node numbers must be unique."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !nodeIdx
          CALL TREE_DESTROY(nodes%NODES_TREE,err,error,*999)
          nodes%NODES_TREE=>newNodesTree
          NULLIFY(newNodesTree)
          DO nodeIdx=1,nodes%NUMBER_OF_NODES
            nodes%NODES(nodeIdx)%GLOBAL_NUMBER=nodeIdx
            nodes%NODES(nodeIdx)%USER_NUMBER=userNumbers(nodeIdx)
          ENDDO !nodesIdx
        ELSE
          localError="The number of specified node user numbers ("// &
            TRIM(NumberToVstring(SIZE(userNumbers,1),"*",err,error))// &
            ") does not match number of nodes ("// &
            TRIM(NumberToVstring(nodes%NUMBER_OF_NODES,"*",err,error))//")."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Nodes is not associated.",err,error,*999)
    ENDIF
    
    EXITS("NodesUserNumbersAllSet")
    RETURN
999 IF(ASSOCIATED(newNodesTree)) CALL TREE_DESTROY(newNodesTree,err,error,*998)
998 ERRORSEXITS("NodesUserNumbersAllSet",err,error)    
    RETURN 1
   
  END SUBROUTINE NodesUserNumbersAllSet
  
  !
  !================================================================================================================================
  !

END MODULE NODE_ROUTINES
