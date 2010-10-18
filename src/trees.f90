!> \file
!> $Id$
!> \author Chris Bradley
!> \brief Implements trees of base types.
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

!> Implements trees of base types.
MODULE TREES

  USE BASE_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup TREES_TreeNodeColourTypes TREES::TreeNodeColourTypes
  !> \brief The colour of the tree nodes
  !> \see TREES
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_BLACK_NODE=0 !<The black colour type for a tree node \see TREES_TreeNodeColourTypes,TREES
  INTEGER(INTG), PARAMETER :: TREE_RED_NODE=1 !<The red colour type for a tree node \see TREES_TreeNodeColourTypes,TREES
  !>@}

  !> \addtogroup TREES_TreeNodeInsertStatus TREES::TreeNodeInsertStatus
  !> \brief The insert status for tree nodes
  !> \see TREES
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_NODE_INSERT_SUCESSFUL=1 !<Successful insert status \see TREES_TreeNodeInsertStatus,TREES
  INTEGER(INTG), PARAMETER :: TREE_NODE_DUPLICATE_KEY=2 !<Duplicate key found for those trees that do not allow duplicate keys \see TREES_TreeNodeInsertStatus,TREES
  !>@}

  !> \addtogroup TREES_TreeInsertTypes TREES::TreeInsertTypes
  !> \brief The insert type for a tree
  !> \see TREES
  !>@{
  INTEGER(INTG), PARAMETER :: TREE_DUPLICATES_ALLOWED_TYPE=1 !<Duplicate keys allowed tree type \see TREES_TreeInsertTypes,TREES
  INTEGER(INTG), PARAMETER :: TREE_NO_DUPLICATES_ALLOWED=2 !<No duplicate keys allowed tree type \see TREES_TreeInsertTypes,TREES
  !>@}

  !Module types

  TYPE TREE_NODE_TYPE
    PRIVATE
    INTEGER(INTG) :: KEY !<The key for the tree node
    INTEGER(INTG) :: VALUE !<The value stored at the tree node
    INTEGER(INTG) :: COLOUR !<The colour of the node for the red-black tree
    TYPE(TREE_NODE_TYPE), POINTER :: LEFT !<The pointer to the left daughter tree node if any
    TYPE(TREE_NODE_TYPE), POINTER :: RIGHT !<The pointer to the right daughter tree node if any
    TYPE(TREE_NODE_TYPE), POINTER :: PARENT !<The pointer to the parent tree node
  END TYPE TREE_NODE_TYPE

  TYPE TREE_TYPE
    PRIVATE
    LOGICAL :: TREE_FINISHED !<Is .TRUE. if the tree has finished being created, .FALSE. if not.
    INTEGER(INTG) :: INSERT_TYPE !<The insert type for duplicate keys for the tree
    INTEGER(INTG) :: NUMBER_IN_TREE !<The number of items currently in the tree
    TYPE(TREE_NODE_TYPE), POINTER :: ROOT !<The pointer to the root of the tree
    TYPE(TREE_NODE_TYPE), POINTER :: NIL !<The pointer to the nil of the tree
  END TYPE TREE_TYPE

  !Module variables
  
  !Interfaces

  PUBLIC TREE_TYPE,TREE_NODE_TYPE

  PUBLIC TREE_NODE_INSERT_SUCESSFUL,TREE_NODE_DUPLICATE_KEY

  PUBLIC TREE_DUPLICATES_ALLOWED_TYPE,TREE_NO_DUPLICATES_ALLOWED
  
  PUBLIC TREE_CREATE_FINISH,TREE_CREATE_START,TREE_DESTROY,TREE_DETACH,TREE_DETACH_AND_DESTROY,TREE_INSERT_TYPE_SET, &
    & TREE_ITEM_DELETE,TREE_ITEM_INSERT,TREE_NODE_KEY_GET,TREE_NODE_VALUE_GET,TREE_NODE_VALUE_SET,TREE_OUTPUT,TREE_SEARCH 

CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Finishes the creation of a tree created with TREE_CREATE_START \see{TREES::TREE_CREATE_START}.
  SUBROUTINE TREE_CREATE_FINISH(TREE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to finish
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        CALL FLAG_ERROR("Tree is already finished",ERR,ERROR,*998)
      ELSE
        !Allocate the nil tree node
        ALLOCATE(TREE%NIL,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NIL tree node",ERR,ERROR,*999)
        CALL TREE_NODE_INITIALISE(TREE,TREE%NIL,ERR,ERROR,*999)
        TREE%NIL%KEY=-99999999 !Set it to something identifiable for debugging
        TREE%NIL%LEFT=>TREE%NIL
        TREE%NIL%RIGHT=>TREE%NIL
        TREE%NIL%PARENT=>TREE%NIL
        !Set the root tree node to NIL        
        TREE%ROOT=>TREE%NIL
        !Finish the tree creation
        TREE%TREE_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("TREE_CREATE_FINISH")
    RETURN
999 CALL TREE_FINALISE(TREE,ERR,ERROR,*998)
998 CALL ERRORS("TREE_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("TREE_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE TREE_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a tree and returns a pointer to the created tree \see{TREES::TREE_CREATE_FINISH}.
  SUBROUTINE TREE_CREATE_START(TREE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variable

    CALL ENTERS("TREE_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(TREE)) THEN
      CALL FLAG_ERROR("Tree is already associated",ERR,ERROR,*998)
    ELSE
      ALLOCATE(TREE,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate tree",ERR,ERROR,*999)
      CALL TREE_INITIALISE(TREE,ERR,ERROR,*999)
      !Set Defaults
      TREE%INSERT_TYPE=TREE_DUPLICATES_ALLOWED_TYPE
    ENDIF

    CALL EXITS("TREE_CREATE_START")
    RETURN
999 CALL TREE_FINALISE(TREE,ERR,ERROR,*998)
998 CALL ERRORS("TREE_CREATE_START",ERR,ERROR)
    CALL EXITS("TREE_CREATE_START")
    RETURN 1
  END SUBROUTINE TREE_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a tree
  SUBROUTINE TREE_DESTROY(TREE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      CALL TREE_FINALISE(TREE,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_DESTROY")
    RETURN
999 CALL ERRORS("TREE_DESTROY",ERR,ERROR)
    CALL EXITS("TREE_DESTROY")
    RETURN 1
  END SUBROUTINE TREE_DESTROY

  !
  !================================================================================================================================
  !

  !>Detaches the tree values and returns them as a pointer to the an array
  SUBROUTINE TREE_DETACH(TREE,NUMBER_IN_TREE,TREE_VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to detach
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_TREE !<On exit, the number in the array that has been detached
    INTEGER(INTG), POINTER :: TREE_VALUES(:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_DETACH",ERR,ERROR,*998)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        IF(ASSOCIATED(TREE_VALUES)) THEN
          CALL FLAG_ERROR("Tree values is already associated.",ERR,ERROR,*998)
        ELSE
          NULLIFY(TREE_VALUES)
          ALLOCATE(TREE_VALUES(TREE%NUMBER_IN_TREE),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate tree values.",ERR,ERROR,*999)
          NUMBER_IN_TREE=0
          CALL TREE_DETACH_IN_ORDER(TREE,TREE%ROOT,NUMBER_IN_TREE,TREE_VALUES,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Tree has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("TREE_DETACH")
    RETURN
999 IF(ASSOCIATED(TREE_VALUES)) DEALLOCATE(TREE_VALUES)
    NUMBER_IN_TREE=0
998 CALL ERRORS("TREE_DETACH",ERR,ERROR)
    CALL EXITS("TREE_DETACH")
    RETURN 1
  END SUBROUTINE TREE_DETACH

  !
  !================================================================================================================================
  !

  !>Detaches the tree values and returns them as a pointer to the an array and then destroys the tree
  SUBROUTINE TREE_DETACH_AND_DESTROY(TREE,NUMBER_IN_TREE,TREE_VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to detach and destroy
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_TREE !<On exit, the number in the array that has been detached
    INTEGER(INTG), POINTER :: TREE_VALUES(:) !<On exit, a pointer to the dettached tree values. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_DETACH_AND_DESTROY",ERR,ERROR,*998)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        IF(ASSOCIATED(TREE_VALUES)) THEN
          CALL FLAG_ERROR("Tree values is associated",ERR,ERROR,*998)
        ELSE
          NULLIFY(TREE_VALUES)
          ALLOCATE(TREE_VALUES(TREE%NUMBER_IN_TREE),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate tree values",ERR,ERROR,*999)
          NUMBER_IN_TREE=0
          CALL TREE_DETACH_IN_ORDER(TREE,TREE%ROOT,NUMBER_IN_TREE,TREE_VALUES,ERR,ERROR,*999)
          CALL TREE_FINALISE(TREE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("TREE_DETACH_AND_DESTROY")
    RETURN
999 IF(ASSOCIATED(TREE_VALUES)) DEALLOCATE(TREE_VALUES)
    NUMBER_IN_TREE=0
998 CALL ERRORS("TREE_DETACH_AND_DESTROY",ERR,ERROR)
    CALL EXITS("TREE_DETACH_AND_DESTROY")
    RETURN 1
  END SUBROUTINE TREE_DETACH_AND_DESTROY

  !
  !================================================================================================================================
  !

  !>Detaches the tree values in order from the specified tree node and adds them to the tree values array
  RECURSIVE SUBROUTINE TREE_DETACH_IN_ORDER(TREE,X,COUNT,TREE_VALUES,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to detach
    TYPE(TREE_NODE_TYPE), POINTER :: X !<A pointer to the specified tree node to detach from
    INTEGER(INTG), INTENT(INOUT) :: COUNT !<The current number in the detached tree values array
    INTEGER(INTG), INTENT(INOUT) :: TREE_VALUES(:) !<The current detached tree values array
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("TREE_DETACH_IN_ORDER",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(.NOT.ASSOCIATED(X,TREE%NIL)) THEN
        CALL TREE_DETACH_IN_ORDER(TREE,X%LEFT,COUNT,TREE_VALUES,ERR,ERROR,*999)
        COUNT=COUNT+1
        IF(COUNT<=SIZE(TREE_VALUES,1)) THEN
          TREE_VALUES(COUNT)=X%VALUE
        ELSE
          LOCAL_ERROR="The current count of the tree values ("//TRIM(NUMBER_TO_VSTRING(COUNT,"*",ERR,ERROR))// &
            & ") is greater than the size of the tree values array ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(TREE_VALUES,1),"*",ERR,ERROR))//")"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
        CALL TREE_DETACH_IN_ORDER(TREE,X%RIGHT,COUNT,TREE_VALUES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_DETACH_IN_ORDER")
    RETURN
999 CALL ERRORS("TREE_DETACH_IN_ORDER",ERR,ERROR)
    CALL EXITS("TREE_DETACH_IN_ORDER")
    RETURN 1
  END SUBROUTINE TREE_DETACH_IN_ORDER

  !
  !================================================================================================================================
  !

  !>Finalises a tree and deallocates all memory
  SUBROUTINE TREE_FINALISE(TREE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      CALL TREE_NODE_FINALISE(TREE,TREE%ROOT,ERR,ERROR,*999)
      IF(ASSOCIATED(TREE%NIL)) DEALLOCATE(TREE%NIL)
      DEALLOCATE(TREE)
    ENDIF

    CALL EXITS("TREE_FINALISE")
    RETURN
999 CALL ERRORS("TREE_FINALISE",ERR,ERROR)
    CALL EXITS("TREE_FINALISE")
    RETURN 1
  END SUBROUTINE TREE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a tree
  SUBROUTINE TREE_INITIALISE(TREE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      TREE%TREE_FINISHED=.FALSE.
      TREE%INSERT_TYPE=0
      TREE%NUMBER_IN_TREE=0
      NULLIFY(TREE%ROOT)
      NULLIFY(TREE%NIL)
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_INITIALISE")
    RETURN
999 CALL ERRORS("TREE_INITIALISE",ERR,ERROR)
    CALL EXITS("TREE_INITIALISE")
    RETURN 1
  END SUBROUTINE TREE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the insert type for a tree
  SUBROUTINE TREE_INSERT_TYPE_SET(TREE,INSERT_TYPE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree
    INTEGER(INTG), INTENT(IN) :: INSERT_TYPE !<The insert type to set \see TREES_TreeInsertTypes,TREES::TreeInsertTypes
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("TREE_INSERT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        CALL FLAG_ERROR("Tree has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(INSERT_TYPE)
        CASE(TREE_DUPLICATES_ALLOWED_TYPE)
          TREE%INSERT_TYPE=TREE_DUPLICATES_ALLOWED_TYPE
        CASE(TREE_NO_DUPLICATES_ALLOWED)
          TREE%INSERT_TYPE=TREE_NO_DUPLICATES_ALLOWED
        CASE DEFAULT
          LOCAL_ERROR="The insert type of "//TRIM(NUMBER_TO_VSTRING(INSERT_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_INSERT_TYPE_SET")
    RETURN
999 CALL ERRORS("TREE_INSERT_TYPE_SET",ERR,ERROR)
    CALL EXITS("TREE_INSERT_TYPE_SET")
    RETURN 1
  END SUBROUTINE TREE_INSERT_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Deletes a tree node specified by a key from a tree 
  SUBROUTINE TREE_ITEM_DELETE(TREE,KEY,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the Red-Black tree to delete from
    INTEGER(INTG), INTENT(IN) :: KEY !<A pointer to the tree node to delete
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    INTEGER(INTG) :: COMPARE_VALUE
    TYPE(TREE_NODE_TYPE), POINTER :: U,V,W,X,Y,Z
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("TREE_ITEM_DELETE",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        !Try and find the key to delete
        Z=>TREE%ROOT
        IF(.NOT.ASSOCIATED(Z,TREE%NIL)) THEN
          COMPARE_VALUE=Z%KEY-KEY
          DO WHILE(COMPARE_VALUE/=0)
            IF(COMPARE_VALUE>0) THEN !Z%KEY > KEY
              Z=>Z%LEFT
            ELSE !Z%KEY < KEY
              Z=>Z%RIGHT
            ENDIF
            IF(ASSOCIATED(Z,TREE%NIL)) THEN
              EXIT
            ELSE
              COMPARE_VALUE=Z%KEY-KEY
            ENDIF
          ENDDO
          IF(COMPARE_VALUE==0) THEN
            !Found the key so delete it
            IF(ASSOCIATED(Z%LEFT,TREE%NIL).OR.ASSOCIATED(Z%RIGHT,TREE%NIL)) THEN
              Y=>Z
            ELSE
              Y=>TREE_SUCCESSOR(TREE,Z,ERR,ERROR)
              IF(ERR/=0) GOTO 999
            ENDIF
            IF(.NOT.ASSOCIATED(Y%LEFT,TREE%NIL)) THEN
              X=>Y%LEFT
            ELSE
              X=>Y%RIGHT
            ENDIF
            X%PARENT=>Y%PARENT
            IF(ASSOCIATED(Y%PARENT,TREE%NIL)) THEN
              TREE%ROOT=>X
            ELSE
              IF(ASSOCIATED(Y,Y%PARENT%LEFT)) THEN
                Y%PARENT%LEFT=>X
              ELSE
                Y%PARENT%RIGHT=>X
              ENDIF
            ENDIF
            IF(Y%COLOUR==TREE_BLACK_NODE) THEN
              !Fixup the delete to ensure the tree has red black properties
              !Note: Due to Fortran restrictions on aliasing pointers in dummy arguments we need to do the fixup and rotations
              !inside this routine rather than call fixup and rotate left and rotate right subroutines.
              DO WHILE(.NOT.ASSOCIATED(X,TREE%ROOT).AND.X%COLOUR==TREE_BLACK_NODE)
                IF(ASSOCIATED(X,X%PARENT%LEFT)) THEN
                  W=>X%PARENT%RIGHT
                  IF(W%COLOUR==TREE_RED_NODE) THEN
                    W%COLOUR=TREE_BLACK_NODE
                    X%PARENT%COLOUR=TREE_RED_NODE
                    !Rotate left on X->Parent
                    U=>X%PARENT
                    V=>U%RIGHT
                    U%RIGHT=>V%LEFT
                    IF(.NOT.ASSOCIATED(V%LEFT,TREE%NIL)) V%LEFT%PARENT=>U
                    V%PARENT=>U%PARENT
                    IF(ASSOCIATED(U%PARENT,TREE%NIL)) THEN
                      TREE%ROOT=>V
                    ELSE
                      IF(ASSOCIATED(U,U%PARENT%LEFT)) THEN
                        U%PARENT%LEFT=>V
                      ELSE
                        U%PARENT%RIGHT=>V
                      ENDIF
                    ENDIF
                    V%LEFT=>U
                    U%PARENT=>V
                    W=>X%PARENT%RIGHT
                  ENDIF
                  IF(W%LEFT%COLOUR==TREE_BLACK_NODE.AND.W%RIGHT%COLOUR==TREE_BLACK_NODE) THEN
                    W%COLOUR=TREE_RED_NODE
                    X=>X%PARENT
                  ELSE
                    IF(W%RIGHT%COLOUR==TREE_BLACK_NODE) THEN
                      W%LEFT%COLOUR=TREE_BLACK_NODE
                      W%COLOUR=TREE_RED_NODE
                      !Rotate right on W
                      U=>W
                      V=>U%LEFT
                      U%LEFT=>V%RIGHT
                      IF(.NOT.ASSOCIATED(V%RIGHT,TREE%NIL)) V%RIGHT%PARENT=>U
                      V%PARENT=>U%PARENT
                      IF(ASSOCIATED(V%PARENT,TREE%NIL)) THEN
                        TREE%ROOT=>V
                      ELSE
                        IF(ASSOCIATED(U,U%PARENT%RIGHT)) THEN
                          U%PARENT%RIGHT=>V
                        ELSE
                          U%PARENT%LEFT=>V
                        ENDIF
                      ENDIF
                      V%RIGHT=>U
                      U%PARENT=>V
                      W=>X%PARENT%RIGHT
                    ENDIF
                    W%COLOUR=X%PARENT%COLOUR
                    X%PARENT%COLOUR=TREE_BLACK_NODE
                    W%RIGHT%COLOUR=TREE_BLACK_NODE
                    !Rotate left on X->Parent
                    U=>X%PARENT
                    V=>U%RIGHT
                    U%RIGHT=>V%LEFT
                    IF(.NOT.ASSOCIATED(V%LEFT,TREE%NIL)) V%LEFT%PARENT=>U
                    V%PARENT=>U%PARENT
                    IF(ASSOCIATED(U%PARENT,TREE%NIL)) THEN
                      TREE%ROOT=>V
                    ELSE
                      IF(ASSOCIATED(U,U%PARENT%LEFT)) THEN
                        U%PARENT%LEFT=>V
                      ELSE
                        U%PARENT%RIGHT=>V
                      ENDIF
                    ENDIF
                    V%LEFT=>U
                    U%PARENT=>V
                    X=>TREE%ROOT
                  ENDIF
                ELSE
                  W=>X%PARENT%LEFT
                  IF(W%COLOUR==TREE_RED_NODE) THEN
                    W%COLOUR=TREE_BLACK_NODE
                    X%PARENT%COLOUR=TREE_RED_NODE
                    !Rotate right on X->Parent
                    U=>X%PARENT
                    V=>U%LEFT
                    U%LEFT=>V%RIGHT
                    IF(.NOT.ASSOCIATED(V%RIGHT,TREE%NIL)) V%RIGHT%PARENT=>U
                    V%PARENT=>U%PARENT
                    IF(ASSOCIATED(V%PARENT,TREE%NIL)) THEN
                      TREE%ROOT=>V
                    ELSE
                      IF(ASSOCIATED(U,U%PARENT%RIGHT)) THEN
                        U%PARENT%RIGHT=>V
                      ELSE
                        U%PARENT%LEFT=>V
                      ENDIF
                    ENDIF
                    V%RIGHT=>U
                    U%PARENT=>V
                    W=>X%PARENT%LEFT
                  ENDIF
                  IF(W%RIGHT%COLOUR==TREE_BLACK_NODE.AND.W%LEFT%COLOUR==TREE_BLACK_NODE) THEN
                    W%COLOUR=TREE_RED_NODE
                    X=>X%PARENT
                  ELSE
                    IF(W%LEFT%COLOUR==TREE_BLACK_NODE) THEN
                      W%RIGHT%COLOUR=TREE_BLACK_NODE
                      W%COLOUR=TREE_RED_NODE
                      !Rotate left on W
                      U=>W
                      V=>U%RIGHT
                      U%RIGHT=>V%LEFT
                      IF(.NOT.ASSOCIATED(V%LEFT,TREE%NIL)) V%LEFT%PARENT=>U
                      V%PARENT=>U%PARENT
                      IF(ASSOCIATED(U%PARENT,TREE%NIL)) THEN
                        TREE%ROOT=>V
                      ELSE
                        IF(ASSOCIATED(U,U%PARENT%LEFT)) THEN
                          U%PARENT%LEFT=>V
                        ELSE
                          U%PARENT%RIGHT=>V
                        ENDIF
                      ENDIF
                      V%LEFT=>U
                      U%PARENT=>V
                      W=>X%PARENT%LEFT
                    ENDIF
                    W%COLOUR=X%PARENT%COLOUR
                    X%PARENT%COLOUR=TREE_BLACK_NODE
                    W%LEFT%COLOUR=TREE_BLACK_NODE
                    !Rotate right on X->Parent
                    U=>X%PARENT
                    V=>U%LEFT
                    U%LEFT=>V%RIGHT
                    IF(.NOT.ASSOCIATED(V%RIGHT,TREE%NIL)) V%RIGHT%PARENT=>U
                    V%PARENT=>U%PARENT
                    IF(ASSOCIATED(V%PARENT,TREE%NIL)) THEN
                      TREE%ROOT=>V
                    ELSE
                      IF(ASSOCIATED(U,U%PARENT%RIGHT)) THEN
                        U%PARENT%RIGHT=>V
                      ELSE
                        U%PARENT%LEFT=>V
                      ENDIF
                    ENDIF
                    V%RIGHT=>U
                    U%PARENT=>V
                    X=>TREE%ROOT
                  ENDIF
                ENDIF
              ENDDO
              X%COLOUR=TREE_BLACK_NODE
            ENDIF
            IF(.NOT.ASSOCIATED(Y,Z)) THEN
              Y%LEFT=>Z%LEFT
              Y%RIGHT=>Z%RIGHT
              Y%PARENT=>Z%PARENT
              Y%COLOUR=Z%COLOUR
              Z%LEFT%PARENT=>Y
              Z%RIGHT%PARENT=>Y
              IF(ASSOCIATED(Z,Z%PARENT%LEFT)) THEN
                Z%PARENT%LEFT=>Y              
              ELSE
                Z%PARENT%RIGHT=>Y
              ENDIF
            ENDIF
            DEALLOCATE(Z)
            TREE%NUMBER_IN_TREE=TREE%NUMBER_IN_TREE-1
          ELSE
            LOCAL_ERROR="Could not find the key "//TRIM(NUMBER_TO_VSTRING(KEY,"*",ERR,ERROR))//" in the tree"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The tree root is NIL. Can not delete the key",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("TREE_ITEM_DELETE")
    RETURN
999 CALL ERRORS("TREE_ITEM_DELETE",ERR,ERROR)
    CALL EXITS("TREE_ITEM_DELETE")
    RETURN 1
  END SUBROUTINE TREE_ITEM_DELETE

  !
  !================================================================================================================================
  !

  !>Inserts a tree node into a red-black tree 
  SUBROUTINE TREE_ITEM_INSERT(TREE,KEY,VALUE,INSERT_STATUS,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the Red-Black tree to insert into
    INTEGER(INTG), INTENT(IN) :: KEY !<The key to insert
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to insert
    INTEGER(INTG), INTENT(OUT) :: INSERT_STATUS !<On exit, the status of the insert \see TREES_TreeNodeInsertStatus,TREES::TreeNodeInsertStatus
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    LOGICAL :: DUPLICATE_KEY
    TYPE(TREE_NODE_TYPE), POINTER :: NEW_TREE_NODE,X,Y,Z

    NULLIFY(NEW_TREE_NODE)
    
    CALL ENTERS("TREE_ITEM_INSERT",ERR,ERROR,*998)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        !Find the position to insert
        Y=>TREE%NIL
        X=>TREE%ROOT
        DUPLICATE_KEY=.FALSE.
        DO WHILE(.NOT.ASSOCIATED(X,TREE%NIL))
          Y=>X
          DUPLICATE_KEY=TREE%INSERT_TYPE==TREE_NO_DUPLICATES_ALLOWED.AND.KEY==X%KEY
          IF(DUPLICATE_KEY) THEN
            EXIT
          ELSE IF(KEY<X%KEY) THEN
            X=>X%LEFT
          ELSE
            X=>X%RIGHT
          ENDIF
        ENDDO
        IF(DUPLICATE_KEY) THEN
          INSERT_STATUS=TREE_NODE_DUPLICATE_KEY
        ELSE
          !Allocate the new tree node and set its key and value
          ALLOCATE(NEW_TREE_NODE,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new tree node",ERR,ERROR,*999)
          CALL TREE_NODE_INITIALISE(TREE,NEW_TREE_NODE,ERR,ERROR,*999)
          NEW_TREE_NODE%KEY=KEY
          NEW_TREE_NODE%VALUE=VALUE
          !Insert the new tree node into the tree
          NEW_TREE_NODE%COLOUR=TREE_RED_NODE
          NEW_TREE_NODE%LEFT=>TREE%NIL
          NEW_TREE_NODE%RIGHT=>TREE%NIL
          NEW_TREE_NODE%PARENT=>Y
          IF(ASSOCIATED(Y,TREE%NIL)) THEN
            TREE%ROOT=>NEW_TREE_NODE
          ELSE
            IF(NEW_TREE_NODE%KEY<Y%KEY) THEN
              Y%LEFT=>NEW_TREE_NODE
            ELSE
              Y%RIGHT=>NEW_TREE_NODE
            ENDIF
          ENDIF
          !Fix up the tree to keep red-black properties
          !Note: Due to Fortran restrictions on aliasing pointers in dummy arguments we need to do the fixup and rotations
          !inside this routine rather than call fixup and rotate left and rotate right subroutines.
          Z=>NEW_TREE_NODE
          DO WHILE(Z%PARENT%COLOUR==TREE_RED_NODE)
            IF(ASSOCIATED(Z%PARENT,Z%PARENT%PARENT%LEFT)) THEN
              Y=>Z%PARENT%PARENT%RIGHT
              IF(Y%COLOUR==TREE_RED_NODE) THEN
                Z%PARENT%COLOUR=TREE_BLACK_NODE
                Y%COLOUR=TREE_BLACK_NODE
                Z%PARENT%PARENT%COLOUR=TREE_RED_NODE
                Z=>Z%PARENT%PARENT
              ELSE
                IF(ASSOCIATED(Z,Z%PARENT%RIGHT)) THEN
                  Z=>Z%PARENT
                  !Rotate the tree left at Z
                  X=>Z                  
                  Y=>X%RIGHT
                  X%RIGHT=>Y%LEFT
                  IF(.NOT.ASSOCIATED(Y%LEFT,TREE%NIL)) Y%LEFT%PARENT=>X
                  Y%PARENT=>X%PARENT
                  IF(ASSOCIATED(X%PARENT,TREE%NIL)) THEN
                    TREE%ROOT=>Y
                  ELSE
                    IF(ASSOCIATED(X,X%PARENT%LEFT)) THEN
                      X%PARENT%LEFT=>Y
                    ELSE
                      X%PARENT%RIGHT=>Y
                    ENDIF
                  ENDIF
                  Y%LEFT=>X
                  X%PARENT=>Y
                ENDIF
                Z%PARENT%COLOUR=TREE_BLACK_NODE
                Z%PARENT%PARENT%COLOUR=TREE_RED_NODE
                !Rotate the tree right at Z->Parent->Parent
                X=>Z%PARENT%PARENT
                Y=>X%LEFT
                X%LEFT=>Y%RIGHT
                IF(.NOT.ASSOCIATED(Y%RIGHT,TREE%NIL)) Y%RIGHT%PARENT=>X
                Y%PARENT=>X%PARENT
                IF(ASSOCIATED(X%PARENT,TREE%NIL)) THEN
                  TREE%ROOT=>Y
                ELSE
                  IF(ASSOCIATED(X,X%PARENT%RIGHT)) THEN
                    X%PARENT%RIGHT=>Y
                  ELSE
                    X%PARENT%LEFT=>Y
                  ENDIF
                ENDIF
                Y%RIGHT=>X
                X%PARENT=>Y
             ENDIF
            ELSE
              Y=>Z%PARENT%PARENT%LEFT
              IF(Y%COLOUR==TREE_RED_NODE) THEN
                Z%PARENT%COLOUR=TREE_BLACK_NODE
                Y%COLOUR=TREE_BLACK_NODE
                Z%PARENT%PARENT%COLOUR=TREE_RED_NODE
                Z=>Z%PARENT%PARENT
              ELSE
                IF(ASSOCIATED(Z,Z%PARENT%LEFT)) THEN
                  Z=>Z%PARENT
                  X=>Z
                  !Rotate the tree right at Z
                  Y=>X%LEFT
                  X%LEFT=>Y%RIGHT
                  IF(.NOT.ASSOCIATED(Y%RIGHT,TREE%NIL)) Y%RIGHT%PARENT=>X
                  Y%PARENT=>X%PARENT
                  IF(ASSOCIATED(X%PARENT,TREE%NIL)) THEN
                    TREE%ROOT=>Y
                  ELSE
                    IF(ASSOCIATED(X,X%PARENT%RIGHT)) THEN
                      X%PARENT%RIGHT=>Y
                    ELSE
                      X%PARENT%LEFT=>Y
                    ENDIF
                  ENDIF
                  Y%RIGHT=>X
                  X%PARENT=>Y
                ENDIF
                Z%PARENT%COLOUR=TREE_BLACK_NODE
                Z%PARENT%PARENT%COLOUR=TREE_RED_NODE
                !Rotate the tree left at Z->Parent->Parent
                X=>Z%PARENT%PARENT
                Y=>X%RIGHT
                X%RIGHT=>Y%LEFT
                IF(.NOT.ASSOCIATED(Y%LEFT,TREE%NIL)) Y%LEFT%PARENT=>X
                Y%PARENT=>X%PARENT
                IF(ASSOCIATED(X%PARENT,TREE%NIL)) THEN
                  TREE%ROOT=>Y
                ELSE
                  IF(ASSOCIATED(X,X%PARENT%LEFT)) THEN
                    X%PARENT%LEFT=>Y
                  ELSE
                    X%PARENT%RIGHT=>Y
                  ENDIF
                ENDIF
                Y%LEFT=>X
                X%PARENT=>Y
              ENDIF
            ENDIF
          ENDDO
          TREE%ROOT%COLOUR=TREE_BLACK_NODE
          !Increment the number in the tree and indicate a successful insertion
          TREE%NUMBER_IN_TREE=TREE%NUMBER_IN_TREE+1
          INSERT_STATUS=TREE_NODE_INSERT_SUCESSFUL
        ENDIF
      ELSE
        CALL FLAG_ERROR("The tree has not been finished",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("TREE_ITEM_INSERT")
    RETURN
999 IF(ASSOCIATED(NEW_TREE_NODE)) DEALLOCATE(NEW_TREE_NODE)
998 CALL ERRORS("TREE_ITEM_INSERT",ERR,ERROR)
    CALL EXITS("TREE_ITEM_INSERT")
    RETURN 1
  END SUBROUTINE TREE_ITEM_INSERT

  !
  !================================================================================================================================
  !

  !>Finalises a tree node and deallocates all memory
  RECURSIVE SUBROUTINE TREE_NODE_FINALISE(TREE,TREE_NODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree containing the tree node to finalise
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE !<A pointer to the tree node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_NODE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(.NOT.ASSOCIATED(TREE_NODE,TREE%NIL)) THEN
        CALL TREE_NODE_FINALISE(TREE,TREE_NODE%LEFT,ERR,ERROR,*999)
        CALL TREE_NODE_FINALISE(TREE,TREE_NODE%RIGHT,ERR,ERROR,*999)
        DEALLOCATE(TREE_NODE)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_NODE_FINALISE")
    RETURN
999 CALL ERRORS("TREE_NODE_FINALISE",ERR,ERROR)
    CALL EXITS("TREE_NODE_FINALISE")
    RETURN 1
  END SUBROUTINE TREE_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a tree node
  SUBROUTINE TREE_NODE_INITIALISE(TREE,TREE_NODE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree containing the tree node to initialise
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE !<A pointer to the tree node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_NODE_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(ASSOCIATED(TREE_NODE)) THEN
        TREE_NODE%KEY=0
        TREE_NODE%VALUE=0
        TREE_NODE%COLOUR=TREE_BLACK_NODE
        NULLIFY(TREE_NODE%LEFT)
        NULLIFY(TREE_NODE%RIGHT)
        NULLIFY(TREE_NODE%PARENT)
      ELSE
        CALL FLAG_ERROR("Tree node is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_NODE_INITIALISE")
    RETURN
999 CALL ERRORS("TREE_NODE_INITIALISE",ERR,ERROR)
    CALL EXITS("TREE_NODE_INITIALISE")
    RETURN 1
  END SUBROUTINE TREE_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the key at a specified tree node
  SUBROUTINE TREE_NODE_KEY_GET(TREE,TREE_NODE,KEY,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree containing the tree node 
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE !<A pointer to the tree node to get the key of
    INTEGER(INTG), INTENT(OUT) :: KEY !<On exit, the key of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_NODE_KEY_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        IF(ASSOCIATED(TREE_NODE)) THEN
          KEY=TREE_NODE%KEY
        ELSE
          CALL FLAG_ERROR("Tree node is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_NODE_KEY_GET")
    RETURN
999 CALL ERRORS("TREE_NODE_KEY_GET",ERR,ERROR)
    CALL EXITS("TREE_NODE_KEY_GET")
    RETURN 1
  END SUBROUTINE TREE_NODE_KEY_GET

  !
  !================================================================================================================================
  !

  !>Gets the value at a specified tree node
  SUBROUTINE TREE_NODE_VALUE_GET(TREE,TREE_NODE,VALUE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree containing the tree node 
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE !<A pointer to the tree node to get the value of
    INTEGER(INTG), INTENT(OUT) :: VALUE !<On exit, the value of the specified tree node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_NODE_VALUE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        IF(ASSOCIATED(TREE_NODE)) THEN
          VALUE=TREE_NODE%VALUE
        ELSE
          CALL FLAG_ERROR("Tree node is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_NODE_VALUE_GET")
    RETURN
999 CALL ERRORS("TREE_NODE_VALUE_GET",ERR,ERROR)
    CALL EXITS("TREE_NODE_VALUE_GET")
    RETURN 1
  END SUBROUTINE TREE_NODE_VALUE_GET

  !
  !================================================================================================================================
  !

  !>Sets the value at a specified tree node
  SUBROUTINE TREE_NODE_VALUE_SET(TREE,TREE_NODE,VALUE,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree containing the tree node 
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE !<A pointer to the tree node to set the value of
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value of the specified tree node to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_NODE_VALUE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        IF(ASSOCIATED(TREE_NODE)) THEN
          TREE_NODE%VALUE=VALUE
        ELSE
          CALL FLAG_ERROR("Tree node is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_NODE_VALUE_SET")
    RETURN
999 CALL ERRORS("TREE_NODE_VALUE_SET",ERR,ERROR)
    CALL EXITS("TREE_NODE_VALUE_SET")
    RETURN 1
  END SUBROUTINE TREE_NODE_VALUE_SET

  !
  !================================================================================================================================
  !

  !>Outputs a tree to the specified output stream ID
  SUBROUTINE TREE_OUTPUT(ID,TREE,ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to search
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_OUTPUT",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        CALL WRITE_STRING(ID,"Tree:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Number of tree nodes = ",TREE%NUMBER_IN_TREE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"Tree insert type = ",TREE%INSERT_TYPE,ERR,ERROR,*999)
        CALL TREE_OUTPUT_IN_ORDER(ID,TREE,TREE%ROOT,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("The tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_OUTPUT")
    RETURN
999 CALL ERRORS("TREE_OUTPUT",ERR,ERROR)
    CALL EXITS("TREE_OUTPUT")
    RETURN 1
  END SUBROUTINE TREE_OUTPUT

  !
  !================================================================================================================================
  !

  !>Outputs a tree in order to the specified output stream ID from the specified tree node
  RECURSIVE SUBROUTINE TREE_OUTPUT_IN_ORDER(ID,TREE,X,ERR,ERROR,*)

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to search
    TYPE(TREE_NODE_TYPE), POINTER :: X !<A pointer to the tree node to output from
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    
    CALL ENTERS("TREE_OUTPUT_IN_ORDER",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(.NOT.ASSOCIATED(X,TREE%NIL)) THEN
        !Output the left subtree first
        CALL TREE_OUTPUT_IN_ORDER(ID,TREE,X%LEFT,ERR,ERROR,*999)
        !Now output the information for this node
        CALL WRITE_STRING(ID,"  Tree Node:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"    Key = ",X%KEY,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"    Value = ",X%VALUE,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(ID,"    Colour = ",X%COLOUR,ERR,ERROR,*999)
        IF(ASSOCIATED(X%LEFT,TREE%NIL)) THEN
          CALL WRITE_STRING(ID,"    Left Key = NIL",ERR,ERROR,*999)
        ELSE
          CALL WRITE_STRING_VALUE(ID,"    Left Key = ",X%LEFT%KEY,ERR,ERROR,*999)
        ENDIF
        IF(ASSOCIATED(X%RIGHT,TREE%NIL)) THEN
          CALL WRITE_STRING(ID,"    Right Key = NIL",ERR,ERROR,*999)
        ELSE
          CALL WRITE_STRING_VALUE(ID,"    Right Key = ",X%RIGHT%KEY,ERR,ERROR,*999)
        ENDIF
        IF(ASSOCIATED(X%PARENT,TREE%NIL)) THEN
          CALL WRITE_STRING(ID,"    Parent Key = NIL",ERR,ERROR,*999)
        ELSE
          CALL WRITE_STRING_VALUE(ID,"    Parent Key = ",X%PARENT%KEY,ERR,ERROR,*999)
        ENDIF
        !Output the right subtree last
        CALL TREE_OUTPUT_IN_ORDER(ID,TREE,X%RIGHT,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_OUTPUT_IN_ORDER")
    RETURN
999 CALL ERRORS("TREE_OUTPUT_IN_ORDER",ERR,ERROR)
    CALL EXITS("TREE_OUTPUT_IN_ORDER")
    RETURN 1
  END SUBROUTINE TREE_OUTPUT_IN_ORDER

  !
  !================================================================================================================================
  !

  !>Returns the predeccessor of a tree at a specified tree node
  FUNCTION TREE_PREDECESSOR(TREE,X,ERR,ERROR)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the Red-Black tree to find the predecessor of
    TYPE(TREE_NODE_TYPE), POINTER :: X !<A pointer to the tree node to return the predecessor of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Function variable
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_PREDECESSOR !<On return the pointer to the predecessor of X or NIL if no predecessor exits
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: Y

    NULLIFY(TREE_PREDECESSOR)
    
    CALL ENTERS("TREE_PREDECESSOR",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(ASSOCIATED(X)) THEN
        Y=>X%LEFT
        IF(ASSOCIATED(Y,TREE%NIL)) THEN
          DO WHILE(.NOT.ASSOCIATED(Y%RIGHT,TREE%NIL))
            Y=>Y%RIGHT
          ENDDO
          TREE_PREDECESSOR=>Y
        ELSE
          Y=>X%PARENT
          DO WHILE(ASSOCIATED(X,Y%LEFT))
            IF(ASSOCIATED(Y,TREE%ROOT)) THEN
              TREE_PREDECESSOR=>TREE%NIL
              EXIT
            ELSE
              X=>Y
              Y=>Y%PARENT
            ENDIF
          ENDDO
          IF(.NOT.ASSOCIATED(TREE_PREDECESSOR)) TREE_PREDECESSOR=>Y
        ENDIF
       ELSE
        CALL FLAG_ERROR("Tree node X is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_PREDECESSOR")
    RETURN
999 CALL ERRORS("TREE_PREDECESSOR",ERR,ERROR)
    CALL EXITS("TREE_PREDECESSOR")
    RETURN 
  END FUNCTION TREE_PREDECESSOR

  !
  !================================================================================================================================
  !

  !>Searches a tree to see if it contains a key
  SUBROUTINE TREE_SEARCH(TREE,KEY,X,ERR,ERROR,*)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the tree to search
    INTEGER(INTG), INTENT(IN) :: KEY !<The key to search for
    TYPE(TREE_NODE_TYPE), POINTER :: X !<On return a pointer to the tree node containing the key. If the key does not exist NULL is returned
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    INTEGER(INTG) :: COMPARE_VALUE
    TYPE(TREE_NODE_TYPE), POINTER :: Y
    
    CALL ENTERS("TREE_SEARCH",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(TREE%TREE_FINISHED) THEN
        IF(ASSOCIATED(X)) THEN
          CALL FLAG_ERROR("The tree node X is already associated",ERR,ERROR,*999)
        ELSE
          NULLIFY(X)
          Y=>TREE%ROOT
          IF(.NOT.ASSOCIATED(Y,TREE%NIL)) THEN
            COMPARE_VALUE=Y%KEY-KEY
            DO WHILE(COMPARE_VALUE/=0)
              IF(COMPARE_VALUE>0) THEN !Y%KEY > KEY
                Y=>Y%LEFT
              ELSE !Y%KEY < KEY
                Y=>Y%RIGHT
              ENDIF
              IF(ASSOCIATED(Y,TREE%NIL)) THEN
                EXIT
              ELSE
                COMPARE_VALUE=Y%KEY-KEY
              ENDIF
            ENDDO
            IF(COMPARE_VALUE==0) X=>Y
          ENDIF          
        ENDIF
      ELSE
        CALL FLAG_ERROR("The tree has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_SEARCH")
    RETURN
999 CALL ERRORS("TREE_SEARCH",ERR,ERROR)
    CALL EXITS("TREE_SEARCH")
    RETURN 1
  END SUBROUTINE TREE_SEARCH

  !
  !================================================================================================================================
  !

  !>Returns the successor of a tree at a specified tree node
  FUNCTION TREE_SUCCESSOR(TREE,X,ERR,ERROR)

    !Argument Variables
    TYPE(TREE_TYPE), POINTER :: TREE !<A pointer to the Red-Black tree to find the successor of
    TYPE(TREE_NODE_TYPE), POINTER :: X !<A pointer to the tree node to return the successor of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Function variable
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_SUCCESSOR !<On return the pointer to the successor of X or NIL if no sucessor exits
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: Y

    NULLIFY(TREE_SUCCESSOR)
    
    CALL ENTERS("TREE_SUCCESSOR",ERR,ERROR,*999)

    IF(ASSOCIATED(TREE)) THEN
      IF(ASSOCIATED(X)) THEN
        Y=>X%RIGHT
        IF(ASSOCIATED(Y,TREE%NIL)) THEN
          DO WHILE(.NOT.ASSOCIATED(Y%LEFT,TREE%NIL))
            Y=>Y%LEFT
          ENDDO
          TREE_SUCCESSOR=>Y
          RETURN
        ELSE
          Y=>X%PARENT
          DO WHILE(ASSOCIATED(X,Y%RIGHT))
            X=>Y
            Y=>Y%PARENT
          ENDDO
          IF(ASSOCIATED(Y,TREE%ROOT)) THEN
            TREE_SUCCESSOR=>TREE%NIL
          ELSE
            TREE_SUCCESSOR=>Y
          ENDIF
        ENDIF
       ELSE
        CALL FLAG_ERROR("Tree node X is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Tree is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("TREE_SUCCESSOR")
    RETURN
999 CALL ERRORS("TREE_SUCCESSOR",ERR,ERROR)
    CALL EXITS("TREE_SUCCESSOR")
    RETURN 
  END FUNCTION TREE_SUCCESSOR

  !
  !================================================================================================================================
  !

END MODULE TREES
