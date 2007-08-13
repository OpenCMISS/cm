!> \file
!> $Id: lists.f90 28 2007-07-27 08:35:14Z cpb $
!> \author Chris Bradley
!> \brief Implements lists of base types.
!> \todo Fix up and have this module use the sorting module for sorts
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

!> Implements lists of base types.
MODULE LISTS

  USE BASE_ROUTINES
  USE CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup LISTS_DataType LISTS::DataType
  !> \brief Data type parameters for a list.
  !> \see LISTS
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_INTG_TYPE=INTEGER_TYPE !<Integer data type for a list \see LISTS_DataType,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SP_TYPE=SINGLE_REAL_TYPE !<Single precision real data type for a list \see LISTS_DataType,LISTS
  INTEGER(INTG), PARAMETER :: LIST_DP_TYPE=DOUBLE_REAL_TYPE !<Double precision real data type for a list \see LISTS_DataType,LISTS
  !>@}
  
  !> \addtogroup LISTS_SortingOrder LISTS::SortingOrder
  !> \brief Sorting order parameters for a list.
  !> \see LISTS
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_UNSORTED_TYPE=1 !<Unsorted list type \see LISTS_SortingOrder,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SORT_ASCENDING_TYPE=2 !<Ascending order for sort \see LISTS_SortingOrder,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SORT_DESCENDING_TYPE=3 !<Descending order for sort \see LISTS_SortingOrder,LISTS
  !>@}

  !> \addtogroup LISTSSortingMethod LISTS::SortingMethod
  !> \brief Sorting method parameters for a list.
  !> \see LISTS
  !>@{
  INTEGER(INTG), PARAMETER :: LIST_BUBBLE_SORT_METHOD=1 !<Bubble sort method \see LISTS_SortingMethod,LISTS
  INTEGER(INTG), PARAMETER :: LIST_SHELL_SORT_METHOD=2 !<Shell sort method \see LISTS_SortingMethod,LISTS
  INTEGER(INTG), PARAMETER :: LIST_HEAP_SORT_METHOD=3 !<Heap sort method \see LISTS_SortingMethod,LISTS
  !>@}

  !Module types

  !>Buffer type to allow arrays of pointers to a list
  TYPE LIST_PTR_TYPE
    TYPE(LIST_TYPE), POINTER :: PTR !<The pointer to the list
  END TYPE LIST_PTR_TYPE

  !>Contains information on a list
  TYPE LIST_TYPE
    LOGICAL :: LIST_FINISHED !<Is .TRUE. if the list has finished being created, .FALSE. if not.
    INTEGER(INTG) :: NUMBER_IN_LIST !<The number of items currently in the list
    INTEGER(INTG) :: INITIAL_SIZE !<The size of the list when it was initially created.
    INTEGER(INTG) :: SIZE !<The current size of the list.
    INTEGER(INTG) :: DATA_TYPE !<The data type of the list \see LISTS_DataType
    INTEGER(INTG) :: SORT_ORDER !<The ordering to be used when sorting the list \see LISTS_SortingOrder
    INTEGER(INTG) :: SORT_METHOD !<The sorting method to be used when sorting the list \see LISTS_SortingMethod
    INTEGER(INTG), POINTER :: LIST_INTG(:) !<A pointer to the integer data for integer lists. 
    REAL(SP), POINTER :: LIST_SP(:) !<A pointer to the single precision data for single precision real lists. 
    REAL(DP), POINTER :: LIST_DP(:) !<A pointer to the double precision data for double precision real lists. 
  END TYPE LIST_TYPE
    
  !Module variables
  
  !Interfaces

  !>Adds an item to the end of a list \see LISTS.
  INTERFACE LIST_ITEM_ADD
    MODULE PROCEDURE LIST_ITEM_ADD_INTG1
    MODULE PROCEDURE LIST_ITEM_ADD_SP1
    MODULE PROCEDURE LIST_ITEM_ADD_DP1
  END INTERFACE !LIST_ITEM_ADD
  
 !>Determines if an item is in a list and returns the position of the item \see LISTS.
 INTERFACE LIST_ITEM_IN_LIST
    MODULE PROCEDURE LIST_ITEM_IN_LIST_INTG1
    MODULE PROCEDURE LIST_ITEM_IN_LIST_SP1
    MODULE PROCEDURE LIST_ITEM_IN_LIST_DP1
  END INTERFACE !LIST_ITEM_IN_LIST

  !>Detaches the list values from a list and returns them as a pointer to a array of base type \see LISTS.
  INTERFACE LIST_DETACH
    MODULE PROCEDURE LIST_DETACH_INTG
    MODULE PROCEDURE LIST_DETACH_SP
    MODULE PROCEDURE LIST_DETACH_DP
  END INTERFACE !LIST_DETACH

  !>Searches a list for a given value and returns the position in the list if the value exists \see LISTS.
  INTERFACE LIST_SEARCH
    MODULE PROCEDURE LIST_SEARCH_INTG_ARRAY
    MODULE PROCEDURE LIST_SEARCH_SP_ARRAY
    MODULE PROCEDURE LIST_SEARCH_DP_ARRAY
  END INTERFACE !LIST_SEARCH

  INTERFACE LIST_SEARCH_LINEAR
    MODULE PROCEDURE LIST_SEARCH_LINEAR_INTG_ARRAY
    MODULE PROCEDURE LIST_SEARCH_LINEAR_SP_ARRAY
    MODULE PROCEDURE LIST_SEARCH_LINEAR_DP_ARRAY
  END INTERFACE !LIST_SEARCH_LINEAR

  INTERFACE LIST_SORT
    MODULE PROCEDURE LIST_SORT_INTG_ARRAY
    MODULE PROCEDURE LIST_SORT_SP_ARRAY
    MODULE PROCEDURE LIST_SORT_DP_ARRAY
  END INTERFACE !LIST_SORT

  INTERFACE LIST_SORT_BUBBLE
    MODULE PROCEDURE LIST_SORT_BUBBLE_INTG_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_SP_ARRAY
    MODULE PROCEDURE LIST_SORT_BUBBLE_DP_ARRAY
  END INTERFACE !LIST_SORT_BUBBLE

  INTERFACE LIST_SORT_HEAP
    MODULE PROCEDURE LIST_SORT_HEAP_INTG_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_SP_ARRAY
    MODULE PROCEDURE LIST_SORT_HEAP_DP_ARRAY
  END INTERFACE !LIST_SORT_HEAP

  INTERFACE LIST_SORT_SHELL
    MODULE PROCEDURE LIST_SORT_SHELL_INTG_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_SP_ARRAY
    MODULE PROCEDURE LIST_SORT_SHELL_DP_ARRAY
  END INTERFACE !LIST_SORT_SHELL

  PUBLIC LIST_TYPE,LIST_PTR_TYPE

  PUBLIC LIST_INTG_TYPE,LIST_SP_TYPE,LIST_DP_TYPE

  PUBLIC LIST_CREATE_FINISH,LIST_CREATE_START,LIST_DATA_TYPE_SET,LIST_DESTROY,LIST_INITIAL_SIZE_SET,LIST_ITEM_ADD, &
    & LIST_ITEM_DELETE,LIST_DETACH,LIST_REMOVE_DUPLICATES

  PUBLIC LIST_SEARCH,LIST_SEARCH_LINEAR
  
  PUBLIC LIST_SORT,LIST_SORT_BUBBLE,LIST_SORT_HEAP,LIST_SORT_SHELL
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a list created with LIST_CREATE_START \see{LISTS::LIST_CREATE_START}.
  SUBROUTINE LIST_CREATE_FINISH(LIST,ERR,ERROR,*)

    !#### Subroutine: LIST_CREATE_FINISH
    !###  Description:
    !###    Finishes the creation of a list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list to finish
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LIST_CREATE_FINISH",ERR,ERROR,*998)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        CALL FLAG_ERROR("List is already finished",ERR,ERROR,*998)
      ELSE
        !Allocate the list
        SELECT CASE(LIST%DATA_TYPE)
        CASE(LIST_INTG_TYPE)
          ALLOCATE(LIST%LIST_INTG(LIST%INITIAL_SIZE),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list integer data",ERR,ERROR,*999)
        CASE(LIST_SP_TYPE)
          ALLOCATE(LIST%LIST_SP(LIST%INITIAL_SIZE),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list single precision data",ERR,ERROR,*999)
        CASE(LIST_DP_TYPE)
          ALLOCATE(LIST%LIST_DP(LIST%INITIAL_SIZE),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list double precision data",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        LIST%LIST_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("LIST_CREATE_FINISH")
    RETURN
999 IF(ASSOCIATED(LIST%LIST_INTG)) DEALLOCATE(LIST%LIST_INTG)
    IF(ASSOCIATED(LIST%LIST_SP)) DEALLOCATE(LIST%LIST_SP)
    IF(ASSOCIATED(LIST%LIST_DP)) DEALLOCATE(LIST%LIST_DP)
998 CALL ERRORS("LIST_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("LIST_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE LIST_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a list and returns a pointer to the created list \see{LISTS::LIST_CREATE_FINISH}.
  SUBROUTINE LIST_CREATE_START(LIST,ERR,ERROR,*)

    !#### Subroutine: LIST_CREATE_START
    !###  Description:
    !###    Starts the creation of a list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list to create. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables

    CALL ENTERS("LIST_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(LIST)) THEN
      CALL FLAG_ERROR("List is already associated",ERR,ERROR,*998)
    ELSE
      ALLOCATE(LIST,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate list",ERR,ERROR,*999)
      CALL LIST_INITIALISE(LIST,ERR,ERROR,*999)
      !Set defaults
      LIST%DATA_TYPE=LIST_INTG_TYPE
      LIST%INITIAL_SIZE=10
      LIST%SORT_ORDER=LIST_SORT_ASCENDING_TYPE
      LIST%SORT_METHOD=LIST_HEAP_SORT_METHOD
    ENDIF

    CALL EXITS("LIST_CREATE_START")
    RETURN
999 CALL LIST_FINALISE(LIST,ERR,ERROR,*998)
998 CALL ERRORS("LIST_CREATE_START",ERR,ERROR)
    CALL EXITS("LIST_CREATE_START")
    RETURN 1
  END SUBROUTINE LIST_CREATE_START

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a list.
  SUBROUTINE LIST_DATA_TYPE_SET(LIST,DATA_TYPE,ERR,ERROR,*)

    !#### Subroutine: LIST_DATA_TYPE_SET
    !###  Description:
    !###    Sets/changes the data type for a list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list 
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type of the list to set. \see{LISTSDataType}
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_DATA_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        CALL FLAG_ERROR("List has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DATA_TYPE)
        CASE(LIST_INTG_TYPE)
          LIST%DATA_TYPE=LIST_INTG_TYPE
        CASE(LIST_SP_TYPE)
          LIST%DATA_TYPE=LIST_SP_TYPE
        CASE(LIST_DP_TYPE)
          LIST%DATA_TYPE=LIST_DP_TYPE
        CASE DEFAULT
          LOCAL_ERROR="The data type of "//TRIM(NUMBER_TO_VSTRING(DATA_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_DATA_TYPE_SET")
    RETURN
999 CALL ERRORS("LIST_DATA_TYPE_SET",ERR,ERROR)
    CALL EXITS("LIST_DATA_TYPE_SET")
    RETURN 1
  END SUBROUTINE LIST_DATA_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Destroys a list.
  SUBROUTINE LIST_DESTROY(LIST,ERR,ERROR,*)

    !#### Subroutine: LIST_DESTROY
    !###  Description:
    !###    Destroys a list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("LIST_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      CALL LIST_FINALISE(LIST,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_DESTROY")
    RETURN
999 CALL ERRORS("LIST_DESTROY",ERR,ERROR)
    CALL EXITS("LIST_DESTROY")
    RETURN 1
  END SUBROUTINE LIST_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a list and deallocates all memory.
  SUBROUTINE LIST_FINALISE(LIST,ERR,ERROR,*)    

    !#### Subroutine: LIST_FINALISE
    !###  Description:
    !###    Finalises a list and deallocates all memory.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("LIST_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(ASSOCIATED(LIST%LIST_INTG)) DEALLOCATE(LIST%LIST_INTG)
      IF(ASSOCIATED(LIST%LIST_SP)) DEALLOCATE(LIST%LIST_SP)
      IF(ASSOCIATED(LIST%LIST_DP)) DEALLOCATE(LIST%LIST_DP)
      DEALLOCATE(LIST)
    ENDIF

    CALL EXITS("LIST_FINALISE")
    RETURN
999 CALL ERRORS("LIST_FINALISE",ERR,ERROR)
    CALL EXITS("LIST_FINALISE")
    RETURN 1
  END SUBROUTINE LIST_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a list and all its components
  SUBROUTINE LIST_INITIALISE(LIST,ERR,ERROR,*)

    !#### Subroutine: LIST_INITIALISE
    !###  Description:
    !###    Initialises a list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("LIST_INTIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      LIST%LIST_FINISHED=.FALSE.
      LIST%NUMBER_IN_LIST=0
      LIST%INITIAL_SIZE=0
      LIST%SIZE=0
      LIST%DATA_TYPE=0
      LIST%SORT_ORDER=0
      LIST%SORT_METHOD=0
      NULLIFY(LIST%LIST_INTG)
      NULLIFY(LIST%LIST_SP)
      NULLIFY(LIST%LIST_DP)
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_INTIALISE")
    RETURN
999 CALL ERRORS("LIST_INTIALISE",ERR,ERROR)
    CALL EXITS("LIST_INITIALISE")
    RETURN 1
  END SUBROUTINE LIST_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the initial size for a list
  SUBROUTINE LIST_INITIAL_SIZE_SET(LIST,INITIAL_SIZE,ERR,ERROR,*)

    !#### Subroutine: LIST_INITIAL_SIZE_SET
    !###  Description:
    !###    Sets/changes the initial size for a list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: INITIAL_SIZE !<The initial size of the list to set. Must be greater than zero.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_INTIAL_SIZE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        CALL FLAG_ERROR("List has been finished",ERR,ERROR,*999)
      ELSE
        IF(INITIAL_SIZE>0) THEN
          LIST%INITIAL_SIZE=INITIAL_SIZE
        ELSE
          LOCAL_ERROR="The initial size of "//TRIM(NUMBER_TO_VSTRING(INITIAL_SIZE,"*",ERR,ERROR))// &
            & " is invalid. The size must be > 0"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_INTIAL_SIZE_SET")
    RETURN
999 CALL ERRORS("LIST_INTIAL_SIZE_SET",ERR,ERROR)
    CALL EXITS("LIST_INITIAL_SIZE_SET")
    RETURN 1
  END SUBROUTINE LIST_INITIAL_SIZE_SET

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_ITEM_ADD
  !###  Description:
  !###    Adds an item to a list.
  !###  Child-subroutines: LIST_ITEM_ADD_INTG1,LIST_ITEM_ADD_SP1,LIST_ITEM_ADD_DP1

  !
  !================================================================================================================================
  !

  !>Adds an item to the end of an integer list. 
  SUBROUTINE LIST_ITEM_ADD_INTG1(LIST,ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_ADD_INTG1
    !###  Description:
    !###    Adds the item ITEM to an integer list.
    !###  Parent-routine: LIST_ITEM_ADD

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list
    INTEGER(INTG), INTENT(IN) :: ITEM !<The item to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    INTEGER(INTG), POINTER :: NEW_LIST(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_LIST)
    
    CALL ENTERS("LIST_ITEM_ADD_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%DATA_TYPE==LIST_INTG_TYPE) THEN
          IF(LIST%NUMBER_IN_LIST==LIST%SIZE) THEN
            !Reallocate
            NEW_SIZE=MAX(2*LIST%NUMBER_IN_LIST,1)
            ALLOCATE(NEW_LIST(NEW_SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new list",ERR,ERROR,*999)
            NEW_LIST(1:LIST%NUMBER_IN_LIST)=LIST%LIST_INTG(1:LIST%NUMBER_IN_LIST)
            IF(ASSOCIATED(LIST%LIST_INTG)) DEALLOCATE(LIST%LIST_INTG)
            LIST%LIST_INTG=>NEW_LIST
            LIST%SIZE=NEW_SIZE
          ENDIF
          LIST%LIST_INTG(LIST%NUMBER_IN_LIST+1)=ITEM
          LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST+1
        ELSE
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not match the integer type of the supplied list item"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The list has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("LIST_ITEM_ADD_INTG1")
    RETURN
999 IF(ASSOCIATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL ERRORS("LIST_ITEM_ADD_INTG1",ERR,ERROR)
    CALL EXITS("LIST_ITEM_ADD_INTG1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_INTG1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a single precision real list. 
  SUBROUTINE LIST_ITEM_ADD_SP1(LIST,ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_ADD_SP1
    !###  Description:
    !###    Adds the item ITEM to a single precision list.
    !###  Parent-routine: LIST_ITEM_ADD

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list
    REAL(SP), INTENT(IN) :: ITEM !<The item to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    REAL(SP), POINTER :: NEW_LIST(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    NULLIFY(NEW_LIST)
    
    CALL ENTERS("LIST_ITEM_ADD_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%DATA_TYPE==LIST_SP_TYPE) THEN
          IF(LIST%NUMBER_IN_LIST==LIST%SIZE) THEN
            !Reallocate
            NEW_SIZE=MAX(2*LIST%NUMBER_IN_LIST,1)
            ALLOCATE(NEW_LIST(NEW_SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new list",ERR,ERROR,*999)
            NEW_LIST(1:LIST%NUMBER_IN_LIST)=LIST%LIST_SP(1:LIST%NUMBER_IN_LIST)
            IF(ASSOCIATED(LIST%LIST_SP)) DEALLOCATE(LIST%LIST_SP)
            LIST%LIST_SP=>NEW_LIST
            LIST%SIZE=NEW_SIZE
          ENDIF
          LIST%LIST_SP(LIST%NUMBER_IN_LIST+1)=ITEM
          LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST+1
        ELSE
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not match the single precision type of the supplied list item"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The list has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("LIST_ITEM_ADD_SP1")
    RETURN
999 IF(ASSOCIATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL ERRORS("LIST_ITEM_ADD_SP1",ERR,ERROR)
    CALL EXITS("LIST_ITEM_ADD_SP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_SP1
  
  !
  !================================================================================================================================
  !

  !>Adds an item to the end of a double precision real list.
  SUBROUTINE LIST_ITEM_ADD_DP1(LIST,ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_ADD_DP1
    !###  Description:
    !###    Adds the item ITEM to a double precision list.
    !###  Parent-routine: LIST_ITEM_ADD

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<A pointer to the list
    REAL(DP), INTENT(IN) :: ITEM !<The item to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_SIZE
    REAL(DP), POINTER :: NEW_LIST(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_LIST)
    
    CALL ENTERS("LIST_ITEM_ADD_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%DATA_TYPE==LIST_DP_TYPE) THEN
          IF(LIST%NUMBER_IN_LIST==LIST%SIZE) THEN
            !Reallocate
            NEW_SIZE=MAX(2*LIST%NUMBER_IN_LIST,1)
            ALLOCATE(NEW_LIST(NEW_SIZE),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new list",ERR,ERROR,*999)
            NEW_LIST(1:LIST%NUMBER_IN_LIST)=LIST%LIST_DP(1:LIST%NUMBER_IN_LIST)
            IF(ASSOCIATED(LIST%LIST_DP)) DEALLOCATE(LIST%LIST_DP)
            LIST%LIST_DP=>NEW_LIST
            LIST%SIZE=NEW_SIZE
          ENDIF
          LIST%LIST_DP(LIST%NUMBER_IN_LIST+1)=ITEM
          LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST+1
        ELSE
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not match the double precision type of the supplied list item"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The list has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
    CALL EXITS("LIST_ITEM_ADD_DP1")
    RETURN
999 IF(ASSOCIATED(NEW_LIST)) DEALLOCATE(NEW_LIST)
    CALL ERRORS("LIST_ITEM_ADD_DP1",ERR,ERROR)
    CALL EXITS("LIST_ITEM_ADD_DP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_ADD_DP1
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_ITEM_IN_LIST
  !###  Description:
  !###    Determines if an item is in a list and returns the position of the item.
  !###  Child-subroutines: LIST_ITEM_IN_LIST_INTG1,LIST_ITEM_IN_LIST_SP1,LIST_ITEM_IN_LIST_DP1

  !
  !================================================================================================================================
  !

  !>Determines if ITEM is in the given integer LIST. If it is LIST_ITEM is the index in the list. If not LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_INTG1(LIST,ITEM,LIST_ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_IN_LIST_INTG1
    !###  Description:
    !###    Determines if ITEM is in the given integer LIST. If it is LIST_ITEM is the index in the list. If not LIST_ITEM is 0.
    !###  Parent-routine: LIST_ITEM_IN_LIST

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: ITEM !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string.
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LIST_ITEM_IN_LIST_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%DATA_TYPE==LIST_INTG_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          CALL LIST_SEARCH_LINEAR(LIST%LIST_INTG(1:LIST%NUMBER_IN_LIST),ITEM,LIST_ITEM,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not match the integer type of the supplied list item"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_ITEM_IN_LIST_INTG1")
    RETURN
999 CALL ERRORS("LIST_ITEM_IN_LIST_INTG1",ERR,ERROR)
    CALL EXITS("LIST_ITEM_IN_LIST_INTG1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_INTG1
  
  !
  !================================================================================================================================
  !

  !> Determines if ITEM is in the given single precision real LIST. If it is LIST_ITEM is the index in the list. If not
  !> LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_SP1(LIST,ITEM,LIST_ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_IN_LIST_SP1
    !###  Description:
    !###    Determines if ITEM is in the given single precision LIST. If it is LIST_ITEM is the index in the list. If not
    !###    LIST_ITEM is 0.
    !###  Parent-routine: LIST_ITEM_IN_LIST

    !Argument Variables    
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    REAL(SP), INTENT(IN) :: ITEM !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.     
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code    
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LIST_ITEM_IN_LIST_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%DATA_TYPE==LIST_SP_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          CALL LIST_SEARCH_LINEAR(LIST%LIST_SP(1:LIST%NUMBER_IN_LIST),ITEM,LIST_ITEM,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not match the single precision type of the supplied list item"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_ITEM_IN_LIST_SP1")
    RETURN
999 CALL ERRORS("LIST_ITEM_IN_LIST_SP1",ERR,ERROR)
    CALL EXITS("LIST_ITEM_IN_LIST_SP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_SP1
  
  !
  !================================================================================================================================
  !

  !> Determines if ITEM is in the given double precision real LIST. If it is LIST_ITEM is the index in the list. If not
  !> LIST_ITEM is 0.
  SUBROUTINE LIST_ITEM_IN_LIST_DP1(LIST,ITEM,LIST_ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_IN_LIST_DP1
    !###  Description:
    !###    Determines if ITEM is in the given double precision LIST. If it is LIST_ITEM is the index in the list. If not
    !###    LIST_ITEM is 0.
    !###  Parent-routine: LIST_ITEM_IN_LIST

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    REAL(DP), INTENT(IN) :: ITEM  !<The item to find.
    INTEGER(INTG), INTENT(OUT) :: LIST_ITEM !<On exit, the position of the item in the list. If the item does not exist then the value of 0 is returned.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LIST_ITEM_IN_LIST_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%DATA_TYPE==LIST_DP_TYPE) THEN
!!TODO: Could search better but requires list to be sorted.
          CALL LIST_SEARCH_LINEAR(LIST%LIST_DP(1:LIST%NUMBER_IN_LIST),ITEM,LIST_ITEM,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
            & " does not match the single precision type of the supplied list item"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_ITEM_IN_LIST_DP1")
    RETURN
999 CALL ERRORS("LIST_ITEM_IN_LIST_DP1",ERR,ERROR)
    CALL EXITS("LIST_ITEM_IN_LIST_DP1")
    RETURN 1
  END SUBROUTINE LIST_ITEM_IN_LIST_DP1

  !
  !================================================================================================================================
  !
  
  !>Deletes the item given by the LIST_ITEM index from the given list.
  SUBROUTINE LIST_ITEM_DELETE(LIST,LIST_ITEM,ERR,ERROR,*)

    !#### Subroutine: LIST_ITEM_DELETE
    !###  Description:
    !###    Deletes the list item given by the list index (or pointer) from the list.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(IN) :: LIST_ITEM !<The position in the list to delete.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_ITEM_DELETE",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST_ITEM>=1.AND.LIST_ITEM<=LIST%NUMBER_IN_LIST) THEN
          SELECT CASE(LIST%DATA_TYPE)
          CASE(LIST_INTG_TYPE)
            LIST%LIST_INTG(LIST_ITEM:LIST%NUMBER_IN_LIST-1)=LIST%LIST_INTG(LIST_ITEM+1:LIST%NUMBER_IN_LIST)
          CASE(LIST_SP_TYPE)
            LIST%LIST_SP(LIST_ITEM:LIST%NUMBER_IN_LIST-1)=LIST%LIST_SP(LIST_ITEM+1:LIST%NUMBER_IN_LIST)
          CASE(LIST_DP_TYPE)
            LIST%LIST_DP(LIST_ITEM:LIST%NUMBER_IN_LIST-1)=LIST%LIST_DP(LIST_ITEM+1:LIST%NUMBER_IN_LIST)
          CASE DEFAULT
            LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST-1
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("LIST_ITEM_DELETE")
    RETURN
999 CALL ERRORS("LIST_ITEM_DELETE",ERR,ERROR)
    CALL EXITS("LIST_ITEM_DELETE")
    RETURN 1
  END SUBROUTINE LIST_ITEM_DELETE
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_DETACH
  !###  Description:
  !###    Detaches the list values from a list and returns them as a pointer to a array of base type. The LIST_VALUES pointer
  !###    must be NULL on entry. It is up to the user to then deallocate the returned list memory.
  !###  Child-subroutines: LIST_DETACH_INTG,LIST_DETACH_SP,LIST_DETACH_DP

  !
  !================================================================================================================================
  !

  !>Detaches the list values from an integer list and returns them as a pointer to a array of base type.
  !>The LIST_VALUES pointer must not be associated on entry. It is up to the user to then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_INTG(LIST,NUMBER_IN_LIST,LIST_VALUES,ERR,ERROR,*)

    !#### Subroutine: LIST_DETACH_INTG
    !###  Description:
    !###    Detaches the list values from a integer list and returns them as a pointer to a array of base type.
    !###    The LIST_VALUES pointer must be NULL on entry. It is up to the user to then deallocate the returned list memory.
    !###  Parent-routine: LIST_DETACH

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    INTEGER(INTG), POINTER :: LIST_VALUES(:) !<On exit, a pointer to the detached list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_DETACH_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(ASSOCIATED(LIST_VALUES)) THEN
          CALL FLAG_ERROR("List values is associated",ERR,ERROR,*999)
        ELSE
          IF(LIST%DATA_TYPE==LIST_INTG_TYPE) THEN
            NUMBER_IN_LIST=LIST%NUMBER_IN_LIST
            !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
            LIST_VALUES=>LIST%LIST_INTG
            LIST%NUMBER_IN_LIST=0
            NULLIFY(LIST%LIST_INTG)
          ELSE
            LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not match the integer type of the supplied list values item"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("LIST_DETACH_INTG")
    RETURN
999 CALL ERRORS("LIST_DETACH_INTG",ERR,ERROR)
    CALL EXITS("LIST_DETACH_INTG")
    RETURN 1
  END SUBROUTINE LIST_DETACH_INTG

  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a single precision real list and returns them as a pointer to a array of base type.
  !>The LIST_VALUES pointer must not be associated on entry. It is up to the user to then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_SP(LIST,NUMBER_IN_LIST,LIST_VALUES,ERR,ERROR,*)

    !#### Subroutine: LIST_DETACH_SP
    !###  Description:
    !###    Detaches the list values from a single precision list and returns them as a pointer to a array of base type.
    !###    The LIST_VALUES pointer must be NULL on entry. It is up to the user to then deallocate the returned list memory.
    !###  Parent-routine: LIST_DETACH

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    REAL(SP), POINTER :: LIST_VALUES(:) !<On exit, a pointer to the detached list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_DETACH_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(ASSOCIATED(LIST_VALUES)) THEN
          CALL FLAG_ERROR("List values is associated",ERR,ERROR,*999)
        ELSE
          IF(LIST%DATA_TYPE==LIST_SP_TYPE) THEN
            NUMBER_IN_LIST=LIST%NUMBER_IN_LIST
            !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
            LIST_VALUES=>LIST%LIST_SP
            LIST%NUMBER_IN_LIST=0
            NULLIFY(LIST%LIST_SP)
          ELSE
            LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not match the single precision type of the supplied list values item"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("LIST_DETACH_SP")
    RETURN
999 CALL ERRORS("LIST_DETACH_SP",ERR,ERROR)
    CALL EXITS("LIST_DETACH_SP")
    RETURN 1
  END SUBROUTINE LIST_DETACH_SP
  !
  !================================================================================================================================
  !
  
  !>Detaches the list values from a double precision real list and returns them as a pointer to a array of base type.
  !>The LIST_VALUES pointer must not be associated on entry. It is up to the user to then deallocate the returned list memory.
  SUBROUTINE LIST_DETACH_DP(LIST,NUMBER_IN_LIST,LIST_VALUES,ERR,ERROR,*)

    !#### Subroutine: LIST_DETACH_DP
    !###  Description:
    !###    Detaches the list values from a double precision list and returns them as a pointer to a array of base type.
    !###    The LIST_VALUES pointer must be NULL on entry. It is up to the user to then deallocate the returned list memory.
    !###  Parent-routine: LIST_DETACH

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: NUMBER_IN_LIST !<On exit, the number in the list that has been detached.
    REAL(DP), POINTER :: LIST_VALUES(:) !<On exit, a pointer to the detached list. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_DETACH_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(ASSOCIATED(LIST_VALUES)) THEN
          CALL FLAG_ERROR("List values is associated",ERR,ERROR,*999)
        ELSE
          IF(LIST%DATA_TYPE==LIST_DP_TYPE) THEN
            NUMBER_IN_LIST=LIST%NUMBER_IN_LIST
            !Note this will return more memory as the list will be bigger. Maybe copy to an array the correct size?
            LIST_VALUES=>LIST%LIST_DP
            LIST%NUMBER_IN_LIST=0
            NULLIFY(LIST%LIST_DP)
          ELSE
            LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))// &
              & " does not match the double precision type of the supplied list values item"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("LIST_DETACH_DP")
    RETURN
999 CALL ERRORS("LIST_DETACH_DP",ERR,ERROR)
    CALL EXITS("LIST_DETACH_DP")
    RETURN 1
  END SUBROUTINE LIST_DETACH_DP

  !
  !================================================================================================================================
  !

  !>Removes duplicate entries from a list. A side effect of this is that the list is sorted.
  SUBROUTINE LIST_REMOVE_DUPLICATES(LIST,ERR,ERROR,*)

    !#### Subroutine: LIST_REMOVE_DUPLICATES
    !###  Description:
    !###    Removes duplicate entries from a list. A side effect of this is that the list is sorted.

    !Argument Variables
    TYPE(LIST_TYPE), POINTER :: LIST !<The pointer to the list
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j,NUMBER_REMOVED
    LOGICAL :: SAME_VALUE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("LIST_REMOVE_DUPLICATES",ERR,ERROR,*999)

    IF(ASSOCIATED(LIST)) THEN
      IF(LIST%LIST_FINISHED) THEN
        IF(LIST%NUMBER_IN_LIST>0) THEN
          SELECT CASE(LIST%DATA_TYPE)
          CASE(LIST_INTG_TYPE)
            CALL LIST_SORT(LIST%LIST_INTG(1:LIST%NUMBER_IN_LIST),ERR,ERROR,*999)
            i=1
            DO WHILE(i<=LIST%NUMBER_IN_LIST)
              !Find the extent of duplicate values if any
              j=i+1
              SAME_VALUE=.TRUE.
              DO WHILE(j<=LIST%NUMBER_IN_LIST.AND.SAME_VALUE)
                IF(LIST%LIST_INTG(j)==LIST%LIST_INTG(i)) THEN
                  j=j+1
                ELSE
                  SAME_VALUE=.FALSE.
                ENDIF
              ENDDO !j
              IF(j>i+1.OR.SAME_VALUE) THEN
                !We have duplicates so remove them
                IF(SAME_VALUE) THEN
                  !Duplicates to the end of the list so just set the number in the list
                  LIST%NUMBER_IN_LIST=i
                ELSE
                  NUMBER_REMOVED=j-i-1
                  LIST%LIST_INTG(i+1:LIST%NUMBER_IN_LIST-NUMBER_REMOVED)=LIST%LIST_INTG(j:LIST%NUMBER_IN_LIST)
                  LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST-NUMBER_REMOVED
                ENDIF
              ENDIF
              i=i+1
            ENDDO !i
          CASE(LIST_SP_TYPE)
           CALL LIST_SORT(LIST%LIST_SP(1:LIST%NUMBER_IN_LIST),ERR,ERROR,*999)
            i=1
            DO WHILE(i<=LIST%NUMBER_IN_LIST)
              !Find the extent of duplicate values if any
              j=i+1
              SAME_VALUE=.TRUE.
              DO WHILE(j<=LIST%NUMBER_IN_LIST.AND.SAME_VALUE)
                IF(ABS(LIST%LIST_SP(j)-LIST%LIST_SP(i))<=ZERO_TOLERANCE_SP) THEN
                  j=j+1
                ELSE
                  SAME_VALUE=.FALSE.
                ENDIF
              ENDDO !j
              IF(j>i+1.OR.SAME_VALUE) THEN
                !We have duplicates so remove them
                IF(SAME_VALUE) THEN
                  !Duplicates to the end of the list so just set the number in the list
                  LIST%NUMBER_IN_LIST=i
                ELSE
                  NUMBER_REMOVED=j-i-1
                  LIST%LIST_SP(i+1:LIST%NUMBER_IN_LIST-NUMBER_REMOVED)=LIST%LIST_SP(j:LIST%NUMBER_IN_LIST)
                  LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST-NUMBER_REMOVED
                ENDIF
              ENDIF
              i=i+1
            ENDDO !i
          CASE(LIST_DP_TYPE)
           CALL LIST_SORT(LIST%LIST_DP(1:LIST%NUMBER_IN_LIST),ERR,ERROR,*999)
            i=1
            DO WHILE(i<=LIST%NUMBER_IN_LIST)
              !Find the extent of duplicate values if any
              j=i+1
              SAME_VALUE=.TRUE.
              DO WHILE(j<=LIST%NUMBER_IN_LIST.AND.SAME_VALUE)
                IF(ABS(LIST%LIST_DP(j)-LIST%LIST_DP(i))<=ZERO_TOLERANCE) THEN
                  j=j+1
                ELSE
                  SAME_VALUE=.FALSE.
                ENDIF
              ENDDO !j
              IF(j>i+1.OR.SAME_VALUE) THEN
                !We have duplicates so remove them
                IF(SAME_VALUE) THEN
                  !Duplicates to the end of the list so just set the number in the list
                  LIST%NUMBER_IN_LIST=i
                ELSE
                  NUMBER_REMOVED=j-i-1
                  LIST%LIST_DP(i+1:LIST%NUMBER_IN_LIST-NUMBER_REMOVED)=LIST%LIST_DP(j:LIST%NUMBER_IN_LIST)
                  LIST%NUMBER_IN_LIST=LIST%NUMBER_IN_LIST-NUMBER_REMOVED
                ENDIF
              ENDIF
              i=i+1
            ENDDO !i
          CASE DEFAULT
            LOCAL_ERROR="The list data type of "//TRIM(NUMBER_TO_VSTRING(LIST%DATA_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ELSE
        CALL FLAG_ERROR("List has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("List is not associated",ERR,ERROR,*999)
    ENDIF
  
    CALL EXITS("LIST_REMOVE_DUPLICATES")
    RETURN
999 CALL ERRORS("LIST_REMOVE_DUPLICATES",ERR,ERROR)
    CALL EXITS("LIST_REMOVE_DUPLICATES")
    RETURN 1
  END SUBROUTINE LIST_REMOVE_DUPLICATES

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_SEARCH
  !###  Description:
  !###    Searches a list.
  !###  Child-subroutines: LIST_SEARCH_INTG,LIST_SEARCH_SP,LIST_SEARCH_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SEARCH_INTG_ARRAY(A,VALUE,POSITION,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SEARCH_INTG_ARRAY
    !###  Description:
    !###    Searches an integer array list A for VALUE. If the search is successful POSITION contains the index of the position
    !###    of VALUE in the list otherwise POSITION is zero.
    !###  Parent-function: LIST_SEARCH
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:)
    INTEGER(INTG), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: POSITION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL ENTERS("LIST_SEARCH_INTG_ARRAY",ERR,ERROR,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,ERR,ERROR,*999)
    
    CALL EXITS("LIST_SEARCH_INTG_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SEARCH_INTG_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SEARCH_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_INTG_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SEARCH_SP_ARRAY(A,VALUE,POSITION,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SEARCH_SP_ARRAY
    !###  Description:
    !###    Searches a single precision real array list A for VALUE. If the search is successful POSITION contains the index of
    !###    the position of VALUE in the list otherwise POSITION is zero.
    !###  Parent-function: LIST_SEARCH
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:)
    REAL(SP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: POSITION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
   
    CALL ENTERS("LIST_SEARCH_SP_ARRAY",ERR,ERROR,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,ERR,ERROR,*999)    
    
    CALL EXITS("LIST_SEARCH_SP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SEARCH_SP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SEARCH_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_SP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SEARCH_DP_ARRAY(A,VALUE,POSITION,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SEARCH_DP_ARRAY
    !###  Description:
    !###    Searches a double precision real array list A for VALUE. If the search is successful POSITION contains the index of
    !###    the position of VALUE in the list otherwise POSITION is zero.
    !###  Parent-function: LIST_SEARCH
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:)
    REAL(DP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: POSITION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL ENTERS("LIST_SEARCH_DP_ARRAY",ERR,ERROR,*999)

    !Default search method is a linear search
    CALL LIST_SEARCH_LINEAR(A,VALUE,POSITION,ERR,ERROR,*999)    
    
    CALL EXITS("LIST_SEARCH_DP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SEARCH_DP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SEARCH_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_DP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_SEARCH_LINEAR
  !###  Description:
  !###    Searches a list using the linear search method.
  !###  Child-subroutines: LIST_SEARCH_LINEAR_INTG,LIST_SEARCH_LINEAR_SP,LIST_SEARCH_LINEAR_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SEARCH_LINEAR_INTG_ARRAY(A,VALUE,POSITION,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SEARCH_LINEAR_INTG_ARRAY
    !###  Description:
    !###    Searches an integer array list A for VALUE using the linear search method. If the search is successful POSITION
    !###    contains the index of the position of VALUE in the list otherwise POSITION is zero.
    !###  Parent-function: LIST_SEARCH_LINEAR
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: A(:)
    INTEGER(INTG), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: POSITION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL ENTERS("LIST_SEARCH_LINEAR_INTG_ARRAY",ERR,ERROR,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(A(i)==VALUE) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    CALL EXITS("LIST_SEARCH_LINEAR_INTG_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SEARCH_LINEAR_INTG_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SEARCH_LINEAR_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_INTG_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SEARCH_LINEAR_SP_ARRAY(A,VALUE,POSITION,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SEARCH_LINEAR_SP_ARRAY
    !###  Description:
    !###    Searches a single precision real array list A for VALUE using the linear search method. If the search is successful
    !###    POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
    !###  Parent-function: LIST_SEARCH_LINEAR
    
    !Argument variables
    REAL(SP), INTENT(IN) :: A(:)
    REAL(SP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: POSITION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL ENTERS("LIST_SEARCH_LINEAR_SP_ARRAY",ERR,ERROR,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(ABS(A(i)-VALUE)<ZERO_TOLERANCE_SP) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    CALL EXITS("LIST_SEARCH_LINEAR_SP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SEARCH_LINEAR_SP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SEARCH_LINEAR_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_SP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SEARCH_LINEAR_DP_ARRAY(A,VALUE,POSITION,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SEARCH_LINEAR_DP_ARRAY
    !###  Description:
    !###    Searches a double precision real array list A for VALUE using the linear search method. If the search is successful
    !###    POSITION contains the index of the position of VALUE in the list otherwise POSITION is zero.
    !###  Parent-function: LIST_SEARCH_LINEAR
    
    !Argument variables
    REAL(DP), INTENT(IN) :: A(:)
    REAL(DP), INTENT(IN) :: VALUE
    INTEGER(INTG), INTENT(OUT) :: POSITION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i
    LOGICAL :: FOUND
    
    CALL ENTERS("LIST_SEARCH_LINEAR_DP_ARRAY",ERR,ERROR,*999)

    FOUND=.FALSE.
    i=1
    DO WHILE(i<=SIZE(A,1).AND..NOT.FOUND)
      IF(ABS(A(i)-VALUE)<ZERO_TOLERANCE) THEN
        FOUND=.TRUE.
      ELSE
        i=i+1
      ENDIF
    ENDDO
    IF(FOUND) THEN
      POSITION=i
    ELSE
      POSITION=0
    ENDIF
    
    CALL EXITS("LIST_SEARCH_LINEAR_DP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SEARCH_LINEAR_DP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SEARCH_LINEAR_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SEARCH_LINEAR_DP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_SORT
  !###  Description:
  !###    Sorts a list into ascending order.
  !###  Child-subroutines: LIST_SORT_INTG,LIST_SORT_SP,LIST_SORT_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_INTG_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_INTG_ARRAY
    !###  Description:
    !###    Sorts an integer array list into ascending order.
    !###  Parent-function: LIST_SORT
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL ENTERS("LIST_SORT_INTG_ARRAY",ERR,ERROR,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,ERR,ERROR,*999)    

    CALL EXITS("LIST_SORT_INTG_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_INTG_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_INTG_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_SP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_SP_ARRAY
    !###  Description:
    !###    Sorts an single precision array list into ascending order.
    !###  Parent-function: LIST_SORT
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    
    CALL ENTERS("LIST_SORT_SP_ARRAY",ERR,ERROR,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,ERR,ERROR,*999)    

    CALL EXITS("LIST_SORT_SP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_SP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_DP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_DP_ARRAY
    !###  Description:
    !###    Sorts an double precision array list into ascending order.
    !###  Parent-function: LIST_SORT
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
     
    CALL ENTERS("LIST_SORT_DP_ARRAY",ERR,ERROR,*999)

    !Default sort method is a heap sort
    CALL LIST_SORT_HEAP(A,ERR,ERROR,*999)    

    CALL EXITS("LIST_SORT_DP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_DP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_DP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_SORT_BUBBLE
  !###  Description:
  !###    Sorts a list into assending order using the bubble sort method.
  !###  Child-subroutines: LIST_SORT_BUBBLE_INTG,LIST_SORT_BUBBLE_SP,LIST_SORT_BUBBLE_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_BUBBLE_INTG_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_BUBBLE_INTG_ARRAY
    !###  Description:
    !###    BUBBLE_SORT_INTG performs a bubble sort on an integer array list
    !###  Parent-function: LIST_SORT_BUBBLE
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,VALUE
    
    CALL ENTERS("LIST_SORT_BUBBLE_INTG_ARRAY",ERR,ERROR,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    CALL EXITS("LIST_SORT_BUBBLE_INTG_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_BUBBLE_INTG_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_BUBBLE_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_INTG_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_BUBBLE_SP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_BUBBLE_SP_ARRAY
    !###  Description:
    !###    BUBBLE_SORT_SP performs a bubble sort on a single precision array list
    !###  Parent-function: LIST_SORT_BUBBLE
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(SP) :: VALUE
    
    CALL ENTERS("LIST_SORT_BUBBLE_SP_ARRAY",ERR,ERROR,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    CALL EXITS("LIST_SORT_BUBBLE_SP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_BUBBLE_SP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_BUBBLE_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_SP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_BUBBLE_DP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_BUBBLE_DP_ARRAY
    !###  Description:
    !###    BUBBLE_SORT_DP performs a bubble sort on a double precision list
    !###  Parent-function: LIST_SORT_BUBBLE
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(DP) :: VALUE
    
    CALL ENTERS("LIST_SORT_BUBBLE_DP_ARRAY",ERR,ERROR,*999)

    IF(SIZE(A,1)>1) THEN
      FLAG=SIZE(A,1)
      DO i=1,SIZE(A,1)
        k=FLAG-1
        FLAG=0
        DO j=1,k
          IF(A(j)>A(j+1)) THEN
            VALUE=A(j)
            A(j)=A(j+1)
            A(j+1)=VALUE
            FLAG=j
          ENDIF
        ENDDO
        IF(FLAG==0) EXIT
      ENDDO
    ENDIF

    CALL EXITS("LIST_SORT_BUBBLE_DP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_BUBBLE_DP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_BUBBLE_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_BUBBLE_DP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_SORT_HEAP
  !###  Description:
  !###    Sorts a list into assending order using the heap sort method.
  !###  Child-subroutines: LIST_SORT_HEAP_INTG_ARRAY,LIST_SORT_HEAP_SP_ARRAY,LIST_SORT_HEAP_DP_ARRAY

  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_HEAP_INTG_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_HEAP_INTG_ARRAY
    !###  Description:
    !###    Sorts an integer array list into assending order using the heap sort method.
    !###  Parent-function: LIST_SORT_HEAP
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L,VALUE
    
    CALL ENTERS("LIST_SORT_HEAP_INTG_ARRAY",ERR,ERROR,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    CALL EXITS("LIST_SORT_HEAP_INTG_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_HEAP_INTG_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_HEAP_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_INTG_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_HEAP_SP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_HEAP_SP_ARRAY
    !###  Description:
    !###    Sorts a real single precision array list into assending order using the heap sort method.
    !###  Parent-function: LIST_SORT_HEAP
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(SP) :: VALUE
    
    CALL ENTERS("LIST_SORT_HEAP_SP_ARRAY",ERR,ERROR,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    CALL EXITS("LIST_SORT_HEAP_SP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_HEAP_SP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_HEAP_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_SP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_HEAP_DP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_HEAP_DP_ARRAY
    !###  Description:
    !###    Sorts a real double precision array list into assending order using the heap sort method.
    !###  Parent-function: HEAP_SORT
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(DP) :: VALUE
    
    CALL ENTERS("LIST_SORT_HEAP_DP_ARRAY",ERR,ERROR,*999)

    IF(SIZE(A,1)>1) THEN      
      L=SIZE(A,1)/2+1
      IVALUE=SIZE(A,1)
      DO 
        IF(L>1) THEN
          L=L-1
          VALUE=A(L)
        ELSE
          VALUE=A(IVALUE)
          A(IVALUE)=A(1)
          IVALUE=IVALUE-1
          IF(IVALUE==1) THEN
            A(1)=VALUE
            EXIT
          ENDIF
        ENDIF
        I=L
        J=L+L
        DO WHILE(J<=IVALUE)
          IF(J<IVALUE) THEN
            IF(A(J)<A(J+1)) J=J+1
          ENDIF
          IF(VALUE<A(J)) THEN
            A(I)=A(J)
            I=J
            J=J+J
          ELSE
            J=IVALUE+1
          ENDIF
        ENDDO
        A(I)=VALUE
      ENDDO
    ENDIF

    CALL EXITS("LIST_SORT_HEAP_DP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_HEAP_DP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_HEAP_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_HEAP_DP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: LIST_SORT_SHELL
  !###  Description:
  !###    Sorts a list into either assending or descending order using the shell sort method.
  !###  Child-subroutines: LIST_SORT_SHELL_INTG_ARRAY,LIST_SORT_SHELL_SP_ARRAY,LIST_SORT_SHELL_DP_ARRAY

  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_SHELL_INTG_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_SHELL_INTG_ARRAY
    !###  Description:
    !###    Sorts an integer array list into either assending or descending order using the shell sort method.
    !###  Parent-function: LIST_SORT_SHELL
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J,VALUE
    
    CALL ENTERS("LIST_SORT_SHELL_INTG_ARRAY",ERR,ERROR,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    CALL EXITS("LIST_SORT_SHELL_INTG_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_SHELL_INTG_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_SHELL_INTG_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_INTG_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_SHELL_SP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_SHELL_SP_ARRAY
    !###  Description:
    !###    Sorts a real single precision array list into either assending or descending order using the shell sort method.
    !###  Parent-function: LIST_SORT_SHELL
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(SP) :: VALUE
    
    CALL ENTERS("LIST_SORT_SHELL_SP_ARRAY",ERR,ERROR,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    CALL EXITS("LIST_SORT_SHELL_SP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_SHELL_SP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_SHELL_SP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_SP_ARRAY
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE LIST_SORT_SHELL_DP_ARRAY(A,ERR,ERROR,*)
  
    !#### Subroutine: LIST_SORT_SHELL_DP_ARRAY
    !###  Description:
    !###    Sorts a real double precision array list into either assending or descending order using the shell sort method.
    !###  Parent-function: LIST_SORT_SHELL
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(DP) :: VALUE
    
    CALL ENTERS("LIST_SORT_SHELL_DP_ARRAY",ERR,ERROR,*999)

    INC=4
    DO WHILE(INC<=SIZE(A,1))
      INC=3*INC+1
    ENDDO
    DO WHILE(INC>1)
      INC=INC/3
      DO i=INC+1,SIZE(A,1)
        VALUE=A(i)
        J=I
        DO WHILE(A(J-INC)>VALUE)
          A(J)=A(J-INC)
          J=J-INC
          IF(J<=INC) EXIT
        ENDDO
        A(J)=VALUE
      ENDDO !i
    ENDDO

    CALL EXITS("LIST_SORT_SHELL_DP_ARRAY")
    RETURN
999 CALL ERRORS("LIST_SORT_SHELL_DP_ARRAY",ERR,ERROR)
    CALL EXITS("LIST_SORT_SHELL_DP_ARRAY")
    RETURN 1
  END SUBROUTINE LIST_SORT_SHELL_DP_ARRAY
  
  !
  !================================================================================================================================
  !

END MODULE LISTS
