!> Only for integer data type for now
MODULE LinkedList_routines

#include "macros.h"  

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  implicit none
  
  private  ! by default

  ! types
  type LinkedListItem
    integer(intg) :: data
    type(LinkedListItem),pointer :: next => NULL()
  end type

  type LinkedList
    type(LinkedListItem),pointer :: root => NULL()
    type(LinkedListItem),pointer :: last => NULL()
  end type

  interface LinkedList_Add
    module procedure LinkedList_Add_Data
    module procedure LinkedList_Add_List
  end interface

  ! public types
  public :: LinkedListItem,LinkedList

  ! public subs
  public :: LinkedList_Add,LinkedList_Destroy,LinkedList_Remove_First,LinkedList_Remove_Last
  public :: LinkedList_is_Empty,LinkedList_to_Array

contains

! -------------------------------------------------------------------

  !> initialises or adds a piece of data to list
  SUBROUTINE LinkedList_Add_Data(list,data,ERR,ERROR,*)

    TYPE(LinkedList),INTENT(INOUT) :: list
    INTEGER(INTG),INTENT(IN) :: data
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(LinkedListItem),pointer :: current

    ENTERS("LinkedList_Add_Data",ERR,ERROR,*999)

    if (associated(list%root)) then
      ! add to the tail end (for now)
      current => list%last
      allocate(current%next)
      current%next%data = data
      list%last => current%next
    else
      allocate(list%root)
      list%root%data = data
      list%last => list%root
    endif

    EXITS("LinkedList_Add_Data")
    RETURN
999 ERRORSEXITS("LinkedList_Add_Data",ERR,ERROR)
    RETURN 1

  End SUBROUTINE LinkedList_Add_Data


! -------------------------------------------------------------------

  !> adds all data from one list to another
  Subroutine LinkedList_Add_List(list,addlist,ERR,ERROR,*)
    type(LinkedList),intent(inout) :: list
    type(LinkedList),intent(in) :: addlist
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    ! local variables
    type(LinkedListItem),pointer :: current

    if (LinkedList_is_Empty(addlist)) return

    current => addlist%root
    do
      call LinkedList_Add_Data(list,current%data,ERR,ERROR,*999)
      if (associated(current%next)) then
        current => current%next
      else
        exit
      endif
    enddo

    EXITS("LinkedList_Add_List")
    RETURN
999 ERRORSEXITS("LinkedList_Add_List",ERR,ERROR)
    RETURN 1

  End Subroutine LinkedList_Add_List

! -------------------------------------------------------------------

  !> removes the first item from list and returns its value in data
  Subroutine LinkedList_Remove_First(list,data,ERR,ERROR,*)
    type(LinkedList),intent(inout) :: list
    integer(intg),intent(out) :: data
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    ! local variables
    type(LinkedListItem),pointer :: next
    
    if (associated(list%root)) then
      data = list%root%data
      next => list%root%next
      deallocate(list%root)
      list%root => next
      if (associated(list%root)) then
        if (.not.associated(list%root%next)) list%last => list%root  ! only one left
      else
        list%last => NULL()
      endif
    else
      write(*,*) ">>> warning: linked list is empty and cannot remove first item"
    endif

  End Subroutine LinkedList_Remove_First

! -------------------------------------------------------------------

  !> removes the first item from list and returns its value in data
  Subroutine LinkedList_Remove_Last(list,data,ERR,ERROR,*)
    type(LinkedList),intent(inout) :: list
    integer(intg),intent(out) :: data
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    ! local variables
    type(LinkedListItem),pointer :: current

    if (.not.associated(list%root)) then
      write(*,*) ">>> warning: linked list is empty and cannot remove last item"
      return
    endif
    current => list%root

    do
      if (associated(current%next)) then
        if (associated(current%next%next)) then
          current => current%next
        else
          ! next one is the last one
          data = current%next%data
          deallocate(current%next)
          current%next => NULL()
          list%last => current
          exit
        endif
      else
        ! there must be only one item in the list?
        data = current%data
        deallocate(list%root)
        list%root => NULL()
        list%last => NULL()
        exit
      endif
    enddo

  End Subroutine LinkedList_Remove_Last

! -------------------------------------------------------------------

  !> will delete and deallocate all items
  Subroutine LinkedList_Destroy(list,ERR,ERROR,*)

    TYPE(LinkedList), INTENT(inout) :: list
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    ! local variables
    TYPE(LinkedListItem), POINTER :: current,next

    if (.not.associated(list%root)) return

    current => list%root
    do
      if (associated(current%next)) then
        next => current%next
        deallocate(current)
        current => next
      else
        deallocate(current)
        exit
      endif
    enddo
    
    list%root => NULL()
    list%last => NULL()

  End SUBROUTINE LinkedList_Destroy

! -------------------------------------------------------------------

  !> returns true if the list is empty
  Function LinkedList_is_Empty(list)
    type(LinkedList),intent(in) :: list
    logical :: LinkedList_is_Empty

    LinkedList_is_Empty = .true.
    if (associated(list%root)) LinkedList_is_Empty = .false.
    
  End Function LinkedList_is_Empty

! -------------------------------------------------------------------

  !> copies out the data to an allocatable array
  Subroutine LinkedList_to_Array(list,array,ERR,ERROR,*)
    type(LinkedList),intent(in) :: list
    integer(INTG),allocatable,intent(out) :: array(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    ! local variables
    integer(INTG) :: i,n
    type(LinkedListItem),pointer :: current

    ! return zero-size array if list is empty
    if (LinkedList_is_Empty(list)) then
      allocate(array(0))
      return
    endif

    ! first traversing to find size
    current => list%root
    n=1
    do
      if (associated(current%next)) then
        n=n+1
        current => current%next
      else
        exit
      endif
    enddo

    ! copy to array
    if (allocated(array)) deallocate(array)
    allocate(array(n),stat=err)
    !IF (ERR/=0) CALL ...
    current => list%root
    do i=1,n
      array(i)=current%data
      current => current%next
    enddo

  End Subroutine LinkedList_to_Array

End Module LinkedList_routines
