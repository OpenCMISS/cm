!> \file
!> \author Chris Bradley
!> \brief This module contains all procedures for sorting. NOTE: THE ROUTINES IN THIS MODULE HAVE NOT BEEN TESTED!!!
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

!> This module contains all procedures for sorting. NOTE: THE ROUTINES IN THIS MODULE HAVE NOT BEEN TESTED!!!
MODULE SORTING

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  INTERFACE BUBBLE_ISORT
    MODULE PROCEDURE BUBBLE_ISORT_INTG  
    MODULE PROCEDURE BUBBLE_ISORT_SP
    MODULE PROCEDURE BUBBLE_ISORT_DP
  END INTERFACE !BUBBLE_ISORT

  INTERFACE BUBBLE_SORT
    MODULE PROCEDURE BUBBLE_SORT_INTG
    MODULE PROCEDURE BUBBLE_SORT_SP
    MODULE PROCEDURE BUBBLE_SORT_DP
  END INTERFACE !BUBBLE_SORT

  INTERFACE HEAP_SORT
    MODULE PROCEDURE HEAP_SORT_INTG
    MODULE PROCEDURE HEAP_SORT_SP
    MODULE PROCEDURE HEAP_SORT_DP
  END INTERFACE !HEAP_SORT

  INTERFACE SHELL_SORT
    MODULE PROCEDURE SHELL_SORT_INTG
    MODULE PROCEDURE SHELL_SORT_SP
    MODULE PROCEDURE SHELL_SORT_DP
  END INTERFACE !SHELL_SORT

  PUBLIC BUBBLE_ISORT
  PUBLIC BUBBLE_SORT,HEAP_SORT,SHELL_SORT

CONTAINS

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: BUBBLE_ISORT
  !###  Description:
  !###    Sorts a list into assending order using the bubble sort method, returning sorting index
  !###  Child-subroutines: BUBBLE_ISORT_INTG,BUBBLE_ISORT_SP,BUBBLE_ISORT_DP
  
  !
  !================================================================================================================================
  !

  SUBROUTINE BUBBLE_ISORT_INTG(A,IND,ERR,ERROR,*)
  
    !#### Subroutine: BUBBLE_ISORT_INTG
    !###  Description:
    !###    BUBBLE_ISORT_INTG performs a bubble sort on an integer list, returning sorting index
    !###  Parent-function: BUBBLE_ISORT
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: IND(:)      
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,VALUE,IVALUE
    
#if DEBUG
    CALL ENTERS("BUBBLE_ISORT_INTG",ERR,ERROR,*999)
#endif

    IF(SIZE(IND,1)==SIZE(A,1)) THEN
      IND(1)=1  
      IF(SIZE(A,1)>1) THEN
        FLAG=SIZE(A,1)
        DO i=1,SIZE(A,1)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(i==1) IND(j+1)=j+1          
            IF(A(j)>A(j+1)) THEN
              VALUE=A(j)
              A(j)=A(j+1)
              A(j+1)=VALUE
              IVALUE=IND(j)
              IND(j)=IND(j+1)
              IND(j+1)=IVALUE              
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Size of input vectors does not match",ERR,ERROR,*999)
    ENDIF      

#if DEBUG
    CALL EXITS("BUBBLE_ISORT_INTG")
#endif
    RETURN
999 CALL ERRORS("BUBBLE_ISORT_INTG",ERR,ERROR)
#if DEBUG
    CALL EXITS("BUBBLE_ISORT_INTG")
#endif
    RETURN 1
  END SUBROUTINE BUBBLE_ISORT_INTG
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE BUBBLE_ISORT_SP(A,IND,ERR,ERROR,*)
  
    !#### Subroutine: BUBBLE_ISORT_SP
    !###  Description:
    !###    BUBBLE_ISORT_SP performs a bubble sort on a single precision list, returning sorting index
    !###  Parent-function: BUBBLE_ISORT
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: IND(:)      
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,IVALUE
    REAL(SP) :: VALUE
    
#if DEBUG
    CALL ENTERS("BUBBLE_ISORT_SP",ERR,ERROR,*999)
#endif

    IF(SIZE(IND,1)==SIZE(A,1)) THEN
      IND(1)=1  
      IF(SIZE(A,1)>1) THEN
        FLAG=SIZE(A,1)
        DO i=1,SIZE(A,1)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(i==1) IND(j+1)=j+1             
            IF(A(j)>A(j+1)) THEN
              VALUE=A(j)
              A(j)=A(j+1)
              A(j+1)=VALUE
              IVALUE=IND(j)
              IND(j)=IND(j+1)
              IND(j+1)=IVALUE              
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Size of input vectors does not match",ERR,ERROR,*999)
    ENDIF      

#if DEBUG
    CALL EXITS("BUBBLE_ISORT_SP")
#endif
    RETURN
999 CALL ERRORS("BUBBLE_ISORT_SP",ERR,ERROR)
#if DEBUG
    CALL EXITS("BUBBLE_ISORT_SP")
#endif
    RETURN 1
  END SUBROUTINE BUBBLE_ISORT_SP
  
  !
  !================================================================================================================================
  !

  SUBROUTINE BUBBLE_ISORT_DP(A,IND,ERR,ERROR,*)
  
    !#### Subroutine: BUBBLE_ISORT_DP
    !###  Description:
    !###    BUBBLE_ISORT_DP performs a bubble sort on a double precision list, returning sorting index
    !###  Parent-function: BUBBLE_ISORT
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: IND(:)    
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,IVALUE
    REAL(DP) :: VALUE
    
#if DEBUG
    CALL ENTERS("BUBBLE_ISORT_DP",ERR,ERROR,*999)
#endif
    
    IF(SIZE(IND,1)==SIZE(A,1)) THEN
      IND(1)=1
      IF(SIZE(A,1)>1) THEN
        FLAG=SIZE(A,1)
        DO i=1,SIZE(A,1)
          k=FLAG-1
          FLAG=0
          DO j=1,k
            IF(i==1) IND(j+1)=j+1          
            IF(A(j)>A(j+1)) THEN
              VALUE=A(j)
              A(j)=A(j+1)
              A(j+1)=VALUE
              IVALUE=IND(j)
              IND(j)=IND(j+1)
              IND(j+1)=IVALUE
              FLAG=j
            ENDIF
          ENDDO
          IF(FLAG==0) EXIT
        ENDDO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Size of input vectors does not match",ERR,ERROR,*999)
    ENDIF

#if DEBUG
    CALL EXITS("BUBBLE_ISORT_DP")
#endif
    RETURN
999 CALL ERRORS("BUBBLE_ISORT_DP",ERR,ERROR)
#if DEBUG
    CALL EXITS("BUBBLE_ISORT_DP")
#endif
    RETURN 1
  END SUBROUTINE BUBBLE_ISORT_DP

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: BUBBLE_SORT
  !###  Description:
  !###    Sorts a list into assending order using the bubble sort method.
  !###  Child-subroutines: BUBBLE_SORT_INTG,BUBBLE_SORT_SP,BUBBLE_SORT_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE BUBBLE_SORT_INTG(A,ERR,ERROR,*)
  
    !#### Subroutine: BUBBLE_SORT_INTG
    !###  Description:
    !###    BUBBLE_SORT_INTG performs a bubble sort on an integer list
    !###  Parent-function: BUBBLE_SORT
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k,VALUE
    
#if DEBUG
    CALL ENTERS("BUBBLE_SORT_INTG",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("BUBBLE_SORT_INTG")
#endif
    RETURN
999 CALL ERRORS("BUBBLE_SORT_INTG",ERR,ERROR)
#if DEBUG
    CALL EXITS("BUBBLE_SORT_INTG")
#endif
    RETURN 1
  END SUBROUTINE BUBBLE_SORT_INTG
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE BUBBLE_SORT_SP(A,ERR,ERROR,*)
  
    !#### Subroutine: BUBBLE_SORT_SP
    !###  Description:
    !###    BUBBLE_SORT_SP performs a bubble sort on a single precision list
    !###  Parent-function: BUBBLE_SORT
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(SP) :: VALUE
    
#if DEBUG
    CALL ENTERS("BUBBLE_SORT_SP",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("BUBBLE_SORT_SP")
#endif
    RETURN
999 CALL ERRORS("BUBBLE_SORT_SP",ERR,ERROR)
#if DEBUG
    CALL EXITS("BUBBLE_SORT_SP")
#endif
    RETURN 1
  END SUBROUTINE BUBBLE_SORT_SP
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE BUBBLE_SORT_DP(A,ERR,ERROR,*)
  
    !#### Subroutine: BUBBLE_SORT_DP
    !###  Description:
    !###    BUBBLE_SORT_DP performs a bubble sort on a double precision list
    !###  Parent-function: BUBBLE_SORT
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: FLAG,i,j,k
    REAL(DP) :: VALUE
    
#if DEBUG
    CALL ENTERS("BUBBLE_SORT_DP",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("BUBBLE_SORT_DP")
#endif
    RETURN
999 CALL ERRORS("BUBBLE_SORT_DP",ERR,ERROR)
#if DEBUG
    CALL EXITS("BUBBLE_SORT_DP")
#endif
    RETURN 1
  END SUBROUTINE BUBBLE_SORT_DP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: HEAP_SORT
  !###  Description:
  !###    Sorts a list into assending order using the heap sort method.
  !###  Child-subroutines: HEAP_SORT_INTG,HEAP_SORT_SP,HEAP_SORT_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE HEAP_SORT_INTG(A,ERR,ERROR,*)
  
    !#### Subroutine: HEAP_SORT_INTG
    !###  Description:
    !###    HEAP_SORT_INTG performs a heap sort on an integer list
    !###  Parent-function: HEAP_SORT
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L,VALUE
    
#if DEBUG
    CALL ENTERS("HEAP_SORT_INTG",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("HEAP_SORT_INTG")
#endif
    RETURN
999 CALL ERRORS("HEAP_SORT_INTG",ERR,ERROR)
#if DEBUG
    CALL EXITS("HEAP_SORT_INTG")
#endif
    RETURN 1
  END SUBROUTINE HEAP_SORT_INTG
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE HEAP_SORT_SP(A,ERR,ERROR,*)
  
    !#### Subroutine: HEAP_SORT_SP
    !###  Description:
    !###    HEAP_SORT_SP performs a heap sort on a single precision list
    !###  Parent-function: HEAP_SORT
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(SP) :: VALUE
    
#if DEBUG
    CALL ENTERS("HEAP_SORT_SP",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("HEAP_SORT_SP")
#endif
    RETURN
999 CALL ERRORS("HEAP_SORT_SP",ERR,ERROR)
#if DEBUG
    CALL EXITS("HEAP_SORT_SP")
#endif
    RETURN 1
  END SUBROUTINE HEAP_SORT_SP
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE HEAP_SORT_DP(A,ERR,ERROR,*)
  
    !#### Subroutine: HEAP_SORT_DP
    !###  Description:
    !###    HEAP_SORT_DP performs a heap sort on a double precision list
    !###  Parent-function: HEAP_SORT
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: I,IVALUE,J,L
    REAL(DP) :: VALUE
    
#if DEBUG
    CALL ENTERS("HEAP_SORT_DP",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("HEAP_SORT_DP")
#endif
    RETURN
999 CALL ERRORS("HEAP_SORT_DP",ERR,ERROR)
#if DEBUG
    CALL EXITS("HEAP_SORT_DP")
#endif
    RETURN 1
  END SUBROUTINE HEAP_SORT_DP
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: SHELL_SORT
  !###  Description:
  !###    Sorts a list into either assending or descending order using the shell sort method.
  !###  Child-subroutines: SHELL_SORT_INTG,SHELL_SORT_SP,SHELL_SORT_DP

  !
  !================================================================================================================================
  !
  
  SUBROUTINE SHELL_SORT_INTG(A,ERR,ERROR,*)
  
    !#### Subroutine: SHELL_SORT_INTG
    !###  Description:
    !###    SHELL_SORT_INTG performs a shell sort on an integer list
    !###  Parent-function: SHELL_SORT
    
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: I,INC,J,VALUE
    
#if DEBUG
    CALL ENTERS("SHELL_SORT_INTG",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("SHELL_SORT_INTG")
#endif
    RETURN
999 CALL ERRORS("SHELL_SORT_INTG",ERR,ERROR)
#if DEBUG
    CALL EXITS("SHELL_SORT_INTG")
#endif
    RETURN 1
  END SUBROUTINE SHELL_SORT_INTG
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE SHELL_SORT_SP(A,ERR,ERROR,*)
  
    !#### Subroutine: SHELL_SORT_SP
    !###  Description:
    !###    SHELL_SORT_SP performs a shell sort on a single precision list
    !###  Parent-function: SHELL_SORT
    
    !Argument variables
    REAL(SP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(SP) :: VALUE
    
#if DEBUG
    CALL ENTERS("SHELL_SORT_SP",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("SHELL_SORT_SP")
#endif
    RETURN
999 CALL ERRORS("SHELL_SORT_SP",ERR,ERROR)
#if DEBUG
    CALL EXITS("SHELL_SORT_SP")
#endif
    RETURN 1
  END SUBROUTINE SHELL_SORT_SP
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE SHELL_SORT_DP(A,ERR,ERROR,*)
  
    !#### Subroutine: SHELL_SORT_DP
    !###  Description:
    !###    SHELL_SORT_DP performs a shell sort on a double precision list
    !###  Parent-function: SHELL_SORT
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: A(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: I,INC,J
    REAL(DP) :: VALUE
    
#if DEBUG
    CALL ENTERS("SHELL_SORT_DP",ERR,ERROR,*999)
#endif

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

#if DEBUG
    CALL EXITS("SHELL_SORT_DP")
#endif
    RETURN
999 CALL ERRORS("SHELL_SORT_DP",ERR,ERROR)
#if DEBUG
    CALL EXITS("SHELL_SORT_DP")
#endif
    RETURN 1
  END SUBROUTINE SHELL_SORT_DP
  
  !
  !================================================================================================================================
  !
  
END MODULE SORTING
