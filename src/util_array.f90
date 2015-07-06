!> \file
!> \author Caton Little
!> \brief This module contains various routines for manipulating arrays
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

!> Implements various dynamic array routines
MODULE UTIL_ARRAY
  USE BASE_ROUTINES
  USE TYPES

  IMPLICIT NONE

  INTERFACE REALLOCATE
    MODULE PROCEDURE REALLOCATE_INT
    MODULE PROCEDURE REALLOCATE_REAL
    MODULE PROCEDURE REALLOCATE_STRING
  END INTERFACE REALLOCATE

  INTERFACE GROW_ARRAY
    MODULE PROCEDURE GROW_ARRAY_INT
    MODULE PROCEDURE GROW_ARRAY_REAL
  END INTERFACE GROW_ARRAY
  
  PUBLIC :: REALLOCATE, GROW_ARRAY

CONTAINS

  !
  !================================================================================================================================
  !
  
  SUBROUTINE REALLOCATE_INT( array, newSize, errorMessage, ERR, ERROR, * )
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER(INTG), INTENT(IN) :: newSize
    CHARACTER(LEN=*), INTENT(IN) :: errorMessage
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
#if DEBUG
    CALL ENTERS("REALLOCATE_INT",ERR,ERROR,*999)
#endif

    IF( ALLOCATED( array ) ) THEN
      DEALLOCATE( array )
    ENDIF
    
    ALLOCATE( array( newSize ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( errorMessage, ERR, ERROR, *999)
    
    array(:) = 0

#if DEBUG
    CALL EXITS("REALLOCATE_INT")
#endif
    RETURN
999 CALL ERRORS("REALLOCATE_INT",ERR,ERROR)
#if DEBUG
    CALL EXITS("REALLOCATE_INT")
#endif
  END SUBROUTINE REALLOCATE_INT

  !
  !================================================================================================================================
  !
  
  SUBROUTINE REALLOCATE_REAL( array, newSize, errorMessage, ERR, ERROR, * )
    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER(INTG), INTENT(IN) :: newSize
    CHARACTER(LEN=*), INTENT(IN) :: errorMessage
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
#if DEBUG
    CALL ENTERS("REALLOCATE_REAL",ERR,ERROR,*999)
#endif

    IF( ALLOCATED( array ) ) THEN
      DEALLOCATE( array )
    ENDIF
    
    ALLOCATE( array( newSize ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( errorMessage, ERR, ERROR, *999)
    
    array(:) = 0

#if DEBUG
    CALL EXITS("REALLOCATE_REAL")
#endif
    RETURN
999 CALL ERRORS("REALLOCATE_REAL",ERR,ERROR)
#if DEBUG
    CALL EXITS("REALLOCATE_REAL")
#endif
  END SUBROUTINE REALLOCATE_REAL

  !
  !================================================================================================================================
  !
  
  SUBROUTINE REALLOCATE_STRING( array, newSize, errorMessage, ERR, ERROR, * )
    TYPE(VARYING_STRING), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER(INTG), INTENT(IN) :: newSize
    CHARACTER(LEN=*), INTENT(IN) :: errorMessage
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: I
    
#if DEBUG
    CALL ENTERS("REALLOCATE_STRING",ERR,ERROR,*999)
#endif

    IF( ALLOCATED( array ) ) THEN
      DO I=1,SIZE(ARRAY,1)
        CALL ERASE(ARRAY(I))
        DEALLOCATE( array )
      ENDDO
    ENDIF
    
    ALLOCATE( array( newSize ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( errorMessage, ERR, ERROR, *999)
    
#if DEBUG
    CALL EXITS("REALLOCATE_STRING")
#endif
    RETURN
999 CALL ERRORS("REALLOCATE_STRING",ERR,ERROR)
#if DEBUG
    CALL EXITS("REALLOCATE_STRING")
#endif
  END SUBROUTINE REALLOCATE_STRING
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE REALLOCATE_2D( array, newSize1, newSize2, errorMessage, ERR, ERROR, * )
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
    INTEGER(INTG), INTENT(IN) :: newSize1
    INTEGER(INTG), INTENT(IN) :: newSize2
    CHARACTER(LEN=*), INTENT(IN) :: errorMessage
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
#if DEBUG
    CALL ENTERS("REALLOCATE_2D",ERR,ERROR,*999)
#endif

    IF( ALLOCATED( array ) ) THEN
      DEALLOCATE( array )
    ENDIF
    
    ALLOCATE( array( newSize1, newSize2 ), STAT = ERR )
    IF( ERR /= 0 ) CALL FLAG_ERROR( errorMessage, ERR, ERROR, *999)
    
    array(:,:) = 0

#if DEBUG
    CALL EXITS("REALLOCATE_2D")
#endif
    RETURN
999 CALL ERRORS("REALLOCATE_2D",ERR,ERROR)
#if DEBUG
    CALL EXITS("REALLOCATE_2D")
#endif
  END SUBROUTINE REALLOCATE_2D

  !
  !================================================================================================================================
  !

  SUBROUTINE GROW_ARRAY_INT( array, delta, errorMessage, ERR, ERROR, * )
    INTEGER(INTG), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER(INTG), INTENT(IN) :: delta
    CHARACTER(LEN=*), INTENT(IN) :: errorMessage
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    INTEGER(INTG), ALLOCATABLE :: tempArray(:)
    INTEGER(INTG) :: oldSize
    
#if DEBUG
    CALL ENTERS("GROW_ARRAY_INT",ERR,ERROR,*999)
#endif

    IF( .NOT.ALLOCATED( array ) ) THEN
      CALL REALLOCATE( array, delta, errorMessage, ERR, ERROR, *999 )
      RETURN
    ENDIF
    
    oldSize = SIZE( array )
    
    CALL REALLOCATE( tempArray, oldSize, errorMessage, ERR, ERROR, *999 )
    
    tempArray(:) = array(:)
    
    CALL REALLOCATE( array, oldSize + delta, errorMessage, ERR, ERROR, *999 )
    
    array(1:oldSize) = tempArray(:)

    DEALLOCATE( tempArray )

#if DEBUG
    CALL EXITS("GROW_ARRAY_INT")
#endif
    RETURN
999 CALL ERRORS("GROW_ARRAY_INT",ERR,ERROR)
#if DEBUG
    CALL EXITS("GROW_ARRAY_INT")
#endif
  END SUBROUTINE GROW_ARRAY_INT
  
  !
  !================================================================================================================================
  !

  SUBROUTINE GROW_ARRAY_REAL( array, delta, errorMessage, ERR, ERROR, * )
    REAL(C_DOUBLE), ALLOCATABLE, INTENT(INOUT) :: array(:)
    INTEGER(INTG), INTENT(IN) :: delta
    CHARACTER(LEN=*), INTENT(IN) :: errorMessage
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    REAL(C_DOUBLE), ALLOCATABLE :: tempArray(:)
    INTEGER(INTG) :: oldSize
    
#if DEBUG
    CALL ENTERS("GROW_ARRAY_REAL",ERR,ERROR,*999)
#endif

    IF( .NOT.ALLOCATED( array ) ) THEN
      CALL REALLOCATE( array, delta, errorMessage, ERR, ERROR, *999 )
      RETURN
    ENDIF
    
    oldSize = SIZE( array )
    
    CALL REALLOCATE( tempArray, oldSize, errorMessage, ERR, ERROR, *999 )
    
    tempArray(:) = array(:)
    
    CALL REALLOCATE( array, oldSize + delta, errorMessage, ERR, ERROR, *999 )
    
    array(1:oldSize) = tempArray(:)

    DEALLOCATE( tempArray )

#if DEBUG
    CALL EXITS("GROW_ARRAY_REAL")
#endif
    RETURN
999 CALL ERRORS("GROW_ARRAY_REAL",ERR,ERROR)
#if DEBUG
    CALL EXITS("GROW_ARRAY_REAL")
#endif
  END SUBROUTINE GROW_ARRAY_REAL
  
  !
  !================================================================================================================================
  !

END MODULE UTIL_ARRAY
