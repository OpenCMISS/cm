!> \file
!> $Id: f90c_f.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module handles calling C from Fortran90 and vise versa. It needs to be linked with the f90c_c.c module.
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

!> This module handles calling C from Fortran90 and vise versa. It needs to be linked with the f90c_c.c module.
MODULE F90C

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE MACHINE_CONSTANTS
  USE ISO_VARYING_STRING

  IMPLICIT NONE

  INTERFACE 
    
    SUBROUTINE CSTRINGLEN(LENGTH, CSTRING)
      !DEC$ ATTRIBUTES C, REFERENCE :: cstringlen
      USE CONSTANTS
      INTEGER(INTG), INTENT(OUT) :: LENGTH
      INTEGER(INTG), INTENT(IN) :: CSTRING(*)
    END SUBROUTINE CSTRINGLEN

    SUBROUTINE PACKCHARACTERS(INTCHAR, POSITION, CSTRING)
      !DEC$ ATTRIBUTES C, REFERENCE :: packcharacters
      USE CONSTANTS
      INTEGER(INTG), INTENT(IN) :: INTCHAR, POSITION
      INTEGER(INTG), INTENT(OUT) :: CSTRING(*)
    END SUBROUTINE PACKCHARACTERS

    SUBROUTINE UNPACKCHARACTERS(INTCHAR, POSITION, CSTRING)
      !DEC$ ATTRIBUTES C, REFERENCE :: unpackcharacters
      USE CONSTANTS
      INTEGER(INTG), INTENT(OUT) :: INTCHAR
      INTEGER(INTG), INTENT(IN) :: POSITION, CSTRING(*)
    END SUBROUTINE UNPACKCHARACTERS
    
  END INTERFACE
  
CONTAINS

  !
  !============================================================================
  !

  FUNCTION CSTRINGLENGTH(CSTRING)
    
    !#### Function: CSTRINGLENGTH
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###     CSTRINGLENGTH returns the length of a c integer string.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: CSTRING(*)
    !Function variable
    INTEGER(INTG) :: CSTRINGLENGTH
    !Local Variables
    INTEGER(INTG) :: LENGTH
    
    CALL CSTRINGLEN(LENGTH, CSTRING)

    CSTRINGLENGTH=LENGTH

    RETURN 
  END FUNCTION CSTRINGLENGTH

  !
  !============================================================================
  !

  SUBROUTINE C2FSTRING(CSTRING,FSTRING,ERR,ERROR,*)

    !#### Subroutine: C2FSTRING
    !###  Description:
    !###    C2FSTRING converts a c integer string to a Fortran
    !###    string.
    
    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: CSTRING(*)
    CHARACTER(LEN=*), INTENT(OUT) :: FSTRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: I, INTCHAR, LENGTH 
    
    CALL ENTERS("C2FSTRING",ERR,ERROR,*999)

    LENGTH=CSTRINGLENGTH(CSTRING)
    IF(LENGTH<=LEN(FSTRING)) THEN
      FSTRING=" "
      DO I=1,LENGTH
        CALL UNPACKCHARACTERS(INTCHAR,I-1,CSTRING)
        FSTRING(I:I)=CHAR(INTCHAR)
      ENDDO!I
    ELSE
      CALL FLAG_ERROR("Fortran string not big enough to hold c string",&
        & ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("C2FSTRING")
    RETURN
999 CALL ERRORS("C2FSTRING",ERR,ERROR)
    CALL EXITS("C2FSTRING")
    RETURN 1
  END SUBROUTINE C2FSTRING

  !
  !============================================================================
  !

  FUNCTION FSTRINGLENGTH(FSTRING)
    
    !#### Function: FSTRINGLENGTH
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    FSTRINGLENGTH returns the length of a fortran string
    !###    with the trailing blanks trimed.

    !Argument Variables
    CHARACTER(LEN=*), INTENT(IN) :: FSTRING
    !Function variable
    INTEGER(INTG) :: FSTRINGLENGTH
    !Local Variables
    
    FSTRINGLENGTH=LEN_TRIM(FSTRING)
    
    RETURN 
  END FUNCTION FSTRINGLENGTH

  !
  !============================================================================
  !

  SUBROUTINE F2CSTRING(CSTRING,FSTRING,ERR,ERROR,*,FSTRINGLEN)

    !#### Subroutine: F2CSTRING
    !###  Description:
    !###    F2CSTRING converts a Fortran90 string into an integer c
    !###    string. If present, the optional argument FSTRINGLEN will
    !###    be used as the string length otherwise the trim length
    !###    will be used.

    !Argument variables    
    INTEGER(INTG), INTENT(OUT) :: CSTRING(:)
    CHARACTER(LEN=*), INTENT(IN) :: FSTRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    INTEGER(INTG), OPTIONAL, INTENT(IN) :: FSTRINGLEN
    !Local Variables
    INTEGER(INTG) :: I,INTCHAR,LENGTH
    
    CALL ENTERS("F2CSTRING",ERR,ERROR,*999)

    IF(PRESENT(FSTRINGLEN)) THEN
      LENGTH=FSTRINGLEN
    ELSE
      LENGTH=FSTRINGLENGTH(FSTRING)
    ENDIF
    IF(LENGTH<=SIZE(CSTRING,1)*INTEGER_SIZE) THEN
      DO I=1,LENGTH
        INTCHAR=ICHAR(FSTRING(I:I))
        CALL PACKCHARACTERS(INTCHAR,I-1,CSTRING)
      ENDDO!I
      INTCHAR=0
      CALL PACKCHARACTERS(INTCHAR,LENGTH,CSTRING)
    ELSE
      CALL FLAG_ERROR("C string is not big enough to hold fortran string",&
        & ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("F2CSTRING")
    RETURN
999 CALL ERRORS("F2CSTRING",ERR,ERROR)
    CALL EXITS("F2CSTRING")
    RETURN 1
  END SUBROUTINE F2CSTRING

  !
  !============================================================================
  !

END MODULE F90C
