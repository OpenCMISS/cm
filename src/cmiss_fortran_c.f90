!> \file
!> \author Chris Bradley
!> \brief This module contains Fortran and C data conversion routines.
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

!> This module contains all Fortran and C data conversion routines.
MODULE CMISS_FORTRAN_C

  USE ISO_C_BINDING
 
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISSC2FString

  PUBLIC CMISSF2CString
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Copys/converts a C string (array of characters) to a Fortran String (length of characters)
  SUBROUTINE CMISSC2FString(Cstring,Fstring)
    !Argument variables
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Cstring(:)
    CHARACTER(LEN=*), INTENT(OUT) :: Fstring
    !Local variables
    INTEGER(C_INT) :: i,LENGTH

    IF(LEN(Fstring)>=SIZE(Cstring,1)-1) THEN
      LENGTH=SIZE(Cstring,1)-1
    ELSE
      LENGTH=LEN(Fstring)
    ENDIF
    Fstring=""
    DO i=1,LENGTH
      IF(Cstring(i)==C_NULL_CHAR) THEN
        EXIT
      ELSE
        Fstring(i:i)=Cstring(i)
      ENDIF
    ENDDO !i
    
    RETURN
    
  END SUBROUTINE CMISSC2FSTRING
   
  !
  !================================================================================================================================
  !

  !>Copys/converts a Fortran String (length of characters) to a C string (array of characters)
  SUBROUTINE CMISSF2CString(Fstring,Cstring)
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: Fstring
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: Cstring(:)
    !Local variables
    INTEGER(C_INT) :: i,LENGTH

    IF(SIZE(Cstring,1)>LEN_TRIM(Fstring)) THEN
      LENGTH=LEN_TRIM(Fstring)
    ELSE
      LENGTH=SIZE(Cstring,1)-1
    ENDIF
    DO i=1,LENGTH     
      Cstring(i)=Fstring(i:i)
    ENDDO !i
    !Null terminate the string
    Cstring(LENGTH+1)=C_NULL_CHAR
    
    RETURN
    
  END SUBROUTINE CMISSF2CSTRING
   
  !
  !================================================================================================================================
  !
 
END MODULE CMISS_FORTRAN_C
