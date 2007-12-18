!> \file
!> $Id: strings.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module contains all string manipulation and transformation routines.
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

!> This module contains all string manipulation and transformation routines.
MODULE STRINGS

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  
  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Interfaces

  INTERFACE IS_ABBREVIATION
    MODULE PROCEDURE IS_ABBREVIATION_C_C
    MODULE PROCEDURE IS_ABBREVIATION_C_VS
    MODULE PROCEDURE IS_ABBREVIATION_VS_C
    MODULE PROCEDURE IS_ABBREVIATION_VS_VS
  END INTERFACE !IS_ABBREVIATION

  INTERFACE LIST_TO_CHARACTER
    MODULE PROCEDURE LIST_TO_CHARACTER_C
    MODULE PROCEDURE LIST_TO_CHARACTER_INTG
    MODULE PROCEDURE LIST_TO_CHARACTER_LINTG
    MODULE PROCEDURE LIST_TO_CHARACTER_L
    MODULE PROCEDURE LIST_TO_CHARACTER_SP
    MODULE PROCEDURE LIST_TO_CHARACTER_DP
  END INTERFACE !LIST_TO_CHARACTER

  INTERFACE NUMBER_TO_CHARACTER
    MODULE PROCEDURE NUMBER_TO_CHARACTER_INTG
    MODULE PROCEDURE NUMBER_TO_CHARACTER_LINTG
    MODULE PROCEDURE NUMBER_TO_CHARACTER_SP
    MODULE PROCEDURE NUMBER_TO_CHARACTER_DP
  END INTERFACE !NUMBER_TO_CHARACTER

  INTERFACE NUMBER_TO_VSTRING
    MODULE PROCEDURE NUMBER_TO_VSTRING_INTG
    MODULE PROCEDURE NUMBER_TO_VSTRING_LINTG
    MODULE PROCEDURE NUMBER_TO_VSTRING_SP
    MODULE PROCEDURE NUMBER_TO_VSTRING_DP
  END INTERFACE !NUMBER_TO_VSTRING

  INTERFACE STRING_TO_DOUBLE
    MODULE PROCEDURE STRING_TO_DOUBLE_C
    MODULE PROCEDURE STRING_TO_DOUBLE_VS
  END INTERFACE !STRING_TO_DOUBLE

  INTERFACE STRING_TO_INTEGER
    MODULE PROCEDURE STRING_TO_INTEGER_C
    MODULE PROCEDURE STRING_TO_INTEGER_VS
  END INTERFACE !STRING_TO_INTEGER

  INTERFACE STRING_TO_LONG_INTEGER
    MODULE PROCEDURE STRING_TO_LONG_INTEGER_C
    MODULE PROCEDURE STRING_TO_LONG_INTEGER_VS
  END INTERFACE !STRING_TO_LONG_INTEGER

  INTERFACE STRING_TO_LOGICAL
    MODULE PROCEDURE STRING_TO_LOGICAL_C
    MODULE PROCEDURE STRING_TO_LOGICAL_VS
  END INTERFACE !STRING_TO_LOGICAL

  INTERFACE STRING_TO_SINGLE
    MODULE PROCEDURE STRING_TO_SINGLE_C
    MODULE PROCEDURE STRING_TO_SINGLE_VS
  END INTERFACE !STRING_TO_SINGLE

  INTERFACE CHARACTER_TO_LOWERCASE
    MODULE PROCEDURE CHARACTER_TO_LOWERCASE_C
    MODULE PROCEDURE CHARACTER_TO_LOWERCASE_VS
  END INTERFACE !CHARACTER_TO_LOWERCASE

  INTERFACE VSTRING_TO_LOWERCASE
    MODULE PROCEDURE VSTRING_TO_LOWERCASE_C
    MODULE PROCEDURE VSTRING_TO_LOWERCASE_VS
  END INTERFACE !VSTRING_TO_LOWERCASE

  INTERFACE CHARACTER_TO_UPPERCASE
    MODULE PROCEDURE CHARACTER_TO_UPPERCASE_C
    MODULE PROCEDURE CHARACTER_TO_UPPERCASE_VS
  END INTERFACE !CHARACTER_TO_UPPERCASE

  INTERFACE VSTRING_TO_UPPERCASE
    MODULE PROCEDURE VSTRING_TO_UPPERCASE_C
    MODULE PROCEDURE VSTRING_TO_UPPERCASE_VS
  END INTERFACE !VSTRING_TO_UPPERCASE

  PUBLIC IS_ABBREVIATION,IS_DIGIT,IS_LETTER,IS_WHITESPACE,LIST_TO_CHARACTER,LOGICAL_TO_CHARACTER,LOGICAL_TO_VSTRING, &
    & NUMBER_TO_CHARACTER,NUMBER_TO_VSTRING,STRING_TO_DOUBLE,STRING_TO_INTEGER,STRING_TO_LONG_INTEGER,STRING_TO_LOGICAL, &
    & STRING_TO_SINGLE,CHARACTER_TO_LOWERCASE,VSTRING_TO_LOWERCASE,CHARACTER_TO_UPPERCASE,VSTRING_TO_UPPERCASE
  
CONTAINS
  
  !
  !================================================================================================================================
  !
 
  !#### Generic-Function: IS_ABBREVIATION
  !###  Type: LOGICAL
  !###  Description:
  !###    Returns .TRUE. if a supplied string is a valid abbreviation of a second supplied string.
  !###  Child-functions: IS_ABBREVIATION_C_C,IS_ABBREVIATION_C_VS,IS_ABBREVIATION_VS_C,IS_ABBREVIATION_VS_VS

  !
  !================================================================================================================================
  !
  
  FUNCTION IS_ABBREVIATION_C_C(SHORT,LONG,MIN_NUM_CHARACTERS)

    !#### Function: IS_ABBREVIATION_C_C
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_ABBREVIATION returns .TRUE. if the character string SHORT is an abbreviation of the character string LONG.
    !###    SHORT must be at least MIN_NUM_CHARACTERS long.
    !###  Parent-function: IS_ABBREVIATION

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: SHORT, LONG
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS
    !Function variable
    LOGICAL :: IS_ABBREVIATION_C_C
    !Local Variables
    INTEGER(INTG) :: noch,NUM_CHARACTERS
    CHARACTER(LEN=LEN(SHORT)) :: UPPER_SHORT
    CHARACTER(LEN=LEN(LONG)) :: UPPER_LONG
    
    IS_ABBREVIATION_C_C=.FALSE.
    UPPER_SHORT=CHARACTER_TO_UPPERCASE(SHORT)
    UPPER_LONG=CHARACTER_TO_UPPERCASE(LONG)
    NUM_CHARACTERS=MIN(LEN(LONG),LEN(SHORT))
    DO noch=MIN_NUM_CHARACTERS,NUM_CHARACTERS
      IF(UPPER_SHORT==UPPER_LONG(:noch)) THEN
          IS_ABBREVIATION_C_C=.TRUE.
          EXIT
      ENDIF
    ENDDO !noch

    RETURN
  END FUNCTION IS_ABBREVIATION_C_C

  !
  !================================================================================================================================
  !
  
  FUNCTION IS_ABBREVIATION_C_VS(SHORT,LONG,MIN_NUM_CHARACTERS)

    !#### Function: IS_ABBREVIATION_C_VS
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_ABBREVIATION returns .TRUE. if the character string SHORT is an abbreviation of the varying string LONG.
    !###    SHORT must be at least MIN_NUM_CHARACTERS long.
    !###  Parent-function: IS_ABBREVIATION

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: SHORT
    TYPE(VARYING_STRING), INTENT(IN) :: LONG
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS
    !Function variable
    LOGICAL :: IS_ABBREVIATION_C_VS
    !Local Variables
    INTEGER(INTG) :: noch,NUM_CHARACTERS
    CHARACTER(LEN=LEN(SHORT)) :: UPPER_SHORT
    TYPE(VARYING_STRING) :: UPPER_LONG
    
    IS_ABBREVIATION_C_VS=.FALSE.
    UPPER_SHORT=CHARACTER_TO_UPPERCASE(SHORT)
    UPPER_LONG=VSTRING_TO_UPPERCASE(LONG)
    NUM_CHARACTERS=MIN(LEN(LONG),LEN(SHORT))
    DO noch=MIN_NUM_CHARACTERS,NUM_CHARACTERS
      IF(UPPER_SHORT==EXTRACT(UPPER_LONG,1,noch)) THEN
          IS_ABBREVIATION_C_VS=.TRUE.
          EXIT
      ENDIF
    ENDDO !noch

    RETURN
  END FUNCTION IS_ABBREVIATION_C_VS

  !
  !================================================================================================================================
  !
  
  FUNCTION IS_ABBREVIATION_VS_C(SHORT,LONG,MIN_NUM_CHARACTERS)

    !#### Function: IS_ABBREVIATION_VS_C
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_ABBREVIATION returns .TRUE. if the varying string SHORT is an abbreviation of the character string LONG.
    !###    SHORT must be at least MIN_NUM_CHARACTERS long.
    !###  Parent-function: IS_ABBREVIATION

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: SHORT
    CHARACTER(LEN=*), INTENT(IN) :: LONG
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS
    !Function variable
    LOGICAL :: IS_ABBREVIATION_VS_C
    !Local Variables
    INTEGER(INTG) :: noch,NUM_CHARACTERS
    TYPE(VARYING_STRING) :: UPPER_SHORT
    CHARACTER(LEN=LEN(LONG)) :: UPPER_LONG
    
    IS_ABBREVIATION_VS_C=.FALSE.
    UPPER_SHORT=VSTRING_TO_UPPERCASE(SHORT)
    UPPER_LONG=CHARACTER_TO_UPPERCASE(LONG)
    NUM_CHARACTERS=MIN(LEN(LONG),LEN(SHORT))
    DO noch=MIN_NUM_CHARACTERS,NUM_CHARACTERS
      IF(UPPER_SHORT==UPPER_LONG(:noch)) THEN
          IS_ABBREVIATION_VS_C=.TRUE.
          EXIT
      ENDIF
    ENDDO !noch

    RETURN
  END FUNCTION IS_ABBREVIATION_VS_C

  !
  !================================================================================================================================
  !
  
  FUNCTION IS_ABBREVIATION_VS_VS(SHORT,LONG,MIN_NUM_CHARACTERS)

    !#### Function: IS_ABBREVIATION_VS_VS
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_ABBREVIATION returns .TRUE. if the varying string SHORT is an abbreviation of the varying string LONG.
    !###    SHORT must be at least MIN_NUM_CHARACTERS long.
    !###  Parent-function: IS_ABBREVIATION

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: SHORT,LONG
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS
    !Function variable
    LOGICAL :: IS_ABBREVIATION_VS_VS
    !Local Variables
    INTEGER(INTG) :: noch,NUM_CHARACTERS
    TYPE(VARYING_STRING) :: UPPER_SHORT,UPPER_LONG
    
    IS_ABBREVIATION_VS_VS=.FALSE.
    UPPER_SHORT=VSTRING_TO_UPPERCASE(SHORT)
    UPPER_LONG=VSTRING_TO_UPPERCASE(LONG)
    NUM_CHARACTERS=MIN(LEN(LONG),LEN(SHORT))
    DO noch=MIN_NUM_CHARACTERS,NUM_CHARACTERS
      IF(UPPER_SHORT==EXTRACT(UPPER_LONG,1,noch)) THEN
          IS_ABBREVIATION_VS_VS=.TRUE.
          EXIT
      ENDIF
    ENDDO !noch

    RETURN
  END FUNCTION IS_ABBREVIATION_VS_VS

  !
  !================================================================================================================================
  !

  FUNCTION IS_DIGIT(CHARAC)

    !#### Function: IS_DIGIT
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_DIGIT returns .TRUE. if the character CHARAC is a digit character (i.e. 0..9)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARAC
    !Function variable
    LOGICAL :: IS_DIGIT
    !Local Variables

    IS_DIGIT=(ICHAR(CHARAC)>=ICHAR("0").AND.ICHAR(CHARAC)<=ICHAR("9"))

    RETURN
  END FUNCTION IS_DIGIT

  !
  !================================================================================================================================
  !

  FUNCTION IS_LETTER(CHARAC)

    !#### Function: IS_LETTER
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_LETTER returns .TRUE. if the character CHARAC is a letter character (i.e. A..Z or a..z)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARAC
    !Function variable
    LOGICAL :: IS_LETTER
    !Local Variables

    IS_LETTER=((ICHAR(CHARAC)>=ICHAR("A").AND.ICHAR(CHARAC)<=ICHAR("Z")).OR.&
	    & (ICHAR(CHARAC)>=ICHAR("a").AND.ICHAR(CHARAC)<=ICHAR("z")))

    RETURN
  END FUNCTION IS_LETTER

  !
  !================================================================================================================================
  !

  FUNCTION IS_LOWERCASE(CHARC)

    !#### Function: IS_LOWERCASE
    !###  Type: LOGICAL
    !###  Description:
    !###    Returns .TRUE. if the supplied character is a lowercase character.

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARC
    !Function variable
    LOGICAL :: IS_LOWERCASE
    !Local Variables

    IF(LGE(CHARC,"a").AND.LLE(CHARC,"z")) THEN
      IS_LOWERCASE=.TRUE.
    ELSE
      IS_LOWERCASE=.FALSE.
    ENDIF

    RETURN
  END FUNCTION IS_LOWERCASE

  !
  !================================================================================================================================
  !

  FUNCTION IS_UPPERCASE(CHARC)

    !#### Function: IS_UPPERCASE
    !###  Type: LOGICAL
    !###  Description:
    !###    Returns .TRUE. if the supplied character is an uppercase character.

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARC
    !Function variable
    LOGICAL :: IS_UPPERCASE
    !Local Variables

    IF(LGE(CHARC,"A").AND.LLE(CHARC,"Z")) THEN
      IS_UPPERCASE=.TRUE.
    ELSE
      IS_UPPERCASE=.FALSE.
    ENDIF

    RETURN
  END FUNCTION IS_UPPERCASE

  !
  !================================================================================================================================
  !

  FUNCTION IS_WHITESPACE(CHARAC)

    !#### Function: IS_WHITESPACE
    !###  Type: LOGICAL
    !###  Description:
    !###    IS_WHITESPACE returns .TRUE. if the character CHARAC is a whitespace character (i.e. space, tabs, etc.)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARAC
    !Function variable
    LOGICAL :: IS_WHITESPACE
    !Local Variables
    
    !!WARNING: Assumes ASCII encoding
    IS_WHITESPACE=(CHARAC==CHAR(32).OR.CHARAC==CHAR(9))

    RETURN
  END FUNCTION IS_WHITESPACE

  !
  !================================================================================================================================
  !
 
  !#### Generic-Function: LIST_TO_CHARACTER
  !###  Type: CHARACTER(LEN=MAXSTRLEN)
  !###  Description:
  !###    Converts a list to its equivalent character string representation.
  !###  Child-functions: LIST_TO_CHARACTER_C,LIST_TO_CHARACTER_INTG,LIST_TO_CHARACTER_LINTG,LIST_TO_CHARACTER_L, &
  !###    LIST_TO_CHARACTER_SP,LIST_TO_CHARACTER_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION LIST_TO_CHARACTER_C(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR,LIST_LENGTHS)
  
    !#### Function: LIST_TO_CHARACTER_C
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a character list to its equivalent character string representation as determined by the supplied format.
    !###    If present, the optional argument LIST_LENGTHS is used for the lengths of each list elements length otherwise the
    !###    trimmed length is used. NOTE: The FORMAT is ignored for this child FUNCTION.
    !###  Parent-function: LIST_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST
    CHARACTER(LEN=*), INTENT(IN) :: LIST(NUMBER_IN_LIST),FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    INTEGER(INTG), OPTIONAL, INTENT(IN) :: LIST_LENGTHS(NUMBER_IN_LIST)
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_C
    !Local variables
    INTEGER(INTG) :: i,POSITION,LENGTH
    
    CALL ENTERS("LIST_TO_CHARACTER_C",ERR,ERROR,*999)

    LIST_TO_CHARACTER_C=""
    IF(NUMBER_IN_LIST>0) THEN
      IF(PRESENT(LIST_LENGTHS)) THEN
        LENGTH=LIST_LENGTHS(1)
        LIST_TO_CHARACTER_C=LIST(1)(1:LENGTH)
        DO i=2,NUMBER_IN_LIST
          IF(LENGTH+LIST_LENGTHS(i)+1<=MAXSTRLEN) THEN
            LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:LENGTH)//","//LIST(i)(1:LIST_LENGTHS(i))
            LENGTH=LENGTH+LIST_LENGTHS(i)+1
          ELSE IF(LENGTH+5<=MAXSTRLEN) THEN
            LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:LENGTH)//",...."
            EXIT
          ELSE
            POSITION=INDEX(LIST_TO_CHARACTER_C(1:MAXSTRLEN-4),",",.TRUE.)
            IF(POSITION/=0) THEN
              LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:POSITION)//"...."
            ELSE
              LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:MAXSTRLEN-5)//",...."
            ENDIF
            EXIT
          ENDIF
        ENDDO !i
      ELSE
        LIST_TO_CHARACTER_C=LIST(1)(1:LEN_TRIM(LIST(1)))
        DO i=2,NUMBER_IN_LIST
          IF(LEN_TRIM(LIST_TO_CHARACTER_C)+LEN_TRIM(LIST(i))+1<=MAXSTRLEN) THEN
            LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:LEN_TRIM(LIST_TO_CHARACTER_C))//","//LIST(i)(1:LEN_TRIM(LIST(i)))
          ELSE IF(LEN_TRIM(LIST_TO_CHARACTER_C)+5<=MAXSTRLEN) THEN
            LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:LEN_TRIM(LIST_TO_CHARACTER_C))//",...."
            EXIT
          ELSE
            POSITION=INDEX(LIST_TO_CHARACTER_C(1:MAXSTRLEN-4),",",.TRUE.)
            IF(POSITION/=0) THEN
              LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:POSITION)//"...."
            ELSE
              LIST_TO_CHARACTER_C=LIST_TO_CHARACTER_C(1:MAXSTRLEN-5)//",...."
            ENDIF
            EXIT
          ENDIF
        ENDDO !i
      ENDIF
    ENDIF

    CALL EXITS("LIST_TO_CHARACTER_C")
    RETURN
999 CALL ERRORS("LIST_TO_CHARACTER_C",ERR,ERROR)
    CALL EXITS("LIST_TO_CHARACTER_C")
    RETURN    
  END FUNCTION LIST_TO_CHARACTER_C
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LIST_TO_CHARACTER_INTG(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !#### Function: LIST_TO_CHARACTER_INTG
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts an integer list to its equivalent character string representation as determined by the supplied format. 
    !###  Parent-function: LIST_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST,LIST(NUMBER_IN_LIST)
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_INTG
    !Local variables
    INTEGER(INTG) :: i,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: LIST_VALUE
    
    CALL ENTERS("LIST_TO_CHARACTER_INTG",ERR,ERROR,*999)

    LIST_TO_CHARACTER_INTG=""
    IF(NUMBER_IN_LIST>0) THEN
      LIST_TO_CHARACTER_INTG=NUMBER_TO_CHARACTER_INTG(LIST(1),FORMAT,ERR,ERROR)
      IF(ERR/=0) GOTO 999
      DO i=2,NUMBER_IN_LIST
        LIST_VALUE=NUMBER_TO_CHARACTER_INTG(LIST(i),FORMAT,ERR,ERROR)
        IF(ERR/=0) GOTO 999
        IF(LEN_TRIM(LIST_TO_CHARACTER_INTG)+LEN_TRIM(LIST_VALUE)+1<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_INTG=LIST_TO_CHARACTER_INTG(1:LEN_TRIM(LIST_TO_CHARACTER_INTG))//","// &
            & LIST_VALUE(1:LEN_TRIM(LIST_VALUE))
        ELSE IF(LEN_TRIM(LIST_TO_CHARACTER_INTG)+5<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_INTG=LIST_TO_CHARACTER_INTG(1:LEN_TRIM(LIST_TO_CHARACTER_INTG))//",...."
          EXIT
        ELSE
          POSITION=INDEX(LIST_TO_CHARACTER_INTG(1:MAXSTRLEN-4),",",.TRUE.)
          IF(POSITION/=0) THEN
            LIST_TO_CHARACTER_INTG=LIST_TO_CHARACTER_INTG(1:POSITION)//"...."
          ELSE
            LIST_TO_CHARACTER_INTG=LIST_TO_CHARACTER_INTG(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    CALL EXITS("LIST_TO_CHARACTER_INTG")
    RETURN
999 CALL ERRORS("LIST_TO_CHARACTER_INTG",ERR,ERROR)
    CALL EXITS("LIST_TO_CHARACTER_INTG")
    RETURN    
  END FUNCTION LIST_TO_CHARACTER_INTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LIST_TO_CHARACTER_LINTG(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !#### Function: LIST_TO_CHARACTER_LINTG
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts an long integer list to its equivalent character string representation as determined by the supplied format. 
    !###  Parent-function: LIST_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST
    INTEGER(LINTG), INTENT(IN) :: LIST(NUMBER_IN_LIST)
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_LINTG
    !Local variables
    INTEGER(INTG) :: i,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: LIST_VALUE
    
    CALL ENTERS("LIST_TO_CHARACTER_LINTG",ERR,ERROR,*999)

    LIST_TO_CHARACTER_LINTG=""
    IF(NUMBER_IN_LIST>0) THEN
      LIST_TO_CHARACTER_LINTG=NUMBER_TO_CHARACTER_LINTG(LIST(1),FORMAT,ERR,ERROR)
      IF(ERR/=0) GOTO 999
      DO i=2,NUMBER_IN_LIST
        LIST_VALUE=NUMBER_TO_CHARACTER_LINTG(LIST(i),FORMAT,ERR,ERROR)
        IF(ERR/=0) GOTO 999
        IF(LEN_TRIM(LIST_TO_CHARACTER_LINTG)+LEN_TRIM(LIST_VALUE)+1<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_LINTG=LIST_TO_CHARACTER_LINTG(1:LEN_TRIM(LIST_TO_CHARACTER_LINTG))//","// &
            & LIST_VALUE(1:LEN_TRIM(LIST_VALUE))
        ELSE IF(LEN_TRIM(LIST_TO_CHARACTER_LINTG)+5<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_LINTG=LIST_TO_CHARACTER_LINTG(1:LEN_TRIM(LIST_TO_CHARACTER_LINTG))//",...."
          EXIT
        ELSE
          POSITION=INDEX(LIST_TO_CHARACTER_LINTG(1:MAXSTRLEN-4),",",.TRUE.)
          IF(POSITION/=0) THEN
            LIST_TO_CHARACTER_LINTG=LIST_TO_CHARACTER_LINTG(1:POSITION)//"...."
          ELSE
            LIST_TO_CHARACTER_LINTG=LIST_TO_CHARACTER_LINTG(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    CALL EXITS("LIST_TO_CHARACTER_LINTG")
    RETURN
999 CALL ERRORS("LIST_TO_CHARACTER_LINTG",ERR,ERROR)
    CALL EXITS("LIST_TO_CHARACTER_LINTG")
    RETURN    
  END FUNCTION LIST_TO_CHARACTER_LINTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LIST_TO_CHARACTER_L(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !#### Function: LIST_TO_CHARACTER_L
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a logical list to its equivalent character string representation as determined by the supplied format string. 
    !###    NOTE: The FORMAT is ignored for this child FUNCTION.
    !###  Parent-function: LIST_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST
    LOGICAL, INTENT(IN) :: LIST(NUMBER_IN_LIST)    
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_L
    !Local variables
    INTEGER(INTG) :: i,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: LIST_VALUE
    
    CALL ENTERS("LIST_TO_CHARACTER_L",ERR,ERROR,*999)

    LIST_TO_CHARACTER_L=""
    IF(NUMBER_IN_LIST>0) THEN
      LIST_TO_CHARACTER_L=LOGICAL_TO_CHARACTER(LIST(1),ERR,ERROR)
      IF(ERR/=0) GOTO 999
      DO i=2,NUMBER_IN_LIST
        LIST_VALUE=LOGICAL_TO_CHARACTER(LIST(i),ERR,ERROR)
        IF(ERR/=0) GOTO 999
        IF(LEN_TRIM(LIST_TO_CHARACTER_L)+LEN_TRIM(LIST_VALUE)+1<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_L=LIST_TO_CHARACTER_L(1:LEN_TRIM(LIST_TO_CHARACTER_L))//","//LIST_VALUE(1:LEN_TRIM(LIST_VALUE))
        ELSE IF(LEN_TRIM(LIST_TO_CHARACTER_L)+5<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_L=LIST_TO_CHARACTER_L(1:LEN_TRIM(LIST_TO_CHARACTER_L))//",...."
          EXIT
        ELSE
          POSITION=INDEX(LIST_TO_CHARACTER_L(1:MAXSTRLEN-4),",",.TRUE.)
          IF(POSITION/=0) THEN
            LIST_TO_CHARACTER_L=LIST_TO_CHARACTER_L(1:POSITION)//"...."
          ELSE
            LIST_TO_CHARACTER_L=LIST_TO_CHARACTER_L(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    CALL EXITS("LIST_TO_CHARACTER_L")
    RETURN
999 CALL ERRORS("LIST_TO_CHARACTER_L",ERR,ERROR)
    CALL EXITS("LIST_TO_CHARACTER_L")
    RETURN    
  END FUNCTION LIST_TO_CHARACTER_L
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LIST_TO_CHARACTER_SP(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !#### Function: LIST_TO_CHARACTER_SP
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a single precision list to its equivalent character string representation as determined by the supplied
    !###    format string.
    !###  Parent-function: LIST_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST
    REAL(SP), INTENT(IN) :: LIST(NUMBER_IN_LIST)    
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_SP
    !Local variables
    INTEGER(INTG) :: i,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: LIST_VALUE
    
    CALL ENTERS("LIST_TO_CHARACTER_SP",ERR,ERROR,*999)

    LIST_TO_CHARACTER_SP=""
    IF(NUMBER_IN_LIST>0) THEN
      LIST_TO_CHARACTER_SP=NUMBER_TO_CHARACTER_SP(LIST(1),FORMAT,ERR,ERROR)
      IF(ERR/=0) GOTO 999
      DO i=2,NUMBER_IN_LIST
        LIST_VALUE=NUMBER_TO_CHARACTER_SP(LIST(i),FORMAT,ERR,ERROR)
        IF(ERR/=0) GOTO 999
        IF(LEN_TRIM(LIST_TO_CHARACTER_SP)+LEN_TRIM(LIST_VALUE)+1<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_SP=LIST_TO_CHARACTER_SP(1:LEN_TRIM(LIST_TO_CHARACTER_SP))//","//LIST_VALUE(1:LEN_TRIM(LIST_VALUE))
        ELSE IF(LEN_TRIM(LIST_TO_CHARACTER_SP)+5<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_SP=LIST_TO_CHARACTER_SP(1:LEN_TRIM(LIST_TO_CHARACTER_SP))//",...."
          EXIT
        ELSE
          POSITION=INDEX(LIST_TO_CHARACTER_SP(1:MAXSTRLEN-4),",",.TRUE.)
          IF(POSITION/=0) THEN
            LIST_TO_CHARACTER_SP=LIST_TO_CHARACTER_SP(1:POSITION)//"...."
          ELSE
            LIST_TO_CHARACTER_SP=LIST_TO_CHARACTER_SP(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    CALL EXITS("LIST_TO_CHARACTER_SP")
    RETURN
999 CALL ERRORS("LIST_TO_CHARACTER_SP",ERR,ERROR)
    CALL EXITS("LIST_TO_CHARACTER_SP")
    RETURN    
  END FUNCTION LIST_TO_CHARACTER_SP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LIST_TO_CHARACTER_DP(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !#### Function: LIST_TO_CHARACTER_DP
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a double precision list to its equivalent character string representation as determined by the supplied
    !###    format string.
    !###  Parent-function: LIST_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST
    REAL(DP), INTENT(IN) :: LIST(NUMBER_IN_LIST)    
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_DP
    !Local variables
    INTEGER(INTG) :: i,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: LIST_VALUE
    
    CALL ENTERS("LIST_TO_CHARACTER_DP",ERR,ERROR,*999)

    LIST_TO_CHARACTER_DP=""
    IF(NUMBER_IN_LIST>0) THEN
      LIST_TO_CHARACTER_DP=NUMBER_TO_CHARACTER_DP(LIST(1),FORMAT,ERR,ERROR)
      IF(ERR/=0) GOTO 999
      DO i=2,NUMBER_IN_LIST
        LIST_VALUE=NUMBER_TO_CHARACTER_DP(LIST(i),FORMAT,ERR,ERROR)
        IF(ERR/=0) GOTO 999
        IF(LEN_TRIM(LIST_TO_CHARACTER_DP)+LEN_TRIM(LIST_VALUE)+1<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_DP=LIST_TO_CHARACTER_DP(1:LEN_TRIM(LIST_TO_CHARACTER_DP))//","//LIST_VALUE(1:LEN_TRIM(LIST_VALUE))
        ELSE IF(LEN_TRIM(LIST_TO_CHARACTER_DP)+5<=MAXSTRLEN) THEN
          LIST_TO_CHARACTER_DP=LIST_TO_CHARACTER_DP(1:LEN_TRIM(LIST_TO_CHARACTER_DP))//",...."
          EXIT
        ELSE
          POSITION=INDEX(LIST_TO_CHARACTER_DP(1:MAXSTRLEN-4),",",.TRUE.)
          IF(POSITION/=0) THEN
            LIST_TO_CHARACTER_DP=LIST_TO_CHARACTER_DP(1:POSITION)//"...."
          ELSE
            LIST_TO_CHARACTER_DP=LIST_TO_CHARACTER_DP(1:MAXSTRLEN-5)//",...."
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

    CALL EXITS("LIST_TO_CHARACTER_DP")
    RETURN
999 CALL ERRORS("LIST_TO_CHARACTER_DP",ERR,ERROR)
    CALL EXITS("LIST_TO_CHARACTER_DP")
    RETURN    
  END FUNCTION LIST_TO_CHARACTER_DP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LOGICAL_TO_CHARACTER(LOGICALVALUE,ERR,ERROR)
  
    !#### Function: LOGICAL_TO_CHARACTER
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a logical value to either a "TRUE" or "FALSE" character string.
    
    !Argument variables
    LOGICAL, INTENT(IN) :: LOGICALVALUE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LOGICAL_TO_CHARACTER
    !Local variables
    
    CALL ENTERS("LOGICAL_TO_CHARACTER",ERR,ERROR,*999)

    IF(LOGICALVALUE) THEN
      LOGICAL_TO_CHARACTER="TRUE"
    ELSE
      LOGICAL_TO_CHARACTER="FALSE"
    ENDIF

    CALL EXITS("LOGICAL_TO_CHARACTER")
    RETURN
999 CALL ERRORS("LOGICAL_TO_CHARACTER",ERR,ERROR)
    CALL EXITS("LOGICAL_TO_CHARACTER")
    RETURN    
  END FUNCTION LOGICAL_TO_CHARACTER
  
  !
  !================================================================================================================================
  !
  
  FUNCTION LOGICAL_TO_VSTRING(LOGICALVALUE,ERR,ERROR)
  
    !#### Function: LOGICAL_TO_CHARACTER
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Converts a logical value to either a "TRUE" or "FALSE" varying string.
    
    !Argument variables
    LOGICAL, INTENT(IN) :: LOGICALVALUE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    TYPE(VARYING_STRING) :: LOGICAL_TO_VSTRING
    !Local variables
    
    CALL ENTERS("LOGICAL_TO_VSTRING",ERR,ERROR,*999)

    IF(LOGICALVALUE) THEN
      LOGICAL_TO_VSTRING="TRUE"
    ELSE
      LOGICAL_TO_VSTRING="FALSE"
    ENDIF

    CALL EXITS("LOGICAL_TO_VSTRING")
    RETURN
999 CALL ERRORS("LOGICAL_TO_VSTRING",ERR,ERROR)
    CALL EXITS("LOGICAL_TO_VSTRING")
    RETURN    
  END FUNCTION LOGICAL_TO_VSTRING
  
  !
  !================================================================================================================================
  !
 
  !#### Generic-Function: NUMBER_TO_CHARACTER
  !###  Type: CHARACTER(LEN=MAXSTRLEN)
  !###  Description:
  !###    Converts a number to its equivalent character string representation.
  !###  Child-functions: NUMBER_TO_CHARACTER_INTG,NUMBER_TO_CHARACTER_LINTG,NUMBER_TO_CHARACTER_SP,NUMBER_TO_CHARACTER_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_CHARACTER_INTG(NUMBER,FORMAT,ERR,ERROR)
  
    !#### Function: NUMBER_TO_CHARACTER_INTG
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts an integer number to its equivalent character string representation as determined by the supplied format.
    !###    The format is of the form of a standard Fortran integer format e.g. "I3".
    !###  Parent-function: NUMBER_TO_CHARACTER
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_INTG
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT
    
    CALL ENTERS("NUMBER_TO_CHARACTER_INTG",ERR,ERROR,*999)

    IF(FORMAT(1:1)=="*") THEN
      LOCAL_FORMAT="(I12)"
    ELSE
      LOCAL_FORMAT="("//FORMAT(1:LEN_TRIM(FORMAT))//")"
    ENDIF
    WRITE(NUMBER_TO_CHARACTER_INTG,LOCAL_FORMAT,ERR=999) NUMBER

    !Trim leading blanks
    NUMBER_TO_CHARACTER_INTG=ADJUSTL(NUMBER_TO_CHARACTER_INTG)

    CALL EXITS("NUMBER_TO_CHARACTER_INTG")
    RETURN
999 CALL FLAG_ERROR("Error converting an integer to a character string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_CHARACTER_INTG",ERR,ERROR)
    CALL EXITS("NUMBER_TO_CHARACTER_INTG")
    RETURN    
  END FUNCTION NUMBER_TO_CHARACTER_INTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_CHARACTER_LINTG(NUMBER,FORMAT,ERR,ERROR)
  
    !#### Function: NUMBER_TO_CHARACTER_LINTG
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a long integer number to its equivalent character string representation as determined by the supplied format.
    !###    The format is of the form of a standard Fortran integer format e.g. "I3".
    !###  Parent-function: NUMBER_TO_CHARACTER
    
    !Argument variables
    INTEGER(LINTG), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_LINTG
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT
    
    CALL ENTERS("NUMBER_TO_CHARACTER_LINTG",ERR,ERROR,*999)

    IF(FORMAT(1:1)=="*") THEN
      LOCAL_FORMAT="(I18)"
    ELSE
      LOCAL_FORMAT="("//FORMAT(1:LEN_TRIM(FORMAT))//")"
    ENDIF
    WRITE(NUMBER_TO_CHARACTER_LINTG,LOCAL_FORMAT,ERR=999) NUMBER

    !Trim leading blanks
    NUMBER_TO_CHARACTER_LINTG=ADJUSTL(NUMBER_TO_CHARACTER_LINTG)

    CALL EXITS("NUMBER_TO_CHARACTER_LINTG")
    RETURN
999 CALL FLAG_ERROR("Error converting a long integer to a character string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_CHARACTER_LINTG",ERR,ERROR)
    CALL EXITS("NUMBER_TO_CHARACTER_LINTG")
    RETURN    
  END FUNCTION NUMBER_TO_CHARACTER_LINTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_CHARACTER_SP(NUMBER, FORMAT, ERR, ERROR)
  
    !#### Function: NUMBER_TO_CHARACTER_SP
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a single precision number to its equivalent character string representation as determined by the supplied
    !###    format string. NOTE: If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to
    !###    maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the
    !###    supplied number in either F8.? or E8.? format depending on its magnitude.
    !###  Parent-function: NUMBER_TO_CHARACTER
    
    !Argument variables
    REAL(SP), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_SP
    !Local variables
    INTEGER(INTG) :: ASTERISK_POS,i0,i1,LENGTH
    CHARACTER(LEN=MAXSTRLEN) :: CI0,CI1
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT
    
    CALL ENTERS("NUMBER_TO_CHARACTER_SP",ERR,ERROR,*999)

    ASTERISK_POS=INDEX(FORMAT,"*")
    LENGTH=LEN_TRIM(FORMAT)
    IF(ASTERISK_POS==1.AND.LENGTH==1) THEN !Free format
      WRITE(NUMBER_TO_CHARACTER_SP,*,ERR=999) NUMBER      
    ELSE IF(ASTERISK_POS>0) THEN !Adjustable format
      CI0=FORMAT(ASTERISK_POS+1:LEN_TRIM(FORMAT))
      READ(CI0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(NUMBER>=0.0_SP) THEN
          IF((NUMBER<10.0_SP**(i0-1)).AND.(NUMBER>=0.1_SP**(MIN(i0-2,5)))) THEN
            IF(NUMBER>1.0_SP) THEN
              i1=i0-2-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ELSE
          IF((-NUMBER<10.0_SP**(i0-2)).AND.(-NUMBER>=0.01_SP**(MIN(i0-2,5)))) THEN
            IF(-NUMBER>=1.0_SP) THEN
              i1=i0-3-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,'(I2)') i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid FORMAT",ERR,ERROR,*999)
        GOTO 999
      ENDIF
    ELSE
      LOCAL_FORMAT='('//FORMAT(1:LEN_TRIM(FORMAT))//')'
      WRITE(NUMBER_TO_CHARACTER_SP,LOCAL_FORMAT,ERR=999) NUMBER
    ENDIF

    !Add an extra zero if required
    IF(NUMBER_TO_CHARACTER_SP(LEN_TRIM(NUMBER_TO_CHARACTER_SP):LEN_TRIM(NUMBER_TO_CHARACTER_SP))==".") &
      & NUMBER_TO_CHARACTER_SP=NUMBER_TO_CHARACTER_SP(1:LEN_TRIM(NUMBER_TO_CHARACTER_SP))//"0"
    !Trim leading blanks
    NUMBER_TO_CHARACTER_SP=ADJUSTL(NUMBER_TO_CHARACTER_SP)

    CALL EXITS("NUMBER_TO_CHARACTER_SP")
    RETURN
999 CALL FLAG_ERROR("Error converting a single precision number to a character string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_CHARACTER_SP",ERR,ERROR)
    CALL EXITS("NUMBER_TO_CHARACTER_SP")
    RETURN    
  END FUNCTION NUMBER_TO_CHARACTER_SP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_CHARACTER_DP(NUMBER, FORMAT, ERR, ERROR)
  
    !#### Function: NUMBER_TO_CHARACTER_DP
    !###  Type: CHARACTER(LEN=MAXSTRLEN)
    !###  Description:
    !###    Converts a double precision number to its equivalent character string representation as determined by the supplied
    !###    format string. Note If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to
    !###    maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the
    !###    supplied number in either F8.? or E8.? format depending on its magnitude.
    !###  Parent-function: NUMBER_TO_CHARACTER
    
    !Argument variables
    REAL(DP), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_DP
    !Local variables
    INTEGER(INTG) :: ASTERISK_POS,i0,i1,LENGTH
    CHARACTER(LEN=2) :: CI0,CI1
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT
    
    CALL ENTERS("NUMBER_TO_CHARACTER_DP",ERR,ERROR,*999)

    ASTERISK_POS=INDEX(FORMAT,"*")
    LENGTH=LEN_TRIM(FORMAT)
    IF(ASTERISK_POS==1.AND.LENGTH==1) THEN !Free format
      WRITE(NUMBER_TO_CHARACTER_DP,*,ERR=999) NUMBER      
    ELSE IF(ASTERISK_POS>0) THEN !Adjustable format
      CI0=FORMAT(ASTERISK_POS+1:LEN_TRIM(FORMAT))
      READ(CI0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(NUMBER>=0.0_DP) THEN
          IF((NUMBER<10.0_DP**(i0-1)).AND.(NUMBER>=0.1_DP**(MIN(i0-2,5)))) THEN
            IF(NUMBER>1.0_DP) THEN
              i1=i0-2-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ELSE
          IF((-NUMBER<10.0_DP**(i0-2)).AND.(-NUMBER>=0.01_DP**(MIN(i0-2,5)))) THEN
            IF(-NUMBER>=1.0_DP) THEN
              i1=i0-3-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,'(I2)') i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid format",ERR,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_FORMAT='('//FORMAT(1:LEN_TRIM(FORMAT))//')'
      WRITE(NUMBER_TO_CHARACTER_DP,LOCAL_FORMAT,ERR=999) NUMBER
    ENDIF

    !Add an extra zero if required
    IF(NUMBER_TO_CHARACTER_DP(LEN_TRIM(NUMBER_TO_CHARACTER_DP):LEN_TRIM(NUMBER_TO_CHARACTER_DP))==".") &
      & NUMBER_TO_CHARACTER_DP=NUMBER_TO_CHARACTER_DP(1:LEN_TRIM(NUMBER_TO_CHARACTER_DP))//"0"
    !Trim leading blanks
    NUMBER_TO_CHARACTER_DP=ADJUSTL(NUMBER_TO_CHARACTER_DP)

    CALL EXITS("NUMBER_TO_CHARACTER_DP")
    RETURN
999 CALL FLAG_ERROR("Error converting double precision number to a character string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_CHARACTER_DP",ERR,ERROR)
    CALL EXITS("NUMBER_TO_CHARACTER_DP")
    RETURN    
  END FUNCTION NUMBER_TO_CHARACTER_DP
  
  !
  !================================================================================================================================
  !
 
  !#### Generic-Function: NUMBER_TO_VSTRING
  !###  Type: TYPE(VARYING_STRING)
  !###  Description:
  !###    Converts a number to its equivalent varying string representation.
  !###  Child-functions: NUMBER_TO_VSTRING_INTG,NUMBER_TO_VSTRING_LINTG,NUMBER_TO_VSTRING_SP,NUMBER_TO_VSTRING_DP

  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_VSTRING_INTG(NUMBER,FORMAT,ERR,ERROR)
  
    !#### Function: NUMBER_TO_VSTRING_INTG
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Converts an integer number to its equivalent varying string representation as determined by the supplied format. The
    !###    format is of the form of a standard Fortran integer format e.g. "I3".
    !###  Parent-function: NUMBER_TO_VSTRING
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_INTG
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT,LOCAL_STRING
    
!!TODO: put back enters,exits.

!    CALL ENTERS("NUMBER_TO_VSTRING_INTG",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING
    
    NUMBER_TO_VSTRING_INTG=""
    
    IF(FORMAT(1:1)=="*") THEN
      LOCAL_FORMAT="(I12)"
    ELSE
      LOCAL_FORMAT="("//FORMAT(1:LEN_TRIM(FORMAT))//")"
    ENDIF
    WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER

    !Trim leading blanks
    NUMBER_TO_VSTRING_INTG=ADJUSTL(LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING)))

!    CALL EXITS("NUMBER_TO_VSTRING_INTG")
    RETURN
999 CALL FLAG_ERROR("Error converting an integer to a varying string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_VSTRING_INTG",ERR,ERROR)
!    CALL EXITS("NUMBER_TO_VSTRING_INTG")
    RETURN    
  END FUNCTION NUMBER_TO_VSTRING_INTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_VSTRING_LINTG(NUMBER,FORMAT,ERR,ERROR)
  
    !#### Function: NUMBER_TO_VSTRING_LINTG
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Converts a long integer number to its equivalent varying string representation as determined by the supplied format.
    !###    The format is of the form of a standard Fortran integer format e.g., "I3".
    !###  Parent-function: NUMBER_TO_VSTRING
    
    !Argument variables
    INTEGER(LINTG), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_LINTG
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT,LOCAL_STRING
    
!!TODO: put back enters,exits.

!    CALL ENTERS("NUMBER_TO_VSTRING_LINTG",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING
    
    NUMBER_TO_VSTRING_LINTG=""
    
    IF(FORMAT(1:1)=="*") THEN
      LOCAL_FORMAT="(I18)"
    ELSE
      LOCAL_FORMAT="("//FORMAT(1:LEN_TRIM(FORMAT))//")"
    ENDIF
    WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER

    !Trim leading blanks
    NUMBER_TO_VSTRING_LINTG=ADJUSTL(LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING)))

!    CALL EXITS("NUMBER_TO_VSTRING_LINTG")
    RETURN
999 CALL FLAG_ERROR("Error converting a long integer to a varying string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_VSTRING_LINTG",ERR,ERROR)
!    CALL EXITS("NUMBER_TO_VSTRING_LINTG")
    RETURN    
  END FUNCTION NUMBER_TO_VSTRING_LINTG
  
  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_VSTRING_SP(NUMBER, FORMAT, ERR, ERROR)
  
    !#### Function: NUMBER_TO_VSTRING_SP
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Converts a single precision number to its equivalent varying string representation as determined by the supplied
    !###    format string. Note If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to
    !###    maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the
    !###    supplied number in either F8.? or E8.? format depending on its magnitude.
    !###  Parent-function: NUMBER_TO_VSTRING
    
    !Argument variables
    REAL(SP), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_SP
    !Local variables
    INTEGER(INTG) :: ASTERISK_POS,i0,i1,LENGTH
    CHARACTER(LEN=MAXSTRLEN) :: CI0,CI1
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT,LOCAL_STRING
    
    CALL ENTERS("NUMBER_TO_VSTRING_SP",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING
    
    NUMBER_TO_VSTRING_SP=""    

    ASTERISK_POS=INDEX(FORMAT,"*")
    LENGTH=LEN_TRIM(FORMAT)
    IF(ASTERISK_POS==1.AND.LENGTH==1) THEN !Free format
      WRITE(LOCAL_STRING,*,ERR=999) NUMBER      
    ELSE IF(ASTERISK_POS>0) THEN !Adjustable format
      CI0=FORMAT(ASTERISK_POS+1:LEN_TRIM(FORMAT))
      READ(CI0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(NUMBER>=0.0_SP) THEN
          IF((NUMBER<10.0_SP**(i0-1)).AND.(NUMBER>=0.1_SP**(MIN(i0-2,5)))) THEN
            IF(NUMBER>1.0_SP) THEN
              i1=i0-2-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ELSE
          IF((-NUMBER<10.0_SP**(i0-2)).AND.(-NUMBER>=0.01_SP**(MIN(i0-2,5)))) THEN
            IF(-NUMBER>=1.0_SP) THEN
              i1=i0-3-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,'(I2)') i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid format",ERR,ERROR,*999)
        GOTO 999
      ENDIF
    ELSE
      LOCAL_FORMAT='('//FORMAT(1:LEN_TRIM(FORMAT))//')'
      WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
    ENDIF

    !Add an extra zero if required
    IF(LOCAL_STRING(LEN_TRIM(LOCAL_STRING):LEN_TRIM(LOCAL_STRING))==".") LOCAL_STRING=LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING))//"0"
    
    !Trim leading blanks
    NUMBER_TO_VSTRING_SP=ADJUSTL(LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING)))

    CALL EXITS("NUMBER_TO_VSTRING_SP")
    RETURN
999 CALL FLAG_ERROR("Error converting a single precision number to a varying string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_VSTRING_SP",ERR,ERROR)
    CALL EXITS("NUMBER_TO_VSTRING_SP")
    RETURN    
  END FUNCTION NUMBER_TO_VSTRING_SP
  
  !
  !================================================================================================================================
  !
  
  FUNCTION NUMBER_TO_VSTRING_DP(NUMBER, FORMAT, ERR, ERROR)
  
    !#### Function: NUMBER_TO_VSTRING_DP
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Converts a double precision number to its equivalent varying string representation as determined by the supplied
    !###    format string. Note If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to
    !###    maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the
    !###    supplied number in either F8.? or E8.? format depending on its magnitude.
    !###  Parent-function: NUMBER_TO_VSTRING
    
    !Argument variables
    REAL(DP), INTENT(IN) :: NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_DP
    !Local variables
    INTEGER(INTG) :: ASTERISK_POS,i0,i1,LENGTH
    CHARACTER(LEN=2) :: CI0,CI1
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_FORMAT,LOCAL_STRING
    
!!TODO: put back enters,exits.

!    CALL ENTERS("NUMBER_TO_VSTRING_DP",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING
    
    NUMBER_TO_VSTRING_DP=""    

    ASTERISK_POS=INDEX(FORMAT,"*")
    LENGTH=LEN_TRIM(FORMAT)
    IF(ASTERISK_POS==1.AND.LENGTH==1) THEN !Free format
      WRITE(LOCAL_STRING,*,ERR=999) NUMBER      
    ELSE IF(ASTERISK_POS>0) THEN !Adjustable format
      CI0=FORMAT(ASTERISK_POS+1:LEN_TRIM(FORMAT))
      READ(CI0,'(BN,I2)') i0
      IF(i0<=MAXSTRLEN) THEN
        IF(NUMBER>=0.0_DP) THEN
          IF((NUMBER<10.0_DP**(i0-1)).AND.(NUMBER>=0.1_DP**(MIN(i0-2,5)))) THEN
            IF(NUMBER>1.0_DP) THEN
              i1=i0-2-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ELSE
          IF((-NUMBER<10.0_DP**(i0-2)).AND.(-NUMBER>=0.01_DP**(MIN(i0-2,5)))) THEN
            IF(-NUMBER>=1.0_DP) THEN
              i1=i0-3-LOG10(NUMBER)
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,'(I2)') i1
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ELSE
              LOCAL_FORMAT="(I2)"
              WRITE(CI1,LOCAL_FORMAT) i0-2
              LOCAL_FORMAT="(F"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
              WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
            ENDIF
          ELSE
            LOCAL_FORMAT="(I2)"
            WRITE(CI1,LOCAL_FORMAT) i0-6
            LOCAL_FORMAT="(E"//CI0(1:LEN_TRIM(CI0))//"."//CI1(1:LEN_TRIM(CI1))//")"
            WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Invalid format",ERR,ERROR,*999)
      ENDIF
    ELSE
      LOCAL_FORMAT='('//FORMAT(1:LEN_TRIM(FORMAT))//')'
      WRITE(LOCAL_STRING,LOCAL_FORMAT,ERR=999) NUMBER
    ENDIF

    !Add an extra zero if required
    IF(LOCAL_STRING(LEN_TRIM(LOCAL_STRING):LEN_TRIM(LOCAL_STRING))==".") LOCAL_STRING=LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING))//"0"

    !!Do you really want to do this???
    !Trim leading blanks
    !NUMBER_TO_VSTRING_DP=ADJUSTL(LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING)))
    NUMBER_TO_VSTRING_DP=LOCAL_STRING(1:LEN_TRIM(LOCAL_STRING))

!    CALL EXITS("NUMBER_TO_VSTRING_DP")
    RETURN
999 CALL FLAG_ERROR("Error converting double precision number to a varying string",ERR,ERROR,*998)
998 CALL ERRORS("NUMBER_TO_VSTRING_DP",ERR,ERROR)
!    CALL EXITS("NUMBER_TO_VSTRING_DP")
    RETURN    
  END FUNCTION NUMBER_TO_VSTRING_DP
  
  !
  !================================================================================================================================
  !

  !#### Generic-Function: STRING_TO_DOUBLE
  !###  Type: REAL(DP)
  !###  Description:
  !###    Converts a string representation of a number to a double precision number.
  !###  Child-functions: STRING_TO_DOUBLE_C,STRING_TO_DOUBLE_VS
  
  !
  !================================================================================================================================
  !

  !> Converts a character string representation of a number to a double precision number.
  FUNCTION STRING_TO_DOUBLE_C(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_DOUBLE_C
    !###  Type: REAL(DP)      
    !###  Description:
    !###    Converts a character string representation of a number to a double precision number.
    !###  Parent-function: STRING_TO_DOUBLE
    
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: STRING_TO_DOUBLE_C
    !Local variables
    
    CALL ENTERS("STRING_TO_DOUBLE_C",ERR,ERROR,*999)

    READ(STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_DOUBLE_C

    CALL EXITS("STRING_TO_DOUBLE_C")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//STRING(1:LEN_TRIM(STRING))//"' to a double real",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_DOUBLE_C",ERR,ERROR)
    CALL EXITS("STRING_TO_DOUBLE_C")
    RETURN    
  END FUNCTION STRING_TO_DOUBLE_C
  
  !
  !================================================================================================================================
  !
  
  FUNCTION STRING_TO_DOUBLE_VS(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_DOUBLE_vs
    !###  Type: REAL(DP)
    !###  Description:
    !###    Converts a varying string representation of a number to a double precision number.
    !###  Parent-function: STRING_TO_DOUBLE
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(DP) :: STRING_TO_DOUBLE_VS
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_STRING
    
    CALL ENTERS("STRING_TO_DOUBLE_VS",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING

    LOCAL_STRING=CHAR(STRING)
    READ(LOCAL_STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_DOUBLE_VS

    CALL EXITS("STRING_TO_DOUBLE_VS")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//CHAR(STRING)//"' to a double real",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_DOUBLE_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_DOUBLE_VS")
    RETURN    
  END FUNCTION STRING_TO_DOUBLE_VS
  
  !
  !================================================================================================================================
  !

  !#### Generic-Function: STRING_TO_INTEGER
  !###  Type: INTEGER(INTG)
  !###  Description:
  !###    Converts a string representation of a number to an integer.
  !###  Child-functions: STRING_TO_INTEGER_C,STRING_TO_INTEGER_VS
  
  !
  !================================================================================================================================
  !

  FUNCTION STRING_TO_INTEGER_C(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_INTEGER_C
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Converts a character string representation of a number to an integer.
    !###  Parent-function: STRING_TO_INTEGER
    
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    INTEGER(INTG) :: STRING_TO_INTEGER_C
    !Local variables
    
    CALL ENTERS("STRING_TO_INTEGER_C",ERR,ERROR,*999)

    READ(STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_INTEGER_C

    CALL EXITS("STRING_TO_INTEGER_C")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//STRING(1:LEN_TRIM(STRING))//"' to an integer",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_INTEGER_C",ERR,ERROR)
    CALL EXITS("STRING_TO_INTEGER_C")
    RETURN 
  END FUNCTION STRING_TO_INTEGER_C
  
  !
  !================================================================================================================================
  !

  FUNCTION STRING_TO_INTEGER_VS(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_INTEGER_VS
    !###  Type: INTEGER(INTG)
    !###  Description:
    !###    Converts a varying string representation of a number to an integer.
    !###  Parent-function: STRING_TO_INTEGER
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    INTEGER(INTG) :: STRING_TO_INTEGER_VS
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_STRING
    
    CALL ENTERS("STRING_TO_INTEGER_VS",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING

    LOCAL_STRING=CHAR(STRING)
    READ(LOCAL_STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_INTEGER_VS

    CALL EXITS("STRING_TO_INTEGER_VS")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//CHAR(STRING)//"' to an integer",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_INTEGER_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_INTEGER_VS")
    RETURN 
  END FUNCTION STRING_TO_INTEGER_VS
  
  !
  !================================================================================================================================
  !

  !#### Generic-Function: STRING_TO_LONG_INTEGER
  !###  Type: INTEGER(INTG)
  !###  Description:
  !###    Converts a string representation of a number to a long integer.
  !###  Child-functions: STRING_TO_LONG_INTEGER_C,STRING_TO_LONG_INTEGER_VS
  
  !
  !================================================================================================================================
  !

  FUNCTION STRING_TO_LONG_INTEGER_C(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_LONG_INTEGER_C
    !###  Type: INTEGER(LINTG)
    !###  Description:
    !###    Converts a character string representation of a number to a long integer.
    !###  Parent-function: STRING_TO_LONG_INTEGER
    
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    INTEGER(LINTG) :: STRING_TO_LONG_INTEGER_C
    !Local variables
    
    CALL ENTERS("STRING_TO_LONG_INTEGER_C",ERR,ERROR,*999)

    READ(STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_LONG_INTEGER_C

    CALL EXITS("STRING_TO_LONG_INTEGER_C")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//STRING(1:LEN_TRIM(STRING))//"' to a long integer",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_LONG_INTEGER_C",ERR,ERROR)
    CALL EXITS("STRING_TO_LONG_INTEGER_C")
    RETURN 
  END FUNCTION STRING_TO_LONG_INTEGER_C
  
  !
  !================================================================================================================================
  !

  FUNCTION STRING_TO_LONG_INTEGER_VS(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_LONG_INTEGER_VS
    !###  Type: INTEGER(LINTG)
    !###  Description:
    !###    Converts a varying string representation of a number to a long integer.
    !###  Parent-function: STRING_TO_INTEGER
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    INTEGER(LINTG) :: STRING_TO_LONG_INTEGER_VS
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_STRING
    
    CALL ENTERS("STRING_TO_LONG_INTEGER_VS",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING

    LOCAL_STRING=CHAR(STRING)
    READ(LOCAL_STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_LONG_INTEGER_VS

    CALL EXITS("STRING_TO_LONG_INTEGER_VS")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//CHAR(STRING)//"' to a long integer",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_LONG_INTEGER_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_LONG_INTEGER_VS")
    RETURN 
  END FUNCTION STRING_TO_LONG_INTEGER_VS
  
  !
  !================================================================================================================================
  !

  !#### Generic-Function: STRING_TO_LOGICAL
  !###  Type: LOGICAL
  !###  Description:
  !###    Converts a string representation of a boolean value (TRUE or FALSE) to a logical.
  !###  Child-functions: STRING_TO_LOGICAL_C,STRING_TO_LOGICAL_VS
  
  !
  !================================================================================================================================
  !
  
  FUNCTION STRING_TO_LOGICAL_C(STRING,ERR,ERROR)
  
    !#### Function: STRING_TO_LOGICAL_C
    !###  Type: LOGICAL
    !###  Description:
    !###    Converts a character string representation of a boolean (TRUE or FALSE) to a logical.
    !###  Parent-function: STRING_TO_LOGICAL
    
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    LOGICAL :: STRING_TO_LOGICAL_C
    !Local variables
    
    CALL ENTERS("STRING_TO_LOGICAL_C",ERR,ERROR,*999)

    READ(STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_LOGICAL_C

    CALL EXITS("STRING_TO_LOGICAL_C")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//STRING(1:LEN_TRIM(STRING))//"' to a logical",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_LOGICAL_C",ERR,ERROR)
    CALL EXITS("STRING_TO_LOGICAL_C")
    RETURN    
  END FUNCTION STRING_TO_LOGICAL_C
  
  !
  !================================================================================================================================
  !
  
  FUNCTION STRING_TO_LOGICAL_VS(STRING,ERR,ERROR)
  
    !#### Function: STRING_TO_LOGICAL_VS
    !###  Type: LOGICAL
    !###  Description:
    !###    Converts a varying string representation of a boolean (TRUE or FALSE) to a logical.
    !###  Parent-function: STRING_TO_LOGICAL
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    LOGICAL :: STRING_TO_LOGICAL_VS
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_STRING
    
    CALL ENTERS("STRING_TO_LOGICAL_VS",ERR,ERROR,*999)

    LOCAL_STRING=CHAR(STRING)
    READ(LOCAL_STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_LOGICAL_VS

    CALL EXITS("STRING_TO_LOGICAL_VS")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//CHAR(STRING)//"' to a logical",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_LOGICAL_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_LOGICAL_VS")
    RETURN    
  END FUNCTION STRING_TO_LOGICAL_VS
  
  !
  !================================================================================================================================
  !

  !#### Generic-Function: STRING_TO_SINGLE
  !###  Type: REAL(SP)
  !###  Description:
  !###    Converts a string representation of a number to a single precision number.
  !###  Child-functions: STRING_TO_SINGLE_C,STRING_TO_SINGLE_VS
  
  !
  !================================================================================================================================
  !

  FUNCTION STRING_TO_SINGLE_C(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_SINGLE_C
    !###  Type: REAL(SP)
    !###  Description:
    !###    Converts a character string representation of a number to a single precision number.
    !###  Parent-function: STRING_TO_SINGLE
    
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(SP) :: STRING_TO_SINGLE_C
    !Local variables
    
    CALL ENTERS("STRING_TO_SINGLE_C",ERR,ERROR,*999)

    READ(STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_SINGLE_C
    
    CALL EXITS("STRING_TO_SINGLE_C")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//STRING(1:LEN_TRIM(STRING))//"' to a single real",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_SINGLE_C",ERR,ERROR)
    CALL EXITS("STRING_TO_SINGLE_C")
    RETURN    
  END FUNCTION STRING_TO_SINGLE_C
    
  !
  !================================================================================================================================
  !

  FUNCTION STRING_TO_SINGLE_VS(STRING, ERR, ERROR)
  
    !#### Function: STRING_TO_SINGLE_VS
    !###  Type: REAL(SP)
    !###  Description:
    !###    Converts a varying string representation of a number to a single precision number.
    !###  Parent-function: STRING_TO_SINGLE
    
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function variable
    REAL(SP) :: STRING_TO_SINGLE_VS
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_STRING
    
    CALL ENTERS("STRING_TO_SINGLE_VS",ERR,ERROR,*999)

!!TODO: remove dependance on LOCAL_STRING

    LOCAL_STRING=CHAR(STRING)
    READ(LOCAL_STRING,*,IOSTAT=ERR,ERR=999) STRING_TO_SINGLE_VS
    
    CALL EXITS("STRING_TO_SINGLE_VS")
    RETURN
999 CALL FLAG_ERROR("Cannot convert '"//CHAR(STRING)//"' to a single real",ERR,ERROR,*998)
998 CALL ERRORS("STRING_TO_SINGLE_VS",ERR,ERROR)
    CALL EXITS("STRING_TO_SINGLE_VS")
    RETURN    
  END FUNCTION STRING_TO_SINGLE_VS
    
  !
  !================================================================================================================================
  !

  !#### Generic-Function: CHARACTER_TO_LOWERCASE
  !###  Type: CHARACTER(LEN=LEN(STRING))
  !###  Description:
  !###    Returns a character string which is the lowercase equivalent of the supplied string.
  !###  Child-functions: CHARACTER_TO_LOWERCASE_C,CHARACTER_TO_LOWERCASE_VS

  !
  !================================================================================================================================
  !
    
  FUNCTION CHARACTER_TO_LOWERCASE_C(STRING)

    !#### Function: CHARACTER_TO_LOWERCASE_C
    !###  Type: CHARACTER(LEN=LEN(STRING))
    !###  Description:
    !###    Returns a character string that is the lowercase equivalent of the supplied character string.
    !###  Parent-function: CHARACTER_TO_LOWERCASE

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_LOWERCASE_C
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: i

    CHARACTER_TO_LOWERCASE_C=STRING
    DO i=1,LEN(STRING)
      IF(IS_UPPERCASE(STRING(i:i))) THEN
        CHARACTER_TO_LOWERCASE_C(i:i)=CHAR(ICHAR(STRING(i:i))+OFFSET)
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION CHARACTER_TO_LOWERCASE_C

  !
  !================================================================================================================================
  !
    
  FUNCTION CHARACTER_TO_LOWERCASE_VS(STRING)

    !#### Function: CHARACTER_TO_LOWERCASE_VS
    !###  Type: CHARACTER(LEN=LEN(STRING))
    !###  Description:
    !###    Returns a character string that is the lowercase equivalent of the supplied varying string.
    !###  Parent-function: CHARACTER_TO_LOWERCASE

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_LOWERCASE_VS
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: i

    CHARACTER_TO_LOWERCASE_VS=CHAR(STRING)
    DO i=1,LEN(STRING)
      IF(IS_UPPERCASE(CHAR(EXTRACT(STRING,i,i)))) THEN
        CHARACTER_TO_LOWERCASE_VS(i:i)=CHAR(ICHAR(EXTRACT(STRING,i,i))+OFFSET)
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION CHARACTER_TO_LOWERCASE_VS

  !
  !================================================================================================================================
  !

  !#### Generic-Function: VSTRING_TO_LOWERCASE
  !###  Type: TYPE(VARYING_STRING)
  !###  Description:
  !###    Returns a varying string which is the lowercase equivalent of the supplied string.
  !###  Child-functions: VSTRING_TO_LOWERCASE_C,VSTRING_TO_LOWERCASE_VS

  !
  !================================================================================================================================
  !
    
  FUNCTION VSTRING_TO_LOWERCASE_C(STRING)

    !#### Function: VSTRING_TO_LOWERCASE_C
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Returns a varying string that is the lowercase equivalent of the supplied character string.
    !###  Parent-function: VSTRING_TO_LOWERCASE

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_LOWERCASE_C
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: i

    VSTRING_TO_LOWERCASE_C=STRING
    DO i=1,LEN(STRING)
      IF(IS_UPPERCASE(STRING(i:i))) THEN
        VSTRING_TO_LOWERCASE_C=INSERT(VSTRING_TO_LOWERCASE_C,i,CHAR(ICHAR(STRING(i:i))+OFFSET))
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION VSTRING_TO_LOWERCASE_C

  !
  !================================================================================================================================
  !
    
  FUNCTION VSTRING_TO_LOWERCASE_VS(STRING)

    !#### Function: VSTRING_TO_LOWERCASE_VS
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Returns a varying string that is the lowercase equivalent of the supplied varying string.
    !###  Parent-function: VSTRING_TO_LOWERCASE

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_LOWERCASE_VS
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("a")-ICHAR("A"))
    INTEGER(INTG) :: i

    VSTRING_TO_LOWERCASE_VS=STRING
    DO i=1,LEN(STRING)
      IF(IS_UPPERCASE(CHAR(EXTRACT(STRING,i,i)))) THEN
        VSTRING_TO_LOWERCASE_VS=INSERT(VSTRING_TO_LOWERCASE_VS,i,CHAR(ICHAR(EXTRACT(STRING,i,i))+OFFSET))
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION VSTRING_TO_LOWERCASE_VS

  !
  !================================================================================================================================
  !

  !#### Generic-Function: CHARACTER_TO_UPPERCASE
  !###  Type: CHARACTER(LEN=LEN(STRING))
  !###  Description:
  !###    Returns a character string which is the uppercase equivalent of the supplied string.
  !###  Child-functions: CHARACTER_TO_UPPERCASE_C,CHARACTER_TO_UPPERCASE_VS

  !
  !================================================================================================================================
  !
  
  FUNCTION CHARACTER_TO_UPPERCASE_C(STRING)

    !#### Function: CHARACTER_TO_UPPERCASE_C
    !###  Type: CHARACTER(LEN=LEN(STRING))
    !###  Description:
    !###    Returns a character string which is uppercase equivalent of the supplied character string.
    !###  Parent-function: CHARACTER_TO_UPPERCASE

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_UPPERCASE_C
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: i

    CHARACTER_TO_UPPERCASE_C=STRING
    DO i=1,LEN(STRING)
      IF(IS_LOWERCASE(STRING(i:i))) THEN
        CHARACTER_TO_UPPERCASE_C(i:i)=CHAR(ICHAR(STRING(i:i))+OFFSET)
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION CHARACTER_TO_UPPERCASE_C

  !
  !================================================================================================================================
  !
  
  FUNCTION CHARACTER_TO_UPPERCASE_VS(STRING)

    !#### Function: CHARACTER_TO_UPPERCASE_VS
    !###  Type: CHARACTER(LEN=LEN(STRING))
    !###  Description:
    !###    Returns a character string which is uppercase equivalent of the supplied varying string.
    !###  Parent-function: CHARACTER_TO_UPPERCASE

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_UPPERCASE_VS
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: i

    CHARACTER_TO_UPPERCASE_VS=CHAR(STRING)
    DO i=1,LEN(STRING)
      IF(IS_LOWERCASE(CHAR(EXTRACT(STRING,i,i)))) THEN
        CHARACTER_TO_UPPERCASE_VS(i:i)=CHAR(ICHAR(EXTRACT(STRING,i,i))+OFFSET)
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION CHARACTER_TO_UPPERCASE_VS

  !
  !================================================================================================================================
  !

  !#### Generic-Function: VSTRING_TO_UPPERCASE
  !###  Type: TYPE(VARYING_STRING)
  !###  Description:
  !###    Returns a varying string which is the uppercase equivalent of the supplied string.
  !###  Child-functions: VSTRING_TO_UPPERCASE_C,VSTRING_TO_UPPERCASE_VS

  !
  !================================================================================================================================
  !
  
  FUNCTION VSTRING_TO_UPPERCASE_C(STRING)

    !#### Function: VSTRING_TO_UPPERCASE_C
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Returns a varying string which is uppercase equivalent of the supplied character string.
    !###  Parent-function: VSTRING_TO_UPPERCASE

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_UPPERCASE_C
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: i

    VSTRING_TO_UPPERCASE_C=STRING
    DO i=1,LEN(STRING)
      IF(IS_LOWERCASE(STRING(i:i))) THEN
        VSTRING_TO_UPPERCASE_C=INSERT(VSTRING_TO_UPPERCASE_C,i,CHAR(ICHAR(STRING(i:i))+OFFSET))
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION VSTRING_TO_UPPERCASE_C

  !
  !================================================================================================================================
  !
  
  FUNCTION VSTRING_TO_UPPERCASE_VS(STRING)

    !#### Function: VSTRING_TO_UPPERCASE_VS
    !###  Type: TYPE(VARYING_STRING)
    !###  Description:
    !###    Returns a varying string which is uppercase equivalent of the supplied varying string.
    !###  Parent-function: VSTRING_TO_UPPERCASE

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_UPPERCASE_VS
    !Local Variables
    INTEGER(INTG), PARAMETER :: OFFSET=(ICHAR("A")-ICHAR("a"))
    INTEGER(INTG) :: i

    VSTRING_TO_UPPERCASE_VS=STRING
    DO i=1,LEN(STRING)
      IF(IS_LOWERCASE(CHAR(EXTRACT(STRING,i,i)))) THEN
        VSTRING_TO_UPPERCASE_VS=INSERT(VSTRING_TO_UPPERCASE_VS,i,CHAR(ICHAR(EXTRACT(STRING,i,i))+OFFSET))
      ENDIF
    ENDDO !i

    RETURN
  END FUNCTION VSTRING_TO_UPPERCASE_VS

  !
  !================================================================================================================================
  !
  
END MODULE STRINGS
