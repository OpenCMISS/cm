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

!>This module contains all string manipulation and transformation routines.
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

  !>Returns .TRUE. if a supplied string is a valid abbreviation of a second supplied string.
  INTERFACE IS_ABBREVIATION
    MODULE PROCEDURE IS_ABBREVIATION_C_C
    MODULE PROCEDURE IS_ABBREVIATION_C_VS
    MODULE PROCEDURE IS_ABBREVIATION_VS_C
    MODULE PROCEDURE IS_ABBREVIATION_VS_VS
  END INTERFACE !IS_ABBREVIATION

  !>Converts a list to its equivalent character string representation.
  INTERFACE LIST_TO_CHARACTER
    MODULE PROCEDURE LIST_TO_CHARACTER_C
    MODULE PROCEDURE LIST_TO_CHARACTER_INTG
    MODULE PROCEDURE LIST_TO_CHARACTER_LINTG
    MODULE PROCEDURE LIST_TO_CHARACTER_L
    MODULE PROCEDURE LIST_TO_CHARACTER_SP
    MODULE PROCEDURE LIST_TO_CHARACTER_DP
  END INTERFACE !LIST_TO_CHARACTER

  !>Converts a number to its equivalent character string representation.
  INTERFACE NUMBER_TO_CHARACTER
    MODULE PROCEDURE NUMBER_TO_CHARACTER_INTG
    MODULE PROCEDURE NUMBER_TO_CHARACTER_LINTG
    MODULE PROCEDURE NUMBER_TO_CHARACTER_SP
    MODULE PROCEDURE NUMBER_TO_CHARACTER_DP
  END INTERFACE !NUMBER_TO_CHARACTER

  !>Converts a number to its equivalent varying string representation.
  INTERFACE NUMBER_TO_VSTRING
    MODULE PROCEDURE NUMBER_TO_VSTRING_INTG
    MODULE PROCEDURE NUMBER_TO_VSTRING_LINTG
    MODULE PROCEDURE NUMBER_TO_VSTRING_SP
    MODULE PROCEDURE NUMBER_TO_VSTRING_DP
  END INTERFACE !NUMBER_TO_VSTRING

  !>Converts a string representation of a number to a double precision number.
  INTERFACE STRING_TO_DOUBLE
    MODULE PROCEDURE STRING_TO_DOUBLE_C
    MODULE PROCEDURE STRING_TO_DOUBLE_VS
  END INTERFACE !STRING_TO_DOUBLE

  !>Converts a string representation of a number to an integer.
  INTERFACE STRING_TO_INTEGER
    MODULE PROCEDURE STRING_TO_INTEGER_C
    MODULE PROCEDURE STRING_TO_INTEGER_VS
  END INTERFACE !STRING_TO_INTEGER

  !>Converts a string representation of a number to a long integer.
  INTERFACE STRING_TO_LONG_INTEGER
    MODULE PROCEDURE STRING_TO_LONG_INTEGER_C
    MODULE PROCEDURE STRING_TO_LONG_INTEGER_VS
  END INTERFACE !STRING_TO_LONG_INTEGER

  !>Converts a string representation of a boolean value (TRUE or FALSE) to a logical.
  INTERFACE STRING_TO_LOGICAL
    MODULE PROCEDURE STRING_TO_LOGICAL_C
    MODULE PROCEDURE STRING_TO_LOGICAL_VS
  END INTERFACE !STRING_TO_LOGICAL

  !>Converts a string representation of a number to a single precision number.
  INTERFACE STRING_TO_SINGLE
    MODULE PROCEDURE STRING_TO_SINGLE_C
    MODULE PROCEDURE STRING_TO_SINGLE_VS
  END INTERFACE !STRING_TO_SINGLE

  !>Returns a character string which is the lowercase equivalent of the supplied string.
  INTERFACE CHARACTER_TO_LOWERCASE
    MODULE PROCEDURE CHARACTER_TO_LOWERCASE_C
    MODULE PROCEDURE CHARACTER_TO_LOWERCASE_VS
  END INTERFACE !CHARACTER_TO_LOWERCASE

  !>Returns a varying string which is the lowercase equivalent of the supplied string.
  INTERFACE VSTRING_TO_LOWERCASE
    MODULE PROCEDURE VSTRING_TO_LOWERCASE_C
    MODULE PROCEDURE VSTRING_TO_LOWERCASE_VS
  END INTERFACE !VSTRING_TO_LOWERCASE

  !>Returns a character string which is the uppercase equivalent of the supplied string.
  INTERFACE CHARACTER_TO_UPPERCASE
    MODULE PROCEDURE CHARACTER_TO_UPPERCASE_C
    MODULE PROCEDURE CHARACTER_TO_UPPERCASE_VS
  END INTERFACE !CHARACTER_TO_UPPERCASE

  !>Returns a varying string which is the uppercase equivalent of the supplied string.
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

  !>IS_ABBREVIATION returns .TRUE. if the character string SHORT is an abbreviation of the character string LONG. SHORT must be at least MIN_NUM_CHARACTERS long.
  FUNCTION IS_ABBREVIATION_C_C(SHORT,LONG,MIN_NUM_CHARACTERS)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: SHORT !<The short form of the string
    CHARACTER(LEN=*), INTENT(IN) :: LONG !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IS_ABBREVIATION_C_C !<On exit, .TRUE. if the short string is an abbreviation
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

  !>IS_ABBREVIATION returns .TRUE. if the character string SHORT is an abbreviation of the varying string LONG. SHORT must be at least MIN_NUM_CHARACTERS long.
  FUNCTION IS_ABBREVIATION_C_VS(SHORT,LONG,MIN_NUM_CHARACTERS)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: SHORT !<The short form of the string
    TYPE(VARYING_STRING), INTENT(IN) :: LONG !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IS_ABBREVIATION_C_VS !<On exit, .TRUE. if the short string is an abbreviation
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

  !>IS_ABBREVIATION returns .TRUE. if the varying string SHORT is an abbreviation of the character string LONG. SHORT must be at least MIN_NUM_CHARACTERS long.
  FUNCTION IS_ABBREVIATION_VS_C(SHORT,LONG,MIN_NUM_CHARACTERS)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: SHORT !<The short form of the string
    CHARACTER(LEN=*), INTENT(IN) :: LONG !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IS_ABBREVIATION_VS_C !<On exit, .TRUE. if the short string is an abbreviation
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

  !>IS_ABBREVIATION returns .TRUE. if the varying string SHORT is an abbreviation of the varying string LONG. SHORT must be at least MIN_NUM_CHARACTERS long.
  FUNCTION IS_ABBREVIATION_VS_VS(SHORT,LONG,MIN_NUM_CHARACTERS)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: SHORT !<The short form of the string
    TYPE(VARYING_STRING), INTENT(IN) :: LONG !<The long form of the string
    INTEGER(INTG), INTENT(IN) :: MIN_NUM_CHARACTERS !<The minimum number of characters to match
    !Function variable
    LOGICAL :: IS_ABBREVIATION_VS_VS !<On exit, .TRUE. if the short string is an abbreviation
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

  !>IS_DIGIT returns .TRUE. if the character CHARAC is a digit character (i.e. 0..9)
  FUNCTION IS_DIGIT(CHARAC)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARAC !<The character to test if it is a digit
    !Function variable
    LOGICAL :: IS_DIGIT !<On exit, .TRUE. if the character is a digit
    !Local Variables

    IS_DIGIT=(ICHAR(CHARAC)>=ICHAR("0").AND.ICHAR(CHARAC)<=ICHAR("9"))

    RETURN
  END FUNCTION IS_DIGIT

  !
  !================================================================================================================================
  !

  !>IS_LETTER returns .TRUE. if the character CHARAC is a letter character (i.e. A..Z or a..z)
  FUNCTION IS_LETTER(CHARAC)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARAC !<The character to test if it is a letter
    !Function variable
    LOGICAL :: IS_LETTER !<On exit, .TRUE. if the character is a letter
    !Local Variables

    IS_LETTER=((ICHAR(CHARAC)>=ICHAR("A").AND.ICHAR(CHARAC)<=ICHAR("Z")).OR.&
	    & (ICHAR(CHARAC)>=ICHAR("a").AND.ICHAR(CHARAC)<=ICHAR("z")))

    RETURN
  END FUNCTION IS_LETTER

  !
  !================================================================================================================================
  !

  !>Returns .TRUE. if the supplied character is a lowercase character.
  FUNCTION IS_LOWERCASE(CHARC)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARC !<The character to test if it is lowercase
    !Function variable
    LOGICAL :: IS_LOWERCASE !<On exit, .TRUE. if the character is lowercase
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

  !>Returns .TRUE. if the supplied character is an uppercase character.
  FUNCTION IS_UPPERCASE(CHARC)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARC !<The character to test if it is uppercase
    !Function variable
    LOGICAL :: IS_UPPERCASE !<On exit, .TRUE. if the character is uppercase
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

  !>IS_WHITESPACE returns .TRUE. if the character CHARAC is a whitespace character (i.e. space, tabs, etc.)
  FUNCTION IS_WHITESPACE(CHARAC)

    !Argument variables
    CHARACTER(LEN=1), INTENT(IN) :: CHARAC !<The character to test if it is whitespace
    !Function variable
    LOGICAL :: IS_WHITESPACE !<On exit, .TRUE. if the character is whitespace
    !Local Variables
    
    !!WARNING: Assumes ASCII encoding
    IS_WHITESPACE=(CHARAC==CHAR(32).OR.CHARAC==CHAR(9))

    RETURN
  END FUNCTION IS_WHITESPACE

  !
  !================================================================================================================================
  !

  !>Converts a character list to its equivalent character string representation as determined by the supplied format. If present, the optional argument LIST_LENGTHS is used for the lengths of each list elements length otherwise the trimmed length is used. NOTE: The FORMAT is ignored for this child FUNCTION.
  FUNCTION LIST_TO_CHARACTER_C(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR,LIST_LENGTHS)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST !<The number of items in the list
    CHARACTER(LEN=*), INTENT(IN) :: LIST(NUMBER_IN_LIST) !<LIST(i). The i'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use. Ignored for character lists.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    INTEGER(INTG), OPTIONAL, INTENT(IN) :: LIST_LENGTHS(NUMBER_IN_LIST) !<LIST_LENGTHS(i). Optional, The length of the i'th list item.
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_C !<On exit, the character equivalent of the list
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
  
  !>Converts an integer list to its equivalent character string representation as determined by the supplied format. 
  FUNCTION LIST_TO_CHARACTER_INTG(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST !<The number of items in the list
    INTEGER(INTG), INTENT(IN) :: LIST(NUMBER_IN_LIST) !<LIST(i). The i'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_INTG !<On exit, the character equivalent of the list
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

  !>Converts an long integer list to its equivalent character string representation as determined by the supplied format. 
  FUNCTION LIST_TO_CHARACTER_LINTG(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST !<The number of items in the list
    INTEGER(LINTG), INTENT(IN) :: LIST(NUMBER_IN_LIST) !<LIST(i). The i'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_LINTG !<On exit, the character equivalent of the list
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

  !>Converts a logical list to its equivalent character string representation as determined by the supplied format string. The FORMAT is ignored for this child FUNCTION.
  FUNCTION LIST_TO_CHARACTER_L(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST !<The number of items in the list
    LOGICAL, INTENT(IN) :: LIST(NUMBER_IN_LIST) !<LIST(i). The i'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use. Ignored for logical lists.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_L !<On exit, the character equivalent of the list
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

  !>Converts a single precision list to its equivalent character string representation as determined by the supplied format string.
  FUNCTION LIST_TO_CHARACTER_SP(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST !<The number of items in the list
    REAL(SP), INTENT(IN) :: LIST(NUMBER_IN_LIST) !<LIST(i). The i'th item in the list 
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_SP !<On exit, the character equivalent of the list
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

  !>Converts a double precision list to its equivalent character string representation as determined by the supplied format string.
  FUNCTION LIST_TO_CHARACTER_DP(NUMBER_IN_LIST,LIST,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER_IN_LIST !<The number of items in the list
    REAL(DP), INTENT(IN) :: LIST(NUMBER_IN_LIST) !<LIST(i). The i'th item in the list
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use for the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LIST_TO_CHARACTER_DP !<On exit, the character equivalent of the list
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

  !>Converts a logical value to either a "TRUE" or "FALSE" character string.
  FUNCTION LOGICAL_TO_CHARACTER(LOGICALVALUE,ERR,ERROR)
  
    !Argument variables
    LOGICAL, INTENT(IN) :: LOGICALVALUE !<The logical value to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: LOGICAL_TO_CHARACTER !<On exit, the character equivalent value
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

  !>Converts a logical value to either a "TRUE" or "FALSE" varying string.
  FUNCTION LOGICAL_TO_VSTRING(LOGICALVALUE,ERR,ERROR)
    
    !Argument variables
    LOGICAL, INTENT(IN) :: LOGICALVALUE !<The logical value to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: LOGICAL_TO_VSTRING !<On exit, the varying string equivalent value
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

  !>Converts an integer number to its equivalent character string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g. "I3".
  FUNCTION NUMBER_TO_CHARACTER_INTG(NUMBER,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_INTG !<On exit, the character equivalent of the number
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

  !>Converts a long integer number to its equivalent character string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g. "I3".
  FUNCTION NUMBER_TO_CHARACTER_LINTG(NUMBER,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(LINTG), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_LINTG !<On exit, the character equivalent of the number
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

  !>Converts a single precision number to its equivalent character string representation as determined by the supplied format string. NOTE: If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NUMBER_TO_CHARACTER_SP(NUMBER, FORMAT, ERR, ERROR)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_SP !<On exit, the character equivalent of the number
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

  !>Converts a double precision number to its equivalent character string representation as determined by the supplied format string. Note If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NUMBER_TO_CHARACTER_DP(NUMBER, FORMAT, ERR, ERROR)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    CHARACTER(LEN=MAXSTRLEN) :: NUMBER_TO_CHARACTER_DP !<On exit, the character equivalent of the number
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

  !>Converts an integer number to its equivalent varying string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g. "I3".
  FUNCTION NUMBER_TO_VSTRING_INTG(NUMBER,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_INTG !<On exit, the varying string equivalent of the number
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

  !>Converts a long integer number to its equivalent varying string representation as determined by the supplied format. The format is of the form of a standard Fortran integer format e.g., "I3".
  FUNCTION NUMBER_TO_VSTRING_LINTG(NUMBER,FORMAT,ERR,ERROR)
  
    !Argument variables
    INTEGER(LINTG), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_LINTG !<On exit, the varying string equivalent of the number
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

  !>Converts a single precision number to its equivalent varying string representation as determined by the supplied format string. Note If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NUMBER_TO_VSTRING_SP(NUMBER, FORMAT, ERR, ERROR)
  
    !Argument variables
    REAL(SP), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_SP !<On exit, the varying string equivalent of the number
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

  !>Converts a double precision number to its equivalent varying string representation as determined by the supplied format string. Note If FORMAT is an asterisk followed by a number between 1 and 32 the format will be chosen to maximise the number of significant digits, e.g., FORMAT="*8" will return a string of 8 characters representing the supplied number in either F8.? or E8.? format depending on its magnitude.
  FUNCTION NUMBER_TO_VSTRING_DP(NUMBER, FORMAT, ERR, ERROR)
      
    !Argument variables
    REAL(DP), INTENT(IN) :: NUMBER !<The number to convert
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT !<The format to use in the conversion
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    TYPE(VARYING_STRING) :: NUMBER_TO_VSTRING_DP !<On exit, the varying string equivalent of the number
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

  !>Converts a character string representation of a number to a double precision number.
  FUNCTION STRING_TO_DOUBLE_C(STRING, ERR, ERROR)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string !<The error string
    !Function variable
    REAL(DP) :: STRING_TO_DOUBLE_C !<On exit, the double precision equivalent of the string
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

  !>Converts a varying string representation of a number to a double precision number.
  FUNCTION STRING_TO_DOUBLE_VS(STRING, ERR, ERROR)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: STRING_TO_DOUBLE_VS !<On exit, the double precision equivalent of the string
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

  !>Converts a character string representation of a number to an integer.
  FUNCTION STRING_TO_INTEGER_C(STRING, ERR, ERROR)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: STRING_TO_INTEGER_C !<On exit, the integer equivalent of the string
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

  !>Converts a varying string representation of a number to an integer.
  FUNCTION STRING_TO_INTEGER_VS(STRING, ERR, ERROR)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: STRING_TO_INTEGER_VS !<On exit, the integer equivalent of the string
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

  !>Converts a character string representation of a number to a long integer.
  FUNCTION STRING_TO_LONG_INTEGER_C(STRING, ERR, ERROR)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(LINTG) :: STRING_TO_LONG_INTEGER_C !<On exit, the long integer equivalent of the string
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

  !>Converts a varying string representation of a number to a long integer.
  FUNCTION STRING_TO_LONG_INTEGER_VS(STRING, ERR, ERROR)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(LINTG) :: STRING_TO_LONG_INTEGER_VS !<On exit, the long integer equivalent of the string
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

  !>Converts a character string representation of a boolean (TRUE or FALSE) to a logical.
  FUNCTION STRING_TO_LOGICAL_C(STRING,ERR,ERROR)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    LOGICAL :: STRING_TO_LOGICAL_C !<On exit, the logical equivalent of the string
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

  !>Converts a varying string representation of a boolean (TRUE or FALSE) to a logical.
  FUNCTION STRING_TO_LOGICAL_VS(STRING,ERR,ERROR)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    LOGICAL :: STRING_TO_LOGICAL_VS !<On exit, the logical equivalent of the string
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

  !>Converts a character string representation of a number to a single precision number.
  FUNCTION STRING_TO_SINGLE_C(STRING, ERR, ERROR)
  
    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(SP) :: STRING_TO_SINGLE_C !<On exit, the single precision equivalent of the string
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

  !>Converts a varying string representation of a number to a single precision number.
  FUNCTION STRING_TO_SINGLE_VS(STRING, ERR, ERROR)
  
    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(SP) :: STRING_TO_SINGLE_VS !<On exit, the single precision equivalent of the string
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
    
  FUNCTION CHARACTER_TO_LOWERCASE_C(STRING)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert to lowercase
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_LOWERCASE_C !<On exit, the lowercase equivalent of the string
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

  !>Returns a character string that is the lowercase equivalent of the supplied varying string.
  FUNCTION CHARACTER_TO_LOWERCASE_VS(STRING)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_LOWERCASE_VS !<On exit, the lowercase equivalent of the string
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

  !>Returns a varying string that is the lowercase equivalent of the supplied character string.
  FUNCTION VSTRING_TO_LOWERCASE_C(STRING)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_LOWERCASE_C !<On exit, the lowercase equivalent of the string
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

  !>Returns a varying string that is the lowercase equivalent of the supplied varying string.
  FUNCTION VSTRING_TO_LOWERCASE_VS(STRING)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_LOWERCASE_VS !<On exit, the lowercase equivalent of the string
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

  !>Returns a character string which is uppercase equivalent of the supplied character string.
  FUNCTION CHARACTER_TO_UPPERCASE_C(STRING)

    !Argument variables 
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_UPPERCASE_C !<On exit, the uppercase equivalent of the string
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

  !>Returns a character string which is uppercase equivalent of the supplied varying string.
  FUNCTION CHARACTER_TO_UPPERCASE_VS(STRING)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    CHARACTER(LEN=LEN(STRING)) :: CHARACTER_TO_UPPERCASE_VS !<On exit, the uppercase equivalent of the string
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

  !>Returns a varying string which is uppercase equivalent of the supplied character string.
  FUNCTION VSTRING_TO_UPPERCASE_C(STRING)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_UPPERCASE_C !<On exit, the uppercase equivalent of the string
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

  !>Returns a varying string which is uppercase equivalent of the supplied varying string.
  FUNCTION VSTRING_TO_UPPERCASE_VS(STRING)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to convert
    !Function variable
    TYPE(VARYING_STRING) :: VSTRING_TO_UPPERCASE_VS !<On exit, the uppercase equivalent of the string
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
