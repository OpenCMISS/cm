  !
  !=================================================================================================================================
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

  END SUBROUTINE CMISSC2FString

  !
  !=================================================================================================================================
  !

  !>Copys/converts a  Fortran String (length of characters) to a C string (array of characters)
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

  END SUBROUTINE CMISSF2CString

  !
  !=================================================================================================================================
  !

  !>Copys/converts a list of C strings (2D array of characters) to an array of Fortran Strings
  SUBROUTINE CMISSC2FStrings(Cstrings,Fstrings)
    !Argument variables
    CHARACTER(LEN=1,KIND=C_CHAR), INTENT(IN) :: Cstrings(:,:)
    CHARACTER(LEN=*), INTENT(INOUT) :: Fstrings(:)
    !Local variables
    INTEGER(C_INT) :: string_idx,i
    INTEGER(C_INT) :: LENGTH=0

    !Cstrings array index order is opposite to C
    IF(LEN(Fstrings(1))>=SIZE(Cstrings,1)-1) THEN
      LENGTH=SIZE(Cstrings,1)-1
    ELSE
      LENGTH=LEN(Fstrings(1))
    ENDIF
    DO string_idx=1,SIZE(Fstrings,1)
      Fstrings(string_idx)=""
      DO i=1,LENGTH
        IF(Cstrings(i,string_idx)==C_NULL_CHAR) THEN
          EXIT
        ELSE
          Fstrings(string_idx)(i:i)=Cstrings(i,string_idx)
        ENDIF
      ENDDO !i
    ENDDO

    RETURN

  END SUBROUTINE CMISSC2FStrings

  !
  !=================================================================================================================================
  !

