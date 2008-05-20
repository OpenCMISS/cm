!> \file
!> $Id: input_output.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module handles all formating and input and output.
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

!> This module handles all formating and input and output.
MODULE INPUT_OUTPUT

  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
  USE ISO_VARYING_STRING
  USE STRINGS

  IMPLICIT NONE

  PRIVATE
  
  !Module parameters

  !> \addtogroup INPUT_OUTPUT_MatrixNameIndexFormat INPUT_OUTPUT::MatrixNameIndexFormat
  !> \brief Output type parameter
  !> \see INPUT_OUTPUT
  !>@{  
  INTEGER(INTG), PARAMETER :: WRITE_STRING_MATRIX_NAME_ONLY=1 !<Write the matrix name with out any indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
  INTEGER(INTG), PARAMETER :: WRITE_STRING_MATRIX_NAME_AND_INDICES=2 !<Write the matrix name together with the matrix indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
  !>@}

  !Module types

  !Interfaces

  !>Write a string to a given output stream
  INTERFACE WRITE_STRING
    MODULE PROCEDURE WRITE_STRING_C
    MODULE PROCEDURE WRITE_STRING_VS
  END INTERFACE !WRITE_STRING

  !>Write a string followed by a value to a given output stream
  INTERFACE WRITE_STRING_VALUE
    MODULE PROCEDURE WRITE_STRING_VALUE_C
    MODULE PROCEDURE WRITE_STRING_VALUE_DP
    MODULE PROCEDURE WRITE_STRING_VALUE_INTG
    MODULE PROCEDURE WRITE_STRING_VALUE_LINTG
    MODULE PROCEDURE WRITE_STRING_VALUE_L
    MODULE PROCEDURE WRITE_STRING_VALUE_SP
    MODULE PROCEDURE WRITE_STRING_VALUE_VS
  END INTERFACE !WRITE_STRING_VALUE

  !>Write a string, value, string then a value to a given output stream
  INTERFACE WRITE_STRING_TWO_VALUE
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_C_C
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_C_DP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_C_INTG
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_C_L
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_C_SP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_C_VS
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_DP_C
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_DP_DP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_DP_INTG
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_DP_L
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_DP_SP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_DP_VS
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_INTG_C
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_INTG_DP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_INTG_INTG
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_INTG_L
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_INTG_SP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_INTG_VS
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_L_C
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_L_DP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_L_INTG
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_L_L
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_L_SP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_L_VS
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_SP_C
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_SP_DP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_SP_INTG
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_SP_L
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_SP_SP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_SP_VS
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_VS_C
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_VS_DP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_VS_INTG
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_VS_L
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_VS_SP
    MODULE PROCEDURE WRITE_STRING_TWO_VALUE_VS_VS
  END INTERFACE !WRITE_STRING_TWO_VALUE

  !>Write a string followed by a value formatted in a particular way to a specified output stream
  INTERFACE WRITE_STRING_FMT_VALUE
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_C
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_DP
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_LINTG
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_L
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_SP
    MODULE PROCEDURE WRITE_STRING_FMT_VALUE_VS
  END INTERFACE !WRITE_STRING_FMT_VALUE
  
  !>Write a string, value, string then a value with the values formatted in a particular way to a given output stream
  INTERFACE WRITE_STRING_FMT_TWO_VALUE
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_C_C
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_C_DP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_C_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_C_L
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_C_SP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_C_VS
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_DP_C
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_DP_DP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_DP_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_DP_L
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_DP_SP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_DP_VS
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_INTG_C
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_INTG_DP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_INTG_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_INTG_L
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_INTG_SP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_INTG_VS
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_L_C
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_L_DP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_L_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_L_L
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_L_SP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_L_VS
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_SP_C
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_SP_DP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_SP_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_SP_L
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_SP_SP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_SP_VS
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_VS_C
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_VS_DP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_VS_INTG
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_VS_L
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_VS_SP
    MODULE PROCEDURE WRITE_STRING_FMT_TWO_VALUE_VS_VS
  END INTERFACE !WRITE_STRING_FMT_TWO_VALUE

  !>Write a string followed by a vector to a specified output stream.
  INTERFACE WRITE_STRING_VECTOR
    MODULE PROCEDURE WRITE_STRING_VECTOR_DP
    MODULE PROCEDURE WRITE_STRING_VECTOR_INTG
    MODULE PROCEDURE WRITE_STRING_VECTOR_LINTG
    MODULE PROCEDURE WRITE_STRING_VECTOR_L
    MODULE PROCEDURE WRITE_STRING_VECTOR_SP
  END INTERFACE !WRITE_STRING_VECTOR
  
  !>Write a string followed by a indexed vector to a specified output stream.
  INTERFACE WRITE_STRING_IDX_VECTOR
    MODULE PROCEDURE WRITE_STRING_IDX_VECTOR_DP
    MODULE PROCEDURE WRITE_STRING_IDX_VECTOR_INTG
    MODULE PROCEDURE WRITE_STRING_IDX_VECTOR_LINTG
    MODULE PROCEDURE WRITE_STRING_IDX_VECTOR_L
    MODULE PROCEDURE WRITE_STRING_IDX_VECTOR_SP
  END INTERFACE !WRITE_STRING_IDX_VECTOR

  !>Write a string followed by a matrix to a specified output stream
  INTERFACE WRITE_STRING_MATRIX
    MODULE PROCEDURE WRITE_STRING_MATRIX_DP
    MODULE PROCEDURE WRITE_STRING_MATRIX_INTG
    MODULE PROCEDURE WRITE_STRING_MATRIX_LINTG
    MODULE PROCEDURE WRITE_STRING_MATRIX_L
    MODULE PROCEDURE WRITE_STRING_MATRIX_SP
  END INTERFACE !WRITE_STRING_MATRIX

  PUBLIC WRITE_STRING_MATRIX_NAME_ONLY,WRITE_STRING_MATRIX_NAME_AND_INDICES
  
  PUBLIC WRITE_STRING,WRITE_STRING_VALUE,WRITE_STRING_TWO_VALUE,WRITE_STRING_FMT_VALUE,WRITE_STRING_FMT_TWO_VALUE, &
    & WRITE_STRING_VECTOR,WRITE_STRING_IDX_VECTOR,WRITE_STRING_MATRIX

  !Module variables

CONTAINS

  !!TODO: put back enters,exits etc.
 
  !
  !================================================================================================================================
  !

  !>Writes the character STRING to the given output stream specified by ID.
  SUBROUTINE WRITE_STRING_C(ID,STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream to write to \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The string to write
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

!    CALL ENTERS("WRITE_STRING_C",ERR,ERROR,*999)
        
    WRITE(OP_STRING,'(A)') STRING
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_C

  !
  !================================================================================================================================
  !

  !>Writes the varying string STRING to the given output stream specified by ID.
  SUBROUTINE WRITE_STRING_VS(ID,STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream to write to \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    TYPE(VARYING_STRING), INTENT(IN) :: STRING !<The string to write
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables

!    CALL ENTERS("WRITE_STRING_VS",ERR,ERROR,*999)
        
    WRITE(OP_STRING,'(A)') CHAR(STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_C(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_VALUE_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//VALUE(1:LEN_TRIM(VALUE))
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_DP(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_VALUE_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_INTG(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_VALUE_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_LINTG(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(LINTG), INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_VALUE_LINTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_LINTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_LINTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_LINTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_LINTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_L(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_VALUE_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_SP(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_VALUE_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. Free format is used to format the value.
  SUBROUTINE WRITE_STRING_VALUE_VS(ID,FIRST_STRING,VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: VALUE !<The value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_VALUE_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_VALUE_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_VALUE_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VALUE_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_VALUE_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_C_C(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_C_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_C_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_C_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_C_DP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_C_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_C_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_C_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_C_INTG(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_C_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_C_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_C_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_C_L(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_C_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_C_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_C_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_C_SP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_C_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_C_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_C_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted varying string SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_C_VS(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_C_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_C_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_C_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_C_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_DP_C(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_DP_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_DP_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_DP_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_DP_DP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_DP_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_DP_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_DP_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_DP_INTG(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_DP_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_DP_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_DP_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_DP_L(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_DP_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_DP_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_DP_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_DP_SP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_DP_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_DP_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_DP_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_DP_VS(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_DP_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_DP_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_DP_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_DP_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_C(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_INTG_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_INTG_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_DP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_INTG_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_INTG_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_INTG(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_INTG_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_INTG_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_L(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_INTG_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_INTG_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_SP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_INTG_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_INTG_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_VS(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_INTG_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_INTG_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_INTG_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_INTG_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_L_C(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_L_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_L_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_L_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_L_DP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_L_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_L_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_L_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_L_INTG(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_L_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    LOCAL_STRING=LOCAL_STRING//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_L_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_L_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_L_L(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_L_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_L_INTG_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_L_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_L_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_L_SP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_L_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    LOCAL_STRING=LOCAL_STRING//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_L_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_L_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_L_VS(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_L_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_L_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_L_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_L_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_SP_C(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_SP_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_SP_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_SP_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_SP_DP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_SP_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_SP_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_SP_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_SP_INTG(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_SP_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_SP_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_SP_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_SP_L(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_SP_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    LOCAL_STRING=LOCAL_STRING//SECOND_STRING//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_SP_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_SP_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_SP_SP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_SP_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_SP_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_SP_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_SP_VS(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_SP_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,"*",ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_SP_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_SP_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_SP_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_VS_C(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_VS_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_VS_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_VS_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_VS_DP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_VS_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_VS_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_VS_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_VS_INTG(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_VS_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_VS_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_VS_INTG

  !
  !================================================================================================================================
  !
  
  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_VS_L(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_VS_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_VS_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_VS_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_VS_SP(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_TWO_VALUE_VS_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,"*",ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_VS_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_VS_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted varying string SECOND_VALUE to the given output stream specified by ID. Free format is used to format both values.
  SUBROUTINE WRITE_STRING_TWO_VALUE_VS_VS(ID,FIRST_STRING,FIRST_VALUE,SECOND_STRING,SECOND_VALUE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_TWO_VALUE_VS_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_TWO_VALUE_VS_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_TWO_VALUE_VS_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_TWO_VALUE_VS_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_C(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING(1:LEN_TRIM(FIRST_STRING))//VALUE(1:LEN_TRIM(VALUE))
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_DP(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,FORMAT_STRING,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_DP
  
  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_INTG(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,FORMAT_STRING,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_LINTG(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(LINTG), INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_LINTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,FORMAT_STRING,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_LINTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_LINTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_LINTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_LINTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_L(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_SP(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(VALUE,FORMAT_STRING,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST STRING followed by a formatted character VALUE to the given output stream specified by ID. FORMAT_STRING is used to format the value.
  SUBROUTINE WRITE_STRING_FMT_VALUE_VS(ID,FIRST_STRING,VALUE,FORMAT_STRING,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: VALUE !<The value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FORMAT_STRING !<The format string to be used to format the value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_VALUE_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_VALUE_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_VALUE_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_VALUE_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_VALUE_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_C(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_C_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_C_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_DP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_C_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_C_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_INTG(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_C_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_C_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_L(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_C_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_C_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_SP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_C_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_C_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted character FIRST_VALUE and the the SECOND_STRING followed by a formatted varying string SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_VS(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_C_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_C_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_C_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_C_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_C(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_DP_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_DP_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_DP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_DP_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_DP_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_INTG(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_DP_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_DP_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_L(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_DP_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_DP_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_SP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_DP_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_DP_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted double precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_VS(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(DP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_DP_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_DP_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_DP_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_DP_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_C(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_INTG_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_INTG_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_DP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_INTG_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_INTG_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_INTG(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_INTG_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_INTG_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_L(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_INTG_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_INTG_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_SP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_INTG_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_INTG_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted integer FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_VS(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)
   
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_INTG_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_INTG_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_INTG_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_INTG_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_C(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_L_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_L_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_DP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_L_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_L_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_INTG(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_L_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    LOCAL_STRING=LOCAL_STRING//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_L_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_L(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_L_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_TWO_L_INTG_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_L_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_SP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_L_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    LOCAL_STRING=LOCAL_STRING//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_L_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted logical FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_VS(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    LOGICAL, INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_L_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//LOGICAL_TO_VSTRING(FIRST_VALUE,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_L_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_L_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_L_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_C(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_SP_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_SP_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_DP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
   
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_SP_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_SP_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_INTG(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_SP_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_SP_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_L(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_SP_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    LOCAL_STRING=LOCAL_STRING//SECOND_STRING//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_SP_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_SP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING,LOCAL_STRING2

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_SP_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    !CPB 21/02/2007 AIX doesn't like concatenating vstrings and reassigning to itself so split this into two steps
    LOCAL_STRING2=LOCAL_STRING//SECOND_STRING
    LOCAL_STRING=LOCAL_STRING2//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_SP_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted single precision FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_VS(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    REAL(SP), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_SP_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//NUMBER_TO_VSTRING(FIRST_VALUE,FIRST_FORMAT,ERR,ERROR)//SECOND_STRING//SECOND_VALUE
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_SP_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_SP_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_SP_VS

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted character SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_C(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_VS_C",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_C")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_VS_C",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_C")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_C

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted double precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_DP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(DP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_VS_DP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_VS_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_DP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_DP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted integer SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_INTG(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE, &
    & SECOND_FORMAT,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    INTEGER(INTG), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_VS_INTG",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_VS_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_INTG")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_INTG

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted logical SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_L(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    LOGICAL, INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_VS_L",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//LOGICAL_TO_VSTRING(SECOND_VALUE,ERR,ERROR)
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_VS_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_L")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_L

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted single precision SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_SP(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    REAL(SP), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

 !   CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_VS_SP",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//NUMBER_TO_VSTRING(SECOND_VALUE,SECOND_FORMAT,ERR,ERROR)
    IF(ERR/=0) GOTO 999
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_VS_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_SP")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_SP

  !
  !================================================================================================================================
  !

  !>Writes the FIRST_STRING followed by a formatted varying string FIRST_VALUE and the the SECOND_STRING followed by a formatted varying string SECOND_VALUE to the given output stream specified by ID. FIRST_FORMAT is used to format the first value and SECOND_FORMAT is used to format the second value.
  SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_VS(ID,FIRST_STRING,FIRST_VALUE,FIRST_FORMAT,SECOND_STRING,SECOND_VALUE,SECOND_FORMAT, &
    & ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_STRING !<The first string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: FIRST_VALUE !<The first value to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used to format the first value
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_STRING !<The second string to be output
    TYPE(VARYING_STRING), INTENT(IN) :: SECOND_VALUE !<The second value to be output
    CHARACTER(LEN=*), INTENT(IN) :: SECOND_FORMAT !<The format string to be used to format the second value
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_STRING

!    CALL ENTERS("WRITE_STRING_FMT_TWO_VALUE_VS_VS",ERR,ERROR,*999)
        
    LOCAL_STRING=FIRST_STRING//FIRST_VALUE//SECOND_STRING//SECOND_VALUE
    WRITE(OP_STRING,'(A)') CHAR(LOCAL_STRING)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
      
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_VS")
    RETURN
999 CALL ERRORS("WRITE_STRING_FMT_TWO_VALUE_VS_VS",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_FMT_TWO_VALUE_VS_VS")
    RETURN 1   
  END SUBROUTINE WRITE_STRING_FMT_TWO_VALUE_VS_VS

  !
  !================================================================================================================================
  !

  !>Writes the given double precision VECTOR to the given output stream specified by ID. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_IDX and LAST_IDX are the extents of the data and DELTA is the NUMBER of indices to skip for each index.
  SUBROUTINE WRITE_STRING_VECTOR_DP(ID,FIRST_IDX,DELTA,LAST_IDX,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT,REPEAT_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_IDX !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: LAST_IDX !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(DP), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    CALL ENTERS("WRITE_STRING_VECTOR_DP",ERR,ERROR,*999)
        
    current=FIRST_IDX
    final=current+(NUMBER_FIRST-1)*DELTA
    IF(final>LAST_IDX) final=LAST_IDX
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR(count),count=current,final,DELTA)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    DO WHILE(final<LAST_IDX) !more stuff to do
      current=final+DELTA
      final=final+NUMBER_REPEAT*DELTA
      IF(final>LAST_IDX) final=LAST_IDX
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR(count),count=current,final,DELTA)
      CALL WRITE_STR(ID,ERR,ERROR,*999)
    ENDDO !final<LAST_IDX

!    CALL EXITS("WRITE_STRING_VECTOR_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_VECTOR_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VECTOR_DP")
    RETURN 1
  END SUBROUTINE WRITE_STRING_VECTOR_DP

  !
  !================================================================================================================================
  !

  !>Writes the given integer VECTOR to the given output stream specified by ID. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_IDX and LAST_IDX are the extents of the data and DELTA is the NUMBER of indices to skip for each index.
  SUBROUTINE WRITE_STRING_VECTOR_INTG(ID,FIRST_IDX,DELTA,LAST_IDX,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT,REPEAT_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_IDX !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: LAST_IDX !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(INTG), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    CALL ENTERS("WRITE_STRING_VECTOR_INTG",ERR,ERROR,*999)
        
    current=FIRST_IDX
    final=current+(NUMBER_FIRST-1)*DELTA
    IF(final>LAST_IDX) final=LAST_IDX
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR(count),count=current,final,DELTA)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    DO WHILE(final<LAST_IDX) !more stuff to do
      current=final+DELTA
      final=final+NUMBER_REPEAT*DELTA
      IF(final>LAST_IDX) final=LAST_IDX
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR(count),count=current,final,DELTA)
      CALL WRITE_STR(ID,ERR,ERROR,*999)
    ENDDO !final<LAST_IDX

!    CALL EXITS("WRITE_STRING_VECTOR_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_VECTOR_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VECTOR_INTG")
    RETURN 1
  END SUBROUTINE WRITE_STRING_VECTOR_INTG

  !
  !================================================================================================================================
  !

  !>Writes the given integer VECTOR to the given output stream specified by ID. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_IDX and LAST_IDX are the extents of the data and DELTA is the NUMBER of indices to skip for each index.
  SUBROUTINE WRITE_STRING_VECTOR_LINTG(ID,FIRST_IDX,DELTA,LAST_IDX,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT,REPEAT_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_IDX !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: LAST_IDX !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(LINTG), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    CALL ENTERS("WRITE_STRING_VECTOR_LINTG",ERR,ERROR,*999)
        
    current=FIRST_IDX
    final=current+(NUMBER_FIRST-1)*DELTA
    IF(final>LAST_IDX) final=LAST_IDX
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR(count),count=current,final,DELTA)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    DO WHILE(final<LAST_IDX) !more stuff to do
      current=final+DELTA
      final=final+NUMBER_REPEAT*DELTA
      IF(final>LAST_IDX) final=LAST_IDX
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR(count),count=current,final,DELTA)
      CALL WRITE_STR(ID,ERR,ERROR,*999)
    ENDDO !final<LAST_IDX

!    CALL EXITS("WRITE_STRING_VECTOR_LINTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_VECTOR_LINTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VECTOR_LINTG")
    RETURN 1
  END SUBROUTINE WRITE_STRING_VECTOR_LINTG

  !
  !================================================================================================================================
  !

  !>Writes the given logical VECTOR to the given output stream specified by ID. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_IDX and LAST_IDX are the extents of the data and DELTA is the NUMBER of indices to skip for each index.
  SUBROUTINE WRITE_STRING_VECTOR_L(ID,FIRST_IDX,DELTA,LAST_IDX,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT,REPEAT_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_IDX !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: LAST_IDX !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    LOGICAL, INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    CALL ENTERS("WRITE_STRING_VECTOR_L",ERR,ERROR,*999)
        
    current=FIRST_IDX
    final=current+(NUMBER_FIRST-1)*DELTA
    IF(final>LAST_IDX) final=LAST_IDX
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR(count),count=current,final,DELTA)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    DO WHILE(final<LAST_IDX) !more stuff to do
      current=final+DELTA
      final=final+NUMBER_REPEAT*DELTA
      IF(final>LAST_IDX) final=LAST_IDX
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR(count),count=current,final,DELTA)
      CALL WRITE_STR(ID,ERR,ERROR,*999)
    ENDDO !final<LAST_IDX

!    CALL EXITS("WRITE_STRING_VECTOR_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_VECTOR_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VECTOR_L")
    RETURN 1
  END SUBROUTINE WRITE_STRING_VECTOR_L

  !
  !================================================================================================================================
  !

  !>Writes the given single precision VECTOR to the given output stream specified by ID. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_IDX and LAST_IDX are the extents of the data and DELTA is the NUMBER of indices to skip for each index.
  SUBROUTINE WRITE_STRING_VECTOR_SP(ID,FIRST_IDX,DELTA,LAST_IDX,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT,REPEAT_FORMAT, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_IDX !<The first index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: LAST_IDX !<The last index of the vector to output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(SP), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,final,count

!    CALL ENTERS("WRITE_STRING_VECTOR_SP",ERR,ERROR,*999)
        
    current=FIRST_IDX
    final=current+(NUMBER_FIRST-1)*DELTA
    IF(final>LAST_IDX) final=LAST_IDX
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR(count),count=current,final,DELTA)
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    DO WHILE(final<LAST_IDX) !more stuff to do
      current=final+DELTA
      final=final+NUMBER_REPEAT*DELTA
      IF(final>LAST_IDX) final=LAST_IDX
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR(count),count=current,final,DELTA)
      CALL WRITE_STR(ID,ERR,ERROR,*999)
    ENDDO !final<LAST_IDX

!    CALL EXITS("WRITE_STRING_VECTOR_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_VECTOR_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_VECTOR_SP")
    RETURN 1
  END SUBROUTINE WRITE_STRING_VECTOR_SP

  !
  !================================================================================================================================
  !

  !>Writes the given indexed double precision VECTOR to the given output stream specified by ID. NUM_INDICES is the number of indices and INDICES(i) contain the indices of the vector to write. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. DELTA is the number of actual indices to skip for each index.
  SUBROUTINE WRITE_STRING_IDX_VECTOR_DP(ID,NUM_INDICES,INDICES,DELTA,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT, &
    & REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: NUM_INDICES !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: INDICES(NUM_INDICES) !<INDICES(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(DP), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,number_to_do

!    CALL ENTERS("WRITE_STRING_IDX_VECTOR_DP",ERR,ERROR,*999)
        
    number_to_do=NUM_INDICES
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=1,MIN(NUMBER_FIRST,NUM_INDICES))
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    number_to_do=NUM_INDICES-NUMBER_FIRST
    current=NUMBER_FIRST+1
    DO WHILE(number_to_do>0) !more stuff to do
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=current,MIN(current+NUMBER_REPEAT-1, &
        & NUM_INDICES))
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      current=current+NUMBER_REPEAT
      number_to_do=number_to_do-NUMBER_REPEAT
    ENDDO !number_to_do > 0

!    CALL EXITS("WRITE_STRING_IDX_VECTOR_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_IDX_VECTOR_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_IDX_VECTOR_DP")
    RETURN 1
  END SUBROUTINE WRITE_STRING_IDX_VECTOR_DP

  !
  !================================================================================================================================
  !

  !>Writes the given indexed integer VECTOR to the given output stream specified by ID. NUM_INDICES is the number of indices and INDICES(i) contain the indices of the vector to write. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. DELTA is the number of actual indices to skip for each index.
  SUBROUTINE WRITE_STRING_IDX_VECTOR_INTG(ID,NUM_INDICES,INDICES,DELTA,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT, &
    & REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: NUM_INDICES !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: INDICES(NUM_INDICES) !<INDICES(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(INTG), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,number_to_do

!    CALL ENTERS("WRITE_STRING_IDX_VECTOR_INTG",ERR,ERROR,*999)
        
    number_to_do=NUM_INDICES
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=1,MIN(NUMBER_FIRST,NUM_INDICES))
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    number_to_do=NUM_INDICES-NUMBER_FIRST
    current=NUMBER_FIRST+1
    DO WHILE(number_to_do>0) !more stuff to do
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=current,MIN(current+NUMBER_REPEAT-1, &
        & NUM_INDICES))
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      current=current+NUMBER_REPEAT
      number_to_do=number_to_do-NUMBER_REPEAT
    ENDDO !number_to_do > 0

!    CALL EXITS("WRITE_STRING_IDX_VECTOR_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_IDX_VECTOR_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_IDX_VECTOR_INTG")
    RETURN 1
  END SUBROUTINE WRITE_STRING_IDX_VECTOR_INTG

  !
  !================================================================================================================================
  !

  !>Writes the given indexed integer VECTOR to the given output stream specified by ID. NUM_INDICES is the number of indices and INDICES(i) contain the indices of the vector to write. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. DELTA is the number of actual indices to skip for each index.
  SUBROUTINE WRITE_STRING_IDX_VECTOR_LINTG(ID,NUM_INDICES,INDICES,DELTA,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT, &
    & REPEAT_FORMAT,ERR,ERROR,*)

    !#### Generic-Subroutine: WRITE_STRING_IDX_VECTOR_LINTG
    !###  Description:
    !###    
    !###  Parent-subroutines: WRITE_STRING_IDX_VECTOR

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: NUM_INDICES !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: INDICES(NUM_INDICES) !<INDICES(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    INTEGER(LINTG), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,number_to_do

!    CALL ENTERS("WRITE_STRING_IDX_VECTOR_LINTG",ERR,ERROR,*999)
        
    number_to_do=NUM_INDICES
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=1,MIN(NUMBER_FIRST,NUM_INDICES))
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    number_to_do=NUM_INDICES-NUMBER_FIRST
    current=NUMBER_FIRST+1
    DO WHILE(number_to_do>0) !more stuff to do
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=current,MIN(current+NUMBER_REPEAT-1, &
        & NUM_INDICES))
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      current=current+NUMBER_REPEAT
      number_to_do=number_to_do-NUMBER_REPEAT
    ENDDO !number_to_do > 0

!    CALL EXITS("WRITE_STRING_IDX_VECTOR_LINTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_IDX_VECTOR_LINTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_IDX_VECTOR_LINTG")
    RETURN 1
  END SUBROUTINE WRITE_STRING_IDX_VECTOR_LINTG

  !
  !================================================================================================================================
  !

  !>Writes the given indexed logical VECTOR to the given output stream specified by ID. NUM_INDICES is the number of indices and INDICES(i) contain the indices of the vector to write. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. DELTA is the number of actual indices to skip for each index.
  SUBROUTINE WRITE_STRING_IDX_VECTOR_L(ID,NUM_INDICES,INDICES,DELTA,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT, &
    & REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: NUM_INDICES !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: INDICES(NUM_INDICES) !<INDICES(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    LOGICAL, INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,number_to_do

!    CALL ENTERS("WRITE_STRING_IDX_VECTOR_L",ERR,ERROR,*999)
        
    number_to_do=NUM_INDICES
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=1,MIN(NUMBER_FIRST,NUM_INDICES))
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    number_to_do=NUM_INDICES-NUMBER_FIRST
    current=NUMBER_FIRST+1
    DO WHILE(number_to_do>0) !more stuff to do
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=current,MIN(current+NUMBER_REPEAT-1, &
        & NUM_INDICES))
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      current=current+NUMBER_REPEAT
      number_to_do=number_to_do-NUMBER_REPEAT
    ENDDO !number_to_do > 0

!    CALL EXITS("WRITE_STRING_IDX_VECTOR_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_IDX_VECTOR_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_IDX_VECTOR_L")
    RETURN 1
  END SUBROUTINE WRITE_STRING_IDX_VECTOR_L

  !
  !================================================================================================================================
  !

  !>Writes the given indexed single precision VECTOR to the given output stream specified by ID. NUM_INDICES is the number of indices and INDICES(i) contain the indices of the vector to write. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. DELTA is the number of actual indices to skip for each index.
  SUBROUTINE WRITE_STRING_IDX_VECTOR_SP(ID,NUM_INDICES,INDICES,DELTA,NUMBER_FIRST,NUMBER_REPEAT,VECTOR,FIRST_FORMAT, &
    & REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: NUM_INDICES !<The number of indices of the vector to output
    INTEGER(INTG), INTENT(IN) :: INDICES(NUM_INDICES) !<INDICES(i). The i'th index of the vector to output
    INTEGER(INTG), INTENT(IN) :: DELTA !<The delta increment to be used when outputing the first through to the last vector index
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of vector elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of vector elements to be output on the second and subsequently repeated lines
    REAL(SP), INTENT(IN) :: VECTOR(:) !<The vector to be output
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current,count,number_to_do

!    CALL ENTERS("WRITE_STRING_IDX_VECTOR_SP",ERR,ERROR,*999)
        
    number_to_do=NUM_INDICES
    WRITE(OP_STRING,FMT=FIRST_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=1,MIN(NUMBER_FIRST,NUM_INDICES))
    CALL WRITE_STR(ID,ERR,ERROR,*999)
    number_to_do=NUM_INDICES-NUMBER_FIRST
    current=NUMBER_FIRST+1
    DO WHILE(number_to_do>0) !more stuff to do
      WRITE(OP_STRING,FMT=REPEAT_FORMAT) (VECTOR((INDICES(count)-1)*DELTA+1),count=current,MIN(current+NUMBER_REPEAT-1, &
        & NUM_INDICES))
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      current=current+NUMBER_REPEAT
      number_to_do=number_to_do-NUMBER_REPEAT
    ENDDO !number_to_do > 0

!    CALL EXITS("WRITE_STRING_IDX_VECTOR_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_IDX_VECTOR_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_IDX_VECTOR_SP")
    RETURN 1
  END SUBROUTINE WRITE_STRING_IDX_VECTOR_SP

  !
  !================================================================================================================================
  !

  !>Writes the given double precision MATRIX to the given output stream specified by ID. The basic output is determined by the flag INDEX_FORMAT_TYPE. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is MATRIX_NAME_FORMAT concatenated named with the FIRST_FORMAT. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is MATRIX_NAME_FORMAT concatenated with ROW_INDEX_FORMAT and concatenated with FIRST_FORMAT. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_ROW/FIRST_COLUMN and LAST_ROW/LAST_COLUMN are the extents of the row/column and DELTA_ROW/DELTA_COLUMN is the NUMBER of indices to skip for each row/column index.
  SUBROUTINE WRITE_STRING_MATRIX_DP(ID,FIRST_ROW,DELTA_ROW,LAST_ROW,FIRST_COLUMN,DELTA_COLUMN,LAST_COLUMN,NUMBER_FIRST, &
    & NUMBER_REPEAT,MATRIX,INDEX_FORMAT_TYPE,MATRIX_NAME_FORMAT,ROW_INDEX_FORMAT,FIRST_FORMAT,REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_ROW !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_ROW !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: LAST_ROW !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_COLUMN !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_COLUMN !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: LAST_COLUMN !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of matrix elements to be output on the second and subsequently repeated lines
    REAL(DP), INTENT(IN) :: MATRIX(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: INDEX_FORMAT_TYPE !<The format type to be used for the matrix name and indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: MATRIX_NAME_FORMAT !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: ROW_INDEX_FORMAT !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current_row,current_column,final_column,count
    CHARACTER(LEN=MAXSTRLEN) :: FORMAT_STR

!    CALL ENTERS("WRITE_STRING_MATRIX_DP",ERR,ERROR,*999)

    IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//FIRST_FORMAT
    ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//ROW_INDEX_FORMAT//FIRST_FORMAT
    ELSE
      CALL FLAG_ERROR("Invalid index format type",ERR,ERROR,*999)
    ENDIF
    DO current_row=FIRST_ROW,LAST_ROW,DELTA_ROW
      current_column=FIRST_COLUMN
      final_column=current_column+(NUMBER_FIRST-1)*DELTA_COLUMN
      IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
      IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) current_row,(MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ENDIF
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      DO WHILE(final_column<LAST_COLUMN) !more stuff to do
        current_column=final_column+DELTA_COLUMN
        final_column=final_column+NUMBER_REPEAT*DELTA_COLUMN
        IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
        WRITE(OP_STRING,FMT=REPEAT_FORMAT) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
        CALL WRITE_STR(ID,ERR,ERROR,*999)
      ENDDO !final_columnn<LAST_COLUMN
    ENDDO !current_row
    
!    CALL EXITS("WRITE_STRING_MATRIX_DP")
    RETURN
999 CALL ERRORS("WRITE_STRING_MATRIX_DP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_MATRIX_DP")
    RETURN 1
  END SUBROUTINE WRITE_STRING_MATRIX_DP

  !
  !================================================================================================================================
  !

  !>Writes the given integer MATRIX to the given output stream specified by ID. The basic output is determined by the flag INDEX_FORMAT_TYPE. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is MATRIX_NAME_FORMAT concatenated named with the FIRST_FORMAT. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is MATRIX_NAME_FORMAT concatenated with ROW_INDEX_FORMAT and concatenated with FIRST_FORMAT. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_ROW/FIRST_COLUMN and LAST_ROW/LAST_COLUMN are the extents of the row/column and DELTA_ROW/DELTA_COLUMN is the NUMBER of indices to skip for each row/column index.
  SUBROUTINE WRITE_STRING_MATRIX_INTG(ID,FIRST_ROW,DELTA_ROW,LAST_ROW,FIRST_COLUMN,DELTA_COLUMN,LAST_COLUMN,NUMBER_FIRST, &
    & NUMBER_REPEAT,MATRIX,INDEX_FORMAT_TYPE,MATRIX_NAME_FORMAT,ROW_INDEX_FORMAT,FIRST_FORMAT,REPEAT_FORMAT,ERR,ERROR,*)

    !#### Generic-Subroutine: WRITE_STRING_MATRIX_INTG
    !###  Description:
    !###    
    !###  Parent-subroutines: WRITE_STRING_MATRIX

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_ROW !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_ROW !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: LAST_ROW !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_COLUMN !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_COLUMN !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: LAST_COLUMN !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of matrix elements to be output on the second and subsequently repeated lines
    INTEGER(INTG), INTENT(IN) :: MATRIX(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: INDEX_FORMAT_TYPE !<The format type to be used for the matrix name and indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: MATRIX_NAME_FORMAT !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: ROW_INDEX_FORMAT !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current_row,current_column,final_column,count
    CHARACTER(LEN=MAXSTRLEN) :: FORMAT_STR

!    CALL ENTERS("WRITE_STRING_MATRIX_INTG",ERR,ERROR,*999)

    IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//FIRST_FORMAT
    ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//ROW_INDEX_FORMAT//FIRST_FORMAT
    ELSE
      CALL FLAG_ERROR("Invalid index format type",ERR,ERROR,*999)
    ENDIF
    DO current_row=FIRST_ROW,LAST_ROW,DELTA_ROW
      current_column=FIRST_COLUMN
      final_column=current_column+(NUMBER_FIRST-1)*DELTA_COLUMN
      IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
      IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) current_row,(MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ENDIF
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      DO WHILE(final_column<LAST_COLUMN) !more stuff to do
        current_column=final_column+DELTA_COLUMN
        final_column=final_column+NUMBER_REPEAT*DELTA_COLUMN
        IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
        WRITE(OP_STRING,FMT=REPEAT_FORMAT) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
        CALL WRITE_STR(ID,ERR,ERROR,*999)
      ENDDO !final_columnn<LAST_COLUMN
    ENDDO !current_row
    
!    CALL EXITS("WRITE_STRING_MATRIX_INTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_MATRIX_INTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_MATRIX_INTG")
    RETURN 1
  END SUBROUTINE WRITE_STRING_MATRIX_INTG

  !
  !================================================================================================================================
  !
  

  !>Writes the given long integer MATRIX to the given output stream specified by ID. The basic output is determined by the flag INDEX_FORMAT_TYPE. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is MATRIX_NAME_FORMAT concatenated named with the FIRST_FORMAT. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is MATRIX_NAME_FORMAT concatenated with ROW_INDEX_FORMAT and concatenated with FIRST_FORMAT. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_ROW/FIRST_COLUMN and LAST_ROW/LAST_COLUMN are the extents of the row/column and DELTA_ROW/DELTA_COLUMN is the NUMBER of indices to skip for each row/column index.
  SUBROUTINE WRITE_STRING_MATRIX_LINTG(ID,FIRST_ROW,DELTA_ROW,LAST_ROW,FIRST_COLUMN,DELTA_COLUMN,LAST_COLUMN,NUMBER_FIRST, &
    & NUMBER_REPEAT,MATRIX,INDEX_FORMAT_TYPE,MATRIX_NAME_FORMAT,ROW_INDEX_FORMAT,FIRST_FORMAT,REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_ROW !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_ROW !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: LAST_ROW !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_COLUMN !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_COLUMN !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: LAST_COLUMN !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of matrix elements to be output on the second and subsequently repeated lines
    INTEGER(LINTG), INTENT(IN) :: MATRIX(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: INDEX_FORMAT_TYPE !<The format type to be used for the matrix name and indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: MATRIX_NAME_FORMAT !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: ROW_INDEX_FORMAT !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current_row,current_column,final_column,count
    CHARACTER(LEN=MAXSTRLEN) :: FORMAT_STR

!    CALL ENTERS("WRITE_STRING_MATRIX_LINTG",ERR,ERROR,*999)

    IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//FIRST_FORMAT
    ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//ROW_INDEX_FORMAT//FIRST_FORMAT
    ELSE
      CALL FLAG_ERROR("Invalid index format type",ERR,ERROR,*999)
    ENDIF
    DO current_row=FIRST_ROW,LAST_ROW,DELTA_ROW
      current_column=FIRST_COLUMN
      final_column=current_column+(NUMBER_FIRST-1)*DELTA_COLUMN
      IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
      IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) current_row,(MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ENDIF
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      DO WHILE(final_column<LAST_COLUMN) !more stuff to do
        current_column=final_column+DELTA_COLUMN
        final_column=final_column+NUMBER_REPEAT*DELTA_COLUMN
        IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
        WRITE(OP_STRING,FMT=REPEAT_FORMAT) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
        CALL WRITE_STR(ID,ERR,ERROR,*999)
      ENDDO !final_columnn<LAST_COLUMN
    ENDDO !current_row
    
!    CALL EXITS("WRITE_STRING_MATRIX_LINTG")
    RETURN
999 CALL ERRORS("WRITE_STRING_MATRIX_LINTG",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_MATRIX_LINTG")
    RETURN 1
  END SUBROUTINE WRITE_STRING_MATRIX_LINTG

  !
  !================================================================================================================================
  !

  !>Writes the given logical MATRIX to the given output stream specified by ID. The basic output is determined by the flag INDEX_FORMAT_TYPE. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is MATRIX_NAME_FORMAT concatenated named with the FIRST_FORMAT. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is MATRIX_NAME_FORMAT concatenated with ROW_INDEX_FORMAT and concatenated with FIRST_FORMAT. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_ROW/FIRST_COLUMN and LAST_ROW/LAST_COLUMN are the extents of the row/column and DELTA_ROW/DELTA_COLUMN is the NUMBER of indices to skip for each row/column index.
  SUBROUTINE WRITE_STRING_MATRIX_L(ID,FIRST_ROW,DELTA_ROW,LAST_ROW,FIRST_COLUMN,DELTA_COLUMN,LAST_COLUMN,NUMBER_FIRST, &
    & NUMBER_REPEAT,MATRIX,INDEX_FORMAT_TYPE,MATRIX_NAME_FORMAT,ROW_INDEX_FORMAT,FIRST_FORMAT,REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_ROW !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_ROW !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: LAST_ROW !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_COLUMN !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_COLUMN !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: LAST_COLUMN !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of matrix elements to be output on the second and subsequently repeated lines
    LOGICAL, INTENT(IN) :: MATRIX(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: INDEX_FORMAT_TYPE !<The format type to be used for the matrix name and indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: MATRIX_NAME_FORMAT !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: ROW_INDEX_FORMAT !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current_row,current_column,final_column,count
    CHARACTER(LEN=MAXSTRLEN) :: FORMAT_STR

!    CALL ENTERS("WRITE_STRING_MATRIX_L",ERR,ERROR,*999)

    IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//FIRST_FORMAT
    ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//ROW_INDEX_FORMAT//FIRST_FORMAT
    ELSE
      CALL FLAG_ERROR("Invalid index format type",ERR,ERROR,*999)
    ENDIF
    DO current_row=FIRST_ROW,LAST_ROW,DELTA_ROW
      current_column=FIRST_COLUMN
      final_column=current_column+(NUMBER_FIRST-1)*DELTA_COLUMN
      IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
      IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) current_row,(MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ENDIF
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      DO WHILE(final_column<LAST_COLUMN) !more stuff to do
        current_column=final_column+DELTA_COLUMN
        final_column=final_column+NUMBER_REPEAT*DELTA_COLUMN
        IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
        WRITE(OP_STRING,FMT=REPEAT_FORMAT) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
        CALL WRITE_STR(ID,ERR,ERROR,*999)
      ENDDO !final_columnn<LAST_COLUMN
    ENDDO !current_row
    
!    CALL EXITS("WRITE_STRING_MATRIX_L")
    RETURN
999 CALL ERRORS("WRITE_STRING_MATRIX_L",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_MATRIX_L")
    RETURN 1
  END SUBROUTINE WRITE_STRING_MATRIX_L

  !
  !================================================================================================================================
  !

  !>Writes the given single precision MATRIX to the given output stream specified by ID. The basic output is determined by the flag INDEX_FORMAT_TYPE. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_ONLY then the first line of output for each row is MATRIX_NAME_FORMAT concatenated named with the FIRST_FORMAT. If INDEX_FORMAT_TYPE is WRITE_STRING_MATRIX_NAME_AND_INDICES then the first line of output for each row is MATRIX_NAME_FORMAT concatenated with ROW_INDEX_FORMAT and concatenated with FIRST_FORMAT. Note that with a WRITE_STRING_MATRIX_NAME_AND_INDICES index format type the row number will be supplied to the format before the matrix data. The FIRST_FORMAT is the format initially used, followed by the REPEAT_FORMAT which is repeated as many times as necessary. NUMBER_FIRST is the number of data items in the FIRST_FORMAT and NUMBER_REPEAT is the number of data items in the REPEAT_FORMAT. FIRST_ROW/FIRST_COLUMN and LAST_ROW/LAST_COLUMN are the extents of the row/column and DELTA_ROW/DELTA_COLUMN is the NUMBER of indices to skip for each row/column index.
  SUBROUTINE WRITE_STRING_MATRIX_SP(ID,FIRST_ROW,DELTA_ROW,LAST_ROW,FIRST_COLUMN,DELTA_COLUMN,LAST_COLUMN,NUMBER_FIRST, &
    & NUMBER_REPEAT,MATRIX,INDEX_FORMAT_TYPE,MATRIX_NAME_FORMAT,ROW_INDEX_FORMAT,FIRST_FORMAT,REPEAT_FORMAT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ID !<The ID of the output stream. An ID of > 9 specifies file output \see BASE_ROUTINES_OutputType,BASE_ROUTINES_FileUnits
    INTEGER(INTG), INTENT(IN) :: FIRST_ROW !<The first row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_ROW !<The delta row increment to be used when outputing the first through to the last matrix row
    INTEGER(INTG), INTENT(IN) :: LAST_ROW !<The last row of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: FIRST_COLUMN !<The first column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: DELTA_COLUMN !<The delta column increate to be used when outputing the first through to the last matrix column
    INTEGER(INTG), INTENT(IN) :: LAST_COLUMN !<The last column of the matrix to be output
    INTEGER(INTG), INTENT(IN) :: NUMBER_FIRST !<The number of matrix elements to be output on the first line
    INTEGER(INTG), INTENT(IN) :: NUMBER_REPEAT !<The number of matrix elements to be output on the second and subsequently repeated lines
    REAL(SP), INTENT(IN) :: MATRIX(:,:) !<The matrix to be output
    INTEGER(INTG), INTENT(IN) :: INDEX_FORMAT_TYPE !<The format type to be used for the matrix name and indices \see INPUT_OUTPUT_MatrixNameIndexFormat,INPUT_OUTPUT::MatrixNameIndexFormat
    CHARACTER(LEN=*), INTENT(IN) :: MATRIX_NAME_FORMAT !<The format string to be used to format the matrix name
    CHARACTER(LEN=*), INTENT(IN) :: ROW_INDEX_FORMAT !<The format string to be used to format the row indices
    CHARACTER(LEN=*), INTENT(IN) :: FIRST_FORMAT !<The format string to be used for the first line of output
    CHARACTER(LEN=*), INTENT(IN) :: REPEAT_FORMAT !<The format type to be used for the second and subsequently repeated lines of output
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) ::  current_row,current_column,final_column,count
    CHARACTER(LEN=MAXSTRLEN) :: FORMAT_STR

!    CALL ENTERS("WRITE_STRING_MATRIX_SP",ERR,ERROR,*999)

    IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//FIRST_FORMAT
    ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
      FORMAT_STR=MATRIX_NAME_FORMAT//ROW_INDEX_FORMAT//FIRST_FORMAT
    ELSE
      CALL FLAG_ERROR("Invalid index format type",ERR,ERROR,*999)
    ENDIF
    DO current_row=FIRST_ROW,LAST_ROW,DELTA_ROW
      current_column=FIRST_COLUMN
      final_column=current_column+(NUMBER_FIRST-1)*DELTA_COLUMN
      IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
      IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_ONLY) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ELSE IF(INDEX_FORMAT_TYPE==WRITE_STRING_MATRIX_NAME_AND_INDICES) THEN
        WRITE(OP_STRING,FMT=FORMAT_STR) current_row,(MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
      ENDIF
      CALL WRITE_STR(ID,ERR,ERROR,*999)
      DO WHILE(final_column<LAST_COLUMN) !more stuff to do
        current_column=final_column+DELTA_COLUMN
        final_column=final_column+NUMBER_REPEAT*DELTA_COLUMN
        IF(final_column>LAST_COLUMN) final_column=LAST_COLUMN
        WRITE(OP_STRING,FMT=REPEAT_FORMAT) (MATRIX(current_row,count),count=current_column,final_column,DELTA_COLUMN)
        CALL WRITE_STR(ID,ERR,ERROR,*999)
      ENDDO !final_columnn<LAST_COLUMN
    ENDDO !current_row
    
!    CALL EXITS("WRITE_STRING_MATRIX_SP")
    RETURN
999 CALL ERRORS("WRITE_STRING_MATRIX_SP",ERR,ERROR)
!    CALL EXITS("WRITE_STRING_MATRIX_SP")
    RETURN 1
  END SUBROUTINE WRITE_STRING_MATRIX_SP

  !
  !================================================================================================================================
  !

END MODULE INPUT_OUTPUT
