!> \file
!> \author Chris Bradley
!> \brief This module handles the reading and writing of binary files.
!> \todo Fix naming convention and update to current code standard.
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


!#### Module: BINARY_FILE 
!###  Description:
!###    This module handles the reading and writing of binary files.
!###  Uses: KINDS,CONSTANTS,MACHINE_CONSTANTS,BASE_ROUTINES,F90C,ISO_VARYING_STRING
!###  Routine: INQUIRE_OPEN_BINARY_FILE
!###  Routine: INQUIRE_EOF_BINARY_FILE
!###  Routine: CLOSE_BINARY_FILE
!###  Routine: CLOSE_CMISS_BINARY_FILE
!###  Routine: OPEN_BINARY_FILE
!###  Routine: OPEN_CMISS_BINARY_FILE
!###  Routine: READ_BINARY_FILE
!###  Routine: READ_BINARY_TAG_HEADER
!###  Routine: RESET_BINARY_NUMBER_TAGS
!###  Routine: SET_BINARY_FILE
!###  Routine: SKIP_CM_BINARY_HEADER
!###  Routine: SKIP_BINARY_FILE
!###  Routine: SKIP_BINARY_TAGS
!###  Routine: WRITE_BINARY_FILE
!###  Routine: WRITE_BINARY_TAG_HEADER

!#### Comment: BINARY FILE FORMAT
!###  Description:
!###    <HTML>
!###    The binary file format used for CMISS is as follows:
!###    Each cmiss binary file has the following header format:
!###    There are three sections in the header, namely:
!###    <UL>
!###    <LI>An identity header section.
!###    <LI>A machine header section.
!###    <LI>A file header section.
!###    </UL>
!###    The format of the identity header section is:
!###    <UL>
!###    <LI>1 identity byte (currently set to the value 0x07).
!###    <LI>1 byte to indicate the revision of the binary file format 
!###        (currently set to the value 0x02).
!###    </UL>
!###    The format of the machine header section is:
!###    <UL>
!###    <LI> If the the revision of the binary file format is 0x00
!###    <UL>
!###    <LI>2 dummy bytes (to be skipped).
!###    </UL>
!###    <LI> If the revision of the binary file format is 0x01
!###    <UL>
!###    <LI>1 byte to indicate the number (in standard base two) of 
!###        following bytes in the machine header section. This is 
!###        11 for this revision and the bytes are:
!###    <LI>1 byte to indicate the machine type which created the
!###        file (0x01 - DEC Alpha, 0x02 - SGI, 0x03 - IBM, 
!###        0x04 - Cray).
!###    <LI>1 byte to indicate the operating system type of the machine
!###        which created the file (redundant - for future development).
!###    <LI>1 byte to indicate the endian ordering of the numbers in
!###        in the file (0x01 - big endian, 0x02 - little endian).
!###        Currently all numbers are converted to big endian format.
!###    <LI>1 byte to indicate the format of standard single precision
!###        numbers (0x01 - IEEE single standard).
!###    <LI>1 byte to indicate the format of standard double precision
!###        numbers (0x01 - IEEE double standard). 
!###    <LI>1 byte to indicate the number of bytes a character type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a integer type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a short integer type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a single precision 
!###        type takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a double precision 
!###        type takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a logical 
!###        type takes up in the file.
!###    </UL>
!###    <LI> If the revision of the binary file format is 0x02
!###    <UL>
!###    </UL>
!###    <LI>1 byte to indicate the number (in standard base two) of 
!###        following bytes in the machine header section. This is 
!###        16 for this revision and the bytes are:
!###    <LI>1 byte to indicate the machine type which created the
!###        file (0x01 - DEC Alpha, 0x02 - SGI, 0x03 - IBM, 
!###        0x04 - Cray).
!###    <LI>1 byte to indicate the operating system type of the machine
!###        which created the file (redundant - for future development).
!###    <LI>1 byte to indicate the endian ordering of the numbers in
!###        in the file (0x01 - big endian, 0x02 - little endian).
!###    <LI>1 byte to indicate the format of characters (0x01 - 
!###        Ascii, 0x02 - Unicode).
!###    <LI>1 byte to indicate the format of integers (0x01 - 
!###        twos complement, 0x02 - Signed magnitude).
!###    <LI>1 byte to indicate the format of standard single precision
!###        numbers (0x01 - IEEE single standard).
!###    <LI>1 byte to indicate the format of standard double precision
!###        numbers (0x01 - IEEE double standard). 
!###    <LI>1 byte to indicate the number of bytes a character type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes an integer type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a short integer type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a long integer type
!###        takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a single precision 
!###        type takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a double precision 
!###        type takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a logical 
!###        type takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a single precision 
!###        complex type takes up in the file.
!###    <LI>1 byte to indicate the number of bytes a double precision 
!###        complex type takes up in the file.
!###    </UL>
!###    </UL>
!###    If the binary file format is below 0x02 then the format of the
!###    file header section is:
!###    <UL>
!###    <LI>Integer to specify the file type. Current binary 
!###        file types are: 1 - Binary matrix file, 2 - Binary time 
!###        series file, 3 - Binary signal file.
!###    <LI>Single precision number to specify the version of the file.
!###    <LI>Integer to specify the number of bytes in the heading.
!###    <LI>The heading (as a string of character bytes).
!###    <LI>Integer to specify how many 'tags' of data are in the file.
!###    </UL>
!###    The rest of the data in the file is made up of `tags' which
!###    contain the actual data. For each tag of data the format is:
!###    <UL>
!###    <LI>Integer to specify the type of tag.
!###    <LI>Integer to specify the number of bytes in the tag 
!###        heading.
!###    <LI>The tag heading (as a string of character bytes).
!###    <LI>Integer to specify the number of bytes in the tag 
!###        (excluding the tag header information).
!###    <LI>The tag data.
!###    </UL>
!###    If the binary file format is 0x02 then the format of the file
!###    header section is:
!###    <UL>
!###    <LI>Integer to specify the file type. Current binary 
!###        file types are: 1 - Binary matrix file, 2 - Binary time 
!###        series file, 3 - Binary signal file, 4 - Binary node file ??,
!###        5 - Binary element file ??
!###    <LI>Three integers to specify the version of the file in the
!###        form xxx.xxx.xxx
!###    <LI>Integer to specify the number of bytes in the heading.
!###    <LI>The heading (as a string of character bytes).
!###    <LI>Integer to specify how many tags are in the next level of
!###        the binary file.
!###    </UL>
!###    The rest of the data in the file is made up of a hierachy of
!###    `tags' which contain the actual data. For each tag of data the
!###    format is:
!###    <UL>
!###    <LI>Integer to specify the type of tag.
!###    <LI>Integer to specify the number of bytes in the tag 
!###        heading.
!###    <LI>The tag heading (as a string of character bytes).
!###    <LI>An integer to specify the number of tags below this tag.
!###        If this number is > 0 the tag is known as a node tag.
!###        If this number is = 0 the tag is known as a leaf tag.
!###    <LI>If the tag is a leaf tag
!###    <UL>
!###    <LI>Integer to specify the number of bytes in the tag 
!###        (excluding the tag header information). NOTE: This restricts
!###        the amount of information within a tag to under 2 GB. Future
!###        revisions of the binary file format will change this quantity
!###        to a long integer (i.e. 64-bits). This will have to wait, however,
!###        until intel based machines can handle 64-bit integers.
!###    <LI>The tag data.
!###    </UL>
!###    </UL>
!###    Note: that if any of the byte numbers are unknown they will have
!###    the value of 0xFF.
!###    </HTML>

!>This module handles the reading and writing of binary files.
MODULE BINARY_FILE

  USE KINDS
  USE CONSTANTS
  USE MACHINE_CONSTANTS
  USE BASE_ROUTINES
  USE F90C
  USE ISO_VARYING_STRING
  
#include "macros.h"

  IMPLICIT NONE

  !PRIVATE

  !Module parameters
  
  INTEGER(INTG), PARAMETER :: MAX_NUM_BINARY_FILES=99

  !Skip file parameters
  INTEGER(INTG), PARAMETER :: FILE_BEGINNING=0
  INTEGER(INTG), PARAMETER :: FILE_CURRENT=1
  INTEGER(INTG), PARAMETER :: FILE_END=2

  !File endian parameters
  INTEGER(INTG), PARAMETER :: FILE_SAME_ENDIAN=0
  INTEGER(INTG), PARAMETER :: FILE_CHANGE_ENDIAN=1
  
  !Binary file parameters
  INTEGER(INTG), PARAMETER :: BINARY_FILE_READABLE=1
  INTEGER(INTG), PARAMETER :: BINARY_FILE_WRITABLE=2

  !CMISS Binary file parameters
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_IDENTITY=7
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_MATRIX_FILE=1
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_HISTORY_FILE=2
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_SIGNAL_FILE=3
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_IDENTITY_HEADER=1
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_MACHINE_HEADER=2
  INTEGER(INTG), PARAMETER :: CMISS_BINARY_FILE_HEADER=3

  !Module types

  TYPE BINARY_FILE_INFO_TYPE
    PRIVATE
    INTEGER(INTG) :: FILE_NUMBER
    INTEGER(INTG) :: BINARY_FILE_REVISION
    INTEGER(INTG) :: MACHINE_TYPE
    INTEGER(INTG) :: OS_TYPE
    INTEGER(INTG) :: ENDIAN_TYPE
    INTEGER(INTG) :: CHAR_FORMAT
    INTEGER(INTG) :: INT_FORMAT
    INTEGER(INTG) :: SP_FORMAT
    INTEGER(INTG) :: DP_FORMAT
    INTEGER(INTG) :: CHARACTER_SIZE
    INTEGER(INTG) :: INTEGER_SIZE
    INTEGER(INTG) :: SINTEGER_SIZE
    INTEGER(INTG) :: LINTEGER_SIZE
    INTEGER(INTG) :: SP_REAL_SIZE
    INTEGER(INTG) :: DP_REAL_SIZE
    INTEGER(INTG) :: LOGICAL_SIZE
    INTEGER(INTG) :: SPC_REAL_SIZE
    INTEGER(INTG) :: DPC_REAL_SIZE
    CHARACTER(LEN=MAXSTRLEN) :: FILE_NAME
    INTEGER(INTG) :: ACCESS_TYPE
  END TYPE BINARY_FILE_INFO_TYPE

  TYPE BINARY_FILE_TYPE
    TYPE(BINARY_FILE_INFO_TYPE), POINTER :: FILE_INFORMATION
  END TYPE BINARY_FILE_TYPE

  TYPE BINARY_TAG_TYPE
    INTEGER(INTG) :: INDEX
    INTEGER(INTG) :: NUM_SUBTAGS
    INTEGER(INTG) :: NUM_BYTES
    INTEGER(INTG) :: NUM_HEADER_BYTES
    CHARACTER(LEN=MAXSTRLEN) :: HEADER
  END TYPE BINARY_TAG_TYPE

  !Module variables

  LOGICAL, SAVE :: BINARY_FILE_USED(MAX_NUM_BINARY_FILES)=.FALSE.

  !Interfaces

  INTERFACE
    
    SUBROUTINE BINARYCLOSEFILE(FILE_NUMBER,ERR,CERROR)
!DEC$ ATTRIBUTES C :: binaryclosefile
      USE CONSTANTS
      INTEGER(INTG), INTENT(IN) :: FILE_NUMBER
      INTEGER(INTG), INTENT(OUT) :: ERR,CERROR(*)
    END SUBROUTINE BINARYCLOSEFILE

    SUBROUTINE BINARYOPENFILE(FILE_NUMBER,CFNAME,CACCESSCODE,ERR,CERROR)
!DEC$ ATTRIBUTES C :: binaryopenfile
      USE CONSTANTS
      INTEGER(INTG), INTENT(IN) :: FILE_NUMBER,CFNAME(*),CACCESSCODE(*)
      INTEGER(INTG), INTENT(OUT) :: ERR,CERROR(*)
    END SUBROUTINE BINARYOPENFILE

!!DEC$ ATTRIBUTES C :: binaryreadfile

    SUBROUTINE BINARYSETFILE(FILE_NUMBER,SET_CODE,ERR,CERROR)
      USE CONSTANTS
!DEC$ ATTRIBUTES C :: binarysetfile
      INTEGER(INTG), INTENT(IN) :: FILE_NUMBER,SET_CODE
      INTEGER(INTG), INTENT(OUT) :: ERR,CERROR(*)
    END SUBROUTINE BINARYSETFILE

    SUBROUTINE BINARYSKIPFILE(FILE_NUMBER,NUMBER_BYTES,ERR,CERROR)
!DEC$ ATTRIBUTES C :: binaryskipfile
     USE CONSTANTS
      INTEGER(INTG), INTENT(IN) :: FILE_NUMBER,NUMBER_BYTES
      INTEGER(INTG), INTENT(OUT) :: ERR,CERROR(*)
    END SUBROUTINE BINARYSKIPFILE

!!    SUBROUTINE BINARYWRITEFILE
!!!DEC$ ATTRIBUTES C :: binarywritefile
!!    END SUBROUTINE BINARYWRITEFILE

    SUBROUTINE ISBINARYFILEOPEN(FILE_NUMBER,RETURNCODE,ERR,CERROR)
!DEC$ ATTRIBUTES C :: isbinaryfileopen
      USE CONSTANTS
      INTEGER(INTG), INTENT(IN) :: FILE_NUMBER
      INTEGER(INTG), INTENT(OUT) :: RETURNCODE, ERR, CERROR(*)
    END SUBROUTINE ISBINARYFILEOPEN

    SUBROUTINE ISENDBINARYFILE(FILE_NUMBER,RETURN_CODE,ERR,CERROR)
!DEC$ ATTRIBUTES C :: isbinaryfileopen
      USE CONSTANTS
      INTEGER(INTG), INTENT(IN) :: FILE_NUMBER
      INTEGER(INTG), INTENT(OUT) :: RETURN_CODE,ERR,CERROR(*)
    END SUBROUTINE ISENDBINARYFILE

  END INTERFACE  

  INTERFACE READ_BINARY_FILE
    MODULE PROCEDURE READ_BINARY_FILE_INTG
    MODULE PROCEDURE READ_BINARY_FILE_INTG1
    MODULE PROCEDURE READ_BINARY_FILE_SINTG
    MODULE PROCEDURE READ_BINARY_FILE_SINTG1
    MODULE PROCEDURE READ_BINARY_FILE_LINTG
    MODULE PROCEDURE READ_BINARY_FILE_LINTG1
    MODULE PROCEDURE READ_BINARY_FILE_SP
    MODULE PROCEDURE READ_BINARY_FILE_SP1
    MODULE PROCEDURE READ_BINARY_FILE_DP
    MODULE PROCEDURE READ_BINARY_FILE_DP1
    MODULE PROCEDURE READ_BINARY_FILE_CHARACTER
    MODULE PROCEDURE READ_BINARY_FILE_LOGICAL
    MODULE PROCEDURE READ_BINARY_FILE_LOGICAL1
    MODULE PROCEDURE READ_BINARY_FILE_SPC
    MODULE PROCEDURE READ_BINARY_FILE_SPC1
    MODULE PROCEDURE READ_BINARY_FILE_DPC
    MODULE PROCEDURE READ_BINARY_FILE_DPC1
  END INTERFACE !READ_BINARY_FILE

  INTERFACE WRITE_BINARY_FILE
    MODULE PROCEDURE WRITE_BINARY_FILE_INTG
    MODULE PROCEDURE WRITE_BINARY_FILE_INTG1
    MODULE PROCEDURE WRITE_BINARY_FILE_SINTG
    MODULE PROCEDURE WRITE_BINARY_FILE_SINTG1
    MODULE PROCEDURE WRITE_BINARY_FILE_LINTG
    MODULE PROCEDURE WRITE_BINARY_FILE_LINTG1
    MODULE PROCEDURE WRITE_BINARY_FILE_SP
    MODULE PROCEDURE WRITE_BINARY_FILE_SP1
    MODULE PROCEDURE WRITE_BINARY_FILE_DP
    MODULE PROCEDURE WRITE_BINARY_FILE_DP1
    MODULE PROCEDURE WRITE_BINARY_FILE_CHARACTER
    MODULE PROCEDURE WRITE_BINARY_FILE_LOGICAL
    MODULE PROCEDURE WRITE_BINARY_FILE_LOGICAL1
    MODULE PROCEDURE WRITE_BINARY_FILE_SPC
    MODULE PROCEDURE WRITE_BINARY_FILE_SPC1
    MODULE PROCEDURE WRITE_BINARY_FILE_DPC
    MODULE PROCEDURE WRITE_BINARY_FILE_DPC1
  END INTERFACE !WRITE_BINARY_FILE

  PUBLIC FILE_BEGINNING,FILE_CURRENT,FILE_END,CMISS_BINARY_MATRIX_FILE,&
    & CMISS_BINARY_HISTORY_FILE,CMISS_BINARY_SIGNAL_FILE,BINARY_FILE_TYPE
  PUBLIC READ_BINARY_FILE,WRITE_BINARY_FILE

CONTAINS

  !
  !============================================================================
  !
  
  FUNCTION INQUIRE_OPEN_BINARY_FILE(FILEID)
    
    !#### Function: INQUIRE_OPEN_BINARY_FILE
    !###  Type: LOGICAL
    !###  Desctiption:
    !###    INQUIRE_OPEN_BINARY_FILE returns .TRUE. if the binary file
    !###    specified by FILEID is OPEN, .FALSE. if not.
    
    !Arguments
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    !Function Variable
    LOGICAL :: INQUIRE_OPEN_BINARY_FILE
    !Local Variables

    INQUIRE_OPEN_BINARY_FILE=ASSOCIATED(FILEID%FILE_INFORMATION)
      
    RETURN
  END FUNCTION INQUIRE_OPEN_BINARY_FILE

  !
  !============================================================================
  !
  
  FUNCTION INQUIRE_EOF_BINARY_FILE(FILEID, ERR, ERROR)
    
    !#### Function: INQUIRE_EOF_BINARY_FILE
    !###  Type: LOGICAL
    !###  Description:
    !###    INQUIRE_EOF_BINARY_FILE returns .TRUE. If the binary file 
    !###    specified by FILEID is at eof, .FALSE. if not.
    
    !Arguments
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function Variable
    LOGICAL :: INQUIRE_EOF_BINARY_FILE
    !Local Variables
    INTEGER(INTG) :: CERROR(100),RETURNCODE
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("INQUIRE_EOF_BINARY_FILE",ERR,ERROR,*999)

    INQUIRE_EOF_BINARY_FILE=.FALSE.
    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL ISENDBINARYFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & RETURNCODE,ERR,CERROR)    
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
      INQUIRE_EOF_BINARY_FILE=(RETURNCODE==1)
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("INQUIRE_EOF_BINARY_FILE")
    RETURN
999 ERRORSEXITS("INQUIRE_EOF_BINARY_FILE",ERR,ERROR)
    RETURN
  END FUNCTION 

  !
  !============================================================================
  !
  
  SUBROUTINE CLOSE_BINARY_FILE(FILEID,ERR,ERROR,*)
    
    !#### Subroutine: CLOSE_BINARY_FILE
    !###  Description:
    !###    CLOSE_BINARY_FILE closes the binary file specified by
    !###    FILEID and deallocates the binary file information.
    
    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(OUT) :: FILEID
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR
    
    ENTERS("CLOSE_BINARY_FILE",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL BINARYCLOSEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER, ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
      !Could be a problem with the file not being closed properly.
      !DEALLOCATE the file information for now and assume that the
      !file has already been closed.
      BINARY_FILE_USED(FILEID%FILE_INFORMATION%FILE_NUMBER)=.FALSE.
      DEALLOCATE(FILEID%FILE_INFORMATION)
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("CLOSE_BINARY_FILE")
    RETURN
999 ERRORSEXITS("CLOSE_BINARY_FILE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLOSE_BINARY_FILE

  !
  !============================================================================
  !
  
  SUBROUTINE CLOSE_CMISS_BINARY_FILE(FILEID,ERR,ERROR,*)
    
    !#### Subroutine: CLOSE_CMISS_BINARY_FILE
    !###  Description:
    !###    CLOSE_CMISS_BINARY_FILE closes the CMISS binary file specified by
    !###    FILEID.
    
    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(INOUT) :: FILEID
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: CERROR(100)
    
    ENTERS("CLOSE_CMISS_BINARY_FILE",ERR,ERROR,*999)

    CALL CLOSE_BINARY_FILE(FILEID,ERR,ERROR,*999)

    EXITS("CLOSE_CMISS_BINARY_FILE")
    RETURN
999 ERRORSEXITS("CLOSE_CMISS_BINARY_FILE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CLOSE_CMISS_BINARY_FILE

  !
  !============================================================================
  !
  
  SUBROUTINE OPEN_BINARY_FILE(FILEID,COMMAND,FILENAME,ERR,ERROR,*)
    
    !#### Subroutine: OPEN_BINARY_FILE
    !###  Description:
    !###    OPEN_BINARY_FILE opens the binary file specified by
    !###    FILEID with the given FILENAME. The file will be opened
    !###    for reading if COMMAND is "READ" and writting if COMMAND
    !###    is "WRITE".
        
    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(OUT) :: FILEID
    CHARACTER(LEN=*), INTENT(IN) :: COMMAND, FILENAME
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: CACCESSCODE(2), CERROR(100), CFNAME(100),&
      & FILENUMBER,FILENUM
    CHARACTER(LEN=6) :: FACCESSCODE
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("OPEN_BINARY_FILE",ERR,ERROR,*999)

    FILENUMBER=0
    DO FILENUM=1,MAX_NUM_BINARY_FILES
      IF(.NOT.BINARY_FILE_USED(FILENUM)) THEN
        FILENUMBER=FILENUM
        EXIT
      ENDIF
    ENDDO
    IF(FILENUMBER/=0) THEN
      ALLOCATE(FILEID%FILE_INFORMATION,STAT=ERR)
      IF(ERR==0) THEN
        BINARY_FILE_USED(FILEID%FILE_INFORMATION%FILE_NUMBER)=.TRUE.
        FILEID%FILE_INFORMATION%FILE_NUMBER=FILENUMBER
        FILEID%FILE_INFORMATION%FILE_NAME=FILENAME
        FILEID%FILE_INFORMATION%ENDIAN_TYPE=MACHINE_ENDIAN
        CALL F2CSTRING(CFNAME,FILENAME,ERR,ERROR,*999)        
        IF(COMMAND(1:4)=="READ") THEN
          FACCESSCODE="rb+"
          FILEID%FILE_INFORMATION%ACCESS_TYPE=BINARY_FILE_READABLE
        ELSE IF(COMMAND(1:5)=="WRITE") THEN
          FACCESSCODE="wb+"
          FILEID%FILE_INFORMATION%ACCESS_TYPE=BINARY_FILE_WRITABLE
        ELSE
          CALL FLAG_ERROR("Invalid command",ERR,ERROR,*999)
        ENDIF
        CALL F2CSTRING(CACCESSCODE,FACCESSCODE,ERR,ERROR,*999)
        CALL BINARYOPENFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & CFNAME,CACCESSCODE,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Could not allocate binary file information",&
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("No free binary files available",ERR,ERROR,*999)
    ENDIF
    
    EXITS("OPEN_BINARY_FILE")
    RETURN
999 IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      BINARY_FILE_USED(FILEID%FILE_INFORMATION%FILE_NUMBER)=.FALSE.
      DEALLOCATE(FILEID%FILE_INFORMATION)
    ENDIF
    ERRORSEXITS("OPEN_BINARY_FILE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE OPEN_BINARY_FILE
  
  !
  !============================================================================
  !
  
  SUBROUTINE OPEN_CMISS_BINARY_FILE(FILEID,FILE_TYPE,NUMBER_TAGS,&
    & VERSION,FILEVERSION,COMMAND,EXTENSION,FILENAME,ERR,ERROR,*)
    
    !#### Subroutine: OPEN_CMISS_BINARY_FILE
    !###  Description:
    !###    OPEN_CMISS_BINARY_FILE opens a CMISS binary file called
    !###    file.bin*** (where *** is the extension) specified by
    !###    FILEID, allocates the file information and reads/writes
    !###    (specified by command) file type, version number, heading
    !###    and number of tags.
    
    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(OUT) :: FILEID
    INTEGER(INTG), INTENT(IN) :: FILE_TYPE
    INTEGER(INTG), INTENT(INOUT) :: NUMBER_TAGS
    INTEGER(INTG), INTENT(IN) :: VERSION(3)
    INTEGER(INTG), INTENT(OUT) :: FILEVERSION(3)
    CHARACTER(LEN=*), INTENT(IN) :: COMMAND, EXTENSION, FILENAME
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: FILE_NUM, FILE_NUMBER, FILETYPE, HEADINGSIZE,&
      & NUMMACHHEADERBYTES
    REAL(SP) :: FVERSION
    CHARACTER(LEN=17) :: CHARDATA
    CHARACTER(LEN=MAXSTRLEN) :: ERROR_STRING,FULL_FILENAME,HEADING,&
      & WARNING_STRING
    
    ENTERS("OPEN_CMISS_BINARY_FILE",ERR,ERROR,*999)

    FULL_FILENAME=FILENAME(1:LEN_TRIM(FILENAME))//".bin"&
      & //EXTENSION(1:LEN_TRIM(EXTENSION))
    CALL OPEN_BINARY_FILE(FILEID, COMMAND,FULL_FILENAME,ERR,ERROR,*999)
    IF(COMMAND(1:4)=="READ") THEN
      !Read identity header section
      CALL READ_BINARY_FILE(FILEID,2,CHARDATA,ERR,ERROR,*999)
      IF(CHARDATA(1:1)/=CHAR(CMISS_BINARY_IDENTITY)) &
        & CALL FLAG_ERROR("Not a CMISS binary file",ERR,ERROR,*999)
      FILEID%FILE_INFORMATION%BINARY_FILE_REVISION=ICHAR(CHARDATA(2:2))
      !Read machine header section
      SELECT CASE(FILEID%FILE_INFORMATION%BINARY_FILE_REVISION)
      CASE(0)
        !Identity format 0 - has 2 extra dummy bytes. Skip them.
        CALL FLAG_WARNING("Old binary file identity format&
          & found. Please update file.",ERR,ERROR,*999)
        CALL SKIP_BINARY_FILE(FILEID,2*CHARACTER_SIZE,ERR,&
          & ERROR,*999)
      CASE(1)
        !Identity format 1 - has machine header section.
        CALL FLAG_WARNING("Old binary file identity format&
          & found. Please update file.",ERR,ERROR,*999)
        CALL READ_BINARY_FILE(FILEID,1,CHARDATA,ERR,ERROR,*999)
        NUMMACHHEADERBYTES=ICHAR(CHARDATA(1:1))
        CALL READ_BINARY_FILE(FILEID,NUMMACHHEADERBYTES,CHARDATA,&
          & ERR,ERROR,*999)
        IF(NUMMACHHEADERBYTES==11) THEN
          FILEID%FILE_INFORMATION%MACHINE_TYPE=ICHAR(CHARDATA(1:1))
          FILEID%FILE_INFORMATION%OS_TYPE=ICHAR(CHARDATA(2:2))
          FILEID%FILE_INFORMATION%ENDIAN_TYPE=ICHAR(CHARDATA(3:3))
          FILEID%FILE_INFORMATION%SP_FORMAT=ICHAR(CHARDATA(4:4))
          FILEID%FILE_INFORMATION%DP_FORMAT=ICHAR(CHARDATA(5:5))
          FILEID%FILE_INFORMATION%CHARACTER_SIZE=ICHAR(CHARDATA(6:6))
          FILEID%FILE_INFORMATION%INTEGER_SIZE=ICHAR(CHARDATA(7:7))
          FILEID%FILE_INFORMATION%SINTEGER_SIZE=ICHAR(CHARDATA(8:8))
          FILEID%FILE_INFORMATION%SP_REAL_SIZE=ICHAR(CHARDATA(9:9))
          FILEID%FILE_INFORMATION%DP_REAL_SIZE=ICHAR(CHARDATA(10:10))
          FILEID%FILE_INFORMATION%LOGICAL_SIZE=ICHAR(CHARDATA(11:11))
        ELSE
          CALL FLAG_ERROR("Invalid number of machine header bytes",&
            & ERR,ERROR,*999)
        ENDIF
      CASE(2)
        !Identity format 2 - has extended machine header section and subtags.
        CALL READ_BINARY_FILE(FILEID,1,CHARDATA,ERR,ERROR,*999)
        NUMMACHHEADERBYTES=ICHAR(CHARDATA(1:1))
        CALL READ_BINARY_FILE(FILEID,NUMMACHHEADERBYTES,CHARDATA,&
          & ERR,ERROR,*999)
        IF(NUMMACHHEADERBYTES==16) THEN
          FILEID%FILE_INFORMATION%MACHINE_TYPE=ICHAR(CHARDATA(1:1))
          FILEID%FILE_INFORMATION%OS_TYPE=ICHAR(CHARDATA(2:2))
          FILEID%FILE_INFORMATION%ENDIAN_TYPE=ICHAR(CHARDATA(3:3))
          FILEID%FILE_INFORMATION%CHAR_FORMAT=ICHAR(CHARDATA(4:4))
          FILEID%FILE_INFORMATION%INT_FORMAT=ICHAR(CHARDATA(5:5))
          FILEID%FILE_INFORMATION%SP_FORMAT=ICHAR(CHARDATA(6:6))
          FILEID%FILE_INFORMATION%DP_FORMAT=ICHAR(CHARDATA(7:7))
          FILEID%FILE_INFORMATION%CHARACTER_SIZE=ICHAR(CHARDATA(8:8))
          FILEID%FILE_INFORMATION%INTEGER_SIZE=ICHAR(CHARDATA(9:9))
          FILEID%FILE_INFORMATION%SINTEGER_SIZE=ICHAR(CHARDATA(10:10))
          FILEID%FILE_INFORMATION%LINTEGER_SIZE=ICHAR(CHARDATA(11:11))
          FILEID%FILE_INFORMATION%SP_REAL_SIZE=ICHAR(CHARDATA(12:12))
          FILEID%FILE_INFORMATION%DP_REAL_SIZE=ICHAR(CHARDATA(13:13))
          FILEID%FILE_INFORMATION%LOGICAL_SIZE=ICHAR(CHARDATA(14:14))
          FILEID%FILE_INFORMATION%SPC_REAL_SIZE=ICHAR(CHARDATA(15:15))
          FILEID%FILE_INFORMATION%DPC_REAL_SIZE=ICHAR(CHARDATA(16:16))
        ELSE
          CALL FLAG_ERROR("Invalid number of machine header bytes",&
            & ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        CALL FLAG_ERROR("Unknown binary file header identity format",&
          & ERR,ERROR,*999)
      END SELECT
      !Read file header section
      !Read file type, version and header
      CALL READ_BINARY_FILE(FILEID,1,FILETYPE,ERR,ERROR,*999)
      IF(FILETYPE/=FILE_TYPE) THEN
        WRITE(ERROR_STRING,'("File has a different file type:",/,&
          & "   File type is     ",I3,/,"   Expected file type is ",I3)')&
          & FILETYPE,FILE_TYPE
        CALL FLAG_ERROR(WARNING_STRING,ERR,ERROR,*999)
      ENDIF
      SELECT CASE(FILEID%FILE_INFORMATION%BINARY_FILE_REVISION)
      CASE(0,1)
        CALL READ_BINARY_FILE(FILEID,1,FVERSION,ERR,ERROR,*999)
        FILEVERSION(1)=INT(FVERSION,INTG)
        FILEVERSION(2)=0
        FILEVERSION(3)=0
        IF(FILEVERSION(1)/=VERSION(1)) THEN
          WRITE(WARNING_STRING,'("File has a different version number:",/,&
            & "   File version is          ",I2,".",I2,".",I2,/,&
            & "   Expected file version is ",I2,".",I2,".",I2)') &
            & FILEVERSION(1),FILEVERSION(2),FILEVERSION(3),&
            & VERSION(1),VERSION(2),VERSION(3)
          CALL FLAG_WARNING(WARNING_STRING,ERR,ERROR,*999)
        ENDIF
      CASE(2)
        CALL READ_BINARY_FILE(FILEID,3,FILEVERSION,ERR,ERROR,*999)
        IF(FILEVERSION(1)/=VERSION(1).OR.&
          & FILEVERSION(2)/=VERSION(2).OR.&
          & FILEVERSION(3)/=VERSION(3)) THEN
          WRITE(WARNING_STRING,'("File has a different version number:",/,&
            & "   File version is          ",I2,".",I2,".",I2,/,&
            & "   Expected file version is ",I2,".",I2,".",I2)') &
            & FILEVERSION(1),FILEVERSION(2),FILEVERSION(3),&
            & VERSION(1),VERSION(2),VERSION(3)
          CALL FLAG_WARNING(WARNING_STRING,ERR,ERROR,*999)
        ENDIF
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid binary file identity format",ERR,ERROR,*999)
      END SELECT
      CALL READ_BINARY_FILE(FILEID,1,HEADINGSIZE,ERR,ERROR,*999)
      IF(HEADINGSIZE==0) THEN
        WRITE(OP_STRING,'("File heading: ")')
      ELSE
        CALL READ_BINARY_FILE(FILEID,HEADINGSIZE,HEADING,ERR,ERROR,*999)
        WRITE(OP_STRING,'("File heading: ",A)') HEADING(1:HEADINGSIZE)
      ENDIF
      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,ERR,ERROR,*999)
      CALL READ_BINARY_FILE(FILEID,1,NUMBER_TAGS,ERR,ERROR,*999)
    ELSE IF(COMMAND(1:5)=="WRITE") THEN
      !Write identity header section
      CHARDATA(1:1)=CHAR(CMISS_BINARY_IDENTITY)
      CHARDATA(2:2)=CHAR(2)
      CALL WRITE_BINARY_FILE(FILEID,2,CHARDATA,ERR,ERROR,*999)
      !Write the machine header section
      CHARDATA(1:1)=CHAR(16)
      CHARDATA(2:2)=CHAR(MACHINE_TYPE)
      CHARDATA(3:3)=CHAR(MACHINE_OS)
      CHARDATA(4:4)=CHAR(MACHINE_ENDIAN)
      CHARDATA(5:5)=CHAR(MACHINE_CHAR_FORMAT)
      CHARDATA(6:6)=CHAR(MACHINE_INT_FORMAT)
      CHARDATA(7:7)=CHAR(MACHINE_SP_FORMAT)
      CHARDATA(8:8)=CHAR(MACHINE_DP_FORMAT)
      CHARDATA(9:9)=CHAR(CHARACTER_SIZE)
      CHARDATA(10:10)=CHAR(INTEGER_SIZE)
      CHARDATA(11:11)=CHAR(SHORT_INTEGER_SIZE)
      CHARDATA(12:12)=CHAR(LONG_INTEGER_SIZE)
      CHARDATA(13:13)=CHAR(SINGLE_REAL_SIZE)
      CHARDATA(14:14)=CHAR(DOUBLE_REAL_SIZE)
      CHARDATA(15:15)=CHAR(LOGICAL_SIZE)
      CHARDATA(16:16)=CHAR(SINGLE_COMPLEX_SIZE)
      CHARDATA(17:17)=CHAR(DOUBLE_COMPLEX_SIZE)
      FILEID%FILE_INFORMATION%MACHINE_TYPE=MACHINE_TYPE
      FILEID%FILE_INFORMATION%OS_TYPE=MACHINE_OS
      FILEID%FILE_INFORMATION%ENDIAN_TYPE=MACHINE_ENDIAN
      FILEID%FILE_INFORMATION%CHAR_FORMAT=MACHINE_CHAR_FORMAT
      FILEID%FILE_INFORMATION%INT_FORMAT=MACHINE_INT_FORMAT
      FILEID%FILE_INFORMATION%SP_FORMAT=MACHINE_SP_FORMAT
      FILEID%FILE_INFORMATION%DP_FORMAT=MACHINE_DP_FORMAT
      FILEID%FILE_INFORMATION%CHARACTER_SIZE=CHARACTER_SIZE
      FILEID%FILE_INFORMATION%INTEGER_SIZE=INTEGER_SIZE
      FILEID%FILE_INFORMATION%SINTEGER_SIZE=SHORT_INTEGER_SIZE
      FILEID%FILE_INFORMATION%LINTEGER_SIZE=LONG_INTEGER_SIZE
      FILEID%FILE_INFORMATION%SP_REAL_SIZE=SINGLE_REAL_SIZE
      FILEID%FILE_INFORMATION%DP_REAL_SIZE=DOUBLE_REAL_SIZE
      FILEID%FILE_INFORMATION%LOGICAL_SIZE=LOGICAL_SIZE
      CALL WRITE_BINARY_FILE(FILEID,12,CHARDATA,ERR,ERROR,*999)
      !Write the file header section
      !Write file type, version and header
      CALL WRITE_BINARY_FILE(FILEID,1,FILE_TYPE,ERR,ERROR,*999)
      CALL WRITE_BINARY_FILE(FILEID,1,VERSION,ERR,ERROR,*999)
      HEADING=" "
      HEADINGSIZE=LEN_TRIM(HEADING)
      CALL WRITE_BINARY_FILE(FILEID,1,HEADINGSIZE,ERR,ERROR,*999)
      CALL WRITE_BINARY_FILE(FILEID,HEADINGSIZE,HEADING,ERR,ERROR,*999)
      !Write number of tags
      CALL WRITE_BINARY_FILE(FILEID,1,NUMBER_TAGS,ERR,ERROR,*999)
    ENDIF

    EXITS("OPEN_CMISS_BINARY_FILE")
    RETURN
999 ERRORSEXITS("OPEN_CMISS_BINARY_FILE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE OPEN_CMISS_BINARY_FILE

  !
  !============================================================================
  !
  
  !#### Generic-Subroutine: READ_BINARY_FILE
  !###  Description:
  !###    READ_BINARY_FILE reads NUM_DATA elements from the binary
  !###    file specified by FILEID into DATA.
  !###  Child-Subroutines: READ_BINARY_FILE_INTG,READ_BINARY_FILE_INTG1,
  !###    READ_BINARY_FILE_SINTG,READ_BINARY_FILE_SINTG1,
  !###    READ_BINARY_FILE_LINTG,READ_BINARY_FILE_LINTG1,
  !###    READ_BINARY_FILE_SP,READ_BINARY_FILE_SP1,
  !###    READ_BINARY_FILE_DP,READ_BINARY_FILE_DP1,
  !###    READ_BINARY_FILE_CHARACTER,READ_BINARY_FILE_LOGICAL,
  !###    READ_BINARY_FILE_LOGICAL1,READ_BINARY_FILE_SPC,
  !###    READ_BINARY_FILE_SPC1,READ_BINARY_FILE_DPC,
  !###    READ_BINARY_FILE_DPC1
  
  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_INTG(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_INTG
    !###  Description:
    !###    READ_BINARY_FILE_INTG reads NUM_DATA integer values from
    !###    the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(INTG), INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        ENDIAN,NUM_DATA,INTEGER_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("READ_BINARY_FILE_INTG")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_INTG",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_INTG
  
  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_INTG1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_INTG1
    !###  Description:
    !###    READ_BINARY_FILE_INTG1 reads 1 integer value from
    !###    the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(INTG), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,INTEGER_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_INTG1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_INTG1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_INTG1
  
  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_SINTG(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_SINTG
    !###  Description:
    !###    READ_BINARY_FILE_SINTG reads NUM_DATA short integer values
    !###    from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(SINTG), INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_SINTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,SHORT_INTEGER_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_SINTG")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_SINTG",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_SINTG

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_SINTG1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_SINTG1
    !###  Description:
    !###    READ_BINARY_FILE_SINTG1 reads 1 short integer value
    !###    from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(SINTG), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_SINTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,SHORT_INTEGER_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("READ_BINARY_FILE_SINTG1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_SINTG1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_SINTG1

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_LINTG(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_LINTG
    !###  Description:
    !###    READ_BINARY_FILE_LINTG reads NUM_DATA long integer values
    !###    from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(LINTG), INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_LINTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,LONG_INTEGER_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_LINTG")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_LINTG",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_LINTG

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_LINTG1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_LINTG1
    !###  Description:
    !###    READ_BINARY_FILE_LINTG1 reads 1 long integer value
    !###    from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(LINTG), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_LINTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,LONG_INTEGER_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("READ_BINARY_FILE_LINTG1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_LINTG1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_LINTG1

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_SP(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_SP
    !###  Description:
    !###    READ_BINARY_FILE_SP reads NUM_DATA single precision real
    !###    values from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(SP), INTENT(OUT) :: DATA(*)
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    INTEGER(INTG), INTENT(OUT) :: ERR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,SINGLE_REAL_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("READ_BINARY_FILE_SP")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_SP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_SP

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_SP1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_SP1
    !###  Description:
    !###    READ_BINARY_FILE_SP1 reads 1 single precision real
    !###    value from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(SP), INTENT(OUT) :: DATA
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    INTEGER(INTG), INTENT(OUT) :: ERR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,SINGLE_REAL_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("READ_BINARY_FILE_SP1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_SP1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_SP1

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_DP(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_DP
    !###  Description:
    !###    READ_BINARY_FILE_DP reads NUM_DATA double precision real
    !###    values from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE
 
    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(DP), INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,DOUBLE_REAL_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
      
    EXITS("READ_BINARY_FILE_DP")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_DP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_DP

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_DP1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_DP1
    !###  Description:
    !###    READ_BINARY_FILE_DP1 reads 1 double precision real
    !###    value from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(DP), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,DOUBLE_REAL_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
      
    EXITS("READ_BINARY_FILE_DP1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_DP1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_DP1

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_CHARACTER(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_CHARACTER
    !###  Description:
    !###    READ_BINARY_FILE_CHARACTER reads NUM_DATA character
    !###    values from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    CHARACTER(LEN=*), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),CSTRING(250),LENGTH
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_CHARACTER",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER, &
        & FILE_SAME_ENDIAN,NUM_DATA,CHARACTER_TYPE,CSTRING,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ELSE
        LENGTH=CSTRINGLENGTH(CSTRING)
        IF(LENGTH==0.AND.NUM_DATA==1) THEN
          !Zero length string and data read i.e. read a single byte
          !of value zero.
          DATA(1:1)=CHAR(0)
        ELSE
          CALL C2FSTRING(CSTRING,DATA,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_CHARACTER")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_CHARACTER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_CHARACTER
  
  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_LOGICAL(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_LOGICAL
    !###  Description:
    !###    READ_BINARY_FILE_LOGICAL reads NUM_DATA logical
    !###    values from the binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    LOGICAL, INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_LOGICAL",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,LOGICAL_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_LOGICAL")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_LOGICAL",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_LOGICAL

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_LOGICAL1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_LOGICAL1
    !###  Description:
    !###    READ_BINARY_FILE_LOGICAL1 reads 1 logical value from the
    !###    binary file specified by FILEID into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    LOGICAL, INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_LOGICAL1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,LOGICAL_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_LOGICAL1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_LOGICAL1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_LOGICAL1

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_SPC(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_SPC
    !###  Description:
    !###    READ_BINARY_FILE_SPC reads NUM_DATA single precision
    !###    complex values from the binary file specified by FILEID
    !###    into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(SPC), INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_SPC",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,SINGLE_COMPLEX_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_SPC")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_SPC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_SPC

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_SPC1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_SPC1
    !###  Description:
    !###    READ_BINARY_FILE_SPC1 reads 1 single precision
    !###    complex value from the binary file specified by FILEID
    !###    into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(SPC), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_SPC1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,SINGLE_COMPLEX_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_SPC1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_SPC1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_SPC1

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_DPC(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_DPC
    !###  Description:
    !###    READ_BINARY_FILE_DPC reads NUM_DATA double precision
    !###    complex values from the binary file specified by FILEID
    !###    into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(DPC), INTENT(OUT) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_DPC",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
        ENDIAN=FILE_CHANGE_ENDIAN
      ELSE
        ENDIAN=FILE_SAME_ENDIAN
      ENDIF
      CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & ENDIAN,NUM_DATA,DOUBLE_COMPLEX_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_DPC")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_DPC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_DPC

  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_FILE_DPC1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_FILE_DPC1
    !###  Description:
    !###    READ_BINARY_FILE_DPC1 reads 1 double precision
    !###    complex value from the binary file specified by FILEID
    !###    into DATA.
    !###  Parent-subroutine: READ_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(DPC), INTENT(OUT) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),ENDIAN
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("READ_BINARY_FILE_DPC1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        IF(FILEID%FILE_INFORMATION%ENDIAN_TYPE/=MACHINE_ENDIAN) THEN
          ENDIAN=FILE_CHANGE_ENDIAN
        ELSE
          ENDIAN=FILE_SAME_ENDIAN
        ENDIF
        CALL BINARYREADFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & ENDIAN,1,DOUBLE_COMPLEX_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_FILE_DPC1")
    RETURN
999 ERRORSEXITS("READ_BINARY_FILE_DPC1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_FILE_DPC1
  
  !
  !============================================================================
  !
  
  SUBROUTINE READ_BINARY_TAG_HEADER(FILEID,TAG,ERR,ERROR,*)

    !#### Subroutine: READ_BINARY_TAG_HEADER
    !###  Description:
    !###    READ_BINARY_TAG_HEADER reads a binary tag header from the
    !###    binary file specified by FILEID.

    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    TYPE(BINARY_TAG_TYPE), INTENT(OUT) :: TAG
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: INTDATA(2)

    ENTERS("READ_BINARY_TAG_HEADER",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL READ_BINARY_FILE(FILEID,2,INTDATA,ERR,ERROR,*999)
      TAG%INDEX=INTDATA(1)
      TAG%NUM_HEADER_BYTES=INTDATA(2)
      IF(TAG%NUM_HEADER_BYTES>MAXSTRLEN) THEN
        CALL FLAG_ERROR("Tag header length greater than maximum &
          &string length",ERR,ERROR,*999)
      ENDIF
      CALL READ_BINARY_FILE(FILEID,TAG%NUM_HEADER_BYTES,TAG%HEADER,ERR,&
        & ERROR,*999)
      SELECT CASE(FILEID%FILE_INFORMATION%BINARY_FILE_REVISION)
      CASE(0,1)
        TAG%NUM_SUBTAGS=0
      CASE(2)
        CALL READ_BINARY_FILE(FILEID,1,TAG%NUM_SUBTAGS,ERR,ERROR,*999)
      CASE DEFAULT
        CALL FLAG_ERROR("Invalid binary file header identity format",&
          & ERR,ERROR,*999)
      END SELECT
      IF(TAG%NUM_SUBTAGS==0) THEN
        CALL READ_BINARY_FILE(FILEID,1,TAG%NUM_BYTES,ERR,ERROR,*999)
      ELSE
        TAG%NUM_BYTES=0
      ENDIF
      IF(TAG%NUM_BYTES<0) CALL FLAG_ERROR("Invalid number of tag bytes",&
        & ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("READ_BINARY_TAG_HEADER")
    RETURN
999 ERRORSEXITS("READ_BINARY_TAG_HEADER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE READ_BINARY_TAG_HEADER

  !
  !============================================================================
  !
  
  SUBROUTINE RESET_BINARY_NUMBER_TAGS(FILEID,NUMBER_TAGS,ERR,ERROR,*)
    
    !#### Subroutine: RESET_BINARY_NUMBER_TAGS
    !###  Description:
    !###    RESET_BINARY_NUMBER_TAGS resets the number of tags in the 
    !###    binary file specified by FILEID to that specified by
    !###    NUMBER_TAGS

!!TODO: THIS PROBABLY NEEDS TO BE REWRITTEN TO ALLOW IT TO RESET THE NUMBER
!!      OF TAGS AT A SPECIFIED LEVEL IN THE TAG HIERACHY

    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUMBER_TAGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) ::  NUMBER_HEADER_BYTES,NUMBER_SKIP_BYTES

    ENTERS("RESET_BINARY_NUMBER_TAGS",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      !Rewind the file to the beginning
      CALL SET_BINARY_FILE(FILEID,FILE_BEGINNING,ERR,ERROR,*999)
      !Skip the identity and machine header sections
      CALL SKIP_CM_BINARY_HEADER(FILEID,CMISS_BINARY_MACHINE_HEADER,&
        & ERR,ERROR,*999)
      !Skip the file type and version
      NUMBER_SKIP_BYTES=INTEGER_SIZE+SINGLE_REAL_SIZE
      CALL SKIP_BINARY_FILE(FILEID,NUMBER_SKIP_BYTES,ERR,ERROR,*999)
      !Skip the file header
      CALL READ_BINARY_FILE(FILEID,1,NUMBER_HEADER_BYTES,ERR,ERROR,*999)
      CALL SKIP_BINARY_FILE(FILEID,NUMBER_HEADER_BYTES,ERR,ERROR,*999)
      !Set the file position
      CALL SET_BINARY_FILE(FILEID,FILE_CURRENT,ERR,ERROR,*999)
      !Write the new number of tags
      CALL WRITE_BINARY_FILE(FILEID,1,NUMBER_TAGS,ERR,ERROR,*999)
      !Go to the end of the file
      CALL SET_BINARY_FILE(FILEID,FILE_END,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("RESET_BINARY_NUMBER_TAGS")
    RETURN
999 ERRORSEXITS("RESET_BINARY_NUMBER_TAGS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE RESET_BINARY_NUMBER_TAGS

  !
  !============================================================================
  !
  
  SUBROUTINE SET_BINARY_FILE(FILEID,SET_CODE,ERR,ERROR,*)

    !#### Subroutine: SET_BINARY_FILE
    !###  Description:
    !###    SET_BINARY_FILE sets the position of the binary file
    !###    specified bye FILEID to that specified by SET_CODE.
    !###    Current SET_CODES are FILE_BEGINNING for the beginning of
    !###    a file, FILE_CURRENT for the current file position and
    !###    FILE_END for the end of a file.

    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: SET_CODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("SET_BINARY_FILE",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL BINARYSETFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & SET_CODE,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("SET_BINARY_FILE")
    RETURN
999 ERRORSEXITS("SET_BINARY_FILE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SET_BINARY_FILE

  !
  !============================================================================
  !
  
  SUBROUTINE SKIP_CM_BINARY_HEADER(FILEID,SKIP,ERR,ERROR,*)

    !#### Subroutine: SKIP_CM_BINARY_HEADER
    !###  Description:
    !###    SKIP_CM_BINARY_HEADER skips CMISS binary header
    !###    information in a binary file specified by FILEID. The
    !###    ammount of information skipping is controlled by SKIP.
    !###    Current values of SKIP are CMISS_BINARY_IDENTITY_HEADER
    !###    to skip the identity header, CMISS_BINARY_MACHINE_HEADER
    !###    to skip the machine header and CMISS_BINARY_FILE_HEADER
    !###    to skip the file header.

    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: SKIP
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: BINARY_FILE_REVISION,INTDATA(1),NUMBER_SKIP_BYTES
    CHARACTER(LEN=11) :: CHARDATA

    ENTERS("SKIP_CM_BINARY_HEADER",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(SKIP/=CMISS_BINARY_IDENTITY_HEADER.OR.&
        & SKIP/=CMISS_BINARY_MACHINE_HEADER.OR.&
        SKIP/=CMISS_BINARY_FILE_HEADER) &
        & CALL FLAG_ERROR("Invalid SKIP code",ERR,ERROR,*999)
      IF(SKIP>=CMISS_BINARY_IDENTITY_HEADER) THEN
        !Skip the identity header section
        CALL READ_BINARY_FILE(FILEID,2,CHARDATA,ERR,ERROR,*999)
        IF(CHARDATA(1:1)/=CHAR(CMISS_BINARY_IDENTITY)) &
          & CALL FLAG_ERROR("Not a CMISS binary file",ERR,ERROR,*999)
      ENDIF
      IF(SKIP>=CMISS_BINARY_MACHINE_HEADER) THEN
        !Skip the machine header section
        BINARY_FILE_REVISION=ICHAR(CHARDATA(2:2))
        SELECT CASE(BINARY_FILE_REVISION)
        CASE(0)
          NUMBER_SKIP_BYTES=3
        CASE(1,2)
          CALL READ_BINARY_FILE(FILEID,1,CHARDATA,ERR,ERROR,*999)
          NUMBER_SKIP_BYTES=ICHAR(CHARDATA(1:1))
        CASE DEFAULT
          CALL FLAG_ERROR("Unknown binary file header identity format",&
            & ERR,ERROR,*999)
        END SELECT
        CALL SKIP_BINARY_FILE(FILEID,NUMBER_SKIP_BYTES,ERR,ERROR,*999)
      ENDIF
      IF(SKIP>=CMISS_BINARY_FILE_HEADER) THEN
        !Skip the file header section
        NUMBER_SKIP_BYTES=INTEGER_SIZE+SINGLE_REAL_SIZE
        CALL SKIP_BINARY_FILE(FILEID,NUMBER_SKIP_BYTES, ERR,ERROR,*999)
        CALL READ_BINARY_FILE(FILEID,1,INTDATA,ERR,ERROR,*999)
        NUMBER_SKIP_BYTES=INTDATA(1)*CHARACTER_SIZE+INTEGER_SIZE
        CALL SKIP_BINARY_FILE(FILEID,NUMBER_SKIP_BYTES,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("SKIP_CM_BINARY_HEADER")
    RETURN
999 ERRORSEXITS("SKIP_CM_BINARY_HEADER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SKIP_CM_BINARY_HEADER

  !
  !============================================================================
  !
  
  SUBROUTINE SKIP_BINARY_FILE(FILEID,NUMBER_BYTES,ERR,ERROR,*)

    !#### Subroutine: SKIP_BINARY_FILE
    !###  Description:
    !###    SKIP_BINARY_FILE skips NUMBER_BYTES in a binary file
    !###    specified bye FILEID.

    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUMBER_BYTES
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("SKIP_BINARY_FILE",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUMBER_BYTES<=0) &
        & CALL FLAG_ERROR("NUMBER_BYTES to skip is <= 0",ERR,ERROR,*999)
      CALL BINARYSKIPFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & NUMBER_BYTES,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("SKIP_BINARY_FILE")
    RETURN
999 ERRORSEXITS("SKIP_BINARY_FILE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SKIP_BINARY_FILE

  !
  !============================================================================
  !
  
  RECURSIVE SUBROUTINE SKIP_BINARY_TAGS(FILEID,TAG,ERR,ERROR,*)

    !#### Subroutine: SKIP_BINARY_TAG
    !###  Description:
    !###    SKIP_BINARY_TAG (recursively) skips a binary tag that has 
    !###    already been read by READ_BINARY_TAG.

    !Argument Variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    TYPE(BINARY_TAG_TYPE), INTENT(IN) :: TAG
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),i
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR
    TYPE(BINARY_TAG_TYPE) :: SUBTAG

    ENTERS("SKIP_BINARY_TAGS",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(TAG%NUM_SUBTAGS>0) THEN
        DO i=1,TAG%NUM_SUBTAGS
          CALL READ_BINARY_TAG_HEADER(FILEID,SUBTAG,ERR,ERROR,*999)
          CALL SKIP_BINARY_TAGS(FILEID,SUBTAG,ERR,ERROR,*999)
        ENDDO !i
      ELSE
        CALL BINARYSKIPFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & TAG%NUM_BYTES,ERR,CERROR)        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("SKIP_BINARY_TAGS")
    RETURN
999 ERRORSEXITS("SKIP_BINARY_TAGS",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SKIP_BINARY_TAGS

  !
  !============================================================================
  !
  
  !#### Generic-Subroutine: WRITE_BINARY_FILE
  !###  Description:
  !###    WRITE_BINARY_FILE writes NUM_DATA elements from the binary
  !###    file specified by FILEID from DATA.
  !###  Child-Subroutines: WRITE_BINARY_FILE_INTG,
  !###    WRITE_BINARY_FILE_INTG1,WRITE_BINARY_FILE_SINTG,
  !###    WRITE_BINARY_FILE_SNT1,WRITE_BINARY_FILE_LINTG,
  !###    WRITE_BINARY_FILE_LINTG1,WRITE_BINARY_FILE_SP,
  !###    WRITE_BINARY_FILE_SP1,WRITE_BINARY_FILE_DP,
  !###    WRITE_BINARY_FILE_DP1,WRITE_BINARY_FILE_CHARACTER,
  !###    WRITE_BINARY_FILE_LOGICAL,WRITE_BINARY_FILE_LOGICAL1,
  !###    WRITE_BINARY_FILE_SPC,WRITE_BINARY_FILE_SPC1,
  !###    WRITE_BINARY_FILE_DPC,WRITE_BINARY_FILE_DPC1
  
  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_INTG(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_INTG
    !###  Description:
    !###    WRITE_BINARY_FILE_INTG writes NUM_DATA integer values to
    !###    the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA, DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN    
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,INTEGER_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("WRITE_BINARY_FILE_INTG")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_INTG",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_INTG
  
  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_INTG1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_INTG1
    !###  Description:
    !###    WRITE_BINARY_FILE_INTG1 writes 1 integer value to
    !###    the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA, DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_INTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN,1,INTEGER_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("WRITE_BINARY_FILE_INTG1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_INTG1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_INTG1
  
  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_SINTG(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_SINTG
    !###  Description:
    !###    WRITE_BINARY_FILE_SINTG writes NUM_DATA short integer
    !###    values to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(SINTG), INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_SINTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,SHORT_INTEGER_TYPE,DATA,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_SINTG")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_SINTG",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_SINTG

  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_SINTG1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_SINTG1
    !###  Description:
    !###    WRITE_BINARY_FILE_SINTG1 writes 1 short integer value to
    !###    the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(SINTG), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_SINTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN,1,SHORT_INTEGER_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_SINTG1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_SINTG1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_SINTG1

  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_LINTG(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_LINTG
    !###  Description:
    !###    WRITE_BINARY_FILE_LINTG writes NUM_DATA long integer
    !###    values to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(LINTG), INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_LINTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,LONG_INTEGER_TYPE,DATA,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_LINTG")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_LINTG",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_LINTG

  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_LINTG1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_LINTG1
    !###  Description:
    !###    WRITE_BINARY_FILE_LINTG1 writes 1 long integer value to
    !###    the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    INTEGER(LINTG), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_LINTG1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN,1,LONG_INTEGER_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_LINTG1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_LINTG1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_LINTG1

  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_SP(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_SP
    !###  Description:
    !###    WRITE_BINARY_FILE_SP writes NUM_DATA single precision
    !###    real values to the binary file specified by FILEID from
    !###    DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(SP), INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,SINGLE_REAL_TYPE,DATA,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("WRITE_BINARY_FILE_SP")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_SP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_SP

  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_FILE_SP1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_SP1
    !###  Description:
    !###    WRITE_BINARY_FILE_SP1 writes 1 single precision real
    !###    value to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(SP), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_SP1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN,1,SINGLE_REAL_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF

    EXITS("WRITE_BINARY_FILE_SP1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_SP1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_SP1

  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_DP(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_DP
    !###  Description:
    !###    WRITE_BINARY_FILE_DP writes NUM_DATA double precision
    !###    real values to the binary file specified by FILEID from
    !###    DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(DP), INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN, NUM_DATA, DOUBLE_REAL_TYPE, DATA, ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
      
    EXITS("WRITE_BINARY_FILE_DP")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_DP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_DP

  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_DP1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_DP1
    !###  Description:
    !###    WRITE_BINARY_FILE_DP1 writes 1 double precision real
    !###    value to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    REAL(DP), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_DP1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN,1,DOUBLE_REAL_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
      
    EXITS("WRITE_BINARY_FILE_DP1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_DP1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_DP1

  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_CHARACTER(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_CHARACTER
    !###  Description:
    !###    WRITE_BINARY_FILE_CHARACTER writes NUM_DATA character
    !###    values to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    CHARACTER(LEN=*), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100),CSTRING(250)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_CHARACTER",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL F2CSTRING(CSTRING,DATA,ERR,ERROR,*999)
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,CHARACTER_TYPE,CSTRING,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_CHARACTER")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_CHARACTER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_CHARACTER
  
  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_LOGICAL(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_LOGICAL
    !###  Description:
    !###    WRITE_BINARY_FILE_LOGICAL writes NUM_DATA logical values
    !###    to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    LOGICAL, INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_LOGICAL",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,LOGICAL_TYPE,DATA,ERR,CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_LOGICAL")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_LOGICAL",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_LOGICAL
  
  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_LOGICAL1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_LOGICAL1
    !###  Description:
    !###    WRITE_BINARY_FILE_LOGICAL1 writes 1 logical value to
    !###    the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    LOGICAL, INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_LOGICAL1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN,1,LOGICAL_TYPE,DATA,ERR,CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_LOGICAL1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_LOGICAL1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_LOGICAL1

  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_SPC(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_SPC
    !###  Description:
    !###    WRITE_BINARY_FILE_SPC writes NUM_DATA single precision
    !###    complex  values to the binary file specified by FILEID
    !###    from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(SPC), INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_SPC",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,SINGLE_COMPLEX_TYPE,DATA,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_SPC")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_SPC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_SPC
  
  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_SPC1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_SPC1
    !###  Description:
    !###    WRITE_BINARY_FILE_SPC1 writes 1 single precision complex
    !###    value to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(SPC), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_SPC1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN      
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data items not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN, 1, SINGLE_COMPLEX_TYPE, DATA, ERR,&
          & CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_SPC1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_SPC1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_SPC1
  
  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_DPC(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_DPC
    !###  Description:
    !###    WRITE_BINARY_FILE_DPC writes NUM_DATA double precision
    !###    complex values to the binary file specified by FILEID
    !###    from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(DPC), INTENT(IN) :: DATA(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_DPC",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
        & FILE_SAME_ENDIAN,NUM_DATA,DOUBLE_COMPLEX_TYPE,DATA,ERR,&
        & CERROR)
      IF(ERR/=0) THEN
        CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
        CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_DPC")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_DPC",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_DPC

  !
  !============================================================================
  !

  SUBROUTINE WRITE_BINARY_FILE_DPC1(FILEID,NUM_DATA,DATA,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_FILE_DPC1
    !###  Description:
    !###    WRITE_BINARY_FILE_ writes 1 double precision complex
    !###    value to the binary file specified by FILEID from DATA.
    !###  Parent-subroutine: WRITE_BINARY_FILE

    !Argument variables 
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    INTEGER(INTG), INTENT(IN) :: NUM_DATA
    COMPLEX(DPC), INTENT(IN) :: DATA
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    ENTERS("WRITE_BINARY_FILE_DPC1",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      IF(NUM_DATA/=1) THEN
        CALL FLAG_ERROR("Number of data item not equal to one",ERR,ERROR,*999)
      ELSE
        CALL BINARYWRITEFILE(FILEID%FILE_INFORMATION%FILE_NUMBER,&
          & FILE_SAME_ENDIAN, 1, DOUBLE_COMPLEX_TYPE, DATA, ERR,&
          & CERROR)
        IF(ERR/=0) THEN
          CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
          CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_FILE_DPC1")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_FILE_DPC1",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_FILE_DPC1

  !
  !============================================================================
  !
  
  SUBROUTINE WRITE_BINARY_TAG_HEADER(FILEID,TAG,ERR,ERROR,*)

    !#### Subroutine: WRITE_BINARY_TAG_HEADER
    !###  Description:
    !###    WRITE_BINARY_TAG_HEADER writes a binary tag header from
    !###    the binary file specified by FILEID.

    !Argument variables
    TYPE(BINARY_FILE_TYPE), INTENT(IN) :: FILEID
    TYPE(BINARY_TAG_TYPE), INTENT(OUT) :: TAG
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: INTDATA(2)

    ENTERS("WRITE_BINARY_TAG_HEADER",ERR,ERROR,*999)

    IF(ASSOCIATED(FILEID%FILE_INFORMATION)) THEN
      INTDATA(1)=TAG%INDEX
      TAG%NUM_HEADER_BYTES=LEN_TRIM(TAG%HEADER)
      INTDATA(2)=TAG%NUM_HEADER_BYTES
      CALL WRITE_BINARY_FILE(FILEID,2,INTDATA,ERR,ERROR,*999)
      CALL WRITE_BINARY_FILE(FILEID,TAG%NUM_HEADER_BYTES,&
        & TAG%HEADER,ERR,ERROR,*999)
      CALL WRITE_BINARY_FILE(FILEID,1,TAG%NUM_SUBTAGS,ERR,ERROR,*999)
      IF(TAG%NUM_SUBTAGS>0) &
        & CALL WRITE_BINARY_FILE(FILEID,1,TAG%NUM_BYTES,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Invalid FILEID",ERR,ERROR,*999)
    ENDIF
    
    EXITS("WRITE_BINARY_TAG_HEADER")
    RETURN
999 ERRORSEXITS("WRITE_BINARY_TAG_HEADER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE WRITE_BINARY_TAG_HEADER

  !
  !============================================================================
  !
  
END MODULE BINARY_FILE
