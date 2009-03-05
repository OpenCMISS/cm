!> \file
!> $Id: distributed_matrix_vector_IO.f90 248 2008-11-28 11:14:17Z chrispbradley $
!> \author Chris Bradley
!> \brief This module handles all distributed matrix vector IO routines.
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


!> This module handles all distributed_matrix_vector_IO.
MODULE DISTRIUBTED_MATRIX_VECTOR_IO
  
  USE BASE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES
  
  IMPLICIT NONE
    

  PRIVATE
  
  !Module types

  !Module variables
  
  INTEGER(INTG) :: FILE_UNIT
  INTEGER(INTG) :: myrank
  TYPE(VARYING_STRING) :: FILE_NAME, FILE_FORMAT, FILE_STATUS

  !Interfaces

  PUBLIC DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE, DISTRIBUTED_MATRIX_VECTOR_IO_OPEN, DISTRIBUTED_VECTOR_IO_READ, DISTRIBUTED_VECTOR_IO_WRITE, DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_START, DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_FINISH, DISTRIBUTED_MATRIX_VECTOR_IO_FILENAME_SET, DISTRIBUTED_MATRIX_VECTOR_IO_FILE_OPEN, DISTRIBUTED_MATRIX_VECTOR_IO_FILE_CLOSE,DISTRIBUTED_MATRIX_VECTOR_IO_MODE_SET,DISTRIBUTED_MATRIX_VECTOR_IO_HEADER_WRITE
     
CONTAINS
 
  !
  !================================================================================================================================
  !
  
  !>Closes a distributed matrix vector IO file. 
  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE(ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE",ERR,ERROR,*999)
    
     
    CALL EXITS("DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE")
    RETURN
999 CALL ERRORS("DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE",ERR,ERROR)
    CALL EXITS("DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE")
    RETURN 1
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_CLOSE

  !
  !================================================================================================================================
  !
  
  !>Opens a distributed matrix vector IO file. 
  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_OPEN(ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DISTRIBUTED_MATRIX_VECTOR_IO_OPEN",ERR,ERROR,*999)
    
     
    CALL EXITS("DISTRIBUTED_MATRIX_VECTOR_IO_OPEN")
    RETURN
999 CALL ERRORS("DISTRIBUTED_MATRIX_VECTOR_IO_OPEN",ERR,ERROR)
    CALL EXITS("DISTRIBUTED_MATRIX_VECTOR_IO_OPEN")
    RETURN 1
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_OPEN

  !
  !================================================================================================================================
  !
  
  !>Reads a distributed vector from a file file. 
  SUBROUTINE DISTRIBUTED_VECTOR_IO_READ(FILE,OFFSET,DISTRIBUTED_VECTOR,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FILE !<Should be pointer to file object
    INTEGER(INTG), INTENT(INOUT) :: OFFSET !Should be large integer
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR !<On exit, a pointer to the distributed vector read. Should not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DISTRIBUTED_VECTOR_IO_READ",ERR,ERROR,*999)
    
     
    CALL EXITS("DISTRIBUTED_VECTOR_IO_READ")
    RETURN
999 CALL ERRORS("DISTRIBUTED_VECTOR_IO_READ",ERR,ERROR)
    CALL EXITS("DISTRIBUTED_VECTOR_IO_READ")
    RETURN 1
  END SUBROUTINE DISTRIBUTED_VECTOR_IO_READ
  
  !
  !================================================================================================================================
  !
  
  !>Reads a distributed vector from a file file. 
  SUBROUTINE DISTRIBUTED_VECTOR_IO_WRITE(FILE,OFFSET,DISTRIBUTED_VECTOR,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: FILE !<Should be pointer ../../../old_cm/examples/AnalyticLaplace/object/x86_64-linux-debug-intel/distriubted_matrix_vector_io.modto file object
    INTEGER(INTG), INTENT(INOUT) :: OFFSET !Should be large integer
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR !<A pointer to the distributed vector to write
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DISTRIBUTED_VECTOR_IO_WRITE",ERR,ERROR,*999)

    IF(ASSOCIATED(DISTRIBUTED_VECTOR)) THEN
    ELSE
      CALL FLAG_ERROR("Distributed vector is not associated.",ERR,ERROR,*999)
    ENDIF    
     
    CALL EXITS("DISTRIBUTED_VECTOR_IO_WRITE")
    RETURN
999 CALL ERRORS("DISTRIBUTED_VECTOR_IO_WRITE",ERR,ERROR)
    CALL EXITS("DISTRIBUTED_VECTOR_IO_WRITE")
    RETURN 1
  END SUBROUTINE DISTRIBUTED_VECTOR_IO_WRITE
  
  !
  !================================================================================================================================
  !





  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_START(MY_COMPUTATIONAL_NODE_NUMBER, FILE_STS,  ERR, ERROR, *)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: MY_COMPUTATIONAL_NODE_NUMBER
    TYPE(VARYING_STRING), INTENT(IN) :: FILE_STS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    
    myrank=MY_COMPUTATIONAL_NODE_NUMBER
    FILE_STATUS=FILE_STS
    ! 
    ! This function initiates the data structure and allocates variables and buffers for the IO to come.
    !  
    !
!    write(*,*) 'DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_START'
    
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_START

  

  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_FINISH(ERR, ERROR, *)

    !Argument variables 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    !
    ! This function initiates the data structure and allocates variables and buffers for the IO to come.
    !  
    !
!    write(*,*)'DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_FINISH'
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_CREATE_FINISH



  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_FILENAME_SET(FILE, ERR, ERROR, *)

    !Argument variables 
    TYPE(VARYING_STRING), INTENT(IN) :: FILE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !
    !
    ! This routine sets the file name in the io file structure and checks for a valid file name.
    !
    !
    CHARACTER*80   :: FILENM, RNR
    
    IF (len(file) > 0) THEN 
       write(RNR,fmt='(I,A)') myrank,'.dta'
       RNR=ADJUSTL(RNR)
       FILENM=FILE//'-'//RNR
!       write(*,*) FILENAME
    ELSE       
!       write(*,*) 'Not a valid file name'
       ERR=1
       ERROR='Not a valid file name'
       CALL ERRORS("Distributed vector is not associated.",ERR,ERROR)       
    ENDIF
    FILE_NAME=FILENM
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_FILENAME_SET




  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_MODE_SET(FILE_FMT, ERR, ERROR, *)

    !Argument variables 
    TYPE(VARYING_STRING), INTENT(IN) :: FILE_FMT
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string    !
    !
    ! This set the mode - which might be simple writing of the fiels local to the current rank in question,
    ! or a more sophisticated method where data is sorted in ways that makes it possible to read it in with a 
    ! different set of ranks or and export format which might use something like parallel HDF5.
    !
    ! There should also be a method called mode_get that will check i the files are written in the mode asked for.
    ! 
    !
    FILE_FORMAT=FILE_FMT
    !
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_MODE_SET



  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_FILE_OPEN(ERR, ERROR, *)

    !Argument variables 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !
    ! This function opens the file and uses tha data structure for this initiated bye previous calls to 
    ! CREATE_START.
    !
    INTEGER :: iostat
    character*9 :: pos
    !
    ! Pick a number
    FILE_UNIT=20
    ! Be sure this unit is not open already.    
    inquire(unit=FILE_UNIT, position=pos, iostat=iostat)
    do while (pos/='UNDEFINED') 
       FILE_UNIT=FILE_UNIT+1
       inquire(unit=FILE_UNIT, position=pos, iostat=iostat)
    enddo
    !
!    write(*,*) 'File open info :', &
!    ' File name :',CHAR(FILE_NAME), &     
!    ' File UNIT :',FILE_UNIT,  &
!    ' File status :',CHAR(FILE_STATUS), &
!    ' File format :',CHAR(FILE_FORMAT) 

    OPEN(UNIT=FILE_UNIT, FILE=CHAR(FILE_NAME), STATUS=CHAR(FILE_STATUS), FORM=CHAR(FILE_FORMAT), iostat=iostat, ERR=998)
    !
    !
    if (iostat/=0) then       
       write(*,*) 'iostat /= 0', iostat
998    CALL ERRORS("DISTRIBUTED_MATRIX_VECTOR_IO_FILE_OPEN",ERR,ERROR)
       CALL EXITS("DISTRIBUTED_MATRIX_VECTOR_IO_FILE_OPEN")
    else
       RETURN 
    endif
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_FILE_OPEN



  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_FILE_CLOSE(ERR, ERROR, *)

    !Argument variables 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string 
    !
    ! This routine closes the file open with DISTRIBUTED_MATRIX_VECTOR_IO_FILE_OPEN.
    !
    INTEGER :: iostat
    !
    !
    !    
    close(UNIT=FILE_UNIT, iostat=iostat)
    !
    if (iostat/=0) then     
       write(*,*) 'iostat /= 0', iostat
       CALL ERRORS("FIELD_IO_FORTRAN_FILE_CLOSE",ERR,ERROR)
       CALL EXITS("FIELD_IO_FORTRAN_FILE_CLOSE")
    endif
    RETURN 
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_FILE_CLOSE



! This is might not correct - do we want one single routine like the below for each one of the
! implemented data types ? string, integer, double ? 
! We do write binary data and hence is is no information about types. One routine for a header and 
! one for binary data.

  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_HEADER_WRITE(ERR, ERROR, *)

    !Argument variables 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string 
    !
    ! Write out a buffer given by the paramaters. This will be a field or any other data type.  
    !
    !
    !
    write(FILE_UNIT) myrank
    ! 
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_HEADER_WRITE




  SUBROUTINE  DISTRIBUTED_MATRIX_VECTOR_IO_DATA_WRITE(ERR, ERROR, *)
  
    !Argument variables 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string 
    !
    ! Write out a buffer given by the paramaters. This will be a field or any other data type.  
    !
    !
    INTEGER :: k=2
    !
    write(FILE_UNIT) k 

    
 
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_DATA_WRITE





  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_HEADER_READ
    !
    !
    ! Read a chunk of data from the file.
    !
    !
    !
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_HEADER_READ






  SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_DATA_READ
    !
    !
    ! Read a chunk of data from the file.
    !
    !
    !
  END SUBROUTINE DISTRIBUTED_MATRIX_VECTOR_IO_DATA_READ




  
 
END MODULE DISTRIUBTED_MATRIX_VECTOR_IO

