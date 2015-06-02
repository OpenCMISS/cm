!> \file
!> \author Chris Bradley
!> \brief This module contains CMISS MPI routines.
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

!> This module contains CMISS MPI routines.
MODULE CMISS_MPI
  
  USE BASE_ROUTINES
  USE CONSTANTS
  USE KINDS
#ifndef NOMPIMOD
  USE MPI
#endif
  USE ISO_VARYING_STRING
  USE STRINGS

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC MPI_ERROR_CHECK
 
CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks to see if an MPI error has occured during an MPI call and flags a CMISS error it if it has.
  SUBROUTINE MPI_ERROR_CHECK(ROUTINE,MPI_ERR_CODE,ERR,ERROR,*)
  
    !Argument Variables
    CHARACTER(LEN=*) :: ROUTINE !<The name of the MPI routine that has just been called.
    INTEGER(INTG), INTENT(IN) :: MPI_ERR_CODE !<The MPI error code returned from the MPI routine.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR, MPI_ERR_STR_LENGTH
    CHARACTER(LEN=MAXSTRLEN) :: MPI_ERR_STR
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MPI_ERROR_CHECK",ERR,ERROR,*999)

    IF(MPI_ERR_CODE/=MPI_SUCCESS) THEN
      CALL MPI_ERROR_STRING(MPI_ERR_CODE,MPI_ERR_STR,MPI_ERR_STR_LENGTH,MPI_IERROR)
      LOCAL_ERROR="MPI error "//TRIM(NUMBER_TO_VSTRING(MPI_ERR_CODE,"*",ERR,ERROR))//" ("// &
        & MPI_ERR_STR(1:MPI_ERR_STR_LENGTH)//") in "//ROUTINE(1:LEN_TRIM(ROUTINE))
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MPI_ERROR_CHECK")
    RETURN
999 CALL ERRORS("MPI_ERROR_CHECK",ERR,ERROR)
    CALL EXITS("MPI_ERROR_CHECK")
    RETURN 1
  END SUBROUTINE MPI_ERROR_CHECK

  !
  !================================================================================================================================
  !
    
END MODULE CMISS_MPI
