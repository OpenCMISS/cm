!> \file
!> $Id: timer_f.f90 27 2007-07-24 16:52:51Z cpb $
!> \author Chris Bradley
!> \brief This module contains routines for timing the program.
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

!> This module contains routines for timing the program.
MODULE TIMER

  USE BASE_ROUTINES
  USE KINDS
  USE F90C
  USE ISO_VARYING_STRING

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module variables
  
  !CPU time parameters
  INTEGER(INTG), PARAMETER :: USER_CPU=1
  INTEGER(INTG), PARAMETER :: SYSTEM_CPU=2
  INTEGER(INTG), PARAMETER :: TOTAL_CPU=3
  
  !Module variables

  !Interfaces

  INTERFACE

    SUBROUTINE CPUTIMER(RETURN_TIME, TIME_TYPE, ERR, CERROR)
      !DEC$ ATTRIBUTES C, reference, alias:'_cputimer' :: cputimer
      USE KINDS
      REAL(DP), INTENT(OUT) :: RETURN_TIME
      INTEGER(INTG), INTENT(IN) :: TIME_TYPE
      INTEGER(INTG), INTENT(OUT) :: ERR
      INTEGER(INTG), INTENT(OUT) :: CERROR(*)
    END SUBROUTINE CPUTIMER

  END INTERFACE

  PUBLIC USER_CPU,SYSTEM_CPU,TOTAL_CPU,CPU_TIMER

CONTAINS

  !
  !============================================================================
  !

  SUBROUTINE CPU_TIMER(TIME_TYPE,TIME,ERR,ERROR,*)

    !#### Subroutine: CPU_TIMER
    !###  Description:
    !###    CPU_TIMER returns the CPU time in TIME(1). TIME_TYPE
    !###    indicates the type of time required i.e. USER_CPU,
    !###    SYSTEM_CPU or TOTAL_CPU

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TIME_TYPE
    REAL(SP), INTENT(OUT) :: TIME(*)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local variables
    REAL(DP) :: RETURN_TIME
    INTEGER(INTG) :: CERROR(100)
    CHARACTER(LEN=MAXSTRLEN) :: DUMMY_ERROR

    CALL ENTERS("CPU_TIMER",ERR, ERROR,*999)
    
    CALL CPUTIMER(RETURN_TIME,TIME_TYPE,ERR,CERROR)
    TIME(1)=REAL(RETURN_TIME,SP)
    IF(ERR/=0) THEN
      CALL C2FSTRING(CERROR,DUMMY_ERROR,ERR,ERROR,*999)
      CALL FLAG_ERROR(DUMMY_ERROR,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("CPU_TIMER")
    RETURN
999 CALL ERRORS("CPU_TIMER",ERR,ERROR)
    CALL EXITS("CPU_TIMER")
    RETURN 1
  END SUBROUTINE CPU_TIMER

END MODULE TIMER
