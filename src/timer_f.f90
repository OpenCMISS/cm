!> \file
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

!> This module contains routines for timing the program.
MODULE TIMER

  USE BASE_ROUTINES
  USE CONSTANTS
  USE BASE_ROUTINES
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module variables
  
  !CPU time parameters
  !> \addtogroup TIMER_TimerType TIMER::TimerType
  !> \brief Timer type parameter
  !> \see TIMER
  !>@{  
  INTEGER(INTG), PARAMETER :: USER_CPU=1 !<User CPU time type \see TIMER_TimerType,TIMER
  INTEGER(INTG), PARAMETER :: SYSTEM_CPU=2 !<System CPU time type \see TIMER_TimerType,TIMER
  INTEGER(INTG), PARAMETER :: TOTAL_CPU=3 !<Total CPU (i.e. User + System) time type \see TIMER_TimerType,TIMER
  !>@}
  
  !Module variables

  !Interfaces

  INTERFACE

    SUBROUTINE CPUTIMER(RETURN_TIME, TIME_TYPE, ERR, CERROR) BIND(C,NAME="CPUTimer")
      USE ISO_C_BINDING
      REAL(C_DOUBLE), INTENT(OUT) :: RETURN_TIME
      INTEGER(C_INT), INTENT(IN) :: TIME_TYPE
      INTEGER(C_INT), INTENT(OUT) :: ERR
      CHARACTER(C_CHAR), INTENT(OUT) :: CERROR(*)
    END SUBROUTINE CPUTIMER

  END INTERFACE

  PUBLIC USER_CPU,SYSTEM_CPU,TOTAL_CPU,CPU_TIMER

CONTAINS

  !
  !============================================================================
  !

  !>CPU_TIMER returns the CPU time in TIME(1). TIME_TYPE indicates the type of time required.
  SUBROUTINE CPU_TIMER(TIME_TYPE,TIME,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: TIME_TYPE !<The required time type \see TIMER_TimerType,TIMER
    REAL(SP), INTENT(OUT) :: TIME(*) !<On return, the requested time.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: RETURN_TIME
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: CERROR

    ENTERS("CPU_TIMER",ERR, ERROR,*999)
    
    CALL CPUTIMER(RETURN_TIME,TIME_TYPE,ERR,CERROR)
    TIME(1)=REAL(RETURN_TIME,SP)
    IF(ERR/=0) CALL FlagError(CERROR,ERR,ERROR,*999)
    
    EXITS("CPU_TIMER")
    RETURN
999 ERRORSEXITS("CPU_TIMER",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CPU_TIMER

  !
  !============================================================================
  !

END MODULE TIMER
