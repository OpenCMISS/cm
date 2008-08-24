!> \file
!> $Id: cmiss.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Chris Bradley
!> \brief The top level cmiss module.
!>
!> \mainpage openCMISS Documentation
!>
!> An open source interactive computer program for Continuum Mechanics, Image analysis, Signal processing and System
!> Identification. Target usage: Bioengineering application of finite element analysis, boundary element and collocation
!> techniques.
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

!> The top level cmiss module.
MODULE CMISS

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE COORDINATE_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MACHINE_CONSTANTS
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(INTG), PARAMETER :: CMISS_MAJOR_VERSION = 0
  INTEGER(INTG), PARAMETER :: CMISS_MINOR_VERSION = 1
  INTEGER(INTG), PARAMETER :: CMISS_BUILD_VERSION = 1
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC CMISS_MAJOR_VERSION,CMISS_MINOR_VERSION,CMISS_BUILD_VERSION
  
  PUBLIC CMISS_WRITE_ERROR,CMISS_FINALISE,CMISS_INITIALISE
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalises CMISS.
  SUBROUTINE CMISS_FINALISE(ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error code
    !Local Variables

    !Finalise the problems
    CALL PROBLEMS_FINALISE(ERR,ERROR,*999)
    !Finalise the regions
    CALL REGIONS_FINALISE(ERR,ERROR,*999)
    !Finalise the coordinate systems
    CALL COORDINATE_SYSTEMS_FINALISE(ERR,ERROR,*999)
    !Finalise bases
    CALL BASES_FINALISE(ERR,ERROR,*999)
    !Finalise computational enviroment
    CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(ERR,ERROR,*999)
    !Initialise the base routines
    CALL BASE_ROUTINES_FINALISE(ERR,ERROR,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE CMISS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises CMISS.
  SUBROUTINE CMISS_INITIALISE(ERR,ERROR,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local Variables

    !Initialise the base routines
    CALL BASE_ROUTINES_INITIALISE(ERR,ERROR,*999)
    !Intialise the computational environment
    CALL COMPUTATIONAL_ENVIRONMENT_INITIALISE(ERR,ERROR,*999)
    !Intialise the bases
    CALL BASES_INITIALISE(ERR,ERROR,*999) !BASES is the pl of basis
    !Initialise the coordinate systems
    CALL COORDINATE_SYSTEMS_INITIALISE(ERR,ERROR,*999)
    !Initialise the regions 
    CALL REGIONS_INITIALISE(ERR,ERROR,*999)
    !Initialise the problems
    CALL PROBLEMS_FINALISE(ERR,ERROR,*999)
    
    RETURN
999 RETURN 1
  END SUBROUTINE CMISS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Writes the error string to screen.
  SUBROUTINE CMISS_WRITE_ERROR(ERR,ERROR)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: INDENT,POSITION
    CHARACTER(LEN=MAXSTRLEN) :: INDENT_STRING=">>"
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
    WRITE(*,'(">>ERROR: ",I5,": ",A)') ERR,CHAR(EXTRACT(ERROR,1,POSITION-1))
    !CPB 20/02/07 aix compiler does not like varying strings so split the remove statement up into two statements
    LOCAL_ERROR=REMOVE(ERROR,1,POSITION)
    ERROR=LOCAL_ERROR
    POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
    INDENT=4
    DO WHILE(POSITION/=0)
      WRITE(*,'(A)') INDENT_STRING(1:INDENT)//CHAR(EXTRACT(ERROR,1,POSITION-1))
      !CPB 20/02/07 aix compiler does not like varying strings so split the remove statement up into two statements
      LOCAL_ERROR=REMOVE(ERROR,1,POSITION)
      ERROR=LOCAL_ERROR
      POSITION=INDEX(ERROR,ERROR_SEPARATOR_CONSTANT)
      INDENT=INDENT+2
    ENDDO
    WRITE(*,'(A)') INDENT_STRING(1:INDENT)//CHAR(ERROR)
    ERR=0
    ERROR=""
    
  END SUBROUTINE CMISS_WRITE_ERROR

  !
  !================================================================================================================================
  !

END MODULE CMISS
