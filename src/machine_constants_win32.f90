!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all machine dependent constants for Win32 systems.
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

!> This module contains all machine dependent constants for Win32 systems.
MODULE MACHINE_CONSTANTS

  USE CONSTANTS
  USE KINDS
  
  IMPLICIT NONE

  !Module parameters

  !Machine constants
  INTEGER(INTG), PARAMETER :: MACHINE_TYPE=PC_COMPUTER
  INTEGER(INTG), PARAMETER :: MACHINE_OS=WINDOWS_OS
  INTEGER(INTG), PARAMETER :: MACHINE_ENDIAN=BIG_ENDIAN_NUMBER
  INTEGER(INTG), PARAMETER :: MACHINE_CHAR_FORMAT=ASCII_CHARACTER
  INTEGER(INTG), PARAMETER :: MACHINE_INT_FORMAT=TWOS_COMPLEMENT_INTEGER
  INTEGER(INTG), PARAMETER :: MACHINE_SP_FORMAT=SPIEEE_NUMBER
  INTEGER(INTG), PARAMETER :: MACHINE_DP_FORMAT=DPIEEE_NUMBER
  INTEGER(INTG), PARAMETER :: INTEGER_SIZE=4
  INTEGER(INTG), PARAMETER :: SHORT_INTEGER_SIZE=2
  INTEGER(INTG), PARAMETER :: LONG_INTEGER_SIZE=8
  INTEGER(INTG), PARAMETER :: SINGLE_REAL_SIZE=4
  INTEGER(INTG), PARAMETER :: DOUBLE_REAL_SIZE=8
  INTEGER(INTG), PARAMETER :: CHARACTER_SIZE=1
  INTEGER(INTG), PARAMETER :: LOGICAL_SIZE=4
  INTEGER(INTG), PARAMETER :: SINGLE_COMPLEX_SIZE=8
  INTEGER(INTG), PARAMETER :: DOUBLE_COMPLEX_SIZE=16

  CHARACTER(LEN=1), PARAMETER :: ERROR_SEPARATOR_CONSTANT=CHAR(0)
  
END MODULE MACHINE_CONSTANTS
