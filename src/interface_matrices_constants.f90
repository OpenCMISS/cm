!> \file
!> \author Chris Bradley
!> \brief This module handles all constants shared across interface condition routines.
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

!> This module defines all constants shared across interface matrices routines.
MODULE INTERFACE_MATRICES_CONSTANTS

  USE KINDS

  IMPLICIT NONE

  !Module parameters

  !> \addtogroup INTERFACE_MATRICES_ROUTINES_InterfaceMatricesTimeDependenceTypes INTERFACE_MATRICES_ROUTINES::InterfaceMatricesTimeDependenceTypes
  !> \brief Interface matrices time dependency types
  !> \see INTERFACE_MATRICES_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_INTERFACE_MATRIX_TYPES=4
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_STATIC=1 !<Interface matrix is of static type \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesTimeDependenceTypes,INTERFACE_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_QUASI_STATIC=2 !<Interface matrix is of quasi-static type \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesTimeDependenceTypes,INTERFACE_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC=3 !<Interface matrix is of first order dynamic type \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesTimeDependenceTypes,INTERFACE_MATRICES_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC=4 !<Interface matrix is of second order dynamic type \see INTERFACE_MATRICES_ROUTINES_InterfaceMatricesTimeDependenceTypes,INTERFACE_MATRICES_ROUTINES
  !>@}

END MODULE INTERFACE_MATRICES_CONSTANTS
