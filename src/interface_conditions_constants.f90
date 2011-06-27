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

!> This module defines all constants shared across interface condition routines.
MODULE INTERFACE_CONDITIONS_CONSTANTS

  USE KINDS

  IMPLICIT NONE

  !Module parameters

  !> \addtogroup INTERFACE_CONDITIONS_CONSTANTS_Methods INTERFACE_CONDITIONS_CONSTANTS::Methods
  !> \brief Interface condition methods.
  !> \see INTERFACE_CONDITIONS_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD=1 !<Lagrange multipliers interface condition method. \see INTERFACE_CONDITIONS_CONSTANTS_Methods,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD=2 !<Augmented Lagrange multiplers interface condition method. \see INTERFACE_CONDITIONS_CONSTANTS_Methods,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_PENALTY_METHOD=3 !<Penalty interface condition method. \see INTERFACE_CONDITIONS_CONSTANTS_Methods,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_POINT_TO_POINT_METHOD=4 !<Point to point interface condition method. \see INTERFACE_CONDITIONS_CONSTANTS_Methods,INTERFACE_CONDITIONS_CONSTANTS
  !>@}

  !> \addtogroup INTERFACE_CONDITIONS_Operators INTERFACE_CONDITIONS_CONSTANTS::Operators
  !> \brief Interface condition operators.
  !> \see INTERFACE_CONDITIONS_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FIELD_CONTINUITY_OPERATOR=1 !<Continuous field operator, i.e., lambda.(u_1-u_2). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FIELD_NORMAL_CONTINUITY_OPERATOR=2 !<Continuous field normal operator, i.e., lambda(u_1.n_1-u_2.n_2). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SOLID_FLUID_OPERATOR=3 !<Solid fluid operator, i.e., lambda.(v_f-du_s/dt). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR=4 !<Solid fluid normal operator, i.e., lambda(v_f.n_f-du_s/dt.n_s). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  !>@}

END MODULE INTERFACE_CONDITIONS_CONSTANTS
