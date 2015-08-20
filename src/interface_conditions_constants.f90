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
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FLS_CONTACT_OPERATOR=3 !<Frictionless contact operator, i.e., lambda.(x_1.n-x_2.n). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR=4 !<Frictionless contact operator, reproject at each newton iteration and has geometric linearisation terms i.e., lambda.(x_1.n-x_2.n). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SOLID_FLUID_OPERATOR=5 !<Solid fluid operator, i.e., lambda.(v_f-du_s/dt). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SOLID_FLUID_NORMAL_OPERATOR=6 !<Solid fluid normal operator, i.e., lambda(v_f.n_f-du_s/dt.n_s). \see INTERFACE_CONDITIONS_CONSTANTS_Operators,INTERFACE_CONDITIONS_CONSTANTS 
  !>@}

  !> \addtogroup INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes INTERFACE_CONDITIONS_CONSTANTS::LinearityTypes
  !> \brief The interface condition linearity types
  !> \see INTERFACE_CONDITIONS_CONSTANTS,OPENCMISS_EquationsLinearityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_INTERFACE_CONDITION_LINEARITY_TYPES=3 !<The number of interface conditions linearity types defined. \see INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_LINEAR=1 !<The interface conditions are linear. \see INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_NONLINEAR=2 !<The interface conditions are non-linear. \see INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_NONLINEAR_BCS=3 !<The interface conditions have non-linear boundary conditions. \see INTERFACE_CONDITIONS_CONSTANTS_LinearityTypes,INTERFACE_CONDITIONS_CONSTANTS
  !>@}

  !> \addtogroup INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes INTERFACE_CONDITIONS_CONSTANTS::TimeDependenceTypes
  !> \brief The interface condition time dependence type parameters
  !> \see INTERFACE_CONDITIONS_CONSTANTS,OPENCMISS_EquationsTimeDependenceTypes
  !>@{
  INTEGER(INTG), PARAMETER :: NUMBER_OF_INTERFACE_CONDITION_TIME_TYPES=5 !<The number of interface conditions time dependence types defined. \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_STATIC=1 !<The interface conditions are static and have no time dependence. \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_QUASISTATIC=2 !<The interface conditions are quasi-static. \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_FIRST_ORDER_DYNAMIC=3 !<The interface conditions are first order dynamic. \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,INTERFACE_CONDITIONS_CONSTANTS
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_SECOND_ORDER_DYNAMIC=4 !<The interface conditions are a second order dynamic. \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_TIME_STEPPING=5 !<The interface conditions are for time stepping. \see INTERFACE_CONDITIONS_CONSTANTS_TimeDependenceTypes,EQUATIONS_ROUTINES
  !>@}
  
  !> \addtogroup INTERFACE_CONDITIONS_IntegrationType INTERFACE_CONDITIONS_CONSTANTS::IntegrationType
  !> \brief Interface condition IntegrationType.
  !> \see INTERFACE_CONDITIONS_CONSTANTS
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_GAUSS_INTEGRATION=1 !<Gauss points integration type, i.e. Loop over element Gauss points and sum up their contribution. \see INTERFACE_CONDITIONS_CONSTANTS_IntegrationType,INTERFACE_CONDITIONS_CONSTANTS 
  INTEGER(INTG), PARAMETER :: INTERFACE_CONDITION_DATA_POINTS_INTEGRATION=2 !< Data points integration type i.e. Loop over data points and  sum up their contribution.\see INTERFACE_CONDITIONS_CONSTANTS_IntegrationType,INTERFACE_CONDITIONS_CONSTANTS 
  !>@}

END MODULE INTERFACE_CONDITIONS_CONSTANTS
