!> \file
!> \author Caton Little
!> \brief This module handles non-IO FieldML logic.
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

!> Types for FieldML

MODULE FIELDML_TYPES

  USE KINDS
  USE TYPES

  IMPLICIT NONE

  !
  !================================================================================================================================
  !
  ! FieldML types

  !<Contains information on the current FieldML parsing state.
  TYPE FIELDML_IO_TYPE
    INTEGER(INTG) :: FML_HANDLE !<The FieldML session handle.
    INTEGER(INTG) :: NODES_HANDLE !<The FieldML global nodes type handle.
    INTEGER(INTG) :: NODES_ARGUMENT_HANDLE !<The FieldML global nodes argument handle.
    INTEGER(INTG) :: MESH_HANDLE !<The FieldML mesh type handle.
    INTEGER(INTG) :: ELEMENTS_HANDLE !<The FieldML mesh elements type handle.
    INTEGER(INTG) :: ELEMENTS_ARGUMENT_HANDLE !<The FieldML mesh elements argument handle.
    INTEGER(INTG) :: XI_HANDLE !<The FieldML mesh chart type handle.
    INTEGER(INTG) :: XI_ARGUMENT_HANDLE !<The FieldML mesh chart argument handle.
    INTEGER(INTG) :: NODE_DOFS_HANDLE !<The FieldML nodal dofs evaluator handle.
    LOGICAL :: IS_OUT !< True if the state is being used for output, false otherwise.
!    INTEGER(C_INT) :: elementDofsHandle !<The FieldML element dofs evaluator handle.
!    INTEGER(C_INT) :: constantDofsHandle !<The FieldML constant dofs evaluator handle.
    TYPE(LIST_TYPE), POINTER :: COMPONENT_HANDLES
    TYPE(LIST_TYPE), POINTER :: BASIS_HANDLES
    TYPE(LIST_TYPE), POINTER :: BASIS_CONNECTIVITY_HANDLES
    TYPE(LIST_TYPE), POINTER :: BASIS_LAYOUT_HANDLES
  END TYPE FIELDML_IO_TYPE

  !Interfaces
  
  PUBLIC :: FIELDML_IO_TYPE

CONTAINS

  !
  !================================================================================================================================
  !

END MODULE FIELDML_TYPES
