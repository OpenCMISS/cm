!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all field related routines.
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

!> This module handles all field related routines.
MODULE FIELD_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE LISTS
  USE MESH_ROUTINES
  USE MPI
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup FIELD_ROUTINES_DependentTypes FIELD_ROUTINES::DependentTypes
  !> \brief Depedent field parameter types
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_INDEPENDENT_TYPE=1 !<Independent field type \see FIELD_ROUTINES_DependentTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEPENDENT_TYPE=2 !<Dependent field type \see FIELD_ROUTINES_DependentTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_DimensionTypes FIELD_ROUTINES::DimensionTypes
  !> \brief Field dimension parameter types
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_SCALAR_DIMENSION_TYPE=1 !<Scalar field \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_VECTOR_DIMENSION_TYPE=2 !<Vector field \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_TENSOR_DIMENSION_TYPE=3 !<Tensor field \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_FieldTypes FIELD_ROUTINES::FieldTypes
  !> \brief Field type parameters
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_GEOMETRIC_TYPE=1 !<Geometric field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_FIBRE_TYPE=2 !<Fibre field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GENERAL_TYPE=3 !<General field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MATERIAL_TYPE=4 !<Material field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_InterpolationTypes FIELD_ROUTINES::InterpolationTypes
  !> \brief Field interpolation parameters
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_CONSTANT_INTERPOLATION=1 !<Constant interpolation. One parameter for the field \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ELEMENT_BASED_INTERPOLATION=2 !<Element based interpolation. Parameters are different in each element \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_NODE_BASED_INTERPOLATION=3 !<Node based interpolation. Parameters are nodal based and a basis function is used \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GRID_POINT_BASED_INTERPOLATION=4 !<Grid point based interpolation. Parameters are different at each grid point \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GAUSS_POINT_BASED_INTERPOLATION=5 !<Gauss point based interpolation. Parameters are different at each Gauss point \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_VariableTypes FIELD_ROUTINES::VariableTypes
  !> \brief Field variable type parameters
  !> \see FIELD_ROUTINES
  !> \todo sort out variable access routines so that you are always accessing by variable type rather than variable number.
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_VARIABLE_TYPES=6 !<Number of different field variable types possible \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
  INTEGER(INTG), PARAMETER :: FIELD_U_VARIABLE_TYPE=1 !<Standard variable type i.e., u \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELUDELN_VARIABLE_TYPE=2 !<Normal derivative variable type i.e., du/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELUDELT_VARIABLE_TYPE=3 !<First time derivative variable type i.e., du/dt \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DEL2UDELT2_VARIABLE_TYPE=4 !<Second type derivative variable type i.e., d^2u/dt^2 \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_V_VARIABLE_TYPE=5 !<Second standard variable type i.e., v \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DELVDELN_VARIABLE_TYPE=6 !<Second normal variable type i.e., dv/dn \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_DofTypes FIELD_ROUTINES::DofTypes
  !> \brief Field dof type parameters
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_CONSTANT_DOF_TYPE=1 !<The dof is from a field variable component with constant interpolation \see FIELD_ROUTINES_DofTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ELEMENT_DOF_TYPE=2 !<The dof is from a field variable component with element based interpolation \see FIELD_ROUTINES_DofTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_NODE_DOF_TYPE=3 !<The dof is from a field variable component with node based interpolation \see FIELD_ROUTINES_DofTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GRID_POINT_DOF_TYPE=4 !<The dof is from a field variable component with grid point based interpolation \see FIELD_ROUTINES_DofTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_GAUSS_POINT_DOF_TYPE=5 !<The dof is from a field variable component with Gauss point based interpolation \see FIELD_ROUTINES_DofTypes,FIELD_ROUTINES
  !>@}

   
  !> \addtogroup FIELD_ROUTINES_DataTypes FIELD_ROUTINES::DataTypes
  !> \brief Field data types
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_INTG_TYPE=1 !<Integer field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_SP_TYPE=2 !<Single precision real field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_DP_TYPE=3 !<Double precision real field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_L_TYPE=4 !<Logical field data type \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_ParameterSetTypes FIELD_ROUTINES::ParameterSetTypes
  !> \brief Field parameter set type parameters \todo make program defined constants negative?
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NUMBER_OF_SET_TYPES=99 !<The maximum number of different parameter sets for a field \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_VALUES_SET_TYPE=1 !<The parameter set corresponding to the field values (at time T+DT for dynamic problems) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_BOUNDARY_CONDITIONS_SET_TYPE=2 !<The parameter set corresponding to the field boundary conditions \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VALUES_SET_TYPE=3 !<The parameter set corresponding to the field initial values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INCREMENTAL_VALUES_SET_TYPE=4 !<The parameter set corresponding to the field incremental values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ANALYTIC_VALUES_SET_TYPE=5 !<The parameter set corresponding to the analytic field values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VALUES_SET_TYPE=6 !<The parameter set corresponding to the previous field values (at time T) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE=7 !<The parameter set corresponding to the mean predicited values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_VELOCITY_VALUES_SET_TYPE=8 !<The parameter set corresponding to the velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_VELOCITY_SET_TYPE=9 !<The parameter set corresponding to the initial velocity values for dynamic problems. This is also the previous velocity values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_VELOCITY_SET_TYPE=9 !<The parameter set corresponding to the previous velocity values (at time T). This is also the initial velocity values for dynamic problems. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE=10 !<The parameter set corresponding to the mean predicited velocity values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ACCELERATION_VALUES_SET_TYPE=11 !<The parameter set corresponding to the acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_INITIAL_ACCELERATION_SET_TYPE=12 !<The parameter set corresponding to the initial acceleration values for dynamic problems. This is also the previous accelearation values \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_PREVIOUS_ACCELERATION_SET_TYPE=12 !<The parameter set corresponding to the previous acceleration values (at time T).This is also the initial acceleration values for dynamic problems. \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE=13 !<The parameter set corresponding to the mean predicited acceleration values (at time T+DT) \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
  !>@}

  !> \addtogroup FIELD_ROUTINES_ScalingTypes FIELD_ROUTINES::ScalingTypes
  !> \brief Field scaling type parameters
  !> \see FIELD_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: FIELD_NO_SCALING=0 !<The field is not scaled \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_UNIT_SCALING=1 !<The field has unit scaling \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ARC_LENGTH_SCALING=2 !<The field has arc length scaling \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_ARITHMETIC_MEAN_SCALING=3 !<The field has arithmetic mean of the arc length scaling \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
  INTEGER(INTG), PARAMETER :: FIELD_HARMONIC_MEAN_SCALING=4 !<The field has geometric mean of the arc length scaling \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  !>Adds the alpha times the parameter set values from one parameter set type to another parameter set type.
  INTERFACE FIELD_PARAMETER_SETS_ADD
    MODULE PROCEDURE FIELD_PARAMETER_SETS_ADD_DP
    MODULE PROCEDURE FIELD_PARAMETER_SETS_ADD_DP1
  END INTERFACE !FIELD_PARAMETER_SETS_ADD

  !>Adds the given value to the given parameter set for the constant of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_ADD_CONSTANT
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_CONSTANT_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_CONSTANT_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_CONSTANT_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_CONSTANT_L
  END INTERFACE !FIELD_PARAMETER_SET_ADD_CONSTANT

  !>Adds the given value to the given parameter set for a particular local dof of the field variable.
  INTERFACE FIELD_PARAMETER_SET_ADD_LOCAL_DOF
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L
  END INTERFACE !FIELD_PARAMETER_SET_ADD_LOCAL_DOF

  !>Adds the given value to the given parameter set for a particular user element of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_ADD_ELEMENT
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_ELEMENT_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_ELEMENT_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_ELEMENT_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_ELEMENT_L
  END INTERFACE !FIELD_PARAMETER_SET_ADD_ELEMENT
  
  !>Adds the given value to the given parameter set for a particular local element of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L
  END INTERFACE !FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT

  !>Adds the given value to the given parameter set for a particular user node and derivative of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_ADD_NODE
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_NODE_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_NODE_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_NODE_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_NODE_L
  END INTERFACE !FIELD_PARAMETER_SET_ADD_NODE
  
  !>Adds the given value to the given parameter set for a particular local node and derivative of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_ADD_LOCAL_NODE
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L
  END INTERFACE !FIELD_PARAMETER_SET_ADD_LOCAL_NODE

  !>Returns a pointer to the specified field parameter set array. The pointer must be restored with a call to FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_RESTORE call. Note: the values can be used for read operations but a FIELD_PARAMETER_SET_UPDATE call must be used to change any values.
  INTERFACE FIELD_PARAMETER_SET_DATA_GET
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_GET_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_GET_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_GET_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_GET_L
  END INTERFACE !FIELD_PARAMETER_SET_DATA_GET
 
  !>Restores the specified field variable parameter set array that was obtained with FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_GET.
  INTERFACE FIELD_PARAMETER_SET_DATA_RESTORE
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_RESTORE_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_RESTORE_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_RESTORE_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_DATA_RESTORE_L
  END INTERFACE !FIELD_PARAMETER_SET_DATA_GET
 
  !>Updates the given parameter set with the given value for the constant of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_UPDATE_CONSTANT
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_CONSTANT_L
  END INTERFACE !FIELD_PARAMETER_SET_UPDATE_CONSTANT

  !>Updates the given parameter set with the given value for a particular local dof of the field variable.
  INTERFACE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L
  END INTERFACE !FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF

  !>Updates the given parameter set with the given value for a particular user element of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_UPDATE_ELEMENT
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_ELEMENT_L
  END INTERFACE !FIELD_PARAMETER_SET_UPDATE_ELEMENT
  
  !>Updates the given parameter set with the given value for a particular local element of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L
  END INTERFACE !FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT

  !>Updates the given parameter set with the given value for a particular user node and derivative of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_UPDATE_NODE
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_NODE_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_NODE_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_NODE_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_NODE_L
  END INTERFACE !FIELD_PARAMETER_SET_UPDATE_NODE
  
  !>Updates the given parameter set with the given value for a particular local node and derivative of the field variable component.
  INTERFACE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP
    MODULE PROCEDURE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L
  END INTERFACE !FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE
  
  PUBLIC FIELD_INDEPENDENT_TYPE,FIELD_DEPENDENT_TYPE

  PUBLIC FIELD_SCALAR_DIMENSION_TYPE,FIELD_VECTOR_DIMENSION_TYPE

  PUBLIC FIELD_GEOMETRIC_TYPE,FIELD_FIBRE_TYPE,FIELD_GENERAL_TYPE,FIELD_MATERIAL_TYPE

  PUBLIC FIELD_CONSTANT_INTERPOLATION,FIELD_ELEMENT_BASED_INTERPOLATION,FIELD_NODE_BASED_INTERPOLATION, &
    & FIELD_GRID_POINT_BASED_INTERPOLATION,FIELD_GAUSS_POINT_BASED_INTERPOLATION

  PUBLIC FIELD_CONSTANT_DOF_TYPE,FIELD_ELEMENT_DOF_TYPE,FIELD_NODE_DOF_TYPE,FIELD_GRID_POINT_DOF_TYPE,FIELD_GAUSS_POINT_DOF_TYPE

  PUBLIC FIELD_NUMBER_OF_VARIABLE_TYPES,FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DELUDELT_VARIABLE_TYPE, &
    & FIELD_DEL2UDELT2_VARIABLE_TYPE

  PUBLIC FIELD_INTG_TYPE,FIELD_SP_TYPE,FIELD_DP_TYPE,FIELD_L_TYPE

  PUBLIC FIELD_VALUES_SET_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,FIELD_INITIAL_VALUES_SET_TYPE,FIELD_INCREMENTAL_VALUES_SET_TYPE, &
    & FIELD_ANALYTIC_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE, &
    & FIELD_VELOCITY_VALUES_SET_TYPE,FIELD_INITIAL_VELOCITY_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
    & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,FIELD_ACCELERATION_VALUES_SET_TYPE,FIELD_INITIAL_ACCELERATION_SET_TYPE, &
    & FIELD_PREVIOUS_ACCELERATION_SET_TYPE,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE

  PUBLIC FIELD_NO_SCALING,FIELD_UNIT_SCALING,FIELD_ARC_LENGTH_SCALING,FIELD_HARMONIC_MEAN_SCALING,FIELD_ARITHMETIC_MEAN_SCALING

  PUBLIC FIELD_COMPONENT_DOF_GET_CONSTANT,FIELD_COMPONENT_DOF_GET_USER_ELEMENT,FIELD_COMPONENT_DOF_GET_USER_NODE

  PUBLIC FIELD_COMPONENT_INTERPOLATION_CHECK,FIELD_COMPONENT_INTERPOLATION_GET,FIELD_COMPONENT_INTERPOLATION_SET, &
    & FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK

  PUBLIC FIELD_COMPONENT_MESH_COMPONENT_CHECK,FIELD_COMPONENT_MESH_COMPONENT_GET,FIELD_COMPONENT_MESH_COMPONENT_SET, &
    & FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK
  
  PUBLIC FIELD_COMPONENT_VALUES_INITIALISE

  PUBLIC FIELD_CREATE_FINISH,FIELD_CREATE_START,FIELD_DESTROY,FIELDS_FINALISE,FIELDS_INITIALISE

  PUBLIC FIELD_DATA_TYPE_CHECK,FIELD_DATA_TYPE_GET,FIELD_DATA_TYPE_SET,FIELD_DATA_TYPE_SET_AND_LOCK

  PUBLIC FIELD_DEPENDENT_TYPE_CHECK,FIELD_DEPENDENT_TYPE_GET,FIELD_DEPENDENT_TYPE_SET,FIELD_DEPENDENT_TYPE_SET_AND_LOCK

  PUBLIC FIELD_DIMENSION_CHECK,FIELD_DIMENSION_GET,FIELD_DIMENSION_SET,FIELD_DIMENSION_SET_AND_LOCK

  PUBLIC FIELD_GEOMETRIC_FIELD_GET,FIELD_GEOMETRIC_FIELD_SET,FIELD_GEOMETRIC_FIELD_SET_AND_LOCK

  PUBLIC FIELD_INTERPOLATE_GAUSS,FIELD_INTERPOLATE_XI

  PUBLIC FIELD_INTERPOLATED_POINT_METRICS_CALCULATE,FIELD_INTERPOLATED_POINT_METRICS_FINALISE, &
    & FIELD_INTERPOLATED_POINT_METRICS_INITIALISE,FIELD_INTERPOLATED_POINT_FINALISE,FIELD_INTERPOLATED_POINT_INITIALISE

  PUBLIC FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET,FIELD_INTERPOLATION_PARAMETERS_FINALISE, &
    & FIELD_INTERPOLATION_PARAMETERS_INITIALISE,FIELD_INTERPOLATION_PARAMETERS_LINE_GET, &
    & FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET,FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET

  PUBLIC FIELD_MESH_DECOMPOSITION_GET,FIELD_MESH_DECOMPOSITION_SET,FIELD_MESH_DECOMPOSITION_SET_AND_LOCK

  PUBLIC FIELD_NUMBER_OF_COMPONENTS_CHECK,FIELD_NUMBER_OF_COMPONENTS_GET,FIELD_NUMBER_OF_COMPONENTS_SET, &
    & FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK

  PUBLIC FIELD_NUMBER_OF_VARIABLES_CHECK,FIELD_NUMBER_OF_VARIABLES_GET,FIELD_NUMBER_OF_VARIABLES_SET, &
    & FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK
  
  PUBLIC FIELD_PARAMETER_SETS_ADD,FIELD_PARAMETER_SETS_COPY
    
  PUBLIC FIELD_PARAMETER_SET_ADD_CONSTANT,FIELD_PARAMETER_SET_ADD_LOCAL_DOF,FIELD_PARAMETER_SET_ADD_ELEMENT, &
    & FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT,FIELD_PARAMETER_SET_ADD_NODE,FIELD_PARAMETER_SET_ADD_LOCAL_NODE, &
    & FIELD_PARAMETER_SET_CREATE,FIELD_PARAMETER_SET_DATA_GET,FIELD_PARAMETER_SET_DATA_RESTORE, &
    & FIELD_PARAMETER_SET_UPDATE_FINISH,FIELD_PARAMETER_SET_UPDATE_START,FIELD_PARAMETER_SET_UPDATE_CONSTANT, &
    & FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF,FIELD_PARAMETER_SET_UPDATE_ELEMENT,FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT, &
    & FIELD_PARAMETER_SET_UPDATE_NODE,FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE,FIELD_PARAMETER_SET_VECTOR_GET

  PUBLIC FIELD_SCALING_TYPE_CHECK,FIELD_SCALING_TYPE_GET,FIELD_SCALING_TYPE_SET,FIELD_SCALING_TYPE_SET_AND_LOCK

  PUBLIC FIELD_TYPE_CHECK,FIELD_TYPE_GET,FIELD_TYPE_SET,FIELD_TYPE_SET_AND_LOCK
  
  PUBLIC FIELD_USER_NUMBER_FIND

  PUBLIC FIELD_VARIABLE_GET
  
  PUBLIC FIELD_VARIABLE_TYPES_CHECK,FIELD_VARIABLE_TYPES_GET,FIELD_VARIABLE_TYPES_SET,FIELD_VARIABLE_TYPES_SET_AND_LOCK
   
CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks the interpolation type for a field variable component.
  SUBROUTINE FIELD_COMPONENT_INTERPOLATION_CHECK(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the interpolation for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type of the field variable component to check \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number of the field variable component to check
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_TYPE !<The interpolation type of the field variable component to check \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_INTERPOLATION_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              SELECT CASE(INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE/=FIELD_CONSTANT_INTERPOLATION) THEN
                  LOCAL_ERROR="Invalid interpolation type. The interpolation type for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " which is not constant interpolation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE/=FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                  LOCAL_ERROR="Invalid interpolation type. The interpolation type for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " which is not element based interpolation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE/=FIELD_NODE_BASED_INTERPOLATION) THEN
                  LOCAL_ERROR="Invalid interpolation type. The interpolation type for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " which is not node based interpolation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE/=FIELD_GRID_POINT_BASED_INTERPOLATION) THEN
                  LOCAL_ERROR="Invalid interpolation type. The interpolation type for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " which is not grid point based interpolation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE/=FIELD_GAUSS_POINT_BASED_INTERPOLATION) THEN
                  LOCAL_ERROR="Invalid interpolation type. The interpolation type for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                    & " which is not Gauss point based interpolation."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE DEFAULT
                LOCAL_ERROR="The specified interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_CHECK")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_INTERPOLATION_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_INTERPOLATION_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the interpolation type for a field variable component identified by a pointer.
  SUBROUTINE FIELD_COMPONENT_INTERPOLATION_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the interpolation for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type of the field variable component to get \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(OUT) :: INTERPOLATION_TYPE !<On return, the interpolation type of the field variable component \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_INTERPOLATION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              INTERPOLATION_TYPE=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_GET")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_INTERPOLATION_GET",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_GET")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_INTERPOLATION_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type for a field variable component.
  SUBROUTINE FIELD_COMPONENT_INTERPOLATION_SET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the interpolation for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type of the field variable component to set \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_TYPE !<The interpolation type to set \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_INTERPOLATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN          
          IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            IF(ANY(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES==VARIABLE_TYPE)) THEN
              IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)) THEN
                IF(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(COMPONENT_NUMBER,VARIABLE_TYPE)) THEN
                  LOCAL_ERROR="The interpolation type has been locked for component number "// &
                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  SELECT CASE(INTERPOLATION_TYPE)                
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE)=INTERPOLATION_TYPE
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE)=INTERPOLATION_TYPE
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE)=INTERPOLATION_TYPE
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE)=INTERPOLATION_TYPE
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE)=INTERPOLATION_TYPE
                  CASE DEFAULT
                    LOCAL_ERROR="The specified interpolation type of "// &
                      & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDIF
              ELSE
                LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE),"*",ERR,ERROR))// &
                  & " components."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_SET")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_INTERPOLATION_SET",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_SET")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_INTERPOLATION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the interpolation type for a field variable component and locks so that no further changes can be made.
  SUBROUTINE FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the interpolation for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type of the field variable component to set \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number of the field variable component to set
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_TYPE !<The interpolation type to set \see FIELD_ROUTINES_InterpolationTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_COMPONENT_INTERPOLATION_SET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(COMPONENT_NUMBER,VARIABLE_TYPE)=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Returns the dof numbers for a field variable component that corresponds to the specified constant
  SUBROUTINE FIELD_COMPONENT_DOF_GET_CONSTANT(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,LOCAL_DOF,GLOBAL_DOF,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the dof for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: LOCAL_DOF !<On exit, the local dof corresponding to the constant
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_DOF !<On exit, the global dof corresponding to the constant
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FIELD_COMPONENT_DOF_GET_CONSTANT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                IF(ASSOCIATED(FIELD_VARIABLE%DOMAIN_MAPPING)) THEN
                  LOCAL_DOF=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                  GLOBAL_DOF=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(LOCAL_DOF)
                ELSE
                  LOCAL_ERROR="The field variable domain mapping is not associated for variable type "// &
                    & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by constant for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by constant for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by constant for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by constant for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                
              CASE DEFAULT
                LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                  & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                  & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_COMPONENT_DOF_GET_CONSTANT")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_DOF_GET_CONSTANT",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_DOF_GET_CONSTANT")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_DOF_GET_CONSTANT

  !
  !================================================================================================================================
  !

  !>Returns the dof numbers for a field component that corresponds to the specified user element.
  SUBROUTINE FIELD_COMPONENT_DOF_GET_USER_ELEMENT(FIELD,VARIABLE_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER,LOCAL_DOF, &
    & GLOBAL_DOF,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the dof for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to get the dof for
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: LOCAL_DOF !<On exit, the local dof corresponding to the user element
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_DOF !<On exit, the global dof corresponding to the user element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FIELD_COMPONENT_DOF_GET_USER_ELEMENT",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user element for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                DECOMPOSITION=>FIELD%DECOMPOSITION
                IF(ASSOCIATED(DECOMPOSITION)) THEN
                  MESH=>DECOMPOSITION%MESH
                  IF(ASSOCIATED(MESH)) THEN
                    MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                    CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                      & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                    IF(MESH_ELEMENT_EXISTS) THEN
                      DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                      IF(ASSOCIATED(DOMAIN)) THEN
                        DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                        IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                          ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                          IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                            CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER,LOCAL_EXISTS, &
                              & DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                            IF(LOCAL_EXISTS) THEN
                              IF(ASSOCIATED(FIELD_VARIABLE%DOMAIN_MAPPING)) THEN
                                LOCAL_DOF=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                  & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                GLOBAL_DOF=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(LOCAL_DOF)
                              ELSE
                                CALL FLAG_ERROR("The field variable domain mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="The specified user element number of "// &
                                & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                & " is not part of this local domain."                                  
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified user element number of "// &
                        & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                        " does not exist in mesh component number "// &
                        & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                        & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user element for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user element for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user element for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                  & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                  & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_COMPONENT_DOF_GET_USER_ELEMENT")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_DOF_GET_USER_ELEMENT",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_DOF_GET_USER_ELEMENT")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_DOF_GET_USER_ELEMENT

  !
  !================================================================================================================================
  !

  !>Returns the dof numbers for a field component that corresponds to the specified user node and derivative.
  SUBROUTINE FIELD_COMPONENT_DOF_GET_USER_NODE(FIELD,VARIABLE_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER,COMPONENT_NUMBER, &
    & LOCAL_DOF,GLOBAL_DOF,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the dof for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the dof for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative number to get the dof for
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to get the dof for
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to get the dof for
    INTEGER(INTG), INTENT(OUT) :: LOCAL_DOF !<On exit, the local dof corresponding to the user node
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_DOF !<On exit, the global dof corresponding to the user node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FIELD_COMPONENT_DOF_GET_USER_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user node for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user node for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                DECOMPOSITION=>FIELD%DECOMPOSITION
                IF(ASSOCIATED(DECOMPOSITION)) THEN
                  MESH=>DECOMPOSITION%MESH
                  IF(ASSOCIATED(MESH)) THEN
                    MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                    CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                      & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                    IF(MESH_NODE_EXISTS) THEN
                      DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                      IF(ASSOCIATED(DOMAIN)) THEN
                        DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                        IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                          NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                          IF(ASSOCIATED(NODES_MAPPING)) THEN
                            CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                              & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                            IF(LOCAL_EXISTS) THEN
                              IF(ASSOCIATED(FIELD_VARIABLE%DOMAIN_MAPPING)) THEN
                                IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                  & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                  LOCAL_DOF=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                    & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                  GLOBAL_DOF=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(LOCAL_DOF)
                                ELSE
                                  LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                    & " is invalid for user node number "// &
                                    & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                    & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                    & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                    & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                    & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("The field variable domain mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="The specified user node number of "// &
                                & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                & " is not part of this local domain."                                  
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The specified user node number of "// &
                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                        " does not exist in mesh component number "// &
                        & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                        & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                ENDIF
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user node for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                LOCAL_ERROR="Can not get the dof by user node for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                  & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                  & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_DOF_GET_USER_NODE")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_DOF_GET_USER_NODE",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_DOF_GET_USER_NODE")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_DOF_GET_USER_NODE

  !
  !================================================================================================================================
  !

  !>Check the mesh component number for a field variable component.
  SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_CHECK(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,MESH_COMPONENT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the mesh component for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to check the field variable component for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to check the field variable component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT !<The mesh component to check for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_MESH_COMPONENT_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER/=MESH_COMPONENT) THEN
                LOCAL_ERROR="Invalid mesh component number. The mesh component number for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                  & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " which is does correspond to the specified mesh component number of "// &
                  & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The field variable type must be > 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_CHECK")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_MESH_COMPONENT_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_CHECK
  
  !
  !================================================================================================================================
  !

  !>Gets the mesh component number for a field variable component.
  SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_GET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,MESH_COMPONENT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the mesh component for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the field variable component for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to get the field variable component for
    INTEGER(INTG), INTENT(OUT) :: MESH_COMPONENT !<On return, the mesh component to get for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_MESH_COMPONENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The field variable type must be > 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_GET")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_MESH_COMPONENT_GET",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_GET")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number for a field variable component.
  SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_SET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to set the mesh component for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component to set for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(MESH_TYPE), POINTER :: MESH   
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        DECOMPOSITION=>FIELD%DECOMPOSITION
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          MESH=>DECOMPOSITION%MESH
          IF(ASSOCIATED(MESH)) THEN
            IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
              IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                IF(ANY(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES==VARIABLE_TYPE)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)) THEN
                    IF(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(COMPONENT_NUMBER,VARIABLE_TYPE)) THEN
                      LOCAL_ERROR="The mesh component has been locked for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ELSE
                      SELECT CASE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE))
                      CASE(FIELD_CONSTANT_INTERPOLATION)
                        LOCAL_ERROR="Can not set a mesh component for field component number "// &
                          & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has constant interpolation."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      CASE(FIELD_ELEMENT_BASED_INTERPOLATION,FIELD_NODE_BASED_INTERPOLATION,FIELD_GRID_POINT_BASED_INTERPOLATION, &
                        & FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                        IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
                          FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(COMPONENT_NUMBER,VARIABLE_TYPE)=MESH_COMPONENT_NUMBER
                        ELSE
                          LOCAL_ERROR="Mesh component number "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid. The component number must be between 1 and "// &
                            & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                            & " for mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      CASE DEFAULT
                        LOCAL_ERROR="The interpolation type "//TRIM(NUMBER_TO_VSTRING(FIELD%CREATE_VALUES_CACHE% &
                          & INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE),"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE),"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field variable type must be > 1 and <= "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF              
            ELSE
              CALL FLAG_ERROR("Field create values cache is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The decomposition mesh is not associated for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The decomposition is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_SET")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_SET")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number for a field variable component and locks it so that no further changes can be made.
  SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to set the mesh component for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component to set for the specified field variable component
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(FIELD,VARIABLE_TYPE,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(COMPONENT_NUMBER,VARIABLE_TYPE)=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Initialises the values of parameter set of a field variable component to a constant value
  SUBROUTINE FIELD_COMPONENT_VALUES_INITIALISE(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the values for 
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to intiialise \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to initialise \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number to initialise
    REAL(DP), INTENT(IN) :: VALUE !<The constant value to initialise the parameter set for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,derivative_idx,field_dof,node_idx,partial_deriv_idx
    REAL(DP), POINTER :: FIELD_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: COMPONENT_DOMAIN
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: FIELD_PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
   
    CALL ENTERS("FIELD_COMPONENT_VALUES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        !Check the variable type
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            !Check the component number
            IF(COMPONENT_NUMBER>0.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              !Check the from set type input
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                FIELD_PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(FIELD_PARAMETER_SET)) THEN
                  !Get the parameters values
                  CALL DISTRIBUTED_VECTOR_DATA_GET(FIELD_PARAMETER_SET%PARAMETERS,FIELD_PARAMETERS,ERR,ERROR,*999)
                  !Set the field components to give a constant value. Note that as the value is constant we can set the ghost dofs
                  !and not worry about updating the field parameter set.
                  SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                  CASE(FIELD_CONSTANT_INTERPOLATION)
                    field_dof=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                    FIELD_PARAMETERS(field_dof)=VALUE
                  CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                    COMPONENT_DOMAIN=>FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN
                    IF(ASSOCIATED(COMPONENT_DOMAIN)) THEN
                      DOMAIN_TOPOLOGY=>COMPONENT_DOMAIN%TOPOLOGY
                      IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DO element_idx=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            field_dof=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                              & ELEMENT_PARAM2DOF_MAP(element_idx)
                            FIELD_PARAMETERS(field_dof)=VALUE
                          ENDDO !element_idx
                        ELSE
                          CALL FLAG_ERROR("Domain topology elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                    ENDIF
                  CASE(FIELD_NODE_BASED_INTERPOLATION)
                    COMPONENT_DOMAIN=>FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN
                    IF(ASSOCIATED(COMPONENT_DOMAIN)) THEN
                      DOMAIN_TOPOLOGY=>COMPONENT_DOMAIN%TOPOLOGY
                      IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                        DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                          DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                            DO derivative_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                              field_dof=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP(derivative_idx,node_idx)
                              partial_deriv_idx=DOMAIN_NODES%NODES(node_idx)%PARTIAL_DERIVATIVE_INDEX(derivative_idx)
                              SELECT CASE(partial_deriv_idx)
                              CASE(NO_PART_DERIV)
                                FIELD_PARAMETERS(field_dof)=VALUE
                              CASE(PART_DERIV_S1)
                                FIELD_PARAMETERS(field_dof)=1.0_DP
                              CASE(PART_DERIV_S1_S1)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE(PART_DERIV_S2)
                                FIELD_PARAMETERS(field_dof)=1.0_DP
                              CASE(PART_DERIV_S2_S2)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE(PART_DERIV_S1_S2)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE(PART_DERIV_S3)
                                FIELD_PARAMETERS(field_dof)=1.0_DP
                              CASE(PART_DERIV_S3_S3)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE(PART_DERIV_S1_S3)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE(PART_DERIV_S2_S3)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE(PART_DERIV_S1_S2_S3)
                                FIELD_PARAMETERS(field_dof)=0.0_DP
                              CASE DEFAULT
                                LOCAL_ERROR="The partial derivative index of "// &
                                  & TRIM(NUMBER_TO_VSTRING(partial_deriv_idx,"*",ERR,ERROR))//" for node number "// &
                                  & TRIM(NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR))//" and derivative number "// &
                                  & TRIM(NUMBER_TO_VSTRING(derivative_idx,"*",ERR,ERROR))//" is invalid."
                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                            ENDDO !derivative_idx
                          ENDDO !node_idx
                        ELSE
                          CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                    ENDIF                    
                  CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                      & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for component number "// &
                      & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                  !Restore the  parameter set 
                  CALL DISTRIBUTED_VECTOR_DATA_RESTORE(FIELD_PARAMETER_SET%PARAMETERS,FIELD_PARAMETERS,ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not be created on variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// & 
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable component number of "// &
                & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" is invalid for a variable type of "//&
                & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" on field number "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//". The number of components must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is not defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The field variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_COMPONENT_VALUES_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_COMPONENT_VALUES_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_COMPONENT_VALUES_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_COMPONENT_VALUES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Checks the data type for a field variable.
  SUBROUTINE FIELD_DATA_TYPE_CHECK(FIELD,VARIABLE_TYPE,DATA_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the data type for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to check \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type of the field variable to check \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DATA_TYPE_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            SELECT CASE(DATA_TYPE)              
            CASE(FIELD_INTG_TYPE)
              IF(FIELD_VARIABLE%DATA_TYPE/=FIELD_INTG_TYPE) THEN
                LOCAL_ERROR="Invalid data type. The data type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " which is not an integer data type."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE(FIELD_SP_TYPE)
              IF(FIELD_VARIABLE%DATA_TYPE/=FIELD_SP_TYPE) THEN
                LOCAL_ERROR="Invalid data type. The data type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " which is not a single precision data type."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE(FIELD_DP_TYPE)
              IF(FIELD_VARIABLE%DATA_TYPE/=FIELD_DP_TYPE) THEN
                LOCAL_ERROR="Invalid data type. The data type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " which is not a double precision data type."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
             CASE(FIELD_L_TYPE)
               IF(FIELD_VARIABLE%DATA_TYPE/=FIELD_L_TYPE) THEN
                 LOCAL_ERROR="Invalid data type. The data type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " which is not a logical data type."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The specified data type of "//TRIM(NUMBER_TO_VSTRING(DATA_TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DATA_TYPE_CHECK")
    RETURN
999 CALL ERRORS("FIELD_DATA_TYPE_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_DATA_TYPE_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_DATA_TYPE_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the data type for a field variable.
  SUBROUTINE FIELD_DATA_TYPE_GET(FIELD,VARIABLE_TYPE,DATA_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the data type for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: DATA_TYPE !<On return, the data type of the field variable \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DATA_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            DATA_TYPE=FIELD_VARIABLE%DATA_TYPE
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DATA_TYPE_GET")
    RETURN
999 CALL ERRORS("FIELD_DATA_TYPE_GET",ERR,ERROR)
    CALL EXITS("FIELD_DATA_TYPE_GET")
    RETURN 1
  END SUBROUTINE FIELD_DATA_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable.
  SUBROUTINE FIELD_DATA_TYPE_SET(FIELD,VARIABLE_TYPE,DATA_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the interpolation for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type of the field variable component to set \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type to set \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DATA_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN          
          IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN            
            IF(ANY(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES==VARIABLE_TYPE)) THEN
              IF(FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED(VARIABLE_TYPE)) THEN
                LOCAL_ERROR="The data type has been locked for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                SELECT CASE(DATA_TYPE)                
                CASE(FIELD_INTG_TYPE)
                  FIELD%CREATE_VALUES_CACHE%DATA_TYPES(VARIABLE_TYPE)=FIELD_INTG_TYPE
                CASE(FIELD_SP_TYPE)
                  FIELD%CREATE_VALUES_CACHE%DATA_TYPES(VARIABLE_TYPE)=FIELD_SP_TYPE
                CASE(FIELD_DP_TYPE)
                  FIELD%CREATE_VALUES_CACHE%DATA_TYPES(VARIABLE_TYPE)=FIELD_DP_TYPE
                CASE(FIELD_L_TYPE)
                  FIELD%CREATE_VALUES_CACHE%DATA_TYPES(VARIABLE_TYPE)=FIELD_L_TYPE
                CASE DEFAULT
                  LOCAL_ERROR="The specified data type of "//TRIM(NUMBER_TO_VSTRING(DATA_TYPE,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Field create values cache is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DATA_TYPE_SET")
    RETURN
999 CALL ERRORS("FIELD_DATA_TYPE_SET",ERR,ERROR)
    CALL EXITS("FIELD_DATA_TYPE_SET")
    RETURN 1
  END SUBROUTINE FIELD_DATA_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the data type for a field variable and locks it so that no further changes can be made.
  SUBROUTINE FIELD_DATA_TYPE_SET_AND_LOCK(FIELD,VARIABLE_TYPE,DATA_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the interpolation for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type of the field variable component to set \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DATA_TYPE !<The data type to set \see FIELD_ROUTINES_DataTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DATA_TYPE_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_DATA_TYPE_SET(FIELD,VARIABLE_TYPE,DATA_TYPE,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED(VARIABLE_TYPE)=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("FIELD_DATA_TYPE_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_DATA_TYPE_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_DATA_TYPE_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_DATA_TYPE_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Finalises a field variable component and deallocates all memory.
  SUBROUTINE FIELD_VARIABLE_COMPONENT_FINALISE(FIELD_VARIABLE_COMPONENT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE) :: FIELD_VARIABLE_COMPONENT !<The field variable component to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_VARIABLE_COMPONENT_FINALISE",ERR,ERROR,*999)

    CALL FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE(FIELD_VARIABLE_COMPONENT,ERR,ERROR,*999)

    CALL EXITS("FIELD_VARIABLE_COMPONENT_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_COMPONENT_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_COMPONENT_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_COMPONENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a field variable component.
  SUBROUTINE FIELD_VARIABLE_COMPONENT_INITIALISE(FIELD_VARIABLE,COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<A pointer to the field variable to initialise the component for
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field component number of the field variable component
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: COMP_NUMBER,DUMMY_ERR,ne,VARIABLE_TYPE
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION    
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_COMPONENT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
      FIELD=>FIELD_VARIABLE%FIELD
      IF(ASSOCIATED(FIELD)) THEN
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          VARIABLE_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          IF(ALLOCATED(FIELD_VARIABLE%COMPONENTS)) THEN
            IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
              FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%COMPONENT_NUMBER=COMPONENT_NUMBER
              FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%FIELD_VARIABLE=>FIELD_VARIABLE
              DECOMPOSITION=>FIELD%DECOMPOSITION
              IF(ASSOCIATED(DECOMPOSITION)) THEN
                MESH=>DECOMPOSITION%MESH
                IF(ASSOCIATED(MESH)) THEN
                  COMP_NUMBER=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(COMPONENT_NUMBER,VARIABLE_TYPE)
                  IF(COMP_NUMBER>0.AND.COMP_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
                    FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER=COMP_NUMBER
                    FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN=>DECOMPOSITION%DOMAIN(COMP_NUMBER)%PTR
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%DOMAIN
                    IF(.NOT.ASSOCIATED(DOMAIN)) THEN
                      LOCAL_ERROR="Field component "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                        & " does not have a domain associated."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The mesh component number of "//TRIM(NUMBER_TO_VSTRING(COMP_NUMBER,"*",ERR,ERROR))// &
                      & " for field component "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. The component number must be between 1 and "// &
                      & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Decomposition mesh is not associated for field number "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="Decomposition is not associated for field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                  
              ENDIF
              FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE= &
                FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(COMPONENT_NUMBER,VARIABLE_TYPE)
              SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS=1
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS=1
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS=-1
                DO ne=1,DOMAIN%TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                  BASIS=>DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
                  IF(BASIS%NUMBER_OF_ELEMENT_PARAMETERS>FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                    & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS) FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                    & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                ENDDO !ne
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                  & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                  & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is invalid."
              END SELECT
              CALL FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER), &
                & ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Field variable components have not been allocated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Field create values cache is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field variable field is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field variable is is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_COMPONENT_INITIALISE")
    RETURN
999 CALL FIELD_VARIABLE_COMPONENT_FINALISE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER),DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_VARIABLE_COMPONENT_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_COMPONENT_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_COMPONENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises a field variable component parameter to dof map and deallocates all memory.
  SUBROUTINE FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE(FIELD_VARIABLE_COMPONENT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE) :: FIELD_VARIABLE_COMPONENT !<The field variable component to finialise the parameter to dof map for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP))  &
      & DEALLOCATE(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP)
    IF(ALLOCATED(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP))  &
      & DEALLOCATE(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP)
    IF(ALLOCATED(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%GRID_POINT_PARAM2DOF_MAP))  &
      & DEALLOCATE(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%GRID_POINT_PARAM2DOF_MAP)
    IF(ALLOCATED(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP))  &
      & DEALLOCATE(FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP)
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_GRID_POINT_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_GAUSS_POINT_PARAMETERS=0

    CALL EXITS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a field variable component parameter to dof map.
  SUBROUTINE FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE(FIELD_VARIABLE_COMPONENT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE) :: FIELD_VARIABLE_COMPONENT !<The field variable component to initialise the parameter to dof map for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE",ERR,ERROR,*999)

    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_GRID_POINT_PARAMETERS=0
    FIELD_VARIABLE_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_GAUSS_POINT_PARAMETERS=0

    CALL EXITS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_COMPONENT_PARAM_TO_DOF_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the field variable components for a field variable and deallocates all memory.
  SUBROUTINE FIELD_VARIABLE_COMPONENTS_FINALISE(FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE) :: FIELD_VARIABLE !<The field variable to finalise the field variable components for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    CALL ENTERS("FIELD_VARIABLE_COMPONENTS_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(FIELD_VARIABLE%COMPONENTS)) THEN
      DO component_idx=1,SIZE(FIELD_VARIABLE%COMPONENTS,1)
        CALL FIELD_VARIABLE_COMPONENT_FINALISE(FIELD_VARIABLE%COMPONENTS(component_idx),ERR,ERROR,*999)
      ENDDO !component_idx
      DEALLOCATE(FIELD_VARIABLE%COMPONENTS)
    ENDIF
    FIELD_VARIABLE%NUMBER_OF_COMPONENTS=0

    CALL EXITS("FIELD_VARIABLE_COMPONENTS_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_COMPONENTS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_COMPONENTS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_COMPONENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the field components.
  SUBROUTINE FIELD_VARIABLE_COMPONENTS_INITIALISE(FIELD,VARIABLE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the field variable components for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to initialise the field variable components for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_COMPONENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(ALLOCATED(FIELD_VARIABLE%COMPONENTS)) THEN
              CALL FLAG_ERROR("Field variable already has allocated components.",ERR,ERROR,*999)
            ELSE
              ALLOCATE(FIELD_VARIABLE%COMPONENTS(FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable components.",ERR,ERROR,*999)
              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                CALL FIELD_VARIABLE_COMPONENT_INITIALISE(FIELD_VARIABLE,component_idx,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field create values cache is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_COMPONENTS_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_COMPONENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_COMPONENTS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_COMPONENTS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a field.
  SUBROUTINE FIELD_CREATE_FINISH(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has already been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        !Check field has a decomposition associated
        IF(ASSOCIATED(FIELD%DECOMPOSITION)) THEN
          !Initialise the components
          CALL FIELD_VARIABLES_INITIALISE(FIELD,ERR,ERROR,*999)
          IF(ASSOCIATED(FIELD%GEOMETRIC_FIELD)) THEN
            CALL FIELD_CREATE_VALUES_CACHE_FINALISE(FIELD%CREATE_VALUES_CACHE,ERR,ERROR,*999)
            FIELD%FIELD_FINISHED=.TRUE.
            !Calculate dof mappings
            CALL FIELD_MAPPINGS_CALCULATE(FIELD,ERR,ERROR,*999)
            !Set up the geometric parameters
            CALL FIELD_GEOMETRIC_PARAMETERS_INITIALISE(FIELD,ERR,ERROR,*999)
            !Initialise the scalings
            CALL FIELD_SCALINGS_INITIALISE(FIELD,ERR,ERROR,*999)
            !Initialise the field parameter sets 
            CALL FIELD_PARAMETER_SETS_INITIALISE(FIELD,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Field does not have a geometric field associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Field does not have a mesh decomposition associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field number : ",FIELD%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number       = ",FIELD%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dependent type      = ",FIELD%DEPENDENT_TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Field type          = ",FIELD%TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of variables = ",FIELD%NUMBER_OF_VARIABLES,ERR,ERROR,*999)
   ENDIF

    CALL EXITS("FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("FIELD_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("FIELD_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE FIELD_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a field defined by a user number in the specified region. 
  !>Default values set for the FIELD's attributes are:
  !>- DEPENDENT_TYPE: 1 (FIELD_INDEPENDENT_TYPE)
  !>- DIMENSION: 2 (FIELD_VECTOR_DIMENSION_TYPE)
  !>- TYPE: 1 (FIELD_GEOMETRIC_TYPE)
  !>- NUMBER_OF_VARIABLES: 1
  !>- GEOMETRIC_FIELD: itself
  !>- SCALINGS%SCALING_TYPE: 3 (FIELD_ARITHMETIC_MEAN_SCALING)
  !>\todo Add in FIELD_INITIALISE
  SUBROUTINE FIELD_CREATE_START(USER_NUMBER,REGION,FIELD,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number for the field
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region in which to create the field
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On return, a pointer to the field being created
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_no
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_TYPE), POINTER :: NEW_FIELD
    TYPE(FIELD_PTR_TYPE), POINTER :: NEW_FIELDS(:)

    NULLIFY(NEW_FIELD)
    NULLIFY(NEW_FIELDS)

    CALL ENTERS("FIELD_CREATE_START",ERR,ERROR,*999)

    NULLIFY(FIELD)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%FIELDS)) THEN
        CALL FIELD_USER_NUMBER_FIND(USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
        IF(ASSOCIATED(FIELD)) THEN
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ELSE
          ALLOCATE(NEW_FIELD,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new field.",ERR,ERROR,*999)
          !Set default field properties
          CALL FIELD_INITIALISE(NEW_FIELD,ERR,ERROR,*999)
          NEW_FIELD%GLOBAL_NUMBER=REGION%FIELDS%NUMBER_OF_FIELDS+1
          NEW_FIELD%USER_NUMBER=USER_NUMBER
          NEW_FIELD%FIELDS=>REGION%FIELDS
          NEW_FIELD%REGION=>REGION
          NEW_FIELD%GEOMETRIC_FIELD=>NEW_FIELD
          NEW_FIELD%NUMBER_OF_VARIABLES=1
          NEW_FIELD%SCALINGS%SCALING_TYPE=FIELD_ARITHMETIC_MEAN_SCALING
          NEW_FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES=0
          CALL FIELD_CREATE_VALUES_CACHE_INITIALISE(NEW_FIELD,ERR,ERROR,*999)
          !Add new field into list of fields in the region
          ALLOCATE(NEW_FIELDS(REGION%FIELDS%NUMBER_OF_FIELDS+1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new fields.",ERR,ERROR,*999)
          DO field_no=1,REGION%FIELDS%NUMBER_OF_FIELDS
            NEW_FIELDS(field_no)%PTR=>REGION%FIELDS%FIELDS(field_no)%PTR
          ENDDO !field_no
          NEW_FIELDS(REGION%FIELDS%NUMBER_OF_FIELDS+1)%PTR=>NEW_FIELD
          IF(ASSOCIATED(REGION%FIELDS%FIELDS)) DEALLOCATE(REGION%FIELDS%FIELDS)
          REGION%FIELDS%FIELDS=>NEW_FIELDS
          REGION%FIELDS%NUMBER_OF_FIELDS=REGION%FIELDS%NUMBER_OF_FIELDS+1
          FIELD=>NEW_FIELD
        ENDIF
      ELSE
        LOCAL_ERROR="The fields on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_FIELD)) DEALLOCATE(NEW_FIELD)
    IF(ASSOCIATED(NEW_FIELDS)) DEALLOCATE(NEW_FIELDS)
998 NULLIFY(FIELD)
    CALL ERRORS("FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("FIELD_CREATE_START")
    RETURN 1
  END SUBROUTINE FIELD_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalise the create values cache for a field.
  SUBROUTINE FIELD_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%DIMENSION)) DEALLOCATE(CREATE_VALUES_CACHE%DIMENSION)
      IF(ALLOCATED(CREATE_VALUES_CACHE%DIMENSION_LOCKED)) DEALLOCATE(CREATE_VALUES_CACHE%DIMENSION_LOCKED)
      IF(ALLOCATED(CREATE_VALUES_CACHE%DATA_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%DATA_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%DATA_TYPES_LOCKED)) DEALLOCATE(CREATE_VALUES_CACHE%DATA_TYPES_LOCKED)
      IF(ALLOCATED(CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS)) DEALLOCATE(CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED)) DEALLOCATE(CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED)
      IF(ALLOCATED(CREATE_VALUES_CACHE%INTERPOLATION_TYPE)) DEALLOCATE(CREATE_VALUES_CACHE%INTERPOLATION_TYPE)
      IF(ALLOCATED(CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED)) DEALLOCATE(CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER)) DEALLOCATE(CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED)) DEALLOCATE(CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
 
    CALL EXITS("FIELD_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the create values cache for a field.
  SUBROUTINE FIELD_CREATE_VALUES_CACHE_INITIALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the create values cache for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,NUMBER_OF_COMPONENTS,component_idx,variable_idx
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("FIELD_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(FIELD%CREATE_VALUES_CACHE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate create values cache.",ERR,ERROR,*999)
        SELECT CASE(FIELD%TYPE)
        CASE(FIELD_GEOMETRIC_TYPE,FIELD_FIBRE_TYPE)
          REGION=>FIELD%REGION
          IF(ASSOCIATED(REGION)) THEN
            COORDINATE_SYSTEM=>REGION%COORDINATE_SYSTEM
            IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
              NUMBER_OF_COMPONENTS=COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            ELSE
              LOCAL_ERROR="Coordinate system is not associated for region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Field region is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_GENERAL_TYPE,FIELD_MATERIAL_TYPE)
          NUMBER_OF_COMPONENTS=1
        CASE DEFAULT
          LOCAL_ERROR="The field type of "//TRIM(NUMBER_TO_VSTRING(FIELD%TYPE,"*",ERR,ERROR))//" is invalid for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache variable types.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%DIMENSION(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache dimension.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache dimension locked.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%DATA_TYPES(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache data type.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache data type locked.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache number of components.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache number of components locked.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache interpolation type.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache interpolation type locked.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache mesh component type.",ERR,ERROR,*999)
        ALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES), &
          & STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocated create values cache mesh component type locked.",ERR,ERROR,*999)
        FIELD%CREATE_VALUES_CACHE%DECOMPOSITION_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%DEPENDENT_TYPE_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%NUMBER_OF_VARIABLES_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%GEOMETRIC_FIELD_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%SCALING_TYPE_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%TYPE_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES=0
        FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%DIMENSION=0
        FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%DATA_TYPES=0
        FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS=0
        FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE=0
        FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED=.FALSE.
        FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER=0
        FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED=.FALSE.
        DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
          FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(variable_idx)=variable_idx
          FIELD%CREATE_VALUES_CACHE%DIMENSION(variable_idx)=FIELD_VECTOR_DIMENSION_TYPE
          FIELD%CREATE_VALUES_CACHE%DATA_TYPES(variable_idx)=FIELD_DP_TYPE
          FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_idx)=NUMBER_OF_COMPONENTS
          DO component_idx=1,NUMBER_OF_COMPONENTS
            FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(component_idx,variable_idx)=FIELD_NODE_BASED_INTERPOLATION
            FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(component_idx,variable_idx)=1
          ENDDO !component_idx
        ENDDO !variable_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL FIELD_CREATE_VALUES_CACHE_FINALISE(FIELD%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Checks the dependent type for a field.
  SUBROUTINE FIELD_DEPENDENT_TYPE_CHECK(FIELD,DEPENDENT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the dependent type for
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_TYPE !<The dependent type to check \see FIELD_ROUTINES_DependentTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DEPENDENT_TYPE_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        SELECT CASE(DEPENDENT_TYPE)
        CASE(FIELD_INDEPENDENT_TYPE)
          IF(FIELD%DEPENDENT_TYPE/=FIELD_INDEPENDENT_TYPE) THEN
            LOCAL_ERROR="Invalid dependent type. The dependent type of field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%DEPENDENT_TYPE,"*",ERR,ERROR))// &
              & " which is not an independent field."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_DEPENDENT_TYPE)
          IF(FIELD%DEPENDENT_TYPE/=FIELD_DEPENDENT_TYPE) THEN
            LOCAL_ERROR="Invalid dependent type. The dependent type of field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%DEPENDENT_TYPE,"*",ERR,ERROR))// &
              & " which is not a dependent field."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The specified dependent type of "//TRIM(NUMBER_TO_VSTRING(DEPENDENT_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DEPENDENT_TYPE_CHECK")
    RETURN
999 CALL ERRORS("FIELD_DEPENDENT_TYPE_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_DEPENDENT_TYPE_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_DEPENDENT_TYPE_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the dependent type for a field.
  SUBROUTINE FIELD_DEPENDENT_TYPE_GET(FIELD,DEPENDENT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the dependent type for
    INTEGER(INTG), INTENT(OUT) :: DEPENDENT_TYPE !<On return, the dependent type to get \see FIELD_ROUTINES_DependentTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DEPENDENT_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        DEPENDENT_TYPE=FIELD%DEPENDENT_TYPE
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DEPENDENT_TYPE_GET")
    RETURN
999 CALL ERRORS("FIELD_DEPENDENT_TYPE_GET",ERR,ERROR)
    CALL EXITS("FIELD_DEPENDENT_TYPE_GET")
    RETURN 1
  END SUBROUTINE FIELD_DEPENDENT_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field.
  SUBROUTINE FIELD_DEPENDENT_TYPE_SET(FIELD,DEPENDENT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set/change the dependent type for
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_TYPE !<The dependent type to set/change \see FIELD_ROUTINES_DependentTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DEPENDENT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%DEPENDENT_TYPE_LOCKED) THEN
            LOCAL_ERROR="The dependent type has been locked for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            SELECT CASE(DEPENDENT_TYPE)
            CASE(FIELD_INDEPENDENT_TYPE)
              FIELD%DEPENDENT_TYPE=FIELD_INDEPENDENT_TYPE
            CASE(FIELD_DEPENDENT_TYPE)
              FIELD%DEPENDENT_TYPE=FIELD_DEPENDENT_TYPE
            CASE DEFAULT
              LOCAL_ERROR="The supplied dependent type of "//TRIM(NUMBER_TO_VSTRING(DEPENDENT_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DEPENDENT_TYPE_SET")
    RETURN
999 CALL ERRORS("FIELD_DEPENDENT_TYPE_SET",ERR,ERROR)
    CALL EXITS("FIELD_DEPENDENT_TYPE_SET")
    RETURN 1
  END SUBROUTINE FIELD_DEPENDENT_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the dependent type for a field and locks so that no further changes are possible.
  SUBROUTINE FIELD_DEPENDENT_TYPE_SET_AND_LOCK(FIELD,DEPENDENT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set/change the dependent type for
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_TYPE !<The dependent type to set/change \see FIELD_ROUTINES_DependentTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DEPENDENT_TYPE_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_DEPENDENT_TYPE_SET(FIELD,DEPENDENT_TYPE,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%DEPENDENT_TYPE_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_DEPENDENT_TYPE_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_DEPENDENT_TYPE_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_DEPENDENT_TYPE_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_DEPENDENT_TYPE_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Destroys a field.
  SUBROUTINE FIELD_DESTROY(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx,field_position,field_position2
    TYPE(FIELD_TYPE), POINTER :: FIELD2,GEOMETRIC_FIELD
    TYPE(FIELD_PTR_TYPE), POINTER :: NEW_FIELDS(:),NEW_FIELDS_USING(:)
    TYPE(REGION_TYPE), POINTER :: REGION

    NULLIFY(NEW_FIELDS)
    NULLIFY(NEW_FIELDS_USING)

    CALL ENTERS("FIELD_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      REGION=>FIELD%REGION
      IF(ASSOCIATED(REGION)) THEN
        field_position=FIELD%GLOBAL_NUMBER
        GEOMETRIC_FIELD=>FIELD%GEOMETRIC_FIELD
        IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
          IF(ASSOCIATED(GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS)) THEN
            !Delete this field from the list of fields using the geometric field.
            field_position2=0
            DO field_idx=1,GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING
              FIELD2=>GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(field_idx)%PTR
              IF(FIELD2%USER_NUMBER==FIELD%USER_NUMBER) THEN
                field_position2=field_idx
                EXIT
              ENDIF
            ENDDO !field_idx
            IF(field_position2/=0) THEN
              ALLOCATE(NEW_FIELDS_USING(GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new fields using.",ERR,ERROR,*999)
              DO field_idx=1,GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING
                IF(field_idx<field_position2) THEN
                  NEW_FIELDS_USING(field_idx)%PTR=>GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(field_idx)%PTR
                ELSE IF(field_idx>field_position2) THEN
                  NEW_FIELDS_USING(field_idx-1)%PTR=>GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(field_idx)%PTR
                ENDIF
              ENDDO !field_idx
              GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING=GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS% &
                & NUMBER_OF_FIELDS_USING-1
              IF(ASSOCIATED(GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING)) &
                & DEALLOCATE(GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING)
              GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING=>NEW_FIELDS_USING
            ELSE
              !??? Error
            ENDIF
          ENDIF
        ENDIF
        CALL FIELD_FINALISE(FIELD,ERR,ERROR,*999)
        IF(REGION%FIELDS%NUMBER_OF_FIELDS>1) THEN
          ALLOCATE(NEW_FIELDS(REGION%FIELDS%NUMBER_OF_FIELDS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new fields.",ERR,ERROR,*999)
          DO field_idx=1,REGION%FIELDS%NUMBER_OF_FIELDS
            IF(field_idx<field_position) THEN
              NEW_FIELDS(field_idx)%PTR=>REGION%FIELDS%FIELDS(field_idx)%PTR
            ELSE IF(field_idx>field_position) THEN
              REGION%FIELDS%FIELDS(field_idx)%PTR%GLOBAL_NUMBER=REGION%FIELDS%FIELDS(field_idx)%PTR%GLOBAL_NUMBER-1
              NEW_FIELDS(field_idx-1)%PTR=>REGION%FIELDS%FIELDS(field_idx)%PTR
            ENDIF
          ENDDO !field_no
          DEALLOCATE(REGION%FIELDS%FIELDS)
          REGION%FIELDS%FIELDS=>NEW_FIELDS
          REGION%FIELDS%NUMBER_OF_FIELDS=REGION%FIELDS%NUMBER_OF_FIELDS-1
        ELSE
          DEALLOCATE(REGION%FIELDS%FIELDS)
          REGION%FIELDS%NUMBER_OF_FIELDS=0
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field region is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_FIELDS)) DEALLOCATE(NEW_FIELDS)
    IF(ASSOCIATED(NEW_FIELDS_USING)) DEALLOCATE(NEW_FIELDS_USING)
    CALL ERRORS("FIELD_DESTROY",ERR,ERROR)
    CALL EXITS("FIELD_DESTROY")
    RETURN 1
  END SUBROUTINE FIELD_DESTROY

  !
  !================================================================================================================================
  !

  !>Checks the field dimension for a field variable.
  SUBROUTINE FIELD_DIMENSION_CHECK(FIELD,VARIABLE_TYPE,DIMENSION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the dimension for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to check \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: DIMENSION_TYPE !<The field dimension to check \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FIELD_DIMENSION_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            SELECT CASE(DIMENSION_TYPE)
            CASE(FIELD_SCALAR_DIMENSION_TYPE)
              IF(FIELD_VARIABLE%DIMENSION/=FIELD_SCALAR_DIMENSION_TYPE) THEN
                LOCAL_ERROR="Invalid dimension type. The dimension type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DIMENSION,"*",ERR,ERROR))// &
                  & " which is not a scalar field."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE(FIELD_VECTOR_DIMENSION_TYPE)
              IF(FIELD_VARIABLE%DIMENSION/=FIELD_VECTOR_DIMENSION_TYPE) THEN
                LOCAL_ERROR="Invalid dimension type. The dimension type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DIMENSION,"*",ERR,ERROR))// &
                  & " which is not a vector field."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            CASE(FIELD_TENSOR_DIMENSION_TYPE) 
              IF(FIELD_VARIABLE%DIMENSION/=FIELD_TENSOR_DIMENSION_TYPE) THEN
                LOCAL_ERROR="Invalid dimension type. The dimension type for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DIMENSION,"*",ERR,ERROR))// &
                  & " which is not a tensor field."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
           CASE DEFAULT
              LOCAL_ERROR="The specified dimension type of "//TRIM(NUMBER_TO_VSTRING(DIMENSION_TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DIMENSION_CHECK")
    RETURN
999 CALL ERRORS("FIELD_DIMENSION_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_DIMENSION_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_DIMENSION_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the field dimension for a field variable.
  SUBROUTINE FIELD_DIMENSION_GET(FIELD,VARIABLE_TYPE,DIMENSION,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the dimension for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(OUT) :: DIMENSION !<On return, the field dimension to get \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FIELD_DIMENSION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            DIMENSION=FIELD_VARIABLE%DIMENSION
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DIMENSION_GET")
    RETURN
999 CALL ERRORS("FIELD_DIMENSION_GET",ERR,ERROR)
    CALL EXITS("FIELD_DIMENSION_GET")
    RETURN 1
  END SUBROUTINE FIELD_DIMENSION_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field dimension for a field variable.
  SUBROUTINE FIELD_DIMENSION_SET(FIELD,VARIABLE_TYPE,FIELD_DIMENSION,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set/change the dimension for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_DIMENSION !<The field dimension to set/change \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS
    INTEGER(INTG), ALLOCATABLE :: OLD_INTERPOLATION_TYPE(:,:),OLD_MESH_COMPONENT_NUMBER(:,:)
    LOGICAL, ALLOCATABLE :: OLD_INTERPOLATION_TYPE_LOCKED(:,:),OLD_MESH_COMPONENT_NUMBER_LOCKED(:,:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DIMENSION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            IF(ANY(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES==VARIABLE_TYPE)) THEN
              IF(FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED(VARIABLE_TYPE)) THEN
                LOCAL_ERROR="The field dimension has been locked for for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                SELECT CASE(FIELD_DIMENSION)
                CASE(FIELD_SCALAR_DIMENSION_TYPE)
                  IF(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)/=1) THEN
                    NUMBER_OF_COMPONENTS=SIZE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE,1)
                    ALLOCATE(OLD_INTERPOLATION_TYPE(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type.",ERR,ERROR,*999)
                    ALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type locked.",ERR,ERROR,*999)
                    ALLOCATE(OLD_MESH_COMPONENT_NUMBER(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old mesh component number.",ERR,ERROR,*999)
                    ALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old mesh component number locked.",ERR,ERROR,*999)
                    OLD_INTERPOLATION_TYPE=FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE
                    OLD_INTERPOLATION_TYPE_LOCKED=FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED
                    OLD_MESH_COMPONENT_NUMBER=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER
                    OLD_MESH_COMPONENT_NUMBER_LOCKED=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED
                    DEALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE)
                    DEALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED)
                    DEALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER)
                    DEALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED)
                    ALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(1,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type.",ERR,ERROR,*999)
                    ALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(1,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type locked.",ERR,ERROR,*999)
                    ALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(1,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh component number.",ERR,ERROR,*999)
                    ALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(1,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh component number locked.",ERR,ERROR,*999)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(1,:)=OLD_INTERPOLATION_TYPE(1,:)
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(1,:)=OLD_INTERPOLATION_TYPE_LOCKED(1,:)
                    FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(1,:)=OLD_MESH_COMPONENT_NUMBER(1,:)
                    FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(1,:)=OLD_MESH_COMPONENT_NUMBER_LOCKED(1,:)
                    DEALLOCATE(OLD_INTERPOLATION_TYPE)
                    DEALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED)
                    DEALLOCATE(OLD_MESH_COMPONENT_NUMBER)
                    DEALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED)
                    FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)=1
                  ENDIF
                  FIELD%CREATE_VALUES_CACHE%DIMENSION(VARIABLE_TYPE)=FIELD_SCALAR_DIMENSION_TYPE
                CASE(FIELD_VECTOR_DIMENSION_TYPE)
                  FIELD%CREATE_VALUES_CACHE%DIMENSION(VARIABLE_TYPE)=FIELD_VECTOR_DIMENSION_TYPE
                CASE(FIELD_TENSOR_DIMENSION_TYPE)
                  FIELD%CREATE_VALUES_CACHE%DIMENSION(VARIABLE_TYPE)=FIELD_TENSOR_DIMENSION_TYPE
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The specified field dimension of "//TRIM(NUMBER_TO_VSTRING(FIELD_DIMENSION,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_DIMENSION_SET")
    RETURN
999 IF(ALLOCATED(OLD_INTERPOLATION_TYPE)) DEALLOCATE(OLD_INTERPOLATION_TYPE)
    IF(ALLOCATED(OLD_INTERPOLATION_TYPE_LOCKED)) DEALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED)
    IF(ALLOCATED(OLD_MESH_COMPONENT_NUMBER)) DEALLOCATE(OLD_MESH_COMPONENT_NUMBER)
    IF(ALLOCATED(OLD_MESH_COMPONENT_NUMBER_LOCKED)) DEALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED)
    CALL ERRORS("FIELD_DIMENSION_SET",ERR,ERROR)
    CALL EXITS("FIELD_DIMENSION_SET")
    RETURN 1
  END SUBROUTINE FIELD_DIMENSION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field dimension for a field variable  and locks so that no further changes can be made.
  SUBROUTINE FIELD_DIMENSION_SET_AND_LOCK(FIELD,VARIABLE_TYPE,FIELD_DIMENSION,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set/change the dimension for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_DIMENSION !<The field dimension to set/change \see FIELD_ROUTINES_DimensionTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_DIMENSION_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_DIMENSION_SET(FIELD,VARIABLE_TYPE,FIELD_DIMENSION,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED(VARIABLE_TYPE)=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_DIMENSION_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_DIMENSION_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_DIMENSION_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_DIMENSION_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Finalises a field and deallocates all memory.
  SUBROUTINE FIELD_FINALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      CALL FIELD_SCALINGS_FINALISE(FIELD,ERR,ERROR,*999)
      CALL FIELD_VARIABLES_FINALISE(FIELD,ERR,ERROR,*999)
      CALL FIELD_CREATE_VALUES_CACHE_FINALISE(FIELD%CREATE_VALUES_CACHE,ERR,ERROR,*999)
      CALL FIELD_GEOMETRIC_PARAMETERS_FINALISE(FIELD%GEOMETRIC_FIELD_PARAMETERS,ERR,ERROR,*999)
      IF(ALLOCATED(FIELD%VARIABLE_TYPE_MAP)) DEALLOCATE(FIELD%VARIABLE_TYPE_MAP)
      DEALLOCATE(FIELD)     
    ENDIF

    CALL EXITS("FIELD_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a field
  SUBROUTINE FIELD_INITIALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,variable_type_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("FIELD_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(FIELD)) THEN
      FIELD%GLOBAL_NUMBER=0
      FIELD%USER_NUMBER=0
      FIELD%FIELD_FINISHED=.FALSE.
      NULLIFY(FIELD%FIELDS)
      NULLIFY(FIELD%REGION)
      FIELD%TYPE=FIELD_GEOMETRIC_TYPE
      FIELD%DEPENDENT_TYPE=FIELD_INDEPENDENT_TYPE
      NULLIFY(FIELD%DECOMPOSITION)
      FIELD%NUMBER_OF_VARIABLES=0
      NULLIFY(FIELD%GEOMETRIC_FIELD)
      NULLIFY(FIELD%GEOMETRIC_FIELD_PARAMETERS)
      NULLIFY(FIELD%CREATE_VALUES_CACHE)
      ALLOCATE(FIELD%VARIABLE_TYPE_MAP(FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable type map.",ERR,ERROR,*999)
      DO variable_type_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
        NULLIFY(FIELD%VARIABLE_TYPE_MAP(variable_type_idx)%PTR)
      ENDDO !variable_type_idx
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_INITIALISE")
    RETURN
999 CALL FIELD_FINALISE(FIELD,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_INITIALISE

  !
  !================================================================================================================================
  !

  !>Interpolates a field at a gauss point to give an interpolated point. PARTIAL_DERIVATIVE_TYPE controls which partial derivatives are evaluated. If it is NO_PART_DERIV then only the field values are interpolated. If it is FIRST_PART_DERIV then the field values and first partial derivatives are interpolated. If it is SECOND_PART_DERIV the the field values and first and second partial derivatives are evaluated. Old CMISS name XEXG, ZEXG
  SUBROUTINE FIELD_INTERPOLATE_GAUSS(PARTIAL_DERIVATIVE_TYPE,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER,INTERPOLATED_POINT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_TYPE !<The partial derivative type of the provided field interpolation
    INTEGER(INTG), INTENT(IN) :: QUADRATURE_SCHEME !<The quadrature scheme of the Gauss points \see BASIS_ROUTINES_QuadratureSchemes,BASIS_ROUTINES
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The number of the Gauss point to interpolate the field at
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<The pointer to the interpolated point which will contain the field interpolation information at the specified Gauss point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ni,nu
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATE_GAUSS",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
      INTERPOLATION_PARAMETERS=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS
      IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
        COORDINATE_SYSTEM=>INTERPOLATION_PARAMETERS%FIELD%REGION%COORDINATE_SYSTEM
        SELECT CASE(PARTIAL_DERIVATIVE_TYPE)
        CASE(NO_PART_DERIV)
          DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)              
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              INTERPOLATED_POINT%VALUES(component_idx,1)=BASIS_INTERPOLATE_GAUSS(INTERPOLATION_PARAMETERS%BASES( &
                & component_idx)%PTR,NO_PART_DERIV,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER,INTERPOLATION_PARAMETERS% &
                & PARAMETERS(:,component_idx),ERR,ERROR)
              IF(ERR/=0) GOTO 999
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS% &
                & FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                & " is invalid for component index "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//"."
            END SELECT
            CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
              & ERR,ERROR,*999)
          ENDDO! component_idx
          INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE=NO_PART_DERIV
        CASE(FIRST_PART_DERIV)
          DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              !Handle the first case of no partial derivative
              INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
              CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                & ERR,ERROR,*999)
              !Now process all the first partial derivatives
              DO ni=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_XI
                nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
                INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                  & ERR,ERROR,*999)
              ENDDO !ni
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              !Handle the first case of no partial derivative
              INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
              CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                & ERR,ERROR,*999)
              !Now process all the first partial derivatives
              DO ni=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_XI
                nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
                INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                  & ERR,ERROR,*999)
              ENDDO !ni
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              !Handle the first case of no partial derivative
              INTERPOLATED_POINT%VALUES(component_idx,1)=BASIS_INTERPOLATE_GAUSS(INTERPOLATION_PARAMETERS%BASES( &
                & component_idx)%PTR,NO_PART_DERIV,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER,INTERPOLATION_PARAMETERS% &
                & PARAMETERS(:,component_idx),ERR,ERROR)
              IF(ERR/=0) GOTO 999
              CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                & ERR,ERROR,*999)
              !Now process all the first partial derivatives
              DO ni=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_XI
                nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
                INTERPOLATED_POINT%VALUES(component_idx,nu)=BASIS_INTERPOLATE_GAUSS(INTERPOLATION_PARAMETERS% &
                  & BASES(component_idx)%PTR,nu,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER, &
                  & INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx),ERR,ERROR)
                IF(ERR/=0) GOTO 999
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                  & ERR,ERROR,*999)
              ENDDO !ni
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS% &
                & FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                & " is invalid for component index "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//"."
            END SELECT
          ENDDO! component_idx
          INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE=FIRST_PART_DERIV
        CASE(SECOND_PART_DERIV)
          DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)
              !Handle the first case of no partial derivative
              INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
              CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                & ERR,ERROR,*999)
              !Now process the rest of partial derivatives
              DO nu=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_PARTIAL_DERIVATIVES
                INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                  & ERR,ERROR,*999)
              ENDDO !nu
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              !Handle the first case of no partial derivative
              INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
              CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                & ERR,ERROR,*999)
              !Now process the rest of partial derivatives
              DO nu=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_PARTIAL_DERIVATIVES
                INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                  & ERR,ERROR,*999)
              ENDDO !nu
            CASE(FIELD_NODE_BASED_INTERPOLATION)              
              DO nu=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_PARTIAL_DERIVATIVES
                INTERPOLATED_POINT%VALUES(component_idx,nu)=BASIS_INTERPOLATE_GAUSS(INTERPOLATION_PARAMETERS% &
                  & BASES(component_idx)%PTR,nu,QUADRATURE_SCHEME,GAUSS_POINT_NUMBER, &
                  & INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx),ERR,ERROR)
                IF(ERR/=0) GOTO 999
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                  & ERR,ERROR,*999)
              ENDDO! nu
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS% &
                & FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                & " is invalid for component index "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//"."
            END SELECT
          ENDDO !component_idx
          INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE=SECOND_PART_DERIV
        CASE DEFAULT
          LOCAL_ERROR="The partial derivative type of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIVATIVE_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Interpolated point interpolation parameters is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolated point is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATE_GAUSS")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATE_GAUSS",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATE_GAUSS")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATE_GAUSS

  !
  !================================================================================================================================
  !

  !>Interpolates a field at a xi location to give an interpolated point. XI is the element location to be interpolated at. PARTIAL_DERIVATIVE_TYPE controls which partial derivatives are evaluated. If it is NO_PART_DERIV then only the field values are interpolated. If it is FIRST_PART_DERIV then the field values and first partial derivatives are interpolated. If it is SECOND_PART_DERIV the the field values and first and second partial derivatives are evaluated. Old CMISS name PXI
  SUBROUTINE FIELD_INTERPOLATE_XI(PARTIAL_DERIVATIVE_TYPE,XI,INTERPOLATED_POINT,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARTIAL_DERIVATIVE_TYPE !<The partial derivative type of the provide field interpolation
    REAL(DP), INTENT(IN) :: XI(:) !<XI(ni). The ni'th Xi coordinate to evaluate the field at
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<The pointer to the interpolated point which will contain the field interpolation information at the specified Xi point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ni,nu
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATE_XI",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
      INTERPOLATION_PARAMETERS=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS
      IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
        !!TODO: Fix this check. You can have less Xi directions than the mesh number of dimensions e.g., interpolating a line
        !IF(SIZE(XI,1)>=INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS) THEN
          COORDINATE_SYSTEM=>INTERPOLATION_PARAMETERS%FIELD%REGION%COORDINATE_SYSTEM
          SELECT CASE(PARTIAL_DERIVATIVE_TYPE)
          CASE(NO_PART_DERIV)
            DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)              
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                INTERPOLATED_POINT%VALUES(component_idx,1)=BASIS_INTERPOLATE_XI(INTERPOLATION_PARAMETERS% &
                  & BASES(component_idx)%PTR,NO_PART_DERIV,XI,INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx),ERR,ERROR)
                IF(ERR/=0) GOTO 999
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for component index "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//"."
              END SELECT
              CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                & ERR,ERROR,*999)
            ENDDO !component_idx
            INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE=NO_PART_DERIV
         CASE(FIRST_PART_DERIV)
            DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                !Handle the first case of no partial derivative
                INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                  & ERR,ERROR,*999)
                !Now process all the first partial derivatives
                DO ni=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_XI
                  nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
                  INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                  CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                    & ERR,ERROR,*999)
                ENDDO !ni
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                !Handle the first case of no partial derivative
                INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                  & ERR,ERROR,*999)
                !Now process all the first partial derivatives
                DO ni=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_XI
                  nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
                  INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                  CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                    & ERR,ERROR,*999)
                ENDDO !ni
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                !Handle the first case of no partial derivative
                INTERPOLATED_POINT%VALUES(component_idx,1)=BASIS_INTERPOLATE_XI(INTERPOLATION_PARAMETERS% &
                  & BASES(component_idx)%PTR,NO_PART_DERIV,XI,INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx),ERR,ERROR)
                IF(ERR/=0) GOTO 999
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                  & ERR,ERROR,*999)
                !Now process all the first partial derivatives
                DO ni=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_XI
                  nu=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni)
                  INTERPOLATED_POINT%VALUES(component_idx,nu)=BASIS_INTERPOLATE_XI(INTERPOLATION_PARAMETERS% &
                    & BASES(component_idx)%PTR,nu,XI,INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx), &
                    & ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                    & ERR,ERROR,*999)
                ENDDO !ni
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for component index "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//"."
              END SELECT
            ENDDO !component_idx
            INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE=FIRST_PART_DERIV
          CASE(SECOND_PART_DERIV)
            DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                !Handle the first case of no partial derivative
                INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                  & ERR,ERROR,*999)
                !Now process the rest of partial derivatives
                DO nu=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_PARTIAL_DERIVATIVES
                  INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                  CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                    & ERR,ERROR,*999)
                ENDDO !nu
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                !Handle the first case of no partial derivative
                INTERPOLATED_POINT%VALUES(component_idx,1)=INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)
                CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,NO_PART_DERIV,INTERPOLATED_POINT%VALUES(component_idx,1), &
                  & ERR,ERROR,*999)
                !Now process the rest of partial derivatives
                DO nu=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_PARTIAL_DERIVATIVES
                  INTERPOLATED_POINT%VALUES(component_idx,nu)=0.0_DP
                  CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                    & ERR,ERROR,*999)
                ENDDO !nu
              CASE(FIELD_NODE_BASED_INTERPOLATION)              
                DO nu=1,INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR%NUMBER_OF_PARTIAL_DERIVATIVES
                  INTERPOLATED_POINT%VALUES(component_idx,nu)=BASIS_INTERPOLATE_XI(INTERPOLATION_PARAMETERS% &
                    & BASES(component_idx)%PTR,nu,XI,INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx), &
                    & ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  CALL COORDINATE_INTERPOLATION_ADJUST(COORDINATE_SYSTEM,nu,INTERPOLATED_POINT%VALUES(component_idx,nu), &
                    & ERR,ERROR,*999)
                ENDDO! nu
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for component index "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//"."
              END SELECT
            ENDDO !component_idx
            INTERPOLATED_POINT%PARTIAL_DERIVATIVE_TYPE=SECOND_PART_DERIV
          CASE DEFAULT
            LOCAL_ERROR="The partial derivative type of "//TRIM(NUMBER_TO_VSTRING(PARTIAL_DERIVATIVE_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !ELSE
        !  LOCAL_ERROR="Invalid number of Xi directions. The supplied Xi has "// &
        !    & TRIM(NUMBER_TO_VSTRING(SIZE(XI,1),"*",ERR,ERROR))//" directions and the required number of directions is "// &
        !    & TRIM(NUMBER_TO_VSTRING(INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS, &
        !    & "*",ERR,ERROR))
        !  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        !ENDIF
      ELSE
        CALL FLAG_ERROR("Interpolated point interpolation parameters is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolated point is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATE_XI")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATE_XI",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATE_XI")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATE_XI

  !
  !================================================================================================================================
  !

  !>Finalises the interpolated point and deallocates all memory.
  SUBROUTINE FIELD_INTERPOLATED_POINT_FINALISE(INTERPOLATED_POINT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<A pointer to the interpolated point to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_INTERPOLATED_POINT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
      IF(ALLOCATED(INTERPOLATED_POINT%VALUES)) DEALLOCATE(INTERPOLATED_POINT%VALUES)
      DEALLOCATE(INTERPOLATED_POINT)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATED_POINT_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATED_POINT_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATED_POINT_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATED_POINT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interpolated point for an interpolation parameters
  SUBROUTINE FIELD_INTERPOLATED_POINT_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<A pointer to the interpolation parameters to initialise the interpolated point for
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<On exit, A pointer to the interpolated point that has been initialised
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,NUMBER_OF_DIMENSIONS
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("FIELD_INTERPOLATED_POINT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
      IF(ASSOCIATED(INTERPOLATION_PARAMETERS%FIELD)) THEN
        IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
          CALL FLAG_ERROR("Interpolated point is already associated.",ERR,ERROR,*998)
        ELSE
          ALLOCATE(INTERPOLATED_POINT,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point",ERR,ERROR,*999)
          INTERPOLATED_POINT%INTERPOLATION_PARAMETERS=>INTERPOLATION_PARAMETERS
          NUMBER_OF_DIMENSIONS=INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS
          INTERPOLATED_POINT%MAX_PARTIAL_DERIVATIVE_INDEX=PARTIAL_DERIVATIVE_MAXIMUM_MAP(NUMBER_OF_DIMENSIONS)
          ALLOCATE(INTERPOLATED_POINT%VALUES(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS, &
            & INTERPOLATED_POINT%MAX_PARTIAL_DERIVATIVE_INDEX),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point values.",ERR,ERROR,*999)
          INTERPOLATED_POINT%VALUES=0.0_DP
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interpolation parameters field is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolation parameters is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATED_POINT_INITIALISE")
    RETURN
999 CALL FIELD_INTERPOLATED_POINT_FINALISE(INTERPOLATED_POINT,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_INTERPOLATED_POINT_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATED_POINT_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATED_POINT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the interpolated point metrics and the associated interpolated point
  SUBROUTINE FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(JACOBIAN_TYPE,INTERPOLATED_POINT_METRICS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: INTERPOLATED_POINT_METRICS !<A pointer to the interpolated point metrics
    INTEGER(INTG), INTENT(IN) :: JACOBIAN_TYPE !<The Jacobian type of the calculation \see COORDINATE_ROUTINES_JacobianTypes,COORDINATE_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS

    CALL ENTERS("FIELD_INTERPOLATED_POINT_METRICS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATED_POINT_METRICS)) THEN
      INTERPOLATED_POINT=>INTERPOLATED_POINT_METRICS%INTERPOLATED_POINT
      INTERPOLATION_PARAMETERS=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS
      FIELD=>INTERPOLATION_PARAMETERS%FIELD
      COORDINATE_SYSTEM=>FIELD%REGION%COORDINATE_SYSTEM
      IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
        CALL COORDINATE_METRICS_CALCULATE(COORDINATE_SYSTEM,JACOBIAN_TYPE,INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("The field is not a geometric or fibre field.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolated point metrics is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATED_POINT_METRICS_CALCULATE")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATED_POINT_METRICS_CALCULATE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATED_POINT_METRICS_CALCULATE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATED_POINT_METRICS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the interpolated point metrics and deallocates all memory.
  SUBROUTINE FIELD_INTERPOLATED_POINT_METRICS_FINALISE(INTERPOLATED_POINT_METRICS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: INTERPOLATED_POINT_METRICS !<A pointer to the interpolated point metrics to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_INTERPOLATED_POINT_METRICS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATED_POINT_METRICS)) THEN
      IF(ALLOCATED(INTERPOLATED_POINT_METRICS%GL)) DEALLOCATE(INTERPOLATED_POINT_METRICS%GL)
      IF(ALLOCATED(INTERPOLATED_POINT_METRICS%GU)) DEALLOCATE(INTERPOLATED_POINT_METRICS%GU)
      IF(ALLOCATED(INTERPOLATED_POINT_METRICS%DX_DXI)) DEALLOCATE(INTERPOLATED_POINT_METRICS%DX_DXI)
      IF(ALLOCATED(INTERPOLATED_POINT_METRICS%DXI_DX)) DEALLOCATE(INTERPOLATED_POINT_METRICS%DXI_DX)
      DEALLOCATE(INTERPOLATED_POINT_METRICS)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATED_POINT_METRICS_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATED_POINT_METRICS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATED_POINT_METRICS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATED_POINT_METRICS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interpolated point metrics for an interpolated point.
  SUBROUTINE FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(INTERPOLATED_POINT,INTERPOLATED_POINT_METRICS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !A pointer to the interpolated point to initliase the interpolated point metrics for
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: INTERPOLATED_POINT_METRICS !<On exit, a pointer to the interpolated point metrics that have been initialised
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_XI_DIMENSIONS,NUMBER_OF_X_DIMENSIONS
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATED_POINT_METRICS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATED_POINT)) THEN
      IF(ASSOCIATED(INTERPOLATED_POINT_METRICS)) THEN
        CALL FLAG_ERROR("Interpolated point metrics is already associated.",ERR,ERROR,*998)
      ELSE
        NUMBER_OF_X_DIMENSIONS=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        NUMBER_OF_XI_DIMENSIONS=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS
        IF(NUMBER_OF_X_DIMENSIONS==SIZE(INTERPOLATED_POINT%VALUES,1)) THEN
          ALLOCATE(INTERPOLATED_POINT_METRICS,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point metrics.",ERR,ERROR,*999)
          ALLOCATE(INTERPOLATED_POINT_METRICS%GL(NUMBER_OF_XI_DIMENSIONS,NUMBER_OF_XI_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point metrics convariant tensor.",ERR,ERROR,*999)
          ALLOCATE(INTERPOLATED_POINT_METRICS%GU(NUMBER_OF_XI_DIMENSIONS,NUMBER_OF_XI_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point metrics contravariant tensor.",ERR,ERROR,*999)
          ALLOCATE(INTERPOLATED_POINT_METRICS%DX_DXI(NUMBER_OF_X_DIMENSIONS,NUMBER_OF_XI_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point metrics dX_dXi.",ERR,ERROR,*999)
          ALLOCATE(INTERPOLATED_POINT_METRICS%DXI_DX(NUMBER_OF_XI_DIMENSIONS,NUMBER_OF_X_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolated point metrics dXi_dX.",ERR,ERROR,*999)
          INTERPOLATED_POINT_METRICS%INTERPOLATED_POINT=>INTERPOLATED_POINT
          INTERPOLATED_POINT_METRICS%NUMBER_OF_X_DIMENSIONS=NUMBER_OF_X_DIMENSIONS
          INTERPOLATED_POINT_METRICS%NUMBER_OF_XI_DIMENSIONS=NUMBER_OF_XI_DIMENSIONS
          INTERPOLATED_POINT_METRICS%GL=0.0_DP
          INTERPOLATED_POINT_METRICS%GU=0.0_DP
          INTERPOLATED_POINT_METRICS%DX_DXI=0.0_DP
          INTERPOLATED_POINT_METRICS%DXI_DX=0.0_DP
          INTERPOLATED_POINT_METRICS%JACOBIAN=0.0_DP
          INTERPOLATED_POINT_METRICS%JACOBIAN_TYPE=0
        ELSE
          LOCAL_ERROR="The number of coordinate dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_X_DIMENSIONS,"*",ERR,ERROR))// &
            & ") does not match the number of components of the interpolated point ("// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(INTERPOLATED_POINT%VALUES,1),"*",ERR,ERROR))//")."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolation point is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATED_POINT_METRICS_INITIALISE")
    RETURN
999 CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(INTERPOLATED_POINT_METRICS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_INTERPOLATED_POINT_METRICS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATED_POINT_METRICS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATED_POINT_METRICS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the interpolation parameters for a particular element. Old CMISS name XPXE, ZPZE
  SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(PARAMETER_SET_NUMBER,ELEMENT_NUMBER,INTERPOLATION_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARAMETER_SET_NUMBER !<The field parameter set number to get the element parameters for
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to get the element parameters for
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<On return, a pointer to the interpolation parameters
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,mk,nk,nn,np,ns,ny,ny2,scaling_idx
    REAL(DP), POINTER :: FIELD_PARAMETER_SET_DATA(:),SCALE_FACTORS(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES_TOPOLOGY
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
      IF(PARAMETER_SET_NUMBER>0.AND.PARAMETER_SET_NUMBER<=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%PARAMETER_SETS% &
        & NUMBER_OF_PARAMETER_SETS) THEN
        PARAMETER_SET=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(PARAMETER_SET_NUMBER)%PTR
        IF(ASSOCIATED(PARAMETER_SET)) THEN
          NULLIFY(FIELD_PARAMETER_SET_DATA)
          CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETER_SET%PARAMETERS,FIELD_PARAMETER_SET_DATA,ERR,ERROR,*999)
          COORDINATE_SYSTEM=>INTERPOLATION_PARAMETERS%FIELD%REGION%COORDINATE_SYSTEM
          DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            ELEMENTS_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
            IF(ELEMENT_NUMBER>0.AND.ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
              BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
              INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR=>BASIS
              SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                ny=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                  & CONSTANT_PARAM2DOF_MAP
                INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
                INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=FIELD_PARAMETER_SET_DATA(ny)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                ny=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                  & ELEMENT_PARAM2DOF_MAP(ELEMENT_NUMBER)
                INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
                INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=FIELD_PARAMETER_SET_DATA(ny)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                ELEMENTS_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
                NODES_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES
                INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                SELECT CASE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALING_TYPE)
                CASE(FIELD_NO_SCALING)
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    np=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nn)
                    DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                      nk=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(mk,nn)
                      ns=BASIS%ELEMENT_PARAMETER_INDEX(mk,nn)
                      ny=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                        & NODE_PARAM2DOF_MAP(nk,np)
                      INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)
                    ENDDO !mk
                  ENDDO !nn
                CASE(FIELD_UNIT_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
                  scaling_idx=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%SCALING_INDEX
                  NULLIFY(SCALE_FACTORS)
                  CALL DISTRIBUTED_VECTOR_DATA_GET(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                    & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    np=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nn)
                    DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                      nk=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(mk,nn)
                      ns=BASIS%ELEMENT_PARAMETER_INDEX(nk,nn)
                      ny=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                        & NODE_PARAM2DOF_MAP(nk,np)
                      ny2=NODES_TOPOLOGY%NODES(np)%DOF_INDEX(nk)
                      !INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)* &
                      !  & INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)%SCALE_FACTORS(ns,ELEMENT_NUMBER)
                      !INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)* &
                      !  & INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)%SCALE_FACTORS(nk,np)
                      INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)*SCALE_FACTORS(ny2)
                      INTERPOLATION_PARAMETERS%SCALE_FACTORS(ns,component_idx)=SCALE_FACTORS(ny2)
                    ENDDO !mk
                  ENDDO !nn
                  CALL DISTRIBUTED_VECTOR_DATA_RESTORE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                    & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
                CASE(FIELD_ARC_LENGTH_SCALING)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The scaling type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%SCALINGS% &
                    & SCALING_TYPE,"*",ERR,ERROR))//" is invalid for field number "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
                  & COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              LOCAL_ERROR="The element number of "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The number must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
                & " for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
                & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !component_idx
          CALL COORDINATE_INTERPOLATION_PARAMETERS_ADJUST(COORDINATE_SYSTEM,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="The field parameter set number of "//TRIM(NUMBER_TO_VSTRING(PARAMETER_SET_NUMBER,"*",ERR,ERROR))// &
            & " has not been created for field number "// &
            & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The field parameter set number of "//TRIM(NUMBER_TO_VSTRING(PARAMETER_SET_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The number must be between 1 and "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
          & PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS,"*",ERR,ERROR))//" for field number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolation parameters is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interpolation parameters:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field number = ",INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field variable number = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & VARIABLE_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parameter set number = ",PARAMETER_SET_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of components = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Component = ",component_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of parameters = ",INTERPOLATION_PARAMETERS% &
          & NUMBER_OF_PARAMETERS(component_idx),ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx),4,4, &
          & INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx),'("      Parameters :",4(X,E13.6))','(18X,4(X,E13.6))', &
          & ERR,ERROR,*999)
      ENDDO !component_idx
    ENDIF

    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET

  !
  !================================================================================================================================
  !

  !>Finalises the interpolation parameters and deallocates all memory
  SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_FINALISE(INTERPOLATION_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<A pointer to the interpolation parameters to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_INTERPOLATION_PARAMETERS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
      IF(ALLOCATED(INTERPOLATION_PARAMETERS%BASES)) DEALLOCATE(INTERPOLATION_PARAMETERS%BASES)
      IF(ALLOCATED(INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS)) DEALLOCATE(INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS)
      IF(ALLOCATED(INTERPOLATION_PARAMETERS%PARAMETERS)) DEALLOCATE(INTERPOLATION_PARAMETERS%PARAMETERS)
      IF(ALLOCATED(INTERPOLATION_PARAMETERS%SCALE_FACTORS)) DEALLOCATE(INTERPOLATION_PARAMETERS%SCALE_FACTORS)      
      DEALLOCATE(INTERPOLATION_PARAMETERS)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATION_PARAMETERS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interpolation parameters for a field variable.
  SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,VARIABLE_TYPE,INTERPOLATION_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the interpolation parameters for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to initialise the interpolation parameters for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<On exit, a pointer to the initialised interpolation parameters.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,DUMMY_ERR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATION_PARAMETERS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
              CALL FLAG_ERROR("Interpolation parameters is already associated.",ERR,ERROR,*998)
            ELSE
              NULLIFY(INTERPOLATION_PARAMETERS)
              ALLOCATE(INTERPOLATION_PARAMETERS,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation parameters.",ERR,ERROR,*999)
              INTERPOLATION_PARAMETERS%FIELD=>FIELD
              INTERPOLATION_PARAMETERS%FIELD_VARIABLE=>FIELD_VARIABLE
              ALLOCATE(INTERPOLATION_PARAMETERS%BASES(FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bases.",ERR,ERROR,*999)
              ALLOCATE(INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type.",ERR,ERROR,*999)
              ALLOCATE(INTERPOLATION_PARAMETERS%PARAMETERS(FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS, &
                & FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate parameters.",ERR,ERROR,*999)
              INTERPOLATION_PARAMETERS%PARAMETERS=0.0_DP
              IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                ALLOCATE(INTERPOLATION_PARAMETERS%SCALE_FACTORS(FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS, &
                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate scale factors.",ERR,ERROR,*999)
                INTERPOLATION_PARAMETERS%SCALE_FACTORS=0.0_DP
              ENDIF
              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                NULLIFY(INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR)
              ENDDO !component_idx
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS=0
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_INITIALISE")
    RETURN
999 CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(INTERPOLATION_PARAMETERS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_INTERPOLATION_PARAMETERS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the interpolation parameters for a particular line. Old CMISS name XPXE, ZPZE
  SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_LINE_GET(PARAMETER_SET_NUMBER,LINE_NUMBER,INTERPOLATION_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: PARAMETER_SET_NUMBER !<The field parameter set number to get the line parameters for
    INTEGER(INTG), INTENT(IN) :: LINE_NUMBER !<The line number to get the line parameters for
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<On return, a pointer to the interpolation parameters
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,mk,nk,nn,np,ns,ny,ny2,scaling_idx
    REAL(DP), POINTER :: FIELD_PARAMETER_SET_DATA(:),SCALE_FACTORS(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(DOMAIN_LINES_TYPE), POINTER :: LINES_TOPOLOGY
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES_TOPOLOGY
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATION_PARAMETERS_LINE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
      IF(PARAMETER_SET_NUMBER>0.AND.PARAMETER_SET_NUMBER<=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%PARAMETER_SETS% &
        & NUMBER_OF_PARAMETER_SETS) THEN
        PARAMETER_SET=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(PARAMETER_SET_NUMBER)%PTR
        IF(ASSOCIATED(PARAMETER_SET)) THEN
          NULLIFY(FIELD_PARAMETER_SET_DATA)
          CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETER_SET%PARAMETERS,FIELD_PARAMETER_SET_DATA,ERR,ERROR,*999)
          COORDINATE_SYSTEM=>INTERPOLATION_PARAMETERS%FIELD%REGION%COORDINATE_SYSTEM
          DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            LINES_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%LINES
            IF(LINE_NUMBER>0.AND.LINE_NUMBER<=LINES_TOPOLOGY%NUMBER_OF_LINES) THEN
              BASIS=>LINES_TOPOLOGY%LINES(LINE_NUMBER)%BASIS
              INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR=>BASIS
              SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)
                NODES_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES
                INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                SELECT CASE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALING_TYPE)
                CASE(FIELD_NO_SCALING)
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    np=LINES_TOPOLOGY%LINES(LINE_NUMBER)%NODES_IN_LINE(nn)
                    DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                      nk=LINES_TOPOLOGY%LINES(LINE_NUMBER)%DERIVATIVES_IN_LINE(mk,nn)
                      ns=BASIS%ELEMENT_PARAMETER_INDEX(mk,nn)
                      ny=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                        & NODE_PARAM2DOF_MAP(nk,np)
                      INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)
                    ENDDO !mk
                  ENDDO !nn
                CASE(FIELD_UNIT_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
                  scaling_idx=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%SCALING_INDEX
                  NULLIFY(SCALE_FACTORS)
                  CALL DISTRIBUTED_VECTOR_DATA_GET(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                    & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    np=LINES_TOPOLOGY%LINES(LINE_NUMBER)%NODES_IN_LINE(nn)
                    DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                      nk=LINES_TOPOLOGY%LINES(LINE_NUMBER)%DERIVATIVES_IN_LINE(mk,nn)
                      ns=BASIS%ELEMENT_PARAMETER_INDEX(mk,nn)
                      ny=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                        & NODE_PARAM2DOF_MAP(nk,np)
                      ny2=NODES_TOPOLOGY%NODES(np)%DOF_INDEX(nk)
                      !INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)* &
                      !  & INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)%SCALE_FACTORS(nk,np)
                      INTERPOLATION_PARAMETERS%PARAMETERS(ns,component_idx)=FIELD_PARAMETER_SET_DATA(ny)*SCALE_FACTORS(ny2)
                      INTERPOLATION_PARAMETERS%SCALE_FACTORS(ns,component_idx)=SCALE_FACTORS(ny2)
                    ENDDO !mk
                  ENDDO !nn
                CASE(FIELD_ARC_LENGTH_SCALING)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The scaling type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%SCALINGS% &
                    & SCALING_TYPE,"*",ERR,ERROR))//" is invalid for field number "// &
                    & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
                  & COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                  & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              LOCAL_ERROR="The line number of "//TRIM(NUMBER_TO_VSTRING(LINE_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The number must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(LINES_TOPOLOGY%NUMBER_OF_LINES,"*",ERR,ERROR))// &
                & " for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
                & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !component_idx
          CALL COORDINATE_INTERPOLATION_PARAMETERS_ADJUST(COORDINATE_SYSTEM,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="The field parameter set number of "//TRIM(NUMBER_TO_VSTRING(PARAMETER_SET_NUMBER,"*",ERR,ERROR))// &
            & " has not been created for field number "// &
            & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The field parameter set number of "//TRIM(NUMBER_TO_VSTRING(PARAMETER_SET_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The number must be between 1 and "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
          & PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS,"*",ERR,ERROR))//" for field number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interpolation parameters is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interpolation parameters:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field number = ",INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field variable number = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & VARIABLE_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parameter set number = ",PARAMETER_SET_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Line number = ",LINE_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of components = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Component = ",component_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of parameters = ",INTERPOLATION_PARAMETERS% &
          & NUMBER_OF_PARAMETERS(component_idx),ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx),4,4, &
          & INTERPOLATION_PARAMETERS%PARAMETERS(:,component_idx),'("      Parameters :",4(X,E13.6))','(18X,4(X,E13.6))', &
          & ERR,ERROR,*999)
      ENDDO !component_idx
    ENDIF

    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_LINE_GET")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATION_PARAMETERS_LINE_GET",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_LINE_GET")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_LINE_GET

  !
  !================================================================================================================================
  !

  !>Gets the interpolation scale factors for a particular element. 
  SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,INTERPOLATION_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to get the element scale factors for
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<On return, a pointer to the interpolation parameters
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,mk,nk,nn,np,ns,ny,scaling_idx
    REAL(DP), POINTER :: SCALE_FACTORS(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
      SELECT CASE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALING_TYPE)
      CASE(FIELD_NO_SCALING)
        CALL FLAG_ERROR("Can not get the scale factors for a field with no scaling.",ERR,ERROR,*999)
      CASE(FIELD_UNIT_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
        DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          ELEMENTS_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
          IF(ELEMENT_NUMBER>0.AND.ELEMENT_NUMBER<=ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS) THEN
            BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
            INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR=>BASIS
            SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)             
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              NODES_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
              scaling_idx=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%SCALING_INDEX
              NULLIFY(SCALE_FACTORS)
              CALL DISTRIBUTED_VECTOR_DATA_GET(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
              DO nn=1,BASIS%NUMBER_OF_NODES
                np=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nn)
                DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                  nk=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(mk,nn)
                  ns=BASIS%ELEMENT_PARAMETER_INDEX(nk,nn)
                  ny=NODES_TOPOLOGY%NODES(np)%DOF_INDEX(nk)
                  INTERPOLATION_PARAMETERS%SCALE_FACTORS(ns,component_idx)=SCALE_FACTORS(ny)
                ENDDO !mk
              ENDDO !nn
              CALL DISTRIBUTED_VECTOR_DATA_RESTORE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE DEFAULT
              LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
              & COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
              & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The element number of "//TRIM(NUMBER_TO_VSTRING(ELEMENT_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS_TOPOLOGY%TOTAL_NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
              & " for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !component_idx
      CASE(FIELD_ARC_LENGTH_SCALING)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The scaling type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%SCALINGS% &
          & SCALING_TYPE,"*",ERR,ERROR))//" is invalid for field number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Interpolation parameters is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interpolation scale factors:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field number = ",INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field variable number = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & VARIABLE_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of components = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Component = ",component_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of parameters = ",INTERPOLATION_PARAMETERS% &
          & NUMBER_OF_PARAMETERS(component_idx),ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx),4,4, &
          & INTERPOLATION_PARAMETERS%SCALE_FACTORS(:,component_idx),'("      Scale factors :",4(X,E13.6))','(21X,4(X,E13.6))', &
          & ERR,ERROR,*999)
      ENDDO !component_idx
    ENDIF

    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET

  !
  !================================================================================================================================
  !

  !>Gets the interpolation scale factors for a particular element. 
  SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET(LINE_NUMBER,INTERPOLATION_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: LINE_NUMBER !<The line number to get the element scale factors for
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS !<On return, a pointer to the interpolation parameters
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,mk,nk,nn,np,ns,ny,scaling_idx
    REAL(DP), POINTER :: SCALE_FACTORS(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_LINES_TYPE), POINTER :: LINES_TOPOLOGY
    TYPE(DOMAIN_NODES_TYPE), POINTER :: NODES_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) THEN
      SELECT CASE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALING_TYPE)
      CASE(FIELD_NO_SCALING)
        CALL FLAG_ERROR("Can not scale factors for a field with no scaling.",ERR,ERROR,*999)
      CASE(FIELD_UNIT_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
        DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          LINES_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%LINES
          IF(LINE_NUMBER>0.AND.LINE_NUMBER<=LINES_TOPOLOGY%NUMBER_OF_LINES) THEN
            BASIS=>LINES_TOPOLOGY%LINES(LINE_NUMBER)%BASIS
            INTERPOLATION_PARAMETERS%BASES(component_idx)%PTR=>BASIS
            SELECT CASE(INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
            CASE(FIELD_CONSTANT_INTERPOLATION)             
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              NODES_TOPOLOGY=>INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%NODES
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
              scaling_idx=INTERPOLATION_PARAMETERS%FIELD_VARIABLE%COMPONENTS(component_idx)%SCALING_INDEX
              NULLIFY(SCALE_FACTORS)
              CALL DISTRIBUTED_VECTOR_DATA_GET(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
              DO nn=1,BASIS%NUMBER_OF_NODES
                np=LINES_TOPOLOGY%LINES(LINE_NUMBER)%NODES_IN_LINE(nn)
                DO mk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                  nk=LINES_TOPOLOGY%LINES(LINE_NUMBER)%DERIVATIVES_IN_LINE(mk,nn)
                  ns=BASIS%ELEMENT_PARAMETER_INDEX(nk,nn)
                  ny=NODES_TOPOLOGY%NODES(np)%DOF_INDEX(nk)
                  INTERPOLATION_PARAMETERS%SCALE_FACTORS(ns,component_idx)=SCALE_FACTORS(ny)
                ENDDO !mk
              ENDDO !nn
              CALL DISTRIBUTED_VECTOR_DATA_RESTORE(INTERPOLATION_PARAMETERS%FIELD%SCALINGS%SCALINGS(scaling_idx)% &
                & SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
            CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
              INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx)=1
              INTERPOLATION_PARAMETERS%PARAMETERS(1,component_idx)=1.0_DP
            CASE DEFAULT
              LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
              & COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
              & TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            LOCAL_ERROR="The line number of "//TRIM(NUMBER_TO_VSTRING(LINE_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(LINES_TOPOLOGY%NUMBER_OF_LINES,"*",ERR,ERROR))// &
              & " for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))//" of field number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ENDDO !component_idx
      CASE(FIELD_ARC_LENGTH_SCALING)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="The scaling type of "//TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%SCALINGS% &
          & SCALING_TYPE,"*",ERR,ERROR))//" is invalid for field number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Interpolation parameters is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Interpolation scale factors:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field number = ",INTERPOLATION_PARAMETERS%FIELD%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Field variable number = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & VARIABLE_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Line number = ",LINE_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of components = ",INTERPOLATION_PARAMETERS%FIELD_VARIABLE% &
        & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      DO component_idx=1,INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Component = ",component_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of parameters = ",INTERPOLATION_PARAMETERS% &
          & NUMBER_OF_PARAMETERS(component_idx),ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,INTERPOLATION_PARAMETERS%NUMBER_OF_PARAMETERS(component_idx),4,4, &
          & INTERPOLATION_PARAMETERS%SCALE_FACTORS(:,component_idx),'("      Scale factors :",4(X,E13.6))','(21X,4(X,E13.6))', &
          & ERR,ERROR,*999)
      ENDDO !component_idx
    ENDIF

    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET")
    RETURN
999 CALL ERRORS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET",ERR,ERROR)
    CALL EXITS("FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET")
    RETURN 1
  END SUBROUTINE FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_LINE_GET

  !
  !================================================================================================================================
  !

  !>Calculates the mappings for a field.
  SUBROUTINE FIELD_MAPPINGS_CALCULATE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to calculate the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx,component_idx,VARIABLE_GLOBAL_DOFS_OFFSET,NUMBER_OF_GLOBAL_VARIABLE_DOFS, &
      & NUMBER_OF_CONSTANT_DOFS,NUMBER_OF_ELEMENT_DOFS,NUMBER_OF_NODE_DOFS,NUMBER_OF_GRID_POINT_DOFS,NUMBER_OF_GAUSS_POINT_DOFS, &
      & NUMBER_OF_LOCAL_VARIABLE_DOFS,TOTAL_NUMBER_OF_VARIABLE_DOFS,NUMBER_OF_DOMAINS,mesh_component_idx,variable_global_ny, &
      & variable_local_ny,domain_idx,domain_no,constant_nyy,element_nyy,node_nyy,grid_point_nyy,Gauss_point_nyy, &
      & MAX_NUMBER_OF_DERIVATIVES,ne,nk,np,ny,NUMBER_OF_COMPUTATIONAL_NODES,my_computational_node_number
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_LOCAL_DOFS_OFFSETS(:)
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING,DOFS_MAPPING,FIELD_VARIABLE_DOFS_MAPPING
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(FIELD_VARIABLE_COMPONENT_TYPE), POINTER :: FIELD_COMPONENT
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_MAPPINGS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      NUMBER_OF_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
      IF(ERR/=0) GOTO 999
      my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      IF(ERR/=0) GOTO 999
      !Calculate the number of global and local degrees of freedom for the field variables and components
      DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
        NUMBER_OF_CONSTANT_DOFS=0
        NUMBER_OF_ELEMENT_DOFS=0
        NUMBER_OF_NODE_DOFS=0
        NUMBER_OF_GRID_POINT_DOFS=0
        NUMBER_OF_GAUSS_POINT_DOFS=0
        NUMBER_OF_LOCAL_VARIABLE_DOFS=0
        TOTAL_NUMBER_OF_VARIABLE_DOFS=0
        NUMBER_OF_GLOBAL_VARIABLE_DOFS=0
        DO component_idx=1,FIELD%VARIABLES(variable_idx)%NUMBER_OF_COMPONENTS
          FIELD_COMPONENT=>FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)
          SELECT CASE(FIELD_COMPONENT%INTERPOLATION_TYPE)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            NUMBER_OF_CONSTANT_DOFS=NUMBER_OF_CONSTANT_DOFS+1
            NUMBER_OF_LOCAL_VARIABLE_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS+1
            TOTAL_NUMBER_OF_VARIABLE_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS+1
            NUMBER_OF_GLOBAL_VARIABLE_DOFS=NUMBER_OF_GLOBAL_VARIABLE_DOFS+1
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            DOMAIN=>FIELD_COMPONENT%DOMAIN
            mesh_component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            MESH=>DOMAIN%MESH
            MESH_TOPOLOGY=>MESH%TOPOLOGY(mesh_component_idx)%PTR
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            NUMBER_OF_ELEMENT_DOFS=NUMBER_OF_ELEMENT_DOFS+DOMAIN_TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
            NUMBER_OF_LOCAL_VARIABLE_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS+DOMAIN_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
            TOTAL_NUMBER_OF_VARIABLE_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS+DOMAIN_TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
            NUMBER_OF_GLOBAL_VARIABLE_DOFS=NUMBER_OF_GLOBAL_VARIABLE_DOFS+MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            DOMAIN=>FIELD_COMPONENT%DOMAIN
            mesh_component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            MESH=>DOMAIN%MESH
            MESH_TOPOLOGY=>MESH%TOPOLOGY(mesh_component_idx)%PTR
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            NUMBER_OF_NODE_DOFS=NUMBER_OF_NODE_DOFS+DOMAIN_TOPOLOGY%DOFS%TOTAL_NUMBER_OF_DOFS
            NUMBER_OF_LOCAL_VARIABLE_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS+DOMAIN_TOPOLOGY%DOFS%NUMBER_OF_DOFS
            TOTAL_NUMBER_OF_VARIABLE_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS+DOMAIN_TOPOLOGY%DOFS%TOTAL_NUMBER_OF_DOFS
            NUMBER_OF_GLOBAL_VARIABLE_DOFS=NUMBER_OF_GLOBAL_VARIABLE_DOFS+MESH_TOPOLOGY%DOFS%NUMBER_OF_DOFS
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interpolation type of "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)%INTERPOLATION_TYPE, &
              & "*",ERR,ERROR))//" is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " of variable type  "//TRIM(NUMBER_TO_VSTRING(FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDDO !component_idx
        FIELD%VARIABLES(variable_idx)%NUMBER_OF_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS
        FIELD%VARIABLES(variable_idx)%TOTAL_NUMBER_OF_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS
        FIELD%VARIABLES(variable_idx)%NUMBER_OF_GLOBAL_DOFS=NUMBER_OF_GLOBAL_VARIABLE_DOFS
        ALLOCATE(FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(2,TOTAL_NUMBER_OF_VARIABLE_DOFS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof to parameter map.",ERR,ERROR,*999)
        FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NUMBER_OF_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS
        IF(NUMBER_OF_CONSTANT_DOFS>0) THEN
          ALLOCATE(FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%CONSTANT_DOF2PARAM_MAP(NUMBER_OF_CONSTANT_DOFS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof to parameter constant map.",ERR,ERROR,*999)
          FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NUMBER_OF_CONSTANT_DOFS=NUMBER_OF_CONSTANT_DOFS
        ENDIF
        IF(NUMBER_OF_ELEMENT_DOFS>0) THEN
          ALLOCATE(FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP(2,NUMBER_OF_ELEMENT_DOFS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof to parameter element map.",ERR,ERROR,*999)
          FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NUMBER_OF_ELEMENT_DOFS=NUMBER_OF_ELEMENT_DOFS
        ENDIF
        IF(NUMBER_OF_NODE_DOFS>0) THEN
          ALLOCATE(FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,NUMBER_OF_NODE_DOFS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof to parameter node map.",ERR,ERROR,*999)
          FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NUMBER_OF_NODE_DOFS=NUMBER_OF_NODE_DOFS
        ENDIF
        IF(NUMBER_OF_GRID_POINT_DOFS>0) THEN
          ALLOCATE(FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%GRID_POINT_DOF2PARAM_MAP(2,NUMBER_OF_GRID_POINT_DOFS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof to parameter grid point map.",ERR,ERROR,*999)
          FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NUMBER_OF_GRID_POINT_DOFS=NUMBER_OF_GRID_POINT_DOFS
        ENDIF
        IF(NUMBER_OF_GAUSS_POINT_DOFS>0) THEN
          ALLOCATE(FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP(3,NUMBER_OF_GAUSS_POINT_DOFS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof to parameter Gauss point map.",ERR,ERROR,*999)
          FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NUMBER_OF_GAUSS_POINT_DOFS=NUMBER_OF_GAUSS_POINT_DOFS
        ENDIF
      ENDDO !variable_idx
      !Allocate the mapping arrays
      DECOMPOSITION=>FIELD%DECOMPOSITION
      ALLOCATE(VARIABLE_LOCAL_DOFS_OFFSETS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable local dofs offsets.",ERR,ERROR,*999)
      !Calculate the local and global numbers and set up the mappings
      DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
        constant_nyy=0
        element_nyy=0
        node_nyy=0
        grid_point_nyy=0
        Gauss_point_nyy=0
        NUMBER_OF_LOCAL_VARIABLE_DOFS=0
        TOTAL_NUMBER_OF_VARIABLE_DOFS=0
        VARIABLE_GLOBAL_DOFS_OFFSET=0
        VARIABLE_LOCAL_DOFS_OFFSETS=0
        FIELD_VARIABLE_DOFS_MAPPING=>FIELD%VARIABLES(variable_idx)%DOMAIN_MAPPING
        IF(ASSOCIATED(FIELD_VARIABLE_DOFS_MAPPING)) THEN
          ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(FIELD%VARIABLES(variable_idx)%NUMBER_OF_GLOBAL_DOFS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable dofs mapping global to local map.",ERR,ERROR,*999)
          FIELD_VARIABLE_DOFS_MAPPING%NUMBER_OF_GLOBAL=FIELD%VARIABLES(variable_idx)%NUMBER_OF_GLOBAL_DOFS
        ENDIF
        DO component_idx=1,FIELD%VARIABLES(variable_idx)%NUMBER_OF_COMPONENTS
          FIELD_COMPONENT=>FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)
          SELECT CASE(FIELD_COMPONENT%INTERPOLATION_TYPE)
          CASE(FIELD_CONSTANT_INTERPOLATION)
            variable_local_ny=1+VARIABLE_LOCAL_DOFS_OFFSETS(my_computational_node_number)
            NUMBER_OF_LOCAL_VARIABLE_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS+1
            TOTAL_NUMBER_OF_VARIABLE_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS+1
            !Allocate and set up global to local domain map for variable mapping
            IF(ASSOCIATED(FIELD_VARIABLE_DOFS_MAPPING)) THEN
              variable_global_ny=1+VARIABLE_GLOBAL_DOFS_OFFSET
              CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny), &
                & ERR,ERROR,*999)
              NUMBER_OF_DOMAINS=NUMBER_OF_COMPUTATIONAL_NODES !Constant is in all domains
              ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_NUMBER(NUMBER_OF_DOMAINS), &
                & STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map local number.",ERR,ERROR,*999)
              ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS), &
                & STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map domain number.",ERR,ERROR,*999)
              ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map domain number.",ERR,ERROR,*999)
              !A constant dof is mapped to all domains.
              FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
              DO domain_idx=1,NUMBER_OF_DOMAINS
                domain_no=domain_idx-1
                FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_NUMBER(domain_idx)= &
                  & 1+VARIABLE_LOCAL_DOFS_OFFSETS(domain_no)
                FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%DOMAIN_NUMBER(domain_idx)=domain_no
                FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_INTERNAL
              ENDDO !domain_idx
            ENDIF
            constant_nyy=constant_nyy+1
            !Setup dof to parameter map
            FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(1,variable_local_ny)=FIELD_CONSTANT_DOF_TYPE
            FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(2,variable_local_ny)=constant_nyy
            FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%CONSTANT_DOF2PARAM_MAP(constant_nyy)=component_idx
            !Setup reverse parameter to dof map
            FIELD_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS=1
            FIELD_COMPONENT%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP=variable_local_ny
            !Adjust the offsets
            VARIABLE_GLOBAL_DOFS_OFFSET=VARIABLE_GLOBAL_DOFS_OFFSET+1
            VARIABLE_LOCAL_DOFS_OFFSETS=VARIABLE_LOCAL_DOFS_OFFSETS+1
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            DOMAIN=>FIELD_COMPONENT%DOMAIN
            ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
            mesh_component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            MESH=>DOMAIN%MESH
            MESH_TOPOLOGY=>MESH%TOPOLOGY(mesh_component_idx)%PTR
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            NUMBER_OF_ELEMENT_DOFS=NUMBER_OF_ELEMENT_DOFS+DOMAIN_TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
            !Allocate parameter to dof map for this field variable component
            DOFS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
            ALLOCATE(FIELD_COMPONENT%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP(DOMAIN_TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS), &
              & STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field component parameter to dof element map.",ERR,ERROR,*999)
            FIELD_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS=DOMAIN_TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
            !Handle global dofs domain mapping
            DO ny=1,ELEMENTS_MAPPING%NUMBER_OF_GLOBAL
              !Handle field mappings
              TOTAL_NUMBER_OF_VARIABLE_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS+1
              !Handle field variable mappings
              IF(ASSOCIATED(FIELD_VARIABLE_DOFS_MAPPING)) THEN
                variable_global_ny=ny+VARIABLE_GLOBAL_DOFS_OFFSET
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(FIELD_VARIABLE_DOFS_MAPPING% &
                  & GLOBAL_TO_LOCAL_MAP(variable_global_ny),ERR,ERROR,*999)
                NUMBER_OF_DOMAINS=DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS
                ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_NUMBER(NUMBER_OF_DOMAINS), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map local number.", &
                  & ERR,ERROR,*999)
                ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map domain number.", &
                  & ERR,ERROR,*999)
                ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_TYPE(NUMBER_OF_DOMAINS), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map domain number.", &
                  & ERR,ERROR,*999)
                FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                DO domain_idx=1,NUMBER_OF_DOMAINS
                  domain_no=DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)
                  FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_NUMBER(domain_idx)= &
                    & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)+VARIABLE_LOCAL_DOFS_OFFSETS(domain_no)
                  FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%DOMAIN_NUMBER(domain_idx)= &
                    & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)
                  FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_TYPE(domain_idx)= &
                    & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)
                ENDDO !domain_idx
              ENDIF
            ENDDO !ny
            !Handle local dofs domain mapping
            DO ne=1,ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL
              variable_local_ny=ne+VARIABLE_LOCAL_DOFS_OFFSETS(my_computational_node_number)
              NUMBER_OF_LOCAL_VARIABLE_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS+1
              element_nyy=element_nyy+1
              !Setup dof to parameter map
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(1,variable_local_ny)=FIELD_ELEMENT_DOF_TYPE
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(2,variable_local_ny)=element_nyy
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP(1,element_nyy)=ne
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP(2,element_nyy)=component_idx
              !Setup reverse parameter to dof map
              FIELD_COMPONENT%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP(ne)=variable_local_ny
            ENDDO !ne
            !Adjust the offsets
            VARIABLE_GLOBAL_DOFS_OFFSET=VARIABLE_GLOBAL_DOFS_OFFSET+ELEMENTS_MAPPING%NUMBER_OF_GLOBAL
            VARIABLE_LOCAL_DOFS_OFFSETS=VARIABLE_LOCAL_DOFS_OFFSETS+ELEMENTS_MAPPING%NUMBER_OF_DOMAIN_LOCAL
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            DOMAIN=>FIELD_COMPONENT%DOMAIN
            mesh_component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            MESH=>DOMAIN%MESH
            MESH_TOPOLOGY=>MESH%TOPOLOGY(mesh_component_idx)%PTR
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            NUMBER_OF_NODE_DOFS=NUMBER_OF_NODE_DOFS+DOMAIN_TOPOLOGY%DOFS%TOTAL_NUMBER_OF_DOFS
            DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
            !Allocate parameter to dof map for this field variable component
            MAX_NUMBER_OF_DERIVATIVES=-1
            DO np=1,DOMAIN_TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
              IF(DOMAIN_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES>MAX_NUMBER_OF_DERIVATIVES) &
                & MAX_NUMBER_OF_DERIVATIVES=DOMAIN_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
            ENDDO !np
            ALLOCATE(FIELD_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(MAX_NUMBER_OF_DERIVATIVES, &
              & DOMAIN_TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field component parameter to dof node map.",ERR,ERROR,*999)
            FIELD_COMPONENT%PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS=DOMAIN_TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
            FIELD_COMPONENT%PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES=MAX_NUMBER_OF_DERIVATIVES
            !Handle global dofs domain mapping
            DO ny=1,DOFS_MAPPING%NUMBER_OF_GLOBAL
              !Handle field mapping
              TOTAL_NUMBER_OF_VARIABLE_DOFS=TOTAL_NUMBER_OF_VARIABLE_DOFS+1
              !Handle variable mapping
              IF(ASSOCIATED(FIELD_VARIABLE_DOFS_MAPPING)) THEN
                variable_global_ny=ny+VARIABLE_GLOBAL_DOFS_OFFSET
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(FIELD_VARIABLE_DOFS_MAPPING% &
                  & GLOBAL_TO_LOCAL_MAP(variable_global_ny),ERR,ERROR,*999)
                NUMBER_OF_DOMAINS=DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS
                ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_NUMBER(NUMBER_OF_DOMAINS), &
                  & STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable dofs global to local map local number.", &
                  & ERR,ERROR,*999)
                ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS), &
                  & STAT=ERR)
                IF(ERR/=0) &
                  & CALL FLAG_ERROR("Could not allocate field variable dofs global to local map domain number.",ERR,ERROR,*999)
                ALLOCATE(FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) &
                  & CALL FLAG_ERROR("Could not allocate field variable dofs global to local map domain number.",ERR,ERROR,*999)
                FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                DO domain_idx=1,NUMBER_OF_DOMAINS
                  domain_no=DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)
                  FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_NUMBER(domain_idx)= &
                    & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)+VARIABLE_LOCAL_DOFS_OFFSETS(domain_no)
                  FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%DOMAIN_NUMBER(domain_idx)= &
                    & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)
                  FIELD_VARIABLE_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(variable_global_ny)%LOCAL_TYPE(domain_idx)= &
                    & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)
                ENDDO !domain_idx
              ENDIF
            ENDDO !ny (global)
            !Handle local dofs domain mapping
            DO ny=1,DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
              variable_local_ny=ny+VARIABLE_LOCAL_DOFS_OFFSETS(my_computational_node_number)
              NUMBER_OF_LOCAL_VARIABLE_DOFS=NUMBER_OF_LOCAL_VARIABLE_DOFS+1
              node_nyy=node_nyy+1
              nk=DOMAIN%TOPOLOGY%DOFS%DOF_INDEX(1,ny)
              np=DOMAIN%TOPOLOGY%DOFS%DOF_INDEX(2,ny)
              !Setup dof to parameter map
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(1,variable_local_ny)=FIELD_NODE_DOF_TYPE
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%DOF_TYPE(2,variable_local_ny)=node_nyy
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(1,node_nyy)=nk
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(2,node_nyy)=np
              FIELD%VARIABLES(variable_idx)%DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP(3,node_nyy)=component_idx
              !Setup reverse parameter to dof map
              FIELD_COMPONENT%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(nk,np)=variable_local_ny
            ENDDO !ny
            !Adjust the offsets
            VARIABLE_GLOBAL_DOFS_OFFSET=VARIABLE_GLOBAL_DOFS_OFFSET+DOFS_MAPPING%NUMBER_OF_GLOBAL
            VARIABLE_LOCAL_DOFS_OFFSETS=VARIABLE_LOCAL_DOFS_OFFSETS+DOFS_MAPPING%NUMBER_OF_DOMAIN_LOCAL
          CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interpolation type of "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)%INTERPOLATION_TYPE, &
              & "*",ERR,ERROR))//" is invalid for component number "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " of variable type "//TRIM(NUMBER_TO_VSTRING(FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDDO !component_idx
        IF(ASSOCIATED(FIELD_VARIABLE_DOFS_MAPPING)) THEN
          CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(FIELD_VARIABLE_DOFS_MAPPING,ERR,ERROR,*999)
        ENDIF
      ENDDO !variable_idx

      IF(ALLOCATED(VARIABLE_LOCAL_DOFS_OFFSETS)) DEALLOCATE(VARIABLE_LOCAL_DOFS_OFFSETS)

    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_MAPPINGS_CALCULATE")
    RETURN
999 IF(ALLOCATED(VARIABLE_LOCAL_DOFS_OFFSETS)) DEALLOCATE(VARIABLE_LOCAL_DOFS_OFFSETS)
    CALL ERRORS("FIELD_MAPPINGS_CALCULATE",ERR,ERROR)
    CALL EXITS("FIELD_MAPPINGS_CALCULATE")
    RETURN 1
  END SUBROUTINE FIELD_MAPPINGS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the dofs to parameters mapping for a field varaible and deallocates all memory.
  SUBROUTINE FIELD_DOF_TO_PARAM_MAP_FINALISE(DOF_TO_PARAM_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE) :: DOF_TO_PARAM_MAP !<The dof to parameter map to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_DOF_TO_PARAM_MAP_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(DOF_TO_PARAM_MAP%DOF_TYPE)) DEALLOCATE(DOF_TO_PARAM_MAP%DOF_TYPE)
    IF(ALLOCATED(DOF_TO_PARAM_MAP%CONSTANT_DOF2PARAM_MAP)) DEALLOCATE(DOF_TO_PARAM_MAP%CONSTANT_DOF2PARAM_MAP)
    IF(ALLOCATED(DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP)) DEALLOCATE(DOF_TO_PARAM_MAP%ELEMENT_DOF2PARAM_MAP)
    IF(ALLOCATED(DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP)) DEALLOCATE(DOF_TO_PARAM_MAP%NODE_DOF2PARAM_MAP)
    IF(ALLOCATED(DOF_TO_PARAM_MAP%GRID_POINT_DOF2PARAM_MAP)) DEALLOCATE(DOF_TO_PARAM_MAP%GRID_POINT_DOF2PARAM_MAP)
    IF(ALLOCATED(DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP)) DEALLOCATE(DOF_TO_PARAM_MAP%GAUSS_POINT_DOF2PARAM_MAP)
    DOF_TO_PARAM_MAP%NUMBER_OF_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_CONSTANT_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_ELEMENT_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_NODE_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_GRID_POINT_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_GAUSS_POINT_DOFS=0

    CALL EXITS("FIELD_DOF_TO_PARAM_MAP_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_DOF_TO_PARAM_MAP_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_DOF_TO_PARAM_FMAP_INALISE")
    RETURN 1
  END SUBROUTINE FIELD_DOF_TO_PARAM_MAP_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the dofs to parameters mappings for a field.
  SUBROUTINE FIELD_DOF_TO_PARAM_MAP_INITIALISE(DOF_TO_PARAM_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_DOF_TO_PARAM_MAP_TYPE) :: DOF_TO_PARAM_MAP !<The dof to parameter map to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_DOF_TO_PARAM_INITIALISE",ERR,ERROR,*999)

    DOF_TO_PARAM_MAP%NUMBER_OF_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_CONSTANT_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_ELEMENT_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_NODE_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_GRID_POINT_DOFS=0
    DOF_TO_PARAM_MAP%NUMBER_OF_GAUSS_POINT_DOFS=0

    CALL EXITS("FIELD_DOF_TO_PARAM_MAP_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_DOF_TO_PARAM_MAP_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_DOF_TO_PARAM_MAP_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_DOF_TO_PARAM_MAP_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the geometric field for a field identified by a pointer.
  SUBROUTINE FIELD_GEOMETRIC_FIELD_GET(FIELD,GEOMETRIC_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the geometric field for
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<On return, a pointer to the geometric field. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_GEOMETRIC_FIELD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
          CALL FLAG_ERROR("Geometric field is already associated.",ERR,ERROR,*999)
        ELSE
          GEOMETRIC_FIELD=>FIELD%GEOMETRIC_FIELD
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_FIELD_GET")
    RETURN
999 CALL ERRORS("FIELD_GEOMETRIC_FIELD_GET",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_FIELD_GET")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_FIELD_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field.
  SUBROUTINE FIELD_GEOMETRIC_FIELD_SET(FIELD,GEOMETRIC_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the geometric field for
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_GEOMETRIC_FIELD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%GEOMETRIC_FIELD_LOCKED) THEN
            LOCAL_ERROR="The geometric field has been locked for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(FIELD%GEOMETRIC_FIELD)) THEN
              LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                & " already has a geometric field associated."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              IF(ASSOCIATED(FIELD%DECOMPOSITION)) THEN
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(GEOMETRIC_FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
                    IF(GEOMETRIC_FIELD%FIELD_FINISHED) THEN
                      IF(FIELD%DECOMPOSITION%MESH%USER_NUMBER==GEOMETRIC_FIELD%DECOMPOSITION%MESH%USER_NUMBER) THEN
                        SELECT CASE(FIELD%TYPE)
                        CASE(FIELD_FIBRE_TYPE,FIELD_GENERAL_TYPE,FIELD_MATERIAL_TYPE)
                          FIELD%GEOMETRIC_FIELD=>GEOMETRIC_FIELD
                        CASE(FIELD_GEOMETRIC_TYPE)
                          CALL FLAG_ERROR("Can not set the geometric field for a geometric field.",ERR,ERROR,*999)
                        CASE DEFAULT
                          LOCAL_ERROR="The field type "//TRIM(NUMBER_TO_VSTRING(FIELD%TYPE,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        LOCAL_ERROR="The specified field is decomposed on mesh user number "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD%DECOMPOSITION%MESH%USER_NUMBER,"*",ERR,ERROR))// &
                          & " and the geometric field is decomposed on mesh user number "// &
                          & TRIM(NUMBER_TO_VSTRING(GEOMETRIC_FIELD%DECOMPOSITION%MESH%USER_NUMBER,"*",ERR,ERROR))// &
                          & ". The two fields must use the same mesh."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("The specified geometric field has not been finished.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The specified geometric field is not a geometric field.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The field does not have a decomposition associated.",ERR,ERROR,*999)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_FIELD_SET")
    RETURN
999 CALL ERRORS("FIELD_GEOMETRIC_FIELD_SET",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_FIELD_SET")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_FIELD_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the geometric field for a field and locks so that no further changes can be made.
  SUBROUTINE FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(FIELD,GEOMETRIC_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the geometric field for
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_GEOMETRIC_FIELD_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_GEOMETRIC_FIELD_SET(FIELD,GEOMETRIC_FIELD,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%GEOMETRIC_FIELD_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_GEOMETRIC_FIELD_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_GEOMETRIC_FIELD_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_FIELD_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_FIELD_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Calculates the geometric parameters (line lengths, areas, volumes, scaling etc.) for a field.
  SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_CALCULATE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<The field to calculate the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
          CALL FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE(FIELD,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN
999 CALL ERRORS("FIELD_GEOMETRIC_PARAMETERS_CALCULATE",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_CALCULATE")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the geometric parameters and deallocates all memory.
  SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_FINALISE(GEOMETRIC_PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_GEOMETRIC_PARAMETERS_TYPE), POINTER :: GEOMETRIC_PARAMETERS !<A pointer to the geometric field parameters to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx
    TYPE(FIELD_TYPE), POINTER :: FIELD2

    CALL ENTERS("FIELD_GEOMETRIC_PARAMETERS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(GEOMETRIC_PARAMETERS)) THEN
      !Nullify the geometric field pointer of those fields using this geometric field.
      DO field_idx=1,GEOMETRIC_PARAMETERS%NUMBER_OF_FIELDS_USING
        FIELD2=>GEOMETRIC_PARAMETERS%FIELDS_USING(field_idx)%PTR
        IF(ASSOCIATED(FIELD2)) NULLIFY(FIELD2%GEOMETRIC_FIELD)
      ENDDO !field_idx
      IF(ASSOCIATED(GEOMETRIC_PARAMETERS%FIELDS_USING)) DEALLOCATE(GEOMETRIC_PARAMETERS%FIELDS_USING)
      IF(ALLOCATED(GEOMETRIC_PARAMETERS%LENGTHS)) DEALLOCATE(GEOMETRIC_PARAMETERS%LENGTHS)
      IF(ALLOCATED(GEOMETRIC_PARAMETERS%AREAS)) DEALLOCATE(GEOMETRIC_PARAMETERS%AREAS)
      IF(ALLOCATED(GEOMETRIC_PARAMETERS%VOLUMES)) DEALLOCATE(GEOMETRIC_PARAMETERS%VOLUMES)
      DEALLOCATE(GEOMETRIC_PARAMETERS)
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_GEOMETRIC_PARAMETERS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the geometric parameters for a geometric field
  SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_INITIALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the geometric parameters for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx
    TYPE(FIELD_PTR_TYPE), POINTER :: NEW_FIELDS_USING(:)

    NULLIFY(NEW_FIELDS_USING)

    CALL ENTERS("FIELD_GEOMETRIC_PARAMETERS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
        !Field is a geometric field
        ALLOCATE(FIELD%GEOMETRIC_FIELD_PARAMETERS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate geometric field parameters.",ERR,ERROR,*999)
        FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_LINES=FIELD%DECOMPOSITION%TOPOLOGY%LINES%NUMBER_OF_LINES
        FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_AREAS=0
        FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_VOLUMES=0
        ALLOCATE(FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_LINES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate lengths.",ERR,ERROR,*999)
        FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS=0.0_DP
        !The field is a geometric field so it must use itself initiallly
        ALLOCATE(FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate fields using.",ERR,ERROR,*999)
        FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(1)%PTR=>FIELD
        FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING=1
      ELSE
        !Field is not a geometric field
        NULLIFY(FIELD%GEOMETRIC_FIELD_PARAMETERS)
        IF(ASSOCIATED(FIELD%GEOMETRIC_FIELD)) THEN
          !Set the geometric field so that it knows that this field is using it
          ALLOCATE(NEW_FIELDS_USING(FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING+1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new fields using.",ERR,ERROR,*999)
          DO field_idx=1,FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING
            NEW_FIELDS_USING(field_idx)%PTR=>FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(field_idx)%PTR
          ENDDO !field_idx
          NEW_FIELDS_USING(FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING+1)%PTR=>FIELD
          FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING=FIELD%GEOMETRIC_FIELD% &
            & GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING+1
          IF(ASSOCIATED(FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING)) &
            & DEALLOCATE(FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING)
          FIELD%GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING=>NEW_FIELDS_USING
        ELSE
          CALL FLAG_ERROR("Field does not have a geometric field associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_INITIALISE")
    RETURN
999 IF(ASSOCIATED(NEW_FIELDS_USING)) DEALLOCATE(NEW_FIELDS_USING)
    CALL ERRORS("FIELD_GEOMETRIC_PARAMETERS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the line lengths from the parameters of a geometric field. Old CMISS name LINSCA
  SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to calculate the line lengths for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,ITERATION_NUMBER,MAXIMUM_DIFFERENCE_LINE,ng,nl
    INTEGER(INTG), PARAMETER :: LINES_MAXIMUM_NUMBER_OF_ITERATIONS=20
    INTEGER(INTG) :: GAUSS_START(4) = (/ 0,1,3,6 /)
    INTEGER(INTG) :: NUMBER_OF_GAUSS_POINTS=4
    REAL(DP) :: LAST_MAXIMUM_LENGTH_DIFFERENCE,LENGTH_DIFFERENCE,MAXIMUM_LENGTH_DIFFERENCE,XI(1),W,DERIV_NORM,LINE_LENGTH, &
      & OLD_LINE_LENGTH
! Doxygen doesn't like this
!    REAL(DP) :: XIG(10) = (/ 0.500000000000000_DP, &
!      &                      0.211324865405187_DP,0.788675134594813_DP, &
!      &                      0.112701665379258_DP,0.500000000000000_DP,0.887298334620742_DP, &
!      &                      0.06943184420297349_DP,0.330009478207572_DP,0.669990521792428_DP,0.930568155797026_DP /)
!    REAL(DP) :: WIG(10) = (/ 1.000000000000000_DP, &
!      &                      0.500000000000000_DP,0.500000000000000_DP, &
!      &                      0.277777777777778_DP,0.444444444444444_DP,0.277777777777778_DP,
!      &                      0.173927422568727_DP,0.326072577431273_DP,0.326072577431273_DP,0.173927422568727_DP /)
    REAL(DP) :: XIG(10),WIG(10)
    REAL(DP), PARAMETER :: LINE_INCREMENT_TOLERANCE=CONVERGENCE_TOLERANCE
    LOGICAL :: ITERATE,UPDATE_FIELDS_USING
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    XIG = (/ 0.500000000000000_DP, &
      &      0.211324865405187_DP,0.788675134594813_DP, &
      &      0.112701665379258_DP,0.500000000000000_DP,0.887298334620742_DP, &
      &      0.06943184420297349_DP,0.330009478207572_DP,0.669990521792428_DP,0.930568155797026_DP /)
    WIG = (/ 1.000000000000000_DP, &
      &      0.500000000000000_DP,0.500000000000000_DP, &
      &      0.277777777777778_DP,0.444444444444444_DP,0.277777777777778_DP, &
      &      0.173927422568727_DP,0.326072577431273_DP,0.326072577431273_DP,0.173927422568727_DP /)

    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)

    CALL ENTERS("FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
          IF(ASSOCIATED(FIELD%GEOMETRIC_FIELD_PARAMETERS)) THEN
            COORDINATE_SYSTEM=>FIELD%REGION%COORDINATE_SYSTEM
            !Iterate to find the line lengths as the line lengths depend on the scaling factors and vise versa.
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(FIELD,FIELD_U_VARIABLE_TYPE,INTERPOLATION_PARAMETERS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINT,ERR,ERROR,*999)
            ITERATE=.TRUE.
            ITERATION_NUMBER=0
            LAST_MAXIMUM_LENGTH_DIFFERENCE=0.0_DP
            DO WHILE(ITERATE.AND.ITERATION_NUMBER<=LINES_MAXIMUM_NUMBER_OF_ITERATIONS)
              MAXIMUM_LENGTH_DIFFERENCE=0.0_DP
              MAXIMUM_DIFFERENCE_LINE=1
              !Loop over the lines
              DO nl=1,FIELD%DECOMPOSITION%TOPOLOGY%LINES%NUMBER_OF_LINES
                CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,nl,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                OLD_LINE_LENGTH=FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(nl)
                LINE_LENGTH=0.0_DP
                !Integrate || dr(xi)/dt || from xi=0 to 1 to determine the arc length.
                DO ng=1,NUMBER_OF_GAUSS_POINTS
                  XI(1)=XIG(GAUSS_START(NUMBER_OF_GAUSS_POINTS)+ng)
                  W=WIG(GAUSS_START(NUMBER_OF_GAUSS_POINTS)+ng)
                  CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
                  CALL COORDINATE_DERIVATIVE_NORM(COORDINATE_SYSTEM,PART_DERIV_S1,INTERPOLATED_POINT,DERIV_NORM, &
                    & ERR,ERROR,*999)
                  LINE_LENGTH=LINE_LENGTH+W*DERIV_NORM
                ENDDO !ng
                FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(nl)=LINE_LENGTH
                LENGTH_DIFFERENCE=ABS(LINE_LENGTH-OLD_LINE_LENGTH)/(1.0_DP+OLD_LINE_LENGTH)
                IF(LENGTH_DIFFERENCE>MAXIMUM_LENGTH_DIFFERENCE) THEN
                  MAXIMUM_LENGTH_DIFFERENCE=LENGTH_DIFFERENCE
                  MAXIMUM_DIFFERENCE_LINE=nl
                ENDIF
              ENDDO !nl
              ITERATE=MAXIMUM_LENGTH_DIFFERENCE>LINE_INCREMENT_TOLERANCE
              IF(ITERATE) THEN
                IF(ITERATION_NUMBER==1) THEN
                  LAST_MAXIMUM_LENGTH_DIFFERENCE=MAXIMUM_LENGTH_DIFFERENCE
                ELSE IF(MAXIMUM_LENGTH_DIFFERENCE<LOOSE_TOLERANCE.AND. &
                  & MAXIMUM_LENGTH_DIFFERENCE>=LAST_MAXIMUM_LENGTH_DIFFERENCE) THEN
                  !Seems to be at a numerical limit
                  ITERATE=.FALSE.
                ELSE
                  LAST_MAXIMUM_LENGTH_DIFFERENCE=MAXIMUM_LENGTH_DIFFERENCE
                ENDIF
              ENDIF
              ITERATION_NUMBER=ITERATION_NUMBER+1
              IF(DIAGNOSTICS2) THEN
                CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Line iteration report:",ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of iterations = ",ITERATION_NUMBER,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Maximum length difference = ",MAXIMUM_LENGTH_DIFFERENCE, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Difference tolerance = ",LINE_INCREMENT_TOLERANCE, &
                  ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Maximum difference line = ",MAXIMUM_DIFFERENCE_LINE, &
                  ERR,ERROR,*999)
              ENDIF
              IF(.NOT.ITERATE.OR.ITERATION_NUMBER==LINES_MAXIMUM_NUMBER_OF_ITERATIONS) THEN
                UPDATE_FIELDS_USING=.TRUE.
              ELSE
                UPDATE_FIELDS_USING=.FALSE.
              ENDIF
              CALL FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE(FIELD,UPDATE_FIELDS_USING,ERR,ERROR,*999)
            ENDDO !iterate
            CALL FIELD_INTERPOLATED_POINT_FINALISE(INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="Geometric parameters are not associated for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is not a geometric field."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Line lengths:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of iterations = ",ITERATION_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Maximum length difference = ",MAXIMUM_LENGTH_DIFFERENCE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Difference tolerance = ",LINE_INCREMENT_TOLERANCE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Maximum difference line = ",MAXIMUM_DIFFERENCE_LINE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",FIELD%DECOMPOSITION%TOPOLOGY%LINES%NUMBER_OF_LINES, &
        & ERR,ERROR,*999)
      DO nl=1,FIELD%DECOMPOSITION%TOPOLOGY%LINES%NUMBER_OF_LINES
        CALL WRITE_STRING_FMT_TWO_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Line ",nl,"(I8)"," length = ",FIELD% &
          & GEOMETRIC_FIELD_PARAMETERS% LENGTHS(nl),"*",ERR,ERROR,*999)
      ENDDO !nl
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE")
    RETURN
999 IF(ASSOCIATED(INTERPOLATED_POINT)) CALL FIELD_INTERPOLATED_POINT_FINALISE(INTERPOLATED_POINT,DUMMY_ERR,DUMMY_ERROR,*999)
    IF(ASSOCIATED(INTERPOLATION_PARAMETERS)) CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(INTERPOLATION_PARAMETERS, &
      & DUMMY_ERR,DUMMY_ERROR,*999)
    CALL ERRORS("FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_LINE_LENGTHS_CALCULATE


  !
  !================================================================================================================================
  !

  !>Finalises the geometric parameters for a field and deallocates all memory.
  SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE(FIELD,UPDATE_FIELDS_USING,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update the scale factors for
    LOGICAL, INTENT(IN) :: UPDATE_FIELDS_USING !<If .TRUE. then update the fields that use this fields geometric parameters.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx,LAST_FIELD_IDX
    TYPE(FIELD_TYPE), POINTER :: FIELD2

    CALL ENTERS("FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
        IF(UPDATE_FIELDS_USING) THEN
          LAST_FIELD_IDX=FIELD%GEOMETRIC_FIELD_PARAMETERS%NUMBER_OF_FIELDS_USING
        ELSE
          LAST_FIELD_IDX=1 !The first field using will be the current field
        ENDIF
        DO field_idx=1,LAST_FIELD_IDX
          FIELD2=>FIELD%GEOMETRIC_FIELD_PARAMETERS%FIELDS_USING(field_idx)%PTR
          CALL FIELD_SCALINGS_CALCULATE(FIELD2,ERR,ERROR,*999)
        ENDDO !field_idx
      ELSE
        CALL FLAG_ERROR("Field is not geometric field.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE")
    RETURN
999 CALL ERRORS("FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE",ERR,ERROR)
    CALL EXITS("FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE")
    RETURN 1
  END SUBROUTINE FIELD_GEOMETRIC_PARAMETERS_SCALE_FACTORS_UPDATE

  !
  !================================================================================================================================
  !

  !>Gets the mesh decomposition for a field.
  SUBROUTINE FIELD_MESH_DECOMPOSITION_GET(FIELD,MESH_DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the decomposition for
    TYPE(DECOMPOSITION_TYPE), POINTER :: MESH_DECOMPOSITION !<On return, a pointer to the mesh decomposition for the field. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("FIELD_MESH_DECOMPOSITION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(ASSOCIATED(MESH_DECOMPOSITION)) THEN
          CALL FLAG_ERROR("Mesh decomposition is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(MESH_DECOMPOSITION)
          MESH_DECOMPOSITION=>FIELD%DECOMPOSITION
          IF(.NOT.ASSOCIATED(MESH_DECOMPOSITION)) CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_MESH_DECOMPOSITION_GET")
    RETURN
999 CALL ERRORS("FIELD_MESH_DECOMPOSITION_GET",ERR,ERROR)
    CALL EXITS("FIELD_MESH_DECOMPOSITION_GET")
    RETURN 1
  END SUBROUTINE FIELD_MESH_DECOMPOSITION_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field.
  SUBROUTINE FIELD_MESH_DECOMPOSITION_SET(FIELD,MESH_DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the decomposition for
    TYPE(DECOMPOSITION_TYPE), POINTER :: MESH_DECOMPOSITION !<A pointer to the mesh decomposition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_MESH_DECOMPOSITION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%DECOMPOSITION_LOCKED) THEN
            LOCAL_ERROR="The mesh decomposition has been locked for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(MESH_DECOMPOSITION)) THEN
              IF(ASSOCIATED(MESH_DECOMPOSITION%MESH)) THEN
                IF(ASSOCIATED(MESH_DECOMPOSITION%MESH%REGION)) THEN
                  IF(ASSOCIATED(FIELD%REGION)) THEN
                    IF(MESH_DECOMPOSITION%MESH%REGION%USER_NUMBER==FIELD%REGION%USER_NUMBER) THEN
                      FIELD%DECOMPOSITION=>MESH_DECOMPOSITION
                    ELSE
                      LOCAL_ERROR="Inconsitent regions. The field is defined on region number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%REGION%USER_NUMBER,"*",ERR,ERROR))// &
                        & " and the mesh decomposition is defined on region number "//&
                        & TRIM(NUMBER_TO_VSTRING(MESH_DECOMPOSITION%MESH%REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Region is not associated for field number "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Region is not associated for the decomposition mesh number "// &
                    & TRIM(NUMBER_TO_VSTRING(MESH_DECOMPOSITION%MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Mesh is not associated for the mesh decomposition.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Mesh decomposition is not assocaited.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_MESH_DECOMPOSITION_SET")
    RETURN
999 CALL ERRORS("FIELD_MESH_DECOMPOSITION_SET",ERR,ERROR)
    CALL EXITS("FIELD_MESH_DECOMPOSITION_SET")
    RETURN 1
  END SUBROUTINE FIELD_MESH_DECOMPOSITION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh decomposition for a field and locks so that no further changes can be made.
  SUBROUTINE FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(FIELD,MESH_DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the decomposition for
    TYPE(DECOMPOSITION_TYPE), POINTER :: MESH_DECOMPOSITION !<A pointer to the mesh decomposition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_MESH_DECOMPOSITION_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_MESH_DECOMPOSITION_SET(FIELD,MESH_DECOMPOSITION,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%DECOMPOSITION_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_MESH_DECOMPOSITION_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_MESH_DECOMPOSITION_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_MESH_DECOMPOSITION_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_MESH_DECOMPOSITION_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Checks the number of field components for a field variable.
  SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_CHECK(FIELD,VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the number of components
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to check \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !The number of components in the field variable to check
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_COMPONENTS_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS/=NUMBER_OF_COMPONENTS) THEN
              LOCAL_ERROR="Invalid number of components. The number components for variable type "// &
                & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                & " which is does correspond to the specified number of components of "// &
                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          ENDIF
        ELSE
          LOCAL_ERROR="The supplied variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The field variable type must be > 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_CHECK")
    RETURN
999 CALL ERRORS("FIELD_NUMBER_OF_COMPONENTS_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the number of field components for a field variable.
  SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_GET(FIELD,VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the number of components
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COMPONENTS !<On return, the number of components in the field variable
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_COMPONENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            NUMBER_OF_COMPONENTS=FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          ENDIF
        ELSE
          LOCAL_ERROR="The supplied variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The field variable type must be > 1 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_GET")
    RETURN
999 CALL ERRORS("FIELD_NUMBER_OF_COMPONENTS_GET",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_GET")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of field components for a field variable.
  SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_SET(FIELD,VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the number of components
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEW_NUMBER_OF_COMPONENTS,OLD_NUMBER_OF_COMPONENTS,variable_idx
    INTEGER(INTG), ALLOCATABLE :: OLD_INTERPOLATION_TYPE(:,:),OLD_MESH_COMPONENT_NUMBER(:,:)
    LOGICAL, ALLOCATABLE :: OLD_INTERPOLATION_TYPE_LOCKED(:,:),OLD_MESH_COMPONENT_NUMBER_LOCKED(:,:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_COMPONENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            IF(ANY(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES==VARIABLE_TYPE)) THEN
              IF(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED(VARIABLE_TYPE)) THEN
                LOCAL_ERROR="The number of components has been locked for variable type "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                SELECT CASE(FIELD%CREATE_VALUES_CACHE%DIMENSION(VARIABLE_TYPE))
                CASE(FIELD_SCALAR_DIMENSION_TYPE)
                  IF(NUMBER_OF_COMPONENTS/=1) THEN
                    LOCAL_ERROR="Scalar fields cannot have "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                CASE(FIELD_VECTOR_DIMENSION_TYPE)
                  IF(NUMBER_OF_COMPONENTS>0) THEN
                    IF(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)/=NUMBER_OF_COMPONENTS) THEN
                      OLD_NUMBER_OF_COMPONENTS=MAXVAL(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS)
                      NEW_NUMBER_OF_COMPONENTS=MAX(OLD_NUMBER_OF_COMPONENTS,NUMBER_OF_COMPONENTS)
                      ALLOCATE(OLD_INTERPOLATION_TYPE(OLD_NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type.",ERR,ERROR,*999)
                      ALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED(OLD_NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type locked.",ERR,ERROR,*999)
                      ALLOCATE(OLD_MESH_COMPONENT_NUMBER(OLD_NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old mesh component number.",ERR,ERROR,*999)
                      ALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED(OLD_NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old mesh component number locked.",ERR,ERROR,*999)
                      OLD_INTERPOLATION_TYPE=FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE
                      OLD_INTERPOLATION_TYPE_LOCKED=FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED
                      OLD_MESH_COMPONENT_NUMBER=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER
                      OLD_MESH_COMPONENT_NUMBER_LOCKED=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED
                      DEALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE)
                      DEALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED)
                      DEALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER)
                      DEALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED)
                      ALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(NEW_NUMBER_OF_COMPONENTS, &
                        & FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type.",ERR,ERROR,*999)
                      ALLOCATE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(NEW_NUMBER_OF_COMPONENTS, &
                        & FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interpolation type locked.",ERR,ERROR,*999)
                      ALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(NEW_NUMBER_OF_COMPONENTS, &
                        & FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh component number.",ERR,ERROR,*999)
                      ALLOCATE(FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(NEW_NUMBER_OF_COMPONENTS, &
                        & FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh component number locked.",ERR,ERROR,*999)
                      FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE=0
                      FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED=.FALSE.
                      FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER=0                   
                      FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED=.FALSE.
                      IF(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)<NUMBER_OF_COMPONENTS) THEN
                        DO variable_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                          FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(1:FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS( &
                            & VARIABLE_TYPE),variable_idx)=OLD_INTERPOLATION_TYPE(1:FIELD%CREATE_VALUES_CACHE% &
                            & NUMBER_OF_COMPONENTS(VARIABLE_TYPE),variable_idx)
                          FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(1:FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS( &
                            & VARIABLE_TYPE),variable_idx)=OLD_INTERPOLATION_TYPE_LOCKED(1:FIELD%CREATE_VALUES_CACHE% &
                            & NUMBER_OF_COMPONENTS(VARIABLE_TYPE),variable_idx)
                          FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS( &
                            & VARIABLE_TYPE)+1:NUMBER_OF_COMPONENTS,variable_idx)=OLD_INTERPOLATION_TYPE(1,variable_idx)
                          FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(1:FIELD%CREATE_VALUES_CACHE% &
                            & NUMBER_OF_COMPONENTS(VARIABLE_TYPE),variable_idx)=OLD_MESH_COMPONENT_NUMBER(1:FIELD% &
                            & CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE),variable_idx)
                          FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(1:FIELD%CREATE_VALUES_CACHE% &
                            & NUMBER_OF_COMPONENTS(VARIABLE_TYPE),variable_idx)=OLD_MESH_COMPONENT_NUMBER_LOCKED(1:FIELD% &
                            & CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE),variable_idx)
                          FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS( &
                            & VARIABLE_TYPE)+1:NUMBER_OF_COMPONENTS,variable_idx)=OLD_MESH_COMPONENT_NUMBER(1,variable_idx)
                        ENDDO !variable_idx
                      ELSE
                        DO variable_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                          FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(1:NUMBER_OF_COMPONENTS,variable_idx)= &
                            & OLD_INTERPOLATION_TYPE(1:NUMBER_OF_COMPONENTS,variable_idx)
                          FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(1:NUMBER_OF_COMPONENTS,variable_idx)= &
                            & OLD_INTERPOLATION_TYPE_LOCKED(1:NUMBER_OF_COMPONENTS,variable_idx)
                          FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(1:NUMBER_OF_COMPONENTS,variable_idx)= &
                            & OLD_MESH_COMPONENT_NUMBER(1:NUMBER_OF_COMPONENTS,variable_idx)
                          FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(1:NUMBER_OF_COMPONENTS,variable_idx)= &
                            & OLD_MESH_COMPONENT_NUMBER_LOCKED(1:NUMBER_OF_COMPONENTS,variable_idx)
                        ENDDO !variable_idx
                      ENDIF
                      FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(VARIABLE_TYPE)=NUMBER_OF_COMPONENTS
                      DEALLOCATE(OLD_INTERPOLATION_TYPE)
                      DEALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED)
                      DEALLOCATE(OLD_MESH_COMPONENT_NUMBER)
                      DEALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Vector fields cannot have "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                CASE(FIELD_TENSOR_DIMENSION_TYPE)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Field dimension "//TRIM(NUMBER_TO_VSTRING(FIELD%CREATE_VALUES_CACHE%DIMENSION( &
                    & VARIABLE_TYPE),"*",ERR,ERROR))//" is not valid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_SET")
    RETURN
999 IF(ALLOCATED(OLD_INTERPOLATION_TYPE)) DEALLOCATE(OLD_INTERPOLATION_TYPE)
    IF(ALLOCATED(OLD_INTERPOLATION_TYPE_LOCKED)) DEALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED)
    IF(ALLOCATED(OLD_MESH_COMPONENT_NUMBER)) DEALLOCATE(OLD_MESH_COMPONENT_NUMBER)
    IF(ALLOCATED(OLD_MESH_COMPONENT_NUMBER_LOCKED)) DEALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED)
    CALL ERRORS("FIELD_NUMBER_OF_COMPONENTS_SET",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_SET")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of field components for a field variable and locks so that no further changes can be made.
  SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(FIELD,VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the number of components
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to be set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_NUMBER_OF_COMPONENTS_SET(FIELD,VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED(VARIABLE_TYPE)=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Checks the number of variables for a field.
  SUBROUTINE FIELD_NUMBER_OF_VARIABLES_CHECK(FIELD,NUMBER_OF_VARIABLES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the number of variables for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_VARIABLES !<The number of variables in the specified field to check
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_VARIABLES_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(FIELD%NUMBER_OF_VARIABLES/=NUMBER_OF_VARIABLES) THEN
          LOCAL_ERROR="Invalid number of variables. The number of variables for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))// &
            & " which is does correspond to the specified number of variables of "// &
            & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_VARIABLES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
       ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_CHECK")
    RETURN
999 CALL ERRORS("FIELD_NUMBER_OF_VARIABLES_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_VARIABLES_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the number of variables for a field.
  SUBROUTINE FIELD_NUMBER_OF_VARIABLES_GET(FIELD,NUMBER_OF_VARIABLES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the number of variables for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_VARIABLES !<On return, the number of variables in the specified field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_VARIABLES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        NUMBER_OF_VARIABLES=FIELD%NUMBER_OF_VARIABLES
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_GET")
    RETURN
999 CALL ERRORS("FIELD_NUMBER_OF_VARIABLES_GET",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_GET")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_VARIABLES_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field.
  SUBROUTINE FIELD_NUMBER_OF_VARIABLES_SET(FIELD,NUMBER_OF_VARIABLES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the number of variables for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_VARIABLES !<The number of variables to set for the field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx,variable_idx2,variable_type
    INTEGER(INTG), ALLOCATABLE :: OLD_VARIABLE_TYPES(:)
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_VARIABLES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_VARIABLES_LOCKED) THEN
            LOCAL_ERROR="The number of variables has been locked field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(NUMBER_OF_VARIABLES>0.AND.NUMBER_OF_VARIABLES<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
              IF(FIELD%NUMBER_OF_VARIABLES/=NUMBER_OF_VARIABLES) THEN
                ALLOCATE(OLD_VARIABLE_TYPES(FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old variable types.",ERR,ERROR,*999)
                OLD_VARIABLE_TYPES=FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES
                DEALLOCATE(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES)
                ALLOCATE(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(NUMBER_OF_VARIABLES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate variable types.",ERR,ERROR,*999)
                FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES=0
                IF(NUMBER_OF_VARIABLES<FIELD%NUMBER_OF_VARIABLES) THEN
                  FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(1:NUMBER_OF_VARIABLES)=OLD_VARIABLE_TYPES(1:NUMBER_OF_VARIABLES)
                  DO variable_idx=NUMBER_OF_VARIABLES+1,FIELD%NUMBER_OF_VARIABLES
                    variable_type=OLD_VARIABLE_TYPES(variable_idx)
                    FIELD%CREATE_VALUES_CACHE%DIMENSION(variable_type)=0
                    FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED(variable_type)=.FALSE.
                    FIELD%CREATE_VALUES_CACHE%DATA_TYPES(variable_type)=0
                    FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED(variable_type)=.FALSE.
                    FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_type)=0
                    FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED(variable_type)=.FALSE.
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(:,variable_type)=0
                    FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(:,variable_type)=.FALSE.
                    FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(:,variable_type)=0
                    FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(:,variable_type)=.FALSE.
                  ENDDO !variable_idx
                ELSE
                  FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(1:FIELD%NUMBER_OF_VARIABLES)= &
                    & OLD_VARIABLE_TYPES(1:FIELD%NUMBER_OF_VARIABLES)
                  DO variable_idx=FIELD%NUMBER_OF_VARIABLES+1,NUMBER_OF_VARIABLES
                    !Find the next available variable type
                    DO variable_type=1,FIELD_NUMBER_OF_VARIABLE_TYPES
                      FOUND=.FALSE.
                      DO variable_idx2=1,FIELD%NUMBER_OF_VARIABLES
                        IF(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(variable_idx2)==variable_type) THEN
                          FOUND=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO !variable_idx2
                      IF(.NOT.FOUND) EXIT
                    ENDDO !variable_type
                    IF(FOUND) THEN
                      CALL FLAG_ERROR("Could not find free variable type???",ERR,ERROR,*999)
                    ELSE
                      FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(variable_idx)=variable_type
                      FIELD%CREATE_VALUES_CACHE%DIMENSION(variable_type)=FIELD%CREATE_VALUES_CACHE%DIMENSION( &
                        & OLD_VARIABLE_TYPES(1))
                      FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED(variable_type)=.FALSE.
                      FIELD%CREATE_VALUES_CACHE%DATA_TYPES(variable_type)=FIELD%CREATE_VALUES_CACHE%DATA_TYPES( &
                        & OLD_VARIABLE_TYPES(1))
                      FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED(variable_type)=.FALSE.
                      FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_type)=FIELD%CREATE_VALUES_CACHE% &
                        & NUMBER_OF_COMPONENTS(OLD_VARIABLE_TYPES(1))
                      FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED(variable_type)=.FALSE.
                      FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(:,variable_type)=FIELD%CREATE_VALUES_CACHE% &
                        INTERPOLATION_TYPE(:,FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(OLD_VARIABLE_TYPES(1)))
                      FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(:,variable_type)=.FALSE.
                      FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(:,variable_type)=FIELD%CREATE_VALUES_CACHE% &
                        MESH_COMPONENT_NUMBER(:,FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(OLD_VARIABLE_TYPES(1)))
                      FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(:,variable_type)=.FALSE.
                    ENDIF
                  ENDDO !variable_idx
                ENDIF
                DEALLOCATE(OLD_VARIABLE_TYPES)
                FIELD%NUMBER_OF_VARIABLES=NUMBER_OF_VARIABLES
              ENDIF
            ELSE
              LOCAL_ERROR="The specified number of variables of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_VARIABLES,"*",ERR,ERROR))// &
                & " is invalid. The number of variables must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_SET")
    RETURN
999 IF(ALLOCATED(OLD_VARIABLE_TYPES)) DEALLOCATE(OLD_VARIABLE_TYPES)
    CALL ERRORS("FIELD_NUMBER_OF_VARIABLES_SET",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_SET")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_VARIABLES_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of variables for a field and locks so that no further changes can be made.
  SUBROUTINE FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(FIELD,NUMBER_OF_VARIABLES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the number of variables for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_VARIABLES !<The number of variables to set for the field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_NUMBER_OF_VARIABLES_SET(FIELD,NUMBER_OF_VARIABLES,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%NUMBER_OF_VARIABLES_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF    
 
    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Adds the alpha times the parameter set values from one parameter set type to another parameter set type \todo make this call distributed vector add???
  SUBROUTINE FIELD_PARAMETER_SETS_ADD_DP(FIELD,VARIABLE_TYPE,ALPHA,FIELD_FROM_SET_TYPE,FIELD_TO_SET_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add the parameter sets for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    REAL(DP), INTENT(IN) :: ALPHA(:) !<The multiplicative factor for the add.
    INTEGER(INTG), INTENT(IN) :: FIELD_FROM_SET_TYPE(:) !<The field parameter set identifier to add the parameters from \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_TO_SET_TYPE !<The field parameter set identifier to add the parameters to \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dof_idx,parameter_set_idx
    REAL(DP) :: VALUE
    TYPE(REAL_DP_PTR_TYPE) :: FIELD_FROM_PARAMETERS(SIZE(FIELD_FROM_SET_TYPE,1))
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: FIELD_FROM_PARAMETER_SET,FIELD_TO_PARAMETER_SET
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE) :: FIELD_FROM_PARAMETER_SETS(SIZE(FIELD_FROM_SET_TYPE,1))
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SETS_ADD_DP",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            !Check the to set type input
            IF(FIELD_TO_SET_TYPE>0.AND.FIELD_TO_SET_TYPE<FIELD_NUMBER_OF_SET_TYPES) THEN
              FIELD_TO_PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_TO_SET_TYPE)%PTR
              IF(ASSOCIATED(FIELD_TO_PARAMETER_SET)) THEN
                IF(SIZE(ALPHA,1)==SIZE(FIELD_FROM_SET_TYPE,1)) THEN
                  DO parameter_set_idx=1,SIZE(FIELD_FROM_SET_TYPE,1)
                    IF(FIELD_FROM_SET_TYPE(parameter_set_idx)>0.AND. &
                      & FIELD_FROM_SET_TYPE(parameter_set_idx)<FIELD_NUMBER_OF_SET_TYPES) THEN
                      FIELD_FROM_PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_TO_SET_TYPE)%PTR
                      IF(ASSOCIATED(FIELD_TO_PARAMETER_SET)) THEN
                        FIELD_FROM_PARAMETER_SETS(parameter_set_idx)%PTR=>FIELD_FROM_PARAMETER_SET
                        NULLIFY(FIELD_FROM_PARAMETERS(parameter_set_idx)%PTR)
                        CALL DISTRIBUTED_VECTOR_DATA_GET(FIELD_FROM_PARAMETER_SETS(parameter_set_idx)%PTR%PARAMETERS, &
                          & FIELD_FROM_PARAMETERS(parameter_set_idx)%PTR,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="The field from set type of "// &
                          & TRIM(NUMBER_TO_VSTRING(FIELD_FROM_SET_TYPE(parameter_set_idx),"*",ERR,ERROR))// &
                          & " in parameter set index "//TRIM(NUMBER_TO_VSTRING(parameter_set_idx,"*",ERR,ERROR))// &
                          & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & "."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="The field from set type of "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_FROM_SET_TYPE(parameter_set_idx),"*",ERR,ERROR))// &
                        & " for parameter set index "//TRIM(NUMBER_TO_VSTRING(parameter_set_idx,"*",ERR,ERROR))// &
                        & " is invalid. The field set TYPE must be between 1 and "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !parameter_set_idx
                  !Do not need to do an update here as each rank already has the values.
                  !Add the field dofs
                  DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                    VALUE=0.0_DP
                    DO parameter_set_idx=1,SIZE(FIELD_FROM_SET_TYPE,1)
                      VALUE=VALUE+ALPHA(parameter_set_idx)*FIELD_FROM_PARAMETERS(parameter_set_idx)%PTR(dof_idx)
                    ENDDO !parameter_set_idx
                    CALL DISTRIBUTED_VECTOR_VALUES_ADD(FIELD_TO_PARAMETER_SET%PARAMETERS,dof_idx,VALUE,ERR,ERROR,*999)
                  ENDDO !dof_idx
                  !Restore the from parameter set transfer
                  DO parameter_set_idx=1,SIZE(FIELD_FROM_SET_TYPE,1)
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(FIELD_FROM_PARAMETER_SETS(parameter_set_idx)%PTR%PARAMETERS, &
                      & FIELD_FROM_PARAMETERS(parameter_set_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !parameter_set_idx
                ELSE
                  LOCAL_ERROR="The size of the alpha array ("//TRIM(NUMBER_TO_VSTRING(SIZE(ALPHA,1),"*",ERR,ERROR))// &
                    & ") does not match the size of the from set type array ("// &
                    & TRIM(NUMBER_TO_VSTRING(SIZE(FIELD_FROM_SET_TYPE,1),"*",ERR,ERROR))//"."
                ENDIF
              ELSE
                LOCAL_ERROR="The field to set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_TO_SET_TYPE,"*",ERR,ERROR))// &
                  & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field to set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_TO_SET_TYPE,"*",ERR,ERROR))// &
                & " is invalid. The field set type must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("FIELD_PARAMETER_SETS_ADD_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SETS_ADD_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SETS_ADD_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SETS_ADD_DP
  
 !
  !================================================================================================================================
  !

  !>Adds the alpha times the parameter set values from one parameter set type to another parameter set type \todo make this call distributed vector add???
  SUBROUTINE FIELD_PARAMETER_SETS_ADD_DP1(FIELD,VARIABLE_TYPE,ALPHA,FIELD_FROM_SET_TYPE,FIELD_TO_SET_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add the parameter sets for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINE
    REAL(DP), INTENT(IN) :: ALPHA !<The multiplicative factor for the add.
    INTEGER(INTG), INTENT(IN) :: FIELD_FROM_SET_TYPE !<The field parameter set identifier to add the parameters from
    INTEGER(INTG), INTENT(IN) :: FIELD_TO_SET_TYPE !<The field parameter set identifier to add the parameters to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_PARAMETER_SETS_ADD_DP1",ERR,ERROR,*999)
    
    CALL FIELD_PARAMETER_SETS_ADD_DP(FIELD,VARIABLE_TYPE,(/ALPHA/),(/FIELD_FROM_SET_TYPE/),FIELD_TO_SET_TYPE,ERR,ERROR,*999)
    
    CALL EXITS("FIELD_PARAMETER_SETS_ADD_DP1")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SETS_ADD_DP1",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SETS_ADD_DP1")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SETS_ADD_DP1
  
  !
  !================================================================================================================================
  !

  !>Copys the parameter set from one parameter set type to another parameter set type \todo make this call distributed vector copy???
  SUBROUTINE FIELD_PARAMETER_SETS_COPY(FIELD,VARIABLE_TYPE,FIELD_FROM_SET_TYPE,FIELD_TO_SET_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to copy the parameters set for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to copy \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_FROM_SET_TYPE !<The field parameter set identifier to copy the parameters from \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_TO_SET_TYPE !<The field parameter set identifier to copy the parameters to \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dof_idx
    REAL(DP) :: VALUE
    REAL(DP), POINTER :: FIELD_FROM_PARAMETERS(:)
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: FIELD_FROM_PARAMETER_SET,FIELD_TO_PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SETS_COPY",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            !Check the from set type input
            IF(FIELD_FROM_SET_TYPE>0.AND.FIELD_FROM_SET_TYPE<FIELD_NUMBER_OF_SET_TYPES) THEN
              FIELD_FROM_PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_FROM_SET_TYPE)%PTR
              IF(ASSOCIATED(FIELD_FROM_PARAMETER_SET)) THEN
                !Check the from set type input
                IF(FIELD_TO_SET_TYPE>0.AND.FIELD_TO_SET_TYPE<FIELD_NUMBER_OF_SET_TYPES) THEN
                  FIELD_TO_PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_TO_SET_TYPE)%PTR
                  !Do not need to do an update here as each rank already has the values.
                  IF(ASSOCIATED(FIELD_TO_PARAMETER_SET)) THEN
                    !Get the from parameter set data
                    CALL DISTRIBUTED_VECTOR_DATA_GET(FIELD_FROM_PARAMETER_SET%PARAMETERS,FIELD_FROM_PARAMETERS,ERR,ERROR,*999)
                    !Do not need to do an update here as each rank already has the values.
                    !Loop over the locals
                    DO dof_idx=1,FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS
                      VALUE=FIELD_FROM_PARAMETERS(dof_idx)
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(FIELD_TO_PARAMETER_SET%PARAMETERS,dof_idx,VALUE,ERR,ERROR,*999)
                    ENDDO !dof_idx
                    !Restore the from parameter set data
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(FIELD_FROM_PARAMETER_SET%PARAMETERS,FIELD_FROM_PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field to set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_TO_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field to set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_TO_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field from set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_FROM_SET_TYPE,"*",ERR,ERROR))// &
                  & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field from set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_FROM_SET_TYPE,"*",ERR,ERROR))// &
                & " is invalid. The field set type must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SETS_COPY")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SETS_COPY",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SETS_COPY")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SETS_COPY

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add to
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_CONSTANT_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                    & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_CONSTANT_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_INTG

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add to
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_CONSTANT_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                    & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_CONSTANT_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_SP

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add to
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_CONSTANT_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                    & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_CONSTANT_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_DP

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add to
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_CONSTANT_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                    & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                    & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_CONSTANT_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_CONSTANT_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_CONSTANT_L

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_INTG

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_SP

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_DP

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_DOF_L

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_ELEMENT_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_ELEMENT_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_INTG

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_ELEMENT_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_ELEMENT_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_SP

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_ELEMENT_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_ELEMENT_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_DP

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_ELEMENT_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_ELEMENT_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_ELEMENT_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_ELEMENT_L

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_INTG

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_SP

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_DP

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_ELEMENT_L

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_NODE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_NODE_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_INTG

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_NODE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_NODE_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_SP

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_NODE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_NODE_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_DP

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_NODE_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_NODE_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_NODE_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_NODE_L

  !
  !================================================================================================================================
  !

  !>Adds the given integer value to the given parameter set for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_INTG

  !
  !================================================================================================================================
  !

  !>Adds the given single precision value to the given parameter set for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    REAL(SP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_SP

  !
  !================================================================================================================================
  !

  !>Adds the given double precision value to the given parameter set for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    REAL(DP), INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_DP

  !
  !================================================================================================================================
  !

  !>Adds the given logical value to the given parameter set for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to add
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to add
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to add
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to add
    LOGICAL, INTENT(IN) :: VALUE !<The value to add
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)                  
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_ADD(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not add element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS( &
                        & COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))//" is invalid for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)              
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_ADD_LOCAL_NODE_L

  !
  !================================================================================================================================
  !

  !>Creates a new parameter set of type set type for a field variable.
  SUBROUTINE FIELD_PARAMETER_SET_CREATE(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to create the parameter set for
    INTEGER(INTG),  INTENT(IN) :: VARIABLE_TYPE !<The variable type to create the parameter set for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,parameter_set_idx
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: NEW_PARAMETER_SET
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE), POINTER :: NEW_PARAMETER_SETS(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR,DUMMY_ERROR

    NULLIFY(NEW_PARAMETER_SET)
    NULLIFY(NEW_PARAMETER_SETS)

    CALL ENTERS("FIELD_PARAMETER_SET_CREATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          !Check the set type input
          IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<FIELD_NUMBER_OF_SET_TYPES) THEN
            !Check if this set type has already been created
            IF(ASSOCIATED(FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR)) THEN
              LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                & " has already been created for variable type of "// &
                & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              ALLOCATE(NEW_PARAMETER_SET,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new parameter set.",ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_INITIALISE(NEW_PARAMETER_SET,ERR,ERROR,*999)
              NEW_PARAMETER_SET%SET_INDEX=FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS+1
              NEW_PARAMETER_SET%SET_TYPE=FIELD_SET_TYPE
              NULLIFY(NEW_PARAMETER_SET%PARAMETERS)
              CALL DISTRIBUTED_VECTOR_CREATE_START(FIELD_VARIABLE%DOMAIN_MAPPING,NEW_PARAMETER_SET%PARAMETERS,ERR,ERROR,*999)
              SELECT CASE(FIELD_VARIABLE%DATA_TYPE)
              CASE(FIELD_INTG_TYPE)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NEW_PARAMETER_SET%PARAMETERS,DISTRIBUTED_MATRIX_VECTOR_INTG_TYPE, &
                  & ERR,ERROR,*999)
              CASE(FIELD_SP_TYPE)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NEW_PARAMETER_SET%PARAMETERS,DISTRIBUTED_MATRIX_VECTOR_SP_TYPE, &
                  & ERR,ERROR,*999)
              CASE(FIELD_DP_TYPE)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NEW_PARAMETER_SET%PARAMETERS,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
                  & ERR,ERROR,*999)
              CASE(FIELD_L_TYPE)
                CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(NEW_PARAMETER_SET%PARAMETERS,DISTRIBUTED_MATRIX_VECTOR_L_TYPE, &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The field data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " is invalid for variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(NEW_PARAMETER_SET%PARAMETERS,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(NEW_PARAMETER_SET%PARAMETERS,0.0_DP,ERR,ERROR,*999)
              !Add the new parameter set to the list of parameter sets
              ALLOCATE(NEW_PARAMETER_SETS(FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new parameter sets.",ERR,ERROR,*999)
              IF(ASSOCIATED(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS)) THEN
                DO parameter_set_idx=1,FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS
                  NEW_PARAMETER_SETS(parameter_set_idx)%PTR=>FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(parameter_set_idx)%PTR
                ENDDO !parameter_set_idx
                DEALLOCATE(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS)
              ENDIF
              NEW_PARAMETER_SETS(FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS+1)%PTR=>NEW_PARAMETER_SET
              ALLOCATE(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS+1), &
                & STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field parameter sets parameter sets.",ERR,ERROR,*999)
              DO parameter_set_idx=1,FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS+1
                FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(parameter_set_idx)%PTR=>NEW_PARAMETER_SETS(parameter_set_idx)%PTR
              ENDDO !parameter_set_idx
              DEALLOCATE(NEW_PARAMETER_SETS)
              FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS=FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS+1
              FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR=>NEW_PARAMETER_SET
            ENDIF
          ELSE
            LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The field set type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
          & " is invalid. The variable type must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_CREATE")
    RETURN
999 IF(ASSOCIATED(NEW_PARAMETER_SET)) THEN
      CALL FIELD_PARAMETER_SET_FINALISE(NEW_PARAMETER_SET,DUMMY_ERR,DUMMY_ERROR,*998)
998   DEALLOCATE(NEW_PARAMETER_SET)
    ENDIF
    IF(ASSOCIATED(NEW_PARAMETER_SETS)) DEALLOCATE(NEW_PARAMETER_SETS)
    CALL ERRORS("FIELD_PARAMETER_SET_CREATE",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_CREATE")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_CREATE

  !
  !================================================================================================================================
  !

  !>Destroys the parameter set of type set type for a field variable and deallocates all memory.
  SUBROUTINE FIELD_PARAMETER_SET_DESTROY(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to destroy a parameter set for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to destroy the parameter set for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: parameter_set_idx,SET_INDEX
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_PARAMETER_SET_PTR_TYPE), POINTER :: NEW_PARAMETER_SETS(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(NEW_PARAMETER_SETS)

    CALL ENTERS("FIELD_PARAMETER_SET_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          !Check the set type input
          IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<FIELD_NUMBER_OF_SET_TYPES) THEN
            !Check if the set type has been created
            PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
            IF(ASSOCIATED(PARAMETER_SET)) THEN
              SET_INDEX=PARAMETER_SET%SET_INDEX
              ALLOCATE(NEW_PARAMETER_SETS(FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new parameter sets",ERR,ERROR,*999)
              DO parameter_set_idx=1,FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS
                IF(parameter_set_idx<SET_INDEX) THEN
                  NEW_PARAMETER_SETS(parameter_set_idx)%PTR=>FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(parameter_set_idx)%PTR
                ELSE IF(parameter_set_idx>SET_INDEX) THEN
                  NEW_PARAMETER_SETS(parameter_set_idx-1)%PTR=>FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(parameter_set_idx)%PTR
                  NEW_PARAMETER_SETS(parameter_set_idx-1)%PTR%SET_INDEX=NEW_PARAMETER_SETS(parameter_set_idx-1)%PTR%SET_INDEX-1
                ENDIF
              ENDDO !parameter_set_idx
              DEALLOCATE(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS)
              FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS=>NEW_PARAMETER_SETS
              FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS=FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS-1
              NULLIFY(FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR)
              CALL FIELD_PARAMETER_SET_FINALISE(PARAMETER_SET,ERR,ERROR,*999)
              DEALLOCATE(PARAMETER_SET)
            ELSE
              LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                & " has not been created for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The field set type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
          & " is invalid. The variable type must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_PARAMETER_SETS)) DEALLOCATE(NEW_PARAMETER_SETS)
    CALL ERRORS("FIELD_PARAMETER_SET_DESTROY",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DESTROY")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the parameter set for a field and deallocates all memory.
  SUBROUTINE FIELD_PARAMETER_SET_FINALISE(FIELD_PARAMETER_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: FIELD_PARAMETER_SET !<A pointer to the field parameter set to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_PARAMETER_SET_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD_PARAMETER_SET)) THEN
      IF(ASSOCIATED(FIELD_PARAMETER_SET%PARAMETERS)) CALL DISTRIBUTED_VECTOR_DESTROY(FIELD_PARAMETER_SET%PARAMETERS,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field integer parameter set array. The pointer must be restored with a call to FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_RESTORE call. Note: the values can be used for read operations but a FIELD_ROUTINES::FIELD_PARAMETER_SET_UPDATE call must be used to change any values.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), POINTER :: PARAMETERS(:) !<On return, a pointer to the field parameter set data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_GET_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(PARAMETERS)) THEN
        CALL FLAG_ERROR("Parameters is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(PARAMETERS)
        IF(FIELD%FIELD_FINISHED) THEN
          IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " does not correspond to the integer data type of the given parameters array."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been finished."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_GET_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_INTG

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field single precision parameter set array. The pointer must be restored with a call to FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_RESTORE call. Note: the values can be used for read operations but a FIELD_ROUTINES::FIELD_PARAMETER_SET_UPDATE call must be used to change any values.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    REAL(SP), POINTER :: PARAMETERS(:) !<On return, a pointer to the field parameter set data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_GET_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(PARAMETERS)) THEN
        CALL FLAG_ERROR("Parameters is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(PARAMETERS)
        IF(FIELD%FIELD_FINISHED) THEN
          IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " does not correspond to the single precision data type of the given parameters array."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been finished."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_GET_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_SP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field double precision parameter set array. The pointer must be restored with a call to FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_RESTORE call. Note: the values can be used for read operations but a FIELD_ROUTINES::FIELD_PARAMETER_SET_UPDATE call must be used to change any values.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    REAL(DP), POINTER :: PARAMETERS(:) !<On return, a pointer to the field parameter set data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_GET_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(PARAMETERS)) THEN
        CALL FLAG_ERROR("Parameters is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(PARAMETERS)
        IF(FIELD%FIELD_FINISHED) THEN
          IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " does not correspond to the double precision data type of the given parameters array."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been finished."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_GET_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_DP

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field logical parameter set array. The pointer must be restored with a call to FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_RESTORE call. Note: the values can be used for read operations but a FIELD_ROUTINES::FIELD_PARAMETER_SET_UPDATE call must be used to change any values.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to get the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    LOGICAL, POINTER :: PARAMETERS(:) !<On return, a pointer to the field parameter set data
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_GET_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(PARAMETERS)) THEN
        CALL FLAG_ERROR("Parameters is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(PARAMETERS)
        IF(FIELD%FIELD_FINISHED) THEN
          IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                  & " does not correspond to the logical data type of the given parameters array."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been finished."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_GET_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_GET_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_GET_L

  !
  !================================================================================================================================
  !

  !>Restores the specified field variable integer parameter set array that was obtained with FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_GET.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to restore the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field varaible type to restore the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), POINTER :: PARAMETERS(:) !<The pointer to the field parameter set data obtained with the parameter set get call
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_RESTORE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(ASSOCIATED(PARAMETERS)) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created on variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Parameters is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given parameters array."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_RESTORE_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_INTG

   !
  !================================================================================================================================
  !

  !>Restores the specified field variable single precision parameter set array that was obtained with FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_GET.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to restore the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field varaible type to restore the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    REAL(SP), POINTER :: PARAMETERS(:) !<The pointer to the field parameter set data obtained with the parameter set get call
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_RESTORE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(ASSOCIATED(PARAMETERS)) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created on variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Parameters is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given parameters array."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_RESTORE_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_SP

   !
  !================================================================================================================================
  !

  !>Restores the specified field variable double precision parameter set array that was obtained with FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_GET.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to restore the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field varaible type to restore the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    REAL(DP), POINTER :: PARAMETERS(:) !<The pointer to the field parameter set data obtained with the parameter set get call
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_RESTORE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(ASSOCIATED(PARAMETERS)) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created on variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Parameters is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given parameters array."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_RESTORE_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_DP

   !
  !================================================================================================================================
  !

  !>Restores the specified field variable logical parameter set array that was obtained with FIELD_ROUTINES::FIELD_PARAMETER_SET_DATA_GET.
  SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,PARAMETERS,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to restore the parameter set from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field varaible type to restore the parameter set data for \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    LOGICAL, POINTER :: PARAMETERS(:) !<The pointer to the field parameter set data obtained with the parameter set get call
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_DATA_RESTORE_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(ASSOCIATED(PARAMETERS)) THEN
                IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                  PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                  IF(ASSOCIATED(PARAMETER_SET)) THEN
                    CALL DISTRIBUTED_VECTOR_DATA_RESTORE(PARAMETER_SET%PARAMETERS,PARAMETERS,ERR,ERROR,*999)
                  ELSE
                    LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                      & " has not been created on variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field parameter set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " is invalid. The field parameter set type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Parameters is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given parameters array."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_DATA_RESTORE_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_DATA_RESTORE_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_DATA_RESTORE_L

  !
  !================================================================================================================================
  !

  !>Initialises the parameter set for a field.
  SUBROUTINE FIELD_PARAMETER_SET_INITIALISE(FIELD_PARAMETER_SET,ERR,ERROR,*)

   !Argument variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: FIELD_PARAMETER_SET !<The field parameter set to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_PARAMETER_SET_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD_PARAMETER_SET)) THEN
      FIELD_PARAMETER_SET%SET_INDEX=0
      FIELD_PARAMETER_SET%SET_TYPE=0
    ELSE
      CALL FLAG_ERROR("Field parameter set is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_INITIALISE

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    INTEGER, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_INTG

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    REAL(SP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_SP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    REAL(DP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_DP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for the constant of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,COMPONENT_NUMBER,VALUE,ERR,ERROR,*)
    
    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    LOGICAL, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      IF(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NUMBER_OF_CONSTANT_PARAMETERS>0) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " does not have any constant parameters."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by constant for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_CONSTANT_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_CONSTANT_L

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to update
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
!!TODO: Allow to specify a global number and then have it all update accordingly???
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_INTG

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to update
    REAL(SP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
!!TODO: Allow to specify a global number and then have it all update accordingly???
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_SP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to update
    REAL(DP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
!!TODO: Allow to specify a global number and then have it all update accordingly???
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_DP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for a particular local dof of the field variable.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DOF_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: DOF_NUMBER !<The dof number to update
    LOGICAL, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_DOF_NUMBER
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L",ERR,ERROR,*999)

!!TODO: Allow multiple dof number and values updates.
    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
!!TODO: Allow to specify a global number and then have it all update accordingly???
                  !Note that dofs are slightly different from other mappings in that all the local dofs are not all at the start.
                  !This is because the dof indicies are from combined field components. Thus need to check that a ghost value is
                  !not being set.
                  IF(DOF_NUMBER>0.AND.DOF_NUMBER<=FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    GLOBAL_DOF_NUMBER=FIELD_VARIABLE%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(DOF_NUMBER)
                    IF(FIELD_VARIABLE%DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_DOF_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,DOF_NUMBER,VALUE,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                        & " is invalid as it is a ghost dof for this domain."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The field dof number of "//TRIM(NUMBER_TO_VSTRING(DOF_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. It must be >0 and <="// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))// &
                      & " for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF_L

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_INTG

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    REAL(SP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_SP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    REAL(DP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_DP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for a particular user element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,USER_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    LOGICAL, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_ELEMENT_NUMBER,MESH_COMPONENT,MESH_GLOBAL_ELEMENT_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_ELEMENT_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_ELEMENT_NUMBER,MESH_ELEMENT_EXISTS, &
                            & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                          IF(MESH_ELEMENT_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                                IF(ASSOCIATED(ELEMENTS_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(ELEMENTS_MAPPING,MESH_GLOBAL_ELEMENT_NUMBER, &
                                    & LOCAL_EXISTS,DOMAIN_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN                                  
                                    ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                      & ELEMENT_PARAM2DOF_MAP(DOMAIN_LOCAL_ELEMENT_NUMBER)
                                    CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                  ELSE
                                    LOCAL_ERROR="The specified user element number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings elements mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user element number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_ELEMENT_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_ELEMENT_L

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_INTG

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    REAL(SP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_SP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    REAL(DP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_DP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for a particular local element of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,LOCAL_ELEMENT_NUMBER,COMPONENT_NUMBER, &
    & VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    INTEGER(INTG), INTENT(IN) :: LOCAL_ELEMENT_NUMBER !<The local element number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component to update
    LOGICAL, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      IF(LOCAL_ELEMENT_NUMBER>0.AND.LOCAL_ELEMENT_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                        & PARAM_TO_DOF_MAP%NUMBER_OF_ELEMENT_PARAMETERS) THEN
                        ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP( &
                          & LOCAL_ELEMENT_NUMBER)                      
                        CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                      ELSE
                        LOCAL_ERROR="Local element number "//TRIM(NUMBER_TO_VSTRING(LOCAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" elements."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has node based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by element for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT_L

  !
  !================================================================================================================================
  !

  !>Finishes the the parameter set update for a field.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_FINISH(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,ERR,ERROR,*)

     !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finish the update for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
              PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
              IF(ASSOCIATED(PARAMETER_SET)) THEN
                CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(PARAMETER_SET%PARAMETERS,ERR,ERROR,*999)
                IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE.AND.FIELD_SET_TYPE==FIELD_VALUES_SET_TYPE) THEN
                  !Geometric field values have changed so update the geometric parmeters (e.g., lines etc.)
                  CALL FIELD_GEOMETRIC_PARAMETERS_CALCULATE(FIELD,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                & " is invalid. The field set type must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_FINISH")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_FINISH",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_FINISH")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_FINISH

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_NODE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_NODE_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_INTG

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    REAL(SP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_NODE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_NODE_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_SP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    REAL(DP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_NODE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_NODE_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_DP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for a particular user node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,USER_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    LOGICAL, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DOMAIN_LOCAL_NODE_NUMBER,MESH_COMPONENT,MESH_GLOBAL_NODE_NUMBER,ny
    LOGICAL :: LOCAL_EXISTS,MESH_NODE_EXISTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_NODE_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      DECOMPOSITION=>FIELD%DECOMPOSITION
                      IF(ASSOCIATED(DECOMPOSITION)) THEN
                        MESH=>DECOMPOSITION%MESH
                        IF(ASSOCIATED(MESH)) THEN
                          MESH_COMPONENT=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%MESH_COMPONENT_NUMBER
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT,USER_NODE_NUMBER,MESH_NODE_EXISTS, &
                            & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                          IF(MESH_NODE_EXISTS) THEN
                            DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT)%PTR
                            IF(ASSOCIATED(DOMAIN)) THEN
                              DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                              IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                                NODES_MAPPING=>DOMAIN_MAPPINGS%NODES
                                IF(ASSOCIATED(NODES_MAPPING)) THEN
                                  CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,MESH_GLOBAL_NODE_NUMBER,LOCAL_EXISTS, &
                                    & DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
                                  IF(LOCAL_EXISTS) THEN
                                    IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                      & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                                      ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,DOMAIN_LOCAL_NODE_NUMBER)
                                      CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                                    ELSE
                                      LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                                        & " is invalid for user node number "// &
                                        & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))//" of component number "// &
                                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has a maximum of "// &
                                        & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                                        & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    LOCAL_ERROR="The specified user node number of "// &
                                      & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                                      & " is not part of this local domain."                                  
                                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("Domain mappings nodes mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain is not associated for this mesh component.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="The specified user node number of "// &
                              & TRIM(NUMBER_TO_VSTRING(USER_NODE_NUMBER,"*",ERR,ERROR))// &
                              " does not exist in mesh component number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT,"*",ERR,ERROR))//" of mesh number "// &
                              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)            
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Field decomposition is not associated.",ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_NODE_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_NODE_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_NODE_L

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given integer value for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    INTEGER(INTG), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_INTG_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the integer data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_INTG

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given single precision value for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    REAL(SP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_SP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the single precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_SP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given double precision value for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    REAL(DP), INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_DP_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the double precision data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_DP

  !
  !================================================================================================================================
  !

  !>Updates the given parameter set with the given logical value for a particular local node and derivative of the field variable component.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER, &
    & COMPONENT_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to update
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to finish the update for \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The node derivative number to update
    INTEGER(INTG), INTENT(IN) :: LOCAL_NODE_NUMBER !<The local node number to update
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The field variable component number to update
    LOGICAL, INTENT(IN) :: VALUE !<The value to update to
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ny
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            IF(FIELD_VARIABLE%DATA_TYPE==FIELD_L_TYPE) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  IF(COMPONENT_NUMBER>=1.AND.COMPONENT_NUMBER<=FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                    SELECT CASE(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE)
                    CASE(FIELD_CONSTANT_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has constant interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has element based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_NODE_BASED_INTERPOLATION)
                      IF(LOCAL_NODE_NUMBER>0.AND.LOCAL_NODE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                        & NUMBER_OF_NODE_PARAMETERS) THEN
                        IF(DERIVATIVE_NUMBER>0.AND.DERIVATIVE_NUMBER<=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES) THEN
                          ny=FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP( &
                            & DERIVATIVE_NUMBER,LOCAL_NODE_NUMBER)
                          CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETER_SET%PARAMETERS,ny,VALUE,ERR,ERROR,*999)
                        ELSE
                          LOCAL_ERROR="Derivative number "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*",ERR,ERROR))// &
                            & " is invalid for local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                            & " of component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                            & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                            & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                            & " which has a maximum of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                            & PARAM_TO_DOF_MAP%MAX_NUMBER_OF_DERIVATIVES,"*",ERR,ERROR))//" derivatives."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Local node number "//TRIM(NUMBER_TO_VSTRING(LOCAL_NODE_NUMBER,"*",ERR,ERROR))// &
                          & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                          & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                          & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(COMPONENT_NUMBER)% &
                          & PARAM_TO_DOF_MAP%NUMBER_OF_NODE_PARAMETERS,"*",ERR,ERROR))//" nodes."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has grid point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                      LOCAL_ERROR="Can not update by node for component number "// &
                        & TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))//" of variable type "// &
                        & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))//" of field number "// &
                        & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has Gauss point based interpolation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The field component interpolation type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE% &
                        & COMPONENTS(COMPONENT_NUMBER)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
                        & " is invalid for component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                        & " of variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                        & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="Component number "//TRIM(NUMBER_TO_VSTRING(COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid for variable type "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                      & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
                      & TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
                      & " components."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable data type of "//TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%DATA_TYPE,"*",ERR,ERROR))// &
                & " does not correspond to the logical data type of the given value."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE_L

  !
  !================================================================================================================================
  !

  !>Starts the the parameter set update for a field variable.
  SUBROUTINE FIELD_PARAMETER_SET_UPDATE_START(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to start the update for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier to update
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_UPDATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
        FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
          IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
            PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
            IF(ASSOCIATED(PARAMETER_SET)) THEN
              CALL DISTRIBUTED_VECTOR_UPDATE_START(PARAMETER_SET%PARAMETERS,ERR,ERROR,*999)
            ELSE
              LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The field set type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
          & " is invalid. The variable type must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_START")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_UPDATE_START",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_UPDATE_START")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_UPDATE_START

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the specified field parameter set distributed vector. 
  SUBROUTINE FIELD_PARAMETER_SET_VECTOR_GET(FIELD,VARIABLE_TYPE,FIELD_SET_TYPE,DISTRIBUTED_VECTOR,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the parameter set vector from
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to update \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES 
    INTEGER(INTG), INTENT(IN) :: FIELD_SET_TYPE !<The field parameter set identifier
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: DISTRIBUTED_VECTOR !<On return, a pointer to the field parameter set distributed vector. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: PARAMETER_SET
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_PARAMETER_SET_VECTOR_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(DISTRIBUTED_VECTOR)) THEN
        CALL FLAG_ERROR("Distributed vector is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(DISTRIBUTED_VECTOR)
        IF(FIELD%FIELD_FINISHED) THEN
          IF(VARIABLE_TYPE>0.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
              IF(FIELD_SET_TYPE>0.AND.FIELD_SET_TYPE<=FIELD_NUMBER_OF_SET_TYPES) THEN
                PARAMETER_SET=>FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_SET_TYPE)%PTR
                IF(ASSOCIATED(PARAMETER_SET)) THEN
                  DISTRIBUTED_VECTOR=>PARAMETER_SET%PARAMETERS
                  IF(.NOT.ASSOCIATED(DISTRIBUTED_VECTOR)) &
                    & CALL FLAG_ERROR("Call parameter set distributed vector is not associated.",ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                    & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The field set type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SET_TYPE,"*",ERR,ERROR))// &
                  & " is invalid. The field set type must be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_SET_TYPES,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been created on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
              & " is invalid. The variable type must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been finished."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SET_VECTOR_GET")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SET_VECTOR_GET",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SET_VECTOR_GET")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SET_VECTOR_GET

  !
  !================================================================================================================================
  !

  !>Finalises the parameter sets for a field and deallocates all memory.
  SUBROUTINE FIELD_PARAMETER_SETS_FINALISE(FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE) :: FIELD_VARIABLE !<A pointer to the field variable to finalise the parameter sets for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: parameter_set_idx

    CALL ENTERS("FIELD_PARAMETER_SETS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE)) DEALLOCATE(FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE)
    IF(ASSOCIATED(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS)) THEN
      DO parameter_set_idx=1,SIZE(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS,1)
        CALL FIELD_PARAMETER_SET_FINALISE(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(parameter_set_idx)%PTR,ERR,ERROR,*999)
      ENDDO !parameter_set_idx
      DEALLOCATE(FIELD_VARIABLE%PARAMETER_SETS%PARAMETER_SETS)
    ENDIF
    FIELD_VARIABLE%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS=0

    CALL EXITS("FIELD_PARAMETER_SETS_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_PARAMETER_SETS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SETS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SETS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the parameter sets for a field.
  SUBROUTINE FIELD_PARAMETER_SETS_INITIALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the variable parameter sets for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,parameter_set_idx,variable_idx
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("FIELD_PARAMETER_SETS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(FIELD)) THEN
      DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES        
        FIELD%VARIABLES(variable_idx)%PARAMETER_SETS%FIELD_VARIABLE=>FIELD%VARIABLES(variable_idx)
        FIELD%VARIABLES(variable_idx)%PARAMETER_SETS%NUMBER_OF_PARAMETER_SETS=0
        NULLIFY(FIELD%VARIABLES(variable_idx)%PARAMETER_SETS%PARAMETER_SETS)
        ALLOCATE(FIELD%VARIABLES(variable_idx)%PARAMETER_SETS%SET_TYPE(FIELD_NUMBER_OF_SET_TYPES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field parameter sets set types.",ERR,ERROR,*999)
        DO parameter_set_idx=1,FIELD_NUMBER_OF_SET_TYPES
          NULLIFY(FIELD%VARIABLES(variable_idx)%PARAMETER_SETS%SET_TYPE(parameter_set_idx)%PTR)
        ENDDO !parameter_set_idx
        !Create a field values parameter set
        CALL FIELD_PARAMETER_SET_CREATE(FIELD,FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
      ENDDO !variable_idx
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_PARAMETER_SETS_INITIALISE")
    RETURN
999 DO variable_idx=1,FIELD_NUMBER_OF_VARIABLE_TYPES
      IF(ASSOCIATED(FIELD%VARIABLE_TYPE_MAP(variable_idx)%PTR)) &
        & CALL FIELD_PARAMETER_SETS_FINALISE(FIELD%VARIABLE_TYPE_MAP(variable_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
    ENDDO !variable_idx
998 CALL ERRORS("FIELD_PARAMETER_SETS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_PARAMETER_SETS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_PARAMETER_SETS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the scaling for a field scaling index and deallocates all memory.
  SUBROUTINE FIELD_SCALING_FINALISE(FIELD,SCALING_INDEX,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finalise the scalings for
    INTEGER(INTG), INTENT(IN) :: SCALING_INDEX !<The scaling index to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(SCALING_INDEX>0.AND.SCALING_INDEX<=FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES) THEN
        !IF(ALLOCATED(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS))  &
        !  & DEALLOCATE(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS)
        CALL DISTRIBUTED_VECTOR_DESTROY(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS,ERR,ERROR,*999)
      ELSE
        LOCAL_ERROR="The scaling index of "//TRIM(NUMBER_TO_VSTRING(SCALING_INDEX,"*",ERR,ERROR))// &
          & " is invalid for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES,"*",ERR,ERROR))// &
          & " scaling indices."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALING_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_SCALING_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_SCALING_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_SCALING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the scalings for a field scaling index corresponding to a mesh component index.
  SUBROUTINE FIELD_SCALING_INITIALISE(FIELD,SCALING_INDEX,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the scaling for
    INTEGER(INTG), INTENT(IN) :: SCALING_INDEX !<The scaling index to initialise
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to initialise for the scaling
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALING_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(SCALING_INDEX>0.AND.SCALING_INDEX<=FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES) THEN
        IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=FIELD%DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS) THEN
          FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
          FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%MAX_NUMBER_OF_ELEMENT_PARAMETERS=FIELD%DECOMPOSITION% &
            & DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS
          FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%MAX_NUMBER_OF_DERIVATIVES=FIELD%DECOMPOSITION% &
            & DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES
          NULLIFY(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS)
          SELECT CASE(FIELD%SCALINGS%SCALING_TYPE)
          CASE(FIELD_NO_SCALING)
            !Do nothing
          CASE(FIELD_UNIT_SCALING,FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
            !ALLOCATE(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)% &
            !  & MAX_NUMBER_OF_DERIVATIVES,FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY% &
            !  & NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
            !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate scale factors",ERR,ERROR,*999)
            !FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS=1.0_DP
            CALL DISTRIBUTED_VECTOR_CREATE_START(FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%DOFS, &
              & FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS,ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS, &
              & DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_CREATE_FINISH(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS,ERR,ERROR,*999)
            IF(FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
              !Initialise the scalings to 1.0 for a geometric field. Other field types will be setup in FIELD_SCALINGS_CALCULATE
              CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(FIELD%SCALINGS%SCALINGS(SCALING_INDEX)%SCALE_FACTORS,1.0_DP,ERR,ERROR,*999)
            ENDIF
          CASE(FIELD_ARC_LENGTH_SCALING)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The scaling type of "//TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " is invalid for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          LOCAL_ERROR="The mesh component number of "//TRIM(NUMBER_TO_VSTRING(SCALING_INDEX,"*",ERR,ERROR))// &
            & " is invalid for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
            & " which is associated with a mesh which has "//TRIM(NUMBER_TO_VSTRING(FIELD%DECOMPOSITION% &
            & MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//" mesh components."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The scaling index of "//TRIM(NUMBER_TO_VSTRING(SCALING_INDEX,"*",ERR,ERROR))// &
          & " is invalid for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " which has "//TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES,"*",ERR,ERROR))// &
          & " scaling indices."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALING_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_SCALING_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_SCALING_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_SCALING_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the scale factors from the geometric field associated with the field.
  SUBROUTINE FIELD_SCALINGS_CALCULATE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to calculate the scalings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,ni,ni1,ni2,nk,nk2,nl,nl2,nlp,np,nu,nu1,nu2,ny,ny1,ny2,ny3,scaling_idx
    REAL(DP) :: LENGTH1,LENGTH2,MEAN_LENGTH,TEMP
    REAL(DP), POINTER :: SCALE_FACTORS(:)
    LOGICAL :: FOUND
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: DECOMPOSITION_LINES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD
    TYPE(FIELD_SCALING_TYPE), POINTER :: FIELD_SCALING
    TYPE(FIELD_SCALINGS_TYPE), POINTER :: FIELD_SCALINGS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALINGS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      FIELD_SCALINGS=>FIELD%SCALINGS
      IF(ASSOCIATED(FIELD_SCALINGS)) THEN
        GEOMETRIC_FIELD=>FIELD%GEOMETRIC_FIELD
        IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
          SELECT CASE(FIELD_SCALINGS%SCALING_TYPE)
          CASE(FIELD_NO_SCALING)
            !Do nothing
          CASE(FIELD_UNIT_SCALING)
            DO scaling_idx=1,FIELD_SCALINGS%NUMBER_OF_SCALING_INDICES
              FIELD_SCALING=>FIELD_SCALINGS%SCALINGS(scaling_idx)
              MESH_COMPONENT_NUMBER=FIELD_SCALING%MESH_COMPONENT_NUMBER
              DOMAIN=>FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
              CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(FIELD_SCALING%SCALE_FACTORS,1.0_DP,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_UPDATE_START(FIELD_SCALING%SCALE_FACTORS,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(FIELD_SCALING%SCALE_FACTORS,ERR,ERROR,*999)
            ENDDO !scaling_idx
          CASE(FIELD_ARC_LENGTH_SCALING)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(FIELD_ARITHMETIC_MEAN_SCALING,FIELD_HARMONIC_MEAN_SCALING)
            DO scaling_idx=1,FIELD_SCALINGS%NUMBER_OF_SCALING_INDICES
              FIELD_SCALING=>FIELD_SCALINGS%SCALINGS(scaling_idx)
              MESH_COMPONENT_NUMBER=FIELD_SCALING%MESH_COMPONENT_NUMBER
              DOMAIN=>FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              DOMAIN_LINES=>DOMAIN%TOPOLOGY%LINES
              DECOMPOSITION_LINES=>FIELD%DECOMPOSITION%TOPOLOGY%LINES
              NULLIFY(SCALE_FACTORS)
              CALL DISTRIBUTED_VECTOR_DATA_GET(FIELD_SCALING%SCALE_FACTORS,SCALE_FACTORS,ERR,ERROR,*999)
              DO np=1,DOMAIN_NODES%NUMBER_OF_NODES
                DO nk=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                  ny=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk)
                  nu=DOMAIN_NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(nk)
                  SELECT CASE(nu)
                  CASE(NO_PART_DERIV)
                    CALL DISTRIBUTED_VECTOR_VALUES_SET(FIELD_SCALING%SCALE_FACTORS,ny,1.0_DP,ERR,ERROR,*999)
                  CASE(PART_DERIV_S1,PART_DERIV_S2,PART_DERIV_S3)
                    IF(nu==PART_DERIV_S1) THEN
                      ni=1
                    ELSE IF(nu==PART_DERIV_S2) THEN
                      ni=2
                    ELSE
                      ni=3
                    ENDIF
                    !Find a line of the correct Xi direction going through this node
                    FOUND=.FALSE.
                    DO nlp=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES
                      nl=DOMAIN_NODES%NODES(np)%NODE_LINES(nlp)
                      IF(DECOMPOSITION_LINES%LINES(nl)%XI_DIRECTION==ni) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nnl
                    IF(FOUND) THEN
                      IF(DOMAIN_LINES%LINES(nl)%NODES_IN_LINE(1)==np) THEN !Current node at the beginning of the line
                        nl2=DECOMPOSITION_LINES%LINES(nl)%ADJACENT_LINES(0)
                      ELSE !Current node at the end of the line
                        nl2=DECOMPOSITION_LINES%LINES(nl)%ADJACENT_LINES(1)
                      ENDIF
                      IF(nl2==0) THEN
                        LENGTH1=GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(nl)
                        MEAN_LENGTH=LENGTH1
                      ELSE
                        LENGTH1=GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(nl)
                        LENGTH2=GEOMETRIC_FIELD%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(nl2)
                        SELECT CASE(FIELD_SCALINGS%SCALING_TYPE)
                        CASE(FIELD_ARITHMETIC_MEAN_SCALING)
                          MEAN_LENGTH=(LENGTH1+LENGTH2)/2.0_DP
                        CASE(FIELD_HARMONIC_MEAN_SCALING)
                          TEMP=LENGTH1*LENGTH2
                          IF(ABS(TEMP)>ZERO_TOLERANCE) THEN
                            MEAN_LENGTH=2.0_DP*TEMP/(LENGTH1+LENGTH2)
                          ELSE
                            MEAN_LENGTH=0.0_DP
                          ENDIF
                        CASE DEFAULT
                          LOCAL_ERROR="The scaling type of "// &
                            & TRIM(NUMBER_TO_VSTRING(FIELD_SCALINGS%SCALING_TYPE,"*",ERR,ERROR))//" is invalid."
                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ENDIF
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(FIELD_SCALING%SCALE_FACTORS,ny,MEAN_LENGTH,ERR,ERROR,*999)
                    ELSE
                      LOCAL_ERROR="Could not find a line in the Xi "//TRIM(NUMBER_TO_VSTRING(ni,"*",ERR,ERROR))// &
                        & " direction going through node number "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(PART_DERIV_S1_S2,PART_DERIV_S1_S3,PART_DERIV_S2_S3,PART_DERIV_S1_S2_S3)
                    IF(nu==PART_DERIV_S1_S2) THEN
                      ni1=1
                      nu1=PART_DERIV_S1
                      ni2=2
                      nu2=PART_DERIV_S2
                    ELSE IF(nu==PART_DERIV_S1_S3) THEN
                      ni1=1
                      nu1=PART_DERIV_S1
                      ni2=3
                      nu2=PART_DERIV_S3
                    ELSE IF(nu==PART_DERIV_S1_S2) THEN
                      ni1=2
                      nu1=PART_DERIV_S2
                      ni2=3
                      nu2=PART_DERIV_S3
                    ELSE
                      ni1=1
                      nu1=PART_DERIV_S1
                      ni2=2
                      nu2=PART_DERIV_S2
                    ENDIF
!!TODO: Shouldn't have to search for the nk directions. Store them somewhere.
                    !Find the first direction nk
                    FOUND=.FALSE.
                    DO nk2=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                      IF(DOMAIN_NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(nk2)==nu1) THEN
                        ny1=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk2)
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nk2
                    IF(FOUND) THEN
                      !Find the second direction nk
                      FOUND=.FALSE.
                      DO nk2=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                        IF(DOMAIN_NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(nk2)==nu2) THEN
                          ny2=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk2)
                          FOUND=.TRUE.
                          EXIT
                        ENDIF
                      ENDDO !nk2
                      IF(FOUND) THEN
                        IF(nu==PART_DERIV_S1_S2_S3) THEN
                          !Find the third direction nk
                          FOUND=.FALSE.
                          DO nk2=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                            IF(DOMAIN_NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(nk2)==PART_DERIV_S3) THEN
                              ny3=DOMAIN_NODES%NODES(np)%DOF_INDEX(nk2)
                              FOUND=.TRUE.
                              EXIT
                            ENDIF
                          ENDDO !nk2
                          IF(FOUND) THEN
                            CALL DISTRIBUTED_VECTOR_VALUES_SET(FIELD_SCALING%SCALE_FACTORS,ny, &
                              SCALE_FACTORS(ny1)*SCALE_FACTORS(ny2)*SCALE_FACTORS(ny3),ERR,ERROR,*999)
                          ELSE
                            LOCAL_ERROR="Could not find the first partial derivative in the s3 direction index for "//&
                              & "local node number "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))//"."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL DISTRIBUTED_VECTOR_VALUES_SET(FIELD_SCALING%SCALE_FACTORS,ny,SCALE_FACTORS(ny1)* &
                            & SCALE_FACTORS(ny2),ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Could not find the first partial derivative in the s"// &
                          & TRIM(NUMBER_TO_VSTRING(ni2,"*",ERR,ERROR))//" direction index for "//&
                          & "local node number "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Could not find the first partial derivative in the s"// &
                        & TRIM(NUMBER_TO_VSTRING(ni1,"*",ERR,ERROR))//" direction index for "//&
                        & "local node number "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))//"."
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(nu,"*",ERR,ERROR))// &
                      & " for derivative number "//TRIM(NUMBER_TO_VSTRING(nk,"*",ERR,ERROR))// &
                      & " of local node number "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))//" is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ENDDO !nk
              ENDDO !np
              CALL DISTRIBUTED_VECTOR_UPDATE_START(FIELD_SCALING%SCALE_FACTORS,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(FIELD_SCALING%SCALE_FACTORS,ERR,ERROR,*999)
            ENDDO !scaling_idx
          CASE DEFAULT
            LOCAL_ERROR="The scaling type of "//TRIM(NUMBER_TO_VSTRING(FIELD_SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Field geometric field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field scalings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALINGS_CALCULATE")
    RETURN
999 CALL ERRORS("FIELD_SCALINGS_CALCULATE",ERR,ERROR)
    CALL EXITS("FIELD_SCALINGS_CALCULATE")
    RETURN 1
  END SUBROUTINE FIELD_SCALINGS_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the scalings for a field and deallocates all memory.
  SUBROUTINE FIELD_SCALINGS_FINALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finalise the scalings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: scaling_idx

    CALL ENTERS("FIELD_SCALINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      DO scaling_idx=1,FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES
        CALL FIELD_SCALING_FINALISE(FIELD,scaling_idx,ERR,ERROR,*999)
      ENDDO !scaling_idx
      IF(ALLOCATED(FIELD%SCALINGS%SCALINGS)) DEALLOCATE(FIELD%SCALINGS%SCALINGS)
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALINGS_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_SCALINGS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_SCALINGS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_SCALINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the scaling parameters sets for a field.
  SUBROUTINE FIELD_SCALINGS_INITIALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the scalings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,DUMMY_ERR,NUMBER_OF_MESH_COMPONENTS,scaling_idx,variable_idx
    INTEGER(INTG), ALLOCATABLE :: MESH_COMPONENTS_MAP(:)
    INTEGER(INTG), POINTER :: MESH_COMPONENTS(:)
    TYPE(LIST_TYPE), POINTER :: MESH_COMPONENTS_LIST
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(MESH_COMPONENTS)
    NULLIFY(MESH_COMPONENTS_LIST)

    CALL ENTERS("FIELD_SCALINGS_INITIALISE",ERR,ERROR,*997)

    IF(ASSOCIATED(FIELD)) THEN
      !Calculate the mesh component numbers involved in the field
      CALL LIST_CREATE_START(MESH_COMPONENTS_LIST,ERR,ERROR,*999)
      CALL LIST_DATA_TYPE_SET(MESH_COMPONENTS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
      CALL LIST_INITIAL_SIZE_SET(MESH_COMPONENTS_LIST,FIELD%DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL LIST_CREATE_FINISH(MESH_COMPONENTS_LIST,ERR,ERROR,*999)
      DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
        DO component_idx=1,FIELD%VARIABLES(variable_idx)%NUMBER_OF_COMPONENTS
          CALL LIST_ITEM_ADD(MESH_COMPONENTS_LIST,FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER, &
            & ERR,ERROR,*999)
        ENDDO !component_idx
      ENDDO !variable_idx
      CALL LIST_REMOVE_DUPLICATES(MESH_COMPONENTS_LIST,ERR,ERROR,*999)
      CALL LIST_DETACH_AND_DESTROY(MESH_COMPONENTS_LIST,NUMBER_OF_MESH_COMPONENTS,MESH_COMPONENTS,ERR,ERROR,*999)
      ALLOCATE(MESH_COMPONENTS_MAP(FIELD%DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh components map.",ERR,ERROR,*999)
      MESH_COMPONENTS_MAP=0
      DO component_idx=1,NUMBER_OF_MESH_COMPONENTS
        MESH_COMPONENTS_MAP(MESH_COMPONENTS(component_idx))=component_idx
      ENDDO !component_idx
      !Allocate the scaling indices and initialise them
      FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES=NUMBER_OF_MESH_COMPONENTS
      ALLOCATE(FIELD%SCALINGS%SCALINGS(FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field scalings.",ERR,ERROR,*999)
      DO scaling_idx=1,FIELD%SCALINGS%NUMBER_OF_SCALING_INDICES
        CALL FIELD_SCALING_INITIALISE(FIELD,scaling_idx,MESH_COMPONENTS(scaling_idx),ERR,ERROR,*999)
      ENDDO !scaling_idx
      !Set the scaling index for all the field variable components
      DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
        DO component_idx=1,FIELD%VARIABLES(variable_idx)%NUMBER_OF_COMPONENTS
          FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)%SCALING_INDEX= &
            & MESH_COMPONENTS_MAP(FIELD%VARIABLES(variable_idx)%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER)
        ENDDO !component_idx
      ENDDO !variable_idx
      DEALLOCATE(MESH_COMPONENTS)
      IF(FIELD%TYPE/=FIELD_GEOMETRIC_TYPE) CALL FIELD_SCALINGS_CALCULATE(FIELD,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*997)
    ENDIF

    CALL EXITS("FIELD_SCALINGS_INITIALISE")
    RETURN
999 IF(ASSOCIATED(MESH_COMPONENTS)) DEALLOCATE(MESH_COMPONENTS)
    IF(ASSOCIATED(MESH_COMPONENTS_LIST)) CALL LIST_DESTROY(MESH_COMPONENTS_LIST,ERR,ERROR,*998)
998 CALL FIELD_SCALINGS_FINALISE(FIELD,DUMMY_ERR,DUMMY_ERROR,*997)
997 CALL ERRORS("FIELD_SCALINGS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_SCALINGS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_SCALINGS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Checks the scaling type for a field.
  SUBROUTINE FIELD_SCALING_TYPE_CHECK(FIELD,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the scaling type for
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE !<The scaling type for the specified field to check \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALING_TYPE_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        SELECT CASE(SCALING_TYPE)
        CASE(FIELD_NO_SCALING)
          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            LOCAL_ERROR="Invalid scaling type. The scaling type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " which is not no scaling."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_UNIT_SCALING)
          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_UNIT_SCALING) THEN
            LOCAL_ERROR="Invalid scaling type. The scaling type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " which is not unit scaling."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_ARC_LENGTH_SCALING)
          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_ARC_LENGTH_SCALING) THEN
            LOCAL_ERROR="Invalid scaling type. The scaling type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " which is not arc length scaling."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_ARITHMETIC_MEAN_SCALING)
          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_ARITHMETIC_MEAN_SCALING) THEN
            LOCAL_ERROR="Invalid scaling type. The scaling type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " which is not arithmetic mean scaling."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_HARMONIC_MEAN_SCALING)
          IF(FIELD%SCALINGS%SCALING_TYPE/=FIELD_HARMONIC_MEAN_SCALING) THEN
            LOCAL_ERROR="Invalid scaling type. The scaling type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%SCALINGS%SCALING_TYPE,"*",ERR,ERROR))// &
              & " which is not harmonic mean scaling."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The specified scaling type of "//TRIM(NUMBER_TO_VSTRING(SCALING_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALING_TYPE_CHECK")
    RETURN
999 CALL ERRORS("FIELD_SCALING_TYPE_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_SCALING_TYPE_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_SCALING_TYPE_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the scaling type for a field.
  SUBROUTINE FIELD_SCALING_TYPE_GET(FIELD,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
     TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the scaling type for
    INTEGER(INTG), INTENT(OUT) :: SCALING_TYPE !<On return, the scaling type for the specified field to get \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALING_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        SCALING_TYPE=FIELD%SCALINGS%SCALING_TYPE
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALING_TYPE_GET")
    RETURN
999 CALL ERRORS("FIELD_SCALING_TYPE_GET",ERR,ERROR)
    CALL EXITS("FIELD_SCALING_TYPE_GET")
    RETURN 1
  END SUBROUTINE FIELD_SCALING_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field.
  SUBROUTINE FIELD_SCALING_TYPE_SET(FIELD,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the scaling type for
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE !<The scaling type to set \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALING_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%SCALING_TYPE_LOCKED) THEN
            LOCAL_ERROR="The field scaling type has been locked for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            SELECT CASE(SCALING_TYPE)
            CASE(FIELD_NO_SCALING)
              FIELD%SCALINGS%SCALING_TYPE=FIELD_NO_SCALING
            CASE(FIELD_UNIT_SCALING)
              FIELD%SCALINGS%SCALING_TYPE=FIELD_UNIT_SCALING
            CASE(FIELD_ARC_LENGTH_SCALING)
              FIELD%SCALINGS%SCALING_TYPE=FIELD_ARC_LENGTH_SCALING
            CASE(FIELD_ARITHMETIC_MEAN_SCALING)
              FIELD%SCALINGS%SCALING_TYPE=FIELD_ARITHMETIC_MEAN_SCALING
            CASE(FIELD_HARMONIC_MEAN_SCALING)
              FIELD%SCALINGS%SCALING_TYPE=FIELD_HARMONIC_MEAN_SCALING
            CASE DEFAULT
              LOCAL_ERROR="The specified scaling type of "//TRIM(NUMBER_TO_VSTRING(SCALING_TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALING_TYPE_SET")
    RETURN
999 CALL ERRORS("FIELD_SCALING_TYPE_SET",ERR,ERROR)
    CALL EXITS("FIELD_SCALING_TYPE_SET")
    RETURN 1
  END SUBROUTINE FIELD_SCALING_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the scaling type for a field and locks it so that not further changes can be made.
  SUBROUTINE FIELD_SCALING_TYPE_SET_AND_LOCK(FIELD,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the scaling type for
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE !<The scaling type to set \see FIELD_ROUTINES_ScalingTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_SCALING_TYPE_SET_AND_LOCK",ERR,ERROR,*999)

    CALL FIELD_SCALING_TYPE_SET(FIELD,SCALING_TYPE,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%SCALING_TYPE_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_SCALING_TYPE_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_SCALING_TYPE_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_SCALING_TYPE_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_SCALING_TYPE_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Checks the field type for a field.
  SUBROUTINE FIELD_TYPE_CHECK(FIELD,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The field type to check \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_TYPE_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        SELECT CASE(TYPE)
        CASE(FIELD_GEOMETRIC_TYPE)
          IF(FIELD%TYPE/=FIELD_GEOMETRIC_TYPE) THEN
            LOCAL_ERROR="Invalid field type. The field type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%TYPE,"*",ERR,ERROR))// &
              & " which is not a geometric field."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_FIBRE_TYPE)
          IF(FIELD%TYPE/=FIELD_FIBRE_TYPE) THEN
            LOCAL_ERROR="Invalid field type. The field type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%TYPE,"*",ERR,ERROR))// &
              & " which is not a fibre field."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_GENERAL_TYPE)
          IF(FIELD%TYPE/=FIELD_GENERAL_TYPE) THEN
            LOCAL_ERROR="Invalid field type. The field type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%TYPE,"*",ERR,ERROR))// &
              & " which is not a general field."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE(FIELD_MATERIAL_TYPE)
          IF(FIELD%TYPE/=FIELD_MATERIAL_TYPE) THEN
            LOCAL_ERROR="Invalid field type. The field type for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%TYPE,"*",ERR,ERROR))// &
              & " which is not a material field."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The specified field type of "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
        END SELECT
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_TYPE_CHECK")
    RETURN
999 CALL ERRORS("FIELD_TYPE_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_TYPE_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_TYPE_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the field type for a field identified by a pointer.
  SUBROUTINE FIELD_TYPE_GET(FIELD,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the field type for the specified field \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        TYPE=FIELD%TYPE
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_TYPE_GET")
    RETURN
999 CALL ERRORS("FIELD_TYPE_GET",ERR,ERROR)
    CALL EXITS("FIELD_TYPE_GET")
    RETURN 1
  END SUBROUTINE FIELD_TYPE_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field type for a field.
  SUBROUTINE FIELD_TYPE_SET(FIELD,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The field type to set \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%TYPE_LOCKED) THEN
            LOCAL_ERROR="The field type has been locked for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            SELECT CASE(TYPE)
            CASE(FIELD_GEOMETRIC_TYPE)
              FIELD%TYPE=FIELD_GEOMETRIC_TYPE
              FIELD%GEOMETRIC_FIELD=>FIELD
            CASE(FIELD_FIBRE_TYPE)
              FIELD%TYPE=FIELD_FIBRE_TYPE
              NULLIFY(FIELD%GEOMETRIC_FIELD)
            CASE(FIELD_GENERAL_TYPE)
              FIELD%TYPE=FIELD_GENERAL_TYPE
              NULLIFY(FIELD%GEOMETRIC_FIELD)
            CASE(FIELD_MATERIAL_TYPE)
              FIELD%TYPE=FIELD_MATERIAL_TYPE
              NULLIFY(FIELD%GEOMETRIC_FIELD)
            CASE DEFAULT
              LOCAL_ERROR="The specified field type of "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_TYPE_SET")
    RETURN
999 CALL ERRORS("FIELD_TYPE_SET",ERR,ERROR)
    CALL EXITS("FIELD_TYPE_SET")
    RETURN 1
  END SUBROUTINE FIELD_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field type for a field and locks it so that no further changes can be made.
  SUBROUTINE FIELD_TYPE_SET_AND_LOCK(FIELD,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The field type to set \see FIELD_ROUTINES_FieldTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_TYPE_SET_AND_LOCK",ERR,ERROR,*999)
    
    CALL FIELD_TYPE_SET(FIELD,TYPE,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%TYPE_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_TYPE_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_TYPE_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_TYPE_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_TYPE_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Finds and returns in FIELD a pointer to the field identified by USER_NUMBER in the given REGION. If no field with that USER_NUMBER exists FIELD is left nullified.
  SUBROUTINE FIELD_USER_NUMBER_FIND(USER_NUMBER,REGION,FIELD,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The field user number to find
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the field
    TYPE(FIELD_TYPE), POINTER :: FIELD !<On exit, a pointer to the field with the given user number. If no field with that user number exists in the region the FIELD is null.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: field_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_USER_NUMBER_FIND",ERR,ERROR,*999)

    NULLIFY(FIELD)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%FIELDS)) THEN
        field_idx=1
        DO WHILE(field_idx<=REGION%FIELDS%NUMBER_OF_FIELDS.AND..NOT.ASSOCIATED(FIELD))
          IF(REGION%FIELDS%FIELDS(field_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            FIELD=>REGION%FIELDS%FIELDS(field_idx)%PTR
          ELSE
            field_idx=field_idx+1
          ENDIF
        ENDDO
      ELSE
        LOCAL_ERROR="The fields on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("FIELD_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("FIELD_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE FIELD_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises a field variable and deallocates all memory.
  SUBROUTINE FIELD_VARIABLE_FINALISE(FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_VARIABLE_TYPE) :: FIELD_VARIABLE !<The field variable to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELD_VARIABLE_FINALISE",ERR,ERROR,*999)

    CALL FIELD_VARIABLE_COMPONENTS_FINALISE(FIELD_VARIABLE,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD_VARIABLE%DOMAIN_MAPPING)) THEN
      CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(FIELD_VARIABLE%DOMAIN_MAPPING,ERR,ERROR,*999)
      DEALLOCATE(FIELD_VARIABLE%DOMAIN_MAPPING)
    ENDIF
    CALL FIELD_DOF_TO_PARAM_MAP_FINALISE(FIELD_VARIABLE%DOF_TO_PARAM_MAP,ERR,ERROR,*999)
    CALL FIELD_PARAMETER_SETS_FINALISE(FIELD_VARIABLE,ERR,ERROR,*999)

    CALL EXITS("FIELD_VARIABLE_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a field variable
  SUBROUTINE FIELD_VARIABLE_GET(FIELD,VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the variable for.
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The type of field variable to set. \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE !<On exit, a pointer to the field variable. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(VARIABLE_TYPE>=1.AND.VARIABLE_TYPE<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
            CALL FLAG_ERROR("Field variable is already associated.",ERR,ERROR,*999)
          ELSE
            NULLIFY(FIELD_VARIABLE)
            FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(VARIABLE_TYPE)%PTR
            IF(.NOT.ASSOCIATED(FIELD_VARIABLE)) THEN
              LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                & " has not been defined on field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The field variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
            & " is invalid. The field variable type must be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_GET")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_GET",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_GET")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_GET

  !
  !================================================================================================================================
  !

  !>Initialises a field variable.
  SUBROUTINE FIELD_VARIABLE_INITIALISE(FIELD,VARIABLE_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the variable for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<The variable number of the field to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,DUMMY_ERR,variable_type
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        IF(VARIABLE_NUMBER>=1.AND.VARIABLE_NUMBER<=FIELD%NUMBER_OF_VARIABLES) THEN
          FIELD%VARIABLES(VARIABLE_NUMBER)%VARIABLE_NUMBER=VARIABLE_NUMBER
          variable_type=FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(VARIABLE_NUMBER)
          IF(variable_type>=1.AND.variable_type<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
            FIELD%VARIABLES(VARIABLE_NUMBER)%VARIABLE_TYPE=FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES(VARIABLE_NUMBER)
          ELSE
            LOCAL_ERROR="A field variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
              & " for variable number "//TRIM(NUMBER_TO_VSTRING(VARIABLE_NUMBER,"*",ERR,ERROR))// &
              & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
          ENDIF
          FIELD%VARIABLE_TYPE_MAP(FIELD%VARIABLES(VARIABLE_NUMBER)%VARIABLE_TYPE)%PTR=>FIELD%VARIABLES(VARIABLE_NUMBER)
          FIELD_VARIABLE=>FIELD%VARIABLE_TYPE_MAP(FIELD%VARIABLES(VARIABLE_NUMBER)%VARIABLE_TYPE)%PTR
          FIELD_VARIABLE%FIELD=>FIELD
          FIELD_VARIABLE%REGION=>FIELD%REGION
          FIELD_VARIABLE%DIMENSION=FIELD%CREATE_VALUES_CACHE%DIMENSION(variable_type)
          FIELD_VARIABLE%DATA_TYPE=FIELD%CREATE_VALUES_CACHE%DATA_TYPES(variable_type)
          IF(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_type)>0) THEN
            FIELD_VARIABLE%NUMBER_OF_COMPONENTS=FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_type)
            CALL FIELD_VARIABLE_COMPONENTS_INITIALISE(FIELD,VARIABLE_NUMBER,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="The number of components of "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_type),"*",ERR,ERROR))// &
              & " for variable type "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
              & " of field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The number must be > 0."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
          FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS=-1
          DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            IF(FIELD_VARIABLE%COMPONENTS(component_idx)%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS> &
              & FIELD_VARIABLE%MAX_NUMBER_OF_INTERPOLATION_PARAMETERS) FIELD_VARIABLE% &
              & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS=FIELD_VARIABLE%COMPONENTS(component_idx)% &
              & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
          ENDDO !component_idx
          FIELD_VARIABLE%NUMBER_OF_DOFS=0
          FIELD_VARIABLE%TOTAL_NUMBER_OF_DOFS=0
          FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS=0
          ALLOCATE(FIELD_VARIABLE%DOMAIN_MAPPING,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate field variable domain mapping.",ERR,ERROR,*999)
          CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(FIELD_VARIABLE%DOMAIN_MAPPING, &
            & FIELD%DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
          CALL FIELD_DOF_TO_PARAM_MAP_INITIALISE(FIELD_VARIABLE%DOF_TO_PARAM_MAP,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="Variable number "//TRIM(NUMBER_TO_VSTRING(VARIABLE_NUMBER,"*",ERR,ERROR))// &
            & " is invalid for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" which has "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))//" variables."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Field create values cache is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_INITIALISE")
    RETURN
999 CALL FIELD_VARIABLE_FINALISE(FIELD_VARIABLE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("FIELD_VARIABLE_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Checks the field variable types for a field.
  SUBROUTINE FIELD_VARIABLE_TYPES_CHECK(FIELD,VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to check the variable types for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). The field variable type for the variable_idx'th field variable to check
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_TYPES_CHECK",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(SIZE(VARIABLE_TYPES,1)>=FIELD%NUMBER_OF_VARIABLES) THEN
          DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
            IF(VARIABLE_TYPES(variable_idx)>=1.AND.VARIABLE_TYPES(variable_idx)<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
              IF(FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE/=VARIABLE_TYPES(variable_idx)) THEN
                LOCAL_ERROR="Invalid variable type. The variable type for variable index number "// &
                  & TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))//" of field number "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                  & TRIM(NUMBER_TO_VSTRING(FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE,"*",ERR,ERROR))// &
                  & " which is does correspond to the specified variable_type of "// &
                  & TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPES(variable_idx),"*",ERR,ERROR))// &
                & " at position number "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                & " is invalid. The variable type must be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                  
            ENDIF
          ENDDO !variable_idx
        ELSE
          LOCAL_ERROR="Invalid variable types. The size of the specified variable types array is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VARIABLE_TYPES,1),"*",ERR,ERROR))//" and it must be >= "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_TYPES_CHECK")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_TYPES_CHECK",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_TYPES_CHECK")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_TYPES_CHECK

  !
  !================================================================================================================================
  !

  !>Gets the field variable types for a field.
  SUBROUTINE FIELD_VARIABLE_TYPES_GET(FIELD,VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to get the variable types for
    INTEGER(INTG), INTENT(OUT) :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). On return, the field variable type variable_idx'th field variable
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_TYPES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        IF(SIZE(VARIABLE_TYPES,1)>=FIELD%NUMBER_OF_VARIABLES) THEN
          VARIABLE_TYPES=0
          DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
            VARIABLE_TYPES(variable_idx)=FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
          ENDDO !variable_idx
        ELSE
          LOCAL_ERROR="Invalid variable types. The size of the specified variable types array is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(VARIABLE_TYPES,1),"*",ERR,ERROR))//" and it must be >= "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has not been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_TYPES_GET")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_TYPES_GET",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_TYPES_GET")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_TYPES_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field variable types for a field.
  SUBROUTINE FIELD_VARIABLE_TYPES_SET(FIELD,VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the type for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). The field variable type for the variable_idx'th field variable to set 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_COMPONENTS,old_variable_type,variable_idx,variable_idx2,variable_type
    INTEGER(INTG) :: OLD_DIMENSION(FIELD_NUMBER_OF_VARIABLE_TYPES),OLD_DATA_TYPES(FIELD_NUMBER_OF_VARIABLE_TYPES), &
      & OLD_NUMBER_OF_COMPONENTS(FIELD_NUMBER_OF_VARIABLE_TYPES)
    INTEGER(INTG), ALLOCATABLE :: OLD_VARIABLE_TYPES(:),OLD_INTERPOLATION_TYPE(:,:),OLD_MESH_COMPONENT_NUMBER(:,:)
    LOGICAL :: OLD_DIMENSION_LOCKED(FIELD_NUMBER_OF_VARIABLE_TYPES),OLD_DATA_TYPES_LOCKED(FIELD_NUMBER_OF_VARIABLE_TYPES), &
      & OLD_NUMBER_OF_COMPONENTS_LOCKED(FIELD_NUMBER_OF_VARIABLE_TYPES)
    LOGICAL, ALLOCATABLE :: OLD_INTERPOLATION_TYPE_LOCKED(:,:),OLD_MESH_COMPONENT_NUMBER_LOCKED(:,:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_TYPES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(FIELD%FIELD_FINISHED) THEN
        LOCAL_ERROR="Field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
          & " has been finished."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
          IF(FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES_LOCKED) THEN
            LOCAL_ERROR="The field variable types has been locked for field number "// &
              & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//" and can not be changed."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(SIZE(VARIABLE_TYPES,1)==FIELD%NUMBER_OF_VARIABLES) THEN
              DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
                variable_type=VARIABLE_TYPES(variable_idx)
                !Check that the variable type is in range
                IF(variable_type<1.OR.VARIABLE_TYPE>FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
                  LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                    & " at position number "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                    & " is invalid. The variable type must be between 1 and "// &
                    & TRIM(NUMBER_TO_VSTRING(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                  
                ENDIF
                !Check that the variable type is not repeated
                DO variable_idx2=variable_idx+1,FIELD%NUMBER_OF_VARIABLES
                  IF(VARIABLE_TYPES(variable_idx2)==variable_type) THEN
                    LOCAL_ERROR="The specified variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                      & " occurs at position number "//TRIM(NUMBER_TO_VSTRING(variable_idx,"*",ERR,ERROR))// &
                      & " and position number "//TRIM(NUMBER_TO_VSTRING(variable_idx2,"*",ERR,ERROR))// &
                      & ". The variable types must be unique."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                    
                  ENDIF
                ENDDO !variable_idx2
              ENDDO !variable_idx
              NUMBER_OF_COMPONENTS=SIZE(FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE,1)
              ALLOCATE(OLD_VARIABLE_TYPES(FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old variable types.",ERR,ERROR,*999)
              ALLOCATE(OLD_INTERPOLATION_TYPE(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type.",ERR,ERROR,*999)
              ALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type locked.",ERR,ERROR,*999)
              ALLOCATE(OLD_MESH_COMPONENT_NUMBER(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old mesh component number.",ERR,ERROR,*999)
              ALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED(NUMBER_OF_COMPONENTS,FIELD_NUMBER_OF_VARIABLE_TYPES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old mesh component number locked.",ERR,ERROR,*999)
              OLD_VARIABLE_TYPES=FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES
              OLD_DIMENSION=FIELD%CREATE_VALUES_CACHE%DIMENSION
              OLD_DIMENSION_LOCKED=FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED
              OLD_DATA_TYPES=FIELD%CREATE_VALUES_CACHE%DATA_TYPES
              OLD_DATA_TYPES_LOCKED=FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED
              OLD_NUMBER_OF_COMPONENTS=FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS
              OLD_NUMBER_OF_COMPONENTS_LOCKED=FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED
              OLD_INTERPOLATION_TYPE=FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE
              OLD_INTERPOLATION_TYPE_LOCKED=FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED
              OLD_MESH_COMPONENT_NUMBER=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER
              OLD_MESH_COMPONENT_NUMBER_LOCKED=FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED
              FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES=0
              FIELD%CREATE_VALUES_CACHE%DIMENSION=0
              FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED=.FALSE.
              FIELD%CREATE_VALUES_CACHE%DATA_TYPES=0
              FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED=.FALSE.
              FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS=0
              FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED=.FALSE.
              FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE=0
              FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED=.FALSE.
              FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER=0
              FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED=.FALSE.
              DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES                
                variable_type=VARIABLE_TYPES(variable_idx)
                old_variable_type=OLD_VARIABLE_TYPES(variable_idx)
                FIELD%CREATE_VALUES_CACHE%DIMENSION(variable_type)=OLD_DIMENSION(old_variable_type)
                FIELD%CREATE_VALUES_CACHE%DIMENSION_LOCKED(variable_type)=OLD_DIMENSION_LOCKED(old_variable_type)
                FIELD%CREATE_VALUES_CACHE%DATA_TYPES(variable_type)=OLD_DATA_TYPES(old_variable_type)
                FIELD%CREATE_VALUES_CACHE%DATA_TYPES_LOCKED(variable_type)=OLD_DATA_TYPES_LOCKED(old_variable_type)
                FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS(variable_type)=OLD_NUMBER_OF_COMPONENTS(old_variable_type)
                FIELD%CREATE_VALUES_CACHE%NUMBER_OF_COMPONENTS_LOCKED(variable_type)=OLD_NUMBER_OF_COMPONENTS_LOCKED( &
                  & old_variable_type)
                FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE(:,variable_type)=OLD_INTERPOLATION_TYPE(:,old_variable_type)
                FIELD%CREATE_VALUES_CACHE%INTERPOLATION_TYPE_LOCKED(:,variable_type)=OLD_INTERPOLATION_TYPE_LOCKED(:, &
                  & old_variable_type)
                FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER(:,variable_type)=OLD_MESH_COMPONENT_NUMBER(:,old_variable_type)
                FIELD%CREATE_VALUES_CACHE%MESH_COMPONENT_NUMBER_LOCKED(:,variable_type)=OLD_MESH_COMPONENT_NUMBER_LOCKED(:, &
                  & old_variable_type)
              ENDDO !variable_idx
              FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES=VARIABLE_TYPES
              DEALLOCATE(OLD_VARIABLE_TYPES)
              DEALLOCATE(OLD_INTERPOLATION_TYPE)
              DEALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED)
              DEALLOCATE(OLD_MESH_COMPONENT_NUMBER)
              DEALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED)              
            ELSE
              LOCAL_ERROR="Invalid variable types. The size of the specified variable types array is "// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(VARIABLE_TYPES,1),"*",ERR,ERROR))// &
                & " and the number of variables for field number "//TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))// &
                & " is "//TRIM(NUMBER_TO_VSTRING(FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="Field create values cache is not associated for field number "// &
            & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_TYPES_SET")
    RETURN
999 IF(ALLOCATED(OLD_VARIABLE_TYPES)) DEALLOCATE(OLD_VARIABLE_TYPES)
    IF(ALLOCATED(OLD_INTERPOLATION_TYPE)) DEALLOCATE(OLD_INTERPOLATION_TYPE)
    IF(ALLOCATED(OLD_INTERPOLATION_TYPE_LOCKED)) DEALLOCATE(OLD_INTERPOLATION_TYPE_LOCKED)
    IF(ALLOCATED(OLD_MESH_COMPONENT_NUMBER)) DEALLOCATE(OLD_MESH_COMPONENT_NUMBER)
    IF(ALLOCATED(OLD_MESH_COMPONENT_NUMBER_LOCKED)) DEALLOCATE(OLD_MESH_COMPONENT_NUMBER_LOCKED)
    CALL ERRORS("FIELD_VARIABLE_TYPES_SET",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_TYPES_SET")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_TYPES_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field variable types for a field and locks it so that no further changes can be made.
  SUBROUTINE FIELD_VARIABLE_TYPES_SET_AND_LOCK(FIELD,VARIABLE_TYPES,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to set the type for
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPES(:) !<VARIABLE_TYPES(variable_idx). The field variable type for the variable_idx'th field variable to set 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FIELD_VARIABLE_TYPES_SET_AND_LOCK",ERR,ERROR,*999)
    
    CALL FIELD_VARIABLE_TYPES_SET(FIELD,VARIABLE_TYPES,ERR,ERROR,*999)
    IF(ASSOCIATED(FIELD)) THEN
      IF(ASSOCIATED(FIELD%CREATE_VALUES_CACHE)) THEN
        FIELD%CREATE_VALUES_CACHE%VARIABLE_TYPES_LOCKED=.TRUE.
      ELSE
        LOCAL_ERROR="Field create values cache is not associated for field number "// &
          & TRIM(NUMBER_TO_VSTRING(FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLE_TYPES_SET_AND_LOCK")
    RETURN
999 CALL ERRORS("FIELD_VARIABLE_TYPES_SET_AND_LOCK",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLE_TYPES_SET_AND_LOCK")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLE_TYPES_SET_AND_LOCK

  !
  !================================================================================================================================
  !

  !>Finalises the field variables for a field and deallocates all memory.
  SUBROUTINE FIELD_VARIABLES_FINALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to finalise the variables for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx

    CALL ENTERS("FIELD_VARIABLES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ALLOCATED(FIELD%VARIABLES)) THEN
        DO variable_idx=1,SIZE(FIELD%VARIABLES,1)
          CALL FIELD_VARIABLE_FINALISE(FIELD%VARIABLES(variable_idx),ERR,ERROR,*999)
        ENDDO !variable_idx
        DEALLOCATE(FIELD%VARIABLES)
      ENDIF
      FIELD%NUMBER_OF_VARIABLES=0
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLES_FINALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLES_FINALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLES_FINALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the field variables.
  SUBROUTINE FIELD_VARIABLES_INITIALISE(FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field to initialise the variables for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_idx

    CALL ENTERS("FIELD_VARIABLES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIELD)) THEN
      IF(ALLOCATED(FIELD%VARIABLES)) THEN
        CALL FLAG_ERROR("Field already has associated variables.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(FIELD%VARIABLES(FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new field variables.",ERR,ERROR,*999)
        DO variable_idx=1,FIELD%NUMBER_OF_VARIABLES
          CALL FIELD_VARIABLE_INITIALISE(FIELD,variable_idx,ERR,ERROR,*999)
        ENDDO !variable_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELD_VARIABLES_INITIALISE")
    RETURN
999 CALL ERRORS("FIELD_VARIABLES_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELD_VARIABLES_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELD_VARIABLES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the fields for the given region and deallocates all memory.
  SUBROUTINE FIELDS_FINALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise the fields for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD

    CALL ENTERS("FIELDS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%FIELDS)) THEN
        DO WHILE(REGION%FIELDS%NUMBER_OF_FIELDS>0)
          FIELD=>REGION%FIELDS%FIELDS(1)%PTR
          CALL FIELD_DESTROY(FIELD,ERR,ERROR,*999)
        ENDDO !field_idx
        DEALLOCATE(REGION%FIELDS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELDS_FINALISE")
    RETURN
999 CALL ERRORS("FIELDS_FINALISE",ERR,ERROR)
    CALL EXITS("FIELDS_FINALISE")
    RETURN 1
  END SUBROUTINE FIELDS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the fields for the given region.
  SUBROUTINE FIELDS_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the fields for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("FIELDS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%FIELDS)) THEN
        CALL FLAG_ERROR("Region already has fields associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(REGION%FIELDS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Region fields could not be allocated.",ERR,ERROR,*999)
        !!TODO: Inherit any fields from the parent region
        REGION%FIELDS%NUMBER_OF_FIELDS=0
        NULLIFY(REGION%FIELDS%FIELDS)
        REGION%FIELDS%REGION=>REGION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FIELDS_INITIALISE")
    RETURN
999 CALL ERRORS("FIELDS_INITIALISE",ERR,ERROR)
    CALL EXITS("FIELDS_INITIALISE")
    RETURN 1
  END SUBROUTINE FIELDS_INITIALISE

END MODULE FIELD_ROUTINES
