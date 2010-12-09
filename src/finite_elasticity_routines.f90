!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all finite elasticity routines.
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
!> Contributor(s): Kumar Mithraratne, Jack Lee, Alice Hung
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

!>This module handles all finite elasticity routines.
MODULE FINITE_ELASTICITY_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES  
  USE COORDINATE_ROUTINES  
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS  
  USE MATRIX_VECTOR
  USE MESH_ROUTINES
  USE MPI
  USE PROBLEM_CONSTANTS
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices FINITE_ELASTICITY_ROUTINES::AnalyticParamIndices
  !> \brief Indices for EQUATIONS_SET_ANALYTIC_TYPE%ANALYTIC_USER_PARAMS
  !> \see FINITE_ELASTICITY_ROUTINES,OPENCMISS_AnalyticParamIndices
  !>@{
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX=1 !<Inner pressure parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX=2 !<Outer pressure parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX=3 !<Lambda parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX=4 !<Tsi parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX=5 !<Inner radius parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX=6 !<Outer radius parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX=7 !<c1 parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  INTEGER(INTG), PARAMETER :: FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX=8 !<c2 parameter index \see FINITE_ELASTICITY_ROUTINES_AnalyticParamIndices, FINITE_ELASTICITY_ROUTINES
  !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX,FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX

  PUBLIC FINITE_ELASTICITY_ANALYTIC_CALCULATE, FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE, &
    & FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SETUP,FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SETUP, &
    & FINITE_ELASTICITY_POST_SOLVE,FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA, &
    & FINITE_ELASTICITY_PRE_SOLVE,FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem
  SUBROUTINE FINITE_ELASTICITY_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: node_idx,component_idx,deriv_idx,variable_idx,dim_idx,local_ny,variable_type
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,user_node,global_node,local_node
    REAL(DP) :: X(3),DEFORMED_X(3),P,VALUE
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN,DOMAIN_PRESSURE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES,DOMAIN_PRESSURE_NODES
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    !BC stuff
    INTEGER(INTG),ALLOCATABLE :: INNER_SURFACE_NODES(:),OUTER_SURFACE_NODES(:),TOP_SURFACE_NODES(:),BOTTOM_SURFACE_NODES(:)
    INTEGER(INTG) :: INNER_NORMAL_XI,OUTER_NORMAL_XI,TOP_NORMAL_XI,BOTTOM_NORMAL_XI,MESH_COMPONENT
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER, DOMAIN_NUMBER, MPI_IERROR
    REAL(DP) :: PIN,POUT,LAMBDA,DEFORMED_Z
    LOGICAL :: X_FIXED,Y_FIXED,NODE_EXISTS, X_OKAY,Y_OKAY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_ANALYTIC_CALCULATE",ERR,ERROR,*999)

    MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            !Get access to geometric coordinates
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            MESH_COMPONENT=GEOMETRIC_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)
            !Assign BC here - it's complicated so separate from analytic calculations
            NULLIFY(BOUNDARY_CONDITIONS)
            CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
            IF(ASSOCIATED(DECOMPOSITION)) THEN
              MESH=>DECOMPOSITION%MESH
              IF(ASSOCIATED(MESH)) THEN
                GENERATED_MESH=>MESH%GENERATED_MESH
                IF(ASSOCIATED(GENERATED_MESH)) THEN
                  NODES_MAPPING=>DECOMPOSITION%DOMAIN(1)%PTR%MAPPINGS%NODES   !HACK - ALL CHECKING INTERMEDIATE SKIPPED
                  IF(ASSOCIATED(NODES_MAPPING)) THEN
                    !Get surfaces (hardcoded): fix two nodes on the bottom face, pressure conditions inside & outside
                    CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,1_INTG, &
                        & INNER_SURFACE_NODES,INNER_NORMAL_XI,ERR,ERROR,*999) !Inner
                    CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,2_INTG, &
                        & OUTER_SURFACE_NODES,OUTER_NORMAL_XI,ERR,ERROR,*999) !Outer
                    CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,3_INTG, &
                        & TOP_SURFACE_NODES,TOP_NORMAL_XI,ERR,ERROR,*999) !Top
                    CALL GENERATED_MESH_SURFACE_GET(GENERATED_MESH,MESH_COMPONENT,4_INTG, &
                        & BOTTOM_SURFACE_NODES,BOTTOM_NORMAL_XI,ERR,ERROR,*999) !Bottom
                    !Set all inner surface nodes to inner pressure (- sign is to make positive P into a compressive force) ?
                    PIN=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX)
                    DO node_idx=1,SIZE(INNER_SURFACE_NODES,1)
                      user_node=INNER_SURFACE_NODES(node_idx)
                      !Need to test if this node is in current decomposition
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                        CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                          & user_node,ABS(INNER_NORMAL_XI),BOUNDARY_CONDITION_PRESSURE_INCREMENTED,PIN,ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                    !Set all outer surface nodes to outer pressure
                    POUT=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX)
                    DO node_idx=1,SIZE(OUTER_SURFACE_NODES,1)
                      user_node=OUTER_SURFACE_NODES(node_idx)
                      !Need to test if this node is in current decomposition
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                        CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                          & user_node,ABS(OUTER_NORMAL_XI),BOUNDARY_CONDITION_PRESSURE_INCREMENTED,POUT,ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                    !Set all top nodes fixed in z plane at lambda*height
                    LAMBDA=EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX)
                    DO node_idx=1,SIZE(TOP_SURFACE_NODES,1)
                      user_node=TOP_SURFACE_NODES(node_idx)
                      !Need to test if this node is in current decomposition
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                        CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,1,user_node,NODE_EXISTS,global_node,ERR,ERROR,*999)
                        IF(.NOT.NODE_EXISTS) CYCLE
                        CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,ERR,ERROR,*999)
                        local_ny=GEOMETRIC_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,local_node)
                        DEFORMED_Z=GEOMETRIC_PARAMETERS(local_ny)*LAMBDA
                        CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,1, &
                          & user_node,ABS(TOP_NORMAL_XI),BOUNDARY_CONDITION_FIXED,DEFORMED_Z,ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                    !Set all bottom nodes fixed in z plane
                    DO node_idx=1,SIZE(BOTTOM_SURFACE_NODES,1)
                      user_node=BOTTOM_SURFACE_NODES(node_idx)
                      !Need to check this node exists in the current domain
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                        CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,1, &
                          & user_node,ABS(BOTTOM_NORMAL_XI),BOUNDARY_CONDITION_FIXED,0.0_DP,ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                    !Set two nodes on the bottom surface to axial displacement only:
                    !Easier for parallel: Fix everything that can be fixed !!!
                    X_FIXED=.FALSE.
                    Y_FIXED=.FALSE.
                    DO node_idx=1,SIZE(BOTTOM_SURFACE_NODES,1)
                      user_node=BOTTOM_SURFACE_NODES(node_idx)
                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,ERR,ERROR,*999)
                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                        CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,1,user_node,NODE_EXISTS,global_node,ERR,ERROR,*999)
                        IF(.NOT.NODE_EXISTS) CYCLE
                        CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,ERR,ERROR,*999)
                        local_ny=GEOMETRIC_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,local_node)
                        X(1)=GEOMETRIC_PARAMETERS(local_ny)
                          CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,1,user_node,NODE_EXISTS,global_node,ERR,ERROR,*999)
                          IF(.NOT.NODE_EXISTS) CYCLE
                          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,ERR,ERROR,*999)
                          local_ny=GEOMETRIC_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,local_node)
                        X(2)=GEOMETRIC_PARAMETERS(local_ny)
                        IF(ABS(X(1))<1E-7_DP) THEN
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,1, &
                            & user_node,1,BOUNDARY_CONDITION_FIXED,0.0_DP,ERR,ERROR,*999)
!                           WRITE(*,*) "COMPUTATIONAL NODE ",MY_COMPUTATIONAL_NODE_NUMBER," user node",user_node, &
!                             & "FIXED IN X DIRECTION"
                          X_FIXED=.TRUE.
                        ENDIF
                        IF(ABS(X(2))<1E-7_DP) THEN
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,1, &
                            & user_node,2,BOUNDARY_CONDITION_FIXED,0.0_DP,ERR,ERROR,*999)
!                           WRITE(*,*) "COMPUTATIONAL NODE ",MY_COMPUTATIONAL_NODE_NUMBER," user node",user_node, &
!                             & "FIXED IN Y DIRECTION"
                          Y_FIXED=.TRUE.
                        ENDIF
                      ENDIF
                    ENDDO
                    !Check it went well
                    CALL MPI_REDUCE(X_FIXED,X_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
                    CALL MPI_REDUCE(Y_FIXED,Y_OKAY,1,MPI_LOGICAL,MPI_LOR,0,MPI_COMM_WORLD,MPI_IERROR)
                    IF(MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
                      IF(.NOT.(X_OKAY.AND.Y_OKAY)) THEN
                        CALL FLAG_ERROR("Could not fix nodes to prevent rigid body motion",ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain nodes mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Generated mesh is not associated. For the Cylinder analytic solution, it must be available "// &
                    & "for automatic boundary condition assignment",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
            ENDIF

            !Now calculate analytic solution
            DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
              IF(variable_idx==1) CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Global number of dofs : ", &
                & FIELD_VARIABLE%NUMBER_OF_GLOBAL_DOFS,ERR,ERROR,*999)
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                component_idx=1 !Assuming components 1..3 use a common mesh component and 4 uses a different one
                IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                  DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                  IF(ASSOCIATED(DOMAIN)) THEN
                    IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        !Also grab the equivalent pointer for pressure component
                        IF(FIELD_VARIABLE%COMPONENTS(4)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                          DOMAIN_PRESSURE=>FIELD_VARIABLE%COMPONENTS(4)%DOMAIN
                          IF(ASSOCIATED(DOMAIN_PRESSURE)) THEN
                            IF(ASSOCIATED(DOMAIN_PRESSURE%TOPOLOGY)) THEN
                              DOMAIN_PRESSURE_NODES=>DOMAIN_PRESSURE%TOPOLOGY%NODES
                                IF(ASSOCIATED(DOMAIN_PRESSURE_NODES)) THEN

                                !Loop over the local nodes excluding the ghosts.
                                DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                  !!TODO \todo We should interpolate the geometric field here and the node position.
                                  DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                    local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                                    X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                  ENDDO !dim_idx
                                  !Loop over the derivatives
                                  DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                    SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                    CASE(EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER)
                                      !Cylinder inflation, extension, torsion
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !Do all components at the same time (r,theta,z)->(x,y,z) & p
                                          CALL FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE(X, &
                                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS,DEFORMED_X,P,ERR,ERROR,*999)
                                        CASE(GLOBAL_DERIV_S1)
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE(GLOBAL_DERIV_S2)
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE(GLOBAL_DERIV_S1_S2)
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                            DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                            & " is invalid."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx))
                                        CASE(NO_GLOBAL_DERIV)
                                          !Not implemented, but don't want to cause an error so do nothing
                                        CASE(GLOBAL_DERIV_S1)
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE(GLOBAL_DERIV_S2)
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
                                        CASE(GLOBAL_DERIV_S1_S2)
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                            DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx),"*",ERR,ERROR))// &
                                            & " is invalid."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR)) &
                                          & //" is invalid."
                                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                      END SELECT
                                    CASE DEFAULT
                                      LOCAL_ERROR="The analytic function type of "// &
                                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                                        & " is invalid."
                                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                    END SELECT
                                    !Set the analytic solution to parameter set
                                    DO component_idx=1,NUMBER_OF_DIMENSIONS
                                      local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                        & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,DEFORMED_X(component_idx),ERR,ERROR,*999)
                                    ENDDO
                                    !Don't forget the pressure component
                                    user_node=DOMAIN_NODES%NODES(node_idx)%USER_NUMBER
                                    CALL MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,DOMAIN_PRESSURE%MESH_COMPONENT_NUMBER,user_node, &
                                      & NODE_EXISTS,global_node,ERR,ERROR,*999)
                                    IF(NODE_EXISTS) THEN
                                      CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node, &
                                        & DOMAIN_PRESSURE%MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*999)
                                      IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                                        !\todo: test the domain node mappings pointer properly
                                        local_node=DOMAIN_PRESSURE%mappings%nodes%global_to_local_map(global_node)%local_number(1)
                                        local_ny=FIELD_VARIABLE%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(deriv_idx, &
                                          & local_node)
                                        !Because p=2.lambda in this particular constitutive law, we'll assign half the 
                                        !hydrostatic pressure to the analytic array
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                        & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,P/2.0_dp,ERR,ERROR,*999)
                                      ENDIF
                                    ENDIF
                                  ENDDO !deriv_idx
                                ENDDO !node_idx

                              ELSE
                                CALL FLAG_ERROR("Domain for pressure topology node is not associated",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Domain for pressure topology is not associated",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("Domain for pressure component is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Non-nodal based interpolation of pressure cannot be used with analytic solutions", &
                            & ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
                ENDIF
                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
              ENDIF

            ENDDO !variable_idx
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
          ENDIF            
        ELSE
          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    
    CALL EXITS("FINITE_ELASTICITY_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_ANALYTIC_CALCULATE

  !
  !================================================================================================================================
  !

  !>Calcualates the analytic solution (deformed coordinates and hydrostatic pressure) for cylinder inflation+extension+torsion problem
  SUBROUTINE FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE(X,ANALYTIC_USER_PARAMS,DEFORMED_X,P,ERR,ERROR,*)
    !Argument variables
    REAL(DP), INTENT(IN) :: X(:)                !<Undeformed coordinates
    REAL(DP), INTENT(IN) :: ANALYTIC_USER_PARAMS(:) !<Array containing the problem parameters
    REAL(DP), INTENT(OUT) :: DEFORMED_X(3)      !<Deformed coordinates
    REAL(DP), INTENT(OUT) :: P                  !<Hydrostatic pressure at the given material coordintae
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    !Local variables
    REAL(DP) :: PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2 !A1=external radius, A2=internal radius
    REAL(DP) :: MU1,MU2,MU,K
    REAL(DP) :: F,F2,DF
    REAL(DP) :: R,THETA ! Undeformed coordinates in radial coordinates
    REAL(DP) :: DEFORMED_R,DEFORMED_THETA
    REAL(DP) :: DELTA,RES
    REAL(DP), PARAMETER :: STEP=1E-5_DP, RELTOL=1E-12_DP
    

    CALL ENTERS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE",ERR,ERROR,*999)

    !Grab problem parameters
    PIN=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX)
    POUT=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX)
    LAMBDA=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX)
    TSI=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX)
    A1=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX) ! external radius
    A2=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX) ! internal radius
    C1=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX)
    C2=ANALYTIC_USER_PARAMS(FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX)

    !Solve for MU1 - Newton's method (\todo: Implement here, or separate out for general use?)
    MU1=1.0_DP  !Initial guess - need a better way!
    DO
      !Calculate f(MU1)
      F=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2)

      !Calculate f'(MU1) by finite differencing
      F2=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1+STEP,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2)
      DF=(F2-F)/STEP

      !Next increment for MU1
      DELTA=-F/DF

      !Ensure that the step actually reduces residual
      F2=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1+DELTA,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2)
      DO
        IF (ABS(F2)<ABS(F).OR.ABS(F2)==0.0_DP) THEN    ! PASS
          MU1=MU1+DELTA
          EXIT
        ELSEIF (DELTA<1E-3_DP) THEN ! FAIL: It's likely that the initial guess is too far away
          CALL ERRORS('FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE failed to converge.',ERR,ERROR)
        ELSE                        ! KEEP GOING
          DELTA=DELTA/2.0_DP
          F2=FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1+DELTA,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2)
        ENDIF
      ENDDO

      !Test for convergence: relative residual
      RES=DELTA/(1.0_DP+MU1)
      IF (RES<RELTOL) EXIT
    ENDDO

    !Calculate MU2
    MU2=SQRT(((A1/A2)**2*(LAMBDA*MU1**2-1.0_DP)+1.0_DP)/LAMBDA)

    !Calculate radius and angle from undeformed coordinates
    R=SQRT(X(1)**2+X(2)**2)
    THETA=ATAN2(X(2),X(1)) ! in radians

    !Calculate deformed coordinates
    K=A1**2*(LAMBDA*MU1**2-1.0_DP)
    MU=SQRT(1.0_DP/LAMBDA*(1.0_DP+K/R**2))
    DEFORMED_R=MU*R
    DEFORMED_THETA=THETA+TSI*LAMBDA*X(3)
    DEFORMED_X(1)=DEFORMED_R*COS(DEFORMED_THETA)
    DEFORMED_X(2)=DEFORMED_R*SIN(DEFORMED_THETA)
    DEFORMED_X(3)=LAMBDA*X(3)
    
    !Calculate pressure
    P=POUT-(C1/LAMBDA+C2*LAMBDA)*(1.0_DP/LAMBDA/MU1**2-R**2/(R**2+K)+LOG(MU**2/MU1**2))+C1*TSI**2*LAMBDA*(R**2-A1**2) &
      & -2.0_DP*(C1/LAMBDA**2/MU**2+C2*(1.0_DP/LAMBDA**2+1.0_DP/MU**2+TSI**2*R**2))

    CALL EXITS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_CYLINDER_ANALYTIC_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual function required to solve for MU1, in the cylinder analytic example
  FUNCTION FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE(MU1,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2)
    !Argument variables
    REAL(DP) :: FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE
    REAL(DP) :: MU1,PIN,POUT,LAMBDA,TSI,A1,A2,C1,C2
    !Local variables
    REAL(DP) :: MU,K

    K=A1**2*(LAMBDA*MU1**2-1.0_DP)
    MU=SQRT(1.0_DP/LAMBDA*(1.0_DP+K/A2**2))

    FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE= &
      &  2.0_DP*(C1/LAMBDA**2/MU**2 + C2*(1.0_DP/LAMBDA**2+1.0_DP/MU**2+TSI**2*A2**2))+ &
      & POUT-(C1/LAMBDA+C2*LAMBDA)*(1.0_DP/LAMBDA/MU1**2-A2**2/(A2**2+K)+2*LOG(MU/MU1))+ &
      & C1*TSI**2*LAMBDA*(A2**2-A1**2)-2.0_DP*(C1/LAMBDA**2/MU**2+C2*(1.0_DP/LAMBDA**2+ &
      & 1.0_DP/MU**2+TSI**2*A2**2))+PIN
  
    RETURN
  END FUNCTION FINITE_ELASTICITY_CYLINDER_ANALYTIC_FUNC_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian matrix entries using finite differencing for a finite elasticity finite element equations set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set 
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD 
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: PARAMETERS
    REAL(DP),POINTER :: DATA(:) ! parameter_set vector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR1
    INTEGER(INTG) :: component_idx,local_ny,derivative_idx,derivative,node_idx,node,column
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: DEPENDENT_COMPONENT_INTERPOLATION_TYPE
    REAL(DP) :: DELTA, ORIG_DEP_VAR ,xnorm

    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
        FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLES(1) ! 'U' variable
        DEPENDENT_NUMBER_OF_COMPONENTS=DEPENDENT_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        PARAMETERS=>DEPENDENT_FIELD%VARIABLES(1)%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS  ! vector of dependent variables, basically
        
        ! make a temporary copy of the unperturbed residuals
        CALL FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999) ! can't we reuse old results?
        ELEMENT_VECTOR1=NONLINEAR_MATRICES%ELEMENT_RESIDUAL

        ! determine step size
        !\todo: will this be robust enough? Try a fraction of the smallest (largest?) entry
!         DELTA=1e-6_DP
        CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETERS,DATA,ERR,ERROR,*999)
        !xnorm=sqrt(sum(DATA**2))/size(DATA)
        !DELTA=(DELTA+xnorm)*DELTA
        DELTA=MAX(MAXVAL(ABS(DATA))*1E-6_DP,1E-9)
        CALL DISTRIBUTED_VECTOR_DATA_RESTORE(PARAMETERS,DATA,ERR,ERROR,*999)

        ! the actual finite differencing algorithm is about 4 lines but since the parameters are all 
        ! distributed out, have to use proper field accessing routines.. 
        ! so let's just loop over component, node/el, derivative
        column=0  ! element jacobian matrix column number
        DO component_idx=1,DEPENDENT_NUMBER_OF_COMPONENTS
          ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
          DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
          SELECT CASE (DEPENDENT_COMPONENT_INTERPOLATION_TYPE)
          CASE (FIELD_NODE_BASED_INTERPOLATION)  !node based
            BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
            DO node_idx=1,BASIS%NUMBER_OF_NODES
              node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
              DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(derivative_idx,node_idx)
                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(derivative,node)
                ! one-sided finite difference
                CALL DISTRIBUTED_VECTOR_VALUES_GET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR+DELTA,ERR,ERROR,*999)
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
                CALL FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
                CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
                column=column+1
                NONLINEAR_MATRICES%JACOBIAN%ELEMENT_JACOBIAN%MATRIX(:,column)= &
                    & (NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR-ELEMENT_VECTOR1%VECTOR)/DELTA
              ENDDO !derivative_idx
            ENDDO !node_idx
          CASE (FIELD_ELEMENT_BASED_INTERPOLATION) ! element based
            local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP(ELEMENT_NUMBER)
            ! one-sided finite difference
            CALL DISTRIBUTED_VECTOR_VALUES_GET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR+DELTA,ERR,ERROR,*999)
            NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
            CALL FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
            CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
            column=column+1
            NONLINEAR_MATRICES%JACOBIAN%ELEMENT_JACOBIAN%MATRIX(:,column)= &
                & (NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR-ELEMENT_VECTOR1%VECTOR)/DELTA
          CASE DEFAULT
            CALL FLAG_ERROR("Unsupported type of interpolation.",ERR,ERROR,*999)
          END SELECT
        ENDDO

        ! put the original residual back in
        NONLINEAR_MATRICES%ELEMENT_RESIDUAL=ELEMENT_VECTOR1

      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual and RHS vectors for a finite elasticity finite element equations set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,COMPONENT_BASIS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME,COMPONENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS, &
      & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(MESH_ELEMENT_TYPE), POINTER :: MESH_ELEMENT
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENT_MAPPING
    !TYPE(VARYING_STRING) :: LOCAL_ERROR   

    INTEGER(INTG) :: component_idx,component_idx2,parameter_idx,gauss_idx,element_dof_idx,FIELD_VAR_TYPE
    INTEGER(INTG) :: idx,dim_idx,global_element_idx
    INTEGER(INTG) :: NDOFS,mh !,mhs
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI,HYDROSTATIC_PRESSURE_COMPONENT
    INTEGER(INTG) :: NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
    INTEGER(INTG) :: DEPENDENT_COMPONENT_INTERPOLATION_TYPE
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS       
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS
    INTEGER(INTG) :: var1 ! Variable number corresponding to 'U' in single physics case
    INTEGER(INTG) :: var2 ! Variable number corresponding to 'DELUDLEN' in single physics case
    REAL(DP) :: DZDNU(3,3),CAUCHY_TENSOR(3,3)
    REAL(DP) :: DFDZ(64,3) !temporary until a proper alternative is found
    REAL(DP) :: GAUSS_WEIGHTS,Jznu,Jxxi
    REAL(DP) :: THICKNESS ! for elastic membrane
    REAL(DP) :: DARCY_MASS_INCREASE,DARCY_VOL_INCREASE,DARCY_RHO_0_F  !coupling with Darcy model

    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS,BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DEPENDENT_BASIS,COMPONENT_BASIS)
    NULLIFY(EQUATIONS,EQUATIONS_MAPPING,EQUATIONS_MATRICES,NONLINEAR_MATRICES)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME,COMPONENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS)
    NULLIFY(MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(DARCY_DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT)
    NULLIFY(MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(DEPENDENT_BASIS_1)
    NULLIFY(DECOMPOSITION,MESH_ELEMENT)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN 
        !Which variables are we working with - find the variable pair used for this equations set
        !\todo: put in checks for all the objects/mappings below (do we want to do this for every element?)
        var1=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE%VARIABLE_NUMBER ! number for 'U'
        var2=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE%VARIABLE_NUMBER ! number for 'DELUDELN'

        !Grab pointers: matrices, fields, decomposition, basis
        BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
        BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(var2)%PTR
        TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_PRESSURE_CONDITIONS+ &
          & BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_PRESSURE_INCREMENTED_CONDITIONS

        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        EQUATIONS_MAPPING =>EQUATIONS%EQUATIONS_MAPPING

        FIBRE_FIELD    =>EQUATIONS%INTERPOLATION%FIBRE_FIELD
        GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
        MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
        DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD

        DECOMPOSITION  =>DEPENDENT_FIELD%DECOMPOSITION
        MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER

        DOMAIN_ELEMENT_MAPPING=>DECOMPOSITION%DOMAIN(1)%PTR%MAPPINGS%ELEMENTS

        DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS       
        DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
        DEPENDENT_NUMBER_OF_COMPONENTS=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS

        NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        NUMBER_OF_XI=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI

        !Stuff used to check if this element is on the mesh boundary
        global_element_idx=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(ELEMENT_NUMBER)
        MESH_ELEMENT=>DECOMPOSITION%MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(global_element_idx)

        !Initialise tensors and matrices
        DZDNU=0.0_DP
        DO idx=1,3
          DZDNU(idx,idx)=1.0_DP
        ENDDO
        CAUCHY_TENSOR=DZDNU ! copy the identity matrix
        DFDZ=0.0_DP ! (parameter_idx,component_idx)

        !Grab interpolation parameters
        FIELD_VAR_TYPE=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE%VARIABLE_TYPE
        DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR
        GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        MATERIALS_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
          DARCY_DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR
        ENDIF

        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & MATERIALS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)

        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        ENDIF

        !Point interpolation pointer
        GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        MATERIALS_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
          DARCY_DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR
        ENDIF

        !SELECT: Compressible or incompressible cases
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) ! 4 dependent components

          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHTS=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
              !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE.OR. &
              & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & DARCY_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            ENDIF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR(DEPENDENT_INTERPOLATED_POINT, &
              & GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DZDNU,Jxxi,ERR,ERROR,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ELEMENT_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,ERR,ERROR,*999)
            ENDIF

            !Calculate Sigma=1/Jznu.FTF', the Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU, &
              & ELEMENT_NUMBER,gauss_idx,ERR,ERROR,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,CAUCHY_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    CAUCHY_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',ERR,ERROR,*999)
            ENDIF

            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
              DARCY_RHO_0_F=1.0E-03 !rather pass it through shared material field
              DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV) 
              DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F
            ENDIF

            !Calculate dPhi/dZ at the gauss point, Phi is the basis function
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,ELEMENT_NUMBER,gauss_idx,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DFDZ,ERR,ERROR,*999)

            !For membrane theory in 3D space, the final equation is multiplied by thickness.
            IF (EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE .AND. NUMBER_OF_DIMENSIONS == 3) THEN
              THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                &         FIELD_VARIABLE%NUMBER_OF_COMPONENTS, 1)
            ELSE
              THICKNESS = 1.0_DP
            ENDIF

            !Now add up the residual terms
            element_dof_idx=0
            DO component_idx=1,NUMBER_OF_DIMENSIONS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                  & ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1
                  DO component_idx2=1,NUMBER_OF_DIMENSIONS
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                      & GAUSS_WEIGHTS*Jxxi*Jznu*THICKNESS*CAUCHY_TENSOR(component_idx,component_idx2)* &
                      & DFDZ(parameter_idx,component_idx2)
                  ENDDO ! component_idx2 (inner component index)
                ENDDO ! parameter_idx (residual vector loop)
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FLAG_ERROR("Finite elasticity with element based interpolation is not implemented.",ERR,ERROR,*999)
              ENDIF
            ENDDO ! component_idx

            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET%SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                COMPONENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                COMPONENT_QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1 
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                      & GAUSS_WEIGHTS*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                  ELSE
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                      & GAUSS_WEIGHTS*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP)
                  ENDIF
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHTS*Jxxi* &
                    & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                ELSE
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHTS*Jxxi* &
                    & (Jznu-1.0_DP)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(MESH_ELEMENT%BOUNDARY_ELEMENT.AND.TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF

        CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)   ! compressible problem (no pressure component)

          !Loop over gauss points and add up residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHTS=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

            !Interpolate fields at the gauss points
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) THEN
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & DARCY_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999) ! 'FIRST_PART_DERIV' required ???
            ENDIF

            !Calculate F=dZ/dNU at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR(DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT, &
              & FIBRE_INTERPOLATED_POINT,NUMBER_OF_DIMENSIONS,NUMBER_OF_XI,DZDNU,Jxxi,ERR,ERROR,*999)

            !Calculate Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU, &
              & ELEMENT_NUMBER,gauss_idx,ERR,ERROR,*999)

            !Calculate dF/DZ at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,ELEMENT_NUMBER,gauss_idx,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DFDZ,ERR,ERROR,*999)

            !Add up the residual terms
            element_dof_idx=0
            DO component_idx=1,DEPENDENT_NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                  & ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS  
                  element_dof_idx=element_dof_idx+1    
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                    & GAUSS_WEIGHTS*Jxxi*Jznu*(CAUCHY_TENSOR(component_idx,1)*DFDZ(parameter_idx,1)+ &
                    & CAUCHY_TENSOR(component_idx,2)*DFDZ(parameter_idx,2)+ &
                    & CAUCHY_TENSOR(component_idx,3)*DFDZ(parameter_idx,3)) 
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FLAG_ERROR("Finite elasticity with element based interpolation is not implemented.",ERR,ERROR,*999)
              ENDIF
            ENDDO !component_idx
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(MESH_ELEMENT%BOUNDARY_ELEMENT.AND.TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    !
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF

        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    !-----------------------------------------------------------------------------------------------------------------------------------
    ! CHECK ELEMENT RESIDUAL VECTOR
    IF(DIAGNOSTICS5) THEN
      IF( ELEMENT_NUMBER == 1 ) THEN
        NDOFS = 0
        FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLES(var1) ! 'U' variable
        DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
          DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          NDOFS = NDOFS + DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
        END DO
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NDOFS: ",NDOFS,ERR,ERROR,*999)
!         CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"NDOFS: ",NDOFS,ERR,ERROR,*999)

              CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number * (Fin.Elast.):",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number (Fin.Elast.): ", &
          & ELEMENT_NUMBER,ERR,ERROR,*999)
!         CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element Vector for element number (Fin.Elast.): ", &
!           & ELEMENT_NUMBER,ERR,ERROR,*999)
!               DO mhs=1,NDOFS
!                 CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"row number = ",mhs,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
            & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(:), &
            & '("",4(X,E13.6))','4(4(X,E13.6))',ERR,ERROR,*999)
!           CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
!             & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(:), &
!             & '("",4(X,E13.6))','4(4(X,E13.6))',ERR,ERROR,*999)
!           CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," ",ERR,ERROR,*999)
!               END DO
      ENDIF
    ENDIF
    !-----------------------------------------------------------------------------------------------------------------------------------

    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN

999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !Evaluates the surface traction (pressure) term of the equilibrium equation
  SUBROUTINE FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: var1 !<'U' variable number in single-physics case
    INTEGER(INTG), INTENT(IN) :: var2 !<'DELUDELN' variable number in single-physics case
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMP_ELEMENT
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMP_FACE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE
    TYPE(BASIS_TYPE), POINTER :: FACE_BASIS,DEPENDENT_BASIS
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: FACE_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FACE_PRESSURE_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FACE_INTERPOLATED_POINT
    INTEGER(INTG) :: FIELD_VAR_U_TYPE,FIELD_VAR_DUDN_TYPE,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: element_face_idx,face_number,normal_component_idx,gauss_idx
    INTEGER(INTG) :: FACE_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: component_idx,component_idx2,element_base_dof_idx,face_node_idx
    INTEGER(INTG) :: node_derivative_idx,element_dof_idx,element_node_idx,parameter_idx
    INTEGER(INTG) :: face_parameter_idx,face_node_derivative_idx
    REAL(DP) :: GAUSS_WEIGHT,PRESSURE_GAUSS,NORMAL_PROJECTION
    REAL(DP) :: DZDXI(3,3),DZDXIT(3,3),GIJL(3,3),GIJU(3,3),G,SQRT_G
    LOGICAL :: NONZERO_PRESSURE

    CALL ENTERS("FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    NULLIFY(EQUATIONS,DEPENDENT_FIELD)
    NULLIFY(NONLINEAR_MATRICES,DECOMPOSITION)
    NULLIFY(DECOMP_ELEMENT,DOMAIN_ELEMENT)
    NULLIFY(DEPENDENT_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(DECOMP_FACE,DOMAIN_FACE)
    NULLIFY(FACE_BASIS,DEPENDENT_BASIS,FACE_QUADRATURE_SCHEME)
    NULLIFY(FACE_PRESSURE_INTERPOLATION_PARAMETERS,FACE_INTERPOLATED_POINT)

    !Grab pointers of interest
    EQUATIONS=>EQUATIONS_SET%EQUATIONS
    NONLINEAR_MATRICES=>EQUATIONS%EQUATIONS_MATRICES%NONLINEAR_MATRICES
    DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
    DECOMPOSITION  =>DEPENDENT_FIELD%DECOMPOSITION
    DECOMP_ELEMENT=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)
    MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
    DOMAIN_ELEMENT=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)

    !Interpolation parameter for metric tensor
    FIELD_VAR_U_TYPE=EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%JACOBIAN_TO_VAR_MAP%VARIABLE%VARIABLE_TYPE
    FIELD_VAR_DUDN_TYPE=EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE
    DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

    !Surface pressure term calculation: Loop over all faces
    DO element_face_idx=1,DEPENDENT_BASIS%NUMBER_OF_LOCAL_FACES
      face_number=DECOMP_ELEMENT%ELEMENT_FACES(element_face_idx)
      DECOMP_FACE=>DECOMPOSITION%TOPOLOGY%FACES%FACES(face_number)

      !Check if it's a boundary face
      IF(DECOMP_FACE%BOUNDARY_FACE) THEN !!temporary until MESH_FACE (or equivalent) is available (decomp face includes ghost faces?)
        !Grab normal xi direction of the face and the other two xi directions
        normal_component_idx=ABS(DECOMP_FACE%XI_DIRECTION)  ! if xi=0, this can be a negative number
!         FACE_COMPONENTS=OTHER_XI_DIRECTIONS3(normal_component_idx,2:3,1)  !Two xi directions for the current face
        !\todo: will FACE_COMPONENTS be a problem with sector elements? Check this.
        !Get pressure interpolation objects (DELUDELN pressure_values_set_type)
        FACE_PRESSURE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_DUDN_TYPE)%PTR
        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,face_number, &
          & FACE_PRESSURE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        FACE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(var2)%PTR

        !Check if nonzero surface pressure is defined on the face
        NONZERO_PRESSURE=.FALSE.
        IF(ANY(ABS(FACE_PRESSURE_INTERPOLATION_PARAMETERS%PARAMETERS(:,normal_component_idx))>ZERO_TOLERANCE)) THEN
          NONZERO_PRESSURE=.TRUE.
        ENDIF

        !Nonzero surface pressure found?
        IF(NONZERO_PRESSURE) THEN
          !Grab some other pointers
          DOMAIN_FACE=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%FACES%FACES(face_number)
          FACE_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%FACES%FACES(face_number)%BASIS       
          FACE_QUADRATURE_SCHEME=>FACE_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          FACE_NUMBER_OF_GAUSS_POINTS=FACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

          !Start integrating
          ! Note: As the code will look for P(appl) in the *normal* component to the face, the
          !       initial assignment of P(appl) will have to be made appropriately during bc assignment
!\todo: hopefully all quadrature stuff will always match up between face basis and local face stuff.
! Annoying issue here that p(appl) is interpolated using the face_basis, while dZdXI has to be evaluated
! using the 3D face interpolation... many variables are shared, probably supposed to be the same but I 
! can't guarantee it and checking every single thing will be a fair bit of overhead
          DO gauss_idx=1,FACE_NUMBER_OF_GAUSS_POINTS 
            !Sort out gauss weight
            GAUSS_WEIGHT=FACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

            !Interpolate p(appl) at gauss point
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx,FACE_INTERPOLATED_POINT, &
              & ERR,ERROR,*999)
            PRESSURE_GAUSS=FACE_INTERPOLATED_POINT%VALUES(normal_component_idx,1)    !Surface pressure at this gauss point
            IF(DECOMP_FACE%XI_DIRECTION<0) PRESSURE_GAUSS=-PRESSURE_GAUSS            !Crucial detail here

            !Interpolate delx_j/delxi_M = dZdxi at the face gauss point
            DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_U_TYPE)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
            DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_U_TYPE)%PTR
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,element_face_idx,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            DZDXI=DEPENDENT_INTERPOLATED_POINT%VALUES(1:3,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1:3)) !(component,derivative)

            !Calculate covariant metric tensor
            CALL MATRIX_TRANSPOSE(DZDXI,DZDXIT,ERR,ERROR,*999)
            CALL MATRIX_PRODUCT(DZDXIT,DZDXI,GIJL,ERR,ERROR,*999) !g_ij = dZdXI' * dZdXI
            CALL INVERT(GIJL,GIJU,G,ERR,ERROR,*999) !g^ij = inv(g_ij), G=DET(GIJL)
            SQRT_G=SQRT(G)

            !Loop over 3 components
            DO component_idx=1,3
              !Calculate g^3M*dZ_j/dxi_M
              NORMAL_PROJECTION=dot_product(GIJU(normal_component_idx,:),DZDXI(component_idx,:))
              IF(ABS(NORMAL_PROJECTION)<ZERO_TOLERANCE) CYCLE !Makes it a bit quicker
              !Update element_base_dof_idx
              element_base_dof_idx=0
              DO component_idx2=1,component_idx-1
                element_base_dof_idx=element_base_dof_idx+DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
              ENDDO
              !Looping here is a bit different to reduce redundancy
              DO face_node_idx=1,FACE_BASIS%NUMBER_OF_NODES !nnf
                element_node_idx=DEPENDENT_BASIS%NODE_NUMBERS_IN_LOCAL_FACE(face_node_idx,element_face_idx) !nn
                DO face_node_derivative_idx=1,FACE_BASIS%NUMBER_OF_DERIVATIVES(face_node_idx) !nkf
                  node_derivative_idx=DEPENDENT_BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(face_node_derivative_idx,element_face_idx)
                  parameter_idx=DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(node_derivative_idx,element_node_idx)
                  face_parameter_idx=FACE_BASIS%ELEMENT_PARAMETER_INDEX(face_node_derivative_idx,face_node_idx)
                  element_dof_idx=element_base_dof_idx+parameter_idx
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &  ! sign: double -'s. p(appl) always opposite to normal
                    & GAUSS_WEIGHT*PRESSURE_GAUSS*NORMAL_PROJECTION* &
                    & FACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(face_parameter_idx,NO_PART_DERIV,gauss_idx)* &
                    & SQRT_G
                ENDDO !node_derivative_idx
              ENDDO !face_node_idx
            ENDDO !componenet_dx
          ENDDO !gauss_idx
        ENDIF !nonzero surface pressure check
      ENDIF !boundary face check
    ENDDO !element_face_idx

    CALL EXITS("FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE")
    RETURN

999 CALL ERRORS("FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the deformation gradient tensor at a given Gauss point
  !> Jznu is not used if IS_2D_ELEMENT_IN_3D_SPACE is true
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR(DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,&
    & FIBRE_INTERPOLATED_POINT,DIMEN,NUMBER_OF_XI,DZDNU,Jxxi,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT, &
      & FIBRE_INTERPOLATED_POINT
    INTEGER(INTG), INTENT(IN) :: DIMEN  !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI !<The number of xi directions
    REAL(DP), INTENT(OUT) :: DZDNU(3,3) !<On return, DZDNU - Deformation Gradient Tensor,
    REAL(DP), INTENT(OUT) :: Jxxi       !<On return, The Jacobian of the transformation from xi to x
    INTEGER(INTG), INTENT(OUT) :: ERR   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: DZDNU_TEMP(DIMEN,DIMEN)
    INTEGER(INTG) :: derivative_idx,component_idx,xi_idx
    REAL(DP) :: DNUDX(DIMEN,DIMEN),DNUDXI(DIMEN,DIMEN)
    REAL(DP) :: DXDNU(DIMEN,DIMEN),DXIDNU(DIMEN,DIMEN)
    REAL(DP) :: DXDXI(DIMEN,DIMEN),DZDXI(DIMEN,DIMEN),Jnuxi

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR",ERR,ERROR,*999)

    DO component_idx=1,DIMEN
      DO xi_idx=1,NUMBER_OF_XI
        derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
        DXDXI(component_idx,xi_idx)=GEOMETRIC_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dx/dxi
        DZDXI(component_idx,xi_idx)=DEPENDENT_INTERPOLATED_POINT%VALUES(component_idx,derivative_idx) !dz/dxi
      ENDDO
    ENDDO

    ! Populate the third vector (orthogonal to first and second vector).  This approach requires less modification to the rest of
    ! the code (cleaner), but will have a propagation of error of 10^-4 compare with using non-square (2 x 3) matrices
    IF (DIMEN == 3 .AND. NUMBER_OF_XI == 2) THEN
        CALL CROSS_PRODUCT(DXDXI(:,1),DXDXI(:,2),DXDXI(:,3),ERR,ERROR,*999)
        DXDXI(:,3) = NORMALISE(DXDXI(:,3),ERR,ERROR)
        CALL CROSS_PRODUCT(DZDXI(:,1),DZDXI(:,2),DZDXI(:,3),ERR,ERROR,*999)
        DZDXI(:,3) = NORMALISE(DZDXI(:,3),ERR,ERROR)
    ENDIF

    CALL COORDINATE_MATERIAL_COORDINATE_SYSTEM_CALCULATE(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,DXDNU,ERR,ERROR,*999)

    CALL MATRIX_TRANSPOSE(DXDNU,DNUDX,ERR,ERROR,*999) !dx/dnu is orthogonal. Therefore transpose is its inverse

    CALL MATRIX_PRODUCT(DNUDX,DXDXI,DNUDXI,ERR,ERROR,*999) !dnu/dxi = dnu/dx * dx/dxi
    CALL INVERT(DNUDXI,DXIDNU,Jnuxi,ERR,ERROR,*999) !dxi/dnu

    CALL MATRIX_PRODUCT(DZDXI,DXIDNU,DZDNU_TEMP,ERR,ERROR,*999) !dz/dnu = dz/dxi * dxi/dnu  (deformation gradient tensor, F)

    Jxxi=DETERMINANT(DXDXI,ERR,ERROR)

    IF (DIMEN == 2) THEN
        DZDNU(1,:) = (/DZDNU_TEMP(1,1),DZDNU_TEMP(1,2),0.0_DP/)
        DZDNU(2,:) = (/DZDNU_TEMP(2,1),DZDNU_TEMP(2,2),0.0_DP/)
        DZDNU(3,:) = (/0.0_DP,0.0_DP,1.0_DP/)
    ELSE
        DZDNU = DZDNU_TEMP
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR

  !
  !================================================================================================================================
  !

  !>Evaluates the Cauchy stress tensor at a given Gauss point
!   SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
!       & MATERIALS_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)

  SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU, &
      & ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DARCY_DEPENDENT_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: CAUCHY_TENSOR(:,:)
    REAL(DP), INTENT(OUT) :: Jznu !Determinant of deformation gradient tensor (AZL)
    REAL(DP), INTENT(IN) :: DZDNU(3,3)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE !<The equation subtype
    INTEGER(INTG) :: i,j,PRESSURE_COMPONENT
    REAL(DP) :: AZL(3,3),AZU(3,3),DZDNUT(3,3),PIOLA_TENSOR(3,3),E(3,3),P,ACTIVTIME
    REAL(DP) :: I1,I2,I3            !Invariants, if needed
    REAL(DP) :: TEMP(3,3),TEMPTERM  !Temporary variables
    REAL(DP) :: PIOLA_TENSOR_ADDITION(3,3)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP), DIMENSION (:), POINTER :: C !Parameters for constitutive laws
    REAL(DP) :: a, B(3,3), Q !Parameters for orthotropic laws
    REAL(DP) :: DARCY_MASS_INCREASE  !coupled elasticity Darcy
    INTEGER(INTG) :: DARCY_MASS_INCREASE_ENTRY !position of mass-increase entry in dependent-variable vector

    !CALL ENTERS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR,*999)
    EQUATIONS_SET_SUBTYPE = EQUATIONS_SET%SUBTYPE

    C => MATERIALS_INTERPOLATED_POINT%VALUES(:,1)

    !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
    !AZU - deformed contravariant tensor; I3 = det(C)
    !E = Green-Lagrange strain tensor = 0.5*(C-I)
    !PIOLA_TENSOR is the second Piola-Kirchoff tensor (PK2 or S)
    !P is the actual hydrostatic pressure, not double it

    CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
    CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)

    PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
    P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,1)

    CALL INVERT(AZL,AZU,I3,ERR,ERROR,*999)
    Jznu=I3**0.5_DP
    E = 0.5_DP*AZL 
    DO i=1,3
      E(i,i)=E(i,i)-0.5_DP
    ENDDO

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+p*(I3-1)
      !Also assumed I3 = det(AZL) = 1.0
      !  Note that because PIOLA = 2.del{W}/del{C}=[...]+2.lambda.J^2.C^{-1}
      !  lambda here is actually half of hydrostatic pressure

      !If subtype is membrane, assume Mooney Rivlin constitutive law
      IF (EQUATIONS_SET_SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
          PIOLA_TENSOR(1,3)=2.0_DP*(       C(2)*(-AZL(3,1))        +P*AZU(1,3))
          PIOLA_TENSOR(2,3)=2.0_DP*(       C(2)*(-AZL(3,2))        +P*AZU(2,3))
          PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
          PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
          PIOLA_TENSOR(3,3)=2.0_DP*(C(1)+C(2)*(AZL(1,1)+AZL(2,2))+P*AZU(3,3))
      ELSE
        ! Membrane Equations
        ! Assume incompressible => I3 = 1 => C33(C11 x C22 - C12*C21) = 1
        AZL(3,3) = 1.0_DP / ((AZL(1,1) * AZL(2,2)) - (AZL(1,2) * AZL (2,1)))
        ! Assume Mooney-Rivlin constitutive relation
        P = -1*((C(1) + C(2) * (AZL(1,1) + AZL(2,2))) * AZL(3,3))
        ! Assume stress normal to the surface is neglible i.e. PIOLA_TENSOR(:,3) = 0,PIOLA_TENSOR(3,:) = 0
        PIOLA_TENSOR(:,3) = 0.0_DP
        PIOLA_TENSOR(3,:) = 0.0_DP
      ENDIF
      PIOLA_TENSOR(1,1)=2.0_DP*(C(1)+C(2)*(AZL(2,2)+AZL(3,3))+P*AZU(1,1))
      PIOLA_TENSOR(1,2)=2.0_DP*(     C(2)*(-AZL(2,1))        +P*AZU(1,2))
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=2.0_DP*(C(1)+C(2)*(AZL(3,3)+AZL(1,1))+P*AZU(2,2))

    CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1/2 (e^Q - 1)
      ! where Q=c2(I1-3)
      ! S = 2*dW/dC + 2pC^-1
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)

      TEMPTERM=C(1)*C(2)*EXP(C(2)*(AZL(1,1)+AZL(2,2)+AZL(3,3)-3.0_DP))
      PIOLA_TENSOR(1,1)=TEMPTERM+2.0_DP*P*AZU(1,1)
      PIOLA_TENSOR(1,2)=2.0_DP*P*AZU(1,2)
      PIOLA_TENSOR(1,3)=2.0_DP*P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=TEMPTERM+2.0_DP*P*AZU(2,2)
      PIOLA_TENSOR(2,3)=2.0_DP*P*AZU(1,2)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=TEMPTERM+2.0_DP*P*AZU(3,3)
    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE) 
      !Form of constitutive model is:
      ! W=c1/2 (e^Q - 1)
      ! where Q=2c2(E11+E22+E33)+c3(E11^2)+c4(E22^2+E33^2+E23^2+E32^2)+c5(E12^2+E21^2+E31^2+E13^2)
      ! with E expressed in fibre coordinates

      TEMPTERM=C(1)*EXP(2.0*C(2)*(E(1,1)+E(2,2)+E(3,3))+C(3)*E(1,1)**2+C(4)*(E(2,2)**2+E(3,3)**2+2.0_DP*E(2,3)**2)+ &
          & C(5)*2.0_DP*(E(1,2)**2+E(1,3)**2))
      PIOLA_TENSOR(1,1)=(C(2)+C(3)*E(1,1))*TEMPTERM+2.0_DP*P*AZU(1,1)
      PIOLA_TENSOR(1,2)=C(5)*E(1,2)*TEMPTERM+2.0_DP*P*AZU(1,2)
      PIOLA_TENSOR(1,3)=C(5)*E(1,3)*TEMPTERM+2.0_DP*P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=(C(2)+C(4)*E(2,2))*TEMPTERM+2.0_DP*P*AZU(2,2)
      PIOLA_TENSOR(2,3)=C(4)*E(2,3)*TEMPTERM+2.0_DP*P*AZU(2,3)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=(C(2)+C(4)*E(3,3))*TEMPTERM+2.0_DP*P*AZU(3,3)
    CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) !added by Robert 2010-01-23
      !Form of constitutive model is:
      ! W=a/2 (e^Q - 1)
      ! where Q=[b_ff 2b_fs 2b_fn b_ss 2b_sn b_nn]'* [E_ff E_fs E_fn E_ss E_sn E_nn].^2;
      ! f,s,n denotes the fibre sheet and sheet-normal direction
      a = MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      B(1,1) = MATERIALS_INTERPOLATED_POINT%VALUES(1+1,1)
      B(1,2) = MATERIALS_INTERPOLATED_POINT%VALUES(1+2,1)
      B(1,3) = MATERIALS_INTERPOLATED_POINT%VALUES(1+3,1)
      B(2,1) = B(1,2);
      B(2,2) = MATERIALS_INTERPOLATED_POINT%VALUES(1+4,1)
      B(2,3) = MATERIALS_INTERPOLATED_POINT%VALUES(1+5,1)
      B(3,1) = B(1,3);
      B(3,2) = B(2,3);
      B(3,3) = MATERIALS_INTERPOLATED_POINT%VALUES(1+6,1)
      !write(*,*) "using costa material a", a,"B(1,1)=",B(1,1),"B(1,2)=",B(1,2),"B(1,3)=",B(1,3)
      Q = 0.0_DP;
      DO i=1,3,1
       DO j=1,3,1
         IF (i==j) THEN
              E(i,j) = 0.5_DP * (AZL(i,j)-1);
         ELSE 
              E(i,j) = 0.5_DP * AZL(i,j);
         ENDIF
         Q = Q + B(i,j) * E(i,j) * E(i,j)
       ENDDO
      ENDDO
      Q = exp(Q);
      DO i=1,3,1
       DO j=1,3,1
         PIOLA_TENSOR(i,j)=a*B(i,j)*E(i,j)*Q + p*AZU(i,j);
         !write(*,*) PIOLA_TENSOR(i,j)
       ENDDO
      ENDDO

      IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
        CALL FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
             & EQUATIONS_SET%EQUATIONS%INTERPOLATION%MATERIALS_FIELD, PIOLA_TENSOR(1,1),E(1,1),         &
             & ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*999)
      ENDIF
    CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+c3*(J-1)^2   (this is actually nearly incompressible)
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)

      PIOLA_TENSOR(1,1)=C(1)+C(2)*(AZL(2,2)+AZL(3,3))
      PIOLA_TENSOR(1,2)=C(2)*(-AZL(2,1))
      PIOLA_TENSOR(1,3)=C(2)*(-AZL(3,1))   
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=C(1)+C(2)*(AZL(3,3)+AZL(1,1))
      PIOLA_TENSOR(2,3)=C(2)*(-AZL(3,2))     
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=C(1)+C(2)*(AZL(1,1)+AZL(2,2))
      PIOLA_TENSOR=PIOLA_TENSOR*2.0_DP

      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  C(1) = ",C(1),ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  C(2) = ",C(2),ERR,ERROR,*999)
        CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
          & 3,3,AZL,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    AZL','(",I1,",:)',' :",3(X,E13.6))', &
          & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      ENDIF

      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) THEN
        C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
        PIOLA_TENSOR=PIOLA_TENSOR+2.0_DP*C(3)*(I3-SQRT(I3))*AZU
      ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
        !Starting point for both models is above compressible form of 2nd PK tensor
        SELECT CASE (EQUATIONS_SET_SUBTYPE)
        CASE (EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) !Nearly incompressible
          C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
          !Adjust for the modified Ciarlet-Geymonat expression: Eq.(22) of the INRIA paper
          ! Question is: What deviation is to be penalized : (J-1) or (J-1-m/rho) ??? Probably the latter !
          ! However, m/rho is a given 'constant' and, upon differentiation, drops out.
          ! But it is important to retain I3 = J^2, since J ~ 1 + m/rho /= 1
          PIOLA_TENSOR=PIOLA_TENSOR+C(3)*(SQRT(I3)-1.0_DP)*AZU
          DARCY_MASS_INCREASE_ENTRY = 5 !fifth entry
          DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(DARCY_MASS_INCREASE_ENTRY,NO_PART_DERIV)
          CALL EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION(AZL,AZU,DARCY_MASS_INCREASE,PIOLA_TENSOR_ADDITION,ERR,ERROR,*999)
          PIOLA_TENSOR = PIOLA_TENSOR + PIOLA_TENSOR_ADDITION
        CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) !Incompressible
          !Constitutive model: W=c1*(I1-3)+c2*(I2-3)+p*(I3-1) 
          ! The term 'p*(I3-1)' gives rise to: '2p I3 AZU'
          ! Retain I3 = J^2, since J ~ 1 + m/rho /= 1 
          IF(DIAGNOSTICS1) THEN
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  I3 = ",I3,ERR,ERROR,*999)
          ENDIF
          !\ToDo: Check: Add the volumetric part (as is done for incompressible Mooney-Rivlin above) and double it ???
          DO i=1,3
           DO j=1,3
             !\todo: Check which form to use. Recall: I3 = (1 + m/rho)^2 here !
!              PIOLA_TENSOR(i,j) = PIOLA_TENSOR(i,j) + 2.0_DP*P*I3*AZU(i,j)
             PIOLA_TENSOR(i,j) = PIOLA_TENSOR(i,j) + 2.0_DP*P*AZU(i,j)
           ENDDO
          ENDDO
          DARCY_MASS_INCREASE_ENTRY = 4 !fourth entry
          DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(DARCY_MASS_INCREASE_ENTRY,NO_PART_DERIV)
          CALL EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION(AZL,AZU,DARCY_MASS_INCREASE,PIOLA_TENSOR_ADDITION,ERR,ERROR,*999)
          PIOLA_TENSOR = PIOLA_TENSOR + PIOLA_TENSOR_ADDITION
        CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_MR_SUBTYPE)
          !Constitutive model: W=C1*(J1-3)+C2*(J2-3)+C3*(J-1)^2+lambda.(J-1-m/rho)
          !J1 and J2 are the modified invariants, adjusted for volume change (J1=I1*J^(-2/3), J2=I2*J^(-4/3))
          !Strictly speaking this law isn't for an incompressible material, but the fourth equation in the elasticity
          !is used to satisfy a subtly different constraint, which is to require the solid portion of the poroelastic
          !material retains its volume. (This law is applied on the whole pororous body).
          
          C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
          C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
          C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)

          !J1 term: del(J1)/del(C)=J^(-2/3)*I-2/3*I_1*J^(-2/3)*C^-1
          TEMPTERM=Jznu**(-2.0_DP/3.0_DP)
          TEMP(1,1)=TEMPTERM
          TEMP(2,2)=TEMPTERM
          TEMP(3,3)=TEMPTERM
          I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
          PIOLA_TENSOR=C(1)* (TEMP-2.0_DP/3.0_DP*I1*TEMPTERM*AZU)

          !J2 term: del(J2)/del(C)=J^(-4/3)*del(I2)/del(C) -4/3*I_2*J^(-4/3)*C^-1
          TEMP=matmul(AZL,AZL)  ! C^2
          I2=0.5_DP*(I1**2-(TEMP(1,1)+TEMP(2,2)+TEMP(3,3)))
          TEMPTERM=Jznu**(-4.0_DP/3.0_DP)
          !TEMP is now del(I2)/del(C)
          TEMP(1,1)=AZL(2,2)+AZL(3,3)
          TEMP(1,2)=-2.0_DP*AZL(1,2)
          TEMP(1,3)=-2.0_DP*AZL(1,3)
          TEMP(2,1)=TEMP(1,2)
          TEMP(2,2)=AZL(1,1)+AZL(3,3)
          TEMP(2,3)=-2.0_DP*AZL(2,3)
          TEMP(3,1)=TEMP(1,3)
          TEMP(3,2)=TEMP(2,3)
          TEMP(3,3)=AZL(1,1)+AZL(2,2)
          PIOLA_TENSOR=PIOLA_TENSOR+C(2)* (TEMPTERM*TEMP-4.0_DP/3.0_DP*I2*TEMPTERM*AZU)
          
          !J (det(F)) term: (2.C3.(J-1)+lambda)*J.C^-1
          PIOLA_TENSOR=PIOLA_TENSOR+(2.0_DP*C(3)*(Jznu-1.0_DP)+P)*Jznu*AZU

          !Don't forget, it's wrt C so there is a factor of 2
          PIOLA_TENSOR=2.0_DP*PIOLA_TENSOR
        END SELECT


        IF(DIAGNOSTICS1) THEN
          CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
            & 3,3,PIOLA_TENSOR-PIOLA_TENSOR_ADDITION,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))','(17X,3(X,E13.6))',ERR,ERROR,*999)
          CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
            & 3,3,PIOLA_TENSOR_ADDITION, &
            & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR_ADDITION','(",I1,",:)',' :",3(X,E13.6))', &
            & '(17X,3(X,E13.6))',ERR,ERROR,*999)
        ENDIF
      ENDIF

    CASE (EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE) ! added by Thomas 2010-04-13
      !Form of the constitutive model is:
      ! W = a/(2*b)*exp[b*(I1-3)] + sum_(i=f,s)[a_i/(2*b_i)*(exp[b_i*(I4i-1)^2]-1)] + a_fs/(2*b_fs)*(exp[b_fs*I8fs^2]-1)
      !Also assumed I3 = det(AZL) = J^2 = 1.0  -  incompressible material
      !Assume directions: fibre f_0=[1 0 0], sheet s_0=[0 1 0], (sheet) normal n_0=[0 0 1]
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1) !a
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1) !b
      C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1) !a_f
      C(4)=MATERIALS_INTERPOLATED_POINT%VALUES(4,1) !a_s
      C(5)=MATERIALS_INTERPOLATED_POINT%VALUES(5,1) !b_f
      C(6)=MATERIALS_INTERPOLATED_POINT%VALUES(6,1) !b_s
      C(7)=MATERIALS_INTERPOLATED_POINT%VALUES(7,1) !a_fs
      C(8)=MATERIALS_INTERPOLATED_POINT%VALUES(8,1) !b_fs
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMPTERM=C(1)*EXP(C(2)*(I1-3.0_DP))
      PIOLA_TENSOR(1,1)=-P*AZU(1,1)+TEMPTERM+2.0_DP*C(3)*(AZL(1,1)-1.0_DP)*EXP(C(5)*(AZL(1,1)-1.0_DP)**2.0_DP)
      PIOLA_TENSOR(1,2)=-P*AZU(1,2)+C(7)*AZL(1,2)*EXP(C(8)*AZL(1,2)**2.0_DP)
      PIOLA_TENSOR(1,3)=-P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=-P*AZU(2,2)+TEMPTERM+2.0_DP*C(4)*(AZL(2,2)-1.0_DP)*EXP(C(6)*(AZL(2,2)-1.0_DP)**2.0_DP)
      PIOLA_TENSOR(2,3)=-P*AZU(2,3)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=-P*AZU(3,3)+TEMPTERM
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for a finite elasticity equation type of an elasticity equation set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    CALL MATRIX_PRODUCT(DZDNU,PIOLA_TENSOR,TEMP,ERR,ERROR,*999)
    CALL MATRIX_PRODUCT(TEMP,DZDNUT,CAUCHY_TENSOR,ERR,ERROR,*999)
    
    CAUCHY_TENSOR=CAUCHY_TENSOR/Jznu

    NULLIFY(C)

    !CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR

  !
  !================================================================================================================================
  !

   ! calculates the current active contraction component using the independent field
  ! Uses a hardcoded tension transient based on GPB+NHS with length-dependence for now
  SUBROUTINE FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION(INDEPENDENT_FIELD,MATERIALS_FIELD,PIOLA_FF,E_FF,&
             & ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(FIELD_TYPE), POINTER, INTENT(IN) :: INDEPENDENT_FIELD, MATERIALS_FIELD
    REAL(DP), INTENT(INOUT) :: PIOLA_FF  !<The (1,1)=(fiber,fiber) component of the stress tensor
    REAL(DP), INTENT(IN)    :: E_FF !<E(1,1)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string    

    INTEGER(INTG)  :: I
    REAL(DP) :: S, LAMBDA, OVERLAP, ISO_TA, TA, ACTIVTIME, TIME, DT, PREV_LAMBDA
    REAL(DP), DIMENSION(1:4) :: QL

    REAL(DP), PARAMETER :: PERIOD = 1000 ! 1 Hz
    REAL(DP), PARAMETER, DIMENSION(28) :: TIMES    =    (/ 0, 20, 30, 40, 60, 80, 100, 120, 150, 160, 170, 175, 180, 190, 200,&
    & 225, 250, 300, 333, 366, 400, 450, 500, 600, 700, 800, 900,1000 /) ! simple tension curve based on GPB/NHS: times

    REAL(DP), PARAMETER, DIMENSION(28) :: TENSIONFRAC = (/ 0.0194, 0.0193, 0.0200, 0.0254, 0.0778, 0.1713, 0.2794, 0.3708,&
    & 0.4472, 0.4578, 0.4624, 0.4627, 0.4618, 0.4567, 0.4478, 0.4121, 0.3614, 0.2326, 0.1471, 0.0920, 0.0681, 0.0526, 0.0438,&
    & 0.0332, 0.0271, 0.0234, 0.0210, 0.0194 /) ! simple isometric tension curve based on GPB/NHS: tension/tref 
    REAL(DP), PARAMETER :: T_REF = 100          ! reference tension
  
    CALL ENTERS("FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION",ERR,ERROR,*999)

    ! Get time, dt, etc from independent field
    CALL FIELD_PARAMETER_SET_GET_CONSTANT(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1, DT,ERR,ERROR,*999)   ! dt
    CALL FIELD_PARAMETER_SET_GET_CONSTANT(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 2, TIME,ERR,ERROR,*999) ! time
    DO I=1,4
      CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER, 2+I, QL(I),ERR,ERROR,*999)  ! Q(1) Q(2) Q(3) Lambda for prev in 3/4/5/6
    END DO

    ! get activation time from material field
    CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,&
      &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER, 1, ACTIVTIME,ERR,ERROR,*999)

    LAMBDA = SQRT(2*E_FF + 1)
    TIME =  MAX( MOD(TIME, PERIOD) - ACTIVTIME, 0.0) ! start activation at this time
   
    I = 1
    DO WHILE (TIMES(I) <= TIME) ! find first I such that times(I) >= time
      I = I+1
    END DO
    S    = (TIME - TIMES(I-1)) /  (TIMES(I) - TIMES(I-1))                     !| linear interpolation of ta/tref
    ISO_TA   = T_REF * (TENSIONFRAC(I-1) * (1-S) + TENSIONFRAC(I) * S)        !/ + multiply by tref
  
    CALL FINITE_ELASTICITY_FMM(TIME,DT,QL(4),LAMBDA,QL,ISO_TA,TA)

    QL(4) = LAMBDA  ! bounds applied in FMM, Qi integrated
    DO I=1,4
      CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER, 6+I, QL(I),ERR,ERROR,*999) ! store Q(1) Q(2) Q(3) Lambda for next in 7/8/9/10
    END DO

    PIOLA_FF = PIOLA_FF + TA

    CALL EXITS("FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION")
    RETURN 1  
  END SUBROUTINE FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION

  !
  !================================================================================================================================
  !

  ! Implements length and velocity dependence. can be used in both weak and strong coupling
  SUBROUTINE FINITE_ELASTICITY_FMM(TIME,DT,PREV_LAMBDA,CURR_LAMBDA,Q123,ISO_TA,TA)
    ! PARAMETERS FROM Niederer Hunter & Smith 2006
    REAL(DP), PARAMETER, DIMENSION(1:3) :: A     = (/-29.0,138.0,129.0/)  ! 'A'
    REAL(DP), PARAMETER, DIMENSION(1:3) :: ALPHA = (/0.03,0.13,0.625/)
    REAL(DP), PARAMETER :: la   = 0.35, BETA_0 = 4.9  ! 'a'

    REAL(DP), INTENT(INOUT), DIMENSION(:) :: Q123
    REAL(DP), INTENT(INOUT) :: CURR_LAMBDA
    REAL(DP), INTENT(IN) :: PREV_LAMBDA, DT, TIME, ISO_TA
    REAL(DP), INTENT(OUT) :: TA

    REAL(DP) :: QFAC, DLAMBDA_DT, Q, OVERLAP
    INTEGER(INTG) :: I

    CURR_LAMBDA = MIN(1.15, MAX(0.8, CURR_LAMBDA))  ! inout -> save this

    IF( TIME - 1e-10 <= 0.0) THEN  ! preload / first step -> update method off
      QFAC = 1.0
    ELSE
      DLAMBDA_DT = (CURR_LAMBDA - PREV_LAMBDA) / DT
      DO I=1,3
        Q123(I) = Q123(I) + DT * (A(I) * DLAMBDA_DT - ALPHA(I) * Q123(I))
      END DO
      Q = Q123(1)+Q123(2)+Q123(3)
      IF(Q < 0.0) THEN 
        QFAC = (la*Q + 1.0) / (1.0 - Q)
      ELSE 
        QFAC = (1.0 + (la+2.0)*Q)/(1.0+Q);
      END IF
    END IF

    OVERLAP= 1.0 + BETA_0 * (CURR_LAMBDA-1.0) 
    TA = OVERLAP * QFAC * ISO_TA  ! length dep * vel dep * isometric tension
  END SUBROUTINE FINITE_ELASTICITY_FMM


  !
  !================================================================================================================================
  !

  !>Evaluates df/dz (derivative of interpolation function wrt deformed coord) matrix at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ(INTERPOLATED_POINT,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,NUMBER_OF_DIMENSIONS, &
    & NUMBER_OF_XI,DFDZ,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT !<Interpolated point for the dependent field
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER !<The gauss point number
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI !<The number of xi directions for the interpolation
    REAL(DP), INTENT(OUT) :: DFDZ(:,:) !<On return, a matrix containing the derivatives of the basis functions wrt the deformed coordinates
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: COMPONENT_BASIS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG) :: NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS,derivative_idx,component_idx,xi_idx,parameter_idx 
    REAL(DP) :: DXIDZ(NUMBER_OF_DIMENSIONS,NUMBER_OF_DIMENSIONS),DZDXI(NUMBER_OF_DIMENSIONS,NUMBER_OF_DIMENSIONS)
    REAL(DP) :: Jzxi,DFDXI(NUMBER_OF_DIMENSIONS,64,NUMBER_OF_XI)!temporary until a proper alternative is found
    
    CALL ENTERS("FINITE_ELASTICITY_GAUSS_DFDZ",ERR,ERROR,*999)

    !Initialise DFDXI array
    DFDXI=0.0_DP  ! DFDXI(component_idx,parameter_idx,xi_idx)
    DFDZ=0.0_DP
    DO component_idx=1,NUMBER_OF_DIMENSIONS !Always 3 spatial coordinates (3D)
      DO xi_idx=1,NUMBER_OF_XI !Thus always 3 element coordinates
        derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx)  !2,4,7      
        DZDXI(component_idx,xi_idx)=INTERPOLATED_POINT%VALUES(component_idx,derivative_idx)  !dz/dxi
      ENDDO
    ENDDO

    ! Populate a 3 x 3 square dzdXi if this is a membrane problem in 3D space
    IF (NUMBER_OF_DIMENSIONS == 3 .AND. NUMBER_OF_XI == 2) THEN
        CALL CROSS_PRODUCT(DZDXI(:,1),DZDXI(:,2),DZDXI(:,3),ERR,ERROR,*999)
        DZDXI(:,3) = NORMALISE(DZDXI(:,3),ERR,ERROR)
    ENDIF

    CALL INVERT(DZDXI,DXIDZ,Jzxi,ERR,ERROR,*999) !dxi/dz

    FIELD=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD
    DO component_idx=1,NUMBER_OF_DIMENSIONS
      COMPONENT_BASIS=>FIELD%VARIABLES(1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(ELEMENT_NUMBER)%BASIS
      QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
      DO parameter_idx=1,size(QUADRATURE_SCHEME%GAUSS_BASIS_FNS,1)
        DO xi_idx=1,NUMBER_OF_XI
          derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx)  !2,4,7
          DFDXI(component_idx,parameter_idx,xi_idx)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,derivative_idx, &
            & GAUSS_POINT_NUMBER)
        ENDDO
      ENDDO
    ENDDO

    DO component_idx=1,NUMBER_OF_DIMENSIONS
      COMPONENT_BASIS=>FIELD%VARIABLES(1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(ELEMENT_NUMBER)%BASIS
!       NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=FIELD%VARIABLES(1)%COMPONENTS(component_idx)% &
!         & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
      NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
      DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
          DO xi_idx=1,NUMBER_OF_XI
            DFDZ(parameter_idx,component_idx)= DFDZ(parameter_idx,component_idx) + &
              & DFDXI(component_idx,parameter_idx,xi_idx) * DXIDZ(xi_idx,component_idx)
        ENDDO
      ENDDO
    ENDDO

    CALL EXITS("FINITE_ELASTICITY_GAUSS_DFDZ")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_DFDZ",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_DFDZ")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_COMPONENTS, &
      & NUMBER_OF_DIMENSIONS, I, NUMBER_OF_DARCY_COMPONENTS
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    LOGICAL :: IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD

    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(EQUATIONS_MATERIALS)
    
    IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD = EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE &
                            & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_MEMBRANE_SUBTYPE &
                            & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE

    NUMBER_OF_DIMENSIONS = EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS

    IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
      NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
    ELSE
      NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS
    ENDIF

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,&
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Default to FEM solution method
            CALL FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !\todo Check dimension of geometric field
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
            & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) THEN
            FIELD_VARIABLE=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
            IF(.NOT.ASSOCIATED(FIELD_VARIABLE%PARAMETER_SETS%SET_TYPE(FIELD_MESH_DISPLACEMENT_SET_TYPE)%PTR)) THEN
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                 & FIELD_MESH_DISPLACEMENT_SET_TYPE, ERR, ERROR, *999)
            ENDIF
          END IF
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          !-----------------------------------------------------------------------
          ! Dependent field setup for single-physics
          !-----------------------------------------------------------------------
          CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,&
              & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
              & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
              & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                !Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                                   
                IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
!kmith :09.06.09 - Do we need this ?      
                !Set the hydrostatic component to that of the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)            
!kmith
                ENDIF
                
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                  
                  IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDIF

                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & ERR,ERROR,*999)
                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%PTR)) THEN
                    LOCAL_ERROR="Variable 2 of type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                      & VARIABLE_TYPE,"*",ERR,ERROR))//" does not have a pressure values set type associated."
!                     write(*,*) char(LOCAL_ERROR) ! this is temporary
!                     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

!kmith :09.06.09 - Hydrostatic pressure could be node-based as well
!                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
!                    & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,4, &
!                    & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!kmith

                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity equation"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          !-------------------------------------------------------------------------------
          ! Shared Dependent field setup for multi-physics: elasticity coupled with Darcy
          !-------------------------------------------------------------------------------
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) 
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+2 !for INRIA model: velocity components, pressure, mass increase
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for standard Darcy: velocity components and pressure
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for Darcy with pressure driven by solid: velocity components and mass increase
                END SELECT

                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)

                !Elasticity: Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx

                IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
!kmith :09.0.06.09 - Do we need this ?      
                  !Set the hydrostatic component to that of the first geometric component
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)            
!kmith
                ENDIF

                !Darcy: Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx

                !Darcy: Default pressure and, if present, mass increase to the first geometric component
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_DARCY_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx

              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity: Set the displacement components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

                  IF (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    !Elasticity: Set the hydrostatic pressure component to node based interpolation
                    !as this is used as the pressure field for the Darcy equations
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ELSE IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
                    !Elasticity: Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDIF

                  !Darcy: Set the velocity, pressure and, if present, mass increase components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELVDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE,&
                  & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/),ERR,ERROR,*999)

                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)

                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & ERR,ERROR,*999)

                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) 
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+2 !for INRIA model: velocity components, pressure, mass increase
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for standard Darcy: velocity components and pressure
                CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                  NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for Darcy with pressure driven by solid: velocity components and mass increase
                END SELECT

                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,NUMBER_OF_DARCY_COMPONENTS, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity:
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                  IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    !If solid hydrostatic pressure is driving Darcy flow, check that pressure uses node based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,4, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDIF

!kmith :09.0.06.09 - Hydrostatic pressure could be node-based as well
!                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
!                    & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,4, &
!                    & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!kmith

                  !Darcy:
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

                CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                  CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                    & " is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              ENDIF
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_INITIAL_VALUES_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_RELATIVE_VELOCITY_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity equation"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            !end: Dependent field setup for elasticity coupled with Darcy
          CASE DEFAULT
            LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT


        !-----------------------------------------------------------------
        ! I N d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)

           IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
            NUMBER_OF_COMPONENTS = 10  ! dt t Q1 Q2 Q3 lambda    prev Q1 Q2 Q3 lambda
             IF(EQUATIONS_SET%SOLUTION_METHOD /= EQUATIONS_SET_FEM_SOLUTION_METHOD .OR. &
               &.NOT. EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
             END IF
             CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
               & INDEPENDENT_FIELD,ERR,ERROR,*999)
             CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)

              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)

             CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
               & GEOMETRIC_FIELD,ERR,ERROR,*999)
             CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
               & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

             DO component_idx=1,2 ! dt t  constant 
              CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
             END DO

             DO component_idx=3,NUMBER_OF_COMPONENTS ! other gauss pt based
               CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
             END DO
  
           ELSE ! coupled Darcy problem


            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                & INDEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_INDEPENDENT_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS !+1 !Include hydrostatic pressure component
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DIMENSIONS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              ENDDO !component_idx

! !kmith :09.06.09 - Do we need this ?      
!               !Set the hydrostatic component to that of the first geometric component
!               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
!                 & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
!               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                 & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
!               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
!                 & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)            
! !kmith

            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Set the displacement components to node based interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
!                 !Set the hydrostatic pressure component to element based interpolation
!                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
!                   & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
!                   & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                & ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx

!kmith :09.06.09 - Hydrostatic pressure could be node-based as well
!                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
!                  & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,4, &
!                  & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!kmith

              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ENDIF
           ENDIF ! else darcy
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
              ! initialize values for active contraction independent field. TODO: actual init for z, trpn, or flag to presolve
              IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,Err,ERROR,*999)
                ENDDO
              ENDIF
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT


        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              SELECT CASE(EQUATIONS_SET%SUBTYPE)
              CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
                & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;
              CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE)
                NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS
              CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 5;
              CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS = 7;
              CASE(EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 3;
              CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)
                NUMBER_OF_COMPONENTS = 8;
              CASE DEFAULT
                LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                  & " is not valid for a finite elasticity equation type of an elasticity equation set class."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)

                IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,2,ERR,ERROR,*999) ! 2nd is activation times
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                  &     (/FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE/), ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,ERR,ERROR,*999) ! just 1 component: activation time
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                   & FIELD_V_VARIABLE_TYPE,1 ,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999) ! gauss pt based interp.

                ELSE
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                ENDIF

                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

                DO I=1,NUMBER_OF_COMPONENTS
                !Default the materials components to the geometric interpolation setup with constant interpolation
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)  ! get 1 = x (?) component
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO

                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                !Set the default values for the materials field
                !Don't bother checking equations types, just default to all componets = 1.0
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER)
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE) THEN
                      !Create analytic field if required
                      !Set analtyic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be a Mooney-Rivlin finite elasticity equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a finite elasticity equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
                  CALL FINITE_ELASTICITY_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set analytic has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Start the equations creation
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              ! sander: Quasistatic / Active contraction. correct location?
              IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
              ELSE
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,0,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              ! set structure and storage types
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,MATRIX_BLOCK_STORAGE_TYPE, &
                    & ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
            CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
            IF(ASSOCIATED(EQUATIONS)) THEN
              IF(EQUATIONS%EQUATIONS_FINISHED) THEN
                CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set equations has not been finished.",ERR,ERROR,*999)               
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a finite elasticity equation type of an elasticity equations set class.
  SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_FIELD_USER_NUMBER
    TYPE(FIELD_TYPE), POINTER :: DUMMY_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DUMMY_FIELD_USER_NUMBER=0
      NULLIFY(DUMMY_FIELD)
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_ELASTICITY_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_FINITE_ELASTICITY_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem.
  SUBROUTINE FINITE_ELASTICITY_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Laplace equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop: default is load increment type now
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
            ! sander - Quasistatic: Change 1/2. worth splitting entire case over in copy/paste?
            IF(PROBLEM%SUBTYPE == PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE) THEN
               CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            ENDIF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equatgions
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            ! sander - Quasistatic: Change 2/2. worth splitting entire case over in copy/paste?
            IF(PROBLEM%SUBTYPE == PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE) THEN
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,ERR,ERROR,*999)
            ELSE 
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            END IF
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity problem."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity problem."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a finite elasticity type .
  SUBROUTINE FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_ELASTICITY_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SUBTYPE      
      CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_ELASTICITY_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE      
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem post solve.
  SUBROUTINE FINITE_ELASTICITY_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: I
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD

    CALL ENTERS("FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
              IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE.AND.SOLVER%GLOBAL_NUMBER==1) THEN
                CALL FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              ! how to check eqn subtype? assume active contraction 
              INDEPENDENT_FIELD => SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                                     & EQUATIONS_SETS(1)%PTR%INDEPENDENT%INDEPENDENT_FIELD
              ! store lambda Q (7-10) in prev lambda Q (3-6)
              DO I=1,4
                CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(INDEPENDENT_FIELD,&
                     & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, I+6, &
                     & INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, I+2,ERR,ERROR,*999)
              END DO
              ! output data
              CALL FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_NO_SUBTYPE)
              !CALL FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a finite elasticity problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_POST_SOLVE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_POST_SOLVE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem post solve output data.
  SUBROUTINE FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP !<A pointer to the control time loop.
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(VARYING_STRING) :: METHOD !,FILE
    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE
    LOGICAL :: EXPORT_FIELD
    INTEGER(INTG) :: CURRENT_LOOP_ITERATION
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER
    INTEGER(INTG) :: equations_set_idx,loop_idx

    CALL ENTERS("FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    !Make sure the equations sets are up to date
                    DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                      METHOD="FORTRAN"
                      EXPORT_FIELD=.TRUE.
                      IF(EXPORT_FIELD) THEN          
                        IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity export fields ... ",ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"STATICSOLUTION",ERR,ERROR,*999)
                        ENDIF
                        CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER, &
                          & "STATICSOLIDSOLUTION",ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                  ENDIF 
                ENDIF
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
               & PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR

                    CONTROL_TIME_LOOP=>CONTROL_LOOP
                    DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
                      IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
                        EXIT
                      ENDIF
                      IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
                        CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
                      ELSE
                        CALL FLAG_ERROR("Could not find a time control loop.",ERR,ERROR,*999)
                      ENDIF
                    ENDDO

                    CURRENT_LOOP_ITERATION=CONTROL_TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_TIME_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                        IF(CURRENT_LOOP_ITERATION<10) THEN
                          WRITE(OUTPUT_FILE,'("S_TIMESTP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("S_TIMESTP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("S_TIMESTP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("S_TIMESTP_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        FILE=OUTPUT_FILE
  !                    FILE="TRANSIENT_OUTPUT"
                        METHOD="FORTRAN"
                        EXPORT_FIELD=.TRUE.
                        IF(EXPORT_FIELD) THEN          
                          IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                            IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity export fields ...",ERR,ERROR,*999)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                            ENDIF
                            CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                              & ERR,ERROR,*999)
                            IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity all fields exported ...",ERR,ERROR,*999)
                            ENDIF
                          ENDIF
                        ENDIF 
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a finite elasticity problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Runs before each time loop for a finite elasticity problem.
  SUBROUTINE FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_SOLID !<A pointer to the solid solver
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_SOLID
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD

    NULLIFY(SOLVER_SOLID)
    NULLIFY(CONTROL_LOOP_SOLID)

    CALL ENTERS("FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
            CONTROL_LOOP_SOLID=>CONTROL_LOOP
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            INDEPENDENT_FIELD=>SOLVER_SOLID%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                                   & EQUATIONS_SETS(1)%PTR%INDEPENDENT%INDEPENDENT_FIELD !?
            ! set component 1 to dt
            CALL FIELD_COMPONENT_VALUES_INITIALISE(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_VALUES_SET_TYPE,1,CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT,ERR,ERROR,*999)
            ! set component 2 to current time.
            CALL FIELD_COMPONENT_VALUES_INITIALISE(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_VALUES_SET_TYPE,2,CONTROL_LOOP%TIME_LOOP%CURRENT_TIME,ERR,ERROR,*999)
          CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)
            ! could do this in one line with problem_solver_get but the dependence on problem_routines causes a circular dependence
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/1,CONTROL_LOOP_NODE/),CONTROL_LOOP_SOLID,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            !--- 3.0 For Standard Elasticity Darcy: Update the boundary conditions of the solid
            CALL FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER_SOLID,ERR,ERROR,*999)
          CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/1,1,CONTROL_LOOP_NODE/),CONTROL_LOOP_SOLID,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            !--- 3.0 For Standard Elasticity Darcy: Update the boundary conditions of the solid
            CALL FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER_SOLID,ERR,ERROR,*999)
          CASE(PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/1,CONTROL_LOOP_NODE/),CONTROL_LOOP_SOLID,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            !--- For PGM: Get the displacement field
            CALL FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT(CONTROL_LOOP,SOLVER_SOLID,ERR,ERROR,*999)
          CASE DEFAULT
            !do nothing
        END SELECT
      ELSE
        CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP")
    RETURN 1

  END SUBROUTINE FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity problem pre-solve.
  SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: solver_matrix_idx
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD

    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              ! do nothing, time values get updated in CONTROL_TIME_LOOOP_PRE_LOOP as there might be 
              ! a load increment loop below the time loop, so we don't want to update times here before
              ! every  solve
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
              IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity pre-solve",ERR,ERROR,*999)
              ENDIF

              !--- Set 'SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.'
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
                  IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                    DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                      SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%PTR
                      IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                        SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                      ELSE
                        CALL FLAG_ERROR("Solver Matrix is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                  ELSE
                    CALL FLAG_ERROR("Solver Matrices is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a finite elasticity problem class."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Read in the displacement field for a PGM elasticity problem
  SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FINITE_ELASTICITY  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_FINITE_ELASTICITY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_FINITE_ELASTICITY  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_FINITE_ELASTICITY !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_FINITE_ELASTICITY !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP

    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    REAL(DP), POINTER :: DUMMY_VALUES2(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP) :: NSUBTRACT

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_FINITE_ELASTICITY
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NDOFS_TO_PRINT, I
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    INTEGER(INTG) :: loop_idx

    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT",ERR,ERROR,*999)

!--- \todo : Do we need for each case a FIELD_PARAMETER_SET_UPDATE_START / FINISH on FIELD_MESH_DISPLACEMENT_SET_TYPE ?

    NULLIFY(SOLVER_FINITE_ELASTICITY)
    NULLIFY(MESH_DISPLACEMENT_VALUES)
    NULLIFY(DUMMY_VALUES2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FLAG_ERROR("Could not find a time control loop.",ERR,ERROR,*999)
        ENDIF
      ENDDO

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              !--- Motion: read in from a file
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_FINITE_ELASTICITY,ERR,ERROR,*999)
                SOLVER_EQUATIONS_FINITE_ELASTICITY=>SOLVER_FINITE_ELASTICITY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_FINITE_ELASTICITY)) THEN
                  SOLVER_MAPPING_FINITE_ELASTICITY=>SOLVER_EQUATIONS_FINITE_ELASTICITY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_FINITE_ELASTICITY)) THEN
                    EQUATIONS_SET_FINITE_ELASTICITY=>SOLVER_MAPPING_FINITE_ELASTICITY%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_FINITE_ELASTICITY)) THEN
                      DEPENDENT_FIELD_FINITE_ELASTICITY=>EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD
                    ELSE
                      CALL FLAG_ERROR("Finite elasticity equations set is not associated.",ERR,ERROR,*999)
                    END IF
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity motion read from a file ... ",ERR,ERROR,*999)

                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_FINITE_ELASTICITY%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)

                    !Copy input to Finite elasticity's dependent field
                    !\todo: Still need to take into account that we are reading in displacement,
                    !       while dependent field is the absolute position of the structure
                    INPUT_TYPE=42
                    INPUT_OPTION=2
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, & 
                      & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP)
                    CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("Finite elasticity solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Finite elasticity solver equations are not associated.",ERR,ERROR,*999)
                END IF
  
               IF(DIAGNOSTICS1) THEN
                 NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                 CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                   & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",4(X,E13.6))','4(4(X,E13.6))', &
                   & ERR,ERROR,*999)
               ENDIF
              ELSE
                ! in case of a solver number different from 3: do nothing ???
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a Finite elasticity equation fluid type of a fluid mechanics problem class."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for finite elasticity pre solve
  SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD, GEOMETRIC_FIELD
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_TIME_LOOP

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: GEOMETRIC_FIELD_VALUES(:) 
    REAL(DP), POINTER :: MESH_POSITION_VALUES(:) 
    REAL(DP), POINTER :: DUMMY_VALUES1(:)

    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE
    INTEGER(INTG) :: dof_number,NUMBER_OF_DOFS,GEOMETRY_NUMBER_OF_DOFS,DEPENDENT_NUMBER_OF_DOFS
    INTEGER(INTG) :: NDOFS_TO_PRINT
    INTEGER(INTG) :: loop_idx

!     TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
!     TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE


    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CONTROL_TIME_LOOP=>CONTROL_LOOP
      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL
        IF(CONTROL_TIME_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_TIME_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
          EXIT
        ENDIF
        IF (ASSOCIATED(CONTROL_LOOP%PARENT_LOOP)) THEN
          CONTROL_TIME_LOOP=>CONTROL_TIME_LOOP%PARENT_LOOP
        ELSE
          CALL FLAG_ERROR("Could not find a time control loop.",ERR,ERROR,*999)
        ENDIF
      ENDDO

      IF(ASSOCIATED(SOLVER)) THEN
        IF(SOLVER%GLOBAL_NUMBER==1) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
              CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
                ! do nothing ???
              CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        SELECT CASE(EQUATIONS_SET%SUBTYPE)
                          CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                            & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
                            & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                            IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity update BCs",ERR,ERROR,*999)
                            ENDIF
                            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                            IF(ASSOCIATED(DEPENDENT_FIELD).AND.ASSOCIATED(GEOMETRIC_FIELD)) THEN
                              BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
                              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                EQUATIONS_MAPPING=>EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING
                                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                                  BOUNDARY_CONDITIONS_VARIABLE=>BOUNDARY_CONDITIONS% & 
                                    & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                                  IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                    IF(DIAGNOSTICS1) THEN
                                      NULLIFY( DUMMY_VALUES1 )
                                      CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                      NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                        & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                                        & '(" DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE (bef) = ",4(X,E13.6))', &
                                        & '4(4(X,E13.6))',ERR,ERROR,*999)
                                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                    ENDIF

                                    ! requires solid dependent field and geometry to be interpolated identically !!!
                                    ! assumes that DOFs for dependent and geometric field are stored in the same order
                                    ! How does this routine take into account the BC value ???
                                    ALPHA = 0.10_DP * sin( 2.0_DP * PI * CURRENT_TIME / 4.0_DP )
                                    CALL FIELD_PARAMETER_SETS_COPY(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                      & FIELD_VALUES_SET_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,ALPHA,ERR,ERROR,*999)

                                    NULLIFY(GEOMETRIC_FIELD_VALUES)
                                    CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                      & FIELD_VALUES_SET_TYPE,GEOMETRIC_FIELD_VALUES,ERR,ERROR,*999)

                                    GEOMETRY_NUMBER_OF_DOFS=GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)% &
                                        & PTR%NUMBER_OF_DOFS
                                    DO dof_number=1,GEOMETRY_NUMBER_OF_DOFS
                                      BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% & 
                                        & GLOBAL_BOUNDARY_CONDITIONS(dof_number)
                                      IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL .OR. &
                                        & BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
                                        !--- To obtain absolute positions, add nodal coordinates on top of mesh displacement
                                        CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(GEOMETRIC_FIELD, & 
                                          & FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,dof_number, & 
                                          & GEOMETRIC_FIELD_VALUES(dof_number),ERR,ERROR,*999)
                                      ELSE
                                        ! do nothing ???
                                      END IF
                                    END DO

                                    NULLIFY(MESH_POSITION_VALUES)
                                    CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                      & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_POSITION_VALUES,ERR,ERROR,*999)

                                    DEPENDENT_NUMBER_OF_DOFS=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)% &
                                        & PTR%NUMBER_OF_DOFS
                                    DO dof_number=1,DEPENDENT_NUMBER_OF_DOFS
                                      BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% & 
                                        & GLOBAL_BOUNDARY_CONDITIONS(dof_number)
                                      IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL .OR. &
                                        & BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN

!---tob
                                      !Update FIELD_BOUNDARY_CONDITIONS_SET_TYPE or FIELD_VALUES_SET_TYPE
                                      !(so it is one or the other, but not both) depending on whether or not load increments are used
                                      IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, & 
                                          & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,dof_number, & 
                                          & MESH_POSITION_VALUES(dof_number),ERR,ERROR,*999)
                                      ELSE
                                        !--- Update the dependent field with the new absolute position
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, & 
                                          & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, & 
                                          & MESH_POSITION_VALUES(dof_number),ERR,ERROR,*999)
                                      ENDIF
!---toe

                                      ELSE
                                        ! do nothing ???
                                      END IF
                                    END DO

                                    IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED) THEN
                                      CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD, &
                                        & FIELD_U_VARIABLE_TYPE, FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
                                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD, &
                                        & FIELD_U_VARIABLE_TYPE, FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
                                    ELSE
                                      CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD, &
                                        & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD, &
                                        & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                    ENDIF

                                    IF(DIAGNOSTICS1) THEN
                                      NULLIFY( DUMMY_VALUES1 )
                                      CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                      NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                                      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                        & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                                        & '(" DEPENDENT_FIELD,FIELD_U_VAR_TYPE,FIELD_VALUES_SET_TYPE (after) = ",4(X,E13.6))', &
                                        & '4(4(X,E13.6))',ERR,ERROR,*999)
                                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                        & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                                  END IF

                                  CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                                    & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                ELSE
                                  CALL FLAG_ERROR("EQUATIONS_MAPPING is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Boundary conditions are not associated.",ERR,ERROR,*999)
                              END IF
                            ELSE
                              CALL FLAG_ERROR("Dependent field and/or geometric field is/are not associated.",ERR,ERROR,*999)
                            END IF
                          CASE DEFAULT
                            ! do nothing ???
!                             LOCAL_ERROR="Equations set subtype " &
!                               & //TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
!                               & " is not valid for a standard elasticity Darcy problem subtype."
!                             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                      ELSE
                        CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                    END IF                
                  ELSE
                    CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
                END IF  
              CASE DEFAULT
                ! do nothing ???
!                 LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
!                   & " is not valid for this problem class."
!               CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          ! do nothing ???
!           CALL FLAG_ERROR("PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS may only be carried out for SOLVER%GLOBAL_NUMBER = 1", &
!             & ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS

  !
  !================================================================================================================================
  !

  !>Evaluates the functions f(J) and f'(J); 
  !>  Eq.(21) in Chapelle, Gerbeau, Sainte-Marie, Vignon-Clementel, Computational Mechanics (2010)
  SUBROUTINE EVALUATE_CHAPELLE_FUNCTION(AZL,Jznu,ffact,dfdJfact,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: AZL(3,3) !<C=F'F
    REAL(DP), INTENT(OUT) :: Jznu !<Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    REAL(DP), INTENT(OUT) :: ffact !<f(Jznu) of the INRIA model
    REAL(DP), INTENT(OUT) :: dfdJfact !<dfdJfact = f'(Jznu) of the INRIA model
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    CALL ENTERS("EVALUATE_CHAPELLE_FUNCTION",ERR,ERROR,*999)

    Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    IF( ABS(Jznu) < 1.0E-10_DP ) THEN
      CALL FLAG_ERROR("EVALUATE_CHAPELLE_FUNCTION: ABS(Jznu) < 1.0E-10_DP",ERR,ERROR,*999)
    END IF

    IF( ABS(Jznu-1.0_DP) > 5.0E-02_DP ) THEN
      !Eq.(21) of the INRIA paper
      ffact = 2.0_DP * (Jznu - 1.0_DP - log(Jznu)) / (Jznu - 1.0_DP)**2.0_DP
      dfdJfact = ( 2.0_DP * (1.0_DP - 1.0_DP/Jznu) * (Jznu - 1.0_DP)**2.0_DP &
                   & - 4.0_DP * (Jznu - 1.0_DP - log(Jznu)) * (Jznu - 1.0_DP) ) / (Jznu - 1.0_DP)**4.0_DP
    ELSE
      ffact = 1.0_DP
      dfdJfact = 0.0_DP
    END IF


    CALL EXITS("EVALUATE_CHAPELLE_FUNCTION")
    RETURN
999 CALL ERRORS("EVALUATE_CHAPELLE_FUNCTION",ERR,ERROR)
    CALL EXITS("EVALUATE_CHAPELLE_FUNCTION")
    RETURN 1
  END SUBROUTINE EVALUATE_CHAPELLE_FUNCTION

  !
  !================================================================================================================================
  !

  !>Evaluates the 2nd Piola-Kirchhoff stress tensor;
  !>  Eq.(13) in Chapelle, Gerbeau, Sainte-Marie, Vignon-Clementel, Computational Mechanics (2010)
  SUBROUTINE EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION(AZL,AZU,DARCY_MASS_INCREASE,PIOLA_TENSOR_ADDITION,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: AZL(3,3) !<C=F'F
    REAL(DP), INTENT(IN) :: AZU(3,3) !<inverse of AZL
    REAL(DP), INTENT(IN) :: DARCY_MASS_INCREASE !<mass increase
    REAL(DP), INTENT(OUT) :: PIOLA_TENSOR_ADDITION(3,3) !<Addition to the 2nd Piola-Kirchhoff tensor
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: Jznu !<Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    REAL(DP) :: ffact !<f(Jznu) of the INRIA model
    REAL(DP) :: dfdJfact !<dfdJfact = f'(Jznu) of the INRIA model
    REAL(DP) :: Mfact, bfact  !<INRIA constitutive law
    REAL(DP) :: DARCY_VOL_INCREASE, DARCY_RHO_0_F
    INTEGER(INTG) :: i,j

    
    CALL ENTERS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION",ERR,ERROR,*999)

    !Parameters settings for coupled elasticity Darcy INRIA model:
    !\ToDo: ensure constants are consistent with Darcy model ! Pass through material parameters ???
    DARCY_RHO_0_F = 1.0E-03_DP
!     Mfact = 2.18E05_DP
    Mfact = 2.18E03_DP
    bfact = 1.0_DP

    DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F

!     write(*,*)'DARCY_VOL_INCREASE(2) = ',DARCY_VOL_INCREASE

    CALL EVALUATE_CHAPELLE_FUNCTION(AZL,Jznu,ffact,dfdJfact,ERR,ERROR,*999)

    DO i=1,3
      DO j=1,3
        PIOLA_TENSOR_ADDITION(i,j) = - Mfact * bfact * DARCY_VOL_INCREASE * (ffact + (Jznu - 1.0_DP) * dfdJfact) * Jznu * AZU(i,j) &
          & + 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * dfdJfact * Jznu * AZU(i,j)
      ENDDO
    ENDDO

!     PIOLA_TENSOR_ADDITION = - Mfact * bfact * DARCY_VOL_INCREASE * (ffact + (Jznu - 1.0_DP) * dfdJfact) * Jznu * AZU &
!       & + 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * dfdJfact * Jznu * AZU

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  DARCY_VOL_INCREASE = ",DARCY_VOL_INCREASE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Jznu = ",Jznu,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ffact = ",ffact,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  dfdJfact = ",dfdJfact,ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,AZU,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    AZU','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,PIOLA_TENSOR_ADDITION, &
        & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR_ADDITION','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
    ENDIF


    CALL EXITS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION")
    RETURN
999 CALL ERRORS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION",ERR,ERROR)
    CALL EXITS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION")
    RETURN 1
  END SUBROUTINE EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION

  !
  !================================================================================================================================
  !

END MODULE FINITE_ELASTICITY_ROUTINES
