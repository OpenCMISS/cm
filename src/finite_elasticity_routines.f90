!> \file
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
!> Contributor(s): Kumar Mithraratne, Jack Lee, Alice Hung, Sander Arens
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
  USE CMISS_PETSC
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
  USE FIELD_IO_ROUTINES
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

  INTEGER(INTG), PARAMETER :: TENSOR_TO_VOIGT(3,3)=RESHAPE([1,4,5,4,2,6,5,6,3],[3,3]) !Converts from rank 2 symmetric tensor indices to Voigt index.
  INTEGER(INTG), PARAMETER :: VOIGT_TO_TENSOR(2,6)=RESHAPE([1,1,2,2,3,3,2,1,3,1,3,2],[2,6]) !Converts from Voigt index to rank 2 symmetric tensor.

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_PIN_IDX,FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_POUT_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_LAMBDA_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_TSI_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_RIN_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_ROUT_IDX, &
    & FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C1_IDX, FINITE_ELASTICITY_ANALYTIC_CYLINDER_PARAM_C2_IDX

  PUBLIC FINITE_ELASTICITY_ANALYTIC_CALCULATE, &
    & FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE,FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SETUP,FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SETUP, &
    & FiniteElasticity_ContactProblemSubtypeSet,FiniteElasticity_ContactProblemSetup, & 
    & FINITE_ELASTICITY_POST_SOLVE,FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA, &
    & FINITE_ELASTICITY_PRE_SOLVE,FINITE_ELASTICITY_CONTROL_TIME_LOOP_PRE_LOOP,FiniteElasticity_ControlLoadIncrementLoopPostLoop, &
    & EVALUATE_CHAPELLE_FUNCTION, GET_DARCY_FINITE_ELASTICITY_PARAMETERS, &
    & FiniteElasticityGaussDeformationGradientTensor,FINITE_ELASTICITY_LOAD_INCREMENT_APPLY, &
    & FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE,FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE

  PUBLIC FiniteElasticityEquationsSet_DerivedVariableCalculate, &
    & FiniteElasticity_StrainInterpolateXi

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem
  SUBROUTINE FINITE_ELASTICITY_ANALYTIC_CALCULATE(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: node_idx,component_idx,deriv_idx,variable_idx,dim_idx,local_ny,variable_type
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,user_node,global_node,local_node
    REAL(DP) :: X(3),DEFORMED_X(3),P
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
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

    NULLIFY(GEOMETRIC_PARAMETERS)

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
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
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
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
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
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,1, &
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
                          CALL MeshTopologyNodeCheckExists(MESH,1,user_node,NODE_EXISTS,global_node,ERR,ERROR,*999)
                          IF(.NOT.NODE_EXISTS) CYCLE
                          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,ERR,ERROR,*999)
                          !Default to version 1 of each node derivative
                          local_ny=GEOMETRIC_VARIABLE%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(local_node)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          DEFORMED_Z=GEOMETRIC_PARAMETERS(local_ny)*LAMBDA
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                            & user_node,ABS(TOP_NORMAL_XI),BOUNDARY_CONDITION_FIXED,DEFORMED_Z,ERR,ERROR,*999)
                        ENDIF
                      ENDDO
                      !Set all bottom nodes fixed in z plane
                      DO node_idx=1,SIZE(BOTTOM_SURFACE_NODES,1)
                        user_node=BOTTOM_SURFACE_NODES(node_idx)
                        !Need to check this node exists in the current domain
                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node,1,DOMAIN_NUMBER,ERR,ERROR,*999)
                        IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                          !Default to version 1 of each node derivative
                          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
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
                          CALL MeshTopologyNodeCheckExists(MESH,1,user_node,NODE_EXISTS,global_node,ERR,ERROR,*999)
                          IF(.NOT.NODE_EXISTS) CYCLE
                          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node,ERR,ERROR,*999)
                          !Default to version 1 of each node derivative
                          local_ny=GEOMETRIC_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(local_node)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          X(1)=GEOMETRIC_PARAMETERS(local_ny)
                            CALL MeshTopologyNodeCheckExists(MESH,1,user_node,NODE_EXISTS,global_node,ERR,ERROR,*999)
                            IF(.NOT.NODE_EXISTS) CYCLE
                            CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(NODES_MAPPING,global_node,NODE_EXISTS,local_node, &
                              & ERR,ERROR,*999)
                            !Default to version 1 of each node derivative
                            local_ny=GEOMETRIC_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(local_node)% &
                            & DERIVATIVES(1)%VERSIONS(1)
                          X(2)=GEOMETRIC_PARAMETERS(local_ny)
                          IF(ABS(X(1))<1E-7_DP) THEN
                            !Default to version 1 of each node derivative
                            CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                              & user_node,1,BOUNDARY_CONDITION_FIXED,0.0_DP,ERR,ERROR,*999)

                            X_FIXED=.TRUE.
                          ENDIF
                          IF(ABS(X(2))<1E-7_DP) THEN
                            !Default to version 1 of each node derivative
                            CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,1, &
                              & user_node,2,BOUNDARY_CONDITION_FIXED,0.0_DP,ERR,ERROR,*999)

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
                    CALL FLAG_ERROR("Generated mesh is not associated. For the Cylinder analytic solution, "// &
                      & "it must be available for automatic boundary condition assignment",ERR,ERROR,*999)
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
                                      !Default to version 1 of each node derivative
                                      local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                        & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                      X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                    ENDDO !dim_idx
                                    !Loop over the derivatives
                                    DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                      SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                      CASE(EQUATIONS_SET_FINITE_ELASTICITY_CYLINDER)
                                        !Cylinder inflation, extension, torsion
                                        SELECT CASE(variable_type)
                                        CASE(FIELD_U_VARIABLE_TYPE)
                                          SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
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
                                              DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                              & ERR,ERROR))//" is invalid."
                                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                          END SELECT
                                        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                          SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
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
                                              DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                                & ERR,ERROR))//" is invalid."
                                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                          END SELECT
                                        CASE DEFAULT
                                          LOCAL_ERROR="The variable type "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR)) &
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
                                        !Default to version 1 of each node derivative
                                        local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                          & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                          & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,DEFORMED_X(component_idx),ERR,ERROR,*999)
                                      ENDDO
                                      !Don't forget the pressure component
                                      user_node=DOMAIN_NODES%NODES(node_idx)%USER_NUMBER
                                      CALL MeshTopologyNodeCheckExists(MESH,DOMAIN_PRESSURE%MESH_COMPONENT_NUMBER,user_node, &
                                        & NODE_EXISTS,global_node,ERR,ERROR,*999)
                                      IF(NODE_EXISTS) THEN
                                        CALL DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,user_node, &
                                          & DOMAIN_PRESSURE%MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*999)
                                        IF(DOMAIN_NUMBER==MY_COMPUTATIONAL_NODE_NUMBER) THEN
                                          !\todo: test the domain node mappings pointer properly
                                          local_node=DOMAIN_PRESSURE%mappings%nodes%global_to_local_map(global_node)%local_number(1)
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(4)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                            & NODES(local_node)%DERIVATIVES(deriv_idx)%VERSIONS(1)
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
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Boundary conditions is not associated.",ERR,ERROR,*999)
            ENDIF
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
        IF (ABS(F2)<ABS(F).OR.ABS(F2)<ZERO_TOLERANCE) THEN    ! PASS
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
  
  !>Evaluates the spatial elasticity and stress tensor in Voigt form at a given Gauss point.
  SUBROUTINE FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,ELASTICITY_TENSOR,STRESS_TENSOR,DZDNU,Jznu,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: ELASTICITY_TENSOR(:,:) !< Rank 4 elasticity tensor in Voigt notation
    REAL(DP), INTENT(OUT) :: STRESS_TENSOR(:) !< Rank 2 stress tensor in Voigt notation
    REAL(DP), INTENT(IN) :: DZDNU(:,:) !< The deformation gradient
    REAL(DP), INTENT(IN) :: Jznu !< The Jacobian
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: PRESSURE_COMPONENT,i,j,dof_idx
    REAL(DP) :: P
    REAL(DP) :: TEMPTERM1,TEMPTERM2
    REAL(DP), PARAMETER :: ONETHIRD=1.0_DP/3.0_DP,TWOTHIRDS=2.0_DP/3.0_DP
    REAL(DP) :: I1,I2,I3            !Invariants, if needed
    REAL(DP), PARAMETER :: DEV_PROJ(6,6)=RESHAPE([TWOTHIRDS,-ONETHIRD,-ONETHIRD,0.0_DP,0.0_DP,0.0_DP, &
                                                  -ONETHIRD,TWOTHIRDS,-ONETHIRD,0.0_DP,0.0_DP,0.0_DP, &
                                                  -ONETHIRD,-ONETHIRD,TWOTHIRDS,0.0_DP,0.0_DP,0.0_DP, &
                                                  0.0_DP,0.0_DP,0.0_DP,1.0_DP,0.0_DP,0.0_DP, &
                                                  0.0_DP,0.0_DP,0.0_DP,0.0_DP,1.0_DP,0.0_DP, &
                                                  0.0_DP,0.0_DP,0.0_DP,0.0_DP,0.0_DP,1.0_DP],[6,6])
    REAL(DP), PARAMETER :: TWOTHIRDS_UNITY(6) = [TWOTHIRDS,TWOTHIRDS,TWOTHIRDS,0.0_DP,0.0_DP,0.0_DP] !Rank 2 unit tensor times 2/3 in Voigt form.
    REAL(DP), PARAMETER :: UNITY_DIAGONAL(6)=[1.0_DP,1.0_DP,1.0_DP,0.5_DP,0.5_DP,0.5_DP] !Diagonal of rank 4 unit tensor in Voigt form.
    REAL(DP), PARAMETER :: UNITY(6) = [1.0_DP,1.0_DP,1.0_DP,0.0_DP,0.0_DP,0.0_DP] !Rank 2 unit tensor in Voigt form.
    REAL(DP), POINTER :: C(:) !Parameters for constitutive laws
    REAL(DP) :: MOD_DZDNU(3,3),MOD_DZDNUT(3,3),AZL(3,3),AZU(3,3)
    REAL(DP) :: TRACE,TWOTHIRDS_TRACE
    REAL(DP) :: B(6),E(6),DQ_DE(6)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR",ERR,ERROR,*999)
    
    NULLIFY(FIELD_VARIABLE,C)

    !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
    !AZU - deformed contravariant tensor; I3 = det(C)
    !E = Green-Lagrange strain tensor = 0.5*(C-I)
    !P is the hydrostatic pressure

    ! From now on, we calculate with the modified deformation tensor (=F*Jznu**(-1/3))
    MOD_DZDNU=DZDNU*Jznu**(-ONETHIRD)
    CALL MATRIX_TRANSPOSE(MOD_DZDNU,MOD_DZDNUT,ERR,ERROR,*999)
    CALL MATRIX_PRODUCT(MOD_DZDNUT,MOD_DZDNU,AZL,ERR,ERROR,*999)

    C=>MATERIALS_INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)

    ELASTICITY_TENSOR=0.0_DP

    SELECT CASE(EQUATIONS_SET%SUBTYPE)
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+p/2*(I3-1)

      ! Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMPTERM1=-2.0_DP*C(2)
      TEMPTERM2=2.0_DP*(C(1)+I1*C(2))
      STRESS_TENSOR(1)=TEMPTERM1*AZL(1,1)+TEMPTERM2
      STRESS_TENSOR(2)=TEMPTERM1*AZL(2,2)+TEMPTERM2
      STRESS_TENSOR(3)=TEMPTERM1*AZL(3,3)+TEMPTERM2
      STRESS_TENSOR(4)=TEMPTERM1*AZL(2,1)
      STRESS_TENSOR(5)=TEMPTERM1*AZL(3,1)
      STRESS_TENSOR(6)=TEMPTERM1*AZL(3,2)
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        !the active stress is stored inside the independent field that has been set up in the user program.
        !for generality we could set up 3 components in independent field for 3 different active stress components
        !!!!! Be aware for modified DZDNU, check if this the right way to do it?
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO i=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(i)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,TEMPTERM1,ERR,ERROR,*999)
          STRESS_TENSOR(i)=STRESS_TENSOR(i)+TEMPTERM1
        ENDDO
      ENDIF

      ! Calculate isochoric fictitious material elasticity tensor (in Voigt form), without the factor Jznu**(-4.0_DP/3.0_DP), as
      ! this will be compensated for in the push-forward with the modified deformation gradient.
      TEMPTERM1=4.0_DP*C(2)
      TEMPTERM2=-2.0_DP*C(2)
      ELASTICITY_TENSOR(2,1)=TEMPTERM1
      ELASTICITY_TENSOR(3,1)=TEMPTERM1
      ELASTICITY_TENSOR(1,2)=TEMPTERM1
      ELASTICITY_TENSOR(3,2)=TEMPTERM1
      ELASTICITY_TENSOR(1,3)=TEMPTERM1
      ELASTICITY_TENSOR(2,3)=TEMPTERM1
      ELASTICITY_TENSOR(4,4)=TEMPTERM2
      ELASTICITY_TENSOR(5,5)=TEMPTERM2
      ELASTICITY_TENSOR(6,6)=TEMPTERM2

      ! Do push-forward of 2nd Piola tensor and the material elasticity tensor. 
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,MOD_DZDNU,Jznu,ERR,ERROR,*999)
      CALL FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR(ELASTICITY_TENSOR,MOD_DZDNU,Jznu,ERR,ERROR,*999)
      
      TRACE=SUM(STRESS_TENSOR(1:3))
      !Calculate isochoric Cauchy tensor (the deviatoric part).
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)-ONETHIRD*TRACE
      TWOTHIRDS_TRACE=TWOTHIRDS*TRACE

      DO i=1,6
        ELASTICITY_TENSOR(i,i)=ELASTICITY_TENSOR(i,i)+TWOTHIRDS_TRACE*UNITY_DIAGONAL(i)
      ENDDO
      ELASTICITY_TENSOR=MATMUL(DEV_PROJ,MATMUL(ELASTICITY_TENSOR,DEV_PROJ))
      DO j=1,6
        DO i=1,6
          ELASTICITY_TENSOR(i,j)=ELASTICITY_TENSOR(i,j)- &
            & (TWOTHIRDS_UNITY(i)*STRESS_TENSOR(j)+TWOTHIRDS_UNITY(j)*STRESS_TENSOR(i))
        ENDDO
      ENDDO
      ! Add volumetric parts.
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)+P
      ELASTICITY_TENSOR(1,1)=ELASTICITY_TENSOR(1,1)-P
      ELASTICITY_TENSOR(2,2)=ELASTICITY_TENSOR(2,2)-P
      ELASTICITY_TENSOR(3,3)=ELASTICITY_TENSOR(3,3)-P
      ELASTICITY_TENSOR(4,4)=ELASTICITY_TENSOR(4,4)-P
      ELASTICITY_TENSOR(5,5)=ELASTICITY_TENSOR(5,5)-P
      ELASTICITY_TENSOR(6,6)=ELASTICITY_TENSOR(6,6)-P
      ELASTICITY_TENSOR(2,1)=ELASTICITY_TENSOR(2,1)+P
      ELASTICITY_TENSOR(3,1)=ELASTICITY_TENSOR(3,1)+P
      ELASTICITY_TENSOR(1,2)=ELASTICITY_TENSOR(1,2)+P
      ELASTICITY_TENSOR(3,2)=ELASTICITY_TENSOR(3,2)+P
      ELASTICITY_TENSOR(1,3)=ELASTICITY_TENSOR(1,3)+P
      ELASTICITY_TENSOR(2,3)=ELASTICITY_TENSOR(2,3)+P

    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      B=[2.0_DP*C(2),2.0_DP*C(3),2.0_DP*C(3),C(4),C(4),C(3)] ![2*b_f,2*b_t,2*b_t,b_ft,b_ft,b_t]
      E=[0.5_DP*(AZL(1,1)-1.0_DP),0.5_DP*(AZL(2,2)-1.0_DP),0.5_DP*(AZL(3,3)-1.0_DP),AZL(2,1),AZL(3,1),AZL(3,2)] !(Modified) strain tensor in Voigt form.
      DQ_DE=B*E
      TEMPTERM1=0.5_DP*C(1)*EXP(0.5_DP*DOT_PRODUCT(E,DQ_DE))
      ! Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      STRESS_TENSOR=TEMPTERM1*DQ_DE
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        !the active stress is stored inside the independent field that has been set up in the user program.
        !for generality we could set up 3 components in independent field for 3 different active stress components
        !!!!! Be aware for modified DZDNU, check if this the right way to do it?
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO i=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(i)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,TEMPTERM1,ERR,ERROR,*999)
          STRESS_TENSOR(i)=STRESS_TENSOR(i)+TEMPTERM1
        ENDDO
      ENDIF

      ! Calculate isochoric fictitious material elasticity tensor (in Voigt form), without the factor Jznu**(-4.0_DP/3.0_DP), as
      ! this will be compensated for in the push-forward with the modified deformation gradient.
      ! First calculate lower part of 6X6 matrix
      DO j=1,6
        DO i=j,6
          ELASTICITY_TENSOR(i,j)=TEMPTERM1*DQ_DE(i)*DQ_DE(j)
        ENDDO
      ENDDO
      DO i=1,6
        ELASTICITY_TENSOR(i,i)=ELASTICITY_TENSOR(i,i)+TEMPTERM1*B(i)
      ENDDO
      ! Then calculate upper part.
      DO j=2,6
        DO i=1,j-1
          ELASTICITY_TENSOR(i,j)=ELASTICITY_TENSOR(j,i)
        ENDDO
      ENDDO

      !Do push-forward of 2nd Piola tensor and the material elasticity tensor. 
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,MOD_DZDNU,Jznu,ERR,ERROR,*999)
      CALL FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR(ELASTICITY_TENSOR,MOD_DZDNU,Jznu,ERR,ERROR,*999)
      
      TRACE=SUM(STRESS_TENSOR(1:3))
      !Calculate isochoric Cauchy tensor (the deviatoric part) and volumetric part (hydrostatic pressure).
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)-ONETHIRD*TRACE

      TWOTHIRDS_TRACE=TWOTHIRDS*TRACE
      DO i=1,6
        ELASTICITY_TENSOR(i,i)=ELASTICITY_TENSOR(i,i)+TWOTHIRDS_TRACE*UNITY_DIAGONAL(i)
      ENDDO
      ELASTICITY_TENSOR=MATMUL(DEV_PROJ,MATMUL(ELASTICITY_TENSOR,DEV_PROJ))
      DO j=1,6
        DO i=1,6
          ELASTICITY_TENSOR(i,j)=ELASTICITY_TENSOR(i,j)- &
            & (TWOTHIRDS_UNITY(i)*STRESS_TENSOR(j)+TWOTHIRDS_UNITY(j)*STRESS_TENSOR(i))
        ENDDO
      ENDDO
      !Add volumetric parts.
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)+P
      ELASTICITY_TENSOR(1,1)=ELASTICITY_TENSOR(1,1)-P
      ELASTICITY_TENSOR(2,2)=ELASTICITY_TENSOR(2,2)-P
      ELASTICITY_TENSOR(3,3)=ELASTICITY_TENSOR(3,3)-P
      ELASTICITY_TENSOR(4,4)=ELASTICITY_TENSOR(4,4)-P
      ELASTICITY_TENSOR(5,5)=ELASTICITY_TENSOR(5,5)-P
      ELASTICITY_TENSOR(6,6)=ELASTICITY_TENSOR(6,6)-P
      ELASTICITY_TENSOR(2,1)=ELASTICITY_TENSOR(2,1)+P
      ELASTICITY_TENSOR(3,1)=ELASTICITY_TENSOR(3,1)+P
      ELASTICITY_TENSOR(1,2)=ELASTICITY_TENSOR(1,2)+P
      ELASTICITY_TENSOR(3,2)=ELASTICITY_TENSOR(3,2)+P
      ELASTICITY_TENSOR(1,3)=ELASTICITY_TENSOR(1,3)+P
      ELASTICITY_TENSOR(2,3)=ELASTICITY_TENSOR(2,3)+P
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for a finite elasticity equation type of an elasticity equation set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    CALL EXITS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix for the given element number for a finite elasticity class finite element equation set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: FIELD_VAR_TYPE,ng,nh,ns,nhs,ni,mh,ms,mhs,mi,oh
    INTEGER(INTG) :: PRESSURE_COMPONENT
    INTEGER(INTG) :: SUM_ELEMENT_PARAMETERS,TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI 
    INTEGER(INTG) :: ELEMENT_BASE_DOF_INDEX(4)
    INTEGER(INTG), PARAMETER :: OFF_DIAG_COMP(3)=[0,1,3],OFF_DIAG_DEP_VAR1(3)=[1,1,2],OFF_DIAG_DEP_VAR2(3)=[2,3,3]
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_ELEMENT_PARAMETERS(4)
    REAL(DP) :: DZDNU(3,3),CAUCHY_TENSOR(3,3),DNUDXI(3,3),DXIDNU(3,3)
    REAL(DP) :: JGW_SUB_MAT(3,3) 
    REAL(DP) :: TEMPVEC(3)
    REAL(DP) :: STRESS_TENSOR(6),ELASTICITY_TENSOR(6,6)
    REAL(DP) :: DPHIDZ(3,64,3)
    REAL(DP) :: JGW_DPHINS_DZ,PHIMS
    REAL(DP) :: Jznu,JGW,SUM1
    TYPE(QUADRATURE_SCHEME_PTR_TYPE) :: QUADRATURE_SCHEMES(4)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERP_POINT,FIBRE_INTERP_POINT, &
      & MATERIALS_INTERP_POINT,DEPENDENT_INTERP_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERP_POINT_METRICS, &
      & DEPENDENT_INTERP_POINT_METRICS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,FIBRE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME

    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR
        IF(JACOBIAN_MATRIX%UPDATE_JACOBIAN) THEN
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          FIBRE_FIELD=>EQUATIONS%INTERPOLATION%FIBRE_FIELD

          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          
          NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
          NUMBER_OF_XI=DEPENDENT_BASIS%NUMBER_OF_XI

          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          NONLINEAR_MAPPING=>EQUATIONS_MAPPING%NONLINEAR_MAPPING
          
          FIELD_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE

          PRESSURE_COMPONENT=FIELD_VARIABLE%NUMBER_OF_COMPONENTS

          BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING% &
            & RHS_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
          TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE)+ &
            & BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
        
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          IF(ASSOCIATED(FIBRE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
          END IF
          
          !Point interpolation pointer
          GEOMETRIC_INTERP_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          GEOMETRIC_INTERP_POINT_METRICS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIBRE_FIELD)) THEN
            FIBRE_INTERP_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          END IF
          MATERIALS_INTERP_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_INTERP_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
          DEPENDENT_INTERP_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR
          
          SUM_ELEMENT_PARAMETERS=0
          !Loop over geometric dependent basis functions.
          DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEMES(nh)%PTR=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            IF(FIELD_VARIABLE%COMPONENTS(nh)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              NUMBER_OF_ELEMENT_PARAMETERS(nh)=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            ELSEIF(FIELD_VARIABLE%COMPONENTS(nh)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
              NUMBER_OF_ELEMENT_PARAMETERS(nh)=1
            ENDIF
            ELEMENT_BASE_DOF_INDEX(nh)=SUM_ELEMENT_PARAMETERS
            SUM_ELEMENT_PARAMETERS=SUM_ELEMENT_PARAMETERS+NUMBER_OF_ELEMENT_PARAMETERS(nh)
          ENDDO !nh

          !Loop over all Gauss points 
          DO ng=1,DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE, &
              & GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE, &
              & DEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & MATERIALS_INTERP_POINT,ERR,ERROR,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & FIBRE_INTERP_POINT,ERR,ERROR,*999)
            ENDIF

            Jznu=DEPENDENT_INTERP_POINT_METRICS%JACOBIAN/GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN 
            JGW=DEPENDENT_INTERP_POINT_METRICS%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            
            !Loop over geometric dependent basis functions.
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                !Loop over derivative directions.
                DO mh=1,NUMBER_OF_DIMENSIONS
                  SUM1=0.0_DP
                  DO ni=1,NUMBER_OF_XI
                    SUM1=SUM1+QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                      & DEPENDENT_INTERP_POINT_METRICS%DXI_DX(ni,mh)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh

            CALL CoordinateMaterialSystemCalculate(GEOMETRIC_INTERP_POINT_METRICS,FIBRE_INTERP_POINT, &
              & DNUDXI,DXIDNU,ERR,ERROR,*999)
            !dZ/dNu = dZ/dXi * dXi/dNu  (deformation gradient tensor, F)
            CALL MATRIX_PRODUCT(DEPENDENT_INTERP_POINT_METRICS%DX_DXI,DXIDNU,DZDNU,ERR,ERROR,*999)

            CALL FINITE_ELASTICITY_GAUSS_ELASTICITY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERP_POINT, &
              & MATERIALS_INTERP_POINT,ELASTICITY_TENSOR,STRESS_TENSOR,DZDNU,Jznu,ELEMENT_NUMBER,ng,ERR,ERROR,*999)

            ! Convert from Voigt form to tensor form.
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO mh=1,NUMBER_OF_DIMENSIONS
                CAUCHY_TENSOR(mh,nh)=STRESS_TENSOR(TENSOR_TO_VOIGT(mh,nh))
              ENDDO
            ENDDO
           
            !First: loop over mh=nh
            !Loop over element columns belonging to geometric dependent variables
            nhs=0
            DO nh=1,NUMBER_OF_DIMENSIONS
              JGW_SUB_MAT=JGW*(ELASTICITY_TENSOR(TENSOR_TO_VOIGT(:,nh),TENSOR_TO_VOIGT(:,nh))+CAUCHY_TENSOR)
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                TEMPVEC=MATMUL(JGW_SUB_MAT,DPHIDZ(:,ns,nh))
                nhs=nhs+1
                mhs=nhs-1
                !Loop over element rows belonging to geometric dependent variables
                DO ms=ns,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  mhs=mhs+1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                    & DOT_PRODUCT(DPHIDZ(:,ms,nh),TEMPVEC)
                ENDDO !ms    
              ENDDO !ns
            ENDDO !nh
         
            !Second: loop over mh>nh
            !Loop over element columns belonging to geometric dependent variables
            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              nhs=ELEMENT_BASE_DOF_INDEX(nh)
              JGW_SUB_MAT=JGW*ELASTICITY_TENSOR(TENSOR_TO_VOIGT(:,mh),TENSOR_TO_VOIGT(:,nh))
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                !Loop over element rows belonging to geometric dependent variables
                TEMPVEC=MATMUL(JGW_SUB_MAT,DPHIDZ(:,ns,nh))
                nhs=nhs+1
                mhs=ELEMENT_BASE_DOF_INDEX(mh)
                DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(mh)
                  mhs=mhs+1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                    & DOT_PRODUCT(DPHIDZ(:,ms,mh),TEMPVEC)
                ENDDO !ms    
              ENDDO !ns
            ENDDO

            !Third: loop over all nh and pressure component
            nhs=0
            IF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  JGW_DPHINS_DZ=JGW*DPHIDZ(nh,ns,nh)
                  nhs=nhs+1
                 !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)
                  DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(PRESSURE_COMPONENT)
                    mhs=mhs+1
                    PHIMS=QUADRATURE_SCHEMES(PRESSURE_COMPONENT)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                      & JGW_DPHINS_DZ*PHIMS
                  ENDDO !ms    
                ENDDO !ns
              ENDDO !nh
            ELSEIF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  JGW_DPHINS_DZ=JGW*DPHIDZ(nh,ns,nh)
                  nhs=nhs+1
                  !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)+1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                    & JGW_DPHINS_DZ
                ENDDO !ns
              ENDDO !nh
            ENDIF
            ! No loop over element columns and rows belonging both to hydrostatic pressure because it is zero.
          ENDDO !ng

          !If symmetric pressure Jacobian uncomment this.
          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
!          IF(DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
!            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
!            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
!          ENDIF

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            ! Following function is necessary, otherwise wrong face scale factors from function call to surface pressure residual are
            ! used.
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER, &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR,ERR,ERROR,*999) 
            nhs=0          
            ! Loop over element columns
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                nhs=nhs+1
                mhs=nhs-1
                ! Loop over element rows
                DO ms=ns,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  mhs=mhs+1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,nh)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                ENDDO !ms
              ENDDO !ns
            ENDDO !nh
            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              nhs=ELEMENT_BASE_DOF_INDEX(nh)
              DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                nhs=nhs+1
                mhs=ELEMENT_BASE_DOF_INDEX(mh)
                !Loop over element rows belonging to geometric dependent variables
                DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(mh)
                  mhs=mhs+1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                ENDDO !ms    
              ENDDO !ns
            ENDDO

            nhs=0
            IF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  nhs=nhs+1
                  !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)
                  DO ms=1,NUMBER_OF_ELEMENT_PARAMETERS(PRESSURE_COMPONENT)
                    mhs=mhs+1
                    JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)* &
                      & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR% &
                      & SCALE_FACTORS(ms,PRESSURE_COMPONENT)* &
                      & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                  ENDDO !ms    
                ENDDO !ns
              ENDDO !nh
            ELSEIF(FIELD_VARIABLE%COMPONENTS(PRESSURE_COMPONENT)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              !Loop over element rows belonging to geometric dependent variables
              DO nh=1,NUMBER_OF_DIMENSIONS
                DO ns=1,NUMBER_OF_ELEMENT_PARAMETERS(nh)
                  nhs=nhs+1
                  !Loop over element rows belonging to hydrostatic pressure
                  mhs=ELEMENT_BASE_DOF_INDEX(PRESSURE_COMPONENT)+1
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)* &
                    & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                ENDDO !ns
              ENDDO !nh
            ENDIF
          ENDIF

          DO nhs=2,SIZE(JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX,2)
            DO mhs=1,nhs-1
              JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(nhs,mhs)
            ENDDO !mhs
          ENDDO !nhs

          !If unsymmetric pressure Jacobian uncomment this.
!         !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Push-forward the rank 4 elasticity tensor.
  SUBROUTINE FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR(ELASTICITY_TENSOR,DZDNU,Jznu,ERR,ERROR,*)
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: ELASTICITY_TENSOR(6,6)
    REAL(DP), INTENT(IN) :: DZDNU(3,3)
    REAL(DP), INTENT(IN) :: Jznu
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j
    INTEGER(INTG), PARAMETER :: idx1(6)=[1,2,3,1,1,2],idx2(6)=[1,2,3,2,3,3]
    REAL(DP) :: t(6,6) 

    CALL ENTERS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR",ERR,ERROR,*999)

    DO j=1,3
      DO i=1,6
        t(i,j)=DZDNU(idx1(i),idx1(j))*DZDNU(idx2(i),idx2(j))
      ENDDO
    END DO
    DO j=4,6
      DO i=1,6
        t(i,j)=DZDNU(idx1(i),idx1(j))*DZDNU(idx2(i),idx2(j))+DZDNU(idx1(i),idx2(j))*DZDNU(idx2(i),idx1(j))
      ENDDO
    END DO

    ELASTICITY_TENSOR=MATMUL(MATMUL(t,ELASTICITY_TENSOR),TRANSPOSE(t))/Jznu

    CALL EXITS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PUSH_ELASTICITY_TENSOR

  !
  !================================================================================================================================
  !

  !>Push-forward the rank 2 Piola stress tensor.
  SUBROUTINE FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,DZDNU,Jznu,ERR,ERROR,*)
    
    !Argument variables
    REAL(DP), INTENT(INOUT) :: STRESS_TENSOR(6)
    REAL(DP), INTENT(IN) :: DZDNU(3,3)
    REAL(DP), INTENT(IN) :: Jznu
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j
    INTEGER(INTG), PARAMETER :: idx1(6)=[1,2,3,1,1,2],idx2(6)=[1,2,3,2,3,3]
    REAL(DP) :: t(6,6) 

    CALL ENTERS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR",ERR,ERROR,*999)

    DO j=1,3
      DO i=1,6
        t(i,j)=DZDNU(idx1(i),idx1(j))*DZDNU(idx2(i),idx2(j))
      ENDDO
    END DO
    DO j=4,6
      DO i=1,6
        t(i,j)=DZDNU(idx1(i),idx1(j))*DZDNU(idx2(i),idx2(j))+DZDNU(idx1(i),idx2(j))*DZDNU(idx2(i),idx1(j))
      ENDDO
    END DO

    STRESS_TENSOR=MATMUL(t,STRESS_TENSOR)/Jznu

    CALL EXITS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PUSH_STRESS_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PUSH_STRESS_TENSOR

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
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,EQUATIONS_SET_FIELD,SOURCE_FIELD
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME,COMPONENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS, &
      & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATION_PARAMETERS,DARCY_MATERIALS_INTERPOLATION_PARAMETERS, &
      & DENSITY_INTERPOLATION_PARAMETERS,INDEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT,SOURCE_INTERPOLATED_POINT, &
      & DENSITY_INTERPOLATED_POINT,INDEPENDENT_INTERPOLATED_POINT,DARCY_MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1,GEOMETRIC_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENT_MAPPING
    LOGICAL :: DARCY_DENSITY,DARCY_DEPENDENT
    INTEGER(INTG) :: component_idx,component_idx2,parameter_idx,gauss_idx,element_dof_idx,FIELD_VAR_TYPE,DARCY_FIELD_VAR_TYPE
    INTEGER(INTG) :: idx,imatrix,Ncompartments
    INTEGER(INTG) :: NDOFS,mh,ms,mhs,mi,nh,ns
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI,HYDROSTATIC_PRESSURE_COMPONENT
    INTEGER(INTG) :: NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
    INTEGER(INTG) :: DEPENDENT_COMPONENT_INTERPOLATION_TYPE
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS       
    INTEGER(INTG) :: MESH_COMPONENT_1,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS
    INTEGER(INTG) :: var1 ! Variable number corresponding to 'U' in single physics case
    INTEGER(INTG) :: var2 ! Variable number corresponding to 'DELUDLEN' in single physics case
    INTEGER(INTG) :: numberOfXDimensions,numberOfZDimensions,numberOfXiDimensions
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    REAL(DP) :: DZDNU(3,3),DZDNUT(3,3),AZL(3,3),AZU(3,3),I3,P,PIOLA_TENSOR(3,3),TEMP(3,3)
    REAL(DP) :: CAUCHY_TENSOR(3,3),JGW_CAUCHY_TENSOR(3,3),STRESS_TENSOR(6)
    REAL(DP) :: DNUDXI(3,3),DXIDNU(3,3)
    REAL(DP) :: DFDZ(64,3,3) !temporary until a proper alternative is found
    REAL(DP) :: DPHIDZ(3,64,3) !temporary until a proper alternative is found
    REAL(DP) :: GAUSS_WEIGHT,Jznu,Jxxi,JGW
    REAL(DP) :: SUM1,TEMPTERM1
    REAL(DP) :: THICKNESS ! for elastic membrane
    REAL(DP) :: DARCY_MASS_INCREASE,DARCY_VOL_INCREASE,DARCY_RHO_0_F,DENSITY  !coupling with Darcy model
    REAL(DP) :: Mfact, bfact, p0fact


    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS,BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DEPENDENT_BASIS,COMPONENT_BASIS)
    NULLIFY(EQUATIONS,EQUATIONS_MAPPING,EQUATIONS_MATRICES,NONLINEAR_MATRICES,RHS_VECTOR)
    NULLIFY(DEPENDENT_FIELD,FIBRE_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD,INDEPENDENT_FIELD)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME,COMPONENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS,SOURCE_INTERPOLATION_PARAMETERS)
    NULLIFY(MATERIALS_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(INDEPENDENT_INTERPOLATION_PARAMETERS,DARCY_MATERIALS_INTERPOLATION_PARAMETERS)
    NULLIFY(DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,DENSITY_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,SOURCE_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)
    NULLIFY(MATERIALS_INTERPOLATED_POINT,DEPENDENT_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(DENSITY_INTERPOLATED_POINT,INDEPENDENT_INTERPOLATED_POINT)
    NULLIFY(DEPENDENT_BASIS_1)
    NULLIFY(DECOMPOSITION)
    NULLIFY(EQUATIONS_SET_FIELD_DATA)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN 
        !Which variables are we working with - find the variable pair used for this equations set
        !\todo: put in checks for all the objects/mappings below (do we want to do this for every element?)
        var1=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%VARIABLE_NUMBER ! number for 'U'
        var2=EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE%VARIABLE_NUMBER ! number for 'DELUDELN'

        !Grab pointers: matrices, fields, decomposition, basis
        !\todo: see if we can separate this residual evaluation from the pressure boundary conditions somehow
        !so that the equations set doesn't need to maintain a pointer to the boundary conditions
        BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING% &
          & RHS_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
        TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS=BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE)+ &
          & BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)

        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
        EQUATIONS_MAPPING =>EQUATIONS%EQUATIONS_MAPPING

        FIBRE_FIELD      =>EQUATIONS%INTERPOLATION%FIBRE_FIELD
        GEOMETRIC_FIELD  =>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
        MATERIALS_FIELD  =>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
        DEPENDENT_FIELD  =>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
        SOURCE_FIELD     =>EQUATIONS%INTERPOLATION%SOURCE_FIELD
        INDEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%INDEPENDENT_FIELD

        DECOMPOSITION    =>DEPENDENT_FIELD%DECOMPOSITION
        MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER

        DOMAIN_ELEMENT_MAPPING=>DECOMPOSITION%DOMAIN(1)%PTR%MAPPINGS%ELEMENTS

        DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS       
        DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
        DEPENDENT_NUMBER_OF_COMPONENTS=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
        GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

        NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        NUMBER_OF_XI=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI

        !Initialise tensors and matrices
        CALL IdentityMatrix(DZDNU,err,error,*999)
        CALL IdentityMatrix(PIOLA_TENSOR,err,error,*999)
        CALL IdentityMatrix(CAUCHY_TENSOR,err,error,*999)
        DFDZ=0.0_DP ! (parameter_idx,component_idx)

        !Set flags for coupled finite elasticity and Darcy problems
        !Check if we need Darcy materials field for Density
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE) THEN
          DARCY_DENSITY=.TRUE.
        ELSE
          DARCY_DENSITY=.FALSE.
        ENDIF
        !Check if we need Darcy dependent field
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE .OR. &
          & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE) THEN
          DARCY_DEPENDENT=.TRUE.
        ELSE
          DARCY_DEPENDENT=.FALSE.
        ENDIF

        !Grab interpolation parameters
        FIELD_VARIABLE=>EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
        FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR
        GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        ENDIF
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          MATERIALS_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        ENDIF
        IF(DARCY_DEPENDENT) THEN
          DARCY_DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR
        ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
          INDEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        ENDIF

        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        END IF
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & MATERIALS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        ENDIF
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        IF(DARCY_DEPENDENT) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & INDEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        ENDIF

        !Point interpolation pointer
        GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        GEOMETRIC_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        END IF
        IF(ASSOCIATED(MATERIALS_FIELD)) THEN
          MATERIALS_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        ENDIF
        DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
        DEPENDENT_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR
        IF(DARCY_DEPENDENT) THEN
          DARCY_DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR
        ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
          INDEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%INDEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        ENDIF

        !SELECT: Compressible or incompressible cases, or poro multicompartment
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
            & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) ! 4 dependent components
          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
            
            !Loop over geometric dependent basis functions.
            DO nh=1,NUMBER_OF_DIMENSIONS
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                !Loop over derivative directions.
                DO mh=1,NUMBER_OF_DIMENSIONS
                  SUM1=0.0_DP
                  DO mi=1,NUMBER_OF_XI
                    SUM1=SUM1+DEPENDENT_INTERPOLATED_POINT_METRICS%DXI_DX(mi,mh)* &
                      & COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(mi),gauss_idx)
                  ENDDO !mi
                  DPHIDZ(mh,ns,nh)=SUM1
                ENDDO !mh
              ENDDO !ns
            ENDDO !nh
            
            CALL CoordinateMaterialSystemCalculate(GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT, &
              & DNUDXI,DXIDNU,ERR,ERROR,*999)
            !dZ/dNu = dZ/dXi * dXi/dNu  (deformation gradient tensor, F)
            CALL MATRIX_PRODUCT(DEPENDENT_INTERPOLATED_POINT_METRICS%DX_DXI,DXIDNU,DZDNU,ERR,ERROR,*999)

            Jznu=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN/GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN
            JGW=DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

           !Calculate the Cauchy stress tensor at the gauss point.
           CALL FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
             & MATERIALS_INTERPOLATED_POINT,STRESS_TENSOR,DZDNU,Jznu,ELEMENT_NUMBER,gauss_idx,ERR,ERROR,*999)
            
            ! Convert from Voigt form to tensor form and multiply with Jacobian and Gauss weight.
            DO nh=1,NUMBER_OF_DIMENSIONS
              DO mh=1,NUMBER_OF_DIMENSIONS
                JGW_CAUCHY_TENSOR(mh,nh)=JGW*STRESS_TENSOR(TENSOR_TO_VOIGT(mh,nh))
              ENDDO
            ENDDO

            !Now add up the residual terms
            mhs=0
            DO mh=1,NUMBER_OF_DIMENSIONS
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)+ &
                  & DOT_PRODUCT(DPHIDZ(:,ms,mh),JGW_CAUCHY_TENSOR(:,mh))
              ENDDO !ms
            ENDDO !mh
            
            JGW=GEOMETRIC_INTERPOLATED_POINT_METRICS%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

            !Hydrostatic pressure component
            MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            COMPONENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            TEMPTERM1=JGW*(Jznu-1.0_DP)
            IF(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1 
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)+ &
                  & TEMPTERM1*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gauss_idx)
              ENDDO
            ELSEIF(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
              mhs=mhs+1
              NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)+TEMPTERM1
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            ! Following function is necessary, otherwise wrong face scale factors from function call to surface pressure residual are
            ! used.
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER, &
              & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999) 
            mhs=0          
            DO mh=1,NUMBER_OF_DIMENSIONS
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              !Loop over residual vector
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)* &
                & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
            IF(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1 
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)=NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(mhs)* &
                & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO
            ENDIF
          ENDIF

        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE,EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE) ! 4 dependent components

          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
              !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
            IF(DARCY_DEPENDENT) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & DARCY_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            ELSE IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & INDEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            ENDIF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticityGaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,Jxxi,ERR,ERROR,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ELEMENT_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,ERR,ERROR,*999)
            ENDIF

            !Calculate Sigma=1/Jznu.FTF', the Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
              & INDEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ELEMENT_NUMBER,gauss_idx,ERR,ERROR,*999)


            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
              !Parameters settings for coupled elasticity Darcy INRIA model:
              CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,ERR,ERROR,*999)
              DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV) 
              DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F
            ENDIF

            !Calculate dPhi/dZ at the gauss point, Phi is the basis function
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,ELEMENT_NUMBER,gauss_idx,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DFDZ,ERR,ERROR,*999)

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(NUMBER_OF_DIMENSIONS == 3) THEN
                THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
              ENDIF
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
                      & GAUSS_WEIGHT*Jxxi*Jznu*THICKNESS*CAUCHY_TENSOR(component_idx,component_idx2)* &
                      & DFDZ(parameter_idx,component_idx2,component_idx)
                  ENDDO ! component_idx2 (inner component index)
                ENDDO ! parameter_idx (residual vector loop)
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FLAG_ERROR("Finite elasticity with element based interpolation is not implemented.",ERR,ERROR,*999)
              ENDIF
            ENDDO ! component_idx

            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
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
                      & GAUSS_WEIGHT*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                  ELSE
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                      & GAUSS_WEIGHT*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP)
                  ENDIF
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHT*Jxxi* &
                    & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                ELSE
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHT*Jxxi* &
                    & (Jznu-1.0_DP)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF

        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE)

          !Loop over gauss points and add residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
            !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            ENDIF
            IF(ASSOCIATED(MATERIALS_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
            ENDIF
            
            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticityGaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,Jxxi,ERR,ERROR,*999)

            CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
            CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)

            HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
            P=DEPENDENT_INTERPOLATED_POINT%VALUES(HYDROSTATIC_PRESSURE_COMPONENT,1)

            CALL INVERT(AZL,AZU,I3,ERR,ERROR,*999)
            Jznu=I3**0.5_DP

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ELEMENT_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,ERR,ERROR,*999)
            ENDIF

            !get the stress field!!!
            IF(NUMBER_OF_DIMENSIONS==3) THEN
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,1,PIOLA_TENSOR(1,1),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,2,PIOLA_TENSOR(1,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,3,PIOLA_TENSOR(1,3),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,4,PIOLA_TENSOR(2,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,5,PIOLA_TENSOR(2,3),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,6,PIOLA_TENSOR(3,3),ERR,ERROR,*999)
              !CellML computes the deviatoric stress. Add the volumetric component!
              PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+2.0_DP*P*AZU(1,1)
              PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+2.0_DP*P*AZU(2,2)
              PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+2.0_DP*P*AZU(3,3)
              PIOLA_TENSOR(1,2)=PIOLA_TENSOR(1,2)+2.0_DP*P*AZU(1,2)
              PIOLA_TENSOR(1,3)=PIOLA_TENSOR(1,3)+2.0_DP*P*AZU(1,3)
              PIOLA_TENSOR(2,3)=PIOLA_TENSOR(2,3)+2.0_DP*P*AZU(2,3)
              PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
              PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
              PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
            ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,1,PIOLA_TENSOR(1,1),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,2,PIOLA_TENSOR(1,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ELEMENT_NUMBER,gauss_idx,3,PIOLA_TENSOR(2,2),ERR,ERROR,*999)
              !CellML computes the deviatoric stress. Add the volumetric component!
              PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+2.0_DP*P*AZU(1,1)
              PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+2.0_DP*P*AZU(2,2)
              PIOLA_TENSOR(1,2)=PIOLA_TENSOR(1,2)+2.0_DP*P*AZU(1,2)
              PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
            ELSE
              CALL FLAG_ERROR("Only 2 and 3 dimensional problems are implemented.",ERR,ERROR,*999)
            ENDIF

            !Compute the CAUCHY stress tensor
            CALL MATRIX_PRODUCT(DZDNU,PIOLA_TENSOR,TEMP,ERR,ERROR,*999)
            CALL MATRIX_PRODUCT(TEMP,DZDNUT,CAUCHY_TENSOR,ERR,ERROR,*999)
            CAUCHY_TENSOR=CAUCHY_TENSOR/Jznu

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,PIOLA_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',ERR,ERROR,*999)

              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,CAUCHY_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    CAUCHY_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',ERR,ERROR,*999)
            ENDIF

            !Calculate dPhi/dZ at the gauss point, Phi is the basis function
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,ELEMENT_NUMBER,gauss_idx,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DFDZ,ERR,ERROR,*999)

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            !!TODO Maybe have the thickness as a component in the equations set field.
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(NUMBER_OF_DIMENSIONS == 3) THEN
                THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
              ENDIF
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
                      & GAUSS_WEIGHT*Jxxi*Jznu*THICKNESS*CAUCHY_TENSOR(component_idx,component_idx2)* &
                      & DFDZ(parameter_idx,component_idx2,component_idx)
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
                      & GAUSS_WEIGHT*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                  ELSE
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                      & GAUSS_WEIGHT*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP)
                  ENDIF
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHT*Jxxi* &
                    & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                ELSE
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHT*Jxxi* &
                    & (Jznu-1.0_DP)
                ENDIF
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF

        ! ---------------------------------------------------------------
        CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
        !keep the multi-compartment case separate for the time being until the formulation has been finalised, then perhaps
        !integrate within the single compartment case
          !Loop over gauss points and add residuals
          EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)

          Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)

          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
              !Interpolate dependent, geometric, fibre and materials fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
             IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ELEMENT_NUMBER,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,ERR,ERROR,*999)
            ENDIF
            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,PIOLA_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',ERR,ERROR,*999)

              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,CAUCHY_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    CAUCHY_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',ERR,ERROR,*999)
            ENDIF

              !Parameters settings for coupled elasticity Darcy INRIA model:
              CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,ERR,ERROR,*999)

              DARCY_MASS_INCREASE = 0.0_DP
              DO imatrix=1,Ncompartments
                DARCY_FIELD_VAR_TYPE=FIELD_V_VARIABLE_TYPE+FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(imatrix-1)
                DARCY_DEPENDENT_INTERPOLATION_PARAMETERS=>&
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(DARCY_FIELD_VAR_TYPE)%PTR

                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
                  & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)

                DARCY_DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(DARCY_FIELD_VAR_TYPE)%PTR
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                  & DARCY_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)               

                DARCY_MASS_INCREASE = DARCY_MASS_INCREASE + DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(4,NO_PART_DERIV)
              ENDDO             

              DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticityGaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,Jxxi,ERR,ERROR,*999)

            !Calculate Sigma=1/Jznu.FTF', the Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
              & INDEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ELEMENT_NUMBER,gauss_idx,ERR,ERROR,*999)

            !Calculate dPhi/dZ at the gauss point, Phi is the basis function
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,ELEMENT_NUMBER,gauss_idx,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DFDZ,ERR,ERROR,*999)

            !For membrane theory in 3D space, the final equation is multiplied by thickness. Default to unit thickness if equation set subtype is not membrane
            THICKNESS = 1.0_DP
            IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              IF(NUMBER_OF_DIMENSIONS == 3) THEN
                THICKNESS = MATERIALS_INTERPOLATED_POINT%VALUES(MATERIALS_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS% &
                  & FIELD_VARIABLE%NUMBER_OF_COMPONENTS,1)
              ENDIF
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
                      & GAUSS_WEIGHT*Jxxi*Jznu*THICKNESS*CAUCHY_TENSOR(component_idx,component_idx2)* &
                      & DFDZ(parameter_idx,component_idx2,component_idx)
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
                    NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                      & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                      & GAUSS_WEIGHT*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)* &
                      & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHT*Jxxi* &
                    & (Jznu-1.0_DP-DARCY_VOL_INCREASE)
              ENDIF
            ENDIF
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    ! 
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF

        CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, & 
            & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
            & EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
          !compressible problem (no pressure component)

          !Loop over gauss points and add up residuals
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)

            !Interpolate fields at the gauss points
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
               CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            END IF
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
            IF(DARCY_DEPENDENT) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & DARCY_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999) ! 'FIRST_PART_DERIV' required ???
            ENDIF

            !Calculate F=dZ/dNU at the gauss point
            CALL FiniteElasticityGaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,Jxxi,ERR,ERROR,*999)

            !Calculate Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
              & INDEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ELEMENT_NUMBER,gauss_idx,ERR,ERROR,*999)

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
                    & GAUSS_WEIGHT*Jxxi*Jznu*(CAUCHY_TENSOR(component_idx,1)*DFDZ(parameter_idx,1,component_idx)+ &
                    & CAUCHY_TENSOR(component_idx,2)*DFDZ(parameter_idx,2,component_idx)+ &
                    & CAUCHY_TENSOR(component_idx,3)*DFDZ(parameter_idx,3,component_idx)) 
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN
                !Will probably never be used
                CALL FLAG_ERROR("Finite elasticity with element based interpolation is not implemented.",ERR,ERROR,*999)
              ENDIF
            ENDDO !component_idx
          ENDDO !gauss_idx

          !Call surface pressure term here: should only be executed if THIS element has surface pressure on it (direct or incremented)
          IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BOUNDARY_ELEMENT.AND. &
            & TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS>0) THEN    !
            CALL FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,var1,var2,ERR,ERROR,*999)
          ENDIF
        END SELECT
        !Gravity loading term
        IF(ASSOCIATED(RHS_VECTOR)) THEN
          IF(ASSOCIATED(SOURCE_FIELD)) THEN
            IF(ASSOCIATED(MATERIALS_FIELD%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR)) THEN
              DENSITY_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
                & DENSITY_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
              DENSITY_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR
              IF(DARCY_DENSITY) THEN
                DARCY_MATERIALS_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS( &
                  & FIELD_U1_VARIABLE_TYPE)%PTR
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
                  & DARCY_MATERIALS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                DARCY_MATERIALS_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U1_VARIABLE_TYPE)%PTR
              ENDIF
              IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                SOURCE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
                  & SOURCE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
                SOURCE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR

                DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
                  GAUSS_WEIGHT=DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)
                  CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & SOURCE_INTERPOLATED_POINT,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx,EQUATIONS%INTERPOLATION% &
                    & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                    & DENSITY_INTERPOLATED_POINT,ERR,ERROR,*999)
                  IF(DARCY_DENSITY) THEN
                    CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                      & DARCY_MATERIALS_INTERPOLATED_POINT,ERR,ERROR,*999)
                    !Account for separate fluid and solid proportions and densities
                    !Total lagrangian density = m_s + m_f = rho^0_s * (1 - phi^0) + rho_f * phi
                    !By assuming solid incompressibility, phi = (J - 1 + phi^0)
                    !\todo: Think about how this fits in with the constitutive relation, and what happens when the solid
                    !isn't incompressible. Can we assume the solid is incompressible if we aren't enforcing that in the
                    !constitutive relation?
                    DENSITY=DENSITY_INTERPOLATED_POINT%VALUES(1,1)*(1.0_DP-DARCY_MATERIALS_INTERPOLATED_POINT%VALUES(8,1)) + &
                      & DARCY_MATERIALS_INTERPOLATED_POINT%VALUES(7,1)*(Jznu-1.0_DP+DARCY_MATERIALS_INTERPOLATED_POINT%VALUES(8,1))
                  ELSE
                    DENSITY=DENSITY_INTERPOLATED_POINT%VALUES(1,1)
                  ENDIF
                  CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
                    & GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  element_dof_idx=0
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    DEPENDENT_BASIS=>DEPENDENT_FIELD%VARIABLES(var1)%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                    DO parameter_idx=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      element_dof_idx=element_dof_idx+1
                      RHS_VECTOR%ELEMENT_VECTOR%VECTOR(element_dof_idx)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(element_dof_idx) + &
                        & DENSITY*SOURCE_INTERPOLATED_POINT%VALUES(component_idx,1) * &
                        & DEPENDENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,NO_PART_DERIV,gauss_idx)*GAUSS_WEIGHT * &
                        & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN
                    ENDDO
                  ENDDO
                ENDDO !gauss_idx
              ENDIF
            ENDIF
          ENDIF
        ELSE
          CALL FLAG_ERROR("RHS vector is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS5) THEN
      !Output element residual vector for first element
      IF(ELEMENT_NUMBER == 1) THEN
        NDOFS = 0
        FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLES(var1) ! 'U' variable
        DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          SELECT CASE(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE)
          CASE(FIELD_NODE_BASED_INTERPOLATION)
            MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            NDOFS = NDOFS + DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"EP: ",DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS,ERR,ERROR,*999)
          CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
            NDOFS = NDOFS + 1
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"EP: ",1,ERR,ERROR,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Interpolation type " &
              & //TRIM(NUMBER_TO_VSTRING(FIELD_VARIABLE%COMPONENTS(mh)%INTERPOLATION_TYPE,"*",ERR,ERROR))// &
              & " is not valid for a finite elasticity equation.",ERR,ERROR,*999)
          END SELECT
        END DO
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NDOFS: ",NDOFS,ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number * (Fin.Elast.):",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number (Fin.Elast.): ", &
          & ELEMENT_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
          & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(:), &
          & '(4(X,E13.6))','4(4(X,E13.6))',ERR,ERROR,*999)
      ENDIF
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN

999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for a finite elasticity finite element equations set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD

    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE)
        DEPENDENT_FIELD=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
        CALL FiniteElasticity_StrainCalculate(EQUATIONS_SET,DEPENDENT_FIELD, &
          & FIELD_U1_VARIABLE_TYPE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
          & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
        !Do nothing ???
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE")
    RETURN

999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_PRE_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for a finite elasticity finite element equations set.
  SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, & 
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, & 
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, & 
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE,&
          & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
        !Do nothing ???
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity equation type of an elasticity equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE")
    RETURN

999 CALL ERRORS("FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_FINITE_ELEMENT_POST_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Calculated an output field for a finite elasticity equations set.
  SUBROUTINE FiniteElasticityEquationsSet_DerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate the output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived field type to calculate. \see EQUATIONS_SET_CONSTANTS_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: derivedVariable

    CALL ENTERS("FiniteElasticityEquationsSet_DerivedVariableCalculate",err,error,*999)

    NULLIFY(derivedVariable)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) THEN
        CALL FLAG_ERROR("Equations set has not been finished.",err,error,*999)
      ELSE
        IF(ASSOCIATED(equationsSet%equations)) THEN
          CALL Equations_DerivedVariableGet(equationsSet%equations,derivedType,derivedVariable,err,error,*999)
          SELECT CASE(derivedType)
          CASE(EQUATIONS_SET_DERIVED_STRAIN)
            CALL FiniteElasticity_StrainCalculate(equationsSet, &
              & derivedVariable%field,derivedVariable%variable_type,err,error,*999)
          CASE(EQUATIONS_SET_DERIVED_STRESS)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Equations set derived field type of "//TRIM(NUMBER_TO_VSTRING(derivedType,"*",err,error))// &
              & " is not valid for a finite elasticity equations set type.",err,error,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Equations set equations are not associated.",err,error,*999)
        END IF
      END IF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    END IF

    CALL EXITS("FiniteElasticityEquationsSet_DerivedVariableCalculate")
    RETURN
999 CALL ERRORS("FiniteElasticityEquationsSet_DerivedVariableCalculate",err,error)
    CALL EXITS("FiniteElasticityEquationsSet_DerivedVariableCalculate")
    RETURN 1
  END SUBROUTINE FiniteElasticityEquationsSet_DerivedVariableCalculate

  !
  !================================================================================================================================
  !

  !>Calculated the strain field for a finite elasticity finite element equations set.
  SUBROUTINE FiniteElasticity_StrainCalculate(equationsSet,strainField,strainFieldVariableType,err,error,*)

    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate strain for
    TYPE(FIELD_TYPE), POINTER, INTENT(INOUT) :: strainField !<The field to store the strain in.
    INTEGER(INTG), INTENT(IN) :: strainFieldVariableType !<The field variable type of the output field to store the strain in.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,FIBRE_FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT, &
      & DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: gauss_idx,i,NUMBER_OF_TIMES,componentIdx
    INTEGER(INTG) :: element_idx,ne
    INTEGER(INTG) :: FIELD_VAR_TYPE
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: var1 ! Variable number corresponding to 'U' in single physics case
    INTEGER(INTG) :: var2 ! Variable number corresponding to 'DELUDLEN' in single physics case
    REAL(DP) :: DZDNU(3,3),E(3,3),AZL(3,3),AZU(3,3),DZDNUT(3,3)
    REAL(DP) :: Jznu,Jxxi,I3
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),SYSTEM_ELAPSED,SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1)

    CALL ENTERS("FiniteElasticity_StrainCalculate",err,error,*999)

    NULLIFY(GEOMETRIC_BASIS,DEPENDENT_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(DEPENDENT_FIELD,GEOMETRIC_FIELD,FIBRE_FIELD)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS,FIBRE_INTERPOLATION_PARAMETERS)
    NULLIFY(DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT)
    NULLIFY(DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS,DEPENDENT_INTERPOLATED_POINT_METRICS)
    NULLIFY(DECOMPOSITION)

    IF(ASSOCIATED(equationsSet)) THEN
      EQUATIONS=>equationsSet%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN 
        NUMBER_OF_DIMENSIONS=equationsSet%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS

        !Check the provided strain field has appropriate components and interpolation
        IF(ASSOCIATED(strainField)) THEN
          CALL FIELD_VARIABLE_TYPE_CHECK(strainField,strainFieldVariableType,err,error,*999)
          SELECT CASE(NUMBER_OF_DIMENSIONS)
          CASE(3)
            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(strainField,strainFieldVariableType,6,err,error,*999)
            DO componentIdx=1,6
              CALL FIELD_COMPONENT_INTERPOLATION_CHECK(strainField,strainFieldVariableType,componentIdx, &
                & FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
            END DO
          CASE(2)
            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(strainField,strainFieldVariableType,3,err,error,*999)
            DO componentIdx=1,3
              CALL FIELD_COMPONENT_INTERPOLATION_CHECK(strainField,strainFieldVariableType,componentIdx, &
                & FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
            END DO
          CASE(1)
            CALL FLAG_ERROR("1D strain calculation not implemented.",err,error,*999)
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid dimension of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
              & " for the equations set.",err,error,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Strain field is not associated.",err,error,*999)
        END IF
      
        !Which variables are we working with - find the variable pair used for this equations set
        !\todo: put in checks for all the objects/mappings below TODO

        var1=equationsSet%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%VARIABLE_NUMBER ! number for 'U'
        var2=equationsSet%EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE%VARIABLE_NUMBER ! number for 'DELUDELN'

        !Grab pointers: fields, decomposition, basis
        NUMBER_OF_DIMENSIONS=equationsSet%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS

        GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
        DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
        FIBRE_FIELD    =>EQUATIONS%INTERPOLATION%FIBRE_FIELD
        DEPENDENT_NUMBER_OF_COMPONENTS=DEPENDENT_FIELD%VARIABLES(var1)%NUMBER_OF_COMPONENTS
        
        DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
        MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER

        !Grab interpolation parameters
        FIELD_VAR_TYPE=equationsSet%EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%VARIABLE_TYPE
        DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR
        GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        IF(ASSOCIATED(FIBRE_FIELD)) THEN
          FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        END IF

        ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & MAPPINGS%ELEMENTS

        NUMBER_OF_TIMES=0
        
        !Loop over the internal elements
        DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
          NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
          ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)

          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
          DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS       
          DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

          NUMBER_OF_XI=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_XI

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
            & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          IF(ASSOCIATED(FIBRE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
              & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          END IF
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
            & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)

          !Point interpolation pointer
          GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          GEOMETRIC_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
          IF(ASSOCIATED(FIBRE_FIELD)) THEN
            FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          END IF
          DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
          DEPENDENT_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR

          !Loop over gauss points
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            !Interpolate dependent, geometric, fibre fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            IF(ASSOCIATED(FIBRE_FIELD)) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)
            END IF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticityGaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,Jxxi,ERR,ERROR,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ne,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,ERR,ERROR,*999)
            ENDIF

            !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
            !E = Green-Lagrange strain tensor = 0.5*(C-I)

            CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
            CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)

            E = 0.5_DP*AZL 
            DO i=1,3 !NUMBER_OF_DIMENSIONS
              E(i,i)=E(i,i)-0.5_DP
            ENDDO

            ! we only want to store the indepent components of the STRAIN FIELD
            IF(NUMBER_OF_DIMENSIONS==3) THEN
              ! 3 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_13, U_22, U_23, U_33 (upper triangular matrix)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,1,E(1,1),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,2,E(1,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,3,E(1,3),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,4,E(2,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,5,E(2,3),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,6,E(3,3),ERR,ERROR,*999)
            ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
              ! 2 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_22 (upper triangular matrix)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,1,E(1,1),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,2,E(1,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,3,E(2,2),ERR,ERROR,*999)
            ELSE !NUMBER_OF_DIMENSIONS
              LOCAL_ERROR="Only 2 dimensional and 3 dimensional problems are implemented at the moment."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF !NUMBER_OF_DIMENSIONS
          ENDDO !gauss_idx
        ENDDO !element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
        !Output timing information if required
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
          CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
          USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
          SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
          ELEMENT_USER_ELAPSED=USER_ELAPSED
          ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for strain field (internal elements) calculation = ", &
            & USER_ELAPSED,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for strain field (internal elements) calculation = ", &
            & SYSTEM_ELAPSED,ERR,ERROR,*999)
        ENDIF !EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT
        !Output timing information if required
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
          CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
          USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
          SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
            & ERR,ERROR,*999)              
        ENDIF !EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT

        !Loop over the boundary and ghost elements
        DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
          NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
          ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)

          DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS       
          DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS

          NUMBER_OF_XI=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_XI

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
            & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
            & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ne, &
            & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)

          !Point interpolation pointer
          GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          GEOMETRIC_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR
          FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
          DEPENDENT_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_TYPE)%PTR

          !Loop over gauss points
          DO gauss_idx=1,DEPENDENT_NUMBER_OF_GAUSS_POINTS
            !Interpolate dependent, geometric, fibre fields
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI,DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & GEOMETRIC_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FIBRE_INTERPOLATED_POINT,ERR,ERROR,*999)

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticityGaussDeformationGradientTensor(DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,Jxxi,ERR,ERROR,*999)

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ne,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",gauss_idx,ERR,ERROR,*999)
            ENDIF

            !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
            !E = Green-Lagrange strain tensor = 0.5*(C-I)

            CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
            CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999)

            E = 0.5_DP*AZL 
            DO i=1,3 !NUMBER_OF_DIMENSIONS ???
              E(i,i)=E(i,i)-0.5_DP
            ENDDO

            ! we only want to store the indepent components of the STRAIN FIELD
            IF(NUMBER_OF_DIMENSIONS==3) THEN
              ! 3 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_13, U_22, U_23, U_33 (upper triangular matrix)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,1,E(1,1),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,2,E(1,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,3,E(1,3),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,4,E(2,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,5,E(2,3),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,6,E(3,3),ERR,ERROR,*999)
            ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
              ! 2 dimensional problem
              ! ORDER OF THE COMPONENTS: U_11, U_12, U_22 (upper triangular matrix)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,1,E(1,1),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,2,E(1,2),ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_GAUSS_POINT(strainField,strainFieldVariableType,FIELD_VALUES_SET_TYPE, &
                & gauss_idx,ne,3,E(2,2),ERR,ERROR,*999)
            ELSE !NUMBER_OF_DIMENSIONS
              LOCAL_ERROR="Only 2 dimensional and 3 dimensional problems are implemented at the moment."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF !NUMBER_OF_DIMENSIONS
          ENDDO !gauss_idx
        ENDDO !element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH                 
        !Output timing information if required
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
          CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
          CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
          USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
          SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
          ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
          ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for strain field (boundary+ghost equations) calculation = ", &
            & USER_ELAPSED,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for strain field (boundary+ghost equations) calculation = ", &
            & SYSTEM_ELAPSED,ERR,ERROR,*999)
          IF(NUMBER_OF_TIMES>0) THEN
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for strain field calculation = ", &
              & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for strain field calculation = ", &
              & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
          ENDIF
        ENDIF !EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FiniteElasticity_StrainCalculate")
    RETURN
999 CALL ERRORS("FiniteElasticity_StrainCalculate",err,error)
    CALL EXITS("FiniteElasticity_StrainCalculate")
    RETURN 1
  END SUBROUTINE FiniteElasticity_StrainCalculate

  !
  !================================================================================================================================
  !

  !>Calculate the Green-Lagrange strain tensor at a given element xi location.
  SUBROUTINE FiniteElasticity_StrainInterpolateXi(equationsSet,userElementNumber,xi,values,err,error,*)
    ! Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate strain for
    INTEGER(INTG), INTENT(IN) :: userElementNumber
    REAL(DP), INTENT(IN) :: xi(:)
    REAL(DP), INTENT(OUT) :: values(6) !<The interpolated strain tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string.
    ! Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint, &
      & fibreInterpolatedPoint,dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricInterpolatedPointMetrics, &
      & dependentInterpolatedPointMetrics
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(BASIS_TYPE), POINTER :: elementBasis
    LOGICAL :: userElementExists,ghostElement
    INTEGER(INTG) :: dependentVarType,meshComponentNumber
    INTEGER(INTG) :: numberOfDimensions,numberOfXi
    INTEGER(INTG) :: localElementNumber,i
    REAL(DP) :: dZdNu(3,3),dZdNuT(3,3),AZL(3,3),E(3,3)
    REAL(DP) :: JXXi

    CALL Enters("FiniteElasticity_StrainInterpolateXi",err,error,*999)

    NULLIFY(equations)
    NULLIFY(dependentField)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(fibreInterpolatedPoint)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(decomposition)
    NULLIFY(decompositionTopology)
    NULLIFY(domainTopology)
    NULLIFY(elementBasis)

    IF(.NOT.ASSOCIATED(equationsSet)) THEN
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF
    equations=>equationsSet%equations
    IF(.NOT.ASSOCIATED(equations)) THEN
      CALL FlagError("Equations set equations is not associated.",err,error,*999)
    END IF

    nonlinearMapping=>equations%equations_mapping%nonlinear_mapping
    IF(.NOT.ASSOCIATED(equations)) THEN
      CALL FlagError("Equations nonlinear mapping is not associated.",err,error,*999)
    END IF
    dependentVarType=nonlinearMapping%residual_variables(1)%ptr%variable_type

    IF(.NOT.ASSOCIATED(equations%interpolation)) THEN
      CALL FlagError("Equations interpolation is not associated.",err,error,*999)
    END IF
    dependentField=>equations%interpolation%dependent_field
    IF(.NOT.ASSOCIATED(dependentField)) THEN
      CALL FlagError("Equations dependent field is not associated.",err,error,*999)
    END IF
    decomposition=>dependentField%decomposition
    IF(.NOT.ASSOCIATED(decomposition)) THEN
      CALL FlagError("Dependent field decomposition is not associated.",err,error,*999)
    END IF
    CALL DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(decomposition,meshComponentNumber,err,error,*999)
    decompositionTopology=>decomposition%topology
    domainTopology=>decomposition%domain(meshComponentNumber)%ptr%topology
    CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decompositionTopology,userElementNumber, &
      & userElementExists,localElementNumber,ghostElement,err,error,*999)
    IF(.NOT.userElementExists) THEN
      CALL FlagError("The specified user element number of "// &
        & TRIM(NumberToVstring(userElementNumber,"*",err,error))// &
        & " does not exist in the decomposition for the dependent field.",err,error,*999)
    END IF
    CALL DomainTopology_ElementBasisGet( &
      & domainTopology,userElementNumber,elementBasis,err,error,*999)

    !Get the interpolation parameters for this element
    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
      & equations%interpolation%geometric_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
    IF(ASSOCIATED(equations%interpolation%fibre_interp_parameters)) THEN
      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
        & equations%interpolation%fibre_interp_parameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
    END IF
    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,localElementNumber, &
      & equations%interpolation%dependent_interp_parameters(dependentVarType)%ptr,err,error,*999)

    !Get interpolated points
    geometricInterpolatedPoint=>equations%interpolation%geometric_interp_point(FIELD_U_VARIABLE_TYPE)%ptr
    IF(ASSOCIATED(equations%interpolation%fibre_interp_point)) THEN
      fibreInterpolatedPoint=>equations%interpolation%fibre_interp_point(FIELD_U_VARIABLE_TYPE)%ptr
    END IF
    dependentInterpolatedPoint=>equations%interpolation%dependent_interp_point(dependentVarType)%ptr

    !Get interpolated point metrics
    geometricInterpolatedPointMetrics=>equations%interpolation% &
      & geometric_interp_point_metrics(FIELD_U_VARIABLE_TYPE)%ptr
    dependentInterpolatedPointMetrics=>equations%interpolation% &
      & dependent_interp_point_metrics(dependentVarType)%ptr


    !Interpolate fields at xi position
    CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,dependentInterpolatedPoint,err,error,*999)
    CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,geometricInterpolatedPoint,err,error,*999)
    IF(ASSOCIATED(fibreInterpolatedPoint)) THEN
      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,xi,fibreInterpolatedPoint,err,error,*999)
    END IF

    ! Calculate field metrics
    CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE( &
      & elementBasis%number_of_xi,geometricInterpolatedPointMetrics,err,error,*999)
    CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE( &
      & elementBasis%number_of_xi,dependentInterpolatedPointMetrics,err,error,*999)

    !Calculate F=dZ/dNU, the deformation gradient tensor at the xi location
    numberOfDimensions=equationsSet%region%coordinate_system%number_of_dimensions
    numberOfXi=elementBasis%number_of_xi
    CALL FiniteElasticityGaussDeformationGradientTensor(dependentInterpolatedPointMetrics, &
      & geometricInterpolatedPointMetrics,fibreInterpolatedPoint, &
      & dZdNu,JXXi,err,error,*999)

    !Calculate E
    CALL MATRIX_TRANSPOSE(dZdNu,dZdNuT,err,error,*999)
    CALL MATRIX_PRODUCT(dZdNuT,dZdNu,AZL,err,error,*999)
    E=0.5_DP*AZL
    DO i=1,3
      E(i,i)=E(i,i)-0.5_DP
    END DO

    !Set output E components
    values(1)=E(1,1)
    values(2)=E(1,2)
    values(3)=E(1,3)
    values(4)=E(2,2)
    values(5)=E(2,3)
    values(6)=E(3,3)

    CALL Exits("FiniteElasticity_StrainInterpolateXi")
    RETURN

999 CALL Errors("FiniteElasticity_StrainInterpolateXi",err,error)
    CALL Exits("FiniteElasticity_StrainInterpolateXi")
    RETURN 1
  END SUBROUTINE FiniteElasticity_StrainInterpolateXi

  !
  !================================================================================================================================
  !

  !Evaluates the Jacobian surface traction (pressure) term of the equilibrium equation. See Rumpel & Schweizerhof for equations. For
  !most physical situations this jacobian should be symmetrical.
  SUBROUTINE FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS
    TYPE(BASIS_PTR_TYPE) :: BASES(3)
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: FACE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: JACOBIAN_MATRIX
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: PRESSURE_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERP_POINT,PRESSURE_INTERP_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: DEPENDENT_INTERP_POINT_METRICS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: DEPENDENT_QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_PTR_TYPE) :: QUADRATURE_SCHEMES(3)
    INTEGER(INTG) :: FACE_NUMBER,NORMAL_COMPONENT
    INTEGER(INTG) :: FIELD_VAR_U_TYPE,FIELD_VAR_DELUDELN_TYPE,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: oh,mh,ms,mhs,nh,ns,nhs,ng,naf 
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_LOCAL_FACES
    INTEGER(INTG) :: SUM_ELEMENT_PARAMETERS,TOTAL_NUMBER_OF_SURFACE_PRESSURE_CONDITIONS
    INTEGER(INTG) :: ELEMENT_BASE_DOF_INDEX(3),NUMBER_OF_FACE_PARAMETERS(3)
    INTEGER(INTG), PARAMETER :: OFF_DIAG_COMP(3)=[0,1,3],OFF_DIAG_DEP_VAR1(3)=[1,1,2],OFF_DIAG_DEP_VAR2(3)=[2,3,3]
    REAL(DP) :: PRESSURE_GAUSS,JGW_PRESSURE,JGW_PRESSURE_W(2)
    REAL(DP) :: TEMPVEC1(4),TEMPVEC2(4)
    LOGICAL :: NONZERO_PRESSURE

    CALL ENTERS("FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    NULLIFY(DEPENDENT_BASIS)
    NULLIFY(DECOMPOSITION)
    NULLIFY(ELEMENT)
    NULLIFY(EQUATIONS,EQUATIONS_MAPPING,EQUATIONS_MATRICES,NONLINEAR_MAPPING,NONLINEAR_MATRICES,JACOBIAN_MATRIX)
    NULLIFY(DEPENDENT_INTERPOLATION_PARAMETERS,PRESSURE_INTERPOLATION_PARAMETERS)
    NULLIFY(DEPENDENT_INTERP_POINT,DEPENDENT_INTERP_POINT_METRICS,PRESSURE_INTERP_POINT)
    NULLIFY(DEPENDENT_FIELD)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(DEPENDENT_QUADRATURE_SCHEME)
    
    NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS 

    EQUATIONS=>EQUATIONS_SET%EQUATIONS
    EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
    NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
    JACOBIAN_MATRIX=>NONLINEAR_MATRICES%JACOBIANS(1)%PTR

    DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
    DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
    MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
    ELEMENT=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)
    NUMBER_OF_LOCAL_FACES=DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
      & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_LOCAL_FACES
    
    FIELD_VARIABLE=>EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
    FIELD_VAR_U_TYPE=EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR%VARIABLE_TYPE
    FIELD_VAR_DELUDELN_TYPE=EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE

    !Surface pressure term calculation: Loop over all faces
    DO naf=1,NUMBER_OF_LOCAL_FACES
      FACE_NUMBER=ELEMENT%ELEMENT_FACES(naf)
      FACE=>DECOMPOSITION%TOPOLOGY%FACES%FACES(FACE_NUMBER)

      !Check if it's a boundary face
      IF(FACE%BOUNDARY_FACE) THEN 
        NORMAL_COMPONENT=ABS(FACE%XI_DIRECTION)  ! This can be a negative number

        PRESSURE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_DELUDELN_TYPE)%PTR
        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,FACE_NUMBER, &
          & PRESSURE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        PRESSURE_INTERP_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_DELUDELN_TYPE)%PTR

        !Check if nonzero surface pressure is defined on the face
        NONZERO_PRESSURE=.FALSE.
        IF(ANY(ABS(PRESSURE_INTERPOLATION_PARAMETERS%PARAMETERS(:,NORMAL_COMPONENT))>ZERO_TOLERANCE)) THEN
          NONZERO_PRESSURE=.TRUE.
        ENDIF
        
        !Nonzero surface pressure found?
        IF(NONZERO_PRESSURE) THEN
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
          DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%FACES%FACES(FACE_NUMBER)%BASIS
          DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_U_TYPE)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DEPENDENT_INTERP_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_U_TYPE)%PTR
          DEPENDENT_INTERP_POINT_METRICS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_U_TYPE)%PTR

          SUM_ELEMENT_PARAMETERS=0
          !Loop over geometric dependent basis functions.
          DO nh=1,NUMBER_OF_DIMENSIONS
            MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
            DEPENDENT_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%FACES%FACES(FACE_NUMBER)%BASIS
            BASES(nh)%PTR=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEMES(nh)%PTR=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
            NUMBER_OF_FACE_PARAMETERS(nh)=DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            ELEMENT_BASE_DOF_INDEX(nh)=SUM_ELEMENT_PARAMETERS
            SUM_ELEMENT_PARAMETERS=SUM_ELEMENT_PARAMETERS+BASES(nh)%PTR%NUMBER_OF_ELEMENT_PARAMETERS
          ENDDO !nh
        
          !Loop over all Gauss points 
          DO ng=1,DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & PRESSURE_INTERP_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & DEPENDENT_INTERP_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
              & DEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)

            PRESSURE_GAUSS=PRESSURE_INTERP_POINT%VALUES(NORMAL_COMPONENT,NO_PART_DERIV)    !Surface pressure at this gauss point
            SELECT CASE(naf)
            CASE(1,3,5) !Local face 2 and 5 are -XI_DIRECTION, because normal is calculated from cross-product g_1 x g_3
              PRESSURE_GAUSS=-PRESSURE_GAUSS
            CASE(2,4,6) !Local face 2 and 5 are -XI_DIRECTION, because normal is calculated from cross-product g_1 x g_3
              !Do nothing
            END SELECT

            JGW_PRESSURE=DEPENDENT_INTERP_POINT_METRICS%JACOBIAN*DEPENDENT_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)*PRESSURE_GAUSS

            !\todo Find a way to correctly multiply with the scale factors if we have an unsymmetric Jacobian.
            !Loop over element columns belonging to geometric dependent variables
            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR1(oh)
              mh=OFF_DIAG_DEP_VAR2(oh)
              JGW_PRESSURE_W(1)=(DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(1,nh)- &
                & DEPENDENT_INTERP_POINT_METRICS%DXI_DX(1,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,nh))*JGW_PRESSURE
              JGW_PRESSURE_W(2)=(DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(2,nh)- &
                & DEPENDENT_INTERP_POINT_METRICS%DXI_DX(2,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,nh))*JGW_PRESSURE
              DO ns=1,NUMBER_OF_FACE_PARAMETERS(nh)
                !Loop over element rows belonging to geometric dependent variables
                nhs=ELEMENT_BASE_DOF_INDEX(nh)+ &
                  & BASES(nh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ns,naf)
                TEMPVEC1(1)=JGW_PRESSURE_W(1)*QUADRATURE_SCHEMES(nh)%PTR% &
                  & GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1),ng)
                TEMPVEC1(2)=JGW_PRESSURE_W(2)*QUADRATURE_SCHEMES(nh)%PTR% &
                  & GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2),ng)
                !TEMPVEC1(3)=-JGW_PRESSURE_W(1)*QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                !TEMPVEC1(4)=-JGW_PRESSURE_W(2)*QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                DO ms=1,NUMBER_OF_FACE_PARAMETERS(mh)
                  mhs=ELEMENT_BASE_DOF_INDEX(mh)+ &
                    & BASES(mh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ms,naf)
                  TEMPVEC2(1)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  TEMPVEC2(2)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                 ! TEMPVEC2(3)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1),ng)
                 ! TEMPVEC2(4)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2),ng)
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                    DOT_PRODUCT(TEMPVEC1(1:2),TEMPVEC2(1:2))* &
                    & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                    & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ns,nh)
                ENDDO !ms    
              ENDDO !ns
            ENDDO !oh
            DO oh=1,OFF_DIAG_COMP(NUMBER_OF_DIMENSIONS)
              nh=OFF_DIAG_DEP_VAR2(oh)
              mh=OFF_DIAG_DEP_VAR1(oh)
              JGW_PRESSURE_W(1)=(DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(1,nh)- &
                & DEPENDENT_INTERP_POINT_METRICS%DXI_DX(1,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,nh))*JGW_PRESSURE
              JGW_PRESSURE_W(2)=(DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(2,nh)- &
                & DEPENDENT_INTERP_POINT_METRICS%DXI_DX(2,mh)*DEPENDENT_INTERP_POINT_METRICS%DXI_DX(3,nh))*JGW_PRESSURE
              DO ns=1,NUMBER_OF_FACE_PARAMETERS(nh)
                !Loop over element rows belonging to geometric dependent variables
                nhs=ELEMENT_BASE_DOF_INDEX(nh)+ &
                  & BASES(nh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ns,naf)
                TEMPVEC1(1)=JGW_PRESSURE_W(1)*QUADRATURE_SCHEMES(nh)%PTR% &
                  & GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1),ng)
                TEMPVEC1(2)=JGW_PRESSURE_W(2)*QUADRATURE_SCHEMES(nh)%PTR% &
                  & GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2),ng)
              !  TEMPVEC1(3)=-JGW_PRESSURE_W(1)*QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
              !  TEMPVEC1(4)=-JGW_PRESSURE_W(2)*QUADRATURE_SCHEMES(nh)%PTR%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                DO ms=1,NUMBER_OF_FACE_PARAMETERS(mh)
                  mhs=ELEMENT_BASE_DOF_INDEX(mh)+ &
                    & BASES(mh)%PTR%ELEMENT_PARAMETERS_IN_LOCAL_FACE(ms,naf)
                  TEMPVEC2(1)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  TEMPVEC2(2)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
               !   TEMPVEC2(3)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(1),ng)
               !   TEMPVEC2(4)=QUADRATURE_SCHEMES(mh)%PTR%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(2),ng)
                  JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)=JACOBIAN_MATRIX%ELEMENT_JACOBIAN%MATRIX(mhs,nhs)+ &
                    DOT_PRODUCT(TEMPVEC1(1:2),TEMPVEC2(1:2))* &
                    & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                    & DEPENDENT_INTERPOLATION_PARAMETERS%SCALE_FACTORS(ns,nh)
                ENDDO !ms    
              ENDDO !ns
            ENDDO !oh
          ENDDO !ng
        ENDIF !Non-zero pressure on face
      ENDIF !Boundary face
    ENDDO !naf

    CALL EXITS("FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE")
    RETURN

999 CALL ERRORS("FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_SURFACE_PRESSURE_JACOBIAN_EVALUATE

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
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_FACE_BASIS,COMPONENT_FACE_BASIS,COMPONENT_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMP_ELEMENT
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMP_FACE
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FACE_DEPENDENT_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: FACE_PRESSURE_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FACE_DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: FACE_DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: FACE_PRESSURE_INTERPOLATED_POINT
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: FACE_QUADRATURE_SCHEME,COMPONENT_FACE_QUADRATURE_SCHEME
    INTEGER(INTG) :: FIELD_VAR_U_TYPE,FIELD_VAR_DUDN_TYPE,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: element_face_idx,face_number,normal_component_idx,gauss_idx
    INTEGER(INTG) :: component_idx,element_base_dof_idx,element_dof_idx,parameter_idx,face_parameter_idx
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_LOCAL_FACES
    REAL(DP) :: PRESSURE_GAUSS,JGW_PRESSURE,JGW_PRESSURE_NORMAL_COMPONENT
    LOGICAL :: NONZERO_PRESSURE

    CALL ENTERS("FINITE_ELASTICITY_SURFACE_PRESSURE_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    NULLIFY(DEPENDENT_FACE_BASIS,COMPONENT_FACE_BASIS,COMPONENT_BASIS)
    NULLIFY(DECOMPOSITION)
    NULLIFY(DECOMP_ELEMENT)
    NULLIFY(DECOMP_FACE)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS,NONLINEAR_MATRICES)
    NULLIFY(DEPENDENT_FIELD,FIELD_VARIABLE)
    NULLIFY(FACE_DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(FACE_DEPENDENT_INTERPOLATED_POINT,FACE_DEPENDENT_INTERPOLATED_POINT_METRICS)
    NULLIFY(FACE_PRESSURE_INTERPOLATION_PARAMETERS,FACE_PRESSURE_INTERPOLATED_POINT)
    NULLIFY(COMPONENT_FACE_QUADRATURE_SCHEME,FACE_QUADRATURE_SCHEME)

    NUMBER_OF_DIMENSIONS=EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS 

    !Grab pointers of interest
    EQUATIONS=>EQUATIONS_SET%EQUATIONS
    NONLINEAR_MATRICES=>EQUATIONS%EQUATIONS_MATRICES%NONLINEAR_MATRICES
    DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
    DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
    MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
    DECOMP_ELEMENT=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)
    
    !Interpolation parameter for metric tensor
    FIELD_VARIABLE=>EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
    FIELD_VAR_U_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
    FIELD_VAR_DUDN_TYPE=EQUATIONS%EQUATIONS_MAPPING%RHS_MAPPING%RHS_VARIABLE_TYPE
    NUMBER_OF_LOCAL_FACES=DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS% &
      & ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_LOCAL_FACES

    !Surface pressure term calculation: Loop over all faces
    DO element_face_idx=1,NUMBER_OF_LOCAL_FACES
      face_number=DECOMP_ELEMENT%ELEMENT_FACES(element_face_idx)
      DECOMP_FACE=>DECOMPOSITION%TOPOLOGY%FACES%FACES(face_number)

      !Check if it's a boundary face
      IF(DECOMP_FACE%BOUNDARY_FACE) THEN !!temporary until MESH_FACE (or equivalent) is available (decomp face includes ghost faces?)
        normal_component_idx=ABS(DECOMP_FACE%XI_DIRECTION)  ! if xi=0, this can be a negative number
        !\todo: will FACE_COMPONENTS be a problem with sector elements? Check this.
        !Get pressure interpolation objects (DELUDELN pressure_values_set_type)
        FACE_PRESSURE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_DUDN_TYPE)%PTR
        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,face_number, &
          & FACE_PRESSURE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
        FACE_PRESSURE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(var2)%PTR

        !Check if nonzero surface pressure is defined on the face
        NONZERO_PRESSURE=.FALSE.
        IF(ANY(ABS(FACE_PRESSURE_INTERPOLATION_PARAMETERS%PARAMETERS(:,normal_component_idx))>ZERO_TOLERANCE)) THEN
          NONZERO_PRESSURE=.TRUE.
        ENDIF

        !Nonzero surface pressure found?
        IF(NONZERO_PRESSURE) THEN
          !Grab some other pointers
          DEPENDENT_FACE_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%FACES%FACES(face_number)%BASIS
          FACE_QUADRATURE_SCHEME=>DEPENDENT_FACE_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          
          FACE_DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_U_TYPE)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & FACE_DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          FACE_DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_U_TYPE)%PTR
          FACE_DEPENDENT_INTERPOLATED_POINT_METRICS=>EQUATIONS%INTERPOLATION% &
            & DEPENDENT_INTERP_POINT_METRICS(FIELD_VAR_U_TYPE)%PTR

          !Start integrating
          ! Note: As the code will look for P(appl) in the *normal* component to the face, the
          !       initial assignment of P(appl) will have to be made appropriately during bc assignment
          DO gauss_idx=1,FACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS 
            !Interpolate p(appl) at gauss point
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FACE_PRESSURE_INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            PRESSURE_GAUSS=FACE_PRESSURE_INTERPOLATED_POINT%VALUES(normal_component_idx,NO_PART_DERIV)    !Surface pressure at this gauss point
            
            SELECT CASE(element_face_idx)
            CASE(1,3,5) !Local face 2 and 5 are -XI_DIRECTION, because normal is calculated from cross-product g_1 x g_3
              PRESSURE_GAUSS=-PRESSURE_GAUSS
            CASE(2,4,6) !Local face 2 and 5 are -XI_DIRECTION, because normal is calculated from cross-product g_1 x g_3
              ! Do nothing
            END SELECT

            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & FACE_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE, &
              & FACE_DEPENDENT_INTERPOLATED_POINT_METRICS,ERR,ERROR,*999)
            
            JGW_PRESSURE=FACE_DEPENDENT_INTERPOLATED_POINT_METRICS%JACOBIAN* &
              & FACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(gauss_idx)*PRESSURE_GAUSS
            element_base_dof_idx=0
            !Loop over 3 components
            DO component_idx=1,NUMBER_OF_DIMENSIONS
              MESH_COMPONENT_NUMBER=FIELD_VARIABLE%COMPONENTS(component_idx)%MESH_COMPONENT_NUMBER
              COMPONENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              COMPONENT_FACE_BASIS=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%FACES%FACES(face_number)%BASIS
              COMPONENT_FACE_QUADRATURE_SCHEME=>COMPONENT_FACE_BASIS% &
                & QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              JGW_PRESSURE_NORMAL_COMPONENT=JGW_PRESSURE*FACE_DEPENDENT_INTERPOLATED_POINT_METRICS%DX_DXI(component_idx,3)
              DO face_parameter_idx=1,COMPONENT_FACE_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                parameter_idx=COMPONENT_BASIS%ELEMENT_PARAMETERS_IN_LOCAL_FACE(face_parameter_idx,element_face_idx)
                element_dof_idx=element_base_dof_idx+parameter_idx
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                  & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ & ! sign: double -'s. p(appl) always opposite to normal'
                  & JGW_PRESSURE_NORMAL_COMPONENT* &
                  & COMPONENT_FACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(face_parameter_idx,NO_PART_DERIV,gauss_idx)
              ENDDO !face_parameter_idx
              !Update element_base_dof_idx
              element_base_dof_idx=element_base_dof_idx+COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
            ENDDO !component_idx
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
  SUBROUTINE FiniteElasticityGaussDeformationGradientTensor(dependentInterpPointMetrics,geometricInterpPointMetrics,&
    & fibreInterpolatedPoint,dZdNu,JXXi,err,error,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: dependentInterpPointMetrics,geometricInterpPointMetrics
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: fibreInterpolatedPoint
    REAL(DP), INTENT(OUT) :: dZdNu(3,3) !<dZdNu(coordinateIdx,coordianteIdx). On return, the deformation gradient tensor
    REAL(DP), INTENT(OUT) :: JXXi       !<On return, The Jacobian of the transformation from Xi to X
    INTEGER(INTG), INTENT(OUT) :: err   !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfXDimensions,numberOfXiDimensions,numberOfZDimensions
    REAL(DP) :: dNuDXi(3,3),dXidNu(3,3),dZdNuTemp(3,3)

    CALL Enters("FiniteElasticityGaussDeformationGradientTensor",err,error,*999)

    IF(ASSOCIATED(dependentInterpPointMetrics)) THEN
      IF(ASSOCIATED(geometricInterpPointMetrics)) THEN
        
        numberOfXDimensions=geometricInterpPointMetrics%NUMBER_OF_X_DIMENSIONS
        numberOfXiDimensions=geometricInterpPointMetrics%NUMBER_OF_XI_DIMENSIONS
        numberOfZDimensions=dependentInterpPointMetrics%NUMBER_OF_X_DIMENSIONS

        CALL CoordinateMaterialSystemCalculate(geometricInterpPointMetrics,fibreInterpolatedPoint, &
          & dNudXi(1:numberOfXDimensions,1:numberOfXiDimensions), &
          & dXidNu(1:numberOfXiDimensions,1:numberOfXDimensions),err,error,*999)
        !dZ/dNu = dZ/dXi * dXi/dNu  (deformation gradient tensor, F)
        CALL MatrixProduct(dependentInterpPointMetrics%DX_DXI(1:numberOfZDimensions,1:numberOfXiDimensions), &
          & dXiDNu(1:numberOfXiDimensions,1:numberOfXDimensions),dZdNuTemp(1:numberOfZDimensions,1:numberOfXDimensions), &
          & err,error,*999)

        JXXi=geometricInterpPointMetrics%JACOBIAN

        IF(numberOfZDimensions == 2) THEN
          dZdNu(1,:) = [dZdNuTemp(1,1),dZdNuTemp(1,2),0.0_DP]
          dZdNu(2,:) = [dZdNuTemp(2,1),dZdNuTemp(2,2),0.0_DP]
          dZdNu(3,:) = [0.0_DP,0.0_DP,1.0_DP]
        ELSE
          dZdNu(1:numberOfZDimensions,1:numberOfXDimensions) = dZdNuTemp(1:numberOfZDimensions,1:numberOfXDimensions)
        ENDIF

        IF(DIAGNOSTICS1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Calculated deformation gradient tensor:",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Z dimensions  = ",numberOfZDimensions,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Number of Xi dimensions = ",numberOfXiDimensions,err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfXDimensions,1,1,numberOfXDimensions, &
            & numberOfXDimensions,numberOfXDimensions,dZdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES, &
            & '("  dZdNu','(",I1,",:)',' :",3(X,E13.6))','(15X,3(X,E13.6))',err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Determinant JXXi = ",JXXi,err,error,*999)
        ENDIF
        
      ELSE
        CALL FlagError("Geometric interpolated point metrics is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Dependent interpolated point metrics is not associated.",err,error,*999)
    ENDIF
    
    CALL Exits("FiniteElasticityGaussDeformationGradientTensor")
    RETURN
999 CALL Errors("FiniteElasticityGaussDeformationGradientTensor",err,error)
    CALL Exits("FiniteElasticityGaussDeformationGradientTensor")
    RETURN 1
  END SUBROUTINE FiniteElasticityGaussDeformationGradientTensor

  !
  !================================================================================================================================
  !

  !>Evaluates the Cauchy stress tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,DARCY_DEPENDENT_INTERPOLATED_POINT, &
      & INDEPENDENT_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set 
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DARCY_DEPENDENT_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INDEPENDENT_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: CAUCHY_TENSOR(:,:)
    REAL(DP), INTENT(OUT) :: Jznu !Determinant of deformation gradient tensor (AZL)
    REAL(DP), INTENT(IN) :: DZDNU(3,3) !Deformation gradient tensor at the Guass point
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE !<The equation subtype
    INTEGER(INTG) :: i,j,PRESSURE_COMPONENT,component_idx,dof_idx
    REAL(DP) :: AZL(3,3),AZU(3,3),DZDNUT(3,3),PIOLA_TENSOR(3,3),E(3,3),P,IDENTITY(3,3),AZLT(3,3),AZUT(3,3)
    REAL(DP) :: AZL_SQUARED(3,3)    
    REAL(DP) :: I1,I2,I3            !Invariants, if needed
    REAL(DP) :: ACTIVE_STRESS_11,ACTIVE_STRESS_22,ACTIVE_STRESS_33 !Active stress to be copied in from independent field.
    REAL(DP) :: TEMP(3,3),TEMPTERM  !Temporary variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP), DIMENSION (:), POINTER :: C !Parameters for constitutive laws
    REAL(DP) :: a, B(3,3), Q !Parameters for orthotropic laws
    REAL(DP) :: ffact,dfdJfact !coupled elasticity Darcy
    INTEGER(INTG) :: DARCY_MASS_INCREASE_ENTRY !position of mass-increase entry in dependent-variable vector
    REAL(DP) :: VALUE,TITIN_VALUE,VAL1,VAL2
    REAL(DP) :: WV_PRIME    

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR,*999)

    NULLIFY(FIELD_VARIABLE)

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
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,E,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    E','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
    ENDIF
    IDENTITY=0.0_DP
    DO i=1,3
      IDENTITY(i,i)=1.0_DP
    ENDDO

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    CASE(EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
    !Form of constitutive model is:
    ! W_hat=c1*(I1_hat-3)+c2*(I2_hat-3)+p*J*C^(-1) + W^v(J)
      ! take W^v(J) = 1/2 * kappa * (J-1)^2
      WV_PRIME = C(3)*(Jznu - 1.0_DP)
      !compute the invariants, I3 a few lines up
      I1 = AZL(1,1) + AZL(2,2) + AZL(3,3)
      CALL MATRIX_PRODUCT(AZL,AZL,AZL_SQUARED,ERR,ERROR,*999)
      I2 = 0.5_DP * (I1**2 - AZL_SQUARED(1,1) - AZL_SQUARED(2,2) - AZL_SQUARED(3,3))
      
      PIOLA_TENSOR=2.0_DP*Jznu**(-2.0_DP/3.0_DP)*((C(1)+C(2)*I1)*IDENTITY-C(2)*AZL &
        & -(C(1)*I1+2.0_DP*C(2)*I2-1.5_DP*WV_PRIME*Jznu**(5.0_DP/3.0_DP))/3.0_DP*AZU) 
   
    CASE(EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
    !Form of constitutive model is:
    ! W_hat=c1*(I1_hat-3)+c2*(I2_hat-3)+p*J*C^(-1)
    
      !compute the invariants, I3 a few lines up
      I1 = AZL(1,1) + AZL(2,2) + AZL(3,3)
      CALL MATRIX_PRODUCT(AZL,AZL,AZL_SQUARED,ERR,ERROR,*999)
      I2 = 0.5_DP * (I1**2 - AZL_SQUARED(1,1) - AZL_SQUARED(2,2) - AZL_SQUARED(3,3))
      
      !compute 2PK
!      PIOLA_TENSOR(1,1) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (C(1) + C(2) * I1 - C(2) * AZL(1,1) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(1,1))
!      PIOLA_TENSOR(1,2) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (-C(2) * AZL(1,2) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(1,2))
!      PIOLA_TENSOR(1,3) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (-C(2) * AZL(1,3) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(1,3))
!      PIOLA_TENSOR(2,1) = PIOLA_TENSOR(1,2)
!      PIOLA_TENSOR(2,2) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (C(1) + C(2) * I1 - C(2) * AZL(2,2) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(2,2))
!      PIOLA_TENSOR(2,3) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (-C(2) * AZL(2,3) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(2,3))
!      PIOLA_TENSOR(3,1) = PIOLA_TENSOR(1,3)
!      PIOLA_TENSOR(3,2) = PIOLA_TENSOR(2,3)
!      PIOLA_TENSOR(3,3) = 2.0_DP * Jznu**(-2.0_DP/3.0_DP) * (C(1) + C(2) * I1 - C(2) * AZL(3,3) &
!                          & - (C(1) * I1 + 2.0_DP * C(2) * I2 - 1.5_DP * P * Jznu**(5.0_DP/3.0_DP)) / 3.0_DP * AZU(3,3))       
      !????
      PIOLA_TENSOR=2.0_DP*Jznu**(-2.0_DP/3.0_DP)*((C(1)+C(2)*I1)*IDENTITY-C(2)*AZL &
        & -(C(1)*I1+2.0_DP*C(2)*I2-1.5_DP*P*Jznu**(5.0_DP/3.0_DP))/3.0_DP*AZU) 
                              
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
      & EQUATIONS_SET_NO_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
      & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+p*(I3-1)
      !Also assumed I3 = det(AZL) = 1.0
      !  Note that because PIOLA = 2.del{W}/del{C}=[...]+2.lambda.J^2.C^{-1}
      !  lambda here is actually half of hydrostatic pressure
      !If subtype is membrane, assume Mooney Rivlin constitutive law
      IF (EQUATIONS_SET_SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
          PIOLA_TENSOR(1,3)=2.0_DP*(C(2)*(-AZL(3,1)))+P*AZU(1,3)
          PIOLA_TENSOR(2,3)=2.0_DP*(C(2)*(-AZL(3,2)))+P*AZU(2,3)
          PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
          PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
          PIOLA_TENSOR(3,3)=2.0_DP*(C(1)+C(2)*(AZL(1,1)+AZL(2,2)))+P*AZU(3,3)
      ELSE
        ! Membrane Equations
        ! Assume incompressible => I3 = 1 => C33(C11 x C22 - C12*C21) = 1
        AZL(3,3) = 1.0_DP / ((AZL(1,1) * AZL(2,2)) - (AZL(1,2) * AZL (2,1)))
        ! Assume Mooney-Rivlin constitutive relation
        P = -1.0_DP*((C(1) + C(2) * (AZL(1,1) + AZL(2,2))) * AZL(3,3))
        ! Assume stress normal to the surface is neglible i.e. PIOLA_TENSOR(:,3) = 0,PIOLA_TENSOR(3,:) = 0
        PIOLA_TENSOR(:,3) = 0.0_DP
        PIOLA_TENSOR(3,:) = 0.0_DP
      ENDIF
      PIOLA_TENSOR(1,1)=2.0_DP*(C(1)+C(2)*(AZL(2,2)+AZL(3,3)))+P*AZU(1,1)
      PIOLA_TENSOR(1,2)=2.0_DP*(     C(2)*(-AZL(2,1)))+P*AZU(1,2)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=2.0_DP*(C(1)+C(2)*(AZL(3,3)+AZL(1,1)))+P*AZU(2,2)


      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE) THEN
      !add active contraction stress value to the trace of the stress tensor - basically adding to hydrostatic pressure.
      !the active stress is stored inside the independent field that has been set up in the user program.
      !for generality we could set up 3 components in independent field for 3 different active stress components
      !1 isotropic value assumed here.
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
          PIOLA_TENSOR(component_idx,component_idx)=PIOLA_TENSOR(component_idx,component_idx)+VALUE
        ENDDO

        ! Following original code doesn't work for parallel stuff, because FIELD_PARAMETER_SET_GET_GAUSS_POINT takes user element
        ! number and ELEMENT_NUMBER is the local domain element number, which is only the same for one domain.

        !CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
        !  &  FIELD_U_VARIABLE_TYPE,&
        !  &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,1,ACTIVE_STRESS_11,ERR,ERROR,*999) ! get the independent field
        !  stress value1111

        !CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
        !  &  FIELD_U_VARIABLE_TYPE,&
        !  &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,2,ACTIVE_STRESS_22,ERR,ERROR,*999) ! get the independent field stress value

        !CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
        !  &  FIELD_U_VARIABLE_TYPE,&
        !  &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,3,ACTIVE_STRESS_33,ERR,ERROR,*999) ! get the independent field stress value

        !PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+ACTIVE_STRESS_11
        !PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+ACTIVE_STRESS_22
        !PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+ACTIVE_STRESS_33
      ENDIF 
     
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
        ! add the active stress component (stored in the independent field) to the 1,1-direction of the 2-PK tensor
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+INDEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
      ELSE IF((EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) .OR. &
        & (EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)) THEN
        !passive anisotropic stiffness -- only in the tension range
        IF(AZL(1,1) > 1.0_DP) THEN
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2)-1)
        ENDIF
        !active stress component
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        dof_idx=FIELD_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
          & ELEMENT_NUMBER)
        CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
        !divide by lambda and multiply by P_max
        VALUE=VALUE/SQRT(AZL(1,1))*C(5)
        
        !HINDAWI paper - force-length relation at the continuum level
!        if((SQRT(AZL(1,1))>0.72_DP).AND.(SQRT(AZL(1,1))<1.68_DP)) then
!          VALUE=VALUE*(-25.0_DP/4.0_DP*AZL(1,1)/1.2_DP/1.2_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.2_DP - 5.25_DP)
!        else
!          VALUE=0.0_DP
!        endif
        
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE
        
        IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
          IF(PIOLA_TENSOR(1,1).GE.0) THEN
            dof_idx=FIELD_VARIABLE%COMPONENTS(2)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(GAUSS_POINT_NUMBER, &
             & ELEMENT_NUMBER)
            CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,dof_idx,TITIN_VALUE,ERR,ERROR,*999)
           !divide by lambda and multiply by P_max and scale by the active stress
           TITIN_VALUE=TITIN_VALUE*VALUE
           PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+TITIN_VALUE
          ENDIF
        ENDIF
      ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE) THEN
        !Additional term for transversely isotropic (fibre-reinforced) materials (Markert, 2005)
        ! W_aniso=c4*(sqrt(I4)^(c5-2)-1/I4)M
        ! with M being the mapping towards the fibre direction, here: I4=C_11
        !C(3)=c4...polynomial coefficient
        !C(4)=c5...power coefficient
        IF(AZL(1,1) > 1.0_DP) THEN ! only in the tension range
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
      ELSE IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE) THEN
        !Isotropic and anisotropic part from above, additionally an active part in fibre direction
        ! W=W_iso+W_aniso+W_act
        !  with W_act=(1/sqrt(I4)*P_max*f*alpha)M
        !C(5)=alpha...activation parameter [0,1]
        IF(AZL(1,1) > 1.0_DP) THEN ! only in the tension range
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2.0_DP)-1.0_DP)
        ENDIF
        IF((SQRT(AZL(1,1))>0.84_DP).AND.(SQRT(AZL(1,1))<1.96_DP)) THEN
          VALUE=(-25.0_DP/4.0_DP*AZL(1,1)/1.4_DP/1.4_DP + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.4_DP - 5.25_DP) !f
          VALUE=VALUE*(1.0_DP/SQRT(AZL(1,1)))*300000000.0_DP*C(5)  ! P_max here as a constant
          PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+VALUE
        ENDIF
      ENDIF


    CASE(EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE)
      !Equations set for transversely isotropic (fibre-reinforced), active contractible bodies consitisting of two materials
      ! The local portion between them is defined by the parameter trans
      ! Material 1 is active contractible, material 2 is only passive
      !W=W_iso+W_aniso+W_act
      ! where the three parts are adopted from above (iso Mooney-Rivlin, aniso Markert, active part)

      !C(1)=c1_m1...Mooney Rivlin parameter material 1
      !C(2)=c2_m1...Mooney Rivlin parameter material 1
      !C(3)=c4_m1...polynomial coefficient (Markert model) material 1
      !C(4)=c5_m1...power coefficient (Markert model) material 1
      !C(5)=c1_m2...Mooney Rivlin parameter material 2
      !C(6)=c2_m2...Mooney Rivlin parameter material 2
      !C(7)=c4_m2...polynomial coefficient (Markert model) material 2
      !C(8)=c5_m2...power coefficient (Markert model) material 2
      !C(9)=alpha...activation parameter [0,1]
      !C(10)=trans...transition parameter [0,1] for the portion between the two materials

      !Weighting the Mooney Rivlin parameters and obtaining resulting c1 and c2
      VAL1=C(1)*C(10)+C(5)*(1.0_DP-C(10))
      VAL2=C(2)*C(10)+C(6)*(1.0_DP-C(10))

      !Mooney-Rivlin for the isotropic part
      PIOLA_TENSOR(1,1)=2.0_DP*(VAL1+VAL2*(AZL(2,2)+AZL(3,3))+P*AZU(1,1))
      PIOLA_TENSOR(1,2)=2.0_DP*(     VAL2*(-AZL(2,1))        +P*AZU(1,2))
      PIOLA_TENSOR(1,3)=2.0_DP*(     VAL2*(-AZL(3,1))        +P*AZU(1,3))
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=2.0_DP*(VAL1+VAL2*(AZL(3,3)+AZL(1,1))+P*AZU(2,2))
      PIOLA_TENSOR(2,3)=2.0_DP*(     VAL2*(-AZL(3,2))        +P*AZU(2,3))
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=2.0_DP*(VAL1+VAL2*(AZL(1,1)+AZL(2,2))+P*AZU(3,3))

      !passive anisotropic part -- only in the tension range (Markert)
        IF(AZL(1,1) > 1.0_DP) THEN
        VAL1=C(3)/AZL(1,1)*(AZL(1,1)**(C(4)/2)-1)
        VAL2=C(7)/AZL(1,1)*(AZL(1,1)**(C(8)/2)-1)
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+(VAL1*C(10)+VAL2*(1.0_DP-C(10)))
        ENDIF

      !active part
        IF((SQRT(AZL(1,1))>0.84_DP).AND.(SQRT(AZL(1,1))<1.96_DP)) THEN
          VALUE=(-25.0_DP/4.0_DP*AZL(1,1)/1.4/1.4 + 25.0_DP/2.0_DP*SQRT(AZL(1,1))/1.4_DP - 5.25_DP)
          VALUE=VALUE*(1/SQRT(AZL(1,1)))*300000000*C(9)*C(10)  ! P_max here as a constant
        ENDIF

    CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1/2 (e^(c2*(I1-3)) - 1)
      ! S = 2*dW/dC + 2pC^-1
      PIOLA_TENSOR=C(1)*C(2)*EXP(C(2)*(AZL(1,1)+AZL(2,2)+AZL(3,3)-3.0_DP))*IDENTITY+2.0_DP*P*AZU
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE)
      !C(1)=Mooney Rivlin parameter
      !C(2)=Mooney Rivlin parameter
      !C(3)=K
      !C(4)=M, Biot modulus
      !C(5)=b, skeleton parameter
      !C(6)=p0, reference pressure

      P=DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV) !Fluid pressure
      CALL MATRIX_TRANSPOSE(AZL,AZLT,ERR,ERROR,*999)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMP=MATMUL(AZL,AZL)
      I2=0.5_DP*(I1**2.0_DP-TEMP(1,1)-TEMP(2,2)-TEMP(3,3))

      CALL EVALUATE_CHAPELLE_FUNCTION(Jznu,ffact,dfdJfact,ERR,ERROR,*999)

      PIOLA_TENSOR=2.0_DP*C(1)*Jznu**(-2.0_DP/3.0_DP)*(IDENTITY-(1.0_DP/3.0_DP)*I1*AZU)
      PIOLA_TENSOR=PIOLA_TENSOR+2.0_DP*C(2)*Jznu**(-4.0_DP/3.0_DP)*(I1*IDENTITY-AZLT-(2.0_DP/3.0_DP)*I2*AZU)
      PIOLA_TENSOR=PIOLA_TENSOR+(C(3)-C(4)*C(5)**2)*(Jznu-1.0_DP)*AZU
      PIOLA_TENSOR=PIOLA_TENSOR-C(5)*(P-C(6))*Jznu*AZU
      PIOLA_TENSOR=PIOLA_TENSOR+0.5_DP*((P-C(6))**2/C(4))*(dfdJfact/(ffact**2))*Jznu*AZU
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
      ! See Holmes MH, Mow VC. The nonlinear characteristics of soft gels and hydrated connective tissues in ultrafiltration.
      ! Journal of Biomechanics. 1990;23(11):1145-1156. DOI: 10.1016/0021-9290(90)90007-P
      ! The form of constitutive relation is:
      ! sigma = sigma^s + sigma^f
      ! sigma^f = -phi^f p I
      ! sigma^s = -phi^s p I + rho_0^s sigma^s_E
      ! sigma^s_E is the effective Cauchy stress obtained by differentiating
      ! the free energy function to get the second Piola-Kirchoff stress tensor:
      ! rho_0^s W^s = c0 exp(c1(I1 - 3) + c2(I2 - 3)) / (I_3^(c1 + 2c2))
      ! Rather than add the "phi^s p I" term to the Cauchy stress, we add it here as "phi^s p J C^-1"
      ! We also set rho_0^s = the solid density * initial solidity, and move the solidity
      ! inside the strain energy density function
      !
      ! c0 = C(1)
      ! c1 = C(2)
      ! c2 = C(3)
      ! phi^s_0 = C(4)

      CALL MATRIX_TRANSPOSE(AZL,AZLT,ERR,ERROR,*999)
      CALL MATRIX_TRANSPOSE(AZU,AZUT,ERR,ERROR,*999)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMP=MATMUL(AZL,AZL)
      I2=0.5_DP*(I1**2.0_DP-TEMP(1,1)-TEMP(2,2)-TEMP(3,3))
      !I3 already defined

      TEMPTERM=2.0_DP*C(4)*C(1)*EXP(C(2)*(I1 - 3.0_DP) + C(3)*(I2 - 3.0_DP)) / (I3**(C(2)+2.0_DP*C(3)))
      PIOLA_TENSOR=C(2)*TEMPTERM*IDENTITY + C(3)*TEMPTERM*(I1*IDENTITY-AZLT) - (C(2)+2.0_DP*C(3))*TEMPTERM*AZUT
      PIOLA_TENSOR=PIOLA_TENSOR - DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)*Jznu*AZU
    
    CASE(EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE)
    ! For of constitutive model is:
    ! W = 0.5lambda*tr(E)^2 + mu*tr(E^2)
    ! S = dW/dE = lambda*tr(E)Identity + 2muE
      PIOLA_TENSOR(1,3)=(2.0_DP*C(2)*E(1,3))+(2.0_DP*P*AZU(1,3))
      PIOLA_TENSOR(2,3)=(2.0_DP*C(2)*E(2,3))+(2.0_DP*P*AZU(2,3))
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=C(1)*(E(1,1)+E(2,2)+E(3,3))+(2.0_DP*E(3,3)*C(2)+(2.0_DP*P*AZU(3,3)))
      
      PIOLA_TENSOR(1,1)=C(1)*(E(1,1)+E(2,2)+E(3,3))+(2.0_DP*E(1,1)*C(2)+(2.0_DP*P*AZU(1,1)))
      PIOLA_TENSOR(1,2)=(2.0_DP*C(2)*E(1,2))+(2.0_DP*P*AZU(1,2))
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=C(1)*(E(1,1)+E(2,2)+E(3,3))+(2.0_DP*E(2,2)*C(2)+(2.0_DP*P*AZU(2,2)))

      CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
        &  FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,1,ACTIVE_STRESS_11,ERR,ERROR,*999) ! get the independent field stress value

      CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
        &  FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,2,ACTIVE_STRESS_22,ERR,ERROR,*999) ! get the independent field stress value

      CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
        &  FIELD_U_VARIABLE_TYPE,&
        &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,3,ACTIVE_STRESS_33,ERR,ERROR,*999) ! get the independent field stress value

      PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+ACTIVE_STRESS_11
      PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+ACTIVE_STRESS_22
      PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+ACTIVE_STRESS_33

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
    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
      ! W=C1*exp*(Q) + p(J-1)
      ! Q=C2*E(1,1)^2 + C3*(E(2,2)^2+E(3,3)^2+2*E(2,3)*E(3,2)) + 2*C4*(E(1,2)*E(2,1)+E(1,3)*E(3,1))
      Q=C(2)*E(1,1)**2 + C(3)*(E(2,2)**2+E(3,3)**2+2.0_DP*E(2,3)**2) + 2.0_DP*C(4)*(E(1,2)**2+E(1,3)**2)
      TEMPTERM=C(1)*exp(Q) ! iso term
      PIOLA_TENSOR(1,1) = C(2) * E(1,1)
      PIOLA_TENSOR(2,2) = C(3) * E(2,2)
      PIOLA_TENSOR(3,3) = C(3) * E(3,3)
      PIOLA_TENSOR(1,2) = C(4) * E(1,2)
      PIOLA_TENSOR(2,1) = PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(1,3) = C(4) * E(1,3)
      PIOLA_TENSOR(3,1) = PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2) = C(3) * E(2,3)
      PIOLA_TENSOR(2,3) = PIOLA_TENSOR(3,2)
      PIOLA_TENSOR = PIOLA_TENSOR * 2.0_DP * TEMPTERM
      ! pressure terms
      !PIOLA_TENSOR = PIOLA_TENSOR + 2.0_DP*p*Jznu*AZU   ! is Jznu required here, or is it omitted everywhere else?
      PIOLA_TENSOR = PIOLA_TENSOR - p*AZU   ! is Jznu required here, or is it omitted everywhere else?
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
      !add active contraction stress value to the trace of the stress tensor - basically adding to hydrostatic pressure.
      !the active stress is stored inside the independent field that has been set up in the user program.
      !for generality we could set up 3 components in independent field for 3 different active stress components
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
          PIOLA_TENSOR(component_idx,component_idx)=PIOLA_TENSOR(component_idx,component_idx)+VALUE
        ENDDO
      ENDIF 
    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE)
      ! W=a*(exp(b(I1-3))-1) + c*(exp(d(alpha-1)^2)-1)
      ! a=C(1), b=C(2), c=C(3), d=C(4)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      PIOLA_TENSOR(1,1)=C(1)*C(2)*EXP(C(2)*(I1-3))+ &
        & C(3)*2.0_DP*(SQRT(AZL(1,1))-1)*C(4)*EXP(C(4)*(SQRT(AZL(1,1))-1)**2)/(2*SQRT(AZL(1,1)))+P*AZU(1,1)
      PIOLA_TENSOR(2,2)=C(1)*C(2)*EXP(C(2)*(I1-3))+P*AZU(2,2)
      PIOLA_TENSOR(3,3)=C(1)*C(2)*EXP(C(2)*(I1-3))+P*AZU(3,3)
      PIOLA_TENSOR(1,2)=P*AZU(1,2)
      PIOLA_TENSOR(1,3)=P*AZU(1,3)
      PIOLA_TENSOR(2,3)=P*AZU(2,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR=PIOLA_TENSOR*2.0_DP
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
       ENDDO
      ENDDO

      IF(EQUATIONS_SET_SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
        CALL FINITE_ELASTICITY_PIOLA_ADD_ACTIVE_CONTRACTION(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
             & EQUATIONS_SET%EQUATIONS%INTERPOLATION%MATERIALS_FIELD, PIOLA_TENSOR(1,1),E(1,1),         &
             & ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*999)
      ENDIF
    CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
      & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
      & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
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

      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE) THEN

        CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          &  FIELD_U_VARIABLE_TYPE,&
          &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,1,ACTIVE_STRESS_11,ERR,ERROR,*999) ! get the independent field stress value

        CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          &  FIELD_U_VARIABLE_TYPE,&
          &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,2,ACTIVE_STRESS_22,ERR,ERROR,*999) ! get the independent field stress value

        CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, & 
          &  FIELD_U_VARIABLE_TYPE,&
          &  FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,3,ACTIVE_STRESS_33,ERR,ERROR,*999) ! get the independent field stress value

        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+ACTIVE_STRESS_11
        PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+ACTIVE_STRESS_22
        PIOLA_TENSOR(3,3)=PIOLA_TENSOR(3,3)+ACTIVE_STRESS_33
      ENDIF
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE .OR. & 
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE) THEN
        C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
        PIOLA_TENSOR=PIOLA_TENSOR+2.0_DP*C(3)*(I3-SQRT(I3))*AZU
      ELSEIF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE.OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
        & EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
        SELECT CASE (EQUATIONS_SET_SUBTYPE)
        CASE (EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE) !Nearly incompressible
          C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
          !Starting point for this models is above compressible form of 2nd PK tensor
          !Adjust for the modified Ciarlet-Geymonat expression: Eq.(22) of the INRIA paper
          ! Question is: What deviation is to be penalized : (J-1) or (J-1-m/rho) ??? Probably the latter !
          ! However, m/rho is a given 'constant' and, upon differentiation, drops out.
          ! But it is important to retain I3 = J^2, since J ~ 1 + m/rho /= 1
          PIOLA_TENSOR=PIOLA_TENSOR+C(3)*(SQRT(I3)-1.0_DP)*AZU
          DARCY_MASS_INCREASE_ENTRY = 5 !fifth entry
        CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
           &  EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) !Incompressible
          !Constitutive model: W=c1*(I1-3)+c2*(I2-3)+p*(I3-1) 
          ! The term 'p*(I3-1)' gives rise to: '2p I3 AZU'
          ! Retain I3 = J^2, since J ~ 1 + m/rho /= 1 
!         CASE (EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_MR_SUBTYPE)
          !Constitutive model: W=C1*(J1-3)+C2*(J2-3)+C3*(J-1)^2+lambda.(J-1-m/rho)
          !J1 and J2 are the modified invariants, adjusted for volume change (J1=I1*J^(-2/3), J2=I2*J^(-4/3))
          !Strictly speaking this law isn't for an incompressible material, but the fourth equation in the elasticity
          !is used to satisfy a subtly different constraint, which is to require the solid portion of the poroelastic
          !material retains its volume. (This law is applied on the whole pororous body).
          
          PIOLA_TENSOR=0.0_DP
          TEMP=0.0_DP
          
          C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
          C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
          C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)

          !J1 term: del(J1)/del(C)=J^(-2/3)*I-2/3*I_1*J^(-2/3)*C^-1
          TEMPTERM=Jznu**(-2.0_DP/3.0_DP)
          TEMP(1,1)=TEMPTERM
          TEMP(2,2)=TEMPTERM
          TEMP(3,3)=TEMPTERM
          I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
          PIOLA_TENSOR=C(1)* (TEMP-1.0_DP/3.0_DP*I1*TEMPTERM*AZU)

          !J2 term: del(J2)/del(C)=J^(-4/3)*del(I2)/del(C) -4/3*I_2*J^(-4/3)*C^-1
          TEMP=MATMUL(AZL,AZL)  ! C^2
          I2=0.5_DP*(I1**2.0_DP-(TEMP(1,1)+TEMP(2,2)+TEMP(3,3)))
          TEMPTERM=Jznu**(-4.0_DP/3.0_DP)
          !TEMP is now del(I2)/del(C)
          TEMP(1,1)=AZL(2,2)+AZL(3,3)
!           TEMP(1,2)=-2.0_DP*AZL(1,2)
          TEMP(1,2)=-1.0_DP*AZL(1,2)
!           TEMP(1,3)=-2.0_DP*AZL(1,3)
          TEMP(1,3)=-1.0_DP*AZL(1,3)
          TEMP(2,1)=TEMP(1,2)
          TEMP(2,2)=AZL(1,1)+AZL(3,3)
!           TEMP(2,3)=-2.0_DP*AZL(2,3)
          TEMP(2,3)=-1.0_DP*AZL(2,3)
          TEMP(3,1)=TEMP(1,3)
          TEMP(3,2)=TEMP(2,3)
          TEMP(3,3)=AZL(1,1)+AZL(2,2)
          PIOLA_TENSOR=PIOLA_TENSOR+C(2)* (TEMPTERM*TEMP-2.0_DP/3.0_DP*I2*TEMPTERM*AZU)
          
          !J (det(F)) term: (2.C3.(J-1)+lambda)*J.C^-1
          PIOLA_TENSOR=PIOLA_TENSOR+(2.0_DP*C(3)*(Jznu-1.0_DP)+P)*Jznu*AZU

          !Don't forget, it's wrt C so there is a factor of 2 - but not for the pressure !!??
          PIOLA_TENSOR=2.0_DP*PIOLA_TENSOR


          DARCY_MASS_INCREASE_ENTRY = 4 !fourth entry

        END SELECT

!         DARCY_MASS_INCREASE = DARCY_DEPENDENT_INTERPOLATED_POINT%VALUES(DARCY_MASS_INCREASE_ENTRY,NO_PART_DERIV)
! 
!         CALL EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION(AZL,AZU,DARCY_MASS_INCREASE,PIOLA_TENSOR_ADDITION,ERR,ERROR,*999)
! 
!         IF(DIAGNOSTICS1) THEN
!           CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
!             & 3,3,PIOLA_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
!             & '(17X,3(X,E13.6))',ERR,ERROR,*999)
!           CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
!             & 3,3,PIOLA_TENSOR_ADDITION, &
!             & WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR_ADDITION','(",I1,",:)',' :",3(X,E13.6))', &
!             & '(17X,3(X,E13.6))',ERR,ERROR,*999)
!         ENDIF
! 
!         PIOLA_TENSOR = PIOLA_TENSOR + PIOLA_TENSOR_ADDITION
      ENDIF

    CASE (EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
        & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE) ! added by Thomas 2010-04-13
      !Form of the constitutive model is:
      ! W = a/(2*b)*exp[b*(I1-3)] + sum_(i=f,s)[H(I4i-1)*a_i/(2*b_i)*(exp[b_i*(I4i-1)^2]-1)] + a_fs/(2*b_fs)*(exp[b_fs*I8fs^2]-1)
      !where H is the Heaviside step function. Fibres only contribute stiffness if in tension.
      !Also assumed I3 = det(AZL) = J^2 = 1.0  -  incompressible material
      !Assume directions: fibre f_0=[1 0 0], sheet s_0=[0 1 0], (sheet) normal n_0=[0 0 1]
      !Based on: Holzapfel, G. A., & Ogden, R. W. (2009). Constitutive modelling of passive myocardium: A structurally based
      !  framework for material characterization. Philosophical Transactions of the Royal Society A: Mathematical, Physical and
      !  Engineering Sciences, 367(1902), 3445–3475. doi:10.1098/rsta.2009.0091
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
      PIOLA_TENSOR(1,1)=-P*AZU(1,1)+TEMPTERM
      IF(AZL(1,1)>1.0_DP) THEN
        PIOLA_TENSOR(1,1)=PIOLA_TENSOR(1,1)+2.0_DP*C(3)*(AZL(1,1)-1.0_DP)*EXP(C(5)*(AZL(1,1)-1.0_DP)**2.0_DP)
      END IF
      PIOLA_TENSOR(1,2)=-P*AZU(1,2)+C(7)*AZL(1,2)*EXP(C(8)*AZL(1,2)**2.0_DP)
      PIOLA_TENSOR(1,3)=-P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=-P*AZU(2,2)+TEMPTERM
      IF(AZL(2,2)>1.0_DP) THEN
        PIOLA_TENSOR(2,2)=PIOLA_TENSOR(2,2)+2.0_DP*C(4)*(AZL(2,2)-1.0_DP)*EXP(C(6)*(AZL(2,2)-1.0_DP)**2.0_DP)
      END IF
      PIOLA_TENSOR(2,3)=-P*AZU(2,3)
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=-P*AZU(3,3)+TEMPTERM
      
      IF(EQUATIONS_SET_SUBTYPE==EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE) THEN
      !add active contraction stress value to the trace of the stress tensor - basically adding to hydrostatic pressure.
      !the active stress is stored inside the independent field that has been set up in the user program.
      !for generality we could set up 3 components in independent field for 3 different active stress components
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,VALUE,ERR,ERROR,*999)
          PIOLA_TENSOR(component_idx,component_idx)=PIOLA_TENSOR(component_idx,component_idx)+VALUE
        ENDDO
      ENDIF 

    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for a finite elasticity equation type of an elasticity equation set class."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    CALL MATRIX_PRODUCT(DZDNU,PIOLA_TENSOR,TEMP,ERR,ERROR,*999)
    CALL MATRIX_PRODUCT(TEMP,DZDNUT,CAUCHY_TENSOR,ERR,ERROR,*999)
    
    CAUCHY_TENSOR=CAUCHY_TENSOR/Jznu
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  ELEMENT_NUMBER = ",ELEMENT_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  gauss_idx = ",GAUSS_POINT_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,PIOLA_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    PIOLA_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
      CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
        & 3,3,CAUCHY_TENSOR,WRITE_STRING_MATRIX_NAME_AND_INDICES,'("    CAUCHY_TENSOR','(",I1,",:)',' :",3(X,E13.6))', &
        & '(17X,3(X,E13.6))',ERR,ERROR,*999)
    ENDIF
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

  !>Evaluates the Cauchy stress tensor at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_STRESS_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,STRESS_TENSOR,DZDNU,Jznu,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: STRESS_TENSOR(:)
    REAL(DP), INTENT(IN) :: DZDNU(3,3) !Deformation gradient tensor at the gauss point
    REAL(DP), INTENT(IN) :: Jznu !Determinant of deformation gradient tensor (AZL)
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER,GAUSS_POINT_NUMBER !<Element/Gauss point number
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,j,PRESSURE_COMPONENT,component_idx,dof_idx
    REAL(DP) :: P
    REAL(DP) :: I1,I2,I3 !Invariants, if needed
    REAL(DP) :: ACTIVE_STRESS_11,ACTIVE_STRESS_22,ACTIVE_STRESS_33 !Active stress to be copied in from independent field.
    REAL(DP) :: TEMPTERM1,TEMPTERM2 !Temporary variables
    REAL(DP) :: ONETHIRD_TRACE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    REAL(DP) :: MOD_DZDNU(3,3),MOD_DZDNUT(3,3),AZL(3,3),AZU(3,3)
    REAL(DP) :: B(6),E(6),DQ_DE(6)
    REAL(DP), POINTER :: C(:) !Parameters for constitutive laws

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR",ERR,ERROR,*999)

    NULLIFY(FIELD_VARIABLE,C)

    !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)
    !AZU - deformed contravariant tensor; I3 = det(C)

    MOD_DZDNU=DZDNU*Jznu**(-1.0_DP/3.0_DP)
    CALL MATRIX_TRANSPOSE(MOD_DZDNU,MOD_DZDNUT,ERR,ERROR,*999)
    CALL MATRIX_PRODUCT(MOD_DZDNUT,MOD_DZDNU,AZL,ERR,ERROR,*999)
    C=>MATERIALS_INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)

    SELECT CASE(EQUATIONS_SET%SUBTYPE)
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      !Form of constitutive model is:
      !W=c1*(I1-3)+c2*(I2-3)+p/2*(I3-1)

      !Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      I1=AZL(1,1)+AZL(2,2)+AZL(3,3)
      TEMPTERM1=-2.0_DP*C(2)
      TEMPTERM2=2.0_DP*(C(1)+I1*C(2))
      STRESS_TENSOR(1)=TEMPTERM1*AZL(1,1)+TEMPTERM2
      STRESS_TENSOR(2)=TEMPTERM1*AZL(2,2)+TEMPTERM2
      STRESS_TENSOR(3)=TEMPTERM1*AZL(3,3)+TEMPTERM2
      STRESS_TENSOR(4)=TEMPTERM1*AZL(2,1)
      STRESS_TENSOR(5)=TEMPTERM1*AZL(3,1)
      STRESS_TENSOR(6)=TEMPTERM1*AZL(3,2)

      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        !the active stress is stored inside the independent field that has been set up in the user program.
        !for generality we could set up 3 components in independent field for 3 different active stress components
        !!!!! Be aware for modified DZDNU, check if this the right way to do it?
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,TEMPTERM1,ERR,ERROR,*999)
          STRESS_TENSOR(component_idx)=STRESS_TENSOR(component_idx)+TEMPTERM1
        ENDDO
      ENDIF

      !Do push-forward of 2nd Piola tensor. 
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,MOD_DZDNU,Jznu,ERR,ERROR,*999)
      !Calculate isochoric Cauchy tensor (the deviatoric part) and add the volumetric part (the hydrostatic pressure).
      ONETHIRD_TRACE=SUM(STRESS_TENSOR(1:3))/3.0_DP
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)-ONETHIRD_TRACE+P
    CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
      PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
      P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,NO_PART_DERIV)
      B=[2.0_DP*C(2),2.0_DP*C(3),2.0_DP*C(3),C(4),C(4),C(3)] ![2*b_f,2*b_t,2*b_t,b_ft,b_ft,b_t]
      E=[0.5_DP*(AZL(1,1)-1.0_DP),0.5_DP*(AZL(2,2)-1.0_DP),0.5_DP*(AZL(3,3)-1.0_DP),AZL(2,1),AZL(3,1),AZL(3,2)] !(Modified) strain tensor in Voigt form.
      DQ_DE=B*E
      TEMPTERM1=0.5_DP*C(1)*EXP(0.5_DP*DOT_PRODUCT(E,DQ_DE))
      ! Calculate isochoric fictitious 2nd Piola tensor (in Voigt form)
      STRESS_TENSOR=TEMPTERM1*DQ_DE
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE) THEN
        !add active contraction stress values
        !the active stress is stored inside the independent field that has been set up in the user program.
        !for generality we could set up 3 components in independent field for 3 different active stress components
        !!!!! Be aware for modified DZDNU, check if this the right way to do it?
        CALL FIELD_VARIABLE_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
        DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          dof_idx=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP% &
            & GAUSS_POINTS(GAUSS_POINT_NUMBER,ELEMENT_NUMBER)
          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,dof_idx,TEMPTERM1,ERR,ERROR,*999)
          STRESS_TENSOR(component_idx)=STRESS_TENSOR(component_idx)+TEMPTERM1
        ENDDO
      ENDIF
      ! Do push-forward of 2nd Piola tensor. 
      CALL FINITE_ELASTICITY_PUSH_STRESS_TENSOR(STRESS_TENSOR,MOD_DZDNU,Jznu,ERR,ERROR,*999)
      !Calculate isochoric Cauchy tensor (the deviatoric part) and add the volumetric part (the hydrostatic pressure).
      ONETHIRD_TRACE=SUM(STRESS_TENSOR(1:3))/3.0_DP
      STRESS_TENSOR(1:3)=STRESS_TENSOR(1:3)-ONETHIRD_TRACE+P
    CASE DEFAULT
      LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
        & " is not valid for a finite elasticity equation type of an elasticity equation set class."
     CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT
  
    CALL EXITS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR")
    RETURN
    999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_STRESS_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_STRESS_TENSOR

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
    REAL(DP) :: S, LAMBDA, ISO_TA, TA, ACTIVTIME, TIME, DT
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
        &  FIELD_VALUES_SET_TYPE,GAUSS_POINT_NUMBER,ELEMENT_NUMBER, 6+I, QL(I),ERR,ERROR,*999) ! store Q(1) Q(2) Q(3) Lambda for next in 7/8/9/10
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
    REAL(DP), INTENT(OUT) :: DFDZ(:,:,:) !<On return, a matrix containing the derivatives of the basis functions wrt the deformed coordinates
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: COMPONENT_BASIS
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    INTEGER(INTG) :: derivative_idx,component_idx1,component_idx2,xi_idx,parameter_idx 
    REAL(DP) :: DXIDZ(NUMBER_OF_DIMENSIONS,NUMBER_OF_DIMENSIONS),DZDXI(NUMBER_OF_DIMENSIONS,NUMBER_OF_DIMENSIONS)
    REAL(DP) :: Jzxi,DFDXI(NUMBER_OF_DIMENSIONS,64,NUMBER_OF_XI)!temporary until a proper alternative is found
    
    CALL ENTERS("FINITE_ELASTICITY_GAUSS_DFDZ",ERR,ERROR,*999)

    !Initialise DFDXI array
    DFDXI=0.0_DP  ! DFDXI(component_idx,parameter_idx,xi_idx)
    DFDZ=0.0_DP
    DO component_idx2=1,NUMBER_OF_DIMENSIONS !Always 3 spatial coordinates (3D)
      DO xi_idx=1,NUMBER_OF_XI !Thus always 3 element coordinates
        derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx)  !2,4,7      
        DZDXI(component_idx2,xi_idx)=INTERPOLATED_POINT%VALUES(component_idx2,derivative_idx)  !dz/dxi
      ENDDO
    ENDDO

    ! Populate a 3 x 3 square dzdXi if this is a membrane problem in 3D space
    IF (NUMBER_OF_DIMENSIONS == 3 .AND. NUMBER_OF_XI == 2) THEN
        CALL CROSS_PRODUCT(DZDXI(:,1),DZDXI(:,2),DZDXI(:,3),ERR,ERROR,*999)
        DZDXI(:,3) = NORMALISE(DZDXI(:,3),ERR,ERROR)
    ENDIF

    CALL INVERT(DZDXI,DXIDZ,Jzxi,ERR,ERROR,*999) !dxi/dz

    FIELD=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD
    DO component_idx1=1,NUMBER_OF_DIMENSIONS
      COMPONENT_BASIS=>FIELD%VARIABLES(1)%COMPONENTS(component_idx1)%DOMAIN%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(ELEMENT_NUMBER)%BASIS
      QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
      DO parameter_idx=1,COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
        DO xi_idx=1,NUMBER_OF_XI
          derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx)
          DFDXI(component_idx1,parameter_idx,xi_idx)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,derivative_idx, &
            & GAUSS_POINT_NUMBER)
        ENDDO
      ENDDO
    ENDDO

    DO component_idx1=1,NUMBER_OF_DIMENSIONS
      COMPONENT_BASIS=>FIELD%VARIABLES(1)%COMPONENTS(component_idx1)%DOMAIN%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(ELEMENT_NUMBER)%BASIS
      DO component_idx2=1,NUMBER_OF_DIMENSIONS
        DO parameter_idx=1,COMPONENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
          DO xi_idx=1,NUMBER_OF_XI
            DFDZ(parameter_idx,component_idx2,component_idx1)=DFDZ(parameter_idx,component_idx2,component_idx1) + &
              & DFDXI(component_idx1,parameter_idx,xi_idx) * DXIDZ(xi_idx,component_idx2)
          ENDDO 
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
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_COMPONENTS, &
      & NUMBER_OF_DIMENSIONS, NUMBER_OF_DARCY_COMPONENTS,GEOMETRIC_COMPONENT_NUMBER,NUMBER_OF_COMPONENTS_2,component_idx, &
      & derivedIdx,varIdx,variableType,NUMBER_OF_FLUID_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    LOGICAL :: IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD
    INTEGER(INTG) :: num_var,Ncompartments           
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS    
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:)

    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(EQUATIONS_MATERIALS)
    NULLIFY(EQUATIONS_EQUATIONS_SET_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_DATA)
    
    IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD = EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE &
      & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE &
      & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_MEMBRANE_SUBTYPE &
      & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE &
      & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE &
      & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE &
      & .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE               

    NUMBER_OF_DIMENSIONS = EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS

    IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
      NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS + 1
    ELSE
      NUMBER_OF_COMPONENTS = NUMBER_OF_DIMENSIONS
    ENDIF

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, & 
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,&
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, & 
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
          & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Default to FEM solution method
            CALL FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
            !setup equations set field to store number of fluid compartments
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES, &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              ENDIF
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)THEN
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1, 1_INTG, ERR, ERROR, *999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 2, 1_INTG, ERR, ERROR, *999)
              ENDIF
            ENDIF 
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !\todo Check dimension of geometric field
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            ! Check whether a fibre field is required, and if so, make sure it has been set
            SELECT CASE(EQUATIONS_SET%SUBTYPE)
              CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
                  & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
                  & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
                  & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
                  & EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
                  & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,&
                  & EQUATIONS_SET_NO_SUBTYPE, &
                  & EQUATIONS_SET_MEMBRANE_SUBTYPE, &
                  & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
                  & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
                  & EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
              ! pass, fibre field isn't required as the constitutive relation is isotropic
            CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, &
                  & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
                  & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &  
                  & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
                  & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
                  & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, &
                  & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
                  & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &                  
                  & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
                  & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
                  & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE, &
                  & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
                  & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
                  & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
                  & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
                  & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
              IF(.NOT.ASSOCIATED(EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)) CALL FLAG_ERROR( &
                & "Finite elascitiy equations require a fibre field.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity equation"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE .OR. &
                & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
              ! Set up mesh displacement and equations set field info for elasticity Darcy problems
              FIELD_VARIABLE=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
              CALL Field_ParameterSetEnsureCreated(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE,ERR,ERROR,*999)
              IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MULTI_COMPARTMENT_DARCY_SUBTYPE .OR. &
                 EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE) THEN
                !Create the equations set field for multi-compartment Darcy
                EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2

                EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD

                IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                    & GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,& 
                    & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)                
                  DO component_idx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  END DO

                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & ERR,ERROR,*999)
                ELSE
                  !Do nothing
                ENDIF
              ENDIF
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            ! do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          !-----------------------------------------------------------------------
          ! Dependent field setup for single-physics
          !-----------------------------------------------------------------------
          CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, & 
              & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
              & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
              & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
              & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, &
              & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
              & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
              & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
              & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
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
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

                !Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                                   
                IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
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
                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
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
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
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
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
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
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity equation"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT

          !-------------------------------------------------------------------------------
          ! Dependent field setup for elasticity evaluated in CellML
          !-------------------------------------------------------------------------------
          CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                IF(NUMBER_OF_DIMENSIONS==3) THEN
                  NUMBER_OF_COMPONENTS_2 = 6
                ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
                  NUMBER_OF_COMPONENTS_2 = 3
                ELSE
                  CALL FLAG_ERROR("Only 2 and 3 dimensional problems are implemented at the moment",ERR,ERROR,*999)
                ENDIF !NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS_2,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS_2,ERR,ERROR,*999)

                !Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                                   
                IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
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

                !Set the stress and strain components to that of the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS_2
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Set the displacement components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                  
                  IF(IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
                    !Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDIF

                  !Set the stress and strain components to gauss point interpolation
                  DO component_idx=1,NUMBER_OF_COMPONENTS_2
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U2_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
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

              ELSE !EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED

                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                IF(NUMBER_OF_DIMENSIONS==3) THEN
                  NUMBER_OF_COMPONENTS_2 = 6
                ELSE IF(NUMBER_OF_DIMENSIONS==2) THEN
                  NUMBER_OF_COMPONENTS_2 = 3
                ELSE
                  CALL FLAG_ERROR("Only 2 and 3 dimensional problems are implemented at the moment",ERR,ERROR,*999)
                ENDIF !NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_2, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_2, &
                  & ERR,ERROR,*999)

                !Check that the pressure values set type is created here?? (second variable is a DELUDELN type, as checked above)
                !\todo: Decide whether these set_types (previous one as well) is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%PTR)) THEN
                    LOCAL_ERROR="Variable 2 of type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%FIELD%VARIABLES(2)% &
                      & VARIABLE_TYPE,"*",ERR,ERROR))//" does not have a pressure values set type associated."
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                  DO component_idx=1,NUMBER_OF_COMPONENTS_2
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,component_idx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,component_idx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
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
              ENDIF !EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED
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

                !Check that the impermeability flag values set type is created here?? 
                !\todo: Decide whether these set_types is to be created by user or automatically..
                IF(.not.ASSOCIATED(EQUATIONS_SET_SETUP%FIELD%VARIABLES(4)%PARAMETER_SETS% &
                  & SET_TYPE(FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE)%PTR)) THEN
                    LOCAL_ERROR="Variable 4 of type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP% &
                      & FIELD%VARIABLES(4)% &
                      & VARIABLE_TYPE,"*",ERR,ERROR))//" does not have an impermeable flag values set type associated."
                ENDIF

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

              CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                 & FIELD_IMPERMEABLE_FLAG_VALUES_SET_TYPE,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity equation"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          !---------------------------------------------------------------------------------------------
          ! Shared Dependent field setup for multi-physics: elasticity coupled with Darcy fluid pressure
          !---------------------------------------------------------------------------------------------
          CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
            NUMBER_OF_DARCY_COMPONENTS=1 !Only solving for the fluid pressure at the moment
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
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
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)

                !Set labels
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"U",ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
                  & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,"V",ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE,"del V/del n", &
                  & ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"x1",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,2,"x2",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,3,"x3",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & "del x1/del n",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2, &
                  & "del x2/del n",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,3, &
                  & "del x3/del n",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"p",ERR,ERROR,*999)
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & "del p/del n",ERR,ERROR,*999)

                !Elasticity: Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
                !Darcy: Default pressure and mass increase to the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
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
                  !Darcy: Set the pressure and mass increase components to node based interpolation
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
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/) &
                    & ,ERR,ERROR,*999)
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
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)
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
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity equation"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          !-------------------------------------------------------------------------------
          ! Shared Dependent field setup for multi-physics: elasticity coupled with multi-compartment Darcy
          !-------------------------------------------------------------------------------
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                !Get the number of Darcy compartments from the equations set field
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                !Set number of variables to be 2+2*Ncompartments
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(2+2*Ncompartments), &
                   & ERR,ERROR,*999)
                ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
                DO num_var=1,Ncompartments+1
                  VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                ENDDO
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for Darcy with pressure driven by solid: velocity components and mass increase

                DO num_var=1,2*Ncompartments+2
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                    & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                    & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                ENDDO

!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
!                   & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
!                   NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for Darcy with pressure driven by solid: velocity components and mass increase
!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
!                   & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)
!                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
!                   & NUMBER_OF_DARCY_COMPONENTS,ERR,ERROR,*999)

                !Elasticity: Default to the geometric interpolation setup
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO !component_idx
 
                !Set the hydrostatic component to that of the first geometric component
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)            
                DO num_var=3,2*Ncompartments+2
                  !Darcy: Default to the geometric interpolation setup
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  ENDDO !component_idx
                  !Darcy: Default pressure and, if present, mass increase to the first geometric component
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDDO
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity: Set the displacement components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

!                   IF (EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE) THEN
                    !Elasticity: Set the hydrostatic pressure component to node based interpolation
                    !as this is used as the pressure field for the Darcy equations
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
!                   ELSE IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN
!                     !Elasticity: Set the hydrostatic pressure component to element based interpolation
!                     CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
!                       & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                     CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
!                       & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                   ENDIF
                  DO num_var=3,2*Ncompartments+2
                    !Darcy: Set the velocity, pressure and, if present, mass increase components to node based interpolation
                    DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                      CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & VARIABLE_TYPES(num_var),component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    ENDDO !component_idx
                  ENDDO
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
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                !Get the number of Darcy compartments from the equations set field
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,(2+2*Ncompartments),ERR,ERROR,*999)
                ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
                DO num_var=1,Ncompartments+1
                  VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                ENDDO
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,ERR,ERROR,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                NUMBER_OF_DARCY_COMPONENTS=NUMBER_OF_DIMENSIONS+1

                DO num_var=1,2*Ncompartments+2
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_VECTOR_DIMENSION_TYPE, &
                    & ERR,ERROR,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),NUMBER_OF_COMPONENTS, &
                    & ERR,ERROR,*999)

                ENDDO

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity:
                 DO component_idx=1,NUMBER_OF_DIMENSIONS
                   CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                     & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                   CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                     & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                 ENDDO !component_idx
                    !If solid hydrostatic pressure is driving Darcy flow, check that pressure uses node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                 CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)

                DO num_var=3,2*Ncompartments+2
                  !Darcy:
                  DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDDO
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
                DEALLOCATE(VARIABLE_TYPES)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
              ENDIF
              EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
              CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,ERR,ERROR,*999)
              Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
              ALLOCATE(VARIABLE_TYPES(2*Ncompartments+2))
              DO num_var=1,Ncompartments+1
                VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
              ENDDO
              DO num_var=3,2*Ncompartments+2
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & FIELD_INITIAL_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & FIELD_RELATIVE_VELOCITY_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_TYPES(num_var), &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ERR,ERROR,*999)
              ENDDO
              DEALLOCATE(VARIABLE_TYPES)
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
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)

           SELECT CASE(EQUATIONS_SET%SUBTYPE)
           ! ACTIVE CONTRACTION
           CASE(EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
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

           !Mooney Rivlin, St Venant Kirchoff and Compressible active contraction subtype
           CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, & 
             & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
             & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
             NUMBER_OF_COMPONENTS = 3 !one contractile stress value for each of the three directions
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
             
             !Set component to be gauss point based 
             DO component_idx=1,NUMBER_OF_COMPONENTS 
               CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
             ENDDO


           ! COUPLED DARCY
           CASE(EQUATIONS_SET_STANDARD_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE,&
              & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_MULTI_COMPARTMENT_DARCY_INRIA_SUBTYPE, &
              & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_MR_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE)
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
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
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
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE/), &
                 & ERR,ERROR,*999)
               CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,ERR, &
                 & ERROR,*999)
               CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                 & ERR,ERROR,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS, &
                 & ERR,ERROR,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,&
                 & ERR,ERROR,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !DO component_idx=1,NUMBER_OF_DIMENSIONS
                 !  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                 !    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                 !  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                 !    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                 !ENDDO !component_idx
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
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           ! BIOELECTRICS COUPLED TO FINITE ELASTICITY
           CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE)
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
               CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,1,ERR,ERROR,*999)
               CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,ERR,ERROR,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,ERR,ERROR,*999)
               !Default to the first component of the geometric interpolation setup
               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
               CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)

               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set to node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
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
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
               CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR, &
                 & ERROR,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                 & ERR,ERROR,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !do/check nothing???
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
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
             IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
               !Create the auto created independent field
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
               CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                 & FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
               IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
               ENDIF
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_DP_TYPE,ERR,ERROR,*999)
               CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & FIELD_INTG_TYPE,ERR,ERROR,*999)
               IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 1,ERR,ERROR,*999)
               ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 2,ERR,ERROR,*999)
               ENDIF
               CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & NUMBER_OF_DIMENSIONS+1,ERR,ERROR,*999)
               !Default to the first component of the geometric interpolation setup
               CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
               IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                 CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & 2,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
               ENDIF
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !Set to node based interpolation
                 CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_U_VARIABLE_TYPE,1,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                 IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_U_VARIABLE_TYPE,2,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999)
                 ENDIF
                 DO component_idx=1,NUMBER_OF_DIMENSIONS
                   CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                     & FIELD_V_VARIABLE_TYPE,component_idx,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                 ENDDO
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
             ELSE !INDEPENDENT_FIELD_AUTO_CREATED
               !Check the user specified field
               CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
               CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
               !Question:Better to leave it up for the user to play around?
               CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,ERR,ERROR,*999)
               CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR, &
                 & ERROR,*999)
               IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR, &
                   & ERROR,*999)
               ENDIF
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
               CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,ERR,ERROR,*999)
               IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                   & ERR,ERROR,*999)
               ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                 CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2, &
                   & ERR,ERROR,*999)
               ENDIF
               CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                 & ERR,ERROR,*999)
               SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
               CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                 !do/check nothing???
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
             ENDIF !INDEPENDENT_FIELD_AUTO_CREATED

           CASE DEFAULT
             LOCAL_ERROR="The equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
               & " is invalid for an independent field of a finite elasticity equation."
             CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
           END SELECT
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

        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              NUMBER_OF_FLUID_COMPONENTS=0
              SELECT CASE(EQUATIONS_SET%SUBTYPE)
              CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
                & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, & 
                & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
                & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
                & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE)
                NUMBER_OF_COMPONENTS = 3;
              CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                NUMBER_OF_COMPONENTS = 5;
              CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE)
                !\todo Currently the number of components for a membrane problem's material field has been set to 3 in 3D space or
                ! 2 in 2D space to work with a Mooney Rivlin material (2 material parameters) and a membrane thickness parameter 
                ! (only if in 3D space). Extra subtypes will need to be added to use other constitutive relations with 
                ! membrane mechanics problems.
                IF (NUMBER_OF_DIMENSIONS==3) THEN
                  NUMBER_OF_COMPONENTS = 3
                ELSE
                  NUMBER_OF_COMPONENTS = 2
                ENDIF
              CASE(EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;
              CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 5;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 4;
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE)
                NUMBER_OF_COMPONENTS = 5;
              CASE(EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE)
                NUMBER_OF_COMPONENTS = 10;
              CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS = 7;
              CASE(EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE,&
                & EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE)
                NUMBER_OF_COMPONENTS = 3;
              CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, & 
                & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE)
                NUMBER_OF_COMPONENTS = 8;
              CASE(EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;  
              CASE(EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE)
                NUMBER_OF_COMPONENTS = 3;                  
              CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
                & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE, &
                & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE)
                NUMBER_OF_COMPONENTS = 4;
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE)
                NUMBER_OF_COMPONENTS = 6;
                NUMBER_OF_FLUID_COMPONENTS=8
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
                NUMBER_OF_COMPONENTS = 4
                NUMBER_OF_FLUID_COMPONENTS=8
              CASE(EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE)
                NUMBER_OF_COMPONENTS = 2;
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
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)  ! get 1 = x (?) component

                !U variable type is constitutive law parameters
                !V variable type has one component, density
                IF(NUMBER_OF_FLUID_COMPONENTS>0) THEN
                  !If coupled with Darcy pressure equation then a shared material field is used and Darcy material parameters are in U1
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,3,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],ERR,ERROR,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,2,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
                ENDIF
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials",ERR,ERROR,*999)

                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Parameters",ERR,ERROR,*999)

                IF(EQUATIONS_SET%SUBTYPE == EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,ERR,ERROR,*999) ! just 1 component: activation time
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                   & FIELD_V_VARIABLE_TYPE,1 ,FIELD_GAUSS_POINT_BASED_INTERPOLATION,ERR,ERROR,*999) ! gauss pt based interp.
                ELSE
                  !Solid density
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                     & 1,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE,"Density",ERR,ERROR,*999)
                ENDIF

                DO component_idx=1,NUMBER_OF_COMPONENTS
                !Default the materials components to the geometric interpolation setup with constant interpolation
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO

                IF(NUMBER_OF_FLUID_COMPONENTS>0) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                     & NUMBER_OF_FLUID_COMPONENTS,ERR,ERROR,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE,"Fluid Parameters", &
                    & ERR,ERROR,*999)
                ENDIF
                DO component_idx=1,NUMBER_OF_FLUID_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                ENDDO

                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_GET(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES)
                CASE(1)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CASE(2)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
                CASE(3)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                      & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Invalid number of variables. The number of variables for field number "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%FIELD%USER_NUMBER,"*",ERR,ERROR))//" is "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%FIELD%NUMBER_OF_VARIABLES,"*",ERR,ERROR))// &
                    & " but should be either 1, 2 or 3"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                IF (EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES>1) THEN
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                      & 1,ERR,ERROR,*999)
                ENDIF
                IF (EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES>2) THEN
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                      & NUMBER_OF_FLUID_COMPONENTS,ERR,ERROR,*999)
                ENDIF
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
                !Initialise density to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
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
          IF(ASSOCIATED(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)) THEN
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
          ELSE
            CALL FLAG_ERROR("Equations set geometrc field is not associated",ERR,ERROR,*999)
          ENDIF
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%SOURCE% &
                & SOURCE_FIELD,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,"Source Field",ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)

              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,"Gravity",ERR,ERROR,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,1,"g1",ERR,ERROR,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,2,"g2",ERR,ERROR,*999)
              IF(NUMBER_OF_COMPONENTS==3) THEN
                CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,3,"g3",ERR,ERROR,*999)
              ENDIF

              DO component_idx=1,NUMBER_OF_COMPONENTS
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%SOURCE%SOURCE_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
              END DO
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
              IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Finish creating the source field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                !Set the default values for the field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                DO component_idx=1,NUMBER_OF_COMPONENTS-1
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,ERR,ERROR,*999)
                ENDDO
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,NUMBER_OF_COMPONENTS,-9.80665_DP,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
            ENDIF
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
              SELECT CASE(EQUATIONS_SET%SUBTYPE)
              CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
                & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
                !Residual vector also depends on the fluid pressure variable
                CALL EQUATIONS_MAPPING_RESIDUAL_VARIABLES_NUMBER_SET(EQUATIONS_MAPPING,2,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET(EQUATIONS_MAPPING, &
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR,ERROR,*999)
              CASE DEFAULT
                !Single residual variable
                CALL EQUATIONS_MAPPING_RESIDUAL_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              END SELECT
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
              SELECT CASE(EQUATIONS_SET%SUBTYPE)
              CASE(EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, & 
                  & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE,EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
                ! Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(EQUATIONS_MATRICES,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & ERR,ERROR,*999)
              CASE DEFAULT
                  ! Do nothing
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
        CASE(EQUATIONS_SET_SETUP_DERIVED_TYPE)
          ! We want to be able to set which derived variables are calculated before finishing the derived
          ! field, so don't create field variables or check the provided field until the finish action.
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%derived%derivedFieldAutoCreated) THEN
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%derived% &
                & derivedField,ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%derived%derivedField,"Derived Field",ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,FIELD_DEPENDENT_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%derived)) THEN
              ALLOCATE(VARIABLE_TYPES(EQUATIONS_SET%derived%numberOfVariables),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derived field variable types.",ERR,ERROR,*999)
              varIdx=0
              DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES
                IF(EQUATIONS_SET%derived%variableTypes(derivedIdx)/=0) THEN
                  varIdx=varIdx+1
                  VARIABLE_TYPES(varIdx)=EQUATIONS_SET%derived%variableTypes(derivedIdx)
                END IF
              END DO
              IF(EQUATIONS_SET%derived%derivedFieldAutoCreated) THEN
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField, &
                  & EQUATIONS_SET%derived%numberOfVariables,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,VARIABLE_TYPES,ERR,ERROR,*999)
                DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES
                  variableType=EQUATIONS_SET%derived%variableTypes(derivedIdx)
                  IF(variableType/=0) THEN
                    CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                      & FIELD_DP_TYPE,ERR,ERROR,*999)
                    SELECT CASE(derivedidx)
                    CASE(EQUATIONS_SET_DERIVED_STRAIN)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Strain",ERR,ERROR,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & 6,ERR,ERROR,*999)
                    CASE(EQUATIONS_SET_DERIVED_STRESS)
                      CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                      CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%derived%derivedField,variableType,"Stress",ERR,ERROR,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & 6,ERR,ERROR,*999)
                    CASE DEFAULT
                      CALL FLAG_ERROR("The specified derived field type of "//TRIM(NUMBER_TO_VSTRING(derivedIdx,"*",ERR,ERROR))// &
                        & " is not supported for a finite elasticity equations set type.",ERR,ERROR,*999)
                    END SELECT
                  END IF
                END DO
                !Finish creating the derived field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%derived%derivedField,ERR,ERROR,*999)
              ELSE
                !Check the user specified derived field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD, &
                  & EQUATIONS_SET%derived%numberOfVariables,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,ERR,ERROR,*999)

                DO derivedIdx=1,EQUATIONS_SET_NUMBER_OF_DERIVED_TYPES
                  variableType=EQUATIONS_SET%derived%variableTypes(derivedIdx)
                  IF(variableType/=0) THEN
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                    SELECT CASE(derivedidx)
                    CASE(EQUATIONS_SET_DERIVED_STRAIN)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & 6,ERR,ERROR,*999)
                    CASE(EQUATIONS_SET_DERIVED_STRESS)
                      CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                      CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET%derived%derivedField,variableType, &
                        & 6,ERR,ERROR,*999)
                    CASE DEFAULT
                      CALL FLAG_ERROR("The specified derived field type of "//TRIM(NUMBER_TO_VSTRING(derivedIdx,"*",ERR,ERROR))// &
                        & " is not supported for a finite elasticity equations set type.",ERR,ERROR,*999)
                    END SELECT
                  END IF
                END DO
              END IF
            ELSE
              CALL FLAG_ERROR("Equations set derived is not associated.",ERR,ERROR,*999)
            ENDIF
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
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, & 
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, & 
          & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, & 
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
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
      CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_MOONEY_RIVLIN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_STVENANT_KIRCHOFF_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE, &
          & EQUATIONS_SET_TRANS_ISOTROPIC_ACTIVE_TRANSITION_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_POLYNOMIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_COMPRESSIBLE_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE,EQUATIONS_SET_ELASTICITY_DARCY_INRIA_MODEL_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NEARLY_INCOMPRESSIBLE_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_ELAST_MULTI_COMP_DARCY_SUBTYPE,EQUATIONS_SET_TRANSVERSE_ISOTROPIC_GUCCIONE_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_HUMPHREY_YIN_SUBTYPE, &
          & EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_HOLZAPFEL_OGDEN_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_GUCCIONE_ACTIVECONTRACTION_SUBTYPE)
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
    TYPE(SOLVER_TYPE), POINTER :: CELLML_SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(CELLML_SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    NULLIFY(CELLML_EQUATIONS)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE)
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

            IF(PROBLEM%SUBTYPE==PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE) THEN
              !Create the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_EVALUATOR_CREATE(SOLVER,CELLML_SOLVER,ERR,ERROR,*999)
              !Link the CellML evaluator solver to the solver
              CALL SOLVER_LINKED_SOLVER_ADD(SOLVER,CELLML_SOLVER,SOLVER_CELLML_EVALUATOR_TYPE,ERR,ERROR,*999)
            ENDIF
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
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Get the CellML evaluator solver
            CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,CELLML_SOLVER,ERR,ERROR,*999)
            !Create the CellML equations
            CALL CELLML_EQUATIONS_CREATE_START(CELLML_SOLVER,CELLML_EQUATIONS, &
              & ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Get the CellML evaluator solver
            CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,CELLML_SOLVER,ERR,ERROR,*999)
            !Get the CellML equations for the CellML evaluator solver
            CALL SOLVER_CELLML_EQUATIONS_GET(CELLML_SOLVER,CELLML_EQUATIONS,ERR,ERROR,*999)
            !Finish the CellML equations creation
            CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity equation."
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

  !>Sets up the finite elasticity problem.
  SUBROUTINE FiniteElasticity_ContactProblemSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Laplace equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: nonlinearSolver,transformationSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FINITE_ELASTICITY_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(nonlinearSolver)
    NULLIFY(transformationSolver)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE)
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
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,ERR,ERROR,*999)
            !Set the first solver to be a geometric transformation solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,transformationSolver,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(transformationSolver,SOLVER_GEOMETRIC_TRANSFORMATION_TYPE,ERR,ERROR,*999)
            !Set the second solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,nonlinearSolver,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(nonlinearSolver,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(nonlinearSolver,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
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
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,nonlinearSolver,ERR,ERROR,*999)
            !Create the solver equatgions
            CALL SOLVER_EQUATIONS_CREATE_START(nonlinearSolver,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,nonlinearSolver,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(nonlinearSolver,SOLVER_EQUATIONS,ERR,ERROR,*999)
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
      CASE(PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE)
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
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,nonlinearSolver,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(nonlinearSolver,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(nonlinearSolver,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
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
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,nonlinearSolver,ERR,ERROR,*999)
            !Create the solver equatgions
            CALL SOLVER_EQUATIONS_CREATE_START(nonlinearSolver,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,nonlinearSolver,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(nonlinearSolver,SOLVER_EQUATIONS,ERR,ERROR,*999)
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
          & " is not valid for a finite elasticity contact type of an elasticity problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FiniteElasticity_ContactProblemSetup")
    RETURN
999 CALL ERRORS("FiniteElasticity_ContactProblemSetup",ERR,ERROR)
    CALL EXITS("FiniteElasticity_ContactProblemSetup")
    RETURN 1
  END SUBROUTINE FiniteElasticity_ContactProblemSetup

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
      CASE(PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE)
        PROBLEM%CLASS=PROBLEM_ELASTICITY_CLASS
        PROBLEM%TYPE=PROBLEM_FINITE_ELASTICITY_TYPE
        PROBLEM%SUBTYPE=PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE
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

  !>Sets/changes the problem subtype for a finite elasticity contact type .
  SUBROUTINE FiniteElasticity_ContactProblemSubtypeSet(problem,problemSubtype,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: problemSubtype !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("FiniteElasticity_ContactProblemSubtypeSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      SELECT CASE(problemSubtype)
      CASE(PROBLEM_NO_SUBTYPE) !Normal finite elasticity problem subject to contact constraint, no extra solvers required        
        problem%CLASS=PROBLEM_ELASTICITY_CLASS
        problem%TYPE=PROBLEM_FINITE_ELASTICITY_TYPE
        problem%SUBTYPE=PROBLEM_NO_SUBTYPE      
      CASE(PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE)        
        problem%CLASS=PROBLEM_ELASTICITY_CLASS
        problem%TYPE=PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE
        problem%SUBTYPE=PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE
      CASE(PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE)        
        problem%CLASS=PROBLEM_ELASTICITY_CLASS
        problem%TYPE=PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE
        problem%SUBTYPE=PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE
      CASE(PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE)        
        problem%CLASS=PROBLEM_ELASTICITY_CLASS
        problem%TYPE=PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE
        problem%SUBTYPE=PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE  
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(problemSubtype,"*",err,error))// &
          & " is not valid for a finite elasticity contact type of an elasticity problem class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("FiniteElasticity_ContactProblemSubtypeSet")
    RETURN
999 CALL ERRORS("FiniteElasticity_ContactProblemSubtypeSet",err,error)
    CALL EXITS("FiniteElasticity_ContactProblemSubtypeSet")
    RETURN 1
  END SUBROUTINE FiniteElasticity_ContactProblemSubtypeSet

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
    INTEGER(INTG) :: I
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD

    CALL ENTERS("FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
          CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
            & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
            !Call divergence test only if finite element loop: THIS IS NOT A PROPER FIX
            IF(CONTROL_LOOP%SUB_LOOP_INDEX==1) THEN
              CALL SOLVER_NONLINEAR_DIVERGENCE_EXIT(SOLVER,ERR,ERROR,*999)
            ENDIF
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
          CASE DEFAULT
            !Check that solver converged
            CALL SOLVER_NONLINEAR_DIVERGENCE_EXIT(SOLVER,ERR,ERROR,*999)
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
    TYPE(CONTROL_LOOP_TYPE), POINTER :: TIME_LOOP !<A pointer to the control time loop.
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
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE,PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE)
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
               & PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
               & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(EQUATIONS_SET%TYPE==EQUATIONS_SET_FINITE_ELASTICITY_TYPE) THEN
                      TIME_LOOP=>CONTROL_LOOP !Initialise time loop (load increment loop on first)
                      !Move up to find outer time loop
                      DO loop_idx=1,CONTROL_LOOP%CONTROL_LOOP_LEVEL-1
                        IF(ASSOCIATED(TIME_LOOP%PARENT_LOOP)) THEN
                          TIME_LOOP=>TIME_LOOP%PARENT_LOOP
                        ELSE
                          CALL FLAG_ERROR("Could not find a time control loop.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO
                      CURRENT_LOOP_ITERATION=TIME_LOOP%TIME_LOOP%ITERATION_NUMBER
                      OUTPUT_ITERATION_NUMBER=TIME_LOOP%TIME_LOOP%OUTPUT_NUMBER

                      !Write out fields at each timestep
                      IF(TIME_LOOP%TIME_LOOP%CURRENT_TIME<=TIME_LOOP%TIME_LOOP%STOP_TIME) THEN
                        WRITE(OUTPUT_FILE,'("S_TIMESTP_",I4.4)') CURRENT_LOOP_ITERATION
                        FILE=OUTPUT_FILE
                        METHOD="FORTRAN"
                        EXPORT_FIELD=.TRUE.
                        IF(EXPORT_FIELD) THEN
                          IF(OUTPUT_ITERATION_NUMBER/=0.AND.MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                            IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity export fields ...",ERR,ERROR,*999)
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                            ENDIF
                            CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                              & ERR,ERROR,*999)
                            IF(SOLVER%OUTPUT_TYPE>=SOLVER_PROGRESS_OUTPUT) THEN
                              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity all fields exported ...",ERR,ERROR,*999)
                            ENDIF
                            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                          ENDIF
                        ENDIF
                      ENDIF !stop_time


                   ENDIF !EQUATIONS_SET_FINITE_ELASTICITY_TYPE
                  ENDDO !equations_set_idx
                ENDIF !Solver_mapping
              ENDIF !Solver_equations
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
          CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            ! could do this in one line with problem_solver_get but the dependence on problem_routines causes a circular dependence
            CALL CONTROL_LOOP_GET(CONTROL_LOOP,(/1,CONTROL_LOOP_NODE/),CONTROL_LOOP_SOLID,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(CONTROL_LOOP_SOLID%SOLVERS,1,SOLVER_SOLID,ERR,ERROR,*999)
            !--- 3.0 For Standard Elasticity Darcy: Update the boundary conditions of the solid
            CALL FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER_SOLID,ERR,ERROR,*999)
          CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
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

  !>Executes after each loop of a control loop for finite elasticity problems, i.e., after each load increment in a load increment loop
  SUBROUTINE FiniteElasticity_ControlLoadIncrementLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(REGION_TYPE), POINTER :: region
    TYPE(FIELDS_TYPE), POINTER :: fields
    INTEGER(INTG) :: solverIdx,equationsSetIdx,incrementIdx,outputNumber
    LOGICAL :: dirExist
    TYPE(VARYING_STRING) :: fileName,method,directory

    CALL ENTERS("FiniteElasticity_ControlLoadIncrementLoopPostLoop",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(controlLoop%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
        incrementIdx=controlLoop%LOAD_INCREMENT_LOOP%ITERATION_NUMBER
        outputNumber=controlLoop%LOAD_INCREMENT_LOOP%OUTPUT_NUMBER
        IF(outputNumber>0) THEN
          IF(MOD(incrementIdx,outputNumber)==0) THEN
            solvers=>controlLoop%SOLVERS
            IF(ASSOCIATED(solvers)) THEN
              DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
                solver=>solvers%SOLVERS(solverIdx)%PTR
                IF(ASSOCIATED(solver)) THEN
                  solverEquations=>SOLVER%SOLVER_EQUATIONS
                  IF(ASSOCIATED(solverEquations)) THEN
                    solverMapping=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING
                    IF(ASSOCIATED(solverMapping)) THEN
                      DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
                        region=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%PTR%REGION
                        NULLIFY(fields)
                        fields=>region%FIELDS
                        directory="results_load/"
                        INQUIRE(FILE=CHAR(directory),EXIST=dirExist)
                        IF(.NOT.dirExist) THEN
                          CALL SYSTEM(CHAR("mkdir "//directory))
                        ENDIF
                        fileName=directory//"mesh"//TRIM(NUMBER_TO_VSTRING(equationsSetIdx,"*",err,error))// &
                          & "_load"//TRIM(NUMBER_TO_VSTRING(incrementIdx,"*",err,error))
                        method="FORTRAN"
                        CALL FIELD_IO_ELEMENTS_EXPORT(fields,fileName,method,err,error,*999)
                        CALL FIELD_IO_NODES_EXPORT(fields,fileName,method,err,error,*999)
                      ENDDO !equationsSetIdx
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO !solverIdx
            ELSE
              CALL FLAG_ERROR("Control loop solvers is not associated.",err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Control loop is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("FiniteElasticity_ControlLoadIncrementLoopPostLoop")
    RETURN
999 CALL ERRORS("FiniteElasticity_ControlLoadIncrementLoopPostLoop",err,error)
    CALL EXITS("FiniteElasticity_ControlLoadIncrementLoopPostLoop")
    RETURN 1
    
  END SUBROUTINE FiniteElasticity_ControlLoadIncrementLoopPostLoop

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
    INTEGER(INTG) :: solver_matrix_idx,equations_set_idx
    TYPE(SOLVER_TYPE), POINTER :: CELLML_SOLVER
    LOGICAL :: VALID_SUBTYPE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    
    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_FINITE_ELASTICITY_CELLML_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  VALID_SUBTYPE=.FALSE.
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTITUTIVE_LAW_IN_CELLML_EVALUATE_SUBTYPE) THEN
                      VALID_SUBTYPE=.TRUE.
                      !compute the strain field
                      DEPENDENT_FIELD=>EQUATIONS_SET%EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
                      CALL FiniteElasticity_StrainCalculate(EQUATIONS_SET,DEPENDENT_FIELD, &
                        & FIELD_U1_VARIABLE_TYPE,ERR,ERROR,*999)
                      !check for a linked CellML solver 
                      CELLML_SOLVER=>SOLVER%NONLINEAR_SOLVER%NEWTON_SOLVER%CELLML_EVALUATOR_SOLVER
                      IF(ASSOCIATED(CELLML_SOLVER)) THEN
                        !evaluate the constiutive equation in CellML
                        CALL SOLVER_SOLVE(CELLML_SOLVER,ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ENDDO !equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  IF(VALID_SUBTYPE .NEQV. .TRUE.) THEN
                    LOCAL_ERROR="The equations set subtype of number "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR, &
                      & ERROR))//"is not valid for a finite elasticity problem subtype of number "//TRIM( &
                      & NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Solver equations is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              ! do nothing, time values get updated in CONTROL_TIME_LOOP_PRE_LOOP as there might be 
              ! a load increment loop below the time loop, so we don't want to update times here before
              ! every  solve
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
              ! do nothing
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
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
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NDOFS_TO_PRINT
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
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
              & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
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
    REAL(DP), POINTER :: DUMMY_VALUES1(:), CURRENT_PRESSURE_VALUES(:)
    REAL(DP), ALLOCATABLE :: NEW_PRESSURE_VALUES(:)

    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE
    INTEGER(INTG) :: dof_number,GEOMETRY_NUMBER_OF_DOFS,DEPENDENT_NUMBER_OF_DOFS
    INTEGER(INTG) :: NDOFS_TO_PRINT
    INTEGER(INTG) :: loop_idx
    INTEGER(INTG) :: SUBITERATION_NUMBER

    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR,*999)


    NULLIFY( CURRENT_PRESSURE_VALUES, DUMMY_VALUES1 )


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
              CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        SELECT CASE(EQUATIONS_SET%SUBTYPE)
                          CASE(EQUATIONS_SET_INCOMPRESSIBLE_ELASTICITY_DRIVEN_DARCY_SUBTYPE)
                            IF(CONTROL_LOOP%sub_loops(1)%ptr%loop_type==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
                              SUBITERATION_NUMBER=CONTROL_LOOP%sub_loops(1)%ptr%while_loop%iteration_number
                              write(*,*)'SUBITERATION_NUMBER = ',SUBITERATION_NUMBER
                            ELSE
                                CALL FLAG_ERROR("Could not find SUBITERATION_NUMBER.",ERR,ERROR,*999)
                            ENDIF

                            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                              BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                EQUATIONS_MAPPING=>EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING
                                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                                  CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VARIABLE, &
                                    & ERR,ERROR,*999)
                                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                                      & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                      IF(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE)>0) THEN
                                        CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                          & FIELD_PRESSURE_VALUES_SET_TYPE,CURRENT_PRESSURE_VALUES,ERR,ERROR,*999)

                                        IF(DIAGNOSTICS1) THEN
                                          NDOFS_TO_PRINT = SIZE(CURRENT_PRESSURE_VALUES,1)
                                          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                            & NDOFS_TO_PRINT,CURRENT_PRESSURE_VALUES, &
                                            & '(" DEP_FIELD,FIELD_U_VAR_TYPE,FIELD_PRESSURE_VAL_SET_TYPE (before) = ",4(X,E13.6))',&
                                            & '4(4(X,E13.6))',ERR,ERROR,*999)
                                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                            & FIELD_PRESSURE_VALUES_SET_TYPE,CURRENT_PRESSURE_VALUES,ERR,ERROR,*999)
                                        ENDIF

                                        DEPENDENT_NUMBER_OF_DOFS=DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_DELUDELN_VARIABLE_TYPE)% &
                                            & PTR%NUMBER_OF_DOFS

                                        ALLOCATE(NEW_PRESSURE_VALUES(DEPENDENT_NUMBER_OF_DOFS))

                                        !Linear increase of cavity pressure: just a test example prototype
                                        !\todo: general time-dependent boundary condition input method?
                                        ALPHA = ( CURRENT_TIME + TIME_INCREMENT ) / CURRENT_TIME
                                        NEW_PRESSURE_VALUES(1:DEPENDENT_NUMBER_OF_DOFS) = ALPHA * &
                                          & CURRENT_PRESSURE_VALUES(1:DEPENDENT_NUMBER_OF_DOFS)

                                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite Elasticity update pressure BCs", &
                                          & ERR,ERROR,*999)
                                        DO dof_number=1,DEPENDENT_NUMBER_OF_DOFS
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD, &
                                            & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE,dof_number, &
                                            & NEW_PRESSURE_VALUES(dof_number),ERR,ERROR,*999)
                                        ENDDO
                                        CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD, &
                                          & FIELD_DELUDELN_VARIABLE_TYPE, FIELD_PRESSURE_VALUES_SET_TYPE,ERR,ERROR,*999)
                                        CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD, &
                                          & FIELD_DELUDELN_VARIABLE_TYPE, FIELD_PRESSURE_VALUES_SET_TYPE,ERR,ERROR,*999)

                                        DEALLOCATE(NEW_PRESSURE_VALUES)

                                        IF(DIAGNOSTICS1) THEN
                                          NULLIFY( DUMMY_VALUES1 )
                                          CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                            & FIELD_PRESSURE_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                          NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                                          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                            & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                                            & '(" DEP_FIELD,FIELD_U_VAR_TYPE,FIELD_PRESSURE_VAL_SET_TYPE (after) = ",4(X,E13.6))', &
                                            & '4(4(X,E13.6))',ERR,ERROR,*999)
                                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                            & FIELD_PRESSURE_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                        ENDIF
                                        CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                                          & FIELD_PRESSURE_VALUES_SET_TYPE,CURRENT_PRESSURE_VALUES,ERR,ERROR,*999)
                                      ENDIF !Pressure_condition_used
                                    ELSE
                                      CALL FLAG_ERROR("Boundary condition variable is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FLAG_ERROR("Dependent field DelUDelN variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FLAG_ERROR("EQUATIONS_MAPPING is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FLAG_ERROR("Boundary conditions are not associated.",ERR,ERROR,*999)
                              END IF
                            ELSE
                              CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
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
                              BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                EQUATIONS_MAPPING=>EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING
                                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                                  CALL FIELD_VARIABLE_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
                                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                                      & BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                      IF(DIAGNOSTICS1) THEN
                                        NULLIFY( DUMMY_VALUES1 )
                                        CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                                          & FIELD_VALUES_SET_TYPE,DUMMY_VALUES1,ERR,ERROR,*999)
                                        NDOFS_TO_PRINT = SIZE(DUMMY_VALUES1,1)
                                        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT, &
                                          & NDOFS_TO_PRINT,DUMMY_VALUES1, &
                                          & '(" DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE (bef) = ",4(X,E13.6))',&
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
                                          & CONDITION_TYPES(dof_number)
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
                                          & CONDITION_TYPES(dof_number)
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
                                    CALL FLAG_ERROR("Dependent field U variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
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

  !>Evaluates the functions f(J) and f\'(J);
  !>  Eq.(21) in Chapelle, Gerbeau, Sainte-Marie, Vignon-Clementel, Computational Mechanics (2010)
  SUBROUTINE EVALUATE_CHAPELLE_FUNCTION(Jznu,ffact,dfdJfact,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(IN) :: Jznu !<Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    REAL(DP), INTENT(OUT) :: ffact !<f(Jznu) of the INRIA model
    REAL(DP), INTENT(OUT) :: dfdJfact !<dfdJfact = f'(Jznu) of the INRIA model
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    CALL ENTERS("EVALUATE_CHAPELLE_FUNCTION",ERR,ERROR,*999)

!     IF( ABS(Jznu-1.0_DP) > 5.0E-02_DP ) THEN
    IF( ABS(Jznu-1.0_DP) > 1.0E-10_DP ) THEN
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
    REAL(DP), INTENT(IN) :: AZL(3,3) !<C=F\'F
    REAL(DP), INTENT(IN) :: AZU(3,3) !<inverse of AZL
    REAL(DP), INTENT(IN) :: DARCY_MASS_INCREASE !<mass increase
    REAL(DP), INTENT(OUT) :: PIOLA_TENSOR_ADDITION(3,3) !<Addition to the 2nd Piola-Kirchhoff tensor
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: Jznu !<Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    REAL(DP) :: ffact !<f(Jznu) of the INRIA model
    REAL(DP) :: dfdJfact !<dfdJfact = f\'(Jznu) of the INRIA model
    REAL(DP) :: Mfact, bfact, p0fact  !<INRIA constitutive law
    REAL(DP) :: DARCY_VOL_INCREASE, DARCY_RHO_0_F
    INTEGER(INTG) :: i,j

    
    CALL ENTERS("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION",ERR,ERROR,*999)

    !Parameters settings for coupled elasticity Darcy INRIA model:
    CALL GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,ERR,ERROR,*999)

    DARCY_VOL_INCREASE = DARCY_MASS_INCREASE / DARCY_RHO_0_F

    Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    IF( ABS(Jznu) < 1.0E-10_DP ) THEN
      CALL FLAG_ERROR("EVALUATE_CHAPELLE_PIOLA_TENSOR_ADDITION: ABS(Jznu) < 1.0E-10_DP",ERR,ERROR,*999)
    END IF

    CALL EVALUATE_CHAPELLE_FUNCTION(Jznu,ffact,dfdJfact,ERR,ERROR,*999)

    DO i=1,3
      DO j=1,3
!         PIOLA_TENSOR_ADDITION(i,j) = - Mfact * bfact * DARCY_VOL_INCREASE * (ffact + (Jznu - 1.0_DP) * dfdJfact) * Jznu * AZU(i,j) &
!           & + 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * dfdJfact * Jznu * AZU(i,j)
        PIOLA_TENSOR_ADDITION(i,j) = 0.5_DP * Mfact * DARCY_VOL_INCREASE**2.0_DP * Jznu * AZU(i,j)
!         PIOLA_TENSOR_ADDITION(i,j) = 0.0_DP
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

  !>Sets some data for the coupled Darcy / finite-elasticity model
  SUBROUTINE GET_DARCY_FINITE_ELASTICITY_PARAMETERS(DARCY_RHO_0_F,Mfact,bfact,p0fact,ERR,ERROR,*)
  
    !Argument variables
    REAL(DP), INTENT(OUT) :: DARCY_RHO_0_F
    REAL(DP), INTENT(OUT) :: Mfact
    REAL(DP), INTENT(OUT) :: bfact
    REAL(DP), INTENT(OUT) :: p0fact
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    
    CALL ENTERS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS",ERR,ERROR,*999)


!   DARCY_RHO_0_F = 1.0E-03_DP
    DARCY_RHO_0_F = 1.0_DP
!   Mfact = 2.18E05_DP
    Mfact = 2.18E00_DP
    bfact = 1.0_DP
    p0fact = 0.0_DP


    CALL EXITS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS")
    RETURN
999 CALL ERRORS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS",ERR,ERROR)
    CALL EXITS("GET_DARCY_FINITE_ELASTICITY_PARAMETERS")
    RETURN 1
  END SUBROUTINE GET_DARCY_FINITE_ELASTICITY_PARAMETERS

  !
  !================================================================================================================================
  !

  !> Apply load increments to the gravity vector
  SUBROUTINE FINITE_ELASTICITY_LOAD_INCREMENT_APPLY(EQUATIONS_SET,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD
    REAL(DP) :: INCREMENT

    CALL ENTERS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD
        IF(ASSOCIATED(SOURCE_FIELD)) THEN
          IF(MAXIMUM_NUMBER_OF_ITERATIONS>1) THEN
            IF(ITERATION_NUMBER==1) THEN
              !Setup initial values parameter set
              CALL Field_ParameterSetEnsureCreated(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SETS_COPY(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & FIELD_INITIAL_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
            ENDIF
            INCREMENT=REAL(ITERATION_NUMBER)/REAL(MAXIMUM_NUMBER_OF_ITERATIONS)
            CALL FIELD_PARAMETER_SETS_COPY(SOURCE_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INITIAL_VALUES_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,INCREMENT,ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_LOAD_INCREMENT_APPLY")
    RETURN 1

  END SUBROUTINE FINITE_ELASTICITY_LOAD_INCREMENT_APPLY

  !
  !================================================================================================================================
  !
  
END MODULE FINITE_ELASTICITY_ROUTINES
