!> \file  
!> \author David Ladd
!> \brief This module handles the characteristic equation routines. These 
!>  equations are often used in concert with 1D fluid modelling to describe
!>  wave propagation phenomena, which is particularly useful for models of
!>  vascular trees. These equations are also often solved using a discontinuous
!>  nodal solution method, rather than FEM.
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
!> Contributor(s): Soroush Safaei
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

!>This module handles all characteristic equation routines.
MODULE CHARACTERISTIC_EQUATION_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS
  USE MATRIX_VECTOR
  USE MESH_ROUTINES
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES


  IMPLICIT NONE

  PRIVATE

  PUBLIC Characteristic_EquationsSet_SubtypeSet
  PUBLIC Characteristic_EquationsSet_SolutionMethodSet
  PUBLIC Characteristic_EquationsSet_Setup
  PUBLIC Characteristic_NodalJacobianEvaluate
  PUBLIC Characteristic_NodalResidualEvaluate

CONTAINS 

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Characteristic equation type of an fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSet_SolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("Characteristic_EquationsSet_SolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)                                
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          equationsSet%SOLUTION_METHOD=EQUATIONS_SET_NODAL_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
          & " is not valid for a Characteristic equation type of a fluid mechanics equations set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("Characteristic_EquationsSet_SolutionMethodSet")
    RETURN
999 CALL ERRORS("Characteristic_EquationsSet_SolutionMethodSet",err,error)
    CALL EXITS("Characteristic_EquationsSet_SolutionMethodSet")
    RETURN 1
  END SUBROUTINE Characteristic_EquationsSet_SolutionMethodSet

!
!================================================================================================================================
!

  !>Sets/changes the equation subtype for a Characteristic type of a fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSet_SubtypeSet(equationsSet,equationsSetSubtype,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: equationsSetSubtype !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Characteristic_EquationsSet_SubtypeSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
        equationsSet%CLASS=EQUATIONS_SET_FLUID_MECHANICS_CLASS
        equationsSet%TYPE=EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE
        equationsSet%SUBTYPE=EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSetSubtype,"*",err,error))// &
          & " is not valid for a Characteristic fluid type of a fluid mechanics equations set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("Characteristic_EquationsSet_SubtypeSet")
    RETURN
999 CALL ERRORS("Characteristic_EquationsSet_SubtypeSet",err,error)
    CALL EXITS("Characteristic_EquationsSet_SubtypeSet")
    RETURN 1
  END SUBROUTINE Characteristic_EquationsSet_SubtypeSet

!
!================================================================================================================================
!

  !>Sets up the Characteristic equations fluid setup.
  SUBROUTINE Characteristic_EquationsSet_Setup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: equationsMaterials
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    INTEGER(INTG) :: numberOfDimensions,componentIdx
    INTEGER(INTG) :: geometricScalingType,geometricMeshComponent,geometricComponentNumber    
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents
    INTEGER(INTG) :: independentFieldNumberOfComponents,independentFieldNumberOfVariables
    INTEGER(INTG) :: materialFieldNumberOfVariables,materialFieldNumberOfComponents
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Characteristic_EquationsSet_Setup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Characteristic_EquationsSet_SolutionMethodSet(equationsSet, &
                & EQUATIONS_SET_NODAL_SOLUTION_METHOD,err,error,*999)
              equationsSet%SOLUTION_METHOD=EQUATIONS_SET_NODAL_SOLUTION_METHOD
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a characteristic equations set."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
            !Do nothing???
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                !define new created field to be dependent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & geometricDecomposition,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 3 (U,DELUDELN,V)=>(Q,A;dQ,dA;W)
                dependentFieldNumberOfVariables=3
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                ! set dimension
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                ! set data type
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                ! number of components for U,DELUDELN=2 (Q,A)
                dependentFieldNumberOfComponents=2
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                ! number of components for V=1 (W)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,1,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup for U,dUdN
                DO componentIdx=1,dependentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                END DO
                !Default to the geometric interpolation setup for V
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                  & FIELD_V_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                !Specify nodal solution method
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  ! (U, dUdN); 2 components (Q,A)
                  DO componentIdx=1,dependentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                  ! V; 1 component (W)
                  componentIdx=1
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType, &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT

              ELSE 
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                dependentFieldNumberOfVariables=3 ! U,dUdN,V
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                !calculate number of components (Q,A) for U and dUdN
                dependentFieldNumberOfComponents=2
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                ! 1 component (W) for V
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & 1,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD, &
                    & "*",err,error))//" is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
              ENDIF

            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
                & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
                & " is invalid for a characteristic equations set."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              independentFieldNumberOfComponents=1 ! normalDirection for wave relative to node
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                !start field creation with name 'INDEPENDENT_FIELD'
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error, & 
                  & *999)
                !define new created field to be independent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & geometricDecomposition,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet% & 
                  & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 1 (1 for U)
                independentFieldNumberOfVariables=1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & independentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                ! characteristic normal direction (normalWave) is +/- 1
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                !calculate number of components with one component for each dimension
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,independentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup
                DO componentIdx=1,independentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                END DO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  DO componentIdx=1,independentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,geometricScalingType, &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT 

              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & independentFieldNumberOfComponents,err,error,*999)
              ENDIF    
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard characteristic equations set"
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a standard characteristic equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d 
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
            !variable X with has Y components
            materialFieldNumberOfVariables=1 ! 1 U-type container variable w/ 12 components
            materialFieldNumberOfComponents=12
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, & 
                    & equationsSet%MATERIALS%MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(equationsMaterials%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition, & 
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSet%MATERIALS%MATERIALS_FIELD, & 
                    & geometricDecomposition,err,error,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,equationsSet%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(equationsMaterials%MATERIALS_FIELD, & 
                    & materialFieldNumberOfVariables,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD, & 
                    & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,materialFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricComponentNumber,err,error,*999)
                  ! Constant Materials
                  DO componentIdx=1,8 !(mu,rho,K,Bs,As,Re,Fr,St)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  ENDDO
                  ! Element-based Materials
                  DO componentIdx=9,12 !(A0,Beta,E,H0)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, & 
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsMaterials%MATERIALS_FIELD,geometricScalingType,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, & 
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialFieldNumberOfComponents,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              END IF
              !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(equationsMaterials%MATERIALS_FIELD,err,error,*999)
                  ! Should be initialized from example file
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*", & 
                & err,error))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*", & 
                & err,error))//" is invalid for characteristic equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN              
                IF(equationsMaterials%MATERIALS_FINISHED) THEN
                  CALL EQUATIONS_CREATE_START(equationsSet,equations,err,error,*999)
                  CALL EQUATIONS_LINEARITY_TYPE_SET(equations,EQUATIONS_NONLINEAR,err,error,*999)
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(equations,EQUATIONS_STATIC,err,error,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set materials has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(equationsSet,equations,err,error,*999)
                CALL EQUATIONS_CREATE_FINISH(equations,err,error,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(equations,equationsMapping,err,error,*999)
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(equationsMapping,1,err,error,*999)
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(equationsMapping,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(equationsMapping,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & err,error,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(equationsMapping,err,error,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(equations,equationsMatrices,err,error,*999)
                ! Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(equationsMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & err,error,*999)
                SELECT CASE(equations%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(equationsMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(equationsMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(equationsMatrices, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(equationsMatrices, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(equationsMatrices, & 
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(equationsMatrices, & 
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(equations%SPARSITY_TYPE,"*",err,error))//" is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
                CALL EQUATIONS_MATRICES_CREATE_FINISH(equationsMatrices,err,error,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FLAG_ERROR("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FLAG_ERROR(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a characteristics equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristics equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a characteristics equation set."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
          & " does not equal a characteristics equation set."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("Characteristic_EquationsSet_Setup")
    RETURN
999 CALL ERRORS("Characteristic_EquationsSet_Setup",err,error)
    CALL EXITS("Characteristic_EquationsSet_Setup")
    RETURN 1
  END SUBROUTINE Characteristic_EquationsSet_Setup

  !
  !================================================================================================================================
  !

  !>Evaluates the residual nodal stiffness matrices and RHS for a characteristic equation nodal equations set.
  SUBROUTINE Characteristic_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to perform the finite nodal calculations on
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The nodal number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: stiffnessMatrix
    TYPE(FIELD_TYPE), POINTER :: materialsField
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_TYPE), POINTER :: independentField
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    REAL(DP), POINTER :: dependentParameters(:)
    REAL(DP), POINTER :: independentParameters(:)
    REAL(DP), POINTER :: materialsParameters(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    INTEGER(INTG) :: numberOfVersions,local_ny,variableType
    INTEGER(INTG) :: derivativeIdx,versionIdx,rowIdx,parameterIdx
    INTEGER(INTG) :: elementNumber,elementIdx,elementNodeIdx,elementNodeNumber,elementNodeVersion,versionElementNumber(3)
    REAL(DP) :: Q_BIF(3),A_BIF(3),A0_PARAM(3),Beta(3),W(3),normalWave(3)
    REAL(DP) :: Bs,As,Fr
    LOGICAL :: updateStiffnessMatrix, updateRhsVector,updateNonlinearResidual
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Characteristic_NodalResidualEvaluate",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(nonlinearMapping)
    NULLIFY(nonlinearMatrices)
    NULLIFY(stiffnessMatrix)
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(materialsField)
    NULLIFY(domain)
    NULLIFY(domainNodes)
    NULLIFY(dependentParameters)
    NULLIFY(independentParameters)
    NULLIFY(materialsParameters)
    NULLIFY(fieldVariable)

    updateStiffnessMatrix=.FALSE.
    updateRhsVector=.FALSE.
    updateNonlinearResidual=.FALSE.

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
      IF(ASSOCIATED(equations)) THEN
        dependentField=>equations%EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          domain=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(domain)) THEN
            domainNodes=>domain%TOPOLOGY%NODES
          ELSE
            CALL FLAG_ERROR("Domain is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent Field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF


    SELECT CASE(equationsSet%SUBTYPE)
    CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
      !Set General and Specific Pointers
      independentField=>equations%INTERPOLATION%INDEPENDENT_FIELD
      materialsField=>equations%INTERPOLATION%MATERIALS_FIELD
      equationsMatrices=>equations%EQUATIONS_MATRICES
      equationsMapping=>equations%EQUATIONS_MAPPING
      linearMatrices=>equationsMatrices%LINEAR_MATRICES
      nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
      stiffnessMatrix=>linearMatrices%MATRICES(1)%PTR
      linearMapping=>equationsMapping%LINEAR_MAPPING
      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
      !Default matrix/vector to 0 and check if update called for
      stiffnessMatrix%NODAL_MATRIX%MATRIX=0.0_DP
      nonlinearMatrices%NODAL_RESIDUAL%VECTOR=0.0_DP
      IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%UPDATE_MATRIX
      IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%UPDATE_RESIDUAL
      ! Get the number of versions at this node
      derivativeIdx=1
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%NUMBER_OF_VERSIONS

      ! Check if this is a standard node or a branching/coupled node (branch/coupled has >1 version)
      ! (if standard, the returned nodal element vector from this routine will be 0)
      IF(numberOfVersions>1) THEN


        !!!--  M A T E R I A L S   P A R A M E T E R S  --!!!
        !Material Values at the node - material field variable 1, components 1-12, type U
        ! First 8 components constant, 9-12 element-based (so will vary with version#)

        ! Element-based material parameters for coupled problems
        IF(dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
           & NODES%NODES(nodeNumber)%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
          elementNumber=dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
           & NODES%NODES(nodeNumber)%SURROUNDING_ELEMENTS(1)
          DO versionIdx=1,numberOfVersions
            parameterIdx=9
            CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
             & elementNumber,parameterIdx,A0_PARAM(versionIdx),err,error,*999)                
            parameterIdx=10
            CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
             & elementNumber,parameterIdx,Beta(versionIdx),err,error,*999)                
          ENDDO
        ! For branching flows, this gets the elements converging on the node and associates them 
        ! with their specific version number: versionElementNumber(version)
        ELSE
          !Loop over the elements surrounding the current node
          DO elementIdx=1,dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
            & NODES%NODES(nodeNumber)%NUMBER_OF_SURROUNDING_ELEMENTS
            ! get the user number for the current element
            elementNumber=dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
             & NODES%NODES(nodeNumber)%SURROUNDING_ELEMENTS(elementIdx)
            ! loop over the nodes on this (surrounding) element
            DO elementNodeIdx=1,equations%INTERPOLATION%GEOMETRIC_FIELD%DECOMPOSITION%MESH%TOPOLOGY(1)% &
              & PTR%ELEMENTS%ELEMENTS(elementNumber)%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
              !get the node number for this node
              elementNodeNumber=equations%INTERPOLATION%GEOMETRIC_FIELD%DECOMPOSITION%MESH%TOPOLOGY(1)% &
                & PTR%ELEMENTS%ELEMENTS(elementNumber)%MESH_ELEMENT_NODES(elementNodeIdx)
              ! check that this node is the same as the current iterative node
              IF(elementNodeNumber==nodeNumber) THEN
                ! Loop over the versions to find the element index that matches the version
                DO versionIdx=1,numberOfVersions
                  ! the version number for the local element node
                  elementNodeVersion=equations%INTERPOLATION%GEOMETRIC_FIELD%DECOMPOSITION%MESH%TOPOLOGY(1)% &
                    & PTR%ELEMENTS%ELEMENTS(elementNumber)% &
                    & USER_ELEMENT_NODE_VERSIONS(derivativeIdx,elementNodeIdx)
                  IF(elementNodeVersion==versionIdx) THEN
                    !Get the element based material parameters for the surrounding elements
                    versionElementNumber(versionIdx)=elementNumber
                    parameterIdx=9
                    CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                     & elementNumber,parameterIdx,A0_PARAM(versionIdx),err,error,*999)                
                    parameterIdx=10
                    CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                     & elementNumber,parameterIdx,Beta(versionIdx),err,error,*999)                
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & materialsParameters,err,error,*999)
        variableType=materialsField%VARIABLES(1)%VARIABLE_TYPE ! variables U
        fieldVariable=>materialsField%VARIABLE_TYPE_MAP(variableType)%PTR
        ! Constant material parameters
        local_ny=fieldVariable%COMPONENTS(4)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
        Bs=materialsParameters(local_ny)
        local_ny=fieldVariable%COMPONENTS(5)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
        As=materialsParameters(local_ny)
        local_ny=fieldVariable%COMPONENTS(7)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
        Fr=materialsParameters(local_ny)
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & materialsParameters,err,error,*999)

        !!!--  D E P E N D E N T   P A R A M E T E R S   ( Q , A ) --!!!
        !Current Q and A Values at the nodes - dependent field variable 1, components 1-2, type U
        CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)
        variableType=dependentField%VARIABLES(1)%VARIABLE_TYPE 
        fieldVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%PTR
        DO versionIdx=1,numberOfVersions
          parameterIdx=1  
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          Q_BIF(versionIdx)=dependentParameters(local_ny)
          parameterIdx=2
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          A_BIF(versionIdx)=dependentParameters(local_ny)
        ENDDO
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)

        !!!-- E X T R A P O L A T E D   C H A R A C T E R I S T I C S  ( W )  --!!!
        CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)
        ! characteristic waves stored in dependent variable 3, component 1, type V
        variableType=dependentField%VARIABLES(3)%VARIABLE_TYPE
        fieldVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%PTR
        parameterIdx=1
        DO versionIdx=1,numberOfVersions
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          W(versionIdx)=dependentParameters(local_ny)
        ENDDO
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)

        !!!-- W a v e   D i r e c t i o n  ( nW ) --!!!
        CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & independentParameters,err,error,*999)
        ! Relative to node, +/- 1, stored in independent field variable 1, component 1, type U
        variableType=independentField%VARIABLES(1)%VARIABLE_TYPE
        fieldVariable=>independentField%VARIABLE_TYPE_MAP(variableType)%PTR
        parameterIdx=1
        DO versionIdx=1,numberOfVersions
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          normalWave(versionIdx)=independentParameters(local_ny)
        ENDDO
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & independentParameters,err,error,*999)

        !!!-- S T I F F N E S S  M A T R I X  --!!!
        IF(updateStiffnessMatrix) THEN
          !Conservation of Mass
          rowIdx=numberOfVersions+1
          DO versionIdx=1,numberOfVersions
            stiffnessMatrix%NODAL_MATRIX%MATRIX(rowIdx,versionIdx)=normalWave(versionIdx)
          ENDDO
        END IF

        !!!-- N O N L I N E A R   V E C T O R --!!!
        IF(updateNonlinearResidual) THEN
          rowIdx=0
          !Characteristics Equations
          DO versionIdx=1,numberOfVersions
            rowIdx=rowIdx + 1
            nonlinearMatrices%NODAL_RESIDUAL%VECTOR(rowIdx)=(Q_BIF(versionIdx)/A_BIF(versionIdx)) &
              & + normalWave(versionIdx)*4.0_DP*((Fr*(Beta(versionIdx)/Bs))**0.5_DP)*(A_BIF(versionIdx)**0.25_DP)- &
              & W(versionIdx)
          ENDDO
          !Continuity of Total Pressure (relative to the first component)
          DO versionIdx=1,numberOfVersions
            rowIdx=rowIdx + 1
            ! will be 0 when versionIdx = 1
            nonlinearMatrices%NODAL_RESIDUAL%VECTOR(rowIdx)=((A_BIF(1)**0.5_DP)-(Beta(versionIdx)/Beta(1))* &
              & (A_BIF(versionIdx)**0.5_DP))-(((A0_PARAM(1)/As)**0.5_DP)-(Beta(versionIdx)/Beta(1))* &
              & ((A0_PARAM(versionIdx)/As)**0.5_DP))+(Bs/(Fr*Beta(1))*0.25_DP*(((Q_BIF(1)/A_BIF(1))**2)- &
              & ((Q_BIF(versionIdx)/A_BIF(versionIdx))**2)))
          ENDDO
        ENDIF
      ENDIF !check for normal node or coupled/branching

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
        & " is not valid for a characteristic equation type of a fluid mechanics equations set class."
      CALL FLAG_ERROR(localError,err,error,*999)
    END SELECT

    
    CALL EXITS("Characteristic_NodalResidualEvaluate")
    RETURN
999 CALL ERRORS("Characteristic_NodalResidualEvaluate",err,error)
    CALL EXITS("Characteristic_NodalResidualEvaluate")
    RETURN 1
  END SUBROUTINE Characteristic_NodalResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian nodal matrix for a Navier-Stokes equation finite nodal equations set.
  SUBROUTINE Characteristic_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to perform the finite nodal calculations on
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The nodal number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(FIELD_TYPE), POINTER :: materialsField
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_TYPE), POINTER :: independentField
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    REAL(DP), POINTER :: dependentParameters(:)
    REAL(DP), POINTER :: independentParameters(:)
    REAL(DP), POINTER :: materialsParameters(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    INTEGER(INTG) :: numberOfVersions,local_ny,variableType
    INTEGER(INTG) :: derivativeIdx,versionIdx,rowIdx,parameterIdx,baseIdx,columnIdx
    REAL(DP) :: Q_BIF(3),A_BIF(3),A0_PARAM(3),Beta(3),W(3),normalWave(3)
    INTEGER(INTG) :: elementNumber,elementIdx,elementNodeIdx,elementNodeNumber,elementNodeVersion,versionElementNumber(3)
    REAL(DP) :: Bs,As,Fr
    LOGICAL :: updateJacobianMatrix
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Characteristic_NodalJacobianEvaluate",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(nonlinearMapping)
    NULLIFY(nonlinearMatrices)
    NULLIFY(jacobianMatrix)
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(materialsField)
    NULLIFY(domain)
    NULLIFY(domainNodes)
    NULLIFY(dependentParameters)
    NULLIFY(independentParameters)
    NULLIFY(materialsParameters)
    NULLIFY(fieldVariable)

    updateJacobianMatrix=.FALSE.

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
      IF(ASSOCIATED(equations)) THEN
        dependentField=>equations%EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          domain=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(domain)) THEN
            domainNodes=>domain%TOPOLOGY%NODES
          ELSE
            CALL FLAG_ERROR("Domain is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent Field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    SELECT CASE(equationsSet%SUBTYPE)
    CASE(EQUATIONS_SET_STATIC_CHARACTERISTIC_SUBTYPE)
      !Set General and Specific Pointers
      independentField=>equations%EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
      materialsField=>equations%INTERPOLATION%MATERIALS_FIELD
      equationsMatrices=>equations%EQUATIONS_MATRICES
      equationsMapping=>equations%EQUATIONS_MAPPING
      linearMatrices=>equationsMatrices%LINEAR_MATRICES
      nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
      linearMapping=>equationsMapping%LINEAR_MAPPING
      jacobianMatrix=>nonlinearMatrices%JACOBIANS(1)%PTR
      !Default matrix/vector to 0 and check if update called for
      jacobianMatrix%NODAL_JACOBIAN%MATRIX=0.0_DP
      IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%UPDATE_JACOBIAN
      ! Set derivative to 1 and get the number of versions at this node
      derivativeIdx=1
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%NUMBER_OF_VERSIONS
      ! Check if this is a standard node or a branching/coupled node (branch/coupled has >1 version)
      ! (if standard, the returned nodal nodal vector from this routine will be 0)
      IF(numberOfVersions>1) THEN


        !!!--  M A T E R I A L S   P A R A M E T E R S  --!!!
        !Material Values at the node - material field variable 1, components 1-12, type U
        ! First 8 components constant, 9-12 element-based (so will vary with version#)

        ! Element-based material parameters for coupled problems
        IF(dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
           & NODES%NODES(nodeNumber)%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
          elementNumber=dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
           & NODES%NODES(nodeNumber)%SURROUNDING_ELEMENTS(1)
          DO versionIdx=1,numberOfVersions
            parameterIdx=9
            CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
             & elementNumber,parameterIdx,A0_PARAM(versionIdx),err,error,*999)                
            parameterIdx=10
            CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
             & elementNumber,parameterIdx,Beta(versionIdx),err,error,*999)                
          ENDDO
        ! For branching flows, this gets the elements converging on the node and associates them 
        ! with their specific version number: versionElementNumber(version)
        ELSE
          !Loop over the elements surrounding the current node
          DO elementIdx=1,dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
            & NODES%NODES(nodeNumber)%NUMBER_OF_SURROUNDING_ELEMENTS
            ! get the user number for the current element
            elementNumber=dependentField%DECOMPOSITION%DOMAIN(1)%PTR%TOPOLOGY% &
             & NODES%NODES(nodeNumber)%SURROUNDING_ELEMENTS(elementIdx)
            ! loop over the nodes on this (surrounding) element
            DO elementNodeIdx=1,equations%INTERPOLATION%GEOMETRIC_FIELD%DECOMPOSITION%MESH%TOPOLOGY(1)% &
              & PTR%ELEMENTS%ELEMENTS(elementNumber)%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
              !get the node number for this node
              elementNodeNumber=equations%INTERPOLATION%GEOMETRIC_FIELD%DECOMPOSITION%MESH%TOPOLOGY(1)% &
                & PTR%ELEMENTS%ELEMENTS(elementNumber)%MESH_ELEMENT_NODES(elementNodeIdx)
              ! check that this node is the same as the current iterative node
              IF(elementNodeNumber==nodeNumber) THEN
                ! Loop over the versions to find the element index that matches the version
                DO versionIdx=1,numberOfVersions
                  ! the version number for the local element node
                  elementNodeVersion=equations%INTERPOLATION%GEOMETRIC_FIELD%DECOMPOSITION%MESH%TOPOLOGY(1)% &
                    & PTR%ELEMENTS%ELEMENTS(elementNumber)% &
                    & USER_ELEMENT_NODE_VERSIONS(derivativeIdx,elementNodeIdx)
                  IF(elementNodeVersion==versionIdx) THEN
                    !Get the element based material parameters for the surrounding elements
                    versionElementNumber(versionIdx)=elementNumber
                    parameterIdx=9
                    CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                     & elementNumber,parameterIdx,A0_PARAM(versionIdx),err,error,*999)                
                    parameterIdx=10
                    CALL FIELD_PARAMETER_SET_GET_ELEMENT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                     & elementNumber,parameterIdx,Beta(versionIdx),err,error,*999)                
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & materialsParameters,err,error,*999)
        variableType=materialsField%VARIABLES(1)%VARIABLE_TYPE ! variables U
        fieldVariable=>materialsField%VARIABLE_TYPE_MAP(variableType)%PTR
        ! Constant material parameters
        local_ny=fieldVariable%COMPONENTS(4)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
        Bs=materialsParameters(local_ny)
        local_ny=fieldVariable%COMPONENTS(5)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
        As=materialsParameters(local_ny)
        local_ny=fieldVariable%COMPONENTS(7)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
        Fr=materialsParameters(local_ny)
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & materialsParameters,err,error,*999)


        !!!--  D E P E N D E N T   P A R A M E T E R S   ( Q , A ) --!!!
        !Current Q and A Values at the nodes - dependent field variable 1, components 1-2, type U
        CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)
        variableType=dependentField%VARIABLES(1)%VARIABLE_TYPE 
        fieldVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%PTR
        DO versionIdx=1,numberOfVersions
          parameterIdx=1  
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          Q_BIF(versionIdx)=dependentParameters(local_ny)
          parameterIdx=2
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          A_BIF(versionIdx)=dependentParameters(local_ny)
        ENDDO
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)


        !!!-- E X T R A P O L A T E D   C H A R A C T E R I S T I C S  ( W )  --!!!
        CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)
        ! characteristic waves stored in dependent variable 3, component 1, type V
        variableType=dependentField%VARIABLES(3)%VARIABLE_TYPE
        fieldVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%PTR
        parameterIdx=1
        DO versionIdx=1,numberOfVersions
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          W(versionIdx)=dependentParameters(local_ny)
        ENDDO
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & dependentParameters,err,error,*999)


        !!!-- W a v e   D i r e c t i o n  ( nW ) --!!!
        CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & independentParameters,err,error,*999)
        ! Relative to node, +/- 1, stored in independent field variable 1, component 1, type U
        variableType=independentField%VARIABLES(1)%VARIABLE_TYPE
        fieldVariable=>independentField%VARIABLE_TYPE_MAP(variableType)%PTR
        parameterIdx=1
        DO versionIdx=1,numberOfVersions
          local_ny=fieldVariable%COMPONENTS(parameterIdx)%PARAM_TO_DOF_MAP% &
           & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          normalWave(versionIdx)=independentParameters(local_ny)
        ENDDO
        CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & independentParameters,err,error,*999)


        !!!--  J A C O B I A N   M A T R I X  --!!!
        IF(updateJacobianMatrix) THEN
          ! Characteristic equations
          columnIdx=0
          DO versionIdx=1,numberOfVersions
            columnIdx=columnIdx+1
            jacobianMatrix%NODAL_JACOBIAN%MATRIX(versionIdx,columnIdx)=1.0_DP/A_BIF(versionIdx)                
          ENDDO
          DO versionIdx=1,numberOfVersions
            columnIdx=columnIdx+1
            jacobianMatrix%NODAL_JACOBIAN%MATRIX(versionIdx,columnIdx)=(-Q_BIF(versionIdx)/(A_BIF(versionIdx)**2)) &
             & +normalWave(versionIdx)*4.0_DP*((Fr*(Beta(versionIdx)/Bs))**(0.5_DP))* &
             & (0.25_DP*(A_BIF(versionIdx)**(-(0.75_DP))))
          ENDDO
          rowIdx=numberOfVersions+1
          ! TODO: check this matrix construction against the theory again
          jacobianMatrix%NODAL_JACOBIAN%MATRIX(rowIdx,:)=0.0_DP !4

          DO baseIdx=2,numberOfVersions ! 2,3
            rowIdx=baseIdx+numberOfVersions ! 5,6
            DO versionIdx=1,numberOfVersions ! 1,2,3
              columnIdx=versionIdx+numberOfVersions 
              IF(versionIdx==1 .OR. versionIdx==baseIdx) THEN
                !Continuity of Total Pressure (dU/dQ)
                jacobianMatrix%NODAL_JACOBIAN%MATRIX(rowIdx,versionIdx)=(Bs/(Fr*Beta(1)))*0.25_DP* &
                 & (normalWave(versionIdx)*2.0_DP*Q_BIF(versionIdx)/(A_BIF(1)**2))

                !columnIdx=versionIdx+numberOfVersions ! 5,4;5,5;6,4;6,6
                !Continuity of Total Pressure (dU/dA) 
                IF(versionIdx==1) THEN
                  jacobianMatrix%NODAL_JACOBIAN%MATRIX(rowIdx,columnIdx)=(0.5_DP*A_BIF(1)**(-(0.5_DP)))+ &
                    & ((Bs/(Fr*Beta(1)))*0.25_DP*(-2.0_DP*(Q_BIF(1)**2)/(A_BIF(1)**3)))
                ELSE
                  jacobianMatrix%NODAL_JACOBIAN%MATRIX(rowIdx,columnIdx)=(-(Beta(versionIdx)/Beta(1))* &
                    & (0.5_DP*A_BIF(versionIdx)**(-(0.5_DP))))+((Bs/(Fr*Beta(1)))*0.25_DP* &
                    & (-2.0_DP*normalWave(versionIdx)*(Q_BIF(versionIdx)**2)/(A_BIF(versionIdx)**3)))
                ENDIF

              ELSE
                jacobianMatrix%NODAL_JACOBIAN%MATRIX(rowIdx,versionIdx)=0.0_DP
                jacobianMatrix%NODAL_JACOBIAN%MATRIX(rowIdx,columnIdx)=0.0_DP
              ENDIF
            ENDDO !versionIdx
          ENDDO !baseidx
        ENDIF !update jacobian
      ENDIF !check for normal or coupled/branch node
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
        & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
      CALL FLAG_ERROR(localError,err,error,*999)
    END SELECT
       
    CALL EXITS("Characteristic_NodalJacobianEvaluate")
    RETURN
999 CALL ERRORS("Characteristic_NodalJacobianEvaluate",err,error)
    CALL EXITS("Characteristic_NodalJacobianEvaluate")
    RETURN 1
  END SUBROUTINE Characteristic_NodalJacobianEvaluate

  !
  !================================================================================================================================
  !

END MODULE CHARACTERISTIC_EQUATION_ROUTINES
