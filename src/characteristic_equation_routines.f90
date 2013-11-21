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

  PUBLIC Characteristic_EquationsSet_SolutionMethodSet
  PUBLIC Characteristic_EquationsSet_SubtypeSet
  PUBLIC Characteristic_EquationsSet_Setup
  PUBLIC Characteristic_NodalResidualEvaluate
  PUBLIC Characteristic_NodalJacobianEvaluate
  PUBLIC Characteristic_PreSolveUpdateBC

CONTAINS 

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Characteristic equation type of an fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSet_SolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: solutionMethod
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("Characteristic_EquationsSet_SolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)                                
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: equationsSetSubtype
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Characteristic_EquationsSet_SubtypeSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
        equationsSet%CLASS=EQUATIONS_SET_FLUID_MECHANICS_CLASS
        equationsSet%TYPE=EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE
        equationsSet%SUBTYPE=EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: equationsMaterials
    INTEGER(INTG) :: numberOfDimensions,componentIdx,geometricScalingType,geometricMeshComponent,geometricComponentNumber
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents,numberComponentsU2
    INTEGER(INTG) :: independentFieldNumberOfComponents,independentFieldNumberOfVariables,numberComponentsV,numberComponentsU1
    INTEGER(INTG) :: materialsFieldNumberOfVariables,materialsFieldNumberOfComponents1,materialsFieldNumberOfComponents2
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Characteristic_EquationsSet_Setup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
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
          CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
            !Do nothing
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
          CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
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
                !set number of variables to 5 (U,DELUDELN,V,U1,U2)=>(Q,A;dQ,dA;W;pCellML;Pressure)
                dependentFieldNumberOfVariables=5
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & err,error,*999)
                ! set dimension
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                ! set data type
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                ! number of components for U,DELUDELN=2 (Q,A)
                dependentFieldNumberOfComponents=2
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                numberComponentsV=2
                numberComponentsU1=1
                numberComponentsU2=1
                ! set number of components for V
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                 & FIELD_V_VARIABLE_TYPE,numberComponentsV,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                 & FIELD_U1_VARIABLE_TYPE,numberComponentsU1,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                 & FIELD_U2_VARIABLE_TYPE,numberComponentsU2,err,error,*999)
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
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                  & FIELD_U1_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                  & FIELD_U2_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
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
                  ! V; 2 components (W1,W2)
                  DO componentIdx=1,numberComponentsV
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                dependentFieldNumberOfVariables=5 ! U,dUdN,V,U1,U2
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE] &
                  & ,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U2_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                !calculate number of components (Q,A) for U and dUdN
                dependentFieldNumberOfComponents=2
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                numberComponentsV=2
                numberComponentsU1=1
                numberComponentsU2=1
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & numberComponentsV,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & numberComponentsU1,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & numberComponentsU2,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U2_VARIABLE_TYPE,1, &
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
          CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              independentFieldNumberOfComponents=2 ! normalDirection for wave relative to node for W1,W2
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
                !set number of variables to 3 (3 for U,V,U1)
                independentFieldNumberOfVariables=3
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & independentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                ! characteristic normal direction (normalWave) is +/- 1
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                !calculate number of components with one component for each dimension
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U_VARIABLE_TYPE,independentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_V_VARIABLE_TYPE,independentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & FIELD_U1_VARIABLE_TYPE,independentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup
                DO componentIdx=1,independentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U1_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                END DO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  DO componentIdx=1,independentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                independentFieldNumberOfVariables=3
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,independentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_U1_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & independentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & independentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & independentFieldNumberOfComponents,err,error,*999)
              ENDIF    
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
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
          CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
            materialsFieldNumberOfVariables=2 ! U type-7 constant / V type-3 variable
            materialsFieldNumberOfComponents1=7
            materialsFieldNumberOfComponents2=3
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
                    & materialsFieldNumberOfVariables,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD, & 
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents1,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents2,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricComponentNumber,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,geometricComponentNumber,err,error,*999)
                  DO componentIdx=1,materialsFieldNumberOfComponents1 !(MU,RHO,K,As,Re,Fr,St)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO
                  DO componentIdx=1,materialsFieldNumberOfComponents2 !(A0,E,H0)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & componentIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsMaterials%MATERIALS_FIELD,geometricScalingType,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,materialsFieldNumberOfVariables,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE], &
                    & err,error,*999)
                  ! U-variable
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, & 
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents1,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents2,err,error,*999)
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
          CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
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
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(equationsMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
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
                  CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(equationsMatrices, & 
                    & [EquationsMatrix_NodalStructure],err,error,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STORAGE_TYPE_SET(equationsMatrices, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EQUATIONS_MATRICES_NONLINEAR_STRUCTURE_TYPE_SET(equationsMatrices, & 
                    & EquationsMatrix_NodalStructure,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: nodeNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: stiffnessMatrix
    TYPE(FIELD_TYPE), POINTER :: materialsField,dependentField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:),materialsParameters1(:)
    REAL(DP) :: Q_BIF(4),A_BIF(4),A0_PARAM(4),E_PARAM(4),H0_PARAM(4),Beta(4),W(2,4),normalWave(2,4),As,Fr,SUM
    INTEGER(INTG) :: derivativeIdx,versionIdx,versionIdx2,componentIdx,rowIdx,columnIdx,componentIdx2,numberOfVersions,local_ny
    LOGICAL :: updateStiffnessMatrix,updateNonlinearResidual,boundaryNode

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
    CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
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
      stiffnessMatrix%NodalMatrix%matrix=0.0_DP
      nonlinearMatrices%NodalResidual%vector=0.0_DP
      IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%UPDATE_MATRIX
      IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%UPDATE_RESIDUAL

      derivativeIdx=1
      normalWave=0.0_DP
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%numberOfVersions
      boundaryNode=dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
        & TOPOLOGY%NODES%NODES(nodeNumber)%BOUNDARY_NODE

      !Get normal wave direction for nodes
      CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & independentParameters,err,error,*999)
      fieldVariable=>independentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
      DO componentIdx=1,2
        DO versionIdx=1,numberOfVersions
          local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
            & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          normalWave(componentIdx,versionIdx)=independentParameters(local_ny)
        ENDDO
      ENDDO
      CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & independentParameters,err,error,*999)

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(ABS(normalWave(1,1))>0 .OR. ABS(normalWave(2,1))>0) THEN
        IF(.NOT. boundaryNode) THEN

          !Get Material Values at the node
          CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999)
          CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters1,err,error,*999)
          fieldVariable=>materialsField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          local_ny=fieldVariable%COMPONENTS(4)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          As=materialsParameters(local_ny)
          local_ny=fieldVariable%COMPONENTS(6)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          Fr=materialsParameters(local_ny)
          fieldVariable=>materialsField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR
          DO versionIdx=1,numberOfVersions
            local_ny=fieldVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            A0_PARAM(versionIdx)=materialsParameters1(local_ny)
            local_ny=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            E_PARAM(versionIdx)=materialsParameters1(local_ny)
            local_ny=fieldVariable%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            H0_PARAM(versionIdx)=materialsParameters1(local_ny)
            Beta(versionIdx)=(4.0_DP*(PI**0.5_DP)*E_PARAM(versionIdx)*H0_PARAM(versionIdx))/ &
              & (3.0_DP*A0_PARAM(versionIdx))            
          ENDDO
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999)
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters1,err,error,*999)

          !Get current Q & A Values at the node
          CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)
          fieldVariable=>dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DO versionIdx=1,numberOfVersions
            local_ny=fieldVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            Q_BIF(versionIdx)=dependentParameters(local_ny)
            local_ny=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            A_BIF(versionIdx)=dependentParameters(local_ny)
          ENDDO
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)

          !Get extrapolated W for the node
          CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)
          fieldVariable=>dependentField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
              W(componentIdx,versionIdx)=dependentParameters(local_ny)
            ENDDO
          ENDDO
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)

          !!!-- S T I F F N E S S  M A T R I X  --!!!
          IF(updateStiffnessMatrix) THEN
            !Conservation of Mass
            rowIdx=numberOfVersions+1
            columnIdx=0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  stiffnessMatrix%NodalMatrix%matrix(rowIdx,columnIdx)=normalWave(componentIdx,versionIdx)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

          !!!-- N O N L I N E A R   V E C T O R --!!!
          IF(updateNonlinearResidual) THEN
            rowIdx=0
            !Characteristics Equations
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  rowIdx=rowIdx+1
                  nonlinearMatrices%NodalResidual%vector(rowIdx)=(Q_BIF(versionIdx)/A_BIF(versionIdx)) &
                    & +normalWave(componentIdx,versionIdx)*4.0_DP*((Fr*(Beta(versionIdx)))**0.5_DP)* &
                    & (A_BIF(versionIdx)**0.25_DP)-W(componentIdx,versionIdx)
                ENDIF
              ENDDO
            ENDDO
            !Continuity of Total Pressure
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  rowIdx=rowIdx+1
                  IF(versionIdx==1) THEN
                    SUM=0.0_DP
                    DO componentIdx2=1,2
                      DO versionIdx2=1,numberOfVersions
                        SUM=SUM+normalWave(componentIdx2,versionIdx2)*Q_BIF(versionIdx2)
                      ENDDO
                    ENDDO
                    nonlinearMatrices%NodalResidual%vector(rowIdx)=SUM
                  ELSE
                    nonlinearMatrices%NodalResidual%vector(rowIdx)=((A_BIF(1)**0.5_DP)-(Beta(versionIdx)/Beta(1))* &
                      & (A_BIF(versionIdx)**0.5_DP))-(((A0_PARAM(1)/As)**0.5_DP)-(Beta(versionIdx)/Beta(1))* &
                      & ((A0_PARAM(versionIdx)/As)**0.5_DP))+(1.0_DP/(Fr*Beta(1))*0.25_DP*(((Q_BIF(1)/A_BIF(1))**2)- &
                      & ((Q_BIF(versionIdx)/A_BIF(versionIdx))**2)))
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        ENDIF
      ENDIF !Find branch nodes

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

  !>Evaluates the Jacobian nodal matrix for a characteristic equation nodal equations set.
  SUBROUTINE Characteristic_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: nodeNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(FIELD_TYPE), POINTER :: materialsField,dependentField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:),materialsParameters1(:)
    REAL(DP) :: Q_BIF(4),A_BIF(4),A0_PARAM(4),E_PARAM(4),H0_PARAM(4),Beta(4),W(2,4),normalWave(2,4),As,Fr
    INTEGER(INTG) :: numberOfVersions,local_ny,startColumn2
    INTEGER(INTG) :: derivativeIdx,versionIdx,rowIdx,columnIdx,columnIdx2,startRow,endRow,componentIdx
    LOGICAL :: updateJacobianMatrix,boundaryNode

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
    CASE(EQUATIONS_SET_Coupled1D0D_CHARACTERISTIC_SUBTYPE)
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
      jacobianMatrix%NodalJacobian%matrix=0.0_DP
      IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%UPDATE_JACOBIAN

      derivativeIdx=1
      normalWave=0.0_DP
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%numberOfVersions
      boundaryNode=dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
        & TOPOLOGY%NODES%NODES(nodeNumber)%BOUNDARY_NODE

      !Get normal wave direction for nodes
      CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & independentParameters,err,error,*999)
      fieldVariable=>independentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
      DO componentIdx=1,2
        DO versionIdx=1,numberOfVersions
          local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
            & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
          normalWave(componentIdx,versionIdx)=independentParameters(local_ny)
        ENDDO
      ENDDO
      CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & independentParameters,err,error,*999)

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(ABS(normalWave(1,1))>0 .OR. ABS(normalWave(2,1))>0) THEN
        IF(.NOT. boundaryNode) THEN

          !Get Material Values at the node
          CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999)
          CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters1,err,error,*999)
          fieldVariable=>materialsField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          local_ny=fieldVariable%COMPONENTS(4)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          As=materialsParameters(local_ny)
          local_ny=fieldVariable%COMPONENTS(6)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
          Fr=materialsParameters(local_ny)
          fieldVariable=>materialsField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR
          DO versionIdx=1,numberOfVersions
            local_ny=fieldVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
              & NODES(nodeNumber)%DERIVATIVES(1)%VERSIONS(versionIdx)
            A0_PARAM(versionIdx)=materialsParameters1(local_ny)
            local_ny=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% & 
              & NODES(nodeNumber)%DERIVATIVES(1)%VERSIONS(versionIdx)
            E_PARAM(versionIdx)=materialsParameters1(local_ny)
            local_ny=fieldVariable%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
              & NODES(nodeNumber)%DERIVATIVES(1)%VERSIONS(versionIdx)
            H0_PARAM(versionIdx)=materialsParameters1(local_ny)
            Beta(versionIdx) = (4.0_DP*(PI**0.5_DP)*E_PARAM(versionIdx)*H0_PARAM(versionIdx))/ &
              & (3.0_DP*A0_PARAM(versionIdx))            
          ENDDO
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999)
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters1,err,error,*999)

          !Get current Q & A Values at the node
          CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)
          fieldVariable=>dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DO versionIdx=1,numberOfVersions
            local_ny=fieldVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            Q_BIF(versionIdx)=dependentParameters(local_ny)
            local_ny=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
              & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
            A_BIF(versionIdx)=dependentParameters(local_ny)
          ENDDO
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)

          !Get extrapolated W for the node
          CALL FIELD_PARAMETER_SET_DATA_GET(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)
          fieldVariable=>dependentField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
              W(componentIdx,versionIdx)=dependentParameters(local_ny)
            ENDDO
          ENDDO
          CALL FIELD_PARAMETER_SET_DATA_RESTORE(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & dependentParameters,err,error,*999)

          !!!--  J A C O B I A N   M A T R I X  --!!!
          IF(updateJacobianMatrix) THEN
            ! Characteristic equations (dW/dU)
            columnIdx=0
            rowIdx=0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  rowIdx=rowIdx+1
                  jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=1.0_DP/A_BIF(versionIdx)
                ENDIF                
              ENDDO
            ENDDO
            rowIdx=0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  rowIdx=rowIdx+1
                  jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=(-Q_BIF(versionIdx)/(A_BIF(versionIdx)**2)) &
                    & +normalWave(componentIdx,versionIdx)*SQRT((Fr*(Beta(versionIdx))))*((A_BIF(versionIdx)**(-0.75_DP)))
                ENDIF
              ENDDO
            ENDDO

            !Continuity of Total Pressure (dP/dU)
            startRow=numberOfVersions+2
            endRow=numberOfVersions*2
            startColumn2=numberOfVersions+1
            DO rowIdx=startRow,endRow
              columnIdx=1
              columnIdx2=startColumn2
              DO componentIdx=1,2
                DO versionIdx=1,numberOfVersions
                  IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                    IF(columnIdx==1) THEN
                      ! dP/dQ
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=(1.0_DP/(2.0_DP*Fr*Beta(1)))* &
                        & (Q_BIF(1)/(A_BIF(1)**2.0_DP))
                      ! dP/dA
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx2)=1.0_DP/(2.0_DP*SQRT(A_BIF(1))) - &
                        & (1.0_DP/(2.0_DP*Fr*Beta(1)))* &
                        & ((Q_BIF(1)**2.0_DP)/(A_BIF(1)**3.0_DP))
                    ELSE IF(columnIdx2==rowIdx) THEN
                      ! dP/dQ
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=(-1.0_DP/(2.0_DP*Fr*Beta(1)))* &
                        & (Q_BIF(versionIdx)/(A_BIF(versionIdx)**2.0_DP))
                      ! dP/dA
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx2)=-Beta(versionIdx)/Beta(1)* &
                        & (1/(2.0_DP*SQRT(A_BIF(versionIdx)))) + (1.0_DP/(2.0_DP*Fr*Beta(1)))* &
                        & (Q_BIF(versionIdx)**2.0_DP)/(A_BIF(versionIdx)**3.0_DP)
                    ELSE
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,versionIdx)=0.0_DP
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=0.0_DP
                    ENDIF
                    columnIdx=columnIdx+1
                    columnIdx2=columnIdx2+1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO

          ENDIF
        ENDIF
      ENDIF !Find branch nodes

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

  !>Extrapolate W for brnach nodes.
  SUBROUTINE Characteristic_PreSolveUpdateBC(solver,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER 
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: dependentBasis,materialsBasis
    TYPE(DOMAIN_TYPE), POINTER :: dependentDomain,materialsDomain
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER ::  dependentField,materialsField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    REAL(DP), POINTER :: independentParameters(:)
    REAL(DP) :: W(2,4),Q_EX(4),A_EX(4),XI(1),A0_PARAM(4),H0_PARAM(4),E_PARAM(4),Beta(4),As,Fr,normalWave(2,4)
    INTEGER(INTG) :: nodeIdx,versionIdx,derivativeIdx,elementIdx,elementNumber,versionElementNumber(4),local_ny
    INTEGER(INTG) :: elementNodeIdx,elementNodeNumber,elementNodeVersion,numberOfVersions,componentIdx,numberOfLocalNodes

    CALL ENTERS("Characteristic_PreSolveUpdateBC",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      solverEquations=>solver%SOLVER_EQUATIONS
      IF(ASSOCIATED(solverEquations)) THEN
        solverMapping=>solverEquations%SOLVER_MAPPING
        IF(ASSOCIATED(solverMapping)) THEN
          equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
          IF(ASSOCIATED(equationsSet)) THEN
            equations=>equationsSet%EQUATIONS
            IF(ASSOCIATED(equations)) THEN
              !Set General and Specific Pointer
              dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
              independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
              materialsField=>equations%INTERPOLATION%MATERIALS_FIELD
              dependentDomain=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
              materialsDomain=>materialsField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR

              numberOfLocalNodes=dependentDomain%TOPOLOGY%NODES%NUMBER_OF_NODES
              derivativeIdx=1

              !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
              DO nodeIdx=1,numberOfLocalNodes
                numberOfVersions=dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions

                !Get node wave direction
                normalWave=0.0_DP
                CALL FIELD_PARAMETER_SET_DATA_GET(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & independentParameters,err,error,*999)
                fieldVariable=>independentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
                DO componentIdx=1,2
                  DO versionIdx=1,numberOfVersions
                    local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                    normalWave(componentIdx,versionIdx)=independentParameters(local_ny)
                  ENDDO
                ENDDO
                CALL FIELD_PARAMETER_SET_DATA_RESTORE(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & independentParameters,err,error,*999)

                !!!-- F i n d   B r a n c h   N o d e s --!!!
                IF(ABS(normalWave(1,1))>0 .OR. ABS(normalWave(2,1))>0) THEN
                  IF(numberOfVersions>1) THEN

                    DO elementIdx=1,dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS
                      elementNumber=dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%SURROUNDING_ELEMENTS(elementIdx)
                      !Loop over the nodes on this (surrounding) element
                      dependentBasis=>dependentDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                      materialsBasis=>materialsDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                      DO elementNodeIdx=1,dependentBasis%NUMBER_OF_NODES
                        elementNodeNumber=dependentDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)% &
                          & ELEMENT_NODES(elementNodeIdx)
                        !Check that this node is the same as the current iterative node
                        IF(elementNodeNumber==nodeIdx) THEN
                          !Loop over the versions to find the element index that matches the version
                          DO versionIdx=1,numberOfVersions
                            !Version number for the local element node
                            elementNodeVersion=dependentDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%&
                              & elementVersions(1,elementNodeIdx)
                            IF(elementNodeVersion==versionIdx) THEN
                              !Get the element based material parameters for the surrounding elements
                              versionElementNumber(versionIdx)=elementNumber
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,4,As,err,error,*999)
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,6,Fr,err,error,*999)
                              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,A0_PARAM(versionIdx),err,error,*999)
                              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,E_PARAM(versionIdx),err,error,*999)   
                              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,3,H0_PARAM(versionIdx),err,error,*999)            
                              Beta(versionIdx) = (4.0_DP*SQRT(PI)*E_PARAM(versionIdx)*H0_PARAM(versionIdx))/ &
                                & (3.0_DP*A0_PARAM(versionIdx))
                            ENDIF
                          ENDDO
                        ENDIF
                      ENDDO
                    ENDDO

                    !Extrapolate Q & A at branch elements
                    DO componentIdx=1,2
                      DO versionIdx=1,numberOfVersions                         
                        IF(ABS(normalWave(componentIdx,versionIdx))> ZERO_TOLERANCE) THEN
                          IF(normalWave(componentIdx,versionIdx)>ZERO_TOLERANCE) THEN
                            !Inlet
                            XI(1)=0.8_DP
                          ELSE
                            !Outlet
                            XI(1)=0.2_DP
                          ENDIF
                          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE, &
                            & versionElementNumber(versionIdx),EQUATIONS%INTERPOLATION% &
                            & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,EQUATIONS%INTERPOLATION% &
                            & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                          Q_EX(versionIdx)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                            & PTR%VALUES(1,NO_PART_DERIV)
                          A_EX(versionIdx)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                            & PTR%VALUES(2,NO_PART_DERIV)
                        ENDIF
                      ENDDO
                    ENDDO

                    !Extrapolate W at branch elements
                    W(:,:)=0.0_DP
                    DO componentIdx=1,2
                      DO versionIdx=1,numberOfVersions
                        IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                          W(componentIdx,versionIdx)=((Q_EX(versionIdx)/A_EX(versionIdx))+ &
                            & normalWave(componentIdx,versionIdx)*4.0_DP*SQRT((Fr*(Beta(versionIdx))))* &
                            & (A_EX(versionIdx)**(0.25_DP)))
                        ENDIF
                      ENDDO
                    ENDDO

                    !Update W at branch elements
                    fieldVariable=>dependentField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR
                    DO componentIdx=1,2
                      DO versionIdx=1,numberOfVersions
                        local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                          & NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,FIELD_V_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,local_ny,W(componentIdx,versionIdx),ERR,ERROR,*999)
                      ENDDO
                    ENDDO
                  ENDIF
                ENDIF !Find branch nodes
              ENDDO !Loop over nodes

            ELSE
              CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Solver mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("Characteristic_PreSolveUpdateBC")
    RETURN
999 CALL ERRORS("Characteristic_PreSolveUpdateBC",ERR,ERROR)
    CALL EXITS("Characteristic_PreSolveUpdateBC")
    RETURN 1
  END SUBROUTINE Characteristic_PreSolveUpdateBC

  !
  !================================================================================================================================
  !

END MODULE CHARACTERISTIC_EQUATION_ROUTINES
