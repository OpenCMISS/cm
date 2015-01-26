!> \file  
!> \author David Ladd
!> \brief This module handles the Stree equation routines. These 
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

!>This module handles all Stree equation routines.
MODULE Stree_EQUATION_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE CMISS_MPI  
  USE COMP_ENVIRONMENT
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

  PUBLIC Stree_EquationsSet_SolutionMethodSet
  PUBLIC Stree_EquationsSet_SubtypeSet
  PUBLIC Stree_EquationsSet_Setup
  PUBLIC Stree_FINITE_ELEMENT_CALCULATE
  PUBLIC Stree_PRE_SOLVE

CONTAINS 

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Stree equation type of an fluid mechanics equations set class.
  SUBROUTINE Stree_EquationsSet_SolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: solutionMethod
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    CALL ENTERS("Stree_EquationsSet_SolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)                                
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          CALL FLAG_ERROR("Not implemented.",err,error,*999)
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
          & " is not valid for a Stree equation type of a fluid mechanics equations set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF
       
    CALL EXITS("Stree_EquationsSet_SolutionMethodSet")
    RETURN
999 CALL ERRORS("Stree_EquationsSet_SolutionMethodSet",err,error)
    CALL EXITS("Stree_EquationsSet_SolutionMethodSet")
    RETURN 1
  END SUBROUTINE Stree_EquationsSet_SolutionMethodSet

!
!================================================================================================================================
!

  !>Sets/changes the equation subtype for a Stree type of a fluid mechanics equations set class.
  SUBROUTINE Stree_EquationsSet_SubtypeSet(equationsSet,equationsSetSubtype,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: equationsSetSubtype
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Stree_EquationsSet_SubtypeSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
        equationsSet%CLASS=EQUATIONS_SET_FLUID_MECHANICS_CLASS
        equationsSet%TYPE=EQUATIONS_SET_Stree_EQUATION_TYPE
        equationsSet%SUBTYPE=EQUATIONS_SET_Stree1D0D_SUBTYPE
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSetSubtype,"*",err,error))// &
          & " is not valid for a Stree fluid type of a fluid mechanics equations set class."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("Stree_EquationsSet_SubtypeSet")
    RETURN
999 CALL ERRORS("Stree_EquationsSet_SubtypeSet",err,error)
    CALL EXITS("Stree_EquationsSet_SubtypeSet")
    RETURN 1
  END SUBROUTINE Stree_EquationsSet_SubtypeSet

!
!================================================================================================================================
!

  !>Sets up the Stree equations fluid setup.
  SUBROUTINE Stree_EquationsSet_Setup(equationsSet,equationsSetSetup,err,error,*)

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
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    INTEGER(INTG) :: I,numberOfDimensions,componentIdx,geometricScalingType,geometricMeshComponent,geometricComponentNumber
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents
    INTEGER(INTG) :: materialsFieldNumberOfVariables,materialsFieldNumberOfComponents
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("Stree_EquationsSet_Setup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Stree_EquationsSet_SolutionMethodSet(equationsSet, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
              equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(equationsSetField,"Equations Set Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSetField,FIELD_GENERAL_TYPE,&
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(equationsSetField,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSetField,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(equationsSetField,FIELD_U_VARIABLE_TYPE,"Stree",ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a Stree equations set."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a Stree equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSetField,geometricDecomposition,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSetField,equationsSet%GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,ERR,ERROR,*999)                
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & geometricScalingType,ERR,ERROR,*999)
              ELSE
                !Do nothing
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Stree equation."
              CALL FLAG_ERROR(localError,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a Stree equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
            dependentFieldNumberOfVariables=2
            dependentFieldNumberOfComponents=1
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
                !set number of variables
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                !set dimension
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                !set data type
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                ! number of components for 1
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup for U,dUdN
                DO I=1,dependentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,I,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,geometricMeshComponent,err,error,*999)
                ENDDO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,dependentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NUMBER_TO_VSTRING(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FLAG_ERROR(localError,err,error,*999)
                END SELECT
              ELSE 
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,dependentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
                & " is invalid for a Stree equations set."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a Stree equations set."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d 
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
            materialsFieldNumberOfVariables=2
            materialsFieldNumberOfComponents=1
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
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricComponentNumber,err,error,*999)
                  DO I=1,materialsFieldNumberOfComponents
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & I,geometricComponentNumber,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & I,geometricComponentNumber,err,error,*999)
                  ENDDO
                  DO I=1,materialsFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsMaterials%MATERIALS_FIELD,geometricScalingType,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,materialsFieldNumberOfVariables,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & materialsFieldNumberOfComponents,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              END IF
              !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(equationsMaterials%MATERIALS_FIELD,err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%ACTION_TYPE,"*", & 
                & err,error))//" for a setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*", & 
                & err,error))//" is invalid for Stree equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a Stree equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN              
                IF(equationsMaterials%MATERIALS_FINISHED) THEN
                  CALL EQUATIONS_CREATE_START(equationsSet,equations,err,error,*999)
                  CALL EQUATIONS_LINEARITY_TYPE_SET(equations,EQUATIONS_LINEAR,err,error,*999)
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(equations,EQUATIONS_STATIC,err,error,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set materials has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(equationsSet,equations,err,error,*999)
                CALL EQUATIONS_CREATE_FINISH(equations,err,error,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(equations,equationsMapping,err,error,*999)
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(equationsMapping,1,err,error,*999)
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(equationsMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(equationsMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(equationsMapping,err,error,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(equations,equationsMatrices,err,error,*999)
                SELECT CASE(equations%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(equationsMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(equationsMatrices, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(equationsMatrices, & 
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
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
                & " is invalid for a Strees equation."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a Strees equation."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a Strees equation set."
          CALL FLAG_ERROR(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
          & " does not equal a Strees equation set."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("Stree_EquationsSet_Setup")
    RETURN
999 CALL ERRORS("Stree_EquationsSet_Setup",err,error)
    CALL EXITS("Stree_EquationsSet_Setup")
    RETURN 1
  END SUBROUTINE Stree_EquationsSet_Setup

  !
  !================================================================================================================================
  !

  !>Evaluates the residual nodal stiffness matrices and RHS for a Stree equation nodal equations set.
  SUBROUTINE STREE_FINITE_ELEMENT_CALCULATE(equationsSet,nodeNumber,err,error,*)

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
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: stiffnessMatrix
    TYPE(FIELD_TYPE), POINTER :: materialsField,dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    REAL(DP), POINTER :: dependentParameters(:),materialsParameters(:)

    CALL ENTERS("STREE_FINITE_ELEMENT_CALCULATE",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(stiffnessMatrix)
    NULLIFY(dependentField)
    NULLIFY(materialsField)
    NULLIFY(domain)
    NULLIFY(domainNodes)
    NULLIFY(dependentParameters)
    NULLIFY(materialsParameters)
    NULLIFY(fieldVariable)

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
    CASE(EQUATIONS_SET_Stree1D0D_SUBTYPE)
      !Set General and Specific Pointers
      equationsMatrices=>equations%EQUATIONS_MATRICES
      equationsMapping=>equations%EQUATIONS_MAPPING
      linearMatrices=>equationsMatrices%LINEAR_MATRICES
      stiffnessMatrix=>linearMatrices%MATRICES(1)%PTR
      linearMapping=>equationsMapping%LINEAR_MAPPING
      stiffnessMatrix%ELEMENT_MATRIX%matrix=0.0_DP

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%SUBTYPE,"*",err,error))// &
        & " is not valid for a Stree equation type of a fluid mechanics equations set class."
      CALL FLAG_ERROR(localError,err,error,*999)
    END SELECT

    CALL EXITS("STREE_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("STREE_FINITE_ELEMENT_CALCULATE",err,error)
    CALL EXITS("STREE_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE STREE_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !      

  !>Evaluates the residual nodal stiffness matrices and RHS for a Stree equation nodal equations set.
  SUBROUTINE Stree_PRE_SOLVE(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,parentLoop,navierstokesLoop
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,navierstokesEquationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations,navierstokesEquations
    TYPE(FIELD_TYPE), POINTER :: equationsSetField,materialsField,geometricField,navierstokesDependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable,dependentVariable,dependentFieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations,navierstokesSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,navierstokesSolverMapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: navierstokesSolver
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeNumber,variableIdx,BOUNDARY_CONDITION_CHECK_VARIABLE,nodeIdx,componentIdx,derivativeIdx,versionIdx
    INTEGER(INTG) :: dependentDof,dependentVariableType,userNodeNumber
    REAL(DP) :: currentTime,timeIncrement,flow

    CALL ENTERS("Stree_PRE_SOLVE",err,error,*999)

    ! Some preliminary sanity checks
    IF(ASSOCIATED(SOLVER)) THEN
      solvers=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,currentTime,timeIncrement,ERR,ERROR,*999)
        parentLoop=>controlLoop%PARENT_LOOP
        navierstokesLoop=>parentLoop%SUB_LOOPS(2)%PTR
        navierstokesSolver=>navierstokesLoop%SOLVERS%SOLVERS(2)%PTR
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          SELECT CASE(controlLoop%PROBLEM%SUBTYPE)
          CASE(PROBLEM_Stree1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_Stree1D0DAdv_NAVIER_STOKES_SUBTYPE)
            solverEquations=>solver%SOLVER_EQUATIONS
            navierstokesSolverEquations=>navierstokesSolver%SOLVER_EQUATIONS
            IF(ASSOCIATED(solverEquations)) THEN
              solverMapping=>solverEquations%SOLVER_MAPPING
              navierstokesSolverMapping=>navierstokesSolverEquations%SOLVER_MAPPING
              IF(ASSOCIATED(solverMapping)) THEN
                equationsSet=>solverMapping%EQUATIONS_SETS(1)%PTR
                navierstokesEquationsSet=>navierstokesSolverMapping%EQUATIONS_SETS(1)%PTR
                IF(ASSOCIATED(equationsSet)) THEN
                  equations=>equationsSet%EQUATIONS
                  navierstokesEquations=>navierstokesEquationsSet%EQUATIONS
                  IF(ASSOCIATED(equations)) THEN
                    materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
                    navierstokesDependentField=>navierstokesEquationsSet%DEPENDENT%DEPENDENT_FIELD
                    IF(.NOT.ASSOCIATED(materialsField)) THEN
                      CALL FLAG_ERROR("Dependent field is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Equations set equations is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FLAG_ERROR("Solver mapping is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FLAG_ERROR("Solver equations is not associated.",err,error,*999)
            END IF
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%SUBTYPE,"*",err,error))// &
              & " is not valid for boundary flux calculation."
            CALL FLAG_ERROR(localError,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FLAG_ERROR("Solvers is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FLAG_ERROR("Solver is not associated.",err,error,*999)
    END IF

    BOUNDARY_CONDITIONS=>navierstokesSolverEquations%BOUNDARY_CONDITIONS
    DO variableIdx=1,navierstokesDependentField%NUMBER_OF_VARIABLES
      dependentVariableType=navierstokesDependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
      NULLIFY(dependentFieldVariable)
      CALL FIELD_VARIABLE_GET(navierstokesDependentField,dependentVariableType,dependentFieldVariable,ERR,ERROR,*999)
      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
        & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
        IF(ASSOCIATED(dependentFieldVariable)) THEN
          DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
            IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
              IF(ASSOCIATED(domain)) THEN
                IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                  domainNodes=>domain%TOPOLOGY%NODES
                  IF(ASSOCIATED(domainNodes)) THEN
                    !Loop over the local nodes excluding the ghosts.
                    DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                      userNodeNumber=domainNodes%NODES(nodeIdx)%USER_NUMBER
                      DO derivativeIdx=1,domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                        DO versionIdx=1,domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                          dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                            & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(dependentDof)
                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FixedStree) THEN  
                            ! Update dependent field value1
                            IF(ASSOCIATED(materialsField)) THEN
                              CALL FIELD_PARAMETER_SET_GET_NODE(navierstokesDependentField,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,flow,ERR,ERROR,*999)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(materialsField,FIELD_V_VARIABLE_TYPE, &
                                & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,int(currentTime)+1,1,flow,err,error,*999)
                            endif
                          ENDIF
                        ENDDO !versionIdx
                      ENDDO !derivativeIdx
                    ENDDO !nodeIdx
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
          ENDDO !componentIdx
        ELSE
          CALL FLAG_ERROR("Dependent field variable is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ENDDO !variableIdx

    CALL EXITS("Stree_PRE_SOLVE")
    RETURN
999 CALL ERRORS("Stree_PRE_SOLVE",err,error)
    CALL EXITS("Stree_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE Stree_PRE_SOLVE

  !
  !================================================================================================================================
  !      

END MODULE Stree_EQUATION_ROUTINES
