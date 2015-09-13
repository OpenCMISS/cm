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

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC Characteristic_EquationsSetSolutionMethodSet
  PUBLIC Characteristic_EquationsSet_SubtypeSet
  PUBLIC Characteristic_EquationsSet_Setup
  PUBLIC Characteristic_NodalResidualEvaluate
  PUBLIC Characteristic_NodalJacobianEvaluate
  PUBLIC Characteristic_Extrapolate
  PUBLIC Characteristic_PrimitiveToCharacteristic

CONTAINS 

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Characteristic equation type of an fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: solutionMethod
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Characteristic_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)                                
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          equationsSet%SOLUTION_METHOD=EQUATIONS_SET_NODAL_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
          & " is not valid for a Characteristic equation type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Characteristic_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Characteristic_EquationsSetSolutionMethodSet",err,error)
    EXITS("Characteristic_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Characteristic_EquationsSetSolutionMethodSet

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

    ENTERS("Characteristic_EquationsSet_SubtypeSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
        equationsSet%CLASS=EQUATIONS_SET_FLUID_MECHANICS_CLASS
        equationsSet%TYPE=EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE
        equationsSet%SUBTYPE=EQUATIONS_SET_CHARACTERISTIC_SUBTYPE
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
          & " is not valid for a Characteristic fluid type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Characteristic_EquationsSet_SubtypeSet")
    RETURN
999 ERRORSEXITS("Characteristic_EquationsSet_SubtypeSet",err,error)
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
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    INTEGER(INTG) :: componentIdx,geometricScalingType,geometricMeshComponent,geometricComponentNumber
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents
    INTEGER(INTG) :: independentFieldNumberOfVariables,independentFieldNumberOfComponents
    INTEGER(INTG) :: materialsFieldNumberOfVariables,materialsFieldNumberOfComponents1,materialsFieldNumberOfComponents2
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_EquationsSet_Setup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Characteristic_EquationsSetSolutionMethodSet(equationsSet, &
                & EQUATIONS_SET_NODAL_SOLUTION_METHOD,err,error,*999)
              equationsSet%SOLUTION_METHOD=EQUATIONS_SET_NODAL_SOLUTION_METHOD
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(equationsSetField,"Equations Set Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSetField,FIELD_GENERAL_TYPE,&
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(equationsSetField, &
                  & 1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSetField,&
                  & [FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_VARIABLE_LABEL_SET(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & "W2Initialise",ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSetField,&
                  & FIELD_U_VARIABLE_TYPE,1,ERR,ERROR,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0_DP,ERR,ERROR,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(equationsSetSetup% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a characteristic equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSetField,&
                  & geometricDecomposition,ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSetField,& 
                  & equationsSet%GEOMETRY%GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,ERR,ERROR,*999)                
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(equationsSetField, &
                  & FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSetField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,geometricScalingType, &
                  & ERR,ERROR,*999)
              ELSE
                !Do nothing
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a characteristic equation."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
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
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                 & FIELD_V_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                 & FIELD_U1_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                 & FIELD_U2_VARIABLE_TYPE,dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & 1,geometricMeshComponent,err,error,*999)
                !Default to the geometric interpolation setup for U,dUdN
                DO componentIdx=1,dependentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U1_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U2_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                END DO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                !Specify nodal solution method
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  ! (U, dUdN); 2 components (Q,A)
                  DO componentIdx=1,dependentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U2_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO
                  CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType, &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
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
                !calculate number of components (Q,A) for U and dUdN
                dependentFieldNumberOfComponents=2
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & dependentFieldNumberOfComponents,err,error,*999)
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
                  localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD, &
                    & "*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
                & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
                & " is invalid for a characteristic equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              independentFieldNumberOfVariables=1   !set number of variables to 1 (W)
              independentFieldNumberOfComponents=2 !normalDirection for wave relative to node for W1,W2
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
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & equationsSet%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & independentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                !characteristic normal direction (normalWave) is +/- 1
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
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
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT 
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(equationsSetSetup%FIELD,independentFieldNumberOfVariables,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & independentFieldNumberOfComponents,err,error,*999)
              ENDIF    
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard characteristic equations set"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a standard characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d 
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            materialsFieldNumberOfVariables=2 ! U type-7 constant / V type-3 variable
            materialsFieldNumberOfComponents1=8
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
                  DO componentIdx=1,materialsFieldNumberOfComponents1 !(MU,RHO,alpha,pressureExternal,LengthScale,TimeScale,MassScale)
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
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
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
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*", & 
                & err,error))//" for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*", & 
                & err,error))//" is invalid for characteristic equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristic equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSet%SUBTYPE)
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%MATERIALS
              IF(ASSOCIATED(equationsMaterials)) THEN              
                IF(equationsMaterials%MATERIALS_FINISHED) THEN
                  CALL EQUATIONS_CREATE_START(equationsSet,equations,err,error,*999)
                  CALL EQUATIONS_LINEARITY_TYPE_SET(equations,EQUATIONS_NONLINEAR,err,error,*999)
                  CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(equations,EQUATIONS_STATIC,err,error,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations materials is not associated.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EQUATIONS_SET_EQUATIONS_GET(equationsSet,equations,err,error,*999)
                CALL EQUATIONS_CREATE_FINISH(equations,err,error,*999)
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(equations,equationsMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(equationsMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(equationsMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
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
                  CALL EquationsMatrices_NonlinearStorageTypeSet(equationsMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(equationsMatrices, & 
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(equationsMatrices, & 
                    & [EQUATIONS_MATRIX_NODAL_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(equationsMatrices, & 
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(equationsMatrices, & 
                    & EQUATIONS_MATRIX_NODAL_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%SPARSITY_TYPE,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EQUATIONS_MATRICES_CREATE_FINISH(equationsMatrices,err,error,*999)
              CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a characteristics equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
              & " is invalid for a characteristics equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a characteristics equation set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
          & " does not equal a characteristics equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Characteristic_EquationsSet_Setup")
    RETURN
999 ERRORSEXITS("Characteristic_EquationsSet_Setup",err,error)
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
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:)
    REAL(DP) :: Q_BIF(4),A_BIF(4),A0_PARAM(4),E_PARAM(4),H0_PARAM(4),Beta(4),W(2,4),normalWave(2,4),SUM,rho
    INTEGER(INTG) :: derivativeIdx,versionIdx,versionIdx2,componentIdx,rowIdx,columnIdx,componentIdx2,numberOfVersions
    LOGICAL :: updateStiffnessMatrix,updateNonlinearResidual,boundaryNode

    ENTERS("Characteristic_NodalResidualEvaluate",err,error,*999)

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
            CALL FlagError("Domain is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Dependent Field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    SELECT CASE(equationsSet%SUBTYPE)
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
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
      DO componentIdx=1,2
        DO versionIdx=1,numberOfVersions
          CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
           & versionIdx,derivativeIdx,nodeNumber,componentIdx,normalWave(componentIdx,versionIdx),err,error,*999)  
        ENDDO
      ENDDO

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(ABS(normalWave(1,1))>0 .OR. ABS(normalWave(2,1))>0) THEN
        IF(.NOT. boundaryNode) THEN

          !Get material constants
          CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
          !Get node-based material parameters
          DO versionIdx=1,numberOfVersions
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,1,A0_PARAM(versionIdx),err,error,*999)  
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,2,E_PARAM(versionIdx),err,error,*999)                
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,3,H0_PARAM(versionIdx),err,error,*999)                
            beta(versionIdx)=(4.0_DP*SQRT(PI)*E_PARAM(versionIdx)*H0_PARAM(versionIdx))/(3.0_DP*A0_PARAM(versionIdx))     
          ENDDO

          DO versionIdx=1,numberOfVersions
            !Get current Q & A Values at the node
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,1,Q_BIF(versionIdx),err,error,*999)                
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,2,A_BIF(versionIdx),err,error,*999)                
            !Set as upwind field values
            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
             & FIELD_UPWIND_VALUES_SET_TYPE,versionIdx,1,nodeNumber,1,Q_BIF(versionIdx),ERR,ERROR,*999)            
            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
             & FIELD_UPWIND_VALUES_SET_TYPE,versionIdx,1,nodeNumber,2,A_BIF(versionIdx),ERR,ERROR,*999)            
            ! If A goes negative during nonlinear iteration, set to A0
            IF (A_BIF(versionIdx) < A0_PARAM(versionIdx)*0.001_DP) A_BIF(versionIdx) = A0_PARAM(versionIdx)*0.001_DP
          ENDDO

          !Get extrapolated W for the node
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,componentIdx,W(componentIdx,versionIdx),err,error,*999)                
            ENDDO
          ENDDO

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
                    & +normalWave(componentIdx,versionIdx)*4.0_DP*SQRT(Beta(versionIdx)/(2.0_DP*rho))* &
                    & (A_BIF(versionIdx)**0.25_DP - A0_PARAM(versionIdx)**0.25_DP)-W(componentIdx,versionIdx)
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
                    nonlinearMatrices%NodalResidual%vector(rowIdx)= &
                      & (rho/2.0_DP*((Q_BIF(1)/A_BIF(1))**2.0_DP) + Beta(1)*(SQRT(A_BIF(1)) - SQRT(A0_PARAM(1)))) - &
                      & (rho/2.0_DP*((Q_BIF(versionIdx)/A_BIF(versionIdx))**2.0_DP) + &
                      & Beta(versionIdx)*(SQRT(A_BIF(versionIdx)) - SQRT(A0_PARAM(versionIdx))))
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        ENDIF
      ENDIF !Find branch nodes

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
        & " is not valid for a characteristic equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Characteristic_NodalResidualEvaluate")
    RETURN
999 ERRORSEXITS("Characteristic_NodalResidualEvaluate",err,error)
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
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:)
    REAL(DP) :: Q_BIF(4),A_BIF(4),A0_PARAM(4),E_PARAM(4),H0_PARAM(4),Beta(4),W(2,4),normalWave(2,4),rho
    INTEGER(INTG) :: numberOfVersions,local_ny,startColumn2
    INTEGER(INTG) :: derivativeIdx,versionIdx,rowIdx,columnIdx,columnIdx2,startRow,endRow,componentIdx
    LOGICAL :: updateJacobianMatrix,boundaryNode

    ENTERS("Characteristic_NodalJacobianEvaluate",err,error,*999)

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
            CALL FlagError("Domain is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Dependent Field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    SELECT CASE(equationsSet%SUBTYPE)
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
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

          !Get material constants
          CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
          !Get node-based material parameters
          DO versionIdx=1,numberOfVersions
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,1,A0_PARAM(versionIdx),err,error,*999)  
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,2,E_PARAM(versionIdx),err,error,*999)                
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,3,H0_PARAM(versionIdx),err,error,*999)                
            beta(versionIdx)=(4.0_DP*SQRT(PI)*E_PARAM(versionIdx)*H0_PARAM(versionIdx))/(3.0_DP*A0_PARAM(versionIdx))     
          ENDDO

          !Get current Q & A Values at the node
          DO versionIdx=1,numberOfVersions
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,1,Q_BIF(versionIdx),err,error,*999)                
            CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & versionIdx,derivativeIdx,nodeNumber,2,A_BIF(versionIdx),err,error,*999)                
            ! If A goes negative during nonlinear iteration, set to A0
            IF (A_BIF(versionIdx) < A0_PARAM(versionIdx)*0.001_DP) A_BIF(versionIdx) = A0_PARAM(versionIdx)*0.001_DP
          ENDDO

          !Get extrapolated W for the node
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,componentIdx,W(componentIdx,versionIdx),err,error,*999)                
            ENDDO
          ENDDO

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
                    & +normalWave(componentIdx,versionIdx)*SQRT(Beta(versionIdx)/(2.0_DP*rho))*(A_BIF(versionIdx)**(-0.75_DP))
                ENDIF
              ENDDO
            ENDDO

            !Conservation of Mass
            rowIdx=numberOfVersions+1
            columnIdx = 0
            DO componentIdx=1,2
              DO versionIdx=1,numberOfVersions
                IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                  columnIdx=columnIdx+1
                  jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=normalWave(componentIdx,versionIdx)
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
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=rho* &
                        & (Q_BIF(1)/(A_BIF(1)**2.0_DP))
                      ! dP/dA
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx2)= &
                        & Beta(1)/(2.0_DP*SQRT(A_BIF(1))) - &
                        & (rho)*((Q_BIF(1)**2.0_DP)/(A_BIF(1)**3.0_DP))
                    ELSE IF(columnIdx2==rowIdx) THEN
                      ! dP/dQ
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=-rho* &
                        & (Q_BIF(versionIdx)/(A_BIF(versionIdx)**2.0_DP))
                      ! dP/dA
                      jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx2)= &
                        & -Beta(versionIdx)/(2.0_DP*SQRT(A_BIF(versionIdx))) + &
                        & (rho)*((Q_BIF(versionIdx)**2.0_DP)/(A_BIF(versionIdx)**3.0_DP))
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
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
        & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Characteristic_NodalJacobianEvaluate")
    RETURN
999 ERRORSEXITS("Characteristic_NodalJacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_NodalJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Extrapolate W for branch nodes and boundaries .
  SUBROUTINE Characteristic_Extrapolate(solver,currentTime,timeIncrement,ERR,ERROR,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER 
    REAL(DP), INTENT(IN) :: currentTime
    REAL(DP), INTENT(IN) :: timeIncrement
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: dependentBasis,materialsBasis
    TYPE(DOMAIN_TYPE), POINTER :: dependentDomain,materialsDomain
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER ::  dependentField,materialsField,independentField,geometricField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    REAL(DP) :: W(2,4),Q_EX(4),A_EX(4),XI(1),A0_PARAM(4),H0_PARAM(4),E_PARAM(4),Beta(4),normalWave(2,4),elementLengths(4)
    REAL(DP) :: A0_EX(4),H0_EX(4),E_EX(4),Beta_EX(4),f(4),l,friction
    REAL(DP) :: QPrevious,APrevious,rho,lambda(4)
    REAL(DP) :: elementLength,extrapolationDistance
    INTEGER(INTG) :: nodeIdx,versionIdx,derivativeIdx,elementIdx,elementNumber,versionElementNumber(4),lineNumber
    INTEGER(INTG) :: elementNodeIdx,elementNodeNumber,elementNodeVersion,numberOfVersions,componentIdx,numberOfLocalNodes
    LOGICAL :: overExtrapolated

    ENTERS("Characteristic_Extrapolate",ERR,ERROR,*999)

    NULLIFY(dependentBasis)
    NULLIFY(materialsBasis)
    NULLIFY(dependentDomain)
    NULLIFY(materialsDomain)
    NULLIFY(equationsSet)
    NULLIFY(equations)
    NULLIFY(geometricField)
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(materialsField)
    NULLIFY(solverEquations)
    NULLIFY(solverMapping)

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
              geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
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

                !Get normal wave direction
                normalWave=0.0_DP
                DO componentIdx=1,2
                  DO versionIdx=1,numberOfVersions
                    CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                     & derivativeIdx,nodeIdx,componentIdx,normalWave(componentIdx,versionIdx),err,error,*999)
                  ENDDO
                ENDDO

                !!!-- F i n d   B r a n c h   a n d   B o u n d a r y    N o d e s --!!!
                IF(ABS(normalWave(1,1)) > ZERO_TOLERANCE .OR. ABS(normalWave(2,1))> ZERO_TOLERANCE) THEN
                  !Get constant material parameters
                  CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)

                  overExtrapolated = .FALSE.
                  !!!-- G e t   E l e m e n t   L e n g t h s --!!!
                  elementLengths = 0.0_DP
                  DO elementIdx=1,dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS
                    elementNumber=dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%SURROUNDING_ELEMENTS(elementIdx)
                    ! Get the line lengths to extrapolate at equidistant points from the branch node
                    lineNumber = geometricField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)% &
                     & ELEMENT_LINES(1)
                    elementLength = geometricField%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(lineNumber)
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
                            versionElementNumber(versionIdx)=elementNumber
                            elementLengths(versionIdx) = elementLength
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO
                  ENDDO

                  !!!-- E x t r a p o l a t e   Q   a n d   A    V a l u e s --!!!
                  ! --------------------------------------------------------------
                  ! Extrapolate along the characteristic curve a distance x - lambda*dt from node location (x) to get 
                  ! values for W(t) from Q,A(t-delta(t)). Note that since the characteristic solver runs before the 
                  ! Navier-Stokes solver, 'previous' values are still in the 'current' field at this time-step as the 
                  ! time integration occurs as part of the Navier-Stokes solution.
                  DO componentIdx=1,2
                    DO versionIdx=1,numberOfVersions                         
                      IF(ABS(normalWave(componentIdx,versionIdx))> ZERO_TOLERANCE) THEN

                        ! Get materials values at node
                        CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                         & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,A0_PARAM(versionIdx),err,error,*999)            
                        CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                         & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,E_PARAM(versionIdx),err,error,*999)            
                        CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE, &
                         & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,3,H0_PARAM(versionIdx),err,error,*999)            
                        Beta(versionIdx) = (4.0_DP*SQRT(PI)*E_PARAM(versionIdx)*H0_PARAM(versionIdx))/ &
                          & (3.0_DP*A0_PARAM(versionIdx))

                        ! Get previous Q,A values at node
                        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,QPrevious,err,error,*999)            
                        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,APrevious,err,error,*999)            

                        ! Calculate wave speed
                        lambda(versionIdx) = QPrevious/APrevious + normalWave(componentIdx,versionIdx)* &
                         & (APrevious**0.25)*SQRT(Beta(versionIdx)/(2.0_DP*rho))
                        ! Check that lambda(1) > 0, lambda(2) < 0
                        IF (lambda(versionIdx)*normalWave(componentIdx,versionIdx) < 0.0_DP) THEN
                          CALL FlagError("Subcritical 1D system violated.",ERR,ERROR,*999)
                        ENDIF

                        ! Calculate extrapolation distance and xi location
                        extrapolationDistance = (timeIncrement)*lambda(versionIdx)
                        !  Convert to xi-space within the element
                        IF((normalWave(componentIdx,versionIdx)>ZERO_TOLERANCE)) THEN
                          ! Parent branch / outlet boundary
                          XI(1)=1.0_DP - extrapolationDistance/(elementLengths(versionIdx))
                        ELSE
                          ! Daughter branch / inlet boundary
                          XI(1)=0.0_DP - extrapolationDistance/(elementLengths(versionIdx))
                        ENDIF
                        IF (XI(1) > 1.0_DP .OR. XI(1) < 0.0_DP) THEN
                          CALL FLAG_WARNING("1D extrapolation location outside of element xi space. Reduce time increment", &
                           & ERR,ERROR,*999)
                          overExtrapolated = .TRUE.
                        ENDIF

                        ! Get Q,A values at extrapolated xi locations
                        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE, &
                         & versionElementNumber(versionIdx),EQUATIONS%INTERPOLATION% &
                         & DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,EQUATIONS%INTERPOLATION% &
                          & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        Q_EX(versionIdx)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                          & PTR%VALUES(1,NO_PART_DERIV)
                        A_EX(versionIdx)=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                          & PTR%VALUES(2,NO_PART_DERIV)
                        ! Get spatially varying material values at extrapolated xi locations
                        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE, &
                          & versionElementNumber(versionIdx),EQUATIONS%INTERPOLATION% &
                          & MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,EQUATIONS%INTERPOLATION% &
                          & MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        A0_EX(versionIdx)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)% &
                          & PTR%VALUES(1,NO_PART_DERIV)
                        E_EX(versionIdx)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)% &
                          & PTR%VALUES(2,NO_PART_DERIV)
                        H0_EX(versionIdx)=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)% &
                          & PTR%VALUES(3,NO_PART_DERIV)
                        Beta_EX(versionIdx) = (4.0_DP*SQRT(PI)*E_EX(versionIdx)*H0_EX(versionIdx))/ &
                          & (3.0_DP*A0_EX(versionIdx))
                        ! Calculate friction term if necessary
                        f(versionIdx) = -Q_EX(versionIdx)/(A_EX(versionIdx)**2.0_DP)
                      ENDIF
                    ENDDO
                  ENDDO

                  !Calculate W
                  W(:,:)=0.0_DP
                  DO componentIdx=1,2
                    DO versionIdx=1,numberOfVersions
                      IF(ABS(normalWave(componentIdx,versionIdx))>ZERO_TOLERANCE) THEN
                        ! W(t+delta(t)) = W_extrap(t)
                        W(componentIdx,versionIdx)= ((Q_EX(versionIdx)/A_EX(versionIdx))+ &
                          & normalWave(componentIdx,versionIdx)*4.0_DP*SQRT(Beta_EX(versionIdx)/(2.0_DP*rho))* &
                          & (A_EX(versionIdx)**(0.25_DP) - (A0_EX(versionIdx))**(0.25_DP)))

                        ! Add friction term if not neglected
                        l = (1.0_DP/(Q_EX(versionIdx)/A_EX(versionIdx) +  &
                          & normalWave(componentIdx,versionIdx)*A_EX(versionIdx)**0.25_DP*SQRT(Beta_EX(versionIdx)/(2.0_DP*rho))))
                        friction = timeIncrement*l*f(versionIdx)
!                        W(componentIdx,versionIdx)= W(componentIdx,versionIdx) + friction

                        ! Check extrapolated wave speed is coherent
                        lambda(versionIdx) = Q_EX(versionIdx)/A_EX(versionIdx) + normalWave(componentIdx,versionIdx)* &
                         & (A_EX(versionIdx)**0.25)*SQRT(Beta(versionIdx)/(2.0_DP*rho))
                        IF (lambda(versionIdx)*normalWave(componentIdx,versionIdx) < -ZERO_TOLERANCE ) THEN
                          CALL FlagError("Subcritical 1D system violated.",ERR,ERROR,*999)
                        ENDIF

                        IF (.NOT. overExtrapolated) THEN
                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_V_VARIABLE_TYPE, &
                           & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,componentIdx,W(componentIdx,versionIdx), &
                           & err,error,*999)
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF ! branch or boundary node
              ENDDO !Loop over nodes

            ELSE
              CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver equations are not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solvers is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("Characteristic_Extrapolate")
    RETURN
999 ERRORSEXITS("Characteristic_Extrapolate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Characteristic_Extrapolate

  !
  !================================================================================================================================
  !

  !>Calculate Characteristic (W) values based on dependent field values
  SUBROUTINE Characteristic_PrimitiveToCharacteristic(equationsSet,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(FIELD_TYPE), POINTER ::  dependentField,materialsField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeNumber,nodeIdx,derivativeIdx,versionIdx,componentIdx,numberOfVersions,dofNumber
    REAL(DP) :: qCurrent(4), aCurrent(4),W(2,4)
    REAL(DP) :: normalWave,A0_PARAM,E_PARAM,H0_PARAM,Beta
    LOGICAL :: boundaryNode

    ENTERS("Characteristic_PrimitiveToCharacteristic",ERR,ERROR,*999)

    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(materialsField)
    NULLIFY(fieldVariable)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%SUBTYPE)
      CASE(EQUATIONS_SET_Coupled1D0D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
        dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
        independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
        materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%SUBTYPE,"*",err,error))// &
          & " is not valid for a call to Characteristic_PrimitiveToCharacteristic"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%NODES

    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
      nodeNumber = domainNodes%NODES(nodeIdx)%local_number
      derivativeIdx = 1
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%numberOfVersions      
      boundaryNode=domainNodes%NODES(nodeNumber)%BOUNDARY_NODE
      !!!-- F i n d    B r a n c h    N o d e s --!!!
      IF(numberOfVersions > 1 .AND. .NOT. boundaryNode) THEN
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)            
            IF(ABS(normalWave)>ZERO_TOLERANCE) THEN
              !Get material parameters
              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,1,A0_PARAM,err,error,*999)  
              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,2,E_PARAM,err,error,*999)                
              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,3,H0_PARAM,err,error,*999)                
              Beta=(4.0_DP*SQRT(PI)*E_PARAM*H0_PARAM)/(3.0_DP*A0_PARAM)     

              ! Get current Q,A values
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,1,qCurrent(versionIdx),err,error,*999)         
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,2,aCurrent(versionIdx),err,error,*999)

              ! Calculate the characteristic based on current Q,A values
              W(componentIdx,versionIdx)= ((qCurrent(versionIdx)/aCurrent(versionIdx))+ &
               & normalWave*4.0_DP*SQRT(((Beta)))*(aCurrent(versionIdx)**(0.25_DP) - (A0_PARAM)**(0.25_DP)))

              !Update W values
              fieldVariable=>dependentField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%PTR
              dofNumber=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
               & NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,FIELD_V_VARIABLE_TYPE, &
               & FIELD_VALUES_SET_TYPE,dofNumber,W(componentIdx,versionIdx),ERR,ERROR,*999)
            ENDIF
          ENDDO
        ENDDO 
      ENDIF ! branch check
    ENDDO ! Loop over nodes

    EXITS("Characteristic_PrimitiveToCharacteristic")
    RETURN
999 ERRORSEXITS("Characteristic_PrimitiveToCharacteristic",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Characteristic_PrimitiveToCharacteristic

  !
  !================================================================================================================================
  !      

END MODULE CHARACTERISTIC_EQUATION_ROUTINES
