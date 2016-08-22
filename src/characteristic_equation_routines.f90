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
  
  PUBLIC Characteristic_EquationsSetSpecificationSet
  
  PUBLIC Characteristic_EquationsSetSetup
  
  PUBLIC Characteristic_NodalResidualEvaluate
  
  PUBLIC Characteristic_NodalJacobianEvaluate
  
  PUBLIC Characteristic_Extrapolate

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
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a characteristic type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(3))
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
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " is not valid for a characteristic type of a fluid mechanics equations set."
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

  !>Sets the equation specification for a Characteristic type of a fluid mechanics equations set class.
  SUBROUTINE Characteristic_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    INTEGER(INTG), INTENT(IN) :: specification(:)
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Characteristic_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a characteristic type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a characteristic type of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_CHARACTERISTIC_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Characteristic_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Characteristic_EquationsSetSpecificationSet",err,error)
    EXITS("Characteristic_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Characteristic_EquationsSetSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Characteristic equations fluid setup.
  SUBROUTINE Characteristic_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

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
    INTEGER(INTG) :: compIdx,geometricScalingType,geometricMeshComponent,geometricComponentNumber
    INTEGER(INTG) :: dependentFieldNumberOfVariables,dependentFieldNumberOfComponents
    INTEGER(INTG) :: independentFieldNumberOfVariables,independentFieldNumberOfComponents
    INTEGER(INTG) :: materialsFieldNumberOfVariables,materialsFieldNumberOfComponents1,materialsFieldNumberOfComponents2
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_EquationsSetSetup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a characteristic type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSet%specification(3))
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
                  & equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(equationsSetField,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSetField,FIELD_GENERAL_TYPE,&
                  & ERR,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(equationsSetField, &
                  & 1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(equationsSetField,&
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & "W2Initialise",err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(equationsSetField,&
                  & FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(equationsSetSetup% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a characteristic equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsEquationsSetField=>equationsSet%EQUATIONS_SET_FIELD
              equationsSetField=>equationsEquationsSetField%EQUATIONS_SET_FIELD_FIELD
              IF(equationsEquationsSetField%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(equationsSetField,&
                  & geometricDecomposition,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(equationsSetField,& 
                  & equationsSet%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricComponentNumber,err,error,*999)                
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(equationsSetField, &
                  & FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSetField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,geometricScalingType, &
                  & ERR,error,*999)
              ELSE
                !Do nothing
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a characteristic equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSet%specification(3))
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
                !set dimension
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
                !set data type
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
                !number of components for U,DELUDELN=2 (Q,A)
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
                DO compIdx=1,dependentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,compIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,compIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_V_VARIABLE_TYPE,compIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U1_VARIABLE_TYPE,compIdx,geometricMeshComponent,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%DEPENDENT%DEPENDENT_FIELD, & 
                    & FIELD_U2_VARIABLE_TYPE,compIdx,geometricMeshComponent,err,error,*999)
                END DO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                !Specify nodal solution method
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  !(U, dUdN); 2 components (Q,A)
                  DO compIdx=1,dependentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U2_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
              localError="The third equations set specification of "// &
                & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
                & " is invalid for a characteristic equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              independentFieldNumberOfVariables=1  !set number of variables to 1 (W)
              independentFieldNumberOfComponents=1 !normalDirection for wave relative to element for W1,W2
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                !start field creation with name 'INDEPENDENT_FIELD'
                CALL FIELD_CREATE_START(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
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
                DO compIdx=1,independentFieldNumberOfComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, & 
                    & FIELD_U_VARIABLE_TYPE,compIdx,geometricMeshComponent,err,error,*999)
                ENDDO
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  DO compIdx=1,independentFieldNumberOfComponents
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !compIdx
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
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a standard characteristic equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d 
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
            materialsFieldNumberOfVariables=2 !U type-7 constant / V type-9 variable
            materialsFieldNumberOfComponents1=9
            materialsFieldNumberOfComponents2=8
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
                  DO compIdx=1,materialsFieldNumberOfComponents1 !(MU,RHO,alpha,Pext,LengthScale,TimeScale,MassScale,G0,Fr)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & compIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  ENDDO
                  DO compIdx=1,materialsFieldNumberOfComponents2 !(A0,E,H,kp,k1,k2,k3,b1)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(equationsMaterials%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & compIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a characteristic equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSet%specification(3))
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
            localError="The third equations set specification of "// &
              & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a characteristics equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a characteristics equation set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
          & " does not equal a characteristics equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Characteristic_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Characteristic_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_EquationsSetSetup

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
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: stiffnessMatrix
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER :: materialsField,dependentField,independentField
    INTEGER(INTG) :: derivIdx,i,compIdx,rowIdx,columnIdx,numberOfVersions
    REAL(DP) :: Q_BIF(4),A_BIF(4),A0(4),E,H,beta(4),W(4),normalWave(4),RHO,kp,k1,k2,k3,b1
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:)
    LOGICAL :: updateStiffnessMatrix,updateNonlinearResidual
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_NodalResidualEvaluate",err,error,*999)

    NULLIFY(domain)
    NULLIFY(domainNodes)
    NULLIFY(dependentField)
    NULLIFY(dependentParameters)
    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(independentField)
    NULLIFY(independentParameters)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(materialsField)
    NULLIFY(materialsParameters)
    NULLIFY(nonlinearMapping)
    NULLIFY(nonlinearMatrices)
    NULLIFY(stiffnessMatrix)
    updateStiffnessMatrix=.FALSE.
    updateNonlinearResidual=.FALSE.
    normalWave=0.0_DP
    derivIdx=1
    compIdx=1

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

    IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
      CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
      CALL FlagError("Equations set specification must have three entries for a characteristic type equations set.", &
        & err,error,*999)
    ENDIF
    SELECT CASE(equationsSet%specification(3))
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
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)%numberOfVersions

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(numberOfVersions>1) THEN
        !Get material constants
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 2,RHO,err,error,*999)

        !Loop over versions
        DO i=1,numberOfVersions
          !Get normal wave direction for versions
          CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,compIdx,normalWave(i),err,error,*999)
          !Get node-based material parameters
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,A0(i),err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,E,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,3,H,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,4,kp,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,5,k1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,6,k2,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,7,k3,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,8,b1,err,error,*999)
          beta(i)=kp*(A0(i)**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

          !Get current Q & A Values at the node
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,Q_BIF(i),err,error,*999)
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,A_BIF(i),err,error,*999)

          !Set as upwind field values
          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
           & FIELD_UPWIND_VALUES_SET_TYPE,i,derivIdx,nodeNumber,1,Q_BIF(i),err,error,*999)
          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
           & FIELD_UPWIND_VALUES_SET_TYPE,i,derivIdx,nodeNumber,2,A_BIF(i),err,error,*999)

          !If A goes negative during nonlinear iteration, set to A0
          IF(A_BIF(i)<A0(i)*0.001_DP) A_BIF(i)=A0(i)*0.001_DP

          !Get extrapolated W for the node
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,compIdx,W(i),err,error,*999)                
        ENDDO

        !!!-- S T I F F N E S S  M A T R I X  --!!!
        IF(updateStiffnessMatrix) THEN
          !Conservation of Mass
          rowIdx=numberOfVersions+1
          columnIdx=0
          DO i=1,numberOfVersions
            columnIdx=columnIdx+1
            stiffnessMatrix%NodalMatrix%matrix(rowIdx,columnIdx)=normalWave(i)
          ENDDO
        ENDIF !update stiffness matrix

        !!!-- N O N L I N E A R   V E C T O R --!!!
        IF(updateNonlinearResidual) THEN
          rowIdx=0
          !Characteristics Equations
          DO i=1,numberOfVersions
            rowIdx=rowIdx+1
            nonlinearMatrices%NodalResidual%vector(rowIdx)= &
              & Q_BIF(i)/A_BIF(i)+normalWave(i)*SQRT(4.0_DP*beta(i)/RHO/b1)* &        
              & ((A_BIF(i)/A0(i))**(b1/2.0_DP)-1.0_DP)-W(i)
          ENDDO

          !Continuity of Total Pressure
          DO i=1,numberOfVersions
            rowIdx=rowIdx+1
            IF(i==1) THEN
              !Do nothing
            ELSE
              nonlinearMatrices%NodalResidual%vector(rowIdx)= &
                & RHO*0.5_DP*(Q_BIF(1)/A_BIF(1))**2-          & !Dynamic pressure in parent vessel  
                & RHO*0.5_DP*(Q_BIF(i)/A_BIF(i))**2+          & !Dynamic pressure in branch vessel
                & beta(1)*((A_BIF(1)/A0(1))**b1-1.0_DP)-      & !Static pressure in parent vessel
                & beta(i)*((A_BIF(i)/A0(i))**b1-1.0_DP)         !Static pressure in branch vessel
            ENDIF
          ENDDO
        ENDIF !update nonlinear vector
      ENDIF !find branch nodes

    CASE DEFAULT
      localError="The third equations set specification of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
        & " is not valid for a characteristic type of a fluid mechanics equations set."
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
    TYPE(EQUATIONS_JACOBIAN_TYPE), POINTER :: jacobianMatrix
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: equationsMapping
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: linearMapping
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: linearMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER :: materialsField,dependentField,independentField
    INTEGER(INTG) :: numberOfVersions,startColumn2,startRow,endRow
    INTEGER(INTG) :: derivIdx,i,rowIdx,columnIdx,columnIdx2,compIdx
    REAL(DP) :: Q_BIF(4),A_BIF(4),A0(4),beta(4),W(4),normalWave(4),E,H,RHO,kp,k1,k2,k3,b1
    REAL(DP), POINTER :: dependentParameters(:),independentParameters(:),materialsParameters(:)
    LOGICAL :: updateJacobianMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("Characteristic_NodalJacobianEvaluate",err,error,*999)

    NULLIFY(domain)
    NULLIFY(domainNodes)
    NULLIFY(dependentField)
    NULLIFY(dependentParameters)
    NULLIFY(equations)
    NULLIFY(equationsMapping)
    NULLIFY(equationsMatrices)
    NULLIFY(independentField)
    NULLIFY(independentParameters)
    NULLIFY(jacobianMatrix)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(materialsField)
    NULLIFY(materialsParameters)
    NULLIFY(nonlinearMapping)
    NULLIFY(nonlinearMatrices)
    updateJacobianMatrix=.FALSE.
    normalWave=0.0_DP
    derivIdx=1
    compIdx=1

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

    IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
      CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    ELSE IF(SIZE(equationsSet%specification,1)/=3) THEN
      CALL FlagError("Equations set specification must have three entries for a characteristic type equations set.", &
        & err,error,*999)
    END IF
    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
      !Set General and Specific Pointers
      equationsMatrices=>equations%EQUATIONS_MATRICES
      equationsMapping=>equations%EQUATIONS_MAPPING
      independentField=>equations%EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
      linearMatrices=>equationsMatrices%LINEAR_MATRICES
      linearMapping=>equationsMapping%LINEAR_MAPPING
      materialsField=>equations%INTERPOLATION%MATERIALS_FIELD
      nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
      nonlinearMapping=>equationsMapping%NONLINEAR_MAPPING
      jacobianMatrix=>nonlinearMatrices%JACOBIANS(1)%PTR
      jacobianMatrix%NodalJacobian%matrix=0.0_DP
      IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%UPDATE_JACOBIAN
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivIdx)%numberOfVersions

      !!!-- F i n d   B r a n c h   N o d e s --!!!
      IF(numberOfVersions>1) THEN
        !Get material constants
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 2,RHO,err,error,*999)

        !Loop over versions
        DO i=1,numberOfVersions
          !Get normal wave direction for nodes
          CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,compIdx,normalWave(i),err,error,*999)
          !Get node-based material parameters
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,A0(i),err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,E,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,3,H,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,4,kp,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,5,k1,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,6,k2,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,7,k3,err,error,*999)
          CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,8,b1,err,error,*999)
          beta(i)=kp*(A0(i)**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

          !Get current Q & A Values at the node
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,1,Q_BIF(i),err,error,*999)
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,2,A_BIF(i),err,error,*999)
          !If A goes negative during nonlinear iteration, set to A0
          IF(A_BIF(i)<A0(i)*0.001_DP) A_BIF(i)=A0(i)*0.001_DP

          !Get extrapolated W for the node
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & i,derivIdx,nodeNumber,compIdx,W(i),err,error,*999)                
        ENDDO

        !!!--  J A C O B I A N   M A T R I X  --!!!
        IF(updateJacobianMatrix) THEN
          !Characteristic equations (dW/dQ)
          columnIdx=0
          rowIdx=0
          DO i=1,numberOfVersions
            columnIdx=columnIdx+1
            rowIdx=rowIdx+1
            jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=1.0_DP/A_BIF(i)
          ENDDO

          !Characteristic equations (dW/dA)
          rowIdx=0
          DO i=1,numberOfVersions
            columnIdx=columnIdx+1
            rowIdx=rowIdx+1
            jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)= &
              & -Q_BIF(i)/(A_BIF(i)**2)+normalWave(i)*SQRT(b1*beta(i)/RHO)/A0(i)* &
              & (A_BIF(i)/A0(i))**(b1/2.0_DP-1.0_DP)
          ENDDO

          !Continuity of Total Pressure (dP/dU)
          startRow=numberOfVersions+2
          endRow=numberOfVersions*2
          startColumn2=numberOfVersions+1
          DO rowIdx=startRow,endRow
            columnIdx=1
            columnIdx2=startColumn2
            DO i=1,numberOfVersions
              IF(columnIdx==1) THEN
                !dP/dQ
                jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)= &
                  & RHO*Q_BIF(1)/(A_BIF(1)**2)
                !dP/dA
                jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx2)=   &
                  & beta(1)/A0(1)*b1*(A_BIF(1)/A0(1))**(b1-1.0_DP) &
                  & -RHO*(Q_BIF(1)**2)/(A_BIF(1)**3)
              ELSE IF (columnIdx2==rowIdx) THEN
                !dP/dQ
                jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)= &
                  & -RHO*Q_BIF(i)/(A_BIF(i)**2)
                !dP/dA
                jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx2)=    &
                  & -beta(i)/A0(i)*b1*(A_BIF(i)/A0(i))**(b1-1.0_DP) &
                  & +RHO*(Q_BIF(i)**2)/(A_BIF(i)**3)
              ELSE
                jacobianMatrix%NodalJacobian%matrix(rowIdx,i)=0.0_DP
                jacobianMatrix%NodalJacobian%matrix(rowIdx,columnIdx)=0.0_DP
              ENDIF
              columnIdx=columnIdx+1
              columnIdx2=columnIdx2+1
            ENDDO
          ENDDO
        ENDIF !update jacobian matrix
      ENDIF !find branch nodes

    CASE DEFAULT
      localError="The third equations set specification of "// &
        & TRIM(NUMBER_TO_VSTRING(equationsSet%specification(3),"*",err,error))// &
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

  !>Extrapolate W for branch nodes and boundaries
  SUBROUTINE Characteristic_Extrapolate(solver,timeIncrement,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    REAL(DP), INTENT(IN) :: timeIncrement
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: dependentBasis,materialsBasis
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: dependentDomain,materialsDomain
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER ::  dependentField,materialsField,independentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    INTEGER(INTG) :: nodeIdx,i,derivIdx,elementIdx,elementNumber,versionElementNumber(4),lineNumber,dependentDof,compIdx,versionIdx
    INTEGER(INTG) :: elementNodeIdx,elementNodeNumber,elementNodeVersion,numberOfVersions,numberOfLocalNodes,boundaryConditionType
    INTEGER(INTG) :: nodeNumber
    REAL(DP) :: QPrevious,APrevious,A0,H,E,beta,kp,k1,k2,k3,b1,Q_EX,A_EX,A0_EX,H_EX,E_EX,beta_EX,RHO
    REAL(DP) :: elementLength,extrapolationDistance,elementLengths(4),lambda(4),normalWave(4),W(4),XI(1)
    LOGICAL :: overExtrapolated,boundaryNode

    ENTERS("Characteristic_Extrapolate",err,error,*999)

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
    W=0.0_DP
    normalWave=0.0_DP
    versionIdx=1
    derivIdx=1
    compIdx=1

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

              !Get the number of local nodes
              domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
                & TOPOLOGY%NODES
              IF(ASSOCIATED(domainNodes)) THEN
                numberOfLocalNodes=domainNodes%NUMBER_OF_NODES
              ELSE
                CALL FlagError("Domain nodes are not associated.",err,error,*999)
              END IF

              !Get constant material parameters
              CALL FIELD_PARAMETER_SET_GET_CONSTANT(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & 2,RHO,err,error,*999)

              !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
              DO nodeIdx=1,numberOfLocalNodes
                nodeNumber=domainNodes%NODES(nodeIdx)%local_number
                numberOfVersions=dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%DERIVATIVES(derivIdx)%numberOfVersions
                boundaryNode=dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
                  & TOPOLOGY%NODES%NODES(nodeIdx)%BOUNDARY_NODE

                !Get the boundary condition type for the dependent field primitive variables (Q,A)
                boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
                NULLIFY(fieldVariable)
                CALL FIELD_VARIABLE_GET(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                dependentDof=fieldVariable%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeNumber)% &
                  & DERIVATIVES(derivIdx)%VERSIONS(versionIdx)
                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable, &
                  & err,error,*999)
                boundaryConditionType=boundaryConditionsVariable%CONDITION_TYPES(dependentDof)     

                !!!-- F i n d   B r a n c h   a n d   B o u n d a r y    N o d e s --!!!
                IF(numberOfVersions>1 .OR. boundaryConditionType>0) THEN
                  overExtrapolated=.FALSE.

                  !!!-- G e t   E l e m e n t   L e n g t h s --!!!
                  elementLengths=0.0_DP
                  DO elementIdx=1,dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%NUMBER_OF_SURROUNDING_ELEMENTS
                    elementNumber=dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%SURROUNDING_ELEMENTS(elementIdx)
                    !Get the line lengths to extrapolate at equidistant points from the branch node
                    lineNumber=geometricField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%ELEMENT_LINES(1)
                    elementLength=geometricField%GEOMETRIC_FIELD_PARAMETERS%LENGTHS(lineNumber)
                    !Loop over the nodes on this (surrounding) element
                    dependentBasis=>dependentDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                    materialsBasis=>materialsDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
                    DO elementNodeIdx=1,dependentBasis%NUMBER_OF_NODES
                      elementNodeNumber=dependentDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%ELEMENT_NODES(elementNodeIdx)
                      !Check that this node is the same as the current iterative node
                      IF(elementNodeNumber==nodeIdx) THEN
                        !Loop over the versions to find the element index that matches the version
                        DO i=1,numberOfVersions
                          !Version number for the local element node
                          elementNodeVersion=dependentDomain%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%&
                            & elementVersions(1,elementNodeIdx)
                          IF(elementNodeVersion==i) THEN
                            versionElementNumber(i)=elementNumber
                            elementLengths(i)=elementLength
                          ENDIF
                        ENDDO
                      ENDIF
                    ENDDO
                  ENDDO

                  !!!-- E x t r a p o l a t e   Q   a n d   A    V a l u e s --!!!
                  ! --------------------------------------------------------------------------------------------------
                  ! Extrapolate along the characteristic curve a distance x - lambda*dt from node location (x) to get 
                  ! values for W(t) from Q,A(t-delta(t)). Note that since the characteristic solver runs before the 
                  ! Navier-Stokes solver, 'previous' values are still in the 'current' field at this time-step as the 
                  ! time integration occurs as part of the Navier-Stokes solution.
                  ! --------------------------------------------------------------------------------------------------
                  DO i=1,numberOfVersions

                    !Get normal wave direction
                    CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,i, &
                      & derivIdx,nodeIdx,compIdx,normalWave(i),err,error,*999)
                    !Get node-based material parameters
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,1,A0,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,2,E,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,3,H,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,4,kp,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,5,k1,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,6,k2,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,7,k3,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & i,derivIdx,nodeIdx,8,b1,err,error,*999)
                    beta=kp*(A0**k1)*(E**k2)*(H**k3) !(kg/m/s2) --> (Pa)

                    !Get previous Q,A values at node
                    CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,i,derivIdx,nodeIdx,1,QPrevious,err,error,*999)
                    CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,i,derivIdx,nodeIdx,2,APrevious,err,error,*999)

                    !Check wave speed is coherent
                    lambda(i)=QPrevious/APrevious+normalWave(i)*SQRT(beta/RHO*b1)*(APrevious/A0)**(b1/2.0_DP)
                    IF(lambda(i)*normalWave(i)<ZERO_TOLERANCE) THEN
                      CALL FlagError("Subcritical 1D system violated.",err,error,*999)
                    ENDIF

                    !Calculate extrapolation distance and xi location
                    extrapolationDistance=timeIncrement*lambda(i)
                    !Convert to xi-space within the element
                    IF(normalWave(i)>ZERO_TOLERANCE) THEN
                      !Parent branch
                      XI(1)=1.0_DP-extrapolationDistance/elementLengths(i)
                    ELSE
                      !Daughter branch
                      XI(1)=0.0_DP-extrapolationDistance/elementLengths(i)
                    ENDIF
                    IF(XI(1)>1.0_DP .OR. XI(1)<0.0_DP) THEN
                      CALL FLAG_WARNING("1D extrapolation location outside of element xi space. Reduce time increment", &
                        & err,error,*999)
                      overExtrapolated=.TRUE.
                    ENDIF

                    !Get Q,A values at extrapolated xi locations
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,versionElementNumber(i), &
                      & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                    CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,EQUATIONS%INTERPOLATION% &
                      & DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                    Q_EX=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                      & PTR%VALUES(1,NO_PART_DERIV)
                    A_EX=EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_U_VARIABLE_TYPE)% &
                      & PTR%VALUES(2,NO_PART_DERIV)

                    !Get spatially varying material values at extrapolated xi locations
                    CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE, &
                      & versionElementNumber(i),EQUATIONS%INTERPOLATION% &
                      & MATERIALS_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR,err,error,*999)
                    CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI,EQUATIONS%INTERPOLATION% &
                      & MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR,err,error,*999)
                    A0_EX=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV)
                    E_EX=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(2,NO_PART_DERIV)
                    H_EX=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR%VALUES(3,NO_PART_DERIV)
                    beta_EX=kp*(A0_EX**k1)*(E_EX**k2)*(H_EX**k3) !(kg/m/s2) --> (Pa)

                    !Calculate W --> W(t+delta(t))=W_extrap(t)
                    W(i)=Q_EX/A_EX+normalWave(i)*SQRT(4.0_DP*beta_EX/RHO/b1)*((A_EX/A0_EX)**(b1/2.0_DP)-1.0_DP)

                    !Check extrapolated wave speed is coherent
                    lambda(i)=Q_EX/A_EX+normalWave(i)*SQRT(beta_EX/RHO*b1)*(A_EX/A0_EX)**(b1/2.0_DP)
                    IF(lambda(i)*normalWave(i)<-ZERO_TOLERANCE) THEN
                      CALL FLAG_ERROR("Subcritical 1D system violated.",err,error,*999)
                    ENDIF

                    IF(.NOT. overExtrapolated) THEN
                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                        & i,derivIdx,nodeIdx,compIdx,W(i),err,error,*999)
                    ENDIF
                  ENDDO
                ENDIF !Branch or boundary node
              ENDDO !Loop over nodes
            ELSE
              CALL FlagError("Equations are not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver equations are not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver mapping is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("Characteristic_Extrapolate")
    RETURN
999 ERRORSEXITS("Characteristic_Extrapolate",err,error)
    RETURN 1
    
  END SUBROUTINE Characteristic_Extrapolate

  !
  !================================================================================================================================
  !      

END MODULE CHARACTERISTIC_EQUATION_ROUTINES
