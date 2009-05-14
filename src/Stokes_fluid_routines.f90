!> \file
!> $Id: Stokes_fluid_routines.f90 372 2009-04-20
!> \author Sebastian Krittian
!> \brief This module handles all Stokes fluid routines.
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
!> The Original Code is openCMISS
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

!>This module handles all Stokes fluid routines.
MODULE STOKES_FLUID_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITION_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_MAPPING_ROUTINES
  USE EQUATIONS_MATRICES_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE


  PUBLIC STOKES_FLUID_EQUATIONS_SET_SETUP,STOKES_FLUID_PROBLEM_SETUP,STOKES_FLUID_FINITE_ELEMENT_CALCULATE,&
    & STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET, STOKES_FLUID_PROBLEM_SUBTYPE_SET


CONTAINS

  !
  !================================================================================================================================
  !
  ! SEBK 20/04/09
  !>Sets up the standard Stokes fluid.
  SUBROUTINE STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type to perform
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type to perform
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NEXT_NUMBER
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD, MATERIALS_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER::DEPENDENT_FIELD_NUMBER_OF_VARIABLES,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER::NUMBER_OF_DIMENSIONS,GEOMETRIC_COMPONENT_NUMBER,NUMBER_OF_MATERIALS_COMPONENTS
    INTEGER::MATERIAL_FIELD_NUMBER_OF_VARIABLES,MATERIAL_FIELD_NUMBER_OF_COMPONENTS,I

    CALL ENTERS("STOKES_FLUID_EQUATION_SET_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_STOKES_SUBTYPE) THEN
      SELECT CASE(SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
            SELECT CASE(ACTION_TYPE)
                CASE(EQUATIONS_SET_SETUP_START_ACTION)
                EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
                CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                !!TODO: Check valid setup
                CASE DEFAULT
                    LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a standard Stokes fluid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
        !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !sebk 30/04/09: define 1 dependent variable with 3(4) components for u,v,(w),p
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! DEPENDENT FIELD

            SELECT CASE(ACTION_TYPE)
                !Set start action
                CASE(EQUATIONS_SET_SETUP_START_ACTION)
                    ! find a user number
                    CALL FIELD_NEXT_NUMBER_FIND(EQUATIONS_SET%REGION,NEXT_NUMBER,ERR,ERROR,*999)
                    !start field creation with name 'DEPENDENT_FIELD'
                    CALL FIELD_CREATE_START(NEXT_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR, &
                    & ERROR,*999)
                    !start creation of a new field
                    CALL FIELD_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                    !define new created field to be dependent
                    CALL FIELD_DEPENDENT_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR, &
                    & ERROR,*999)
                    !look for decomposition rule already defined
                    CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,&
                    & ERR,ERROR,*999)
                    !apply decomposition rule found on new created field
                    CALL FIELD_MESH_DECOMPOSITION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION,ERR,&
                    & ERROR,*999)
                    !point new field to geometric field
                    CALL FIELD_GEOMETRIC_FIELD_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                    & ERR,ERROR,*999)
                    ! set number of variables to 2 (1 for U and one for DELUDELN)
                    DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                    CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,&
                    & ERR,ERROR,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                    !calculate number of components with one component for each dimension and one for pressure
                    DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                    CALL FIELD_NUMBER_OF_COMPONENTS_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,&
                    DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                    !set first i components to velocity field (2) and the (i+1)th component to the pressure field (3)
                    DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    IF(I<DEPENDENT_FIELD_NUMBER_OF_COMPONENTS) THEN
                          ! set velocity components
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,&
                          & I,2,ERR,ERROR,*999)
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,&
                          FIELD_DELUDELN_VARIABLE_TYPE,I,2,ERR,ERROR,*999)
                    ELSE
                          ! set pressure component
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,&
                          FIELD_U_VARIABLE_TYPE,I,3,ERR,ERROR,*999)
                          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,&
                          FIELD_DELUDELN_VARIABLE_TYPE,I,3,ERR,ERROR,*999)
                          END IF
                    END DO

                    SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                    !Specify fem solution method
                        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                            DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,&
                                & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,&
                                FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                                END DO
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,&
                            & ERR,ERROR,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,&
                            & ERR,ERROR,*999)
                    !Other solutions not defined yet
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
                            LOCAL_ERROR="The solution method of "&
                            & //TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// " is invalid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                        END SELECT
                !Specify finish action
                CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                   CALL FIELD_CREATE_FINISH(EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)

                CASE DEFAULT
                    LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                    & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                    & " is invalid for a standard Stokes fluid"
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
             END SELECT


        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !sebk 02/05/09: define 1 dependent variable with 1 (2) components for viscosity (density)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! MATERIALS FIELD

          !variable X with has Y components, here Y represents viscosity only
          MATERIAL_FIELD_NUMBER_OF_VARIABLES=1	!X
          MATERIAL_FIELD_NUMBER_OF_COMPONENTS=1	!Y
          SELECT CASE(ACTION_TYPE)

              !Specify start action
              CASE(EQUATIONS_SET_SETUP_START_ACTION)

              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                  IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                  ! find a user number
                  CALL FIELD_NEXT_NUMBER_FIND(EQUATIONS_SET%REGION,NEXT_NUMBER,ERR,ERROR,*999)
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(NEXT_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  CALL FIELD_TYPE_SET(NEXT_NUMBER,EQUATIONS_SET%REGION,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,&
                  & ERR,ERROR,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(NEXT_NUMBER,EQUATIONS_SET%REGION,MATERIAL_FIELD_NUMBER_OF_VARIABLES,ERR &
                  &,ERROR,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET(NEXT_NUMBER,EQUATIONS_SET%REGION,MATERIAL_FIELD_NUMBER_OF_COMPONENTS,ERR &
                  &,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)

                ELSE
                  CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
                END IF

              !Specify start action
              CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)

              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%REGION,EQUATIONS_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
                  !Set the default values for the materials field
                  !First set the mu values to 0.001
                  DO I=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                     CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,FIELD_VALUES_SET_TYPE,0.001_DP,ERR,ERROR,*999)
                  ENDDO !component_idx
                  ELSE
                 CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
                ENDIF

              CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for Stokes equation."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT


        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          !TO DO: INCLUDE GRAVITY AS SOURCE TYPE
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT


        CASE(EQUATIONS_SET_SETUP_FIXED_CONDITIONS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT


        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%FIXED_CONDITIONS)) THEN
              IF(EQUATIONS_SET%FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
                CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
                CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
                CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_STATIC,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set fixed conditions has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set fixed conditions is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_U_VARIABLE_TYPE/), &
                & ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%SPARSITY_TYPE)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                  & ERR,ERROR,*999)
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
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
     END SELECT


        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Stokes fluid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Stokes fluid subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP")



    RETURN
999 CALL ERRORS("STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP")



    RETURN 1
  END SUBROUTINE STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP

! END SEBK 20/04/09
  !
  !================================================================================================================================
  !

! SEBK 10/05/09

  !>Calculates the element stiffness matrices and RHS for a Stokes fluid finite element equations set.
  SUBROUTINE STOKES_FLUID_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ng,mh,mhs,mi,ms,nh,nhs,ni,ns,MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3),PGN,PGM,MU_PARAM
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("STOKES_FLUID_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_STANDARD_STOKES_SUBTYPE)

!!TODO: move these and scale factor adjustment out once generalised Stokes is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD

          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD

          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE

          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)

          !Define MU_PARAM
          MU_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(1,NO_PART_DERIV)

          !Loop over gauss points	2^DIM... bei 3D also ng=1,8
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)

            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & MATERIALS_INTERP_POINT,ERR,ERROR,*999)

            !Calculate RWG.
!!TODO: Think about symmetric problems.

            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0

             DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS1=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS


! sebk 11/05/09 input to determine number of parameters per component

              DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS

                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN


                  !Loop over element columns

                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS

! sebk 11/05/9 input to determine number of parameters per component
                MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS2=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS

                   DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS	! also das gleiche nur in die andere richtung
                      nhs=nhs+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MOMENTUM EQUATION FOR VELOCITIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                     IF (nh==mh.AND.nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN

! DEPENDENT_BASIS??????
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      ENDDO !ni

                      SUM=0.0_DP

! DEPENDENT_BASIS??????
                      DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          SUM=SUM+MU_PARAM*PGMSI(mi)*PGNSI(ni)*EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%GU(mi,ni)
                        ENDDO !ni
                      ENDDO !mi


! CHANGED TO MINUS SUM!!!!!

                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)-SUM*RWG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MOMENTUM EQUATION FOR PRESSURE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                   ELSE IF (nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS.AND.mh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN


                      SUM=0.0_DP
                      PGN=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

! DEPENDENT_BASIS??????
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)

! NOTE mh INDEX in DXI_DX
                        SUM=SUM+PGN*PGMSI(ni)*EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,mh)
                      ENDDO !ni

                     EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MASS EQUATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    ELSE IF (mh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS.AND.nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
!!!FIX UP THIS BIT !!!!!
!!!
!!! CAN THIS PART BE OPTIMISED BY USING THE TRANSPOSE???
!!!
                      SUM=0.0_DP
                      PGM=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

! DEPENDENT_BASIS??????
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)

! NOTE nh INDEX in DXI_DX
                        SUM=SUM+PGM*PGNSI(ni)*EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nh)
                      ENDDO !ni

                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET ALL OTHER COMPONENTS TO ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    ELSE

                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=0.0_DP

                    END IF

                    ENDDO !ns


                  ENDDO !nh
                ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NO RIGHT HAND SIDE FOR THIS CASE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEFINE SCALING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!FIX UP THIS BIT !!!!!
!!!
!!! TO DO FIX THIS BIT FOR SCALING!
!!!
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
              & DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                        & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF

        CASE(EQUATIONS_SET_GENERALISED_STOKES_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("STOKES_FLUID_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE STOKES_FLUID_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

! SEBK 20/04/09

  !>Sets up the Stokes fluid type of a fluid mechanics equations set class.
  SUBROUTINE STOKES_FLUID_EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Stokes fluid on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("STOKES_FLUID_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_STOKES_SUBTYPE)
        CALL STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_STOKES_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Stokes fluid type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("STOKES_FLUID_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE STOKES_FLUID_EQUATIONS_SET_SETUP


! END SEBK
  !
  !================================================================================================================================
  !

! SEBK 20/04/09

  !>Sets/changes the equation subtype for a Stokes fluid type of a fluid mechanics equations set class.
  SUBROUTINE STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_STOKES_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_FLUID_MECHANICS_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_STOKES_FLUID_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_STANDARD_STOKES_SUBTYPE
        CALL STOKES_FLUID_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE, &
          & EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_STOKES_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Stokes fluid type of a fluid mechanics equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE STOKES_FLUID_EQUATIONS_SET_SUBTYPE_SET

! END SEBK 20/04/09
  !
  !================================================================================================================================
  !

  !>Sets up the Stokes problem.
  SUBROUTINE STOKES_FLUID_PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Stokes fluid on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("STOKES_FLUID_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_STANDARD_STOKES_SUBTYPE)
        CALL STOKES_FLUID_PROBLEM_STANDARD_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_STOKES_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Stokes fluid type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("STOKES_FLUID_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE STOKES_FLUID_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a Stokes fluid type .
  SUBROUTINE STOKES_FLUID_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("STOKES_FLUID_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_STOKES_SUBTYPE)
        PROBLEM%CLASS=PROBLEM_FLUID_MECHANICS_CLASS
        PROBLEM%TYPE=PROBLEM_STOKES_FLUID_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_STOKES_SUBTYPE
        CALL STOKES_FLUID_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP_INITIAL_TYPE,PROBLEM_SETUP_START_ACTION, &
          & ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_STOKES_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Stokes fluid type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("STOKES_FLUID_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE STOKES_FLUID_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the standard Stokes fluid problem.
  SUBROUTINE STOKES_FLUID_PROBLEM_STANDARD_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type to perform
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type to perform
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(PROBLEM_EQUATIONS_ADD_TYPE), POINTER :: EQUATIONS_TO_ADD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("STOKES_FLUID_PROBLEM_STANDARD_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_STANDARD_STOKES_SUBTYPE) THEN
        SELECT CASE(SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
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
          CASE(PROBLEM_SETUP_DO_ACTION)
            EQUATIONS_TO_ADD=>PROBLEM%EQUATIONS_TO_ADD
            IF(ASSOCIATED(EQUATIONS_TO_ADD)) THEN
              EQUATIONS_SET=>EQUATIONS_TO_ADD%EQUATIONS_SET_TO_ADD
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
               !Check the equations set is from a standard Stokes fluid
                IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_FLUID_MECHANICS_CLASS.AND. &
                  & EQUATIONS_SET%TYPE==EQUATIONS_SET_STOKES_FLUID_TYPE.AND. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_STOKES_SUBTYPE) THEN
                  !Get the control loop
                  CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
                  CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,EQUATIONS_TO_ADD%CONTROL_LOOP_IDENTIFIER,CONTROL_LOOP,ERR,ERROR,*999)
                  !Get the solver equations
                  CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
                  CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
                  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
                  !Add in the equations set
                  CALL SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_TO_ADD%EQUATIONS_SET_TO_ADD, &
                    & EQUATIONS_TO_ADD%EQUATIONS_SET_ADDED_INDEX,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The equations set to add is not a standard Stokes fluid set.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set to add is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Problem equations to add is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Stokes fluid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Stokes fluid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Stokes fluid subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("STOKES_FLUID_PROBLEM_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("STOKES_FLUID_PROBLEM_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("STOKES_FLUID_PROBLEM_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE STOKES_FLUID_PROBLEM_STANDARD_SETUP

  !
  !================================================================================================================================
  !
END MODULE STOKES_FLUID_ROUTINES

