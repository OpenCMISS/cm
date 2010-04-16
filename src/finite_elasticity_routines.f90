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
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
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
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATHS  
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE SOLVER_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FINITE_ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE,FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SETUP,FINITE_ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & FINITE_ELASTICITY_EQUATIONS_SET_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SUBTYPE_SET,FINITE_ELASTICITY_PROBLEM_SETUP, &
    & FINITE_ELASTICITY_POST_SOLVE, FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA, &
    & FINITE_ELASTICITY_PRE_SOLVE    

CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a finite elasticity finite element equations set.
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

        ! determine step size: will this be robust enough?
        DELTA=1e-4_dp
        CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETERS,DATA,ERR,ERROR,*999)
        xnorm=sqrt(sum(DATA**2))/size(DATA)
        DELTA=(DELTA+xnorm)*DELTA
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
    !TYPE(VARYING_STRING) :: LOCAL_ERROR   

    INTEGER(INTG) :: component_idx,component_idx2,parameter_idx,gauss_idx,element_dof_idx,FIELD_VAR_TYPE

    INTEGER(INTG) :: idx
    INTEGER(INTG) :: NDOFS,mh,mhs
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NUMBER_OF_XI,HYDROSTATIC_PRESSURE_COMPONENT
    INTEGER(INTG) :: NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
    INTEGER(INTG) :: DEPENDENT_COMPONENT_INTERPOLATION_TYPE
    INTEGER(INTG) :: DEPENDENT_NUMBER_OF_GAUSS_POINTS       
    INTEGER(INTG) :: MESH_COMPONENT_1
    REAL(DP) :: DZDNU(3,3),CAUCHY_TENSOR(3,3)
    REAL(DP) :: DFDZ(64,3) !temporary until a proper alternative is found
    REAL(DP) :: GAUSS_WEIGHTS,Jznu,Jxxi
    REAL(DP) :: THICKNESS ! for elastic membrane
    
    CALL ENTERS("FINITE_ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN 
        !Grab pointers: matrices, fields, basis
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        EQUATIONS_MAPPING =>EQUATIONS%EQUATIONS_MAPPING

        FIBRE_FIELD    =>EQUATIONS%INTERPOLATION%FIBRE_FIELD
        GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
        MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
        DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD

        DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS       
        DEPENDENT_QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        DEPENDENT_NUMBER_OF_GAUSS_POINTS=DEPENDENT_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
        DEPENDENT_NUMBER_OF_COMPONENTS=DEPENDENT_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS

        NUMBER_OF_DIMENSIONS = EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
        NUMBER_OF_XI = DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI

        !Initialise tensors and matrices
        DZDNU=0.0_DP
        DO idx=1,3
          DZDNU(idx,idx)=1.0_DP
        ENDDO
        CAUCHY_TENSOR=DZDNU ! copy the identity matrix
        DFDZ=0.0_DP ! (parameter_idx,component_idx)
  
        !Grab interpolation parameters
        FIELD_VAR_TYPE=FIELD_U_VARIABLE_TYPE !Future-proofing for coupled problems (applies only to dependent field)
        GEOMETRIC_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        FIBRE_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        MATERIALS_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR
        DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_VAR_TYPE)%PTR
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE) THEN  
          DARCY_DEPENDENT_INTERPOLATION_PARAMETERS=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS(FIELD_V_VARIABLE_TYPE)%PTR
        END IF

        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & GEOMETRIC_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & FIBRE_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & MATERIALS_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
          & DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE) THEN  
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DARCY_DEPENDENT_INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        END IF

        !Point interpolation pointer
        GEOMETRIC_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        FIBRE_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        MATERIALS_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%PTR
        DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_VAR_TYPE)%PTR
        IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE) THEN  
          DARCY_DEPENDENT_INTERPOLATED_POINT=>EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT(FIELD_V_VARIABLE_TYPE)%PTR
        END IF

        !SELECT: Compressible or incompressible cases
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE, EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE, &
          & EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE, EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, &
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE) ! 4 dependent components

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
          IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE) THEN  
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
              & DARCY_DEPENDENT_INTERPOLATED_POINT,ERR,ERROR,*999) ! 'FIRST_PART_DERIV' required ???
          END IF

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR(DEPENDENT_INTERPOLATED_POINT, &
              & GEOMETRIC_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT,NUMBER_OF_DIMENSIONS, &
              & NUMBER_OF_XI,DZDNU,Jxxi,ERR,ERROR,*999)

            !Calculate Sigma=1/Jznu.FTF', the Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
              & MATERIALS_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ERR,ERROR,*999)

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
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(component_idx)% &
                  & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
                  element_dof_idx=element_dof_idx+1
                    DO component_idx2=1,NUMBER_OF_DIMENSIONS
                      NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                        & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                        & GAUSS_WEIGHTS*Jxxi*Jznu*THICKNESS*CAUCHY_TENSOR(component_idx,component_idx2)* &
                        & DFDZ(parameter_idx,component_idx2)
                    ENDDO ! component_idx2 (inner component index)
                ENDDO ! parameter_idx (residual vector loop)
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based - probably not required
                element_dof_idx=element_dof_idx+1
                !TODO:  element-based interpolation will probably never be used...
              ENDIF
            ENDDO ! component_idx

            !Hydrostatic pressure component (skip for membrane problems)
            IF (EQUATIONS_SET%SUBTYPE /= EQUATIONS_SET_MEMBRANE_SUBTYPE) THEN
              HYDROSTATIC_PRESSURE_COMPONENT=DEPENDENT_FIELD%VARIABLES(1)%NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                COMPONENT_BASIS=>DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)%DOMAIN% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                COMPONENT_QUADRATURE_SCHEME=>COMPONENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_FIELD%VARIABLES(1)% &
                COMPONENTS(HYDROSTATIC_PRESSURE_COMPONENT)% &
                  & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS 
                  element_dof_idx=element_dof_idx+1 
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                    & GAUSS_WEIGHTS*Jxxi*COMPONENT_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(parameter_idx,1,gauss_idx)*(Jznu-1.0_DP)
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based
                element_dof_idx=element_dof_idx+1
                NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                  & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+GAUSS_WEIGHTS*Jxxi*(Jznu-1.0_DP)      
              ENDIF
            ENDIF
          ENDDO !gauss_idx

        CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE)   ! compressible problem (no pressure component)

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

            !Calculate F=dZ/dNU at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR(DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT, &
              & FIBRE_INTERPOLATED_POINT,NUMBER_OF_DIMENSIONS,NUMBER_OF_XI,DZDNU,Jxxi,ERR,ERROR,*999)
            !Caculate Cauchy stress tensor at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
            & MATERIALS_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ERR,ERROR,*999)
            !Calculate dF/DZ at the gauss point
            CALL FINITE_ELASTICITY_GAUSS_DFDZ(DEPENDENT_INTERPOLATED_POINT,ELEMENT_NUMBER,gauss_idx,NUMBER_OF_DIMENSIONS, &
            & NUMBER_OF_XI,DFDZ,ERR,ERROR,*999)

            !Add up the residual terms
            element_dof_idx=0
            DO component_idx=1,DEPENDENT_NUMBER_OF_COMPONENTS
              DEPENDENT_COMPONENT_INTERPOLATION_TYPE=DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(component_idx)%INTERPOLATION_TYPE
              IF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN !node based
                NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=DEPENDENT_FIELD%VARIABLES(1)%COMPONENTS(component_idx)% &
                  & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
                DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS  
                  element_dof_idx=element_dof_idx+1    
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)= &
                    & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(element_dof_idx)+ &
                    & GAUSS_WEIGHTS*Jxxi*Jznu*(CAUCHY_TENSOR(component_idx,1)*DFDZ(parameter_idx,1)+ &
                    & CAUCHY_TENSOR(component_idx,2)*DFDZ(parameter_idx,2)+ &
                    & CAUCHY_TENSOR(component_idx,3)*DFDZ(parameter_idx,3)) 
                ENDDO
              ELSEIF(DEPENDENT_COMPONENT_INTERPOLATION_TYPE==FIELD_ELEMENT_BASED_INTERPOLATION) THEN !element based - probably not required
                element_dof_idx=element_dof_idx+1
                !TODO:  element-based interpolation will probably never be used...
              ENDIF
            ENDDO !component_idx
          ENDDO !gauss_idx
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
        FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLES(1) ! 'U' variable
        DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
          MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
          DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          NDOFS = NDOFS + DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
        END DO
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"NDOFS: ",NDOFS,ERR,ERROR,*999)

!               CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number * (Fin.Elast.):",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Element Vector for element number (Fin.Elast.): ", &
          & ELEMENT_NUMBER,ERR,ERROR,*999)
!               DO mhs=1,NDOFS
!                 CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"row number = ",mhs,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS,NDOFS,NDOFS,&
            & NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR(:), &
            & '("",4(X,E13.6))','4(4(X,E13.6))',ERR,ERROR,*999)
          CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE," ",ERR,ERROR,*999)
!               END DO
      ENDIF
    ENDiF
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

  !>Evaluates the deformation gradient tensor at a given Gauss point
  !> Jznu is not used if IS_2D_ELEMENT_IN_3D_SPACE is true
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DEFORMATION_GRADIENT_TENSOR(DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT,&
    & FIBRE_INTERPOLATED_POINT,DIMEN,NUMBER_OF_XI,DZDNU,Jxxi,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,GEOMETRIC_INTERPOLATED_POINT, &
      & FIBRE_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: DZDNU(3,3) !DZDNU - Deformation Gradient Tensor,
    INTEGER(INTG), INTENT(IN) :: DIMEN
    REAL(DP) :: DZDNU_TEMP(DIMEN,DIMEN),Jxxi
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
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
  SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR(EQUATIONS_SET,DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,CAUCHY_TENSOR,Jznu,DZDNU,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<A pointer to the equations set 
    REAL(DP), INTENT(IN) :: DZDNU(3,3)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: DEPENDENT_INTERPOLATED_POINT,MATERIALS_INTERPOLATED_POINT
    REAL(DP), INTENT(OUT) :: CAUCHY_TENSOR(:,:)
    REAL(DP), INTENT(OUT) :: Jznu !Deformation Gradient Tensor
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: EQUATIONS_SET_SUBTYPE !<The equation subtype
    INTEGER(INTG) :: i,j,PRESSURE_COMPONENT
    REAL(DP) :: AZL(3,3),AZU(3,3),DZDNUT(3,3),PIOLA_TENSOR(3,3),E(3,3),TEMP(3,3),TEMPTERM,I3,P,I1
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP), DIMENSION (:), POINTER :: C !Parameters for constitutive laws
    REAL(DP) :: a, B(3,3), Q !Parameters for orthotropic laws

    CALL ENTERS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR,*999)
    EQUATIONS_SET_SUBTYPE = EQUATIONS_SET%SUBTYPE

    C => MATERIALS_INTERPOLATED_POINT%VALUES(:,1)

    CALL MATRIX_TRANSPOSE(DZDNU,DZDNUT,ERR,ERROR,*999)
    CALL MATRIX_PRODUCT(DZDNUT,DZDNU,AZL,ERR,ERROR,*999) !AZL = F'*F (deformed covariant or right cauchy deformation tensor, C)

    PRESSURE_COMPONENT=DEPENDENT_INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD_VARIABLE%NUMBER_OF_COMPONENTS
    P=DEPENDENT_INTERPOLATED_POINT%VALUES(PRESSURE_COMPONENT,1)

    CALL INVERT(AZL,AZU,I3,ERR,ERROR,*999) !AZU - deformed contravariant tensor; I3 = det(C)
    E = 0.5_DP*AZL !Green-Lagrange strain tensor, E
    DO i=1,3
      E(i,i)=E(i,i)-0.5_DP
    ENDDO

    !PIOLA_TENSOR is the second Piola-Kirchoff tensor (PK2 or S)
    !p is the actual hydrostatic pressure, not double it

    SELECT CASE(EQUATIONS_SET_SUBTYPE)
    CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+p*(I3-1)
      !Also assumed I3 = det(AZL) = 1.0

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
        PIOLA_TENSOR(1,2)=2.0_DP*(       C(2)*(-AZL(2,1))        +P*AZU(1,2))
        PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
        PIOLA_TENSOR(2,2)=2.0_DP*(C(1)+C(2)*(AZL(3,3)+AZL(1,1))+P*AZU(2,2))

    CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1/2 (e^Q - 1)
      ! where Q=c2(I1-3)
      ! S = 2*dW/dC + 2pC^-1
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)

      TEMPTERM=C(1)*C(2)*EXP(AZL(1,1)+AZL(2,2)+AZL(3,3)-3.0_DP)
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
        CALL FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
          & PIOLA_TENSOR(1,1),ERR,ERROR,*999)
      END IF
    CASE (EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE)
      !Form of constitutive model is:
      ! W=c1*(I1-3)+c2*(I2-3)+c3*(J-1)^2   (this is actually nearly incompressible)
      C(1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
      C(2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
      C(3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
      PIOLA_TENSOR(1,1)=C(1)+C(2)*(AZL(2,2)+AZL(3,3))
      PIOLA_TENSOR(1,2)=C(2)*(-AZL(2,1))
      PIOLA_TENSOR(1,3)=C(2)*(-AZL(3,1))   
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=C(1)+C(2)*(AZL(3,3)+AZL(1,1))
      PIOLA_TENSOR(2,3)=C(2)*(-AZL(3,2))     
      PIOLA_TENSOR(3,1)=PIOLA_TENSOR(1,3)
      PIOLA_TENSOR(3,2)=PIOLA_TENSOR(2,3)
      PIOLA_TENSOR(3,3)=C(1)+C(2)*(AZL(1,1)+AZL(2,2))
      PIOLA_TENSOR=PIOLA_TENSOR+C(3)*(I3-SQRT(I3))*AZU
      PIOLA_TENSOR=PIOLA_TENSOR*2.0_DP
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
      TEMPTERM=C(1)*EXP(C(2)*(I1-3))
      PIOLA_TENSOR(1,1)=-P*AZU(1,1)+TEMPTERM+2.0*C(3)*(AZL(1,1)-1)*EXP(C(5)*(AZL(1,1)-1)**2)
      PIOLA_TENSOR(1,2)=-P*AZU(1,2)+C(7)*AZL(1,2)*EXP(C(8)*AZL(1,2)**2)
      PIOLA_TENSOR(1,3)=-P*AZU(1,3)
      PIOLA_TENSOR(2,1)=PIOLA_TENSOR(1,2)
      PIOLA_TENSOR(2,2)=-P*AZU(2,2)+TEMPTERM+2.0*C(4)*(AZL(2,2)-1)*EXP(C(6)*(AZL(2,2)-1)**2)
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
    
    Jznu=DETERMINANT(AZL,ERR,ERROR)**0.5_DP
    CAUCHY_TENSOR=CAUCHY_TENSOR/Jznu

    NULLIFY(C)

    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_TENSOR

  !
  !================================================================================================================================
  !

  ! calculates the current active contraction component using the independent field
  ! TODO: actual active contraction law. current: Ta = Kt
  ! TODO: which element/gp are we in?
  SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION(INDEPENDENT_FIELD,CAUCHY_FF,ERR,ERROR,*)
    !Argument variables
    REAL(DP), INTENT(INOUT) :: CAUCHY_FF !< the (1,1)=(fiber,fiber) component of the stress tensor
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string    

    REAL(DP), PARAMETER :: K = 2  ! kPa/ms, rate of increase of active tension.  Ta = Kt
    REAL(DP) :: TIME
  
    CALL ENTERS("FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION",ERR,ERROR,*999)
    ! TODO: which element/gp are we in? : here 1/1 since it's constant anyway
    CALL FIELD_PARAMETER_SET_GET_GAUSS_POINT(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
      & 2, TIME ,ERR,ERROR,*999)
!    WRITE (*,*) 'ADDING ',K * TIME,' EXTRA FORCE'
    CAUCHY_FF = CAUCHY_FF + K * TIME

    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION")
    RETURN 1  
  END SUBROUTINE FINITE_ELASTICITY_GAUSS_CAUCHY_ADD_ACTIVE_CONTRACTION

  !
  !================================================================================================================================
  !

  !>Evaluates df/dz (derivative of interpolation function wrt deformed coord) matrix at a given Gauss point
  SUBROUTINE FINITE_ELASTICITY_GAUSS_DFDZ(INTERPOLATED_POINT,ELEMENT_NUMBER,GAUSS_POINT_NUMBER,NUMBER_OF_DIMENSIONS, &
    & NUMBER_OF_XI,DFDZ,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT       
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER       
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT_NUMBER   
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_XI
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    REAL(DP) :: DFDZ(:,:)  
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
      NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=FIELD%VARIABLES(1)%COMPONENTS(component_idx)% &
        & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
      DO parameter_idx=1,NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS
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
      NUMBER_OF_FIELD_COMPONENT_INTERPOLATION_PARAMETERS=FIELD%VARIABLES(1)%COMPONENTS(component_idx)% &
        & MAX_NUMBER_OF_INTERPOLATION_PARAMETERS
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
      & NUMBER_OF_DIMENSIONS
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    LOGICAL :: IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD

    CALL ENTERS("FINITE_ELASTICITY_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD = EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE &
                              .AND. EQUATIONS_SET%SUBTYPE/=EQUATIONS_SET_MEMBRANE_SUBTYPE

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
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
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
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SUBTYPE)
          !-----------------------------------------------------------------------
          ! Dependent field setup for single-physics
          !-----------------------------------------------------------------------
          CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE, &
              & EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE,&
              & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE,&
              & EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE, EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_MEMBRANE_SUBTYPE)
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
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
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
          CASE(EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE) 
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
                NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1 !for Darcy: velocity components and pressure
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,ERR,ERROR,*999)

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
                DO component_idx=1,NUMBER_OF_DIMENSIONS+1
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
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

                  IF (IS_HYDROSTATIC_PRESSURE_DEPENDENT_FIELD) THEN                           
                    !Elasticity: Set the hydrostatic pressure component to element based interpolation
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_COMPONENTS,FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDIF

                  !Darcy: Set the velocity and pressure components to node based interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS+1
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
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                  & ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS+1,ERR,ERROR,*999)

                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  !Elasticity:
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx

!kmith :09.0.06.09 - Hydrostatic pressure could be node-based as well
!                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
!                    & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,4, &
!                    & FIELD_ELEMENT_BASED_INTERPOLATION,ERR,ERROR,*999)
!kmith

                  !Darcy:
                  DO component_idx=1,NUMBER_OF_DIMENSIONS+1
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
             NUMBER_OF_COMPONENTS = 17  ! this number doesn't mean anything yet
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

             DO component_idx=1,NUMBER_OF_COMPONENTS
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
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
                  & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 5,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 7,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 3,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 8,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                    & " is not valid for a finite elasticity equation type of an elasticity equation set class."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
                !Default the materials components to the geometric interpolation setup with constant interpolation
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
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
                SELECT CASE(EQUATIONS_SET%SUBTYPE)
                CASE(EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE,EQUATIONS_SET_NO_SUBTYPE, &
                  & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_MEMBRANE_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,3,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_TRANSVERSE_ISOTROPIC_EXPONENTIAL_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,5,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_COSTA_SUBTYPE,EQUATIONS_SET_ACTIVECONTRACTION_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,7,ERR,ERROR,*999)        
                CASE(EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,3,ERR,ERROR,*999)
                CASE(EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,8,ERR,ERROR,*999)
                CASE DEFAULT
                  LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                    & " is not valid for a finite elasticity equation type of an elasticity equation set class."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
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
                !Default to Mooney-Rivlin for now and set the c10 and c01 constants to be 2.0 and 3.0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,2.0_DP,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,3.0_DP,ERR,ERROR,*999)
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
              !Do nothing
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
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
              END IF
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
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
          & EQUATIONS_SET_ORTHOTROPIC_MATERIAL_HOLZAPFEL_OGDEN_SUBTYPE)        
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
          & EQUATIONS_SET_INCOMPRESSIBLE_FINITE_ELASTICITY_DARCY_SUBTYPE, &
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

    CALL ENTERS("FINITE_ELASTICITY_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==3) THEN
                CALL FINITE_ELASTICITY_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              END IF
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
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(VARYING_STRING) :: METHOD !,FILE
    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE
    LOGICAL :: EXPORT_FIELD
    INTEGER(INTG) :: CURRENT_LOOP_ITERATION
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER
    INTEGER(INTG) :: equations_set_idx

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
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now exporting fields... ",ERR,ERROR,*999)
                        CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER, &
                          & "STATICSOLIDSOLUTION",ERR,ERROR,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"STATICSOLUTION",ERR,ERROR,*999)
                        CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                      ENDIF
                    ENDDO
                  ENDIF 
                ENDIF
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR

                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
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
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now exporting fields...",ERR,ERROR,*999)
                            CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
                              & ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"All fields exported...",ERR,ERROR,*999)
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
              WRITE(*,*) "t=",CONTROL_LOOP%TIME_LOOP%CURRENT_TIME
              ! how to check eqn subtype?
              
              INDEPENDENT_FIELD => SOLVER%SOLVERS%SOLVERS(1)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                                     & EQUATIONS_SETS(1)%PTR%INDEPENDENT%INDEPENDENT_FIELD !?
              ! set component 2 to current time.
              CALL FIELD_COMPONENT_VALUES_INITIALISE(INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,2,CONTROL_LOOP%TIME_LOOP%CURRENT_TIME,ERR,ERROR,*999)

            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite elasticity pre solve ... ",ERR,ERROR,*999)

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

              IF(SOLVER%GLOBAL_NUMBER==3) THEN
                !--- 3.1 Get the Darcy pressure
                CALL FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !--- 3.2 For PGM: Get the displacement field
                IF(CONTROL_LOOP%PROBLEM%SUBTYPE==PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE) THEN
                  CALL FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                END IF
              END IF
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

  !>Copy Darcy pressure into independent field of finite elasticity
  SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_FINITE_ELASTICITY, SOLVER_DARCY  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD_FINITE_ELASTICITY, DEPENDENT_FIELD_DARCY
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_FINITE_ELASTICITY, SOLVER_EQUATIONS_DARCY  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_FINITE_ELASTICITY, SOLVER_MAPPING_DARCY !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_FINITE_ELASTICITY, EQUATIONS_SET_DARCY !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    REAL(DP), POINTER :: DUMMY_VALUES2(:)

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_INDEPENDENT_FIELD_FINITE_ELASTICITY,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DARCY
    INTEGER(INTG) :: NDOFS_TO_PRINT


    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_FINITE_ELASTICITY)
      NULLIFY(SOLVER_DARCY)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SUBTYPE,PROBLEM_QUASISTATIC_FINITE_ELASTICITY_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==3) THEN
                !--- Get the independent field of the finite elasticity equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Copying Darcy pressure into independent field of finite elasticity ... ", &
                  & ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,3,SOLVER_FINITE_ELASTICITY,ERR,ERROR,*999)
                SOLVER_EQUATIONS_FINITE_ELASTICITY=>SOLVER_FINITE_ELASTICITY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_FINITE_ELASTICITY)) THEN
                  SOLVER_MAPPING_FINITE_ELASTICITY=>SOLVER_EQUATIONS_FINITE_ELASTICITY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_FINITE_ELASTICITY)) THEN
                    EQUATIONS_SET_FINITE_ELASTICITY=>SOLVER_MAPPING_FINITE_ELASTICITY%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_FINITE_ELASTICITY)) THEN
                      INDEPENDENT_FIELD_FINITE_ELASTICITY=>EQUATIONS_SET_FINITE_ELASTICITY%INDEPENDENT%INDEPENDENT_FIELD
                      IF(ASSOCIATED(INDEPENDENT_FIELD_FINITE_ELASTICITY)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(INDEPENDENT_FIELD_FINITE_ELASTICITY, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_INDEPENDENT_FIELD_FINITE_ELASTICITY,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("INDEPENDENT_FIELD_FINITE_ELASTICITY is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Finite elasticity equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Finite elasticity solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Finite elasticity solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Get the dependent field for the Darcy equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DARCY,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DARCY=>SOLVER_DARCY%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DARCY)) THEN
                  SOLVER_MAPPING_DARCY=>SOLVER_EQUATIONS_DARCY%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DARCY)) THEN
                    EQUATIONS_SET_DARCY=>SOLVER_MAPPING_DARCY%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DARCY)) THEN
                      DEPENDENT_FIELD_DARCY=>EQUATIONS_SET_DARCY%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DARCY)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DARCY, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DARCY,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DARCY is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Darcy equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Darcy solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Darcy solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the result from Darcy's dependent field (pressure) to finite-elasticity's independent field
                CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DARCY, & 
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,INDEPENDENT_FIELD_FINITE_ELASTICITY, & 
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,ERR,ERROR,*999)

                IF(DIAGNOSTICS3) THEN
                  NULLIFY( DUMMY_VALUES2 )
                  CALL FIELD_PARAMETER_SET_DATA_GET(INDEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
                  NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
                  CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
                    & '(" INDEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
                    & '4(4(X,E13.6))',ERR,ERROR,*999)
                  CALL FIELD_PARAMETER_SET_DATA_RESTORE(INDEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
                ENDIF

              ELSE  
                ! do nothing ???
              END IF
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

    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE")
    RETURN
999 CALL ERRORS("FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE",ERR,ERROR)
    CALL EXITS("FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE")
    RETURN 1
  END SUBROUTINE FINITE_ELASTICITY_PRE_SOLVE_GET_DARCY_PRESSURE

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

    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)
    REAL(DP), POINTER :: DUMMY_VALUES2(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP) :: NSUBTRACT

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_FINITE_ELASTICITY
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,NDOFS_TO_PRINT, I
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION


    CALL ENTERS("FINITE_ELASTICITY_PRE_SOLVE_GET_SOLID_DISPLACEMENT",ERR,ERROR,*999)

!--- \todo : Do we need for each case a FIELD_PARAMETER_SET_UPDATE_START / FINISH on FIELD_MESH_DISPLACEMENT_SET_TYPE ?

    NULLIFY(SOLVER_FINITE_ELASTICITY)
    NULLIFY(MESH_DISPLACEMENT_VALUES)
    NULLIFY(DUMMY_VALUES2)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
              !--- Motion: read in from a file
              IF(SOLVER%GLOBAL_NUMBER==3) THEN
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,3,SOLVER_FINITE_ELASTICITY,ERR,ERROR,*999)
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
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Motion read in from a file ... ",ERR,ERROR,*999)

                    CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET_FINITE_ELASTICITY%GEOMETRY%GEOMETRIC_FIELD, & 
                      & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)

                    !Copy input to Finite elasticity's dependent field
                    INPUT_TYPE=42
                    INPUT_OPTION=2
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET_FINITE_ELASTICITY%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,MESH_DISPLACEMENT_VALUES,ERR,ERROR,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,MESH_DISPLACEMENT_VALUES, & 
                      & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CURRENT_TIME)
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


END MODULE FINITE_ELASTICITY_ROUTINES
