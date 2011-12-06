!> \file
!> \author Chris Bradley
!> \brief This module handles general finite element routines
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

!> This module handles all general finite element routines
MODULE FiniteElementRoutines

  USE BASE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EquationsSet_FiniteElementJacobianEvaluateFD

CONTAINS

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix entries using finite differencing for a general finite element equations set.
  SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: NONLINEAR_MAPPING
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: PARAMETERS
    REAL(DP),POINTER :: DATA(:) ! parameter_set vector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ROW_VARIABLE,COLUMN_VARIABLE
    TYPE(ELEMENT_VECTOR_TYPE) :: ELEMENT_VECTOR1
    INTEGER(INTG) :: component_idx,local_ny,version,derivative_idx,derivative,node_idx,node,column,jacobian_idx
    INTEGER(INTG) :: COMPONENT_INTERPOLATION_TYPE
    REAL(DP) :: DELTA, ORIG_DEP_VAR

    CALL ENTERS("EquationsSet_FiniteElementJacobianEvaluateFD",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
        NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
        NONLINEAR_MAPPING=>EQUATIONS%EQUATIONS_MAPPING%NONLINEAR_MAPPING
        ! The first residual variable is always the row variable, which is the variable the
        ! residual is calculated for
        ROW_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(1)%PTR
        ! make a temporary copy of the unperturbed residuals
        CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999) ! can't we reuse old results?
        ELEMENT_VECTOR1=NONLINEAR_MATRICES%ELEMENT_RESIDUAL
        DO jacobian_idx=1,NONLINEAR_MATRICES%NUMBER_OF_JACOBIANS
          ! For coupled nonlinear problems there will be multiple Jacobians
          ! For this equations set, we calculate the residual for the row variable
          ! while pertubing parameters from the column variable.
          ! For non coupled problems these two variables will be the same
          COLUMN_VARIABLE=>NONLINEAR_MAPPING%RESIDUAL_VARIABLES(jacobian_idx)%PTR
          PARAMETERS=>COLUMN_VARIABLE%PARAMETER_SETS%PARAMETER_SETS(1)%PTR%PARAMETERS  ! vector of dependent variables, basically
          ! determine step size
          CALL DISTRIBUTED_VECTOR_DATA_GET(PARAMETERS,DATA,ERR,ERROR,*999)
          DELTA=MAX(MAXVAL(ABS(DATA))*1E-6_DP,1E-9)
          CALL DISTRIBUTED_VECTOR_DATA_RESTORE(PARAMETERS,DATA,ERR,ERROR,*999)
          ! the actual finite differencing algorithm is about 4 lines but since the parameters are all
          ! distributed out, have to use proper field accessing routines..
          ! so let's just loop over component, node/el, derivative
          column=0  ! element jacobian matrix column number
          DO component_idx=1,COLUMN_VARIABLE%NUMBER_OF_COMPONENTS
            ELEMENTS_TOPOLOGY=>COLUMN_VARIABLE%COMPONENTS(component_idx)%DOMAIN%TOPOLOGY%ELEMENTS
            COMPONENT_INTERPOLATION_TYPE=COLUMN_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE
            SELECT CASE (COMPONENT_INTERPOLATION_TYPE)
            CASE (FIELD_NODE_BASED_INTERPOLATION)
              BASIS=>ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS
              DO node_idx=1,BASIS%NUMBER_OF_NODES
                node=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(node_idx)
                DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(node_idx)
                  derivative=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(1,derivative_idx,node_idx)
                  version=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_DERIVATIVES(2,derivative_idx,node_idx)
                  local_ny=COLUMN_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & DERIVATIVES(derivative)%VERSIONS(version)
                  ! one-sided finite difference
                  CALL DISTRIBUTED_VECTOR_VALUES_GET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR+DELTA,ERR,ERROR,*999)
                  NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
                  CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
                  column=column+1
                  NONLINEAR_MATRICES%JACOBIANS(jacobian_idx)%PTR%ELEMENT_JACOBIAN%MATRIX(:,column)= &
                      & (NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR-ELEMENT_VECTOR1%VECTOR)/DELTA
                ENDDO !derivative_idx
              ENDDO !node_idx
            CASE (FIELD_ELEMENT_BASED_INTERPOLATION)
              local_ny=COLUMN_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(ELEMENT_NUMBER)
              ! one-sided finite difference
              CALL DISTRIBUTED_VECTOR_VALUES_GET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR+DELTA,ERR,ERROR,*999)
              NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_VALUES_SET(PARAMETERS,local_ny,ORIG_DEP_VAR,ERR,ERROR,*999)
              column=column+1
              NONLINEAR_MATRICES%JACOBIANS(jacobian_idx)%PTR%ELEMENT_JACOBIAN%MATRIX(:,column)= &
                  & (NONLINEAR_MATRICES%ELEMENT_RESIDUAL%VECTOR-ELEMENT_VECTOR1%VECTOR)/DELTA
            CASE DEFAULT
              CALL FLAG_ERROR("Unsupported type of interpolation.",ERR,ERROR,*999)
            END SELECT
          ENDDO

          ! put the original residual back in
          NONLINEAR_MATRICES%ELEMENT_RESIDUAL=ELEMENT_VECTOR1
        ENDDO
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EquationsSet_FiniteElementJacobianEvaluateFD")
    RETURN
999 CALL ERRORS("EquationsSet_FiniteElementJacobianEvaluateFD",ERR,ERROR)
    CALL EXITS("EquationsSet_FiniteElementJacobianEvaluateFD")
    RETURN 1
  END SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD

  !
  !================================================================================================================================
  !

END MODULE FiniteElementRoutines
