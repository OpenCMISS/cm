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
  SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet  !<A pointer to the equations set to evaluate the element Jacobian for
    INTEGER(INTG), INTENT(IN) :: elementNumber  !<The element number to calculate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err  !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: equations
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: equationsMatrices
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: nonlinearMatrices
    TYPE(EQUATIONS_MAPPING_NONLINEAR_TYPE), POINTER :: nonlinearMapping
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elementsTopology
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: parameters
    REAL(DP),POINTER :: columnData(:)  ! parameter set vector
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowVariable,columnVariable
    TYPE(ELEMENT_VECTOR_TYPE) :: elementVector
    INTEGER(INTG) :: componentIdx,localNy,version,derivativeIdx,derivative,nodeIdx,node,column,jacobianIdx
    INTEGER(INTG) :: componentInterpolationType
    REAL(DP) :: delta,origDepVar

    CALL ENTERS("EquationsSet_FiniteElementJacobianEvaluateFD",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
      IF(ASSOCIATED(equations)) THEN
        equationsMatrices=>equations%EQUATIONS_MATRICES
        nonlinearMatrices=>equationsMatrices%NONLINEAR_MATRICES
        nonlinearMapping=>equations%EQUATIONS_MAPPING%NONLINEAR_MAPPING
        ! The first residual variable is always the row variable, which is the variable the
        ! residual is calculated for
        rowVariable=>nonlinearMapping%RESIDUAL_VARIABLES(1)%PTR
        ! make a temporary copy of the unperturbed residuals
        ! can't reuse old results as they will be for another element
        CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
        elementVector=nonlinearMatrices%ELEMENT_RESIDUAL
        DO jacobianIdx=1,nonlinearMatrices%NUMBER_OF_JACOBIANS
          ! For coupled nonlinear problems there will be multiple Jacobians
          ! For this equations set, we calculate the residual for the row variable
          ! while pertubing parameters from the column variable.
          ! For non coupled problems these two variables will be the same
          columnVariable=>nonlinearMapping%RESIDUAL_VARIABLES(jacobianIdx)%PTR
          parameters=>columnVariable%PARAMETER_SETS%PARAMETER_SETS(FIELD_VALUES_SET_TYPE)%PTR%PARAMETERS  ! vector of dependent variables, basically
          ! determine step size
          CALL DISTRIBUTED_VECTOR_DATA_GET(parameters,columnData,err,error,*999)
          delta=MAX(MAXVAL(ABS(columnData))*1E-6_DP,1E-9)
          CALL DISTRIBUTED_VECTOR_DATA_RESTORE(parameters,columnData,err,error,*999)
          ! the actual finite differencing algorithm is about 4 lines but since the parameters are all
          ! distributed out, have to use proper field accessing routines..
          ! so let's just loop over component, node/el, derivative
          column=0  ! element jacobian matrix column number
          DO componentIdx=1,columnVariable%NUMBER_OF_COMPONENTS
            elementsTopology=>columnVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%ELEMENTS
            componentInterpolationType=columnVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE
            SELECT CASE (componentInterpolationType)
            CASE (FIELD_NODE_BASED_INTERPOLATION)
              basis=>elementsTopology%ELEMENTS(elementNumber)%BASIS
              DO nodeIdx=1,basis%NUMBER_OF_NODES
                node=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_NODES(nodeIdx)
                DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
                  derivative=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_DERIVATIVES(1,derivativeIdx,nodeIdx)
                  version=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_DERIVATIVES(2,derivativeIdx,nodeIdx)
                  localNy=columnVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                    & DERIVATIVES(derivative)%VERSIONS(version)
                  ! one-sided finite difference
                  CALL DISTRIBUTED_VECTOR_VALUES_GET(parameters,localNy,origDepVar,err,error,*999)
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar+delta,err,error,*999)
                  nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
                  CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
                  CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar,err,error,*999)
                  column=column+1
                  nonlinearMatrices%JACOBIANS(jacobianIdx)%PTR%ELEMENT_JACOBIAN%MATRIX(:,column)= &
                      & (nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR-elementVector%VECTOR)/delta
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            CASE (FIELD_ELEMENT_BASED_INTERPOLATION)
              localNy=columnVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(elementNumber)
              ! one-sided finite difference
              CALL DISTRIBUTED_VECTOR_VALUES_GET(parameters,localNy,origDepVar,err,error,*999)
              CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar+delta,err,error,*999)
              nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR=0.0_DP ! must remember to flush existing results, otherwise they're added
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
              CALL DISTRIBUTED_VECTOR_VALUES_SET(parameters,localNy,origDepVar,err,error,*999)
              column=column+1
              nonlinearMatrices%JACOBIANS(jacobianIdx)%PTR%ELEMENT_JACOBIAN%MATRIX(:,column)= &
                  & (nonlinearMatrices%ELEMENT_RESIDUAL%VECTOR-elementVector%VECTOR)/delta
            CASE DEFAULT
              CALL FLAG_ERROR("Unsupported type of interpolation.",err,error,*999)
            END SELECT
          ENDDO

          ! put the original residual back in
          nonlinearMatrices%ELEMENT_RESIDUAL=elementVector
        ENDDO
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("EquationsSet_FiniteElementJacobianEvaluateFD")
    RETURN
999 CALL ERRORS("EquationsSet_FiniteElementJacobianEvaluateFD",err,error)
    CALL EXITS("EquationsSet_FiniteElementJacobianEvaluateFD")
    RETURN 1
  END SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD

  !
  !================================================================================================================================
  !

END MODULE FiniteElementRoutines
