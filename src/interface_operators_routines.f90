!> \file
!> \author Chris Bradley
!> \brief This module contains all interface conditions operators routines.
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
!> Contributor(s): Xiani (Nancy) Yan, Thiranja Prasad Babarenda Gamage
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

!>This module contains all interface conditions routines. 
MODULE INTERFACE_OPERATORS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CONSTANTS
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_EQUATIONS_ROUTINES
  USE INTERFACE_MAPPING_ROUTINES
  USE INTERFACE_MATRICES_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE STRINGS
  USE TIMER
  USE TYPES

#include "macros.h"  

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FieldContinuity_FiniteElementCalculate
  
  PUBLIC FrictionlessContact_FiniteElementCalculate

  PUBLIC SolidFluidOperator_FiniteElementCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for field continuity operator 
  SUBROUTINE FieldContinuity_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: GaussPoint, rowComponentIdx, rowIdx, rowParameterIdx, colComponentIdx, colIdx, colParameterIdx
    INTEGER(INTG) :: rowMeshComponentNumber,derivativeIdx,derivative,localElementNode,interfaceNode,interfaceDerivative
    INTEGER(INTG) :: coupledMeshElementNumber,coupledMeshIdx,coupledMeshVariableType,lagrangeVariableType
    INTEGER(INTG) :: connectedLine,decompositionLineNumber,localLineNodeIdx,connectedFace,decompositionFaceNumber,localFaceNodeIdx
    REAL(DP) :: XI(3),rwg,PGMSI,PGNSI,matrixCoefficient
    TYPE(BASIS_TYPE), POINTER :: interfaceDependentBasis,coupledMeshBasis,interfaceGeometricBasis, &
      & interfacePenaltyBasis,interfaceConnectivityBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: interfaceQuadratureScheme
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField,interfaceDependentField,interfaceGeometricField, &
      & interfacePenaltyField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: interfaceMatrixVariable,lagrangeVariable
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), POINTER :: interfaceInterpolation
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity
    TYPE(DOMAIN_LINE_TYPE), POINTER :: coupledMeshDomainLine
    TYPE(DOMAIN_FACE_TYPE), POINTER :: coupledMeshDomainFace
    TYPE(VARYING_STRING) :: localError

    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(BASIS_TYPE), POINTER :: coupledMeshDependentBasis
    TYPE(FIELD_TYPE), POINTER :: coupledMeshGeometricField
    INTEGER(INTG) :: meshComponentNumber,numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements
    INTEGER(INTG) :: dataPointIdx,localElementNumber,matrixElementIdx
    INTEGER(INTG) :: matrixCoefficients(2),interfaceelementnumber

    ENTERS("FieldContinuity_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FLAG_error("Interface condition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE_EQUATIONS)) CALL FLAG_error("Interface equations is not associated." &
      & ,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE)) CALL FLAG_error("Interface is not associated.",err,error,*999)

    interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS

    SELECT CASE(interfaceCondition%METHOD)

    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FLAG_error("Not implemented.",err,error,*999)

    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)

      SELECT CASE(interfaceCondition%integrationType)

      CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
        !Pointers to interface variables (columns of interface element matrix)
        interfaceInterpolation=>interfaceEquations%INTERPOLATION%INTERFACE_INTERPOLATION
        interfaceGeometricField=>interfaceInterpolation%GEOMETRIC_FIELD
        interfaceDependentField=>interfaceInterpolation%DEPENDENT_FIELD
        interfaceGeometricBasis=>interfaceGeometricField%DECOMPOSITION%DOMAIN(interfaceGeometricField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        interfaceDependentBasis=>interfaceDependentField%DECOMPOSITION%DOMAIN(interfaceDependentField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        SELECT CASE(interfaceCondition%METHOD)
        CASE(INTERFACE_CONDITION_PENALTY_METHOD)
          interfacePenaltyField=>interfaceInterpolation%PENALTY_FIELD
          interfacePenaltyBasis=>interfacePenaltyField%DECOMPOSITION%DOMAIN(interfacePenaltyField% &
            & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,interfaceInterpolation% &
            & PENALTY_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        ENDSELECT
        !Integrate using the interface quadrature scheme
        interfaceQuadratureScheme=>interfaceGeometricBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        lagrangeVariable=>interfaceEquations%INTERFACE_MAPPING%LAGRANGE_VARIABLE
        lagrangeVariableType=lagrangeVariable%VARIABLE_TYPE
        !Get element interpolation parameters from the first geometric interpolation set (to get Jacobian for interface surface integral)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,interfaceInterpolation% &
          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
        matrixCoefficient=1.0_DP
        DO coupledMeshIdx=1,interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
            IF(coupledMeshIdx>1) THEN
              matrixCoefficient=-1.0_DP
            ENDIF 
            !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
            coupledMeshDependentField=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_FIELD
            elementConnectivity=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(elementNumber,coupledMeshIdx)
            coupledMeshElementNumber=elementConnectivity%COUPLED_MESH_ELEMENT_NUMBER
            interfaceMatrixVariable=> &
              & interfaceEquations%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(coupledMeshIdx)%VARIABLE
            coupledMeshVariableType=interfaceMatrixVariable%VARIABLE_TYPE
            interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
            interfaceConnectivityBasis=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%BASIS

            !coupledMeshDependentInterpolation=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
            !  & DEPENDENT_INTERPOLATION

            !Loop over gauss points
            DO GaussPoint=1,interfaceQuadratureScheme%NUMBER_OF_GAUSS
              !CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint, &
              !  & coupledMeshDependentInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR, &
              !  & err,error,*999)
              !CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
              !  & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)

              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              rwg=interfaceInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                & JACOBIAN*interfaceQuadratureScheme%GAUSS_WEIGHTS(GaussPoint)
              IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                  & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint,interfaceInterpolation% &
                  & PENALTY_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  !Loop over the Lagrange variable matrix rows
                  DO rowParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(rowParameterIdx,NO_PART_DERIV,GaussPoint)
                    rowIdx=rowIdx+1
                    interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,rowIdx)- &
                      & (1.0_DP/interfaceInterpolation%PENALTY_INTERPOLATION(1)% &
                      & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1))*PGNSI**2.0_DP*rwg
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ELSE
                !\todo defaults to first mesh component, generalise
                !XI=InterfaceOperators_InterfToCoupledMeshGaussTransform( &
                !  & elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
                XI(1:interfaceDependentBasis%NUMBER_OF_XI)=InterfaceOperators_InterfToCoupledMeshGaussTransform( &
                  & elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
                !XI=interfaceCondition%interface%pointsConnectivity%pointsConnectivity(GaussPoint,coupledMeshIdx)%xi
                ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS

                  SELECT CASE(interfaceDependentBasis%NUMBER_OF_XI)

                  CASE(1) !1D interface (line)
                    connectedLine = elementConnectivity%CONNECTED_LINE
                    decompositionLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                    coupledMeshDomainLine=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                      & LINES%LINES(decompositionLineNumber)
                    DO localLineNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                      localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                      DO derivativeIdx=1,coupledMeshDomainLine%BASIS%NUMBER_OF_DERIVATIVES(localLineNodeIdx)
                        derivative=coupledMeshBasis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                        derivative=coupledMeshDomainLine%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                        rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV,XI,err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                          DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                            PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                            colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                              & PGNSI*PGMSI*rwg*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !localLineNodeIdx

                  CASE(2) !2D interface (face)

                    SELECT CASE(coupledMeshBasis%NUMBER_OF_XI)

                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNode=1,coupledMeshBasis%NUMBER_OF_NODES
                        DO derivative=1,coupledMeshBasis%NUMBER_OF_DERIVATIVES(localElementNode)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivative
                      ENDDO !localElementNode

                    CASE(3) !Coupled Mesh has 3 xi directions
                      connectedFace = elementConnectivity%CONNECTED_FACE
                      decompositionFaceNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                      coupledMeshDomainFace=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                        & FACES%FACES(decompositionFaceNumber)
                      DO localFaceNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                        localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                        DO derivativeIdx=1,coupledMeshDomainFace%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                          derivative=coupledMeshBasis% &
                            & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivativeIdx
                      ENDDO !FaceNodeIdx

                    END SELECT !coupledMeshBasis%NUMBER_OF_XI

                  END SELECT !interfaceDependentBasis%NUMBER_OF_XI

                ENDDO !rowComponentIdx
              ENDIF
            ENDDO !GaussPoint

            !Scale factor adjustment
            !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
            !\todo update looping of variables/components for non-zero matrix elements as done above 
            IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
              & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber, &
                  & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  !Loop over element Lagrange variable rows
                  DO rowParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    interfaceElementMatrix%MATRIX(rowIdx,rowIdx)=interfaceElementMatrix%MATRIX(rowIdx,rowIdx) * &
                      & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)% &
                      & INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)**2
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ELSE
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber, &
                  & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx) * &
                        & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)% &
                        & INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR%SCALE_FACTORS(colParameterIdx,colComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
              !Scale factor adjustment for the row dependent variable
              IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(coupledMeshElementNumber, &
                  & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                  & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,interfaceMatrixVariable%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)* &
                        & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                        & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR% &
                        & SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! coupledMeshIdx

      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)

        matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
        matrixCoefficients(2)=-1;
        interfaceElementNumber = elementNumber!todo simplify
        interface=>interfaceCondition%INTERFACE
        pointsConnectivity=>interface%pointsConnectivity
        numberOfInterfaceMeshXi=pointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
        IF(ASSOCIATED(pointsConnectivity)) THEN
          decompositionElementData=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD%DECOMPOSITION%TOPOLOGY%dataPoints% &
            & elementDataPoint(interfaceElementNumber)
         
          !Calculate PGSMI, update interface matrices with PGSMI, and update scale factors
          DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
            IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
              numberOfMatrixCoupledElements=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                & numberOfCoupledElements
              numberOfCoupledMeshXi=interface%COUPLED_MESHES(coupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
              coupledMeshGeometricField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                & GEOMETRY%GEOMETRIC_FIELD
              coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                & DEPENDENT%DEPENDENT_FIELD

              numberOfCoupledMeshGeoComp=coupledMeshGeometricField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
              interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
              !mesh component number is the same for all geometric components in elasticity problems
              meshComponentNumber=coupledMeshDependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)% &
                & MESH_COMPONENT_NUMBER
              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                  & coupledMeshElementNumber

                !Calculate the element index (non-conforming element) for this interface matrix
                matrixElementIdx=1
                DO WHILE ((localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                    & elementNumbers(matrixElementIdx)).AND.(matrixElementIdx/= &
                    & pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)%numberOfCoupledElements))
                  matrixElementIdx=matrixElementIdx+1
                ENDDO   
                xi(1:numberOfCoupledMeshXi)=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                  & xi(1:numberOfCoupledMeshXi)
                !Calculate PGSMI for each data point component
                coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                  DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    PGMSI=BASIS_EVALUATE_XI(coupledMeshDependentBasis,rowParameterIdx,NO_PART_DERIV, &
                      & xi(1:numberOfCoupledMeshXi),ERR,ERROR)*matrixCoefficients(coupledMeshIdx)
                    rowIdx=rowParameterIdx+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                    colIdx=dataPointIdx+decompositionElementData%numberOfProjectedData*(rowComponentIdx-1)
                    interfaceElementMatrix%MATRIX(rowIdx,colIdx)=PGMSI !Update interface element matrix with contact point contribution
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDDO !dataPointIdx

              !scale factor update
              IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                  localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                    & coupledMeshElementNumber
                  CALL Field_InterpolationParametersScaleFactorsElementGet(localElementNumber,interfaceEquations% &
                    & INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                    & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                  !Calculate the element index (non-conforming element) for this interface matrix
                  matrixElementIdx=1
                  DO WHILE ((localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                      & elementNumbers(matrixElementIdx)).AND.(matrixElementIdx/= &
                      & pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)%numberOfCoupledElements))
                    matrixElementIdx=matrixElementIdx+1
                  ENDDO
                  coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                  DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                    DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                      rowIdx=rowParameterIdx+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                      colIdx=dataPointIdx+decompositionElementData%numberOfProjectedData*(rowComponentIdx-1)
                      interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)* &
                        & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                        & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                        & PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                    ENDDO !rowParameterIdx
                  ENDDO !rowComponentIdx
                ENDDO !dataPointIdx
              ENDIF !.NOT. FIELD_NO_SCALING

            ENDIF !UPDATE_MATRIX
          ENDDO !coupledMeshIdx
        ELSE
          CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
        ENDIF

      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType

    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FLAG_error("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
        & " is not valid."
      CALL FLAG_error(localError,err,error,*999)
    END SELECT

    EXITS("FieldContinuity_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FieldContinuity_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FieldContinuity_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for frictionless contact operator
  SUBROUTINE FrictionlessContact_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The interface element number to calcualte the interface element matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    TYPE(INTERFACE_TYPE), POINTER :: interface !<A pointer to the interface 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField,penaltyField
    TYPE(INTERFACE_PENALTY_TYPE), POINTER :: interfacePenalty
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointsMetrics(:)
    TYPE(BASIS_TYPE), POINTER :: coupledMeshDependentBasis
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: penaltyMatrix
    INTEGER(INTG) :: meshComponentNumber,numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements,localDof
    INTEGER(INTG) :: dataPointIdx,coupledMeshIdx,xiIdx,localElementNumber,localFaceLineNumber,matrixElementIdx,rowComponentIdx, &
      & rowParameterIdx,rowIdx,colIdx,componentIdx,globalDataPointNumber
    INTEGER(INTG) :: matrixCoefficients(2)
    REAL(DP) :: PGMSI,contactStiffness
    REAL(DP) :: positionPoint(3),normalPoint(3),tangentsPoint(3,3),xi(3)
    LOGICAL :: reverseNormal
    REAL(DP), ALLOCATABLE :: gaps(:),gapsComponents(:,:),normals(:,:)
    LOGICAL, ALLOCATABLE :: orthogonallyProjected(:)
    
    
    TYPE(VARYING_STRING) :: localError

    ENTERS("FrictionlessContact_FiniteElementCalculate",err,error,*999)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
      IF(ASSOCIATED(interfaceEquations)) THEN
        interface=>interfaceCondition%INTERFACE
        IF(ASSOCIATED(interface)) THEN
          SELECT CASE(interfaceCondition%METHOD)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            SELECT CASE(interfaceCondition%integrationType)
            CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
              CALL FlagError("Mesh connectivity is not implemented for frictionless contact.",err,error,*999)
            CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
              matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
              matrixCoefficients(2)=-1;
              pointsConnectivity=>interface%pointsConnectivity
              numberOfInterfaceMeshXi=pointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
              IF(ASSOCIATED(pointsConnectivity)) THEN
                decompositionElementData=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD%DECOMPOSITION%TOPOLOGY%dataPoints% &
                  & elementDataPoint(interfaceElementNumber)
                !###################################################################################################################
                
                !Test if datapoints were orthogonally projected.  
                !\todo: Allow the user to choose to only include orthogonally projected points or not (check is commented when populating element matrix below).  
                ALLOCATE(orthogonallyProjected(decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate orthogonal projected logicals.",err,error,*999)
                orthogonallyProjected=.TRUE. !Initialise orthogonal projected logicals
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                    & DEPENDENT%DEPENDENT_FIELD
                  !mesh component number is the same for all geometric components in elasticity problems
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                    DO xiIdx=1,SIZE(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%reducedXi,1)
                      IF(ABS(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%reducedXi(xiIdx)) &
                          & < ZERO_TOLERANCE) THEN
                        orthogonallyProjected(dataPointIdx)=.FALSE.
                      ENDIF
                    ENDDO !xiIdx
                  ENDDO !dataPointIdx
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Allocate memory for local allocatable variables
                ALLOCATE(gaps(decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate gaps.",err,error,*999)
                gaps=0.0_DP !Initialise gap functions
                ALLOCATE(gapsComponents(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate component gaps.",err,error,*999)
                gapsComponents=0.0_DP !Initialise gap functions
                ALLOCATE(normals(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FlagError("Could not allocate normals.",err,error,*999)
                normals=0.0_DP !Initialise gap functions

                !Calculate Gap for each data point 
                !\todo: This is only required if only penetration is to penalized (ie seperation of meshes allowed.)
                ! If a no seperation condition is also required then calculation of the gap is not required.
                ! Need to allow user to choose which type of problem to solve.
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDependentField=>interfaceCondition%DEPENDENT%FIELD_VARIABLES(coupledMeshIdx)%PTR%FIELD
                  numberOfCoupledMeshGeoComp=coupledMeshDependentField%GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)% &
                    & PTR%NUMBER_OF_COMPONENTS
                  NULLIFY(interpolatedPoints)
                  NULLIFY(interpolationParameters)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(coupledMeshDependentField,interpolationParameters,err,error, &
                    & *999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
                    & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                  interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                    !Only interpolate if orthogonally projected
                    !\todo: Allow the user to choose to only include orthogonally projected points or not (currenlty commented out).  
                    !IF(orthogonallyProjected(dataPointIdx)) THEN
                      localElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                        & coupledMeshElementNumber
                      localFaceLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)% &
                        & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                        & elementLineFaceNumber)
                      SELECT CASE(numberOfInterfaceMeshXi) !Use face/line interpolation parameters for normal calculation
                      CASE(1)
                        CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      CASE(2)
                        SELECT CASE(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                            & elementLineFaceNumber)
                        CASE(1,3,5)
                          reverseNormal=.FALSE.
                        CASE(2,4,6)
                          reverseNormal=.TRUE.
                        END SELECT
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                      END SELECT
                      ! Determine the gap. 
                      ! \todo: Note that FIELD_INTERPOLATE_XI(FIRST_PART_DERIV by default calculates NO_PART_DERIV too
                      ! and is used because the FIRST_PART_DERIV is also need for the surface normal calculation. However the
                      ! normal is only calculated for one of the coupled bodies so unnecessary computation. Need to generalize
                      ! FIELD_INTERPOLATE_XI to allow the user to specify which PART_DERIV to calculate.  
                      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,pointsConnectivity%pointsConnectivity(globalDataPointNumber, &
                        & coupledMeshIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                      gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx)=gapsComponents(1:numberOfCoupledMeshGeoComp, &
                        & dataPointIdx)+interpolatedPoint%VALUES(1:numberOfCoupledMeshGeoComp,NO_PART_DERIV)* &
                        & matrixCoefficients(coupledMeshIdx) !Calculate 3 components gap function for each contact point
                      !Calculate surface normal (use 2nd coupled mesh surface normal)
                      !\todo: Allow the user to choose which surface normal to calculate or alternatively allow for a weighted average of the two.  
                      IF (coupledMeshIdx==2) THEN
                        CALL Field_InterpolatedPointsMetricsInitialise(interpolatedPoints,interpolatedPointsMetrics, &
                          & err,error,*999)
                        CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(numberOfCoupledMeshGeoComp,interpolatedPointsMetrics &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                        CALL Field_PositionNormalTangentsCalculateIntPtMetric(interpolatedPointsMetrics &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,reverseNormal,positionPoint,normalPoint,tangentsPoint,err,error,*999)
                        normals(1:numberOfCoupledMeshGeoComp,dataPointIdx)=normalPoint(1:numberOfCoupledMeshGeoComp)
                        CALL Field_InterpolatedPointsMetricsFinalise(interpolatedPointsMetrics,err,error,*999)
                      ENDIF !coupledMeshIdx==1
                    !ENDIF !orthogonallyProjected(dataPointIdx)
                  ENDDO !dataPointIdx
                  CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoints,err,error,*999)
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Calcualte 1 component gap
                DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                  gaps(dataPointIdx)=DOT_PRODUCT(gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx), &
                    & normals(1:numberOfCoupledMeshGeoComp,dataPointIdx))
                ENDDO !dataPointIdx
                
                !###################################################################################################################
                
                !Calculate PGSMI, update interface matrices with PGSMI, and update scale factors
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
                    numberOfMatrixCoupledElements=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                      & numberOfCoupledElements
                    numberOfCoupledMeshXi=interface%COUPLED_MESHES(coupledMeshIdx)%PTR%NUMBER_OF_DIMENSIONS
                    numberOfCoupledMeshGeoComp=interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR%GEOMETRY% &
                      & GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_COMPONENTS
                    coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                      & DEPENDENT%DEPENDENT_FIELD
                    interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
                    DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                      globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                      !\todo: Allow the user to choose gap tolerance or default to zero tolerance (currently commented out).  
                      !IF(gaps(dataPointIdx)>1.0E-10) THEN !Only add contact point contribution if the gap is a penetration
                        localElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                          & coupledMeshElementNumber
                        !Calculate the element index (non-conforming element) for this interface matrix
                        matrixElementIdx=1
                        DO WHILE (localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                            & elementNumbers(matrixElementIdx))
                          matrixElementIdx=matrixElementIdx+1
                        ENDDO
                        CALL Field_InterpolationParametersScaleFactorsElementGet(localElementNumber,interfaceEquations% &
                          & INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                          & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                        xi(1:numberOfCoupledMeshXi)=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                          & xi(1:numberOfCoupledMeshXi)                  
                        DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                          meshComponentNumber=coupledMeshDependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                            & COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                          !Calculate PGSMI for each data point component
                          coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(meshComponentNumber)%PTR% &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                          !\todo: Loop over the number of coupled mesh dependent basis element parameters on the contact face to save a bit of computation.
                          DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                            PGMSI=BASIS_EVALUATE_XI(coupledMeshDependentBasis,rowParameterIdx,NO_PART_DERIV, &
                              & xi(1:numberOfCoupledMeshXi),ERR,ERROR)*normals(rowComponentIdx,dataPointIdx)* &
                              & matrixCoefficients(coupledMeshIdx)
                            rowIdx=coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*numberOfMatrixCoupledElements* &
                              & (rowComponentIdx-1)+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS* &
                              & (matrixElementIdx-1)+rowParameterIdx
                            colIdx=dataPointIdx
                            !Update interface element matrix with contact point contribution
                            !\todo: Seperate multiplication of scale factors if required.  
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=PGMSI*interfaceEquations%INTERPOLATION% &
                              & VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                              & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                          ENDDO !rowParameterIdx
                        ENDDO !rowComponentIdx
                      !ENDIF !gaps(dataPointIdx)>ZERO_TOLERANCE
                    ENDDO !dataPointIdx

                    IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%FIRST_ASSEMBLY) &
                      & interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%FIRST_ASSEMBLY=.FALSE.
                  ENDIF !UPDATE_MATRIX
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Deallocate memory
                IF(ALLOCATED(orthogonallyProjected)) DEALLOCATE(orthogonallyProjected)
                IF(ALLOCATED(gapsComponents)) DEALLOCATE(gapsComponents)
                IF(ALLOCATED(gaps)) DEALLOCATE(gaps)
                
                !###################################################################################################################
                
                !Calculate penalty matrix if required
                IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD) THEN
                  interfacePenalty=>interfaceCondition%PENALTY
                  IF(ASSOCIATED(interfacePenalty)) THEN
                    penaltyField=>interfacePenalty%PENALTY_FIELD
                    IF(ASSOCIATED(penaltyField)) THEN
                      penaltyMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(interfaceEquations% &
                        & INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES)%PTR
                      IF(ASSOCIATED(penaltyMatrix)) THEN
                        IF(penaltyMatrix%FIRST_ASSEMBLY .AND. penaltyMatrix%UPDATE_MATRIX) THEN
                          DO componentIdx=1,penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%NUMBER_OF_COMPONENTS
                            SELECT CASE(penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR% &
                              & COMPONENTS(componentIdx)%INTERPOLATION_TYPE)
                            CASE(FIELD_CONSTANT_INTERPOLATION)
                              CALL FIELD_PARAMETER_SET_GET_CONSTANT(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & componentIdx,contactStiffness,err,error,*999)
                              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                                penaltyMatrix%ELEMENT_MATRIX%MATRIX(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                              ENDDO !dataPointIdx
                            CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                              localDof=penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS(componentIdx)% &
                                & PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(interfaceElementNumber)
                              CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                & localDof,contactStiffness,err,error,*999)
                              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                                penaltyMatrix%ELEMENT_MATRIX%MATRIX(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                              ENDDO !dataPointIdx
                            CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                                localDof=penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS(componentIdx)% &
                                  & PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointIdx)
                                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & localDof,contactStiffness,err,error,*999)
                                penaltyMatrix%ELEMENT_MATRIX%MATRIX(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                              ENDDO !dataPointIdx
                            CASE DEFAULT
                              localError="The interpolation type for component number "// &
                                & TRIM(NUMBER_TO_VSTRING(componentIdx,"*",err,error))// &
                                & " of variable type "//TRIM(NUMBER_TO_VSTRING(FIELD_U_VARIABLE_TYPE,"*",err,error))// &
                                & " of field number "//TRIM(NUMBER_TO_VSTRING(penaltyField%USER_NUMBER,"*",err,error))//" is "// &
                                & TRIM(NUMBER_TO_VSTRING(penaltyField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS &
                                & (componentIdx)%INTERPOLATION_TYPE,"*", err,error))// " which is invalid for penalty field."
                              CALL FlagError(localError,err,error,*999)
                            END SELECT
                          ENDDO !componentIdx
                        ENDIF              
                      ELSE
                        CALL FlagError("Interface penalty matrix is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Interface penalty field is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Interface penalty is not associated.",err,error,*999)
                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("Interface points connectivity is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
                & "*",err,error))// " is not valid."
              CALL FlagError(localError,err,error,*999)
            END SELECT !interfaceCondition%integrationType
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
              & " is not valid."
            CALL FlagError(localError,err,error,*999)
          END SELECT !interfaceCondition%METHOD
        ELSE
          CALL FlagError("Interface is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Interface equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Interface condition is not associated.",err,error,*999)
    ENDIF

    EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FrictionlessContact_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FrictionlessContact_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the element stiffness matrices for the given element number for fluid solid operator 
  !Note First interface matrix must be fluid equations set's interface matrix
  !SUBROUTINE FluidSolidOperator_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)

  !
  !================================================================================================================================
  !
  
  !>Calculates the element stiffness matrices for the given element number for solid fluid operator 
  !Note First interface matrix must be solid equations set's interface matrix
  SUBROUTINE SolidFluidOperator_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)
  
    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: interfaceEquations !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: GaussPoint, rowComponentIdx, rowIdx, rowParameterIdx, colComponentIdx, colIdx, colParameterIdx
    INTEGER(INTG) :: rowMeshComponentNumber,derivativeIdx,derivative,localElementNode,interfaceNode,interfaceDerivative
    INTEGER(INTG) :: coupledMeshElementNumber,coupledMeshIdx,coupledMeshVariableType,lagrangeVariableType
    INTEGER(INTG) :: connectedLine,decompositionLineNumber,localLineNodeIdx,connectedFace,decompositionFaceNumber,localFaceNodeIdx
    REAL(DP) :: XI(3),rwg,PGMSI,PGNSI,matrixCoefficient
    TYPE(BASIS_TYPE), POINTER :: interfaceDependentBasis,coupledMeshBasis,interfaceGeometricBasis, &
      & interfaceConnectivityBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: interfaceQuadratureScheme
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField,interfaceDependentField,interfaceGeometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: interfaceMatrixVariable,lagrangeVariable
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), POINTER :: interfaceInterpolation
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity
    TYPE(DOMAIN_LINE_TYPE), POINTER :: coupledMeshDomainLine
    TYPE(DOMAIN_FACE_TYPE), POINTER :: coupledMeshDomainFace
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolidFluidOperator_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE_EQUATIONS)) CALL FlagError("Interface equations is not associated." &
      & ,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition%INTERFACE)) CALL FlagError("Interface is not associated.",err,error,*999)

    interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS

    !===============================================================================================================================
    !Select Interface method
    SELECT CASE(interfaceCondition%METHOD)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
      !=============================================================================================================================
      !Select Integration type
      SELECT CASE(interfaceCondition%integrationType)
      CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
        !Pointers to interface variables (columns of interface element matrix)
        interfaceInterpolation=>interfaceEquations%INTERPOLATION%INTERFACE_INTERPOLATION
        interfaceGeometricField=>interfaceInterpolation%GEOMETRIC_FIELD
        interfaceDependentField=>interfaceInterpolation%DEPENDENT_FIELD
        interfaceGeometricBasis=>interfaceGeometricField%DECOMPOSITION%DOMAIN(interfaceGeometricField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        interfaceDependentBasis=>interfaceDependentField%DECOMPOSITION%DOMAIN(interfaceDependentField% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
        !Integrate using the interface quadrature scheme
        interfaceQuadratureScheme=>interfaceGeometricBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
        lagrangeVariable=>interfaceEquations%INTERFACE_MAPPING%LAGRANGE_VARIABLE
        !lagrangeVariableNumberOfComponents=>interfaceEquations%INTERFACE_MAPPING%LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
        lagrangeVariableType=lagrangeVariable%VARIABLE_TYPE
        !Get element interpolation parameters from the first geometric interpolation set (to get Jacobian for interface surface integral)
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,interfaceInterpolation% &
          & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
        matrixCoefficient=1.0_DP
        !Loop over interface matrices (1st solid, 2nd fluid)
        DO coupledMeshIdx=1,interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
          IF(interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%UPDATE_MATRIX) THEN
            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
            IF(coupledMeshIdx>1) THEN
              matrixCoefficient=-1.0_DP
            ENDIF 
            !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
            coupledMeshDependentField=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_FIELD
            elementConnectivity=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%ELEMENT_CONNECTIVITY(elementNumber,coupledMeshIdx)
            coupledMeshElementNumber=elementConnectivity%COUPLED_MESH_ELEMENT_NUMBER
            interfaceMatrixVariable=> &
              & interfaceEquations%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(coupledMeshIdx)%VARIABLE
            coupledMeshVariableType=interfaceMatrixVariable%VARIABLE_TYPE
            interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
            interfaceConnectivityBasis=>interfaceCondition%INTERFACE%MESH_CONNECTIVITY%BASIS

            !coupledMeshDependentInterpolation=>interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
            !  & DEPENDENT_INTERPOLATION
            
            !=======================================================================================================================
            !Loop over gauss points
            DO GaussPoint=1,interfaceQuadratureScheme%NUMBER_OF_GAUSS
              !CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint, &
              !  & coupledMeshDependentInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR, &
              !  & err,error,*999)
              !CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
              !  & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              !=====================================================================================================================
              !Interpolates field at given gauss point, includes first partial derivatives
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,GaussPoint,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              !Calculates the interpolated point metrics and the associated interpolated point
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(interfaceGeometricBasis%NUMBER_OF_XI,interfaceInterpolation% &
                & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
              !=====================================================================================================================
              ! R W G = GAUSSWEIGTHS * JACOBIAN
              rwg=interfaceInterpolation%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
                & JACOBIAN*interfaceQuadratureScheme%GAUSS_WEIGHTS(GaussPoint)
              IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                  & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
                CALL FlagError("Not implemented.",err,error,*999)
              ELSE
                !===================================================================================================================
                !\todo defaults to first mesh component, generalise
                !TODO Originally XI=...
                XI(1:SIZE(elementConnectivity%XI,1))=InterfaceOperators_InterfToCoupledMeshGaussTransform( &
                  & elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
                ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS

                  SELECT CASE(interfaceDependentBasis%NUMBER_OF_XI)

                  CASE(1) !1D interface (line)
                    connectedLine=elementConnectivity%CONNECTED_LINE
                    decompositionLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                    coupledMeshDomainLine=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                      & LINES%LINES(decompositionLineNumber)
                    DO localLineNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                      localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                      DO derivativeIdx=1,coupledMeshDomainLine%BASIS%NUMBER_OF_DERIVATIVES(localLineNodeIdx)
                      !???????????????????
                        derivative=coupledMeshBasis%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                        derivative=coupledMeshDomainLine%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                      !???????????????????
                        rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        !===========================================================================================================
                        ! P G M S I - this represents the D E P E N D E N T _ F I E L D S (solid, fluid)
                        !Evaluates the appropriate partial derivative index at position XI for the basis
                        PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV,XI,err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                          DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                            !=======================================================================================================
                            ! P G N S I - this represents the L A M B D A
                            PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                            colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                              & PGNSI*PGMSI*rwg*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !localLineNodeIdx

                  CASE(2) !2D interface (face)

                    SELECT CASE(coupledMeshBasis%NUMBER_OF_XI)

                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNode=1,coupledMeshBasis%NUMBER_OF_NODES
                        DO derivative=1,coupledMeshBasis%NUMBER_OF_DERIVATIVES(localElementNode)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          !=========================================================================================================
                          ! P G M S I
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              !=====================================================================================================
                              ! P G N S I
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivative
                      ENDDO !localElementNode

                    CASE(3) !Coupled Mesh has 3 xi directions
                      connectedFace = elementConnectivity%CONNECTED_FACE
                      decompositionFaceNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY% &
                        & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                      coupledMeshDomainFace=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% &
                        & FACES%FACES(decompositionFaceNumber)
                      DO localFaceNodeIdx=1,coupledMeshBasis%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                        localElementNode=coupledMeshBasis%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                        DO derivativeIdx=1,coupledMeshDomainFace%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                          derivative=coupledMeshBasis% &
                            & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                          rowParameterIdx=coupledMeshBasis%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                          !=========================================================================================================
                          ! P G M S I
                          PGMSI=BASIS_EVALUATE_XI(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & XI(1:coupledMeshBasis%NUMBER_OF_XI),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%NUMBER_OF_NODES
                            DO interfaceDerivative=1,interfaceDependentBasis%NUMBER_OF_DERIVATIVES(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                              !=====================================================================================================
                              ! P G N S I
                              PGNSI=interfaceQuadratureScheme%GAUSS_BASIS_FNS(colParameterIdx,NO_PART_DERIV,GaussPoint)
                              colIdx=colParameterIdx+interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)+ &
                                & PGNSI*PGMSI*rwg*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivativeIdx
                      ENDDO !FaceNodeIdx

                    END SELECT !coupledMeshBasis%NUMBER_OF_XI

                  END SELECT !interfaceDependentBasis%NUMBER_OF_XI

                ENDDO !rowComponentIdx
              ENDIF
            ENDDO !GaussPoint

            !Scale factor adjustment
            !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
            !\todo update looping of variables/components for non-zero matrix elements as done above 
            IF(interfaceCondition%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
              & coupledMeshIdx==interfaceEquations%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
              CALL FlagError("Not implemented.",err,error,*999)
            ELSE
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber, &
                  & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                  rowMeshComponentNumber=interfaceMatrixVariable%COMPONENTS(rowComponentIdx)%MESH_COMPONENT_NUMBER
                  coupledMeshBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN(rowMeshComponentNumber)%PTR%TOPOLOGY% & 
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx) * &
                        & interfaceInterpolation%DEPENDENT_INTERPOLATION(1)% &
                        & INTERPOLATION_PARAMETERS(lagrangeVariableType)%PTR%SCALE_FACTORS(colParameterIdx,colComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
              !Scale factor adjustment for the row dependent variable
              IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(coupledMeshElementNumber, &
                  & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                  & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,interfaceMatrixVariable%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%NUMBER_OF_ELEMENT_PARAMETERS
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%NUMBER_OF_COMPONENTS
                      DO colParameterIdx=1,interfaceDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        colIdx=colIdx+1
                        interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx,colIdx)* &
                        & interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                        & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(coupledMeshVariableType)%PTR% &
                        & SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! coupledMeshIdx
        
      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType

    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("SolidFluidOperator_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("SolidFluidOperator_FiniteElementCalculate",err,error)
    RETURN 1
  
  END SUBROUTINE SolidFluidOperator_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  FUNCTION InterfaceOperators_InterfToCoupledMeshGaussTransform(elementConnectivity,interfaceConnectivityBasis,GaussPoint,err,error)
  
    !Argument variables
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: elementConnectivity !<A pointer to the element connectivity
    TYPE(BASIS_TYPE), POINTER :: interfaceConnectivityBasis !<A pointer to the interface mesh connectivity basis
    INTEGER(INTG), INTENT(IN) :: GaussPoint !< Index to the gauss point which needs to be transformed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: InterfaceOperators_InterfToCoupledMeshGaussTransform(SIZE(elementConnectivity%XI,1))
    !Local Variables
    INTEGER(INTG) :: rowParameterIdx

    ENTERS("InterfaceOperators_InterfToCoupledMeshGaussTransform",err,error,*999)
    
    InterfaceOperators_InterfToCoupledMeshGaussTransform=0.0_DP
    DO rowParameterIdx = 1,interfaceConnectivityBasis%NUMBER_OF_ELEMENT_PARAMETERS
      InterfaceOperators_InterfToCoupledMeshGaussTransform(:)= InterfaceOperators_InterfToCoupledMeshGaussTransform(:) + &
        & interfaceConnectivityBasis%QUADRATURE%QUADRATURE_SCHEME_MAP &
        & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_BASIS_FNS(rowParameterIdx,NO_PART_DERIV,GaussPoint) * &
        & elementConnectivity%XI(:,1,rowParameterIdx)
    ENDDO
     
    EXITS("InterfaceOperators_InterfToCoupledMeshGaussTransform")
    RETURN
999 ERRORS("InterfaceOperators_InterfToCoupledMeshGaussTransform",err,error)
    EXITS("InterfaceOperators_InterfToCoupledMeshGaussTransform")
    RETURN
    
  END FUNCTION InterfaceOperators_InterfToCoupledMeshGaussTransform

END MODULE INTERFACE_OPERATORS_ROUTINES
