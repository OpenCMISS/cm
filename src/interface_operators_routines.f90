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

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FieldContinuity_FiniteElementCalculate
  
  PUBLIC FrictionlessContact_FiniteElementCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for field continuity operator 
  SUBROUTINE FieldContinuity_FiniteElementCalculate(INTERFACE_CONDITION,ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: ELEMENT !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ng, mh, mhs, ms, nh, nhs, ns, mhc, dir,derivativeIdx,derivative,localNode,localNodeIdx,localElementNode
    INTEGER(INTG) :: connectedLine,lineNodeIdx,decompositionLineNumber,localLineNode,localLineNodeIdx
    INTEGER(INTG) :: connectedFace,faceNodeIdx,decompositionFaceNumber,localFaceNode,localFaceNodeIdx
    INTEGER(INTG) :: interfaceNode,interfaceDerivative,coupledMeshElementNumber
    REAL(DP) :: XI(3),RWG,PGMSI,PGNSI,MATRIX_COEFFICIENT
    INTEGER(INTG) :: interface_matrix_idx
    INTEGER(INTG) :: INTERFACE_MATRIX_VARIABLE_TYPE,LAGRANGE_VARIABLE_TYPE
    TYPE(BASIS_TYPE), POINTER :: INTERFACE_DEPENDENT_BASIS,COUPLED_MESH_BASIS,INTERFACE_GEOMETRIC_BASIS, &
      & INTERFACE_PENALTY_BASIS,INTERFACE_CONNECTIVITY_BASIS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: INTERFACE_MATRIX_VARIABLE,LAGRANGE_VARIABLE
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX,ELEMENT_MATRIX,INTERFACE_ELEMENT_MATRIX
    TYPE(INTERFACE_EQUATIONS_DOMAIN_INTERPOLATION_TYPE), POINTER :: INTERFACE_INTERPOLATION
    TYPE(FIELD_TYPE), POINTER :: INTERFACE_MATRIX_DEPENDENT_FIELD,INTERFACE_DEPENDENT_FIELD,INTERFACE_GEOMETRIC_FIELD, &
      & INTERFACE_PENALTY_FIELD
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: INTERFACE_QUADRATURE_SCHEME
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY
    TYPE(DOMAIN_LINE_TYPE), POINTER :: COUPLED_MESH_DOMAIN_LINE
    TYPE(DOMAIN_FACE_TYPE), POINTER :: COUPLED_MESH_DOMAIN_FACE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FieldContinuity_FiniteElementCalculate",ERR,ERROR,*999)

    IF(.NOT.ASSOCIATED(INTERFACE_CONDITION)) CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) CALL FLAG_ERROR("Interface equations is not associated." &
      & ,ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTERFACE_CONDITION%INTERFACE)) CALL FLAG_ERROR("Interface is not associated.",err,error,*999)

    INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
    SELECT CASE(INTERFACE_CONDITION%METHOD)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      !Pointers to interface variables (columns of interface element matrix)
      INTERFACE_INTERPOLATION=>INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION
      INTERFACE_GEOMETRIC_FIELD=>INTERFACE_INTERPOLATION%GEOMETRIC_FIELD
      INTERFACE_GEOMETRIC_BASIS=>INTERFACE_GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(INTERFACE_GEOMETRIC_FIELD% &
        & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT)%BASIS
      INTERFACE_DEPENDENT_FIELD=>INTERFACE_INTERPOLATION%DEPENDENT_FIELD
      INTERFACE_DEPENDENT_BASIS=>INTERFACE_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(INTERFACE_DEPENDENT_FIELD% &
        & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT)%BASIS
      SELECT CASE(INTERFACE_CONDITION%METHOD)
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        INTERFACE_PENALTY_FIELD=>INTERFACE_INTERPOLATION%PENALTY_FIELD
        INTERFACE_PENALTY_BASIS=>INTERFACE_PENALTY_FIELD%DECOMPOSITION%DOMAIN(INTERFACE_PENALTY_FIELD% &
          & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT)%BASIS
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT,INTERFACE_INTERPOLATION% &
          & PENALTY_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
      ENDSELECT
      !Integrate using the interface quadrature scheme
      INTERFACE_QUADRATURE_SCHEME=>INTERFACE_GEOMETRIC_BASIS%QUADRATURE% &
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
      LAGRANGE_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING%LAGRANGE_VARIABLE
      LAGRANGE_VARIABLE_TYPE=LAGRANGE_VARIABLE%VARIABLE_TYPE
      !Get element interpolation parameters from the first geometric interpolation set (to get Jacobian for interface surface integral)
      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT,INTERFACE_INTERPOLATION% &
        & GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
      !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
      MATRIX_COEFFICIENT=1.0_DP
      DO interface_matrix_idx=1,INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
        IF(INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%UPDATE_MATRIX) THEN
          !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
          IF(interface_matrix_idx>1) THEN
            MATRIX_COEFFICIENT=-1.0_DP
          ENDIF
          !Pointers to the interface_matrix_idx'th coupled mesh variables (rows of interface element matrix)
          INTERFACE_MATRIX_DEPENDENT_FIELD=>INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
            & DEPENDENT_FIELD
          INTERFACE_MATRIX_VARIABLE=>INTERFACE_EQUATIONS%INTERFACE_MAPPING% & 
            & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
          INTERFACE_MATRIX_VARIABLE_TYPE=INTERFACE_MATRIX_VARIABLE%VARIABLE_TYPE
          INTERFACE_ELEMENT_MATRIX=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(interface_matrix_idx)%PTR%ELEMENT_MATRIX
          ELEMENT_CONNECTIVITY=>INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY% &
            & ELEMENT_CONNECTIVITY(ELEMENT,interface_matrix_idx)
          coupledMeshElementNumber=ELEMENT_CONNECTIVITY%COUPLED_MESH_ELEMENT_NUMBER
          INTERFACE_CONNECTIVITY_BASIS=>INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY%BASIS

          !Loop over gauss points
          DO ng=1,INTERFACE_QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,INTERFACE_INTERPOLATION% &
              & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(INTERFACE_GEOMETRIC_BASIS%NUMBER_OF_XI,INTERFACE_INTERPOLATION% &
              & GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
            RWG=INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR% &
              & JACOBIAN*INTERFACE_QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
                & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,INTERFACE_INTERPOLATION% &
                & PENALTY_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              mhs=0
              DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over the Lagrange variable matrix rows
                DO ms=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  mhs=mhs+1
                  INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs)- &
                    & (1.0_DP/INTERFACE_INTERPOLATION%PENALTY_INTERPOLATION(1)% &
                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,1))*PGNSI**2.0_DP*RWG
                ENDDO !ms
              ENDDO !mh
            ELSE
              !\todo defaults to first mesh component, Generalise
              XI(1:INTERFACE_DEPENDENT_BASIS%NUMBER_OF_XI)=INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM( &
                & ELEMENT_CONNECTIVITY,INTERFACE_CONNECTIVITY_BASIS,ng,ERR,ERROR)
              ! Loop over number of Lagrange variable components as not all components in the dependent field variable may be coupled
              !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
              DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                mhc=INTERFACE_MATRIX_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% & 
                  & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                SELECT CASE(INTERFACE_DEPENDENT_BASIS%NUMBER_OF_XI)
                CASE(1) !1D interface (line)
                  connectedLine = ELEMENT_CONNECTIVITY%CONNECTED_LINE
                  decompositionLineNumber=INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY% &
                    & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_LINES(connectedLine)
                  COUPLED_MESH_DOMAIN_LINE=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% &
                    & LINES%LINES(decompositionLineNumber)
                  DO localLineNodeIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(connectedLine)
                    localElementNode=COUPLED_MESH_BASIS%NODE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                    DO derivativeIdx=1,COUPLED_MESH_DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(localLineNodeIdx)
                      derivative=COUPLED_MESH_BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(localLineNodeIdx,connectedLine)
                      derivative=COUPLED_MESH_DOMAIN_LINE%DERIVATIVES_IN_LINE(1,derivativeIdx,localLineNodeIdx)
                      ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                      IF (mh==4) THEN
                        PGMSI=1.0_DP
                      ELSE
                        PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,ms,NO_PART_DERIV,XI,ERR,ERROR)
                      ENDIF
                      mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                      DO interfaceNode=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_NODES
                        DO interfaceDerivative=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_DERIVATIVES(interfaceNode)
                          !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                          ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                          PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                          nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                          !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                          INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                            & PGNSI*PGMSI*RWG*MATRIX_COEFFICIENT
                        ENDDO !interfaceDerivative
                      ENDDO !interfaceNode
                    ENDDO !derivativeIdx
                  ENDDO !lineNodeIdx
                CASE(2) !2D interface (face)
                  SELECT CASE(COUPLED_MESH_BASIS%NUMBER_OF_XI)
                  CASE(2) !Coupled Mesh has 2 xi directions
                    DO localElementNode=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES
                      DO derivative=1,COUPLED_MESH_BASIS%NUMBER_OF_DERIVATIVES(localElementNode)
                        ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        IF (mh==4) THEN
                          PGMSI=1.0_DP
                        ELSE
                          PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,ms,NO_PART_DERIV, &
                            & XI(1:COUPLED_MESH_BASIS%NUMBER_OF_XI),ERR,ERROR)
                        ENDIF
                        mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                        DO interfaceNode=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_NODES
                          DO interfaceDerivative=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_DERIVATIVES(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                            PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                            nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                            INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                              & PGNSI*PGMSI*RWG*MATRIX_COEFFICIENT
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivative
                    ENDDO !localElementNode
                  CASE(3) !Coupled Mesh has 3 xi directions
                    connectedFace = ELEMENT_CONNECTIVITY%CONNECTED_FACE
                    decompositionFaceNumber=INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%TOPOLOGY% &
                      & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%ELEMENT_FACES(connectedFace)
                    COUPLED_MESH_DOMAIN_FACE=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% &
                      & FACES%FACES(decompositionFaceNumber)
                    DO localFaceNodeIdx=1,COUPLED_MESH_BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(connectedFace)
                      localElementNode=COUPLED_MESH_BASIS%NODE_NUMBERS_IN_LOCAL_FACE(localFaceNodeIdx,connectedFace)
                      DO derivativeIdx=1,COUPLED_MESH_DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(localFaceNodeIdx)
                        derivative=COUPLED_MESH_BASIS% &
                          & DERIVATIVE_NUMBERS_IN_LOCAL_FACE(derivativeIdx,localFaceNodeIdx,connectedFace)
                        ms=COUPLED_MESH_BASIS%ELEMENT_PARAMETER_INDEX(derivative,localElementNode)
                        IF (mh==4) THEN
                          PGMSI=1.0_DP
                        ELSE
                          PGMSI=BASIS_EVALUATE_XI(COUPLED_MESH_BASIS,ms,NO_PART_DERIV, &
                            & XI(1:COUPLED_MESH_BASIS%NUMBER_OF_XI),ERR,ERROR)
                        ENDIF
                        mhs=ms+COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                        DO interfaceNode=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_NODES
                          DO interfaceDerivative=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_DERIVATIVES(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            ns=INTERFACE_DEPENDENT_BASIS%ELEMENT_PARAMETER_INDEX(interfaceDerivative,interfaceNode)
                            PGNSI=INTERFACE_QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                            nhs=ns+INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS*(mh-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
                            INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                              & PGNSI*PGMSI*RWG*MATRIX_COEFFICIENT
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !FaceNodeIdx
                  END SELECT
                END SELECT
              ENDDO !mh
            ENDIF
          ENDDO !ng

          !Scale factor adjustment
          !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix contribution to the residual for non-linear problems
          !\todo update looping of variables/components for non-zero matrix elements as done above 
          IF(INTERFACE_CONDITION%METHOD==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
            & interface_matrix_idx==INTERFACE_EQUATIONS%INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES) THEN
            !Scale factor adjustment for the Lagrange Variable (columns)
            IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT, &
                & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR, &
                & ERR,ERROR,*999)
              mhs=0
              !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
              !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
              DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element Lagrange variable rows
                DO ms=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,mhs) * &
                    & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                    & INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR%SCALE_FACTORS(ms,mh)**2
                ENDDO !ms
              ENDDO !mh
            ENDIF
          ELSE
            !Scale factor adjustment for the Lagrange Variable (columns)
            IF(INTERFACE_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEMENT, &
                & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR, &
                & ERR,ERROR,*999)
              mhs=0
              !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
              !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable component numbers. Generalise ordering
              DO mh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                mhc=INTERFACE_MATRIX_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                COUPLED_MESH_BASIS=>INTERFACE_MATRIX_DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY% & 
                  & ELEMENTS%ELEMENTS(coupledMeshElementNumber)%BASIS
                !Loop over element rows
                DO ms=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs) * &
                      & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
                      & INTERPOLATION_PARAMETERS(LAGRANGE_VARIABLE_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh
            ENDIF
            !Scale factor adjustment for the row dependent variable
            IF(INTERFACE_MATRIX_DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(coupledMeshElementNumber, &
                & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
              mhs=0
              DO mh=1,INTERFACE_MATRIX_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,COUPLED_MESH_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,LAGRANGE_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,INTERFACE_DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                      & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(interface_matrix_idx)% &
                      & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(INTERFACE_MATRIX_VARIABLE_TYPE)%PTR% &
                      & SCALE_FACTORS(ms,mh)
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh
            ENDIF
          ENDIF
        ENDIF
      ENDDO ! interface_matrix_idx
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE DEFAULT
      LOCAL_ERROR="Interface condition method "//TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
        & " is not valid."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    END SELECT

    CALL EXITS("FieldContinuity_FiniteElementCalculate")
    RETURN
999 CALL ERRORS("FieldContinuity_FiniteElementCalculate",ERR,ERROR)
    CALL EXITS("FieldContinuity_FiniteElementCalculate")
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
    TYPE(FIELD_TYPE), POINTER :: coupledMeshDependentField
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointsMetrics(:)
    TYPE(BASIS_TYPE), POINTER :: coupledMeshDependentBasis
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: interfaceElementMatrix
    INTEGER(INTG) :: meshComponentNumber,numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements
    INTEGER(INTG) :: dataPointIdx,coupledMeshIdx,xiIdx,localElementNumber,localFaceLineNumber,matrixElementIdx,rowComponentIdx, &
      & rowParameterIdx,rowIdx,colIdx
    INTEGER(INTG) :: matrixCoefficients(2)
    REAL(DP) :: PGMSI
    REAL(DP) :: positionPoint(3),normalPoint(3),tangentsPoint(3,3),xi(3)
    REAL(DP), ALLOCATABLE :: gaps(:),gapsComponents(:,:),normals(:,:)
    LOGICAL, ALLOCATABLE :: orthogonallyProjected(:)
    
    
    TYPE(VARYING_STRING) :: localError

    CALL ENTERS("FrictionlessContact_FiniteElementCalculate",err,error,*999)
    
    IF(ASSOCIATED(interfaceCondition)) THEN
      interfaceEquations=>interfaceCondition%INTERFACE_EQUATIONS
      IF(ASSOCIATED(interfaceEquations)) THEN
        interface=>interfaceCondition%INTERFACE
        IF(ASSOCIATED(interface)) THEN
          SELECT CASE(interfaceCondition%METHOD)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
            SELECT CASE(interfaceCondition%integrationType)
            CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
              CALL FLAG_ERROR("Mesh connectivity is not implemented for frictionless contact.",err,error,*999)
            CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
              matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
              matrixCoefficients(2)=-1;
              pointsConnectivity=>interface%pointsConnectivity
              numberOfInterfaceMeshXi=pointsConnectivity%interfaceMesh%NUMBER_OF_DIMENSIONS
              IF(ASSOCIATED(pointsConnectivity)) THEN
                decompositionElementData=>interfaceCondition%LAGRANGE%LAGRANGE_FIELD%DECOMPOSITION%TOPOLOGY%dataPoints% &
                  & elementDataPoint(interfaceElementNumber)
                !###################################################################################################################
                
                !Test for orthogonal projected
                ALLOCATE(orthogonallyProjected(decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate orthogonal projected logicals.",err,error,*999)
                orthogonallyProjected=.TRUE. !Initialise orthogonal projected logicals
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                    & DEPENDENT%DEPENDENT_FIELD
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    IF(ALL(pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%reducedXi == 0.0_DP)) THEN
                      orthogonallyProjected(dataPointIdx)=.FALSE.
                    ENDIF
                  ENDDO !dataPointIdx
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Allocate memory for local allocatable variables
                ALLOCATE(gaps(decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate gaps.",err,error,*999)
                gaps=0.0_DP !Initialise gap functions
                ALLOCATE(gapsComponents(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate component gaps.",err,error,*999)
                gapsComponents=0.0_DP !Initialise gap functions
                ALLOCATE(normals(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate normals.",err,error,*999)
                normals=0.0_DP !Initialise gap functions
                
                !Calculate Gap for each data point 
                DO coupledMeshIdx=1,interface%NUMBER_OF_COUPLED_MESHES
                  coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                    & DEPENDENT%DEPENDENT_FIELD
                  NULLIFY(interpolatedPoints)
                  NULLIFY(interpolationParameters)
                  CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(coupledMeshDependentField,interpolationParameters,err,error,*999)
                  CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999)
                  interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%PTR
                  numberOfCoupledMeshGeoComp=coupledMeshDependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    !Only interpolate if orthogonally projected
                    IF(orthogonallyProjected(dataPointIdx)) THEN
                      localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%coupledMeshElementNumber
                      localFaceLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)% &
                        & ELEMENT_FACES(pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%elementLineFaceNumber)
                      SELECT CASE(numberOfInterfaceMeshXi) !Use face/line interpolation parameters for normal calculation
                      CASE(1)
                        CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                      CASE(2)
                        CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                          & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                      END SELECT
                      CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,pointsConnectivity%pointsConnectivity(dataPointIdx, &
                        & coupledMeshIdx)%reducedXi(:),interpolatedPoint,err,error,*999) !Interpolate contact data points on each surface
                      gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx)=gapsComponents(1:numberOfCoupledMeshGeoComp, &
                        & dataPointIdx)+interpolatedPoint%VALUES(1:numberOfCoupledMeshGeoComp,NO_PART_DERIV)* &
                        & matrixCoefficients(coupledMeshIdx) !Calculate 3 components gap function for each contact point
                      !Calculate surface normal (use 1st coupled mesh surface normal)
                      IF (coupledMeshIdx==1) THEN
                        CALL FIELD_INTERPOLATED_POINTS_METRICS_INITIALISE(interpolatedPoints,interpolatedPointsMetrics, &
                          & err,error,*999)
                        CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(numberOfCoupledMeshGeoComp,interpolatedPointsMetrics &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                        CALL FIELD_POSITION_NORMAL_TANGENTS_CALCULATE_INT_PT_METRIC(interpolatedPointsMetrics &
                          & (FIELD_U_VARIABLE_TYPE)%PTR,positionPoint,normalPoint,tangentsPoint,err,error,*999)
                        normals(1:numberOfCoupledMeshGeoComp,dataPointIdx)=normalPoint(1:numberOfCoupledMeshGeoComp)
                        CALL FIELD_INTERPOLATED_POINTS_METRICS_FINALISE(interpolatedPointsMetrics,err,error,*999)
                      ENDIF !coupledMeshIdx==1
                    ENDIF !orthogonallyProjected(dataPointIdx)
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
                    numberOfCoupledMeshGeoComp=coupledMeshDependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)%NUMBER_OF_COMPONENTS
                    coupledMeshDependentField=>interfaceCondition%DEPENDENT%EQUATIONS_SETS(coupledMeshIdx)%PTR% &
                      & DEPENDENT%DEPENDENT_FIELD
                    interfaceElementMatrix=>interfaceEquations%INTERFACE_MATRICES%MATRICES(coupledMeshIdx)%PTR%ELEMENT_MATRIX
                    !mesh component number is the same for all geometric components in elasticity problems
                    meshComponentNumber=coupledMeshDependentField%VARIABLES(FIELD_U_VARIABLE_TYPE)%COMPONENTS(1)% &
                      & MESH_COMPONENT_NUMBER
                    DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                      IF(gaps(dataPointIdx)>0.0_dp) THEN !Only add contact point contribution if the gap is a penetration
                        localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                          & coupledMeshElementNumber
                        !Calculate the element index (non-conforming element) for this interface matrix
                        matrixElementIdx=1
                        DO WHILE (localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                            & elementNumbers(matrixElementIdx))
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
                              & xi(1:numberOfCoupledMeshXi),ERR,ERROR)*normals(rowComponentIdx,dataPointIdx)* &
                              & matrixCoefficients(coupledMeshIdx)
                            rowIdx=coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS*numberOfMatrixCoupledElements* &
                              & (rowComponentIdx-1)+coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS* &
                              & (matrixElementIdx-1)+rowParameterIdx
                            colIdx=dataPointIdx
                            interfaceElementMatrix%MATRIX(rowIdx,colIdx)=PGMSI !Update interface element matrix with contact point contribution
                          ENDDO !rowParameterIdx
                        ENDDO !rowComponentIdx
                      ENDIF !gaps(dataPointIdx)>0.0_dp
                    ENDDO !dataPointIdx
                    !scale factor update
                    IF(coupledMeshDependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                      rowIdx=0
                      DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                        DO matrixElementIdx=1,numberOfMatrixCoupledElements
                          localElementNumber=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                          & elementNumbers(matrixElementIdx)
                          coupledMeshDependentBasis=>coupledMeshDependentField%DECOMPOSITION%DOMAIN &
                            & (meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
                          CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(localElementNumber,interfaceEquations% &
                            & INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)%DEPENDENT_INTERPOLATION(1)% &
                            & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR,ERR,ERROR,*999)
                          DO rowParameterIdx=1,coupledMeshDependentBasis%NUMBER_OF_ELEMENT_PARAMETERS 
                            rowIdx=rowIdx+1
                            DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                              colIdx=dataPointIdx
                              interfaceElementMatrix%MATRIX(rowIdx,colIdx)=interfaceElementMatrix%MATRIX(rowIdx, &
                                & colIdx)*interfaceEquations%INTERPOLATION%VARIABLE_INTERPOLATION(coupledMeshIdx)% &
                                & DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)% &
                                & PTR%SCALE_FACTORS(rowParameterIdx,rowComponentIdx)
                            ENDDO !dataPointIdx
                          ENDDO !rowParameterIdx 
                        ENDDO !coupledMeshElementIdx
                      ENDDO !rowComponentIdx
                    ENDIF !.NOT. FIELD_NO_SCALING
                  ENDIF !UPDATE_MATRIX
                ENDDO !coupledMeshIdx
                
                !###################################################################################################################
                
                !Deallocate memory
                IF(ALLOCATED(orthogonallyProjected)) DEALLOCATE(orthogonallyProjected)
                IF(ALLOCATED(gapsComponents)) DEALLOCATE(gapsComponents)
                IF(ALLOCATED(gaps)) DEALLOCATE(gaps)
              ELSE
                CALL FLAG_ERROR("Interface points connectivity is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="Interface condition integration type "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%integrationType, &
                & "*",err,error))// " is not valid."
              CALL FLAG_ERROR(localError,err,error,*999)
            END SELECT !interfaceCondition%integrationType
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="Interface condition method "//TRIM(NUMBER_TO_VSTRING(interfaceCondition%METHOD,"*",err,error))// &
              & " is not valid."
            CALL FLAG_ERROR(localError,err,error,*999)
          END SELECT !interfaceCondition%METHOD
        ELSE
          CALL FLAG_ERROR("Interface is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN
999 CALL ERRORS("FrictionlessContact_FiniteElementCalculate",err,error)
    CALL EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE FrictionlessContact_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  FUNCTION INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(ELEMENT_CONNECTIVITY,INTERFACE_CONNECTIVITY_BASIS,GAUSS_POINT,ERR,ERROR)
  
    !Argument variables
    TYPE(INTERFACE_ELEMENT_CONNECTIVITY_TYPE), POINTER :: ELEMENT_CONNECTIVITY !<A pointer to the element connectivity
    TYPE(BASIS_TYPE), POINTER :: INTERFACE_CONNECTIVITY_BASIS !<A pointer to the interface mesh connectivity basis
    INTEGER(INTG), INTENT(IN) :: GAUSS_POINT !< Index to the gauss point which needs to be transformed
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(SIZE(ELEMENT_CONNECTIVITY%XI,1))
    !Local Variables
    INTEGER(INTG) :: ms

    CALL ENTERS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM",ERR,ERROR,*999)
    
    INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM=0.0_DP
    DO ms = 1,INTERFACE_CONNECTIVITY_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
      INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:)= INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM(:) + &
        & INTERFACE_CONNECTIVITY_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP &
        & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,GAUSS_POINT) * ELEMENT_CONNECTIVITY%XI(:,1,ms)
    ENDDO
     
    CALL EXITS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM")
    RETURN
999 CALL ERRORS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM",ERR,ERROR)
    CALL EXITS("INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM")
    RETURN
    
  END FUNCTION INTERFACE_TO_COUPLED_MESH_GAUSSPOINT_TRANSFORM

END MODULE INTERFACE_OPERATORS_ROUTINES
