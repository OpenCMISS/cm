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

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for field continuity operator \todo: This subroutine will be completely replace by Prasad's tight coupling code
  SUBROUTINE FieldContinuity_FiniteElementCalculate(INTERFACE_CONDITION,ELEM,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: ELEM !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ng, mh, mhs, ms, nh, nhs, ns, nhc, mhc, dir
    REAL(DP) :: XI(3)
    REAL(DP) :: RWG, PGMSI, PGNSI
    INTEGER(INTG) :: MID
    INTEGER(INTG) :: FVAR_TYPE_RW, FVAR_TYPE_CL
    TYPE(BASIS_TYPE), POINTER :: COLBAS,ROWBAS,GEOMBAS
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FVAR_RW, FVAR_CL
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX, EQMAT
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(FIELD_TYPE), POINTER :: ROWVAR, COLVAR, GEOMVAR
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADSCHEME, QUADSCHEMECOL, QUADSCHEMEROW
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("FieldContinuity_FiniteElementCalculate",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTERFACE_CONDITION)) CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    IF(.NOT.ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) CALL FLAG_ERROR("Interface equations is not associated." &
      & ,ERR,ERROR,*999)

    INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
    SELECT CASE(INTERFACE_CONDITION%METHOD)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
      GEOMVAR=>INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%GEOMETRIC_FIELD
      COLVAR =>INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%DEPENDENT_FIELD
      FVAR_CL=>INTERFACE_EQUATIONS%INTERFACE_MAPPING%LAGRANGE_VARIABLE
      FVAR_TYPE_CL=FVAR_CL%VARIABLE_TYPE
      DO MID = 1,INTERFACE_CONDITION%INTERFACE%NUMBER_OF_COUPLED_MESHES
        IF (INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(MID)%PTR%UPDATE_MATRIX) THEN
          ROWVAR=>INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(MID)%DEPENDENT_FIELD
          EQMAT =>INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(MID)%PTR%ELEMENT_MATRIX
          FVAR_RW=>INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(MID)%VARIABLE
          FVAR_TYPE_RW=FVAR_RW%VARIABLE_TYPE
          GEOMBAS=>GEOMVAR%DECOMPOSITION%DOMAIN(GEOMVAR%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
             & TOPOLOGY%ELEMENTS%ELEMENTS(ELEM)%BASIS
          INTERFACE_MATRIX=>INTERFACE_EQUATIONS%INTERFACE_MATRICES%MATRICES(MID)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEM,INTERFACE_EQUATIONS%INTERPOLATION% &  
             & INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%PTR, &
             & ERR,ERROR,*999)
          QUADSCHEME=>GEOMBAS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          !Loop over gauss points
          DO ng=1,QUADSCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,INTERFACE_EQUATIONS%INTERPOLATION% &
              & INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%PTR, &
              & ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMBAS%NUMBER_OF_XI,INTERFACE_EQUATIONS%INTERPOLATION% &  
              & INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)%INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR, &
              & ERR,ERROR,*999)
            RWG=INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%GEOMETRIC_INTERPOLATION(1)% &
              & INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN * QUADSCHEME%GAUSS_WEIGHTS(ng)
            dir=SIZE(INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY%ELEMENTS_CONNECTIVITY(ELEM,MID)%XI,1)
            XI(1:dir)=INTERFACE_TRANSFORM_GPT(INTERFACE_CONDITION%INTERFACE%MESH_CONNECTIVITY,ELEM,MID,ng,ERR,ERROR)
            mhs=0
            DO mh=1,FVAR_RW%NUMBER_OF_COMPONENTS
              mhc=FVAR_RW%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              ROWBAS=>ROWVAR%DECOMPOSITION%DOMAIN(mhc)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEM)%BASIS
              QUADSCHEMEROW=>ROWBAS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
              !Loop over element rows
              DO ms=1,ROWBAS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1; nhs=0; PGMSI=BASIS_EVALUATE_XI(ROWBAS,mhs,NO_PART_DERIV,XI,ERR,ERROR)
                !Loop over element columns
                DO nh=1,FVAR_CL%NUMBER_OF_COMPONENTS
                  nhc=FVAR_CL%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                  COLBAS=>COLVAR%DECOMPOSITION%DOMAIN(nhc)%PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEM)%BASIS
                  QUADSCHEMECOL=>COLBAS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                  DO ns=1,COLBAS%NUMBER_OF_ELEMENT_PARAMETERS
                    nhs=nhs+1; PGNSI=QUADSCHEMECOL%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                    INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+PGNSI*PGMSI*RWG
                  ENDDO !ns
                ENDDO !nh
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
          !Scale factor adjustment for the Column Variable
          IF(COLVAR%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEM,INTERFACE_EQUATIONS%INTERPOLATION% &  
              & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FVAR_TYPE_CL)%PTR,ERR,ERROR,*999)
            mhs=0
            DO mh=1,FVAR_RW%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,ROWBAS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1; nhs=0
                !Loop over element columns
                DO nh=1,FVAR_CL%NUMBER_OF_COMPONENTS
                  DO ns=1,COLBAS%NUMBER_OF_ELEMENT_PARAMETERS
                    nhs=nhs+1
!                     INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs) * &
!                     & INTERFACE_EQUATIONS%INTERPOLATION%INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)% &
!                     & INTERPOLATION_PARAMETERS(FVAR_TYPE_CL)%PTR%SCALE_FACTORS(ms,mh) * INTERFACE_EQUATIONS%INTERPOLATION% &
!                     & INTERFACE_INTERPOLATION%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FVAR_TYPE_CL)%PTR% &
!                     & SCALE_FACTORS(ns,nh)
                  ENDDO !ns
                ENDDO !nh
              ENDDO !ms
            ENDDO !mh
          ENDIF
          !Scale factor adjustment for the Column Variable
          IF(ROWVAR%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEM_GET(ELEM,INTERFACE_EQUATIONS%INTERPOLATION% &  
              & VARIABLE_INTERPOLATION(MID)%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FVAR_TYPE_RW)%PTR,ERR,ERROR,*999)
            mhs=0
            DO mh=1,FVAR_RW%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,ROWBAS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1; nhs=0
                !Loop over element columns
                DO nh=1,FVAR_CL%NUMBER_OF_COMPONENTS
                  DO ns=1,COLBAS%NUMBER_OF_ELEMENT_PARAMETERS
                    nhs=nhs+1
!                     INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=INTERFACE_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
!                     & INTERFACE_EQUATIONS%INTERPOLATION%VARIABLE_INTERPOLATION(MID)%DEPENDENT_INTERPOLATION(1)% & 
!                     & INTERPOLATION_PARAMETERS(FVAR_TYPE_RW)%PTR%SCALE_FACTORS(ms,mh) * INTERFACE_EQUATIONS%INTERPOLATION% &
!                     & VARIABLE_INTERPOLATION(MID)%DEPENDENT_INTERPOLATION(1)%INTERPOLATION_PARAMETERS(FVAR_TYPE_RW)%PTR% &
!                     & SCALE_FACTORS(ns,nh)
                  ENDDO !ns
                ENDDO !nh
              ENDDO !ms
            ENDDO !mh
          ENDIF
        ENDIF ! UPDATE MATRIX
      ENDDO ! MID
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
    CASE(INTERFACE_CONDITION_PENALTY_METHOD)
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

  FUNCTION INTERFACE_TRANSFORM_GPT(INTERFACE_MESH_CONNECTIVITY,ELEMENT_NUMBER,DOMAIN_NUMBER,GAUSSPT,ERR,ERROR)
  
    !Argument variables
    TYPE(INTERFACE_MESH_CONNECTIVITY_TYPE), POINTER :: INTERFACE_MESH_CONNECTIVITY !<A pointer to the interface meshes connectivity to finish creating
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !< Index to the element number
    INTEGER(INTG), INTENT(IN) :: DOMAIN_NUMBER !< Index to the domain number
    INTEGER(INTG), INTENT(IN) :: GAUSSPT !< Index to the element coupled list
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    REAL(DP) :: INTERFACE_TRANSFORM_GPT(SIZE(INTERFACE_MESH_CONNECTIVITY%ELEMENTS_CONNECTIVITY &
      & (ELEMENT_NUMBER,DOMAIN_NUMBER)%XI,1))
    !Local Variables
    INTEGER(INTG) :: ms

    CALL ENTERS("INTERFACE_TRANSFORM_GPT",ERR,ERROR,*999)
    
    INTERFACE_TRANSFORM_GPT=0.0_DP
    IF(.NOT.ASSOCIATED(INTERFACE_MESH_CONNECTIVITY)) CALL FLAG_ERROR("Mesh Connectivity is not associated.",ERR,ERROR,*999)
    IF(.NOT.ALLOCATED(INTERFACE_MESH_CONNECTIVITY%ELEMENTS_CONNECTIVITY(ELEMENT_NUMBER,DOMAIN_NUMBER)%XI)) THEN
      CALL FLAG_ERROR("Coupled Mesh xi array not allocated.",ERR,ERROR,*999)
    END IF

    DO ms = 1,INTERFACE_MESH_CONNECTIVITY%BASIS%NUMBER_OF_ELEMENT_PARAMETERS
      INTERFACE_TRANSFORM_GPT(:)= INTERFACE_TRANSFORM_GPT(:) + INTERFACE_MESH_CONNECTIVITY%BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP &
         & (BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,GAUSSPT) * INTERFACE_MESH_CONNECTIVITY% & 
         & ELEMENTS_CONNECTIVITY(ELEMENT_NUMBER,DOMAIN_NUMBER)%XI(:,1,ms)
    END DO
     
    CALL EXITS("INTERFACE_TRANSFORM_GPT")
    RETURN
999 CALL ERRORS("INTERFACE_TRANSFORM_GPT",ERR,ERROR)
    CALL EXITS("INTERFACE_TRANSFORM_GPT")
    RETURN
    
  END FUNCTION INTERFACE_TRANSFORM_GPT
  
  !
  !================================================================================================================================
  !

END MODULE INTERFACE_OPERATORS_ROUTINES
