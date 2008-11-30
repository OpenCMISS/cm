!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all Laplace equations routines.
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

!>This module handles all Laplace equations routines.
MODULE LAPLACE_EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITION_ROUTINES
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
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
  USE SOLUTION_MAPPING_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

!!MERGE: move

  PUBLIC LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE,LAPLACE_EQUATION_EQUATIONS_SET_SETUP, &
    & LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET,LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET,LAPLACE_EQUATION_PROBLEM_SETUP
  
CONTAINS

  !
  !================================================================================================================================
  !

!!MERGE: check

  !>Calculates the element stiffness matrices and RHS for a Laplace equation finite element equations set. \todo for regular mesh only
  SUBROUTINE LAPLACE_EQUATION_ANALYTIC_CALCULATE(FIELD,ANALYTIC_FUNCTION,GEOMETRIC_PARAMETERS,DERIVATIVE_NUMBER,NODE_NUMBER, &
    & COMPONENT_NUMBER,VARIABLE_NUMBER,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: FIELD 
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION 
    REAL(DP), INTENT(IN) :: GEOMETRIC_PARAMETERS(:)  
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER 
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER
    REAL(DP), INTENT(OUT)     :: VALUE
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS, NUMBER_OF_DERIVATIVES,NUMBER_OF_NODES
    REAL(DP) :: x,y,z
    INTEGER(INTG), ALLOCATABLE :: NODE_PARAM2DOF_MAP_X(:,:,:),NODE_PARAM2DOF_MAP_Y(:,:,:),NODE_PARAM2DOF_MAP_Z(:,:,:)
    
    CALL ENTERS("LAPLACE_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(FIELD)) THEN
      NUMBER_OF_DIMENSIONS=FIELD%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
      NUMBER_OF_NODES=FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES%NUMBER_OF_NODES
      NUMBER_OF_DERIVATIVES=FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%DOMAIN%TOPOLOGY%NODES%NODES(NODE_NUMBER)% & 
        & NUMBER_OF_DERIVATIVES
      ALLOCATE(NODE_PARAM2DOF_MAP_X(NUMBER_OF_DERIVATIVES,NUMBER_OF_NODES,FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
      ALLOCATE(NODE_PARAM2DOF_MAP_Y(NUMBER_OF_DERIVATIVES,NUMBER_OF_NODES,FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
      ALLOCATE(NODE_PARAM2DOF_MAP_Z(NUMBER_OF_DERIVATIVES,NUMBER_OF_NODES,FIELD%NUMBER_OF_VARIABLES),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old interpolation type",ERR,ERROR,*999)
      NODE_PARAM2DOF_MAP_X=FIELD%GEOMETRIC_FIELD%VARIABLES(1)%COMPONENTS(1)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP
      NODE_PARAM2DOF_MAP_Y=FIELD%GEOMETRIC_FIELD%VARIABLES(1)%COMPONENTS(2)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP
      x = GEOMETRIC_PARAMETERS(NODE_PARAM2DOF_MAP_X(1,NODE_NUMBER,1))
      y = GEOMETRIC_PARAMETERS(NODE_PARAM2DOF_MAP_Y(1,NODE_NUMBER,1))
      
      IF(NUMBER_OF_DIMENSIONS==2) THEN
        SELECT CASE(VARIABLE_NUMBER)
        CASE(FIELD_STANDARD_VARIABLE_TYPE)
	      SELECT CASE(ANALYTIC_FUNCTION)
	      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
!TODO Analytical calculation for regular mesh only, need to implement generic analytical calculation,
!i.e. du/ds1=(du/dxi1)(dxi1/ds1)=(du/dx)*(dx/dxi1)*(dxi1/ds1),du/ds2=(du/dxi2)(dxi2/ds2)=(du/dx)*(dx/dxi2)*(dxi2/ds2)	      
	        SELECT CASE(DERIVATIVE_NUMBER)
            CASE(NO_GLOBAL_DERIV)
              VALUE=x**2+2*x*y-y**2
            CASE(GLOBAL_DERIV_S1)
              VALUE=2*x+2*y
            CASE(GLOBAL_DERIV_S2)
              VALUE=2*x-2*y
            CASE(GLOBAL_DERIV_S1_S2)
              VALUE=2.0_DP
            CASE DEFAULT
              CALL FLAG_ERROR("The derivativehas not been implemented.",ERR,ERROR,*999)
            END SELECT
	      CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
!TODO Analytical calculation for regular mesh only, need to implement generic analytical calculation,
!i.e. du/ds1=(du/dxi1)(dxi1/ds1)=(du/dx)*(dx/dxi1)*(dxi1/ds1),du/ds2=(du/dxi2)(dxi2/ds2)=(du/dx)*(dx/dxi2)*(dxi2/ds2)
	        SELECT CASE(DERIVATIVE_NUMBER)
            CASE(NO_GLOBAL_DERIV)
	          VALUE=cos(x)*cosh(y)
	        CASE(GLOBAL_DERIV_S1)
	          VALUE=-sin(x)*cosh(y)
	        CASE(GLOBAL_DERIV_S2)
	          VALUE=cos(x)*sinh(y)
	        CASE(GLOBAL_DERIV_S1_S2)
              VALUE=-sin(x)*sinh(y)
	        CASE DEFAULT
              CALL FLAG_ERROR("The derivativehas not been implemented.",ERR,ERROR,*999)
	        END SELECT
	      CASE DEFAULT
            CALL FLAG_ERROR("The equation is not implemented.",ERR,ERROR,*999)
	      END SELECT
	    CASE(FIELD_NORMAL_VARIABLE_TYPE)
	      ! TODO fill the calculation
	      SELECT CASE(ANALYTIC_FUNCTION)
          CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)            
            VALUE=(2*x+2*y)*x+(2*x-2*y)*y
          CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2) 
            VALUE=-sin(x)*cosh(y)*x+cos(x)*sinh(y)*y
          CASE DEFAULT
            CALL FLAG_ERROR("The equation is not implemented.",ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          CALL FLAG_ERROR("This variable type is either not implemented or not valid.",ERR,ERROR,*999)
        END SELECT  
	  ELSE IF(NUMBER_OF_DIMENSIONS==3)THEN
	    NODE_PARAM2DOF_MAP_Z=FIELD%GEOMETRIC_FIELD%VARIABLES(1)%COMPONENTS(3)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP
	    z = GEOMETRIC_PARAMETERS(NODE_PARAM2DOF_MAP_Z(1,NODE_NUMBER,1))
	    SELECT CASE(VARIABLE_NUMBER)
        CASE(FIELD_STANDARD_VARIABLE_TYPE)
          SELECT CASE(ANALYTIC_FUNCTION)
          CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
!TODO Analytical calculation for regular mesh only, need to implement generic analytical calculation            
            SELECT CASE(DERIVATIVE_NUMBER)
            CASE(NO_GLOBAL_DERIV)
              VALUE=x**2-2*y**2+z**2
            CASE(GLOBAL_DERIV_S1)
              VALUE=2*x
            CASE(GLOBAL_DERIV_S2)
              VALUE=-4*y
            CASE(GLOBAL_DERIV_S1_S2)
              VALUE=0.0_DP
            CASE(GLOBAL_DERIV_S3)
              VALUE=2*z
            CASE(GLOBAL_DERIV_S1_S3)
              VALUE=0.0_DP
            CASE(GLOBAL_DERIV_S2_S3)
              VALUE=0.0_DP
            CASE(GLOBAL_DERIV_S1_S2_S3)
              VALUE=0.0_DP
            CASE DEFAULT
              CALL FLAG_ERROR("The derivativehas not been implemented.",ERR,ERROR,*999)
            END SELECT            
          CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2) 
!TODO Analytical calculation for regular mesh only, need to implement generic analytical calculation 
            SELECT CASE(DERIVATIVE_NUMBER)
            CASE(NO_GLOBAL_DERIV)
              VALUE=cos(x)*cosh(y)*z
            CASE(GLOBAL_DERIV_S1)
              VALUE=-sin(x)*cosh(y)*z
            CASE(GLOBAL_DERIV_S2)
              VALUE=cos(x)*sinh(y)*z
            CASE(GLOBAL_DERIV_S1_S2)
              VALUE=-sin(x)*sinh(y)*z
            CASE(GLOBAL_DERIV_S3)
              VALUE=cos(x)*cosh(y)
            CASE(GLOBAL_DERIV_S1_S3)
              VALUE=-sin(x)*cosh(y)
            CASE(GLOBAL_DERIV_S2_S3)
              VALUE=cos(x)*sinh(y) 
            CASE(GLOBAL_DERIV_S1_S2_S3)
              VALUE=-sin(x)*sinh(y)
            CASE DEFAULT
              CALL FLAG_ERROR("The derivativehas not been implemented.",ERR,ERROR,*999)
            END SELECT            
          CASE DEFAULT
            CALL FLAG_ERROR("The equation is not implemented.",ERR,ERROR,*999)
          END SELECT
	    CASE(FIELD_NORMAL_VARIABLE_TYPE)
	      !TODO fill the calculation
	      SELECT CASE(ANALYTIC_FUNCTION)
          CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)            
            VALUE=2*x**2-4*y**2+2*z**2
          CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2) 
            VALUE=-sin(x)*cosh(y)*z*x+cos(x)*sinh(y)*z*y+cos(x)*cosh(y)*z
          CASE DEFAULT
            CALL FLAG_ERROR("The equation is not implemented.",ERR,ERROR,*999)
          END SELECT
	    CASE DEFAULT
          CALL FLAG_ERROR("This variable type is either not implemented or not valid.",ERR,ERROR,*999)
        END SELECT
	  ENDIF 
    ELSE
      CALL FLAG_ERROR("The field is not associated.",ERR,ERROR,*999)
    ENDIF
      
    CALL EXITS("LAPLACE_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_ANALYTIC_CALCULATE
  
  !
  !================================================================================================================================
  !

!!MERGE: what is this
  
  !>Calculates the element stiffness matrices and RHS for a Laplace equation finite element equations set.
  SUBROUTINE LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE(EQUATIONS_SET,GEOMETRIC_PARAMETERS,NODE_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET 
    REAL(DP), INTENT(IN) :: GEOMETRIC_PARAMETERS(:)
    INTEGER(INTG), INTENT(IN) :: NODE_TYPE 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD 
    INTEGER(INTG) :: ANALYTIC_FUNCTION 
    INTEGER(INTG) :: var_idx,comp_idx,node_idx,node_number,dev_idx,NODE_TYPE_COUNT
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    REAL(DP) :: VALUE
    
    
    CALL ENTERS("LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE",ERR,ERROR,*999)
    
    FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
    ANALYTIC_FUNCTION=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION
    IF(ASSOCIATED(FIELD)) THEN
      DO var_idx=1,FIELD%NUMBER_OF_VARIABLES
        DO comp_idx=1,FIELD%VARIABLES(var_idx)%NUMBER_OF_COMPONENTS
          DOMAIN_NODES=>FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%TOPOLOGY%NODES
          IF(ASSOCIATED(DOMAIN_NODES)) THEN
	        NODES_MAPPING=>FIELD%VARIABLES(1)%COMPONENTS(comp_idx)%DOMAIN%MAPPINGS%NODES
	        SELECT CASE(NODE_TYPE)
            CASE(DOMAIN_LOCAL_INTERNAL)            
              NODE_TYPE_COUNT=NODES_MAPPING%NUMBER_OF_INTERNAL
            CASE(DOMAIN_LOCAL_BOUNDARY) 
              NODE_TYPE_COUNT=NODES_MAPPING%NUMBER_OF_BOUNDARY
            CASE DEFAULT
              CALL FLAG_ERROR("Invalid node type.",ERR,ERROR,*999)
            END SELECT
	        DO node_idx=1,NODE_TYPE_COUNT
	          SELECT CASE(NODE_TYPE)
              CASE(DOMAIN_LOCAL_INTERNAL)            
                node_number=NODES_MAPPING%INTERNAL_LIST(node_idx)
              CASE(DOMAIN_LOCAL_BOUNDARY) 
                node_number=NODES_MAPPING%BOUNDARY_LIST(node_idx)
              CASE DEFAULT
                CALL FLAG_ERROR("Invalid node type.",ERR,ERROR,*999)
              END SELECT
	          
	          DO dev_idx=1,DOMAIN_NODES%NODES(node_number)%NUMBER_OF_DERIVATIVES
	            CALL LAPLACE_EQUATION_ANALYTIC_CALCULATE(FIELD,ANALYTIC_FUNCTION,GEOMETRIC_PARAMETERS,dev_idx,node_number,comp_idx, &
	              & var_idx,VALUE,ERR,ERROR,*999)	            
	            CALL FIELD_PARAMETER_SET_UPDATE_NODE(FIELD,FIELD_ANALYTIC_SET_TYPE,dev_idx,node_number,comp_idx,var_idx,VALUE, &
	              & ERR,ERROR,*999)   
	          ENDDO ! dev_idx
            ENDDO ! node_idx
          ELSE
            CALL FLAG_ERROR("Domain nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ENDDO ! comp_idx
      ENDDO ! var_idx
    ELSE
      CALL FLAG_ERROR("The field is not associated.",ERR,ERROR,*999) 
    ENDIF
    
    CALL EXITS("LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Laplace equation finite element equations set.
  SUBROUTINE LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ng,mh,mhs,mi,ms,nh,nhs,ni,ns
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Laplace is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
          EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(1)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
          FIELD_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(1)%VARIABLE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,EQUATIONS%INTERPOLATION% &
              & GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%JACOBIAN*QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(EQUATIONS_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      ENDDO !ni
                      SUM=0.0_DP
                      DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          SUM=SUM+PGMSI(mi)*PGNSI(ni)*EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%GU(mi,ni)
                        ENDDO !ni
                      ENDDO !mi
                      EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=EQUATIONS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
          
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_SCALE_FACTORS_ELEMENT_GET(ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
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
       
        CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a Laplace equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

        
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the Laplace equation type of a classical field equations set class.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Laplace equation on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Laplace equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a Laplace equation type of a classical field equations set class.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE
        CALL LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE, &
          & EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Laplace equation type of a classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the standard Laplace equation.
  SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type to perform
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type to perform
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NEXT_NUMBER
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    
    CALL ENTERS("LAPLACE_EQUATION_EQUATION_SET_STANDARD_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE) THEN
        SELECT CASE(SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_SET%LINEARITY=EQUATIONS_SET_LINEAR
            EQUATIONS_SET%TIME_TYPE=EQUATIONS_SET_STATIC
            EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
!!TODO: maybe given negative user numbers to openCMISS generated fields???
            CALL FIELD_NEXT_NUMBER_FIND(EQUATIONS_SET%REGION,NEXT_NUMBER,ERR,ERROR,*999)
            CALL FIELD_CREATE_START(NEXT_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
            CALL FIELD_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
            CALL FIELD_DEPENDENT_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
            CALL FIELD_MESH_DECOMPOSITION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD% &
              & DECOMPOSITION,ERR,ERROR,*999)
            CALL FIELD_GEOMETRIC_FIELD_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
              & ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,1,ERR,ERROR,*999)
            !Default to the geometric interpolation setup
            CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_STANDARD_VARIABLE_TYPE,1, &
              & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)%COMPONENTS(1)% &
              & MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
            CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_NORMAL_VARIABLE_TYPE,1, &
              & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD%VARIABLES(FIELD_STANDARD_VARIABLE_TYPE)%COMPONENTS(1)% &
              & MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_STANDARD_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_NORMAL_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
              CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD% &
                & SCALINGS%SCALING_TYPE,ERR,ERROR,*999)
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
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL FIELD_CREATE_FINISH(EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,FIELD_ANALYTIC_SET_TYPE,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_ANALYTIC_SET_TYPE,ERR,ERROR,*999)
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              
              CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD%GEOMETRIC_FIELD,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
                & ERR,ERROR,*999)
              
              ! Set up boundary nodes 
              CALL LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE(EQUATIONS_SET,GEOMETRIC_PARAMETERS,DOMAIN_LOCAL_BOUNDARY, &
                & ERR,ERROR,*999)
              
	          CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_ANALYTIC_SET_TYPE,ERR,ERROR,*999)
	          
	          ! Set up internal nodes
	          CALL LAPLACE_EQUATION_ANALYTIC_PARAMETER_SET_UPDATE(EQUATIONS_SET,GEOMETRIC_PARAMETERS,DOMAIN_LOCAL_INTERNAL, &
	            & ERR,ERROR,*999)
	          
	          CALL BOUNDARY_CONDITION_PARAMETER_SET_UPDATE_FROM_ANALYTIC_VALUE(EQUATIONS_SET,ERR,ERROR,*999)
             
              CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_ANALYTIC_SET_TYPE,ERR,ERROR,*999)
              EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
            EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED=.TRUE.
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
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
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%FIXED_CONDITIONS)) THEN
              IF(EQUATIONS_SET%FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
                !Do nothing
                !?Initialise problem solution???
              ELSE
                CALL FLAG_ERROR("Equations set fixed conditions has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set fixed conditions is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              EQUATIONS=>EQUATIONS_SET%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                !Create the equations mapping.
                CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_NUMBER_SET(EQUATIONS_MAPPING,1,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_LINEAR_MATRICES_VARIABLE_TYPES_SET(EQUATIONS_MAPPING,(/FIELD_STANDARD_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_NORMAL_VARIABLE_TYPE,ERR,ERROR,*999)
                CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
                !Create the equations matrices
                CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
                IF(EQUATIONS%SPARSITY_TYPE==EQUATIONS_MATRICES_FULL_MATRICES) THEN
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_BLOCK_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                ELSE IF(EQUATIONS%SPARSITY_TYPE==EQUATIONS_MATRICES_SPARSE_MATRICES) THEN
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES,(/MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_LINEAR_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES,(/EQUATIONS_MATRIX_FEM_STRUCTURE/), &
                    & ERR,ERROR,*999)
                ELSE
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                CALL EQUATIONS_MATRICES_CREATE_FINISH(EQUATIONS_MATRICES,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
              ENDIF
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
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_EQUATIONS_SET_STANDARD_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the Laplace solution.
  SUBROUTINE LAPLACE_EQUATION_PROBLEM_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Laplace equation on.
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)
        CALL LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Laplace equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a Laplace equation type .
  SUBROUTINE LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_LAPLACE_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_STANDARD_LAPLACE_SUBTYPE     
        CALL LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP_INITIAL_TYPE,PROBLEM_SETUP_START_ACTION, &
          & ERR,ERROR,*999)
      CASE(PROBLEM_GENERALISED_LAPLACE_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a Laplace equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the standard Laplace equations solution.
  SUBROUTINE LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The setup type to perform
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The action type to perform
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(SOLUTION_TYPE), POINTER :: SOLUTION
    TYPE(SOLUTION_MAPPING_TYPE), POINTER :: SOLUTION_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_STANDARD_LAPLACE_SUBTYPE) THEN
        SELECT CASE(SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !DO nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            PROBLEM%NUMBER_OF_SOLUTIONS=1
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
!!TODO:
          CASE(PROBLEM_SETUP_FINISH_ACTION)
!!TODO:
          CASE(PROBLEM_SETUP_DO_ACTION)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a standard Laplace equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLUTION_TYPE)
          SOLUTION=>PROBLEM%SOLUTIONS(1)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SELECT CASE(ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              CALL SOLUTION_MAPPING_CREATE_START(SOLUTION,SOLUTION_MAPPING,ERR,ERROR,*999)
              CALL SOLUTION_MAPPING_SOLVER_MATRICES_NUMBER_SET(SOLUTION_MAPPING,1,ERR,ERROR,*999)
           CASE(PROBLEM_SETUP_FINISH_ACTION)
             SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING
             IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
               CALL SOLUTION_MAPPING_CREATE_FINISH(SOLUTION_MAPPING,ERR,ERROR,*999)
             ELSE
               CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
             ENDIF
           CASE(PROBLEM_SETUP_DO_ACTION)
             EQUATIONS_SET=>SOLUTION%EQUATIONS_SET_TO_ADD
             IF(ASSOCIATED(EQUATIONS_SET)) THEN
               !Check the equations set is from a standard Laplace equation
               IF(EQUATIONS_SET%CLASS==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                 & EQUATIONS_SET%TYPE==EQUATIONS_SET_LAPLACE_EQUATION_TYPE.AND. &
                 & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE) THEN
                 SOLUTION_MAPPING=>SOLUTION%SOLUTION_MAPPING            
                 IF(ASSOCIATED(SOLUTION_MAPPING)) THEN
                   CALL SOLUTION_MAPPING_EQUATIONS_SET_ADD(SOLUTION_MAPPING,SOLUTION%EQUATIONS_SET_TO_ADD,SOLUTION% &
                     & EQUATIONS_SET_ADDED_INDEX,ERR,ERROR,*999)
                   CALL SOLUTION_MAPPING_EQUATIONS_VARIABLES_TO_SOLVER_MATRIX_SET(SOLUTION_MAPPING,1,SOLUTION% &
                     & EQUATIONS_SET_ADDED_INDEX,(/FIELD_STANDARD_VARIABLE_TYPE/),ERR,ERROR,*999)
                 ELSE
                   CALL FLAG_ERROR("Solution mapping is not associated.",ERR,ERROR,*999)
                 ENDIF
               ELSE
                 CALL FLAG_ERROR("The equations set to add is not a standard Laplace equations set.",ERR,ERROR,*999)
               ENDIF
             ELSE
               CALL FLAG_ERROR("Equations set to add is not associated.",ERR,ERROR,*999)
             ENDIF
           CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Laplace equation."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem solution is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(PROBLEM_SETUP_SOLVER_TYPE)
          SOLUTION=>PROBLEM%SOLUTIONS(1)%PTR
          IF(ASSOCIATED(SOLUTION)) THEN
            SELECT CASE(ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              CALL SOLVER_CREATE_START(SOLUTION,SOLVER_LINEAR_TYPE,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
              CALL SOLVER_SPARSITY_TYPE_SET(SOLVER,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              SOLVER=>SOLUTION%SOLVER
              IF(ASSOCIATED(SOLVER)) THEN                
                CALL SOLVER_CREATE_FINISH(SOLVER,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Solution solver is not associated.",ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_SETUP_DO_ACTION)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Laplace equation."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Problem solution is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a standard Laplace equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a standard Laplace equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP")
    RETURN
999 CALL ERRORS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP",ERR,ERROR)
    CALL EXITS("LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP")
    RETURN 1
  END SUBROUTINE LAPLACE_EQUATION_PROBLEM_STANDARD_SETUP

  !
  !================================================================================================================================
  !
 
END MODULE LAPLACE_EQUATIONS_ROUTINES
