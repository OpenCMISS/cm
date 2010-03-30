!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all diffusion equation routines.
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

!>This module handles all diffusion equation routines.
MODULE DIFFUSION_EQUATION_ROUTINES

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
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MATRIX_VECTOR
  USE PROBLEM_CONSTANTS
  USE STRINGS
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  USE FLUID_MECHANICS_IO_ROUTINES

  IMPLICIT NONE

  PRIVATE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS
  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC DIFFUSION_EQUATION_EQUATIONS_SET_SETUP,DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET, &
    & DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET,DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE, &
    & DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET,DIFFUSION_EQUATION_PROBLEM_SETUP, &
    & DIFFUSION_EQUATION_PRE_SOLVE, DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE, &
    & DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION
  
CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !Calculates a two-dimensional unsteady heat equation solution (diffusion coefficient is 1)
  SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type
    REAL(DP) :: VALUE,X(3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    !TYPE(VARYING_STRING) :: LOCAL_ERROR    
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    !THESE ARE TEMPORARY VARIABLES - they need to be replace by constant field values and the current simulation time
    REAL(DP) :: CURRENT_TIME
    CURRENT_TIME = 1.0_DP

    CALL ENTERS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & ERR,ERROR,*999)            
             NULLIFY(BOUNDARY_CONDITIONS)
             CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
             DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                    IF(ASSOCIATED(DOMAIN)) THEN
                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                          !Loop over the local nodes excluding the ghosts.
                          DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                            !!TODO \todo We should interpolate the geometric field here and the node position.
                            DO dim_idx=1,NUMBER_OF_DIMENSIONS
                              local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                              X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                            ENDDO !dim_idx
                            !Loop over the derivatives
                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                              ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                              GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx)
                              CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,CURRENT_TIME,variable_type, & 
                                & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,ERR,ERROR,*999)
                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                  !If we are a boundary node then set the analytic value on the boundary
                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDDO !deriv_idx
                          ENDDO !node_idx
                        ELSE
                          CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
                  ENDIF
                ENDDO !component_idx
                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                  & ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
              ENDIF

             ENDDO !variable_idx
             CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
             CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
          ENDIF            
        ELSE
          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_CALCULATE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_CALCULATE


  !
  !================================================================================================================================
  !
  !>Calculates the various analytic solutions given X and time, can be called from within analytic calculate or elsewhere if needed
  SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X,CURRENT_TIME,VARIABLE_TYPE, & 
    & GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    REAL(DP), INTENT(OUT) :: VALUE
    REAL(DP), INTENT(IN) :: CURRENT_TIME
    REAL(DP), INTENT(IN), DIMENSION(3) :: X
    !Local variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP) :: k,phi,A,B,C,D
    INTEGER(INTG) :: variable_type,GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    !TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    !TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES

    !These are parameters for the analytical solution
    k = 1.0_DP !this is a time decay constant for the exponential term
    phi = 0.785398163397_DP !pi/4 - this sets the orientation of the solution relative to the axes
    !These are parameters for the 3D analytical solution with a linear source
    A = -0.25_DP
    B = 0.5_DP   
    C = 0.5_DP
    D = 0.5_DP   

    CALL ENTERS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS",ERR,ERROR,*999)

     SELECT CASE(ANALYTIC_FUNCTION_TYPE)
     CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
       !u=exp(-kt)*sin(sqrt(k)*(x*cos(phi)+y*sin(phi)))
       SELECT CASE(variable_type)
       CASE(FIELD_U_VARIABLE_TYPE)
         SELECT CASE(GLOBAL_DERIV_INDEX)
         CASE(NO_GLOBAL_DERIV)
           VALUE=EXP(-k*CURRENT_TIME)*SIN((SQRT(k))*(X(1)*COS(phi)+X(2)*SIN(phi)))!Need to specify time, k and phi!
         CASE(GLOBAL_DERIV_S1)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(GLOBAL_DERIV_S2)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(GLOBAL_DERIV_S1_S2)
           CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
         CASE DEFAULT
           LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
            GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
            & " is invalid."
           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
         END SELECT
       CASE(FIELD_DELUDELN_VARIABLE_TYPE)
         SELECT CASE(GLOBAL_DERIV_INDEX)
         CASE(NO_GLOBAL_DERIV)
           VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
         CASE(GLOBAL_DERIV_S1)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(GLOBAL_DERIV_S2)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
         CASE(GLOBAL_DERIV_S1_S2)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE DEFAULT
           LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
            GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
            & " is invalid."
           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
         END SELECT
       CASE DEFAULT
         LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
          & " is invalid."
         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
       END SELECT
     CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
       !u=exp(-kt)*sin(sqrt(k)*(x*cos(phi)+y*sin(phi)))
       SELECT CASE(variable_type)
       CASE(FIELD_U_VARIABLE_TYPE)
         SELECT CASE(GLOBAL_DERIV_INDEX)
         CASE(NO_GLOBAL_DERIV)
           VALUE=EXP(A*CURRENT_TIME)*EXP(B*X(1))*EXP(C*X(2))*EXP(D*X(3))!Need to specify time, k and phi!
         CASE(GLOBAL_DERIV_S1)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(GLOBAL_DERIV_S2)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(GLOBAL_DERIV_S1_S2)
           CALL FLAG_ERROR("Not implmented.",ERR,ERROR,*999)
         CASE DEFAULT
           LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
            GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
            & " is invalid."
           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
         END SELECT
       CASE(FIELD_DELUDELN_VARIABLE_TYPE)
         SELECT CASE(GLOBAL_DERIV_INDEX)
         CASE(NO_GLOBAL_DERIV)
           VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
         CASE(GLOBAL_DERIV_S1)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE(GLOBAL_DERIV_S2)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                                    
         CASE(GLOBAL_DERIV_S1_S2)
           CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
         CASE DEFAULT
           LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
            GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
            & " is invalid."
           CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
         END SELECT
       CASE DEFAULT
         LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
          & " is invalid."
         CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
       END SELECT                                                             
     CASE DEFAULT
       LOCAL_ERROR="The analytic function type of "// &
        & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
        & " is invalid."
       CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
     END SELECT

    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS
  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a diffusion equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
         CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field equation set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)        
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
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)        
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
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)        
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
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)        
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
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)        
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
          & " is not valid for a diffusion equation type of an classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation subtype for a diffusion equation type of a classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET(EQUATIONS_SET,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the equation subtype for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET_SUBTYPE)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
        EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
        EQUATIONS_SET%TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TYPE
        EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field equations set class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the linear diffusion equation.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS, &
      & NUMBER_OF_MATERIALS_COMPONENTS, NUMBER_OF_SOURCE_COMPONENTS
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SOURCE
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_EQUATION_SET_LINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(EQUATIONS_MAPPING)
    NULLIFY(EQUATIONS_MATRICES)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL DIFFUSION_EQUATION_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,ERR,ERROR,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",ERR,ERROR,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,ERR,ERROR,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE/),ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & ERR,ERROR,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
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
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                & ERR,ERROR,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                & ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
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
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation"
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
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",ERR,ERROR,*999)
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
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ELSE
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
                ENDIF
                 !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx
                !Default the source materials components to the first component geometric interpolation with constant interpolation
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !components_idx
                ENDIF
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
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                   & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & ERR,ERROR,*999)
                ELSE
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                    & ERR,ERROR,*999)
                ENDIF
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
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS             
                ELSEIF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Now set the linear source values to 1.0
                  WRITE(*,'("Setting components for linear source :")')
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !This has been pasted from the materials field section - need to modify so that it is a source field
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Create the auto created source field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SOURCE% &
                  & SOURCE_FIELD,ERR,ERROR,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,"Source Field",ERR,ERROR,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,(/FIELD_U_VARIABLE_TYPE/), &
                  & ERR,ERROR,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                NUMBER_OF_SOURCE_COMPONENTS=1
                !Set the number of source components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_SOURCE_COMPONENTS,ERR,ERROR,*999)
                !Default the source components to the geometric interpolation setup with constant interpolation
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                  DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  ENDDO !component_idx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE/),ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Finish creating the source field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
                !Set the default values for the source field
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                  NUMBER_OF_SOURCE_COMPONENTS=1
                ELSE
                  NUMBER_OF_SOURCE_COMPONENTS=0
                ENDIF
                !Now set the source values to 1.0
                DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
                    !Check there is no source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE) THEN
                      !Check that domain is 2D
                      IF(NUMBER_OF_DIMENSIONS/=2) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 2 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analtyic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an no source diffusion equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
                    !Check there is no source
                    IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                      !Check that domain is 3D
                      IF(NUMBER_OF_DIMENSIONS/=3) THEN
                        LOCAL_ERROR="The number of geometric dimensions of "// &
                          & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                          & " is invalid. The analytic function type of "// &
                          & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                          & " requires that there be 3 geometric dimensions."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                      !Create analytic field if required
                      !Set analtyic function type
                      EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                        & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1
                    ELSE
                      LOCAL_ERROR="The equations set subtype of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                        & " requires that the equations set subtype be an no source diffusion equation."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
                      & " is invalid for a linear diffusion equation."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
                  CALL DIFFUSION_EQUATION_ANALYTIC_CALCULATE(EQUATIONS_SET,ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Equations set analytic has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations set dependent has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_LINEARITY_TYPE_SET(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
              CALL EQUATIONS_TIME_DEPENDENCE_TYPE_SET(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations
              CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL EQUATIONS_CREATE_FINISH(EQUATIONS,ERR,ERROR,*999)
              !Create the equations mapping.
              CALL EQUATIONS_MAPPING_CREATE_START(EQUATIONS,EQUATIONS_MAPPING,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_MATRICES_SET(EQUATIONS_MAPPING,.TRUE.,.TRUE.,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_DYNAMIC_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EQUATIONS_MAPPING_RHS_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              IF(EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                   & EQUATIONS_SET%SUBTYPE==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN              
                CALL EQUATIONS_MAPPING_SOURCE_VARIABLE_TYPE_SET(EQUATIONS_MAPPING,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              ENDIF
              CALL EQUATIONS_MAPPING_CREATE_FINISH(EQUATIONS_MAPPING,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EQUATIONS_MATRICES_CREATE_START(EQUATIONS,EQUATIONS_MATRICES,ERR,ERROR,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%LUMPING_TYPE==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EQUATIONS_MATRICES_DYNAMIC_LUMPING_TYPE_SET(EQUATIONS_MATRICES, &
                  & (/EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED/),ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                  & (/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE/),ERR,ERROR,*999)
                CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                  (/EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE/),ERR,ERROR,*999)
              ELSE
                SELECT CASE(EQUATIONS%SPARSITY_TYPE)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EQUATIONS_MATRICES_LINEAR_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & (/DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE/),ERR,ERROR,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STORAGE_TYPE_SET(EQUATIONS_MATRICES, &
                    & (/DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE/), &
                    & ERR,ERROR,*999)
                  CALL EQUATIONS_MATRICES_DYNAMIC_STRUCTURE_TYPE_SET(EQUATIONS_MATRICES, &
                    (/EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE/),ERR,ERROR,*999)                  
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(EQUATIONS%SPARSITY_TYPE,"*",ERR,ERROR))//" is invalid."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END SELECT
              ENDIF
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EQUATIONS_SET_EQUATIONS_GET(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
            IF(EQUATIONS%EQUATIONS_FINISHED) THEN
              CALL BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Equations set equations has not been finished.",ERR,ERROR,*999)               
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_GET(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
            CALL BOUNDARY_CONDITIONS_CREATE_FINISH(BOUNDARY_CONDITIONS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
          & " is not a linear diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_LINEAR_SETUP

  !
  !================================================================================================================================
  !
  
!   SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)
!     !Argument variables
!     TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
!     TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
!     INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
! 
!     CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!     CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
!     RETURN
! 999 CALL ERRORS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP",ERR,ERROR)
!     CALL EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
!     RETURN 1
! 
!   END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a diffusion equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM%SUBTYPE)
      CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
!        CALL DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the diffusion problem pre-solve.
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local variables
    LOGICAL :: UPDATE_MATERIALS
    LOGICAL :: UPDATE_BOUNDARY_CONDITIONS

    UPDATE_MATERIALS = .FALSE.
    UPDATE_BOUNDARY_CONDITIONS = .TRUE.

    IF( UPDATE_MATERIALS ) THEN
      !CALL DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_MATERIALS_FIELD(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
    ENDIF

    IF( UPDATE_BOUNDARY_CONDITIONS ) THEN
      !CALL DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !
  !>Within the diffusion pre-solve, update the boundary conditions
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
!    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
!    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
!    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
!    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
!    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS!,BOUNDARY_CONDITION_CHECK_VARIABLE

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP) :: VALUE,X(3) !<The value to add
!     REAL(DP) :: k_xx, k_yy, k_zz
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,variable_idx
    INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX
!    INTEGER(INTG) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
!    INTEGER(INTG) :: DERIVATIVE_NUMBER !<The node derivative number
!    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number
!    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of (geometry) nodes
!    INTEGER(INTG) :: LOCAL_NODE_NUMBER
!    INTEGER(INTG) :: EQUATIONS_SET_IDX
!    INTEGER(INTG) :: equations_row_number

    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,ERR,ERROR,*999)
       !write(*,*)'CURRENT_TIME = ',CURRENT_TIME
       !write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                            NULLIFY(GEOMETRIC_VARIABLE)
                            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
                              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                             DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
                                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                    IF(ASSOCIATED(DOMAIN)) THEN
                                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                          !Loop over the local nodes excluding the ghosts.
                                          DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                            !!TODO \todo We should interpolate the geometric field here and the node position.
                                            DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                              local_ny= & 
                                          & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                                              X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                            ENDDO !dim_idx
                                            !Loop over the derivatives
                                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                              ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                              GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx)
                                              CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, & 
                                                & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
                                                & ANALYTIC_FUNCTION_TYPE,ERR,ERROR,*999)
                                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
!                                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
!                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                                  !If we are a boundary node then set the analytic value on the boundary
!                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
!                                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
!                                                ENDIF
!                                              ENDIF
                                            ENDDO !deriv_idx
                                          ENDDO !node_idx
                                        ELSE
                                          CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !component_idx
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                              ELSE
                                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                              ENDIF

                             ENDDO !variable_idx
                             CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
                              & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                          ENDIF            
                        ELSE
                          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                  END IF                
                ELSE
                  CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
                END IF  
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !do nothing?! 
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
            !do nothing?! 
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
                              & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                            NULLIFY(GEOMETRIC_VARIABLE)
                            CALL FIELD_VARIABLE_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,ERR,ERROR,*999)
                            CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
                              & GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                             DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                              FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
                                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                    IF(ASSOCIATED(DOMAIN)) THEN
                                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                        DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                          !Loop over the local nodes excluding the ghosts.
                                          DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                            !!TODO \todo We should interpolate the geometric field here and the node position.
                                            DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                              local_ny= & 
                                          & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
                                              X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                            ENDDO !dim_idx
                                            !Loop over the derivatives
                                            DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                              ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                              GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(deriv_idx)
                                              CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, & 
                                                & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
                                                & ANALYTIC_FUNCTION_TYPE,ERR,ERROR,*999)
                                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,ERR,ERROR,*999)
!                                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
!                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                                  !If we are a boundary node then set the analytic value on the boundary
!                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
!                                                    & BOUNDARY_CONDITION_FIXED,VALUE,ERR,ERROR,*999)
!                                                ENDIF
!                                              ENDIF
                                            ENDDO !deriv_idx
                                          ENDDO !node_idx
                                        ELSE
                                          CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FLAG_ERROR("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !component_idx
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type, &
                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                              ELSE
                                CALL FLAG_ERROR("Field variable is not associated.",ERR,ERROR,*999)
                              ENDIF

                             ENDDO !variable_idx
                             CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
                              & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,ERR,ERROR,*999)
                          ELSE
                            CALL FLAG_ERROR("Equations set geometric field is not associated.",ERR,ERROR,*999)
                          ENDIF            
                        ELSE
                          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations are not associated.",ERR,ERROR,*999)
                  END IF                
                ELSE
                  CALL FLAG_ERROR("Solver equations are not associated.",ERR,ERROR,*999)
                END IF  
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
              ! do nothing ???
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS
  !   
  !================================================================================================================================
  !
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION_ONE !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION_ONE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION_ONE !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DIFFUSION_ONE !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_DIFFUSION_ONE !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    REAL(DP), POINTER :: DUMMY_VALUES2(:)

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
    INTEGER(INTG) :: NDOFS_TO_PRINT, I

    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_DIFFUSION_ONE)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion-one equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Store value diffusion-one dependent field at time, t ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_DIFFUSION_ONE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      DEPENDENT_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_ONE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION_ONE, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DIFFUSION_ONE is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion-one equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion-one solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion-one solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,ERR,ERROR,*999)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIELDMESHDISPLACEMENTTYPE needs to be changed to appropriate type for this problem
                  END DO


!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',ERR,ERROR,*999)
!                   CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                 ENDIF

              ELSE  
                ! do nothing ???
!                 CALL FLAG_ERROR("DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT may only be carried out for SOLVER%GLOBAL_NUMBER = 2", &
!                   & ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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



    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLUTION    
  !   
  !================================================================================================================================
  !
  SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION_ONE, SOLVER_DIFFUSION_TWO  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION_TWO, SOURCE_FIELD_DIFFUSION_ONE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION_ONE, SOLVER_EQUATIONS_DIFFUSION_TWO  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DIFFUSION_ONE, SOLVER_MAPPING_DIFFUSION_TWO !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_DIFFUSION_ONE, EQUATIONS_SET_DIFFUSION_TWO !<A pointer to the equations set
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    REAL(DP), POINTER :: DUMMY_VALUES2(:)

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE
    INTEGER(INTG) :: NDOFS_TO_PRINT, I


    CALL ENTERS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_DIFFUSION_ONE)
      NULLIFY(SOLVER_DIFFUSION_TWO)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          SELECT CASE(CONTROL_LOOP%PROBLEM%SUBTYPE)
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion_two equations
                CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Update diffusion-one source field ... ",ERR,ERROR,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION_TWO,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_TWO=>SOLVER_DIFFUSION_TWO%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_TWO)) THEN
                  SOLVER_MAPPING_DIFFUSION_TWO=>SOLVER_EQUATIONS_DIFFUSION_TWO%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_TWO)) THEN
                    EQUATIONS_SET_DIFFUSION_TWO=>SOLVER_MAPPING_DIFFUSION_TWO%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_TWO)) THEN
                      DEPENDENT_FIELD_DIFFUSION_TWO=>EQUATIONS_SET_DIFFUSION_TWO%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_TWO)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(DEPENDENT_FIELD_DIFFUSION_TWO, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("DEPENDENT_FIELD_DIFFUSION_TWO is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion-two equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion-two solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion-two solver equations are not associated.",ERR,ERROR,*999)
                END IF


                !--- Get the source field for the diffusion_one equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_DIFFUSION_ONE,ERR,ERROR,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%PTR
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      SOURCE_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%SOURCE%SOURCE_FIELD
                      IF(ASSOCIATED(SOURCE_FIELD_DIFFUSION_ONE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(SOURCE_FIELD_DIFFUSION_ONE, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE,ERR,ERROR,*999)
                      ELSE
                        CALL FLAG_ERROR("SOURCE_FIELD_DIFFUSION_ONE is not associated.",ERR,ERROR,*999)
                      END IF
                    ELSE
                      CALL FLAG_ERROR("Diffusion-one equations set is not associated.",ERR,ERROR,*999)
                    END IF
                  ELSE
                    CALL FLAG_ERROR("Diffusion-one solver mapping is not associated.",ERR,ERROR,*999)
                  END IF
                ELSE
                  CALL FLAG_ERROR("Diffusion-one solver equations are not associated.",ERR,ERROR,*999)
                END IF

                !--- Copy the result from diffusion-two's dependent field to diffusion-one's source field
                IF(NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE
                    CALL FIELD_PARAMETERS_TO_FIELD_PARAMETERS_COMPONENT_COPY(DEPENDENT_FIELD_DIFFUSION_TWO, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,ERR,ERROR,*999)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIELDMESHDISPLACEMENTTYPE needs to be changed to appropriate type for this problem
                  END DO
                ELSE
                  LOCAL_ERROR="Number of components of diffusion-two dependent field "// &
                    & "is not consistent with diffusion-one-equation source field."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                END IF

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',ERR,ERROR,*999)
!                   CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,ERR,ERROR,*999)
!                 ENDIF

              ELSE  
                ! do nothing ???
!                 CALL FLAG_ERROR("DARCY_EQUATION_PRE_SOLVE_GET_SOLID_DISPLACEMENT may only be carried out for SOLVER%GLOBAL_NUMBER = 2", &
!                   & ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
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



    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE
  !   
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) mh,mhs,ms,ng,nh,nhs,ni,nj,ns
    REAL(DP) :: C_PARAM,K_PARAM,RWG,SUM,PGMJ(3),PGNJ(3),A_PARAM
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_DYNAMIC_TYPE), POINTER :: DYNAMIC_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: DAMPING_MATRIX,STIFFNESS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SUBTYPE)
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
          !Loop over gauss points
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
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
          
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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
          !MOST OF THIS CODE WILL BE SAME AS WITHOUT A SOURCE-NEED TO GET THE SOURCE CONTRIBUTION INTO THE RHS_VECTOR
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD !should this be source field, or independent field?
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
          !Loop over gauss points
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
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF              
              ENDDO !ms
            ENDDO !mh
             IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                !C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+1, &
                 ! & NO_PART_DERIV)
                 C_PARAM=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT%VALUES(1, NO_PART_DERIV)
                 mhs=0
                 DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                 !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                   DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                       & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                   ENDDO !ms
                 ENDDO !mh
             ENDIF
             IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP 
          ENDDO !ng
          
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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)
                IF(SOURCE_VECTOR%UPDATE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
          DEPENDENT_FIELD=>EQUATIONS%INTERPOLATION%DEPENDENT_FIELD
          GEOMETRIC_FIELD=>EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD
          MATERIALS_FIELD=>EQUATIONS%INTERPOLATION%MATERIALS_FIELD
          SOURCE_FIELD=>EQUATIONS%INTERPOLATION%SOURCE_FIELD !should this be source field, or independent field?
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
          STIFFNESS_MATRIX=>DYNAMIC_MATRICES%MATRICES(1)%PTR
          DAMPING_MATRIX=>DYNAMIC_MATRICES%MATRICES(2)%PTR
          RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
          SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
          EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
          DYNAMIC_MAPPING=>EQUATIONS_MAPPING%DYNAMIC_MAPPING
          FIELD_VARIABLE=>DYNAMIC_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(1)%VARIABLE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%PTR
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,EQUATIONS%INTERPOLATION% &
            & SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
          !Loop over gauss points
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
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        A_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS, &
                          & NO_PART_DERIV)
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+SUM*RWG- &
                          & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF              
              ENDDO !ms
            ENDDO !mh
             IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                !C_PARAM=EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT%VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+1, &
                 ! & NO_PART_DERIV)
                 C_PARAM=EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT%VALUES(1, NO_PART_DERIV)
                 mhs=0
                 DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                 !DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                   DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)+ &
                       & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                   ENDDO !ms
                 ENDDO !mh
             ENDIF
             IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=0.0_DP 
          ENDDO !ng
          
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
                IF(STIFFNESS_MATRIX%UPDATE_MATRIX.OR.DAMPING_MATRIX%UPDATE_MATRIX) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(STIFFNESS_MATRIX%UPDATE_MATRIX) THEN
                        STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=STIFFNESS_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(DAMPING_MATRIX%UPDATE_MATRIX) THEN
                        DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)=DAMPING_MATRIX%ELEMENT_MATRIX%MATRIX(mhs,nhs)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)* &
                          & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(RHS_VECTOR%UPDATE_VECTOR) RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=RHS_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)
                IF(SOURCE_VECTOR%UPDATE_VECTOR) SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)=SOURCE_VECTOR%ELEMENT_VECTOR%VECTOR(mhs)* &
                  & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not calculate finite element stiffness matrices for a nonlinear source.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
          CALL FLAG_ERROR("Can not calculate finite element stiffness matrices for a nonlinear source.",ERR,ERROR,*999)  
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SUBTYPE,"*",ERR,ERROR))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !INSERT CODE HERE TO DEAL WITH THE NON-LINEAR SOURCE TERMS - DECIDE HOW TO SOLVE THEM, USE JACOBIAN AS FOR NON-LINEAR POISSON EXAMPLE?

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a diffusion equation type.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET(PROBLEM,PROBLEM_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to set the problem subtype for
    INTEGER(INTG), INTENT(IN) :: PROBLEM_SUBTYPE !<The problem subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SUBTYPE)
      CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE  
      CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE
      CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)        
        PROBLEM%CLASS=PROBLEM_CLASSICAL_FIELD_CLASS
        PROBLEM%TYPE=PROBLEM_DIFFUSION_EQUATION_TYPE
        PROBLEM%SUBTYPE=PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SUBTYPE,"*",ERR,ERROR))// &
          & " is not valid for a diffusion equation type of a classical field problem class."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SUBTYPE_SET

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(PROBLEM%SUBTYPE==PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM%SUBTYPE==PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
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
            !Set the solver to be a first order dynamic solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICHOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a linear diffusion equation."
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
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
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
              & " is invalid for a linear diffusion equation."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a linear diffusion equation."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SUBTYPE,"*",ERR,ERROR))// &
          & " does not equal a linear diffusion equation subtype."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN
999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",ERR,ERROR)
    CALL EXITS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP
  
  !
  !================================================================================================================================
  !
    !>Sets up the diffusion equations.
!   SUBROUTINE DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
! 
!     !Argument variables
!     TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
!     TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
!     INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
!     CALL EXITS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
!     RETURN
! 999 CALL ERRORS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",ERR,ERROR)
!     CALL EXITS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
!     RETURN 1
! 
!   END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP

END MODULE DIFFUSION_EQUATION_ROUTINES
