!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all equations set routines.
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

!> This module handles all equations set routines.
MODULE EQUATIONS_SET_ROUTINES

  USE BASE_ROUTINES
  USE BIOELECTRIC_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CLASSICAL_FIELD_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE ELASTICITY_ROUTINES
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE FLUID_MECHANICS_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE MPI
  USE MULTI_PHYSICS_ROUTINES
  USE NODE_ROUTINES
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EQUATIONS_SET_ANALYTIC_CREATE_START,EQUATIONS_SET_ANALYTIC_CREATE_FINISH,EQUATIONS_SET_ANALYTIC_DESTROY

  PUBLIC EQUATIONS_SET_BACKSUBSTITUTE

  PUBLIC EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC,EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH, &
    & EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START,EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY
  
  PUBLIC EQUATIONS_SET_CREATE_START,EQUATIONS_SET_CREATE_FINISH,EQUATIONS_SET_DESTROY,EQUATIONS_SETS_INITIALISE, &
    & EQUATIONS_SETS_FINALISE

  PUBLIC EQUATIONS_SET_EQUATIONS_CREATE_FINISH,EQUATIONS_SET_EQUATIONS_CREATE_START,EQUATIONS_SET_EQUATIONS_DESTROY
  
  PUBLIC EQUATIONS_SET_MATERIALS_CREATE_START,EQUATIONS_SET_MATERIALS_CREATE_FINISH,EQUATIONS_SET_MATERIALS_DESTROY

  PUBLIC EQUATIONS_SET_DEPENDENT_CREATE_START,EQUATIONS_SET_DEPENDENT_CREATE_FINISH,EQUATIONS_SET_DEPENDENT_DESTROY

  PUBLIC EQUATIONS_SET_INDEPENDENT_CREATE_START,EQUATIONS_SET_INDEPENDENT_CREATE_FINISH,EQUATIONS_SET_INDEPENDENT_DESTROY
 
  PUBLIC EQUATIONS_SET_JACOBIAN_EVALUATE,EQUATIONS_SET_RESIDUAL_EVALUATE

  PUBLIC EQUATIONS_SET_SOLUTION_METHOD_GET,EQUATIONS_SET_SOLUTION_METHOD_SET
  
  PUBLIC EQUATIONS_SET_SOURCE_CREATE_START,EQUATIONS_SET_SOURCE_CREATE_FINISH,EQUATIONS_SET_SOURCE_DESTROY

  PUBLIC EQUATIONS_SET_SPECIFICATION_GET,EQUATIONS_SET_SPECIFICATION_SET

  PUBLIC EQUATIONS_SET_ASSEMBLE

  PUBLIC EQUATIONS_SET_USER_NUMBER_FIND
  
CONTAINS

  !
  !================================================================================================================================
  !
      
  !>Finish the creation of a analytic solution for equations set. \see OPENCMISS::CMISSEquationsSetAnalyticCreateFinish
  SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create the analytic for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD

    CALL ENTERS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
          CALL FLAG_ERROR("Equations set analytic has already been finished.",ERR,ERROR,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
          IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=ANALYTIC_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>ANALYTIC_FIELD
          ENDIF
          !Finish the equations set specific analytic setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
          !Finish the analytic creation
          EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a analytic solution for a equations set. \see OPENCMISS::CMISSEquationsSetAnalyticCreateStart
  SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_START(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,ANALYTIC_FIELD_USER_NUMBER,ANALYTIC_FIELD, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of an analytic for.
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The analytic function type to setup \see EQUATIONS_SET_CONSTANTS_AnalyticFunctionTypes,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FIELD_USER_NUMBER !<The user specified analytic field number
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD !<If associated on entry, a pointer to the user created analytic field which has the same user number as the specified analytic field user number. If not associated on entry, on exit, a pointer to the created analytic field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,ANALYTIC_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_ANALYTIC_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        CALL FLAG_ERROR("The equations set analytic is already associated.",ERR,ERROR,*998)        
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
            !Check the analytic field has been finished
            IF(ANALYTIC_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(ANALYTIC_FIELD_USER_NUMBER/=ANALYTIC_FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified analytic field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified analytic field of "// &
                  & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              ANALYTIC_FIELD_REGION=>ANALYTIC_FIELD%REGION
              IF(ASSOCIATED(ANALYTIC_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(ANALYTIC_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified analytic field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified analtyic field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,ANALYTIC_FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified analtyic field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The geometric field is not associated for the specified equations set.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified analtyic field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified analtyic field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            CALL FIELD_USER_NUMBER_FIND(ANALYTIC_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified analytic field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Initialise the equations set analytic
          CALL EQUATIONS_SET_ANALYTIC_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(ANALYTIC_FIELD)) EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=ANALYTIC_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>ANALYTIC_FIELD
          EQUATIONS_SET_SETUP_INFO%ANALYTIC_FUNCTION_TYPE=ANALYTIC_FUNCTION_TYPE
          !Start the equations set specific analytic setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Set pointers
          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
            ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
          ELSE
            EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD=>ANALYTIC_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ANALYTIC_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_ANALYTIC_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ANALYTIC_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroy the analytic solution for an equations set. \see OPENCMISS::CMISSEquationsSetAnalyticDestroy
  SUBROUTINE EQUATIONS_SET_ANALYTIC_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the analytic solutins for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_ANALYTIC_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN        
        CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set analytic is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ANALYTIC_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ANALYTIC_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ANALYTIC_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the analytic solution for an equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET_ANALYTIC,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_SET_ANALYTIC !<A pointer to the equations set analytic to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_ANALYTIC_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET_ANALYTIC)) THEN        
      DEALLOCATE(EQUATIONS_SET_ANALYTIC)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ANALYTIC_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ANALYTIC_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ANALYTIC_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the analytic solution for an equations set.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("EQUATIONS_SET_ANALYTIC_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        CALL FLAG_ERROR("Analytic is already associated for this equations set.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%ANALYTIC,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set analytic.",ERR,ERROR,*999)
        EQUATIONS_SET%ANALYTIC%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED=.FALSE.
        EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ANALYTIC_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_ANALYTIC_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ANALYTIC_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_INITIALISE

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an equations set.
  SUBROUTINE EQUATIONS_SET_ASSEMBLE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%EQUATIONS_FINISHED) THEN
          SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
          CASE(EQUATIONS_STATIC)
            SELECT CASE(EQUATIONS%LINEARITY)
            CASE(EQUATIONS_LINEAR)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_NONLINEAR)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_NONLINEAR_BCS)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations linearity of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_QUASISTATIC)
! chrm, 17/09/09
            SELECT CASE(EQUATIONS%LINEARITY)
            CASE(EQUATIONS_LINEAR)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_NONLINEAR)
                CALL EQUATIONS_SET_ASSEMBLE_QUASISTATIC_NONLINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*999)
            CASE(EQUATIONS_NONLINEAR_BCS)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations linearity of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
            SELECT CASE(EQUATIONS%LINEARITY)
            CASE(EQUATIONS_LINEAR)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_NONLINEAR)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)             
            CASE(EQUATIONS_NONLINEAR_BCS)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set linearity of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_TIME_STEPPING)
            CALL FLAG_ERROR("Time stepping equations are not assembled.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations time dependence type of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ASSEMBLE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ASSEMBLE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ASSEMBLE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ASSEMBLE

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a dynamic linear equations set using the finite element method.
  SUBROUTINE EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    
    CALL ENTERS("EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ASSEMBLE_DYNAMIC_LINEAR_FEM

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a linear static equations set using the finite element method.
  SUBROUTINE EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    
!#ifdef TAUPROF
!    CHARACTER(28) :: CVAR
!    INTEGER :: PHASE(2) = (/ 0, 0 /)
!    SAVE PHASE
!#endif

    CALL ENTERS("EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Initialise the matrices and rhs vector
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_VALUES_INITIALISE()")
#endif
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_VALUES_INITIALISE()")
#endif
            !Assemble the elements
            !Allocate the element matrices 
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_ELEMENT_INITIALISE()")
#endif
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_ELEMENT_INITIALISE()")
#endif
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements

#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("Internal Elements Loop")
#endif
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
!#ifdef TAUPROF
!              WRITE (CVAR,'(a23,i3)') 'Internal Elements Loop ',element_idx
!              CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
!              CALL TAU_PHASE_START(PHASE)
!#endif
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
!#ifdef TAUPROF
!              CALL TAU_PHASE_STOP(PHASE)
!#endif
            ENDDO !element_idx
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("Internal Elements Loop")
#endif

            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("Boundary and Ghost Elements Loop")
#endif
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("Boundary and Ghost Elements Loop")
#endif
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_START("EQUATIONS_MATRICES_ELEMENT_FINALISE()")
#endif
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
#ifdef TAUPROF
            CALL TAU_STATIC_PHASE_STOP("EQUATIONS_MATRICES_ELEMENT_FINALISE()")
#endif
            !Output equations matrices and vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix, residuals and rhs for a nonlinear static equations set using the finite element method.
  SUBROUTINE EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    
    CALL ENTERS("EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
             !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM

  ! sander, 26/03/10
  !>Assembles the equations stiffness matrix, residuals and rhs for a nonlinear quasistatic equations set using the finite element method.
  !> currently the same as the static nonlinear case
  SUBROUTINE EQUATIONS_SET_ASSEMBLE_QUASISTATIC_NONLINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_NONLINEAR_FEM",ERR,ERROR,*999)
    ! currently no difference
    CALL EQUATIONS_SET_ASSEMBLE_STATIC_NONLINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*999)
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_NONLINEAR_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_NONLINEAR_FEM")
    RETURN 1
  END  SUBROUTINE EQUATIONS_SET_ASSEMBLE_QUASISTATIC_NONLINEAR_FEM


  !
  !================================================================================================================================
  !

! chrm, 17/09/09
  
  !>Assembles the equations stiffness matrix and rhs for a linear quasistatic equations set using the finite element method.
  SUBROUTINE EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    
    CALL ENTERS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ASSEMBLE_QUASISTATIC_LINEAR_FEM

  !
  !================================================================================================================================
  !

  !>Backsubstitutes with an equations set to calculate unknown right hand side vectors
  SUBROUTINE EQUATIONS_SET_BACKSUBSTITUTE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to backsubstitute
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,equations_matrix_idx,equations_row_number, &
      & EQUATIONS_STORAGE_TYPE,rhs_boundary_condition,rhs_global_dof,rhs_variable_dof,RHS_VARIABLE_TYPE,variable_dof,VARIABLE_TYPE
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(DP) :: DEPENDENT_VALUE,MATRIX_VALUE,RHS_VALUE,SOURCE_VALUE
    REAL(DP), POINTER :: DEPENDENT_PARAMETERS(:),EQUATIONS_MATRIX_DATA(:),SOURCE_VECTOR_DATA(:)
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: RHS_BOUNDARY_CONDITIONS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOMAIN_MAPPING,RHS_DOMAIN_MAPPING
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOURCE_DISTRIBUTED_VECTOR
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,RHS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_BACKSUBSTITUTE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
            IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
              DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
              IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
                !CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
              ELSE
                LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
                IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                  EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                  IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                    LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                    IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                      RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                      SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                      IF(ASSOCIATED(RHS_MAPPING)) THEN
                        BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                            SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
                            IF(ASSOCIATED(SOURCE_VECTOR)) THEN
                              SOURCE_DISTRIBUTED_VECTOR=>SOURCE_VECTOR%VECTOR
                              IF(ASSOCIATED(SOURCE_DISTRIBUTED_VECTOR)) THEN
                                CALL DISTRIBUTED_VECTOR_DATA_GET(SOURCE_DISTRIBUTED_VECTOR,SOURCE_VECTOR_DATA,ERR,ERROR,*999)
                              ELSE
                                CALL FLAG_ERROR("Source distributed vector is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("Source vector is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ENDIF
                          RHS_VARIABLE=>RHS_MAPPING%RHS_VARIABLE
                          IF(ASSOCIATED(RHS_VARIABLE)) THEN
                            RHS_VARIABLE_TYPE=RHS_VARIABLE%VARIABLE_TYPE
                            RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
                            IF(ASSOCIATED(RHS_DOMAIN_MAPPING)) THEN
                              RHS_BOUNDARY_CONDITIONS=>BOUNDARY_CONDITIONS%BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP( &
                                & RHS_VARIABLE_TYPE)%PTR
                              IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN
                                !Loop over the equations matrices
                                DO equations_matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                                  DEPENDENT_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)%VARIABLE
                                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                                    VARIABLE_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
                                    !Get the dependent field variable parameters
                                    CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                      & DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                                    EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(equations_matrix_idx)%PTR
                                    IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                                      COLUMN_DOMAIN_MAPPING=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VAR_MAPS(equations_matrix_idx)% &
                                        & COLUMN_DOFS_MAPPING
                                      IF(ASSOCIATED(COLUMN_DOMAIN_MAPPING)) THEN
                                        EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                        IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                                          CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                            & EQUATIONS_STORAGE_TYPE,ERR,ERROR,*999)
                                          CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                            & ERR,ERROR,*999)
                                          SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                          CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                            !Loop over the non ghosted rows in the equations set
                                            DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                              RHS_VALUE=0.0_DP
                                              rhs_variable_dof=RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(equations_row_number)
                                              rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                                              rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS( &
                                                & rhs_global_dof)
                                              SELECT CASE(rhs_boundary_condition)
                                              CASE(BOUNDARY_CONDITION_NOT_FIXED,BOUNDARY_CONDITION_FREE_WALL)
                                                !Back substitute
                                                !Loop over the local columns of the equations matrix
                                                DO equations_column_idx=1,COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                                  equations_column_number=COLUMN_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP( &
                                                    & equations_column_idx)
                                                  variable_dof=equations_column_idx
                                                  MATRIX_VALUE=EQUATIONS_MATRIX_DATA(equations_row_number+ &
                                                    & (equations_column_number-1)*EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)
                                                  DEPENDENT_VALUE=DEPENDENT_PARAMETERS(variable_dof)
                                                  RHS_VALUE=RHS_VALUE+MATRIX_VALUE*DEPENDENT_VALUE
                                                ENDDO !equations_column_idx
                                              CASE(BOUNDARY_CONDITION_FIXED)
                                                !Do nothing
                                              CASE(BOUNDARY_CONDITION_MIXED)
                                                !Robin or is it Cauchy??? boundary conditions
                                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                              CASE DEFAULT
                                                LOCAL_ERROR="The RHS variable boundary condition of "// &
                                                  & TRIM(NUMBER_TO_VSTRING(rhs_boundary_condition,"*",ERR,ERROR))// &
                                                  & " for RHS variable dof number "// &
                                                  & TRIM(NUMBER_TO_VSTRING(rhs_variable_dof,"*",ERR,ERROR))//" is invalid."
                                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                              END SELECT
                                              IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                                SOURCE_VALUE=SOURCE_VECTOR_DATA(equations_row_number)
                                                RHS_VALUE=RHS_VALUE-SOURCE_VALUE
                                              ENDIF
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,RHS_VARIABLE_TYPE, &
                                                & FIELD_VALUES_SET_TYPE,rhs_variable_dof,RHS_VALUE,ERR,ERROR,*999)
                                            ENDDO !equations_row_number
                                          CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                          CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                          CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                            CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX, &
                                              & ROW_INDICES,COLUMN_INDICES,ERR,ERROR,*999)
                                            !Loop over the non-ghosted rows in the equations set
                                            DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                              RHS_VALUE=0.0_DP
                                              rhs_variable_dof=RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(equations_row_number)
                                              rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                                              rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS( &
                                                & rhs_global_dof)
                                              SELECT CASE(rhs_boundary_condition)
                                              CASE(BOUNDARY_CONDITION_NOT_FIXED,BOUNDARY_CONDITION_FREE_WALL)
                                                !Back substitute
                                                !Loop over the local columns of the equations matrix
                                                DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                                  ROW_INDICES(equations_row_number+1)-1
                                                  equations_column_number=COLUMN_INDICES(equations_column_idx)
                                                  variable_dof=equations_column_idx-ROW_INDICES(equations_row_number)+1
                                                  MATRIX_VALUE=EQUATIONS_MATRIX_DATA(equations_column_idx)
                                                  DEPENDENT_VALUE=DEPENDENT_PARAMETERS(variable_dof)
                                                  RHS_VALUE=RHS_VALUE+MATRIX_VALUE*DEPENDENT_VALUE
                                                ENDDO !equations_column_idx
                                              CASE(BOUNDARY_CONDITION_FIXED)
                                                !Do nothing
                                              CASE(BOUNDARY_CONDITION_MIXED)
                                                !Robin or is it Cauchy??? boundary conditions
                                                CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                              CASE DEFAULT
                                                LOCAL_ERROR="The global boundary condition of "// &
                                                  & TRIM(NUMBER_TO_VSTRING(rhs_boundary_condition,"*",ERR,ERROR))// &
                                                  & " for RHS variable dof number "// &
                                                  & TRIM(NUMBER_TO_VSTRING(rhs_variable_dof,"*",ERR,ERROR))//" is invalid."
                                                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                              END SELECT
                                              IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                                SOURCE_VALUE=SOURCE_VECTOR_DATA(equations_row_number)
                                                RHS_VALUE=RHS_VALUE-SOURCE_VALUE
                                              ENDIF
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,RHS_VARIABLE_TYPE, &
                                                & FIELD_VALUES_SET_TYPE,rhs_variable_dof,RHS_VALUE,ERR,ERROR,*999)
                                            ENDDO !equations_row_number
                                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                        
                                          CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                          CASE DEFAULT
                                            LOCAL_ERROR="The matrix storage type of "// &
                                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_STORAGE_TYPE,"*",ERR,ERROR))//" is invalid."
                                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                          END SELECT
                                          CALL DISTRIBUTED_MATRIX_DATA_RESTORE(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                            & ERR,ERROR,*999)
                                        ELSE
                                          CALL FLAG_ERROR("Equations matrix distributed matrix is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FLAG_ERROR("Equations column domain mapping is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FLAG_ERROR("Equations equations matrix is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                    !Restore the dependent field variable parameters
                                    CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                      & DEPENDENT_PARAMETERS,ERR,ERROR,*999)
                                  ELSE
                                    CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !equations_matrix_idx
                                !Start the update of the field parameters
                                CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,RHS_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & ERR,ERROR,*999)
                                !Finish the update of the field parameters
                                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,RHS_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & ERR,ERROR,*999)
                              ELSE
                                CALL FLAG_ERROR("RHS boundary conditions variable is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FLAG_ERROR("RHS variable domain mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FLAG_ERROR("RHS variable is not associated.",ERR,ERROR,*999)
                          ENDIF
                          IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                            CALL DISTRIBUTED_VECTOR_DATA_RESTORE(SOURCE_DISTRIBUTED_VECTOR,SOURCE_VECTOR_DATA,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Equations set boundary conditions are not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Equations mapping RHS mappings is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Equations mapping linear mapping is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE            
        CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
          
    CALL EXITS("EQUATIONS_SET_BACKSUBSTITUTE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_BACKSUBSTITUTE",ERR,ERROR)    
    CALL EXITS("EQUATIONS_SET_BACKSUBSTITUTE")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_SET_BACKSUBSTITUTE
  
  !
  !================================================================================================================================
  !

  !>Set boundary conditions for an equation set according to the analytic equations. \see OPENCMISS::CMISSEquationsSetBoundaryConditionsAnalytic
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the analyticboundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
     
    CALL ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      !Initialise the setup
      CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
      EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_GENERATE_ACTION
      !Finish the equations specific solution setup.
      CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the setup
      CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC

  !
  !================================================================================================================================
  !

  !>Finish the creation of boundary conditions for an equation set. \see OPENCMISS::CMISSEquationsSetBoundaryConditionsCreateFinish
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create the boundary conditions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    
    CALL ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      !Initialise the setup
      CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE
      EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
      !Finish the equations specific solution setup.
      CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the setup
      CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of boundary conditions for an equations set. \see OPENCMISS::CMISSEquationsSetBoundaryConditionsCreateStart
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<On return, a pointer to the created boundary conditions. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
 
    CALL ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
        CALL FLAG_ERROR("Boundary conditions is already associated.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_SET%BOUNDARY_CONDITIONS)) THEN
          CALL FLAG_ERROR("The equations set boundary conditions is already associated.",ERR,ERROR,*999)        
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)



          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_BOUNDARY_CONDITIONS_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          !Start the equations set specific solution setup




          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)



          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Return the pointer
          BOUNDARY_CONDITIONS=>EQUATIONS_SET%BOUNDARY_CONDITIONS



        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the BOUNDARY conditions for an equations set and deallocate all memory. \see OPENCMISS::CMISSEquationsSetBoundaryConditionsDestroy
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%BOUNDARY_CONDITIONS)) THEN
        CALL BOUNDARY_CONDITIONS_DESTROY(EQUATIONS_SET%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set boundary conditions is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equation set on a region. \see OPENCMISS::CMISSEquationsSetCreateStart
  SUBROUTINE EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    
    CALL ENTERS("EQUATIONS_SET_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        CALL FLAG_ERROR("Equations set has already been finished.",ERR,ERROR,*999)
      ELSE            
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INITIAL_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
        !Finish the equations set specific setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_GEOMETRY_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
        !Finish the equations set specific geometry setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
        !Finish the equations set creation
        EQUATIONS_SET%EQUATIONS_SET_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("EQUATIONS_SET_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE EQUATIONS_SET_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an equations set defined by USER_NUMBER in the region identified by REGION. \see OPENCMISS::CMISSEquationsSetCreateStart
  !>Default values set for the EQUATIONS_SET's attributes are:
  !>- CLASS: 4 (EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
  !>- TYPE: 1 (EQUATIONS_SET_LAPLACE_EQUATION_TYPE)
  !>- SUBTYPE: 1 (EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
  !>- LINEARITY: 1 (EQUATIONS_SET_LINEAR)
  !>- TIME_DEPENDENCE: 1 (EQUATIONS_SET_STATIC)
  !>- SOLUTION_METHOD: 1 (EQUATIONS_SET_FEM_SOLUTION_METHOD)
  !>- GEOMETRY 
  !>- MATERIALS 
  !>- SOURCE 
  !>- DEPENDENT
  !>- ANALYTIC
  !>- FIXED_CONDITIONS 
  !>- EQUATIONS 
  SUBROUTINE EQUATIONS_SET_CREATE_START(USER_NUMBER,REGION,GEOM_FIBRE_FIELD,EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the equations set
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the equations set on
    TYPE(FIELD_TYPE), POINTER :: GEOM_FIBRE_FIELD !<A pointer to the either the geometry or, in appropriate, the fibre field for the equation set
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<On return, a pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: NEW_EQUATIONS_SET
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    NULLIFY(NEW_EQUATIONS_SET)
    NULLIFY(NEW_EQUATIONS_SETS)

    CALL ENTERS("EQUATIONS_SET_CREATE_START",ERR,ERROR,*997)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        CALL EQUATIONS_SET_USER_NUMBER_FIND(USER_NUMBER,REGION,NEW_EQUATIONS_SET,ERR,ERROR,*997)
        IF(ASSOCIATED(NEW_EQUATIONS_SET)) THEN
          LOCAL_ERROR="Equations set user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
        ELSE
          NULLIFY(NEW_EQUATIONS_SET)
          IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
            IF(GEOM_FIBRE_FIELD%FIELD_FINISHED) THEN
              IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.GEOM_FIBRE_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
                !Allocate the new equtions set
                ALLOCATE(NEW_EQUATIONS_SET,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations set.",ERR,ERROR,*999)
                !Initalise equations set
                CALL EQUATIONS_SET_INITIALISE(NEW_EQUATIONS_SET,ERR,ERROR,*999)
                !Set default equations set values
                NEW_EQUATIONS_SET%USER_NUMBER=USER_NUMBER
                NEW_EQUATIONS_SET%GLOBAL_NUMBER=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1
                NEW_EQUATIONS_SET%EQUATIONS_SETS=>REGION%EQUATIONS_SETS
                NEW_EQUATIONS_SET%REGION=>REGION
                !Default to a standardised Laplace.
                NEW_EQUATIONS_SET%CLASS=EQUATIONS_SET_CLASSICAL_FIELD_CLASS
                NEW_EQUATIONS_SET%TYPE=EQUATIONS_SET_LAPLACE_EQUATION_TYPE
                NEW_EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE
                NEW_EQUATIONS_SET%EQUATIONS_SET_FINISHED=.FALSE.
                !Initialise the setup
                CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
                EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INITIAL_TYPE
                EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
                !Start equations set specific setup
                CALL EQUATIONS_SET_SETUP(NEW_EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
                !Set up the equations set geometric fields
                CALL EQUATIONS_SET_GEOMETRY_INITIALISE(NEW_EQUATIONS_SET,ERR,ERROR,*999)
                IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
                  NEW_EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD
                  NULLIFY(NEW_EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)
                ELSE
                  NEW_EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD%GEOMETRIC_FIELD
                  NEW_EQUATIONS_SET%GEOMETRY%FIBRE_FIELD=>GEOM_FIBRE_FIELD
                ENDIF
                EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_GEOMETRY_TYPE
                EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
                EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=GEOM_FIBRE_FIELD%USER_NUMBER
                EQUATIONS_SET_SETUP_INFO%FIELD=>GEOM_FIBRE_FIELD
                !Set up equations set specific geometry
                CALL EQUATIONS_SET_SETUP(NEW_EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)                                
                !Finalise the setup
                CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
                !Add new equations set into list of equations set in the region
                ALLOCATE(NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets.",ERR,ERROR,*999)
                DO equations_set_idx=1,REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
                  NEW_EQUATIONS_SETS(equations_set_idx)%PTR=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
                ENDDO !equations_set_idx
                NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1)%PTR=>NEW_EQUATIONS_SET
                IF(ASSOCIATED(REGION%EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
                REGION%EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
                REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1
                EQUATIONS_SET=>NEW_EQUATIONS_SET
              ELSE
                CALL FLAG_ERROR("The specified geometric field is not a geometric or fibre field.",ERR,ERROR,*997)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified geometric field is not finished.",ERR,ERROR,*997)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The specified geometric field is not associated.",ERR,ERROR,*997)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The equations sets on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*997)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_EQUATIONS_SET)) CALL EQUATIONS_SET_FINALISE(NEW_EQUATIONS_SET,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
997 CALL ERRORS("EQUATIONS_SET_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_CREATE_START")
    RETURN 1   
  END SUBROUTINE EQUATIONS_SET_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroys an equations set identified by a user number on the give region and deallocates all memory. \see OPENCMISS::CMISSEquationsSetDestroy
  SUBROUTINE EQUATIONS_SET_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the equations set to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<The region of the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,equations_set_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)

    NULLIFY(NEW_EQUATIONS_SETS)

    CALL ENTERS("EQUATIONS_SET_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        
        !Find the equations set identified by the user number
        FOUND=.FALSE.
        equations_set_position=0
        DO WHILE(equations_set_position<REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS.AND..NOT.FOUND)
          equations_set_position=equations_set_position+1
          IF(REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          EQUATIONS_SET=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_position)%PTR
          
          !Destroy all the equations set components
          CALL EQUATIONS_SET_FINALISE(EQUATIONS_SET,ERR,ERROR,*999)
          
          !Remove the equations set from the list of equations set
          IF(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>1) THEN
            ALLOCATE(NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets.",ERR,ERROR,*999)
            DO equations_set_idx=1,REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
              IF(equations_set_idx<equations_set_position) THEN
                NEW_EQUATIONS_SETS(equations_set_idx)%PTR=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
              ELSE IF(equations_set_idx>equations_set_position) THEN
                REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR%GLOBAL_NUMBER=REGION%EQUATIONS_SETS% &
                  & EQUATIONS_SETS(equations_set_idx)%PTR%GLOBAL_NUMBER-1
                NEW_EQUATIONS_SETS(equations_set_idx-1)%PTR=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
              ENDIF
            ENDDO !equations_set_idx
            IF(ASSOCIATED(REGION%EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
            REGION%EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
            REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1
          ELSE
            DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
            REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Equations set number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The equations sets on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("EQUATIONS_SET_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
998 CALL ERRORS("EQUATIONS_SET_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DESTROY_NUMBER")
    RETURN 1   
  END SUBROUTINE EQUATIONS_SET_DESTROY_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Destroys an equations set identified by a pointer and deallocates all memory. \see OPENCMISS::CMISSEquationsSetDestroy
  SUBROUTINE EQUATIONS_SET_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,equations_set_position
    TYPE(EQUATIONS_SETS_TYPE), POINTER :: EQUATIONS_SETS
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)

    NULLIFY(NEW_EQUATIONS_SETS)

    CALL ENTERS("EQUATIONS_SET_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SETS=>EQUATIONS_SET%EQUATIONS_SETS
      IF(ASSOCIATED(EQUATIONS_SETS)) THEN
        equations_set_position=EQUATIONS_SET%GLOBAL_NUMBER

        !Destroy all the equations set components
        CALL EQUATIONS_SET_FINALISE(EQUATIONS_SET,ERR,ERROR,*999)
        
        !Remove the equations set from the list of equations set
        IF(EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>1) THEN
          ALLOCATE(NEW_EQUATIONS_SETS(EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets.",ERR,ERROR,*999)
          DO equations_set_idx=1,EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
            IF(equations_set_idx<equations_set_position) THEN
              NEW_EQUATIONS_SETS(equations_set_idx)%PTR=>EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
            ELSE IF(equations_set_idx>equations_set_position) THEN
              EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR%GLOBAL_NUMBER=EQUATIONS_SETS% &
                & EQUATIONS_SETS(equations_set_idx)%PTR%GLOBAL_NUMBER-1
              NEW_EQUATIONS_SETS(equations_set_idx-1)%PTR=>EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
            ENDIF
          ENDDO !equations_set_idx
          IF(ASSOCIATED(EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(EQUATIONS_SETS%EQUATIONS_SETS)
          EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
          EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1
        ELSE
          DEALLOCATE(EQUATIONS_SETS%EQUATIONS_SETS)
          EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Equations set equations set is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("EQUATIONS_SET_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
998 CALL ERRORS("EQUATIONS_SET_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DESTROY")
    RETURN 1   
  END SUBROUTINE EQUATIONS_SET_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalise the equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_FINALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      CALL EQUATIONS_SET_GEOMETRY_FINALISE(EQUATIONS_SET%GEOMETRY,ERR,ERROR,*999)
      CALL EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET%DEPENDENT,ERR,ERROR,*999)
      CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,ERR,ERROR,*999)
      CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,ERR,ERROR,*999)
      CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,ERR,ERROR,*999)
      CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) CALL EQUATIONS_DESTROY(EQUATIONS_SET%EQUATIONS,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET%BOUNDARY_CONDITIONS)) CALL BOUNDARY_CONDITIONS_DESTROY(EQUATIONS_SET%BOUNDARY_CONDITIONS, &
        & ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_SET)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FINALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a finite element equations set.
  SUBROUTINE EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(ELEMENT_VECTOR_TYPE), POINTER :: ELEMENT_VECTOR
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_DYNAMIC_TYPE), POINTER :: DYNAMIC_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%CLASS)
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
       CALL FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
        CALL BIOELECTRIC_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
       CALL MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",ERR,ERROR,*999)          
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
            DYNAMIC_MATRICES=>EQUATIONS_MATRICES%DYNAMIC_MATRICES
            IF(ASSOCIATED(DYNAMIC_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Dynamic matrices:",ERR,ERROR,*999)                        
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",DYNAMIC_MATRICES% &
                & NUMBER_OF_DYNAMIC_MATRICES,ERR,ERROR,*999)
              DO matrix_idx=1,DYNAMIC_MATRICES%NUMBER_OF_DYNAMIC_MATRICES
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR% &
                  & UPDATE_MATRIX,ERR,ERROR,*999)
                IF(DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                  ELEMENT_MATRIX=>DYNAMIC_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                    & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                    & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                    & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                    & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !matrix_idx
            ENDIF
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            IF(ASSOCIATED(LINEAR_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Linear matrices:",ERR,ERROR,*999)                        
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",LINEAR_MATRICES% &
                & NUMBER_OF_LINEAR_MATRICES,ERR,ERROR,*999)
              DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",LINEAR_MATRICES%MATRICES(matrix_idx)%PTR% &
                  & UPDATE_MATRIX,ERR,ERROR,*999)
                IF(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                  ELEMENT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                    & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                    & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                    & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                    & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !matrix_idx
            ENDIF
            RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element RHS vector :",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",RHS_VECTOR%UPDATE_VECTOR,ERR,ERROR,*999)
              IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                ELEMENT_VECTOR=>RHS_VECTOR%ELEMENT_VECTOR
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDIF
            SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
            IF(ASSOCIATED(SOURCE_VECTOR)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element source vector :",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",SOURCE_VECTOR%UPDATE_VECTOR,ERR,ERROR,*999)
              IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                ELEMENT_VECTOR=>SOURCE_VECTOR%ELEMENT_VECTOR
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equation matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF    

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE()")
#endif
       
    CALL EXITS("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian for the given element number for a finite element equations set.
  SUBROUTINE EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%CLASS)
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element Jacobian matrix:",ERR,ERROR,*999)          
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element Jacobian:",ERR,ERROR,*999)          
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update Jacobian = ",NONLINEAR_MATRICES%JACOBIAN%UPDATE_JACOBIAN, &
                & ERR,ERROR,*999)
              IF(NONLINEAR_MATRICES%JACOBIAN%UPDATE_JACOBIAN) THEN
                ELEMENT_MATRIX=>NONLINEAR_MATRICES%JACOBIAN%ELEMENT_JACOBIAN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                  & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                  & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                  & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                  & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                  & '(16X,8(X,E13.6))',ERR,ERROR,*999)
!!TODO: Write out the element residual???
              ENDIF
            ELSE
              CALL FLAG_ERROR("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equation matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF    
       
    CALL EXITS("EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vector for the given element number for a finite element equations set.
  SUBROUTINE EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(ELEMENT_VECTOR_TYPE), POINTER :: ELEMENT_VECTOR
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_NONLINEAR_TYPE), POINTER :: NONLINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%CLASS)
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element residual matrices and vectors:",ERR,ERROR,*999)          
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
            LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
            IF(ASSOCIATED(LINEAR_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Linear matrices:",ERR,ERROR,*999)                        
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",LINEAR_MATRICES% &
                & NUMBER_OF_LINEAR_MATRICES,ERR,ERROR,*999)
              DO matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",LINEAR_MATRICES%MATRICES(matrix_idx)%PTR% &
                  & UPDATE_MATRIX,ERR,ERROR,*999)
                IF(LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                  ELEMENT_MATRIX=>LINEAR_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_MATRIX%NUMBER_OF_ROWS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of columns = ",ELEMENT_MATRIX%NUMBER_OF_COLUMNS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_MATRIX%MAX_NUMBER_OF_ROWS, &
                    & ERR,ERROR,*999)
                  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",ELEMENT_MATRIX% &
                    & MAX_NUMBER_OF_COLUMNS,ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,8,8,ELEMENT_MATRIX%ROW_DOFS, &
                    & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX% &
                    & COLUMN_DOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                  CALL WRITE_STRING_MATRIX(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_MATRIX%NUMBER_OF_ROWS,1,1,ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS,8,8,ELEMENT_MATRIX%MATRIX(1:ELEMENT_MATRIX%NUMBER_OF_ROWS,1:ELEMENT_MATRIX% &
                    & NUMBER_OF_COLUMNS),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
                    & '(16X,8(X,E13.6))',ERR,ERROR,*999)
                ENDIF
              ENDDO !matrix_idx
            ENDIF
            NONLINEAR_MATRICES=>EQUATIONS_MATRICES%NONLINEAR_MATRICES
            IF(ASSOCIATED(NONLINEAR_MATRICES)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element residual vector:",ERR,ERROR,*999)                        
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",NONLINEAR_MATRICES%UPDATE_RESIDUAL,ERR,ERROR,*999)
              IF(NONLINEAR_MATRICES%UPDATE_RESIDUAL) THEN
                ELEMENT_VECTOR=>NONLINEAR_MATRICES%ELEMENT_RESIDUAL
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDIF
            RHS_VECTOR=>EQUATIONS_MATRICES%RHS_VECTOR
            IF(ASSOCIATED(RHS_VECTOR)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element RHS vector :",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",RHS_VECTOR%UPDATE_VECTOR,ERR,ERROR,*999)
              IF(RHS_VECTOR%UPDATE_VECTOR) THEN
                ELEMENT_VECTOR=>RHS_VECTOR%ELEMENT_VECTOR
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDIF
            SOURCE_VECTOR=>EQUATIONS_MATRICES%SOURCE_VECTOR
            IF(ASSOCIATED(SOURCE_VECTOR)) THEN
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Element source vector :",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update vector = ",SOURCE_VECTOR%UPDATE_VECTOR,ERR,ERROR,*999)
              IF(SOURCE_VECTOR%UPDATE_VECTOR) THEN
                ELEMENT_VECTOR=>SOURCE_VECTOR%ELEMENT_VECTOR
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Number of rows = ",ELEMENT_VECTOR%NUMBER_OF_ROWS,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",ELEMENT_VECTOR%MAX_NUMBER_OF_ROWS, &
                  & ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%ROW_DOFS, &
                  & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',ERR,ERROR,*999)
                CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,ELEMENT_VECTOR%NUMBER_OF_ROWS,8,8,ELEMENT_VECTOR%VECTOR, &
                  & '("  Vector(:)    :",8(X,E13.6))','(16X,8(X,E13.6))',ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equation matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF    
       
    CALL EXITS("EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE")
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Finish the creation of independent variables for an equations set. \see OPENCMISS::CMISSEquationsSetIndependentCreateFinish
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the independent field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD

    CALL ENTERS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FINISHED) THEN
          CALL FLAG_ERROR("Equations set independent field has already been finished.",ERR,ERROR,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INDEPENDENT_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
          IF(ASSOCIATED(INDEPENDENT_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=INDEPENDENT_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>INDEPENDENT_FIELD
            !Finish equations set specific startup
            CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set independent independent field is not associated.",ERR,ERROR,*999)
          ENDIF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
          !Finish independent creation
          EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set independent is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of independent variables for an equations set. \see OPENCMISS::CMISSEquationsSetIndependentCreateStart
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_START(EQUATIONS_SET,INDEPENDENT_FIELD_USER_NUMBER,INDEPENDENT_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(IN) :: INDEPENDENT_FIELD_USER_NUMBER !<The user specified independent field number
    TYPE(FIELD_TYPE), POINTER :: INDEPENDENT_FIELD !<If associated on entry, a pointer to the user created independent field which has the same user number as the specified independent field user number. If not associated on entry, on exit, a pointer to the created independent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,INDEPENDENT_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_INDEPENDENT_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        CALL FLAG_ERROR("The equations set independent is already associated",ERR,ERROR,*998)        
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(INDEPENDENT_FIELD)) THEN
            !Check the independent field has been finished
            IF(INDEPENDENT_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(INDEPENDENT_FIELD_USER_NUMBER/=INDEPENDENT_FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified independent field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(INDEPENDENT_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified independent field of "// &
                  & TRIM(NUMBER_TO_VSTRING(INDEPENDENT_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              INDEPENDENT_FIELD_REGION=>INDEPENDENT_FIELD%REGION
              IF(ASSOCIATED(INDEPENDENT_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(INDEPENDENT_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified independent field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(INDEPENDENT_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified independent field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,INDEPENDENT_FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified independent field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The geometric field is not associated for the specified equations set.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified independent field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified independent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            CALL FIELD_USER_NUMBER_FIND(INDEPENDENT_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified independent field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(INDEPENDENT_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Initialise the equations set independent
          CALL EQUATIONS_SET_INDEPENDENT_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(INDEPENDENT_FIELD)) EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INDEPENDENT_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=INDEPENDENT_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>INDEPENDENT_FIELD
          !Start equations set specific startup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Set pointers
          IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN            
            INDEPENDENT_FIELD=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
          ELSE
            EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD=>INDEPENDENT_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equation set region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_INDEPENDENT_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the independent field for an equations set. \see OPENCMISS::CMISSEquationsSetIndependentDestroy
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the independent field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_INDEPENDENT_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set indpendent is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_INDEPENDENT_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the independent field for an equations set.
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET_INDEPENDENT,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_INDEPENDENT_TYPE), POINTER :: EQUATIONS_SET_INDEPENDENT !<A pointer to the equations set independent to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_INDEPENDENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET_INDEPENDENT)) THEN
      DEALLOCATE(EQUATIONS_SET_INDEPENDENT)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_INDEPENDENT_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the independent field for an equations set.
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the independent for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_SET_INDEPENDENT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        CALL FLAG_ERROR("Independent field is already associated for these equations sets.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%INDEPENDENT,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set independent field.",ERR,ERROR,*999)
        EQUATIONS_SET%INDEPENDENT%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FINISHED=.FALSE.
        EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_INDEPENDENT_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_INDEPENDENT_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations set.
  SUBROUTINE EQUATIONS_SET_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<The pointer to the equations set to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SET%USER_NUMBER=0
      EQUATIONS_SET%GLOBAL_NUMBER=0
      EQUATIONS_SET%EQUATIONS_SET_FINISHED=.FALSE.
      NULLIFY(EQUATIONS_SET%EQUATIONS_SETS)
      NULLIFY(EQUATIONS_SET%REGION)
      EQUATIONS_SET%CLASS=EQUATIONS_SET_NO_CLASS
      EQUATIONS_SET%TYPE=EQUATIONS_SET_NO_TYPE
      EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SUBTYPE
      EQUATIONS_SET%SOLUTION_METHOD=0
      CALL EQUATIONS_SET_GEOMETRY_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
      CALL EQUATIONS_SET_DEPENDENT_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
      NULLIFY(EQUATIONS_SET%INDEPENDENT)
      NULLIFY(EQUATIONS_SET%MATERIALS)
      NULLIFY(EQUATIONS_SET%SOURCE)
      NULLIFY(EQUATIONS_SET%ANALYTIC)
      NULLIFY(EQUATIONS_SET%EQUATIONS)
      NULLIFY(EQUATIONS_SET%BOUNDARY_CONDITIONS)
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise the geometry for an equations set
  SUBROUTINE EQUATIONS_SET_GEOMETRY_FINALISE(EQUATIONS_SET_GEOMETRY,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_GEOMETRY_TYPE) :: EQUATIONS_SET_GEOMETRY !<A pointer to the equations set geometry to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_GEOMETRY_FINALISE",ERR,ERROR,*999)
    
    NULLIFY(EQUATIONS_SET_GEOMETRY%GEOMETRIC_FIELD)
    NULLIFY(EQUATIONS_SET_GEOMETRY%FIBRE_FIELD)
       
    CALL EXITS("EQUATIONS_SET_GEOMETRY_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_GEOMETRY_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_GEOMETRY_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_GEOMETRY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the geometry for an equation set
  SUBROUTINE EQUATIONS_SET_GEOMETRY_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the geometry for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_GEOMETRY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SET%GEOMETRY%EQUATIONS_SET=>EQUATIONS_SET
      NULLIFY(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)
      NULLIFY(EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_GEOMETRY_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_GEOMETRY_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_GEOMETRY_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_GEOMETRY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of materials for an equations set. \see OPENCMISS::CMISSEquationsSetMaterialsCreateFinish
  SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the materials field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD

    CALL ENTERS("EQUATIONS_SET_MATERIALS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Equations set materials has already been finished.",ERR,ERROR,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_MATERIALS_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          IF(ASSOCIATED(MATERIALS_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=MATERIALS_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>MATERIALS_FIELD
            !Finish equations set specific startup
            CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set materials materials field is not associated.",ERR,ERROR,*999)
          ENDIF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
          !Finish materials creation
          EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of materials for a problem. \see OPENCMISS::CMISSEquationsSetMaterialsCreateStart
  SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,MATERIALS_FIELD_USER_NUMBER,MATERIALS_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(IN) :: MATERIALS_FIELD_USER_NUMBER !<The user specified materials field number
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD !<If associated on entry, a pointer to the user created materials field which has the same user number as the specified materials field user number. If not associated on entry, on exit, a pointer to the created materials field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,MATERIALS_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_MATERIALS_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL FLAG_ERROR("The equations set materials is already associated",ERR,ERROR,*998)        
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(MATERIALS_FIELD)) THEN
            !Check the materials field has been finished
            IF(MATERIALS_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(MATERIALS_FIELD_USER_NUMBER/=MATERIALS_FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified materials field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(MATERIALS_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified materials field of "// &
                  & TRIM(NUMBER_TO_VSTRING(MATERIALS_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              MATERIALS_FIELD_REGION=>MATERIALS_FIELD%REGION
              IF(ASSOCIATED(MATERIALS_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(MATERIALS_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified materials field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(MATERIALS_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified materials field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,MATERIALS_FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified materials field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The geometric field is not associated for the specified equations set.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified materials field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified materials field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            CALL FIELD_USER_NUMBER_FIND(MATERIALS_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified materials field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(MATERIALS_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Initialise the equations set materials
          CALL EQUATIONS_SET_MATERIALS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(MATERIALS_FIELD)) EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_MATERIALS_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=MATERIALS_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>MATERIALS_FIELD
          !Start equations set specific startup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Set pointers
          IF(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN            
            MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          ELSE
            EQUATIONS_SET%MATERIALS%MATERIALS_FIELD=>MATERIALS_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equation set region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_MATERIALS_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the materials for an equations set. \see OPENCMISS::CMISSEquationsSetMaterialsDestroy
  SUBROUTINE EQUATIONS_SET_MATERIALS_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the materials for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_MATERIALS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the materials for an equations set.
  SUBROUTINE EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET_MATERIALS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_SET_MATERIALS !<A pointer to the equations set materials to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_MATERIALS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET_MATERIALS)) THEN
      DEALLOCATE(EQUATIONS_SET_MATERIALS)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the materials for an equations set.
  SUBROUTINE EQUATIONS_SET_MATERIALS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the materials for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_SET_MATERIALS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL FLAG_ERROR("Materials is already associated for these equations sets.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%MATERIALS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set materials.",ERR,ERROR,*999)
        EQUATIONS_SET%MATERIALS%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED=.FALSE.
        EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_MATERIALS_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_INITIALISE

  !
  !
  !================================================================================================================================
  !

  !>Finish the creation of a dependent variables for an equations set. \see OPENCMISS::CMISSEquationsSetDependentCreateFinish
  SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)
    
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD

    CALL ENTERS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Equations set dependent has already been finished",ERR,ERROR,*999)
      ELSE
        !Initialise the setup
        CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_DEPENDENT_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=DEPENDENT_FIELD%USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>DEPENDENT_FIELD
          !Finish equations set specific setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Equations set dependent dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
        !Finish the equations set creation
        EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of dependent variables for an equations set. \see OPENCMISS::CMISSEquationsSetDependentCreateStart
  SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,DEPENDENT_FIELD_USER_NUMBER,DEPENDENT_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a dependent field on
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_FIELD_USER_NUMBER !<The user specified dependent field number
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<If associated on entry, a pointer to the user created dependent field which has the same user number as the specified dependent field user number. If not associated on entry, on exit, a pointer to the created dependent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,DEPENDENT_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_SET_DEPENDENT_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("The equations set dependent has been finished.",ERR,ERROR,*999)
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
            !Check the dependent field has been finished
            IF(DEPENDENT_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(DEPENDENT_FIELD_USER_NUMBER/=DEPENDENT_FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified dependent field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified dependent field of "// &
                  & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              DEPENDENT_FIELD_REGION=>DEPENDENT_FIELD%REGION
              IF(ASSOCIATED(DEPENDENT_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(DEPENDENT_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified dependent field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified dependent field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,DEPENDENT_FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified dependent field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The geometric field is not associated for the specified equations set.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified dependent field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified dependent field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            CALL FIELD_USER_NUMBER_FIND(DEPENDENT_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified dependent field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(DEPENDENT_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & " has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
            EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.TRUE.
          ENDIF
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_DEPENDENT_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=DEPENDENT_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>DEPENDENT_FIELD
          !Start the equations set specfic solution setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Set pointers
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          ELSE
            EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD=>DEPENDENT_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equation set region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations_set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_DEPENDENT_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET%DEPENDENT,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_DEPENDENT_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Destroy the dependent variables for an equations set. \see OPENCMISS::CMISSEquationsSetDependentDestroy
  SUBROUTINE EQUATIONS_SET_DEPENDENT_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<The pointer to the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_DEPENDENT_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      CALL EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET%DEPENDENT,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_DEPENDENT_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_DESTROY
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET_DEPENDENT,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_DEPENDENT_TYPE) :: EQUATIONS_SET_DEPENDENT !<The pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_DEPENDENT_FINALISE",ERR,ERROR,*999)

    NULLIFY(EQUATIONS_SET_DEPENDENT%EQUATIONS_SET)
    EQUATIONS_SET_DEPENDENT%DEPENDENT_FINISHED=.FALSE.
    EQUATIONS_SET_DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.FALSE.
    NULLIFY(EQUATIONS_SET_DEPENDENT%DEPENDENT_FIELD)
    
    CALL EXITS("EQUATIONS_SET_DEPENDENT_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the dependent variables for a equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_DEPENDENT_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SET%DEPENDENT%EQUATIONS_SET=>EQUATIONS_SET
      EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED=.FALSE.
      EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.FALSE.
      NULLIFY(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_DEPENDENT_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets up the specifices for an equation set.
  SUBROUTINE EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the setup on
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP_INFO !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%CLASS)
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
        CALL BIOELECTRIC_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SETUP")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SETUP",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SETUP")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

 !>Finish the creation of equations for the equations set. \see OPENCMISS::CMISSEquationsSetEquationsCreateFinish
  SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      !Initialise the setup
      CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_EQUATIONS_TYPE
      EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
      !Finish the equations specific solution setup.
      CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
      !Finalise the setup
      CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of equations for the equation set. \see CMISSEquationsSetEquationsCreateStart
  !>Default values set for the EQUATIONS's attributes are:
  !>- OUTPUT_TYPE: 0 (EQUATIONS_SET_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (EQUATIONS_SET_SPARSE_MATRICES)
  !>- NONLINEAR_JACOBIAN_TYPE: 0
  !>- INTERPOLATION: null
  !>- LINEAR_DATA: null 
  !>- NONLINEAR_DATA: null
  !>- TIME_DATA: null
  !>- EQUATIONS_MAPPING:  
  !>- EQUATIONS_MATRICES:  
  SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create equations for
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS)) THEN
        CALL FLAG_ERROR("Equations is already associated.",ERR,ERROR,*999)
      ELSE
        !Initialise the setup
        CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_EQUATIONS_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
        !Start the equations set specific solution setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        !Return the pointer
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATION_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the equations for an equations set. \see OPENCMISS::CMISSEquationsSetEquationsDestroy
  SUBROUTINE EQUATIONS_SET_EQUATIONS_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) THEN        
        CALL EQUATIONS_FINALISE(EQUATIONS_SET%EQUATIONS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear equations set.
  SUBROUTINE EQUATIONS_SET_JACOBIAN_EVALUATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_JACOBIAN_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%EQUATIONS_FINISHED) THEN
          SELECT CASE(EQUATIONS%LINEARITY)
          CASE(EQUATIONS_LINEAR)            
            CALL FLAG_ERROR("Can not evaluate a Jacobian for linear equations.",ERR,ERROR,*999)
          CASE(EQUATIONS_NONLINEAR)
            SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
            CASE(EQUATIONS_STATIC)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_QUASISTATIC)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
! sebk 15/09/09
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_TIME_STEPPING)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations linearity of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_JACOBIAN_EVALUATE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_JACOBIAN_EVALUATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_JACOBIAN_EVALUATE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_JACOBIAN_EVALUATE

 !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an static equations set using the finite element method
  SUBROUTINE EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
  
    CALL ENTERS("EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
!!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
             !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
             !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF            
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_JACOBIAN_EVALUATE_STATIC_FEM

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an dynamic equations set using the finite element method
  SUBROUTINE EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
  
    CALL ENTERS("EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
!!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
             !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
             !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF            
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_JACOBIAN_EVALUATE_DYNAMIC_FEM

  !
  !================================================================================================================================
  !


  !>Evaluates the residual for an equations set.
  SUBROUTINE EQUATIONS_SET_RESIDUAL_EVALUATE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_RESIDUAL_EVALUATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%EQUATIONS_FINISHED) THEN
          SELECT CASE(EQUATIONS%LINEARITY)
          CASE(EQUATIONS_LINEAR)            
            CALL FLAG_ERROR("Can not evaluate a residual for linear equations.",ERR,ERROR,*999)
          CASE(EQUATIONS_NONLINEAR)
            SELECT CASE(EQUATIONS%TIME_DEPENDENCE)
            CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC) ! quasistatic handled like static
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
! sebk 19/08/09
!|
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM(EQUATIONS_SET,ERR,ERROR,*999)
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
                LOCAL_ERROR="The equations set solution method  of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
!|
! sebk 19/08/09
            CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_TIME_STEPPING)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set time dependence type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS%TIME_DEPENDENCE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations linearity of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS%LINEARITY,"*",ERR,ERROR))//" is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_RESIDUAL_EVALUATE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_RESIDUAL_EVALUATE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_RESIDUAL_EVALUATE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_RESIDUAL_EVALUATE

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an dynamic equations set using the finite element method
  SUBROUTINE EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
 
    CALL ENTERS("EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
             !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_RESIDUAL_EVALUATE_DYNAMIC_FEM

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an static equations set using the finite element method
  SUBROUTINE EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
 
    CALL ENTERS("EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !!Do we need to transfer parameter sets???
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for equations setup and initialisation = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for equations setup and initialisation = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              ELEMENT_USER_ELAPSED=0.0_SP
              ELEMENT_SYSTEM_ELAPSED=0.0_SP
            ENDIF
            NUMBER_OF_TIMES=0
            !Loop over the internal elements
            DO element_idx=ELEMENTS_MAPPING%INTERNAL_START,ELEMENTS_MAPPING%INTERNAL_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
              SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
              ELEMENT_USER_ELAPSED=USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
             ENDIF
             !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME4,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME4,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME4(1)-USER_TIME3(1)
              SYSTEM_ELAPSED=SYSTEM_TIME4(1)-SYSTEM_TIME3(1)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for parameter transfer completion = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for parameter transfer completion = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)              
            ENDIF
            !Loop over the boundary and ghost elements
            DO element_idx=ELEMENTS_MAPPING%BOUNDARY_START,ELEMENTS_MAPPING%GHOST_FINISH
              ne=ELEMENTS_MAPPING%DOMAIN_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME5,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME5,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME5(1)-USER_TIME4(1)
              SYSTEM_ELAPSED=SYSTEM_TIME5(1)-SYSTEM_TIME4(1)
              ELEMENT_USER_ELAPSED=ELEMENT_USER_ELAPSED+USER_ELAPSED
              ELEMENT_SYSTEM_ELAPSED=ELEMENT_SYSTEM_ELAPSED+USER_ELAPSED
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for boundary+ghost equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for boundary+ghost equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              IF(NUMBER_OF_TIMES>0) THEN
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element user time for equations assembly = ", &
                  & ELEMENT_USER_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Average element system time for equations assembly = ", &
                  & ELEMENT_SYSTEM_ELAPSED/NUMBER_OF_TIMES,ERR,ERROR,*999)
              ENDIF
            ENDIF
            !Finalise the element matrices
            CALL EQUATIONS_MATRICES_ELEMENT_FINALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME6,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME6,ERR,ERROR,*999)
              USER_ELAPSED=USER_TIME6(1)-USER_TIME1(1)
              SYSTEM_ELAPSED=SYSTEM_TIME6(1)-SYSTEM_TIME1(1)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total user time for equations assembly = ",USER_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Total system time for equations assembly = ",SYSTEM_ELAPSED, &
                & ERR,ERROR,*999)
              CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"***",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations matrices is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Dependent field is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_RESIDUAL_EVALUATE_STATIC_FEM

  !
  !================================================================================================================================
  !

  !>Finalises the equations set setup and deallocates all memory
  SUBROUTINE EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(OUT) :: EQUATIONS_SET_SETUP_INFO !<The equations set setup to be finalised
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_SETUP_FINALISE",ERR,ERROR,*999)

    EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=0
    EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=0
    EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=0
    NULLIFY(EQUATIONS_SET_SETUP_INFO%FIELD)
    EQUATIONS_SET_SETUP_INFO%ANALYTIC_FUNCTION_TYPE=0
    
    CALL EXITS("EQUATIONS_SET_SETUP_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SETUP_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SETUP_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SETUP_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialise the equations set setup.
  SUBROUTINE EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(OUT) :: EQUATIONS_SET_SETUP_INFO !<The equations set setup to be initialised
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_SETUP_INITIALISE",ERR,ERROR,*999)

    EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=0
    EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=0
    EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=0
    NULLIFY(EQUATIONS_SET_SETUP_INFO%FIELD)
    EQUATIONS_SET_SETUP_INFO%ANALYTIC_FUNCTION_TYPE=0
    
    CALL EXITS("EQUATIONS_SET_SETUP_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SETUP_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SETUP_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SETUP_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for an equations set. \see OPENCMISS::CMISSEquationsSetSolutionMethodSet
  SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The equations set solution method to set \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        CALL FLAG_ERROR("Equations set has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(EQUATIONS_SET%CLASS)
        CASE(EQUATIONS_SET_ELASTICITY_CLASS)
          CALL ELASTICITY_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
          CALL BIOELECTRIC_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SOLUTION_METHOD_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the solution method for an equations set. \see OPENCMISS::CMISSEquationsSetSolutionMethodGet
  SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_GET(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to get the solution method for
    INTEGER(INTG), INTENT(OUT) :: SOLUTION_METHOD !<On return, the equations set solution method \see EQUATIONS_SET_CONSTANTS_SolutionMethods,EQUATIONS_SET_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_SOLUTION_METHOD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        SOLUTION_METHOD=EQUATIONS_SET%SOLUTION_METHOD
      ELSE
        CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_SOLUTION_METHOD_GET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SOLUTION_METHOD_GET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOLUTION_METHOD_GET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_GET
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of a source for an equation set. \see OPENCMISS::CMISSEquationsSetSourceCreateFinish
  SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a souce for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD

    CALL ENTERS("EQUATIONS_SET_SOURCE_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        IF(EQUATIONS_SET%SOURCE%SOURCE_FINISHED) THEN
          CALL FLAG_ERROR("Equations set source has already been finished.",ERR,ERROR,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_SOURCE_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
          IF(ASSOCIATED(SOURCE_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=SOURCE_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>SOURCE_FIELD
            !Finish the equation set specific source setup
            CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Equations set source source field is not associated.",ERR,ERROR,*999)
          ENDIF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
          !Finish the source creation
          EQUATIONS_SET%SOURCE%SOURCE_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set source is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SOURCE_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SOURCE_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOURCE_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a source for an equations set. \see OPENCMISS::CMISSEquationsSetSourceCreateStart
  SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_START(EQUATIONS_SET,SOURCE_FIELD_USER_NUMBER,SOURCE_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a source for
    INTEGER(INTG), INTENT(IN) :: SOURCE_FIELD_USER_NUMBER !<The user specified source field number
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD !<If associated on entry, a pointer to the user created source field which has the same user number as the specified source field user number. If not associated on entry, on exit, a pointer to the created source field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,SOURCE_FIELD_REGION
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_SOURCE_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL FLAG_ERROR("The equations set source is already associated.",ERR,ERROR,*998)        
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(SOURCE_FIELD)) THEN
            !Check the source field has been finished
            IF(SOURCE_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(SOURCE_FIELD_USER_NUMBER/=SOURCE_FIELD%USER_NUMBER) THEN
                LOCAL_ERROR="The specified source field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOURCE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not match the user number of the specified source field of "// &
                  & TRIM(NUMBER_TO_VSTRING(SOURCE_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
              SOURCE_FIELD_REGION=>SOURCE_FIELD%REGION
              IF(ASSOCIATED(SOURCE_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(SOURCE_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  LOCAL_ERROR="Invalid region setup. The specified source field has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(SOURCE_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                !Check the specified source field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,SOURCE_FIELD%DECOMPOSITION)) THEN
                    CALL FLAG_ERROR("The specified source field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The geometric field is not associated for the specified equations set.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified source field region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified source field has not been finished.",ERR,ERROR,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            CALL FIELD_USER_NUMBER_FIND(SOURCE_FIELD_USER_NUMBER,REGION,FIELD,ERR,ERROR,*999)
            IF(ASSOCIATED(FIELD)) THEN
              LOCAL_ERROR="The specified source field user number of "// &
                & TRIM(NUMBER_TO_VSTRING(SOURCE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
          !Initialise the equations set source
          CALL EQUATIONS_SET_SOURCE_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(SOURCE_FIELD)) EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_SOURCE_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=SOURCE_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>SOURCE_FIELD
          !Start the equation set specific source setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
          !Set pointers
          IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN            
            SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
          ELSE
            EQUATIONS_SET%SOURCE%SOURCE_FIELD=>SOURCE_FIELD
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equation set region is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SOURCE_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_SOURCE_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOURCE_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the source for an equations set. \see OPENCMISS::CMISSEquationsSetSourceDestroy
  SUBROUTINE EQUATIONS_SET_SOURCE_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the source for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_SOURCE_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SOURCE_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SOURCE_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOURCE_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the source for a equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET_SOURCE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SET_SOURCE !<A pointer to the equations set source to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_SOURCE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET_SOURCE)) THEN
       DEALLOCATE(EQUATIONS_SET_SOURCE)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SOURCE_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SOURCE_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOURCE_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the source for an equations set.
  SUBROUTINE EQUATIONS_SET_SOURCE_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the source field for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    
    CALL ENTERS("EQUATIONS_SET_SOURCE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL FLAG_ERROR("Source is already associated for this equations set.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%SOURCE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations set source.",ERR,ERROR,*999)
        EQUATIONS_SET%SOURCE%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%SOURCE%SOURCE_FINISHED=.FALSE.
        EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%SOURCE%SOURCE_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SOURCE_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_SOURCE_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOURCE_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the equations set specification i.e., equations set class, type and subtype for an equations set. \see OPENCMISS::CMISSEquationsSetSpecificationGet
  SUBROUTINE EQUATIONS_SET_SPECIFICATION_GET(EQUATIONS_SET,EQUATIONS_SET_CLASS,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to get the specification for
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_CLASS !<On return, the equations set class.
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_TYPE_ !<On return, the equations set type.
    INTEGER(INTG), INTENT(OUT) :: EQUATIONS_SET_SUBTYPE !<On return, the equations set subtype.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_SPECIFICATION_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        EQUATIONS_SET_CLASS=EQUATIONS_SET%CLASS
        EQUATIONS_SET_TYPE_=EQUATIONS_SET%TYPE
        EQUATIONS_SET_SUBTYPE=EQUATIONS_SET%SUBTYPE
      ELSE
        CALL FLAG_ERROR("Equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_GET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SPECIFICATION_GET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_GET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SPECIFICATION_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the equations set specification i.e., equations set class, type and subtype for an equations set. \see OPENCMISS::CMISSEquationsSetSpecificationSet
  SUBROUTINE EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_CLASS,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_CLASS !<The equations set class to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_TYPE_ !<The equations set type to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equations set subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_SPECIFICATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        CALL FLAG_ERROR("Equations set has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(EQUATIONS_SET_CLASS)
        CASE(EQUATIONS_SET_ELASTICITY_CLASS)
          CALL ELASTICITY_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
          CALL FLUID_MECHANICS_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
          CALL BIOELECTRIC_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
          CALL MULTI_PHYSICS_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_CLASS,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
        !Initialise the setup
        CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INITIAL_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
        !Peform the initial equations set setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,ERR,ERROR,*999)          
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SPECIFICATION_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SPECIFICATION_SET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in EQUATIONS_SET a pointer to the equations set identified by USER_NUMBER in the given REGION. If no equations set with that USER_NUMBER exists EQUATIONS_SET is left nullified.
  SUBROUTINE EQUATIONS_SET_USER_NUMBER_FIND(USER_NUMBER,REGION,EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find the equation set
    TYPE(REGION_TYPE), POINTER :: REGION !<The region to find the equations set in
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<On return, a pointer to the equations set if an equations set with the specified user number exists in the given region. If no equation set with the specified number exists a NULL pointer is returned. The pointer must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL FLAG_ERROR("Equations set is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(EQUATIONS_SET)
        IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
          equations_set_idx=1
          DO WHILE(equations_set_idx<=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS.AND..NOT.ASSOCIATED(EQUATIONS_SET))
            IF(REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              EQUATIONS_SET=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
            ELSE
              equations_set_idx=equations_set_idx+1
            ENDIF
          ENDDO
        ELSE
          LOCAL_ERROR="The equations sets on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises all equations sets on a region and deallocates all memory.
  SUBROUTINE EQUATIONS_SETS_FINALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise the problems for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SETS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        DO WHILE(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>0)
          CALL EQUATIONS_SET_DESTROY(REGION%EQUATIONS_SETS%EQUATIONS_SETS(1)%PTR,ERR,ERROR,*999)
        ENDDO !problem_idx
        DEALLOCATE(REGION%EQUATIONS_SETS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EQUATIONS_SETS_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SETS_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SETS_FINALISE")
    RETURN 1   
  END SUBROUTINE EQUATIONS_SETS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all equations sets on a region.
  SUBROUTINE EQUATIONS_SETS_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the equations sets for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        CALL FLAG_ERROR("Region already has associated equations sets",ERR,ERROR,*998)
      ELSE
!!TODO: Inherit any equations sets from the parent region???
        ALLOCATE(REGION%EQUATIONS_SETS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region equations sets",ERR,ERROR,*999)
        REGION%EQUATIONS_SETS%REGION=>REGION
        REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=0
        NULLIFY(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("EQUATIONS_SETS_INITIALISE")
    RETURN
999 IF(ASSOCIATED(REGION%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS)
998 CALL ERRORS("EQUATIONS_SETS_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SETS_INITIALISE")
    RETURN 1   
  END SUBROUTINE EQUATIONS_SETS_INITIALISE
  
  !
  !================================================================================================================================
  !
  
END MODULE EQUATIONS_SET_ROUTINES
