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

!> This module handles all equations set routines.
MODULE EQUATIONS_SET_ROUTINES

  USE BASE_ROUTINES
  USE CLASSICAL_FIELD_ROUTINES
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE ELASTICITY_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_MATRICES_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MATRIX_VECTOR
  USE MPI
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

  INTERFACE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF
    MODULE PROCEDURE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS
    MODULE PROCEDURE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1
  END INTERFACE !EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF
  
  INTERFACE EQUATIONS_SET_SPECIFICATION_SET
    MODULE PROCEDURE EQUATIONS_SET_SPECIFICATION_SET_NUMBER
    MODULE PROCEDURE EQUATIONS_SET_SPECIFICATION_SET_PTR
  END INTERFACE !EQUATIONS_SET_SPECIFICATION_SET

  PUBLIC EQUATIONS_SET_BACKSUBSTITUTE
  
  PUBLIC EQUATIONS_SET_CREATE_START,EQUATIONS_SET_CREATE_FINISH,EQUATIONS_SET_DESTROY,EQUATIONS_SETS_INITIALISE, &
    & EQUATIONS_SETS_FINALISE

  PUBLIC EQUATIONS_SET_EQUATIONS_CREATE_START,EQUATIONS_SET_EQUATIONS_CREATE_FINISH, &
    & EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET,EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET

  PUBLIC EQUATIONS_SET_FIXED_CONDITIONS_APPLY,EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START, &
    & EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH,EQUATIONS_SET_FIXED_CONDITIONS_DESTROY, &
    & EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF,EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE

  PUBLIC EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET,EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET, &
    & EQUATIONS_SET_MATERIALS_CREATE_START,EQUATIONS_SET_MATERIALS_CREATE_FINISH,EQUATIONS_SET_MATERIALS_DESTROY, &
    & EQUATIONS_SET_MATERIALS_SCALING_SET

  PUBLIC EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET,EQUATIONS_SET_DEPENDENT_CREATE_START, &
    & EQUATIONS_SET_DEPENDENT_CREATE_FINISH,EQUATIONS_SET_DEPENDENT_DESTROY,EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET, &
    & EQUATIONS_SET_DEPENDENT_SCALING_SET
 
  PUBLIC EQUATIONS_SET_ANALYTIC_CREATE_START,EQUATIONS_SET_ANALYTIC_CREATE_FINISH,EQUATIONS_SET_ANALYITIC_FUNCTION_SET

  PUBLIC EQUATIONS_SET_JACOBIAN_EVALUATE,EQUATIONS_SET_RESIDUAL_EVALUATE
  
  PUBLIC EQUATIONS_SET_SOURCE_CREATE_START,EQUATIONS_SET_SOURCE_CREATE_FINISH,EQUATIONS_SET_SOURCE_DESTROY, &
    & EQUATIONS_SET_SOURCE_SCALING_SET

  PUBLIC EQUATIONS_SET_ASSEMBLE
  
  PUBLIC EQUATIONS_SET_SPECIFICATION_SET
  
CONTAINS

  !
  !================================================================================================================================
  !
      
  !>Finish the creation of a analytic solution for equations set.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create the analytic for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
          CALL FLAG_ERROR("Equations set analytic has already been finished",ERR,ERROR,*999)
        ELSE
          !Finish the equations set specific analytic setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_ANALYTIC_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION, &
            & ERR,ERROR,*999)
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

  !>Start the creation of a analytic solution for a equations set.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_START(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of an analytic for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_SET_ANALYTIC_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        CALL FLAG_ERROR("The equations set analytic is already associated.",ERR,ERROR,*998)        
      ELSE
        !Initialise the equations set analytic
        CALL EQUATIONS_SET_ANALYTIC_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
        !Start the equations set specific analytic setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_ANALYTIC_TYPE,EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
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
  
  !>Sets/changes the analytic equation number.
  SUBROUTINE EQUATIONS_SET_ANALYITIC_FUNCTION_SET(EQUATIONS_SET,ANALYTIC_FUNCTION,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the number
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION !<The analytic function to set 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_ANALYITIC_FUNCTION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
        CALL FLAG_ERROR("Equations set analytic has been finished",ERR,ERROR,*999)
      ELSE
        EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION=ANALYTIC_FUNCTION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_ANALYITIC_FUNCTION_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_ANALYITIC_FUNCTION_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_ANALYITIC_FUNCTION_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYITIC_FUNCTION_SET

  !
  !================================================================================================================================
  !

  !>Destroy the analytic solution for an equations set.
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
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
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
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
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
          SELECT CASE(EQUATIONS_SET%TIME_TYPE)
          CASE(EQUATIONS_SET_STATIC)
            SELECT CASE(EQUATIONS_SET%LINEARITY)
            CASE(EQUATIONS_SET_LINEAR)
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
            CASE(EQUATIONS_SET_NONLINEAR)
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
            CASE(EQUATIONS_SET_NONLINEAR_BCS)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set linearity of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%LINEARITY,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_DYNAMIC)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_QUASISTATIC)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations set time type of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TIME_TYPE,"*",ERR,ERROR))//" is invalid."
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
    
    CALL ENTERS("EQUATIONS_SET_ASSEMBLE_STATIC_LINEAR_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        EQUATIONS=>EQUATIONS_SET%EQUATIONS
        IF(ASSOCIATED(EQUATIONS)) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Start the transfer of the solution values that have been set as part of the boundary conditions
            CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Problem interpolation setup
            CALL EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS,ERR,ERROR,*999)
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL
              ne=ELEMENTS_MAPPING%INTERNAL_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finish the transfer of the solution values.
            !CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
!!TODO: sort out combining boundary and ghost list
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY
              ne=ELEMENTS_MAPPING%BOUNDARY_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_GHOST
              ne=ELEMENTS_MAPPING%GHOST_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx          
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finalise the problem interpolation
            CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !Start the transfer of the solution values that have been set as part of the boundary conditions
            CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Problem interpolation setup
            CALL EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS,ERR,ERROR,*999)
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL
              ne=ELEMENTS_MAPPING%INTERNAL_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finish the transfer of the solution values.
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
!!TODO: sort out combining boundary and ghost list
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY
              ne=ELEMENTS_MAPPING%BOUNDARY_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_GHOST
              ne=ELEMENTS_MAPPING%GHOST_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx          
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finalise the problem interpolation
            CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
      & EQUATIONS_STORAGE_TYPE,field_dof,rhs_boundary_condition,rhs_field_dof,rhs_variable_dof,variable_dof
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(DP) :: DEPENDENT_VALUE,MATRIX_VALUE,RHS_VALUE,SOURCE_VALUE
    REAL(DP), POINTER :: DEPENDENT_PARAMETERS(:),EQUATIONS_MATRIX_DATA(:),SOURCE_VECTOR_DATA(:)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOMAIN_MAPPING
    TYPE(DISTRIBUTED_MATRIX_TYPE), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX
    TYPE(DISTRIBUTED_VECTOR_TYPE), POINTER :: SOURCE_DISTRIBUTED_VECTOR
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_MAPPING_TYPE), POINTER :: EQUATIONS_MAPPING
    TYPE(EQUATIONS_MAPPING_LINEAR_TYPE), POINTER :: LINEAR_MAPPING
    TYPE(EQUATIONS_MAPPING_RHS_TYPE), POINTER :: RHS_MAPPING
    TYPE(EQUATIONS_MAPPING_SOURCE_TYPE), POINTER :: SOURCE_MAPPING
    TYPE(EQUATIONS_MATRICES_TYPE), POINTER :: EQUATIONS_MATRICES
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(EQUATIONS_MATRIX_TYPE), POINTER :: EQUATIONS_MATRIX
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,RHS_VARIABLE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_BACKSUBSTITUTE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          !Get the dependent field parameters
          CALL FIELD_PARAMETER_SET_GET(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
          EQUATIONS=>EQUATIONS_SET%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
            IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
              LINEAR_MATRICES=>EQUATIONS_MATRICES%LINEAR_MATRICES
              IF(ASSOCIATED(LINEAR_MATRICES)) THEN
                EQUATIONS_MAPPING=>EQUATIONS%EQUATIONS_MAPPING
                IF(ASSOCIATED(EQUATIONS_MAPPING)) THEN
                  LINEAR_MAPPING=>EQUATIONS_MAPPING%LINEAR_MAPPING
                  IF(ASSOCIATED(LINEAR_MAPPING)) THEN
                    RHS_MAPPING=>EQUATIONS_MAPPING%RHS_MAPPING
                    SOURCE_MAPPING=>EQUATIONS_MAPPING%SOURCE_MAPPING
                    IF(ASSOCIATED(RHS_MAPPING)) THEN
                      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
                      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
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
                          !Loop over the equations matrices
                          DO equations_matrix_idx=1,LINEAR_MATRICES%NUMBER_OF_LINEAR_MATRICES
                            DEPENDENT_VARIABLE=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(equations_matrix_idx)%VARIABLE
                            IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                              EQUATIONS_MATRIX=>LINEAR_MATRICES%MATRICES(equations_matrix_idx)%PTR
                              IF(ASSOCIATED(EQUATIONS_MATRIX)) THEN
                                COLUMN_DOMAIN_MAPPING=>LINEAR_MAPPING%EQUATIONS_MATRIX_TO_VARIABLE_MAPS(equations_matrix_idx)% &
                                  & COLUMN_DOFS_MAPPING
                                IF(ASSOCIATED(COLUMN_DOMAIN_MAPPING)) THEN
                                  EQUATIONS_DISTRIBUTED_MATRIX=>EQUATIONS_MATRIX%MATRIX
                                  IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                                    CALL DISTRIBUTED_MATRIX_STORAGE_TYPE_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_STORAGE_TYPE, &
                                      & ERR,ERROR,*999)
                                    CALL DISTRIBUTED_MATRIX_DATA_GET(EQUATIONS_DISTRIBUTED_MATRIX,EQUATIONS_MATRIX_DATA, &
                                      & ERR,ERROR,*999)
                                    SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                    CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)                                    
                                      !Loop over the non ghosted rows in the equations set
                                      DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                        RHS_VALUE=0.0_DP
                                        rhs_variable_dof=RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(equations_row_number)
                                        rhs_field_dof=RHS_VARIABLE%DOF_LIST(rhs_variable_dof)
                                        rhs_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(rhs_field_dof)
                                        SELECT CASE(rhs_boundary_condition)
                                        CASE(EQUATIONS_SET_NOT_FIXED)
                                          !Back substitute
                                          !Loop over the local columns of the equations matrix
                                          DO equations_column_idx=1,COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                            equations_column_number=COLUMN_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(equations_column_idx)
                                            variable_dof=equations_column_idx
                                            field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                            MATRIX_VALUE=EQUATIONS_MATRIX_DATA(equations_row_number+(equations_column_number-1)* &
                                              & EQUATIONS_MATRICES%TOTAL_NUMBER_OF_ROWS)
                                            DEPENDENT_VALUE=DEPENDENT_PARAMETERS(field_dof)                                        
                                            RHS_VALUE=RHS_VALUE+MATRIX_VALUE*DEPENDENT_VALUE
                                          ENDDO !equations_column_idx
                                        CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                                          !Do nothing
                                        CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                                          !Robin or is it Cauchy??? boundary conditions
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global boundary condition of "// &
                                            & TRIM(NUMBER_TO_VSTRING(rhs_boundary_condition,"*",ERR,ERROR))// &
                                            & " for RHS field dof number "//TRIM(NUMBER_TO_VSTRING(rhs_field_dof,"*",ERR,ERROR))// &
                                            & " is invalid."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                          SOURCE_VALUE=SOURCE_VECTOR_DATA(equations_row_number)
                                          RHS_VALUE=RHS_VALUE-SOURCE_VALUE
                                        ENDIF
                                        CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,rhs_field_dof, &
                                          & RHS_VALUE,ERR,ERROR,*999)
                                      ENDDO !equations_row_number
                                    CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)                      
                                    CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                    CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                      CALL DISTRIBUTED_MATRIX_STORAGE_LOCATIONS_GET(EQUATIONS_DISTRIBUTED_MATRIX,ROW_INDICES, &
                                        & COLUMN_INDICES,ERR,ERROR,*999)
                                      !Loop over the non-ghosted rows in the equations set
                                      DO equations_row_number=1,EQUATIONS_MAPPING%NUMBER_OF_ROWS
                                        RHS_VALUE=0.0_DP
                                        rhs_variable_dof=RHS_MAPPING%EQUATIONS_ROW_TO_RHS_DOF_MAP(equations_row_number)
                                        rhs_field_dof=RHS_VARIABLE%DOF_LIST(rhs_variable_dof)
                                        rhs_boundary_condition=FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(rhs_field_dof)
                                        SELECT CASE(rhs_boundary_condition)
                                        CASE(EQUATIONS_SET_NOT_FIXED)
                                          !Back substitute
                                          !Loop over the local columns of the equations matrix                                      
                                          DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                            ROW_INDICES(equations_row_number+1)-1
                                            equations_column_number=COLUMN_INDICES(equations_column_idx)
                                            variable_dof=equations_column_idx-ROW_INDICES(equations_row_number)+1
                                            field_dof=DEPENDENT_VARIABLE%DOF_LIST(variable_dof)
                                            MATRIX_VALUE=EQUATIONS_MATRIX_DATA(equations_column_idx)
                                            DEPENDENT_VALUE=DEPENDENT_PARAMETERS(field_dof)
                                            RHS_VALUE=RHS_VALUE+MATRIX_VALUE*DEPENDENT_VALUE
                                          ENDDO !equations_column_idx
                                        CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                                          !Do nothing
                                        CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                                          !Robin or is it Cauchy??? boundary conditions
                                          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global boundary condition of "// &
                                            & TRIM(NUMBER_TO_VSTRING(rhs_boundary_condition,"*",ERR,ERROR))// &
                                            & " for RHS field dof number "//TRIM(NUMBER_TO_VSTRING(rhs_field_dof,"*",ERR,ERROR))// &
                                            & " is invalid."
                                          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                                          SOURCE_VALUE=SOURCE_VECTOR_DATA(equations_row_number)
                                          RHS_VALUE=RHS_VALUE-SOURCE_VALUE
                                        ENDIF
                                        CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,rhs_field_dof, &
                                          & RHS_VALUE,ERR,ERROR,*999)
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
                            ELSE
                              CALL FLAG_ERROR("Dependent variable is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ENDDO !equations_matrix_idx
                        ELSE
                          CALL FLAG_ERROR("RHS variable is not associated.",ERR,ERROR,*999)
                        ENDIF
                        IF(ASSOCIATED(SOURCE_MAPPING)) THEN
                          CALL DISTRIBUTED_VECTOR_DATA_RESTORE(SOURCE_DISTRIBUTED_VECTOR,SOURCE_VECTOR_DATA,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Fixed conditions are not associated.",ERR,ERROR,*999)
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
            ELSE
              CALL FLAG_ERROR("Equations matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
          ENDIF
          !Restore the dependent field parameters
          CALL FIELD_PARAMETER_SET_RESTORE(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,DEPENDENT_PARAMETERS,ERR,ERROR,*999)
          !Start the update of the field parameters
          CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
          !Finish the update of the field parameters
          CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
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

  !>Finishes the process of creating an equation set on a region.
  SUBROUTINE EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        CALL FLAG_ERROR("Equations set has already been finished.",ERR,ERROR,*999)
      ELSE            
        !Finish the equations set specific setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION,ERR,ERROR,*999)
        !Finish the equations set specific geometry setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_GEOMETRY_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION,ERR,ERROR,*999)
        !Finish the equations set creation
        EQUATIONS_SET%EQUATIONS_SET_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
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

  !>Starts the process of creating an equations set defined by USER_NUMBER in the region identified by REGION.    
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
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    NULLIFY(NEW_EQUATIONS_SET)
    NULLIFY(NEW_EQUATIONS_SETS)

    CALL ENTERS("EQUATIONS_SET_CREATE_START",ERR,ERROR,*997)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        CALL EQUATIONS_SET_USER_NUMBER_FIND(USER_NUMBER,REGION,NEW_EQUATIONS_SET,ERR,ERROR,*997)
        IF(ASSOCIATED(NEW_EQUATIONS_SET)) THEN
          LOCAL_ERROR="Equations set user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
        ELSE
          NULLIFY(NEW_EQUATIONS_SET)
          IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
            IF(GEOM_FIBRE_FIELD%FIELD_FINISHED) THEN
              IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.GEOM_FIBRE_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
                !Allocate the new equtions set
                ALLOCATE(NEW_EQUATIONS_SET,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations set",ERR,ERROR,*999)
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
                !Start equations set specific setup
                CALL EQUATIONS_SET_SETUP(NEW_EQUATIONS_SET,EQUATIONS_SET_SETUP_INITIAL_TYPE,EQUATIONS_SET_SETUP_START_ACTION, &
                  & ERR,ERROR,*999)
                !Set up the equations set geometric fields
                CALL EQUATIONS_SET_GEOMETRY_INITIALISE(NEW_EQUATIONS_SET,ERR,ERROR,*999)
                IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
                  NEW_EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD
                  NULLIFY(NEW_EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)
                ELSE
                  NEW_EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD%GEOMETRIC_FIELD
                  NEW_EQUATIONS_SET%GEOMETRY%FIBRE_FIELD=>GEOM_FIBRE_FIELD
                ENDIF
                !Set up equations set specific geometry
                CALL EQUATIONS_SET_SETUP(NEW_EQUATIONS_SET,EQUATIONS_SET_SETUP_GEOMETRY_TYPE,EQUATIONS_SET_SETUP_START_ACTION, &
                  & ERR,ERROR,*999)                                
                !Add new equations set into list of equations set in the region
                ALLOCATE(NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets",ERR,ERROR,*999)
                DO equations_set_idx=1,REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
                  NEW_EQUATIONS_SETS(equations_set_idx)%PTR=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%PTR
                ENDDO !equations_set_idx
                NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1)%PTR=>NEW_EQUATIONS_SET
                IF(ASSOCIATED(REGION%EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
                REGION%EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
                REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1
                EQUATIONS_SET=>NEW_EQUATIONS_SET
              ELSE
                CALL FLAG_ERROR("The specified geometric field is not a geometric or fibre field",ERR,ERROR,*997)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The specified geometric field is not finished",ERR,ERROR,*997)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The specified geometric field is not associated",ERR,ERROR,*997)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The equations sets on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*997)
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

  !>Destroys an equations set identified by a user number on the give region and deallocates all memory.
  SUBROUTINE EQUATIONS_SET_DESTROY(USER_NUMBER,REGION,ERR,ERROR,*)

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

    CALL ENTERS("EQUATIONS_SET_DESTROY",ERR,ERROR,*999)

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
      CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,ERR,ERROR,*999)
      CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,ERR,ERROR,*999)
      CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,ERR,ERROR,*999)
      CALL EQUATIONS_SET_FIXED_CONDITIONS_FINALISE(EQUATIONS_SET%FIXED_CONDITIONS,ERR,ERROR,*999)
      CALL EQUATIONS_SET_EQUATIONS_FINALISE(EQUATIONS_SET%EQUATIONS,ERR,ERROR,*999)
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
    TYPE(EQUATIONS_MATRICES_LINEAR_TYPE), POINTER :: LINEAR_MATRICES
    TYPE(EQUATIONS_MATRICES_RHS_TYPE), POINTER :: RHS_VECTOR
    TYPE(EQUATIONS_MATRICES_SOURCE_TYPE), POINTER :: SOURCE_VECTOR
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("EQUATIONS_SET_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%CLASS)
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_ELEMENT_MATRIX_OUTPUT) THEN
          EQUATIONS_MATRICES=>EQUATIONS%EQUATIONS_MATRICES
          IF(ASSOCIATED(EQUATIONS_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",ERR,ERROR,*999)          
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
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_ELEMENT_MATRIX_OUTPUT) THEN
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
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_ELEMENT_MATRIX_OUTPUT) THEN
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
  
  !>Finalises the interpolation information for equations and deallocates all memory
  SUBROUTINE EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS_INTERPOLATION,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_INTERPOLATION_TYPE), POINTER :: EQUATIONS_INTERPOLATION !<A pointer to the equations interpolation to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_INTERPOLATION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_INTERPOLATION)) THEN
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%FIBRE_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(EQUATIONS_INTERPOLATION%SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_FINALISE(EQUATIONS_INTERPOLATION%GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_FINALISE(EQUATIONS_INTERPOLATION%DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_FINALISE(EQUATIONS_INTERPOLATION%FIBRE_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_FINALISE(EQUATIONS_INTERPOLATION%MATERIALS_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_FINALISE(EQUATIONS_INTERPOLATION%SOURCE_INTERP_POINT,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(EQUATIONS_INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(EQUATIONS_INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
      CALL FIELD_INTERPOLATED_POINT_METRICS_FINALISE(EQUATIONS_INTERPOLATION%FIBRE_INTERP_POINT_METRICS,ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_INTERPOLATION)
    ENDIF
       
    CALL EXITS("EQUATIONS_INTERPOLATION_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_INTERPOLATION_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_INTERPOLATION_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_INTERPOLATION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interpolation information for equations
  SUBROUTINE EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<The pointer to the equations to initialise the interpolation for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    
    CALL ENTERS("EQUATIONS_INTERPOLATION_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      EQUATIONS_SET=>EQUATIONS%EQUATIONS_SET
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        IF(ASSOCIATED(EQUATIONS%INTERPOLATION)) THEN
          CALL FLAG_ERROR("Interpolation is already associated for these equations.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(EQUATIONS%INTERPOLATION,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations interpolation",ERR,ERROR,*999)
          EQUATIONS%INTERPOLATION%EQUATIONS=>EQUATIONS
          NULLIFY(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS)
          NULLIFY(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT)
          NULLIFY(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS)
          NULLIFY(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS)
          NULLIFY(EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT_METRICS)
          
          EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          EQUATIONS%INTERPOLATION%FIBRE_FIELD=>EQUATIONS_SET%GEOMETRY%FIBRE_FIELD
          EQUATIONS%INTERPOLATION%DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
            EQUATIONS%INTERPOLATION%MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          ELSE
            NULLIFY(EQUATIONS%INTERPOLATION%MATERIALS_FIELD)
          ENDIF
          IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
            EQUATIONS%INTERPOLATION%SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
          ELSE
            NULLIFY(EQUATIONS%INTERPOLATION%SOURCE_FIELD)
          ENDIF
          
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%GEOMETRIC_FIELD, &
            & FIELD_STANDARD_VARIABLE_TYPE,EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_INITIALISE(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_PARAMETERS, &
            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT, &
            & EQUATIONS%INTERPOLATION%GEOMETRIC_INTERP_POINT_METRICS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%DEPENDENT_FIELD, &
            & FIELD_STANDARD_VARIABLE_TYPE,EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATED_POINT_INITIALISE(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_PARAMETERS, &
            & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT,ERR,ERROR,*999)
          IF(EQUATIONS%INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR. &
            & EQUATIONS%INTERPOLATION%DEPENDENT_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
            CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT, &
              & EQUATIONS%INTERPOLATION%DEPENDENT_INTERP_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%FIBRE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%FIBRE_FIELD, &
              & FIELD_STANDARD_VARIABLE_TYPE,EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_INITIALISE(EQUATIONS%INTERPOLATION%FIBRE_INTERP_PARAMETERS,  &
              &  EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_INITIALISE(EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT,  &
              &  EQUATIONS%INTERPOLATION%FIBRE_INTERP_POINT_METRICS,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%MATERIALS_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%MATERIALS_FIELD, &
              & FIELD_STANDARD_VARIABLE_TYPE,EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_INITIALISE(EQUATIONS%INTERPOLATION%MATERIALS_INTERP_PARAMETERS,  &
              & EQUATIONS%INTERPOLATION%MATERIALS_INTERP_POINT,ERR,ERROR,*999)
          ENDIF
          IF(ASSOCIATED(EQUATIONS%INTERPOLATION%SOURCE_FIELD)) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(EQUATIONS%INTERPOLATION%SOURCE_FIELD, &
              & FIELD_STANDARD_VARIABLE_TYPE,EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS,ERR,ERROR,*999)
            CALL FIELD_INTERPOLATED_POINT_INITIALISE(EQUATIONS%INTERPOLATION%SOURCE_INTERP_PARAMETERS, &
              & EQUATIONS%INTERPOLATION%SOURCE_INTERP_POINT,ERR,ERROR,*999)
          ENDIF
          
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations equation set is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_INTERPOLATION_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_INTERPOLATION_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_INTERPOLATION_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_INTERPOLATION_INITIALISE
  
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
      NULLIFY(EQUATIONS_SET%EQUATIONS_SETS)
      NULLIFY(EQUATIONS_SET%REGION)
      EQUATIONS_SET%CLASS=EQUATIONS_SET_NO_CLASS
      EQUATIONS_SET%TYPE=EQUATIONS_SET_NO_TYPE
      EQUATIONS_SET%SUBTYPE=EQUATIONS_SET_NO_SUBTYPE
      EQUATIONS_SET%LINEARITY=0
      EQUATIONS_SET%TIME_TYPE=0
      EQUATIONS_SET%SOLUTION_METHOD=0
      CALL EQUATIONS_SET_GEOMETRY_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
      CALL EQUATIONS_SET_DEPENDENT_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
      NULLIFY(EQUATIONS_SET%MATERIALS)
      NULLIFY(EQUATIONS_SET%SOURCE)
      NULLIFY(EQUATIONS_SET%ANALYTIC)
      NULLIFY(EQUATIONS_SET%FIXED_CONDITIONS)
      NULLIFY(EQUATIONS_SET%EQUATIONS)
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

  !>Applies the fixed conditions in an equation set to the dependent field in an equations set
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_APPLY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to apply the fixed conditions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS

    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_APPLY",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              CALL FIELD_PARAMETER_SET_COPY(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            ELSE
              CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations dependent field has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The equations set fixed conditions have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set fixed conditions is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_APPLY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_APPLY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_APPLY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_APPLY

  !
  !================================================================================================================================
  !

  !>Finish the creation of fixed conditions for an equation set.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create the fixed conditions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_computational_nodes,MPI_IERROR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS

    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
            DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
            IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
              !Start the transfer of the boundary conditions array. Note that the acutal boundary condition values will be
              !transferred in the assemble routines.
              CALL DISTRIBUTED_VECTOR_UPDATE_START(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              IF(number_computational_nodes>1) THEN
                !Transfer all the fixed conditions to all the computational nodes. At the moment just use an MPI_ALLREDUCE as the
                !dofs belonging to each computational node are not continuous (different field components) which prevents the
                !straightforward use of MPI_ALLGATHERV. The ALLREDUCE will have more transfers but we will see how bad it is later.
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS,DEPENDENT_DOFS_MAPPING% &
                  & NUMBER_OF_GLOBAL,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
              ENDIF
              !Finish the transfer of the boundary conditions
              CALL DISTRIBUTED_VECTOR_UPDATE_FINISH(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              !Finish equations set specific setting up 
              CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_FIXED_CONDITIONS_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION, &
                & ERR,ERROR,*999)
              !Apply the fixed conditions
              !Finish the fixed conditions
              FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED=.TRUE.
            ELSE
              CALL FLAG_ERROR("Dependent field dofs mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Dependent field has not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set fixed conditions is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Fixed conditions:",ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DEPENDENT_DOFS_MAPPING%NUMBER_OF_GLOBAL,8,8, &
        & FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS,'("  Global BCs:",8(X,I8))','(13X,8(X,I8))',ERR,ERROR,*999)      
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of fixed conditions for a problem.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the fixed conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD

    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%FIXED_CONDITIONS)) THEN
        CALL FLAG_ERROR("The equations set fixed conditions is already associated.",ERR,ERROR,*999)        
      ELSE
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          IF(DEPENDENT_FIELD%FIELD_FINISHED) THEN
            DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
            IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
              CALL EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
              ALLOCATE(EQUATIONS_SET%FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(DEPENDENT_DOFS_MAPPING%NUMBER_OF_GLOBAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global fixed conditions.",ERR,ERROR,*999)
              EQUATIONS_SET%FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS=EQUATIONS_SET_NOT_FIXED
              CALL DISTRIBUTED_VECTOR_CREATE_START(DEPENDENT_DOFS_MAPPING,EQUATIONS_SET%FIXED_CONDITIONS%BOUNDARY_CONDITIONS, &
                & ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_DATA_TYPE_SET(EQUATIONS_SET%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,MATRIX_VECTOR_INTG_TYPE, &
                & ERR,ERROR,*999)
              CALL DISTRIBUTED_VECTOR_CREATE_FINISH(EQUATIONS_SET%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
              !Initialise boundary conditions
              CALL DISTRIBUTED_VECTOR_ALL_VALUES_SET(EQUATIONS_SET%FIXED_CONDITIONS%BOUNDARY_CONDITIONS,EQUATIONS_SET_NOT_FIXED, &
                & ERR,ERROR,*999)
              !Perform equations set specific setup
              CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_FIXED_CONDITIONS_TYPE,EQUATIONS_SET_SETUP_START_ACTION, &
                & ERR,ERROR,*999)            
            ELSE
              CALL FLAG_ERROR("Dependent field dofs mapping is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The dependent field has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Equations set dependent field is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the fixed conditions for an equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_DESTROY(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the fixed conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%FIXED_CONDITIONS)) THEN
        CALL EQUATIONS_SET_FIXED_CONDITIONS_FINALISE(EQUATIONS_SET%FIXED_CONDITIONS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Equations set fixed conditions is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_DESTROY")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_DESTROY",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_DESTROY")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the fixed conditions for an equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_FINALISE(EQUATIONS_SET_FIXED_CONDITIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: EQUATIONS_SET_FIXED_CONDITIONS !<A pointer to the equations set fixed conditions to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET_FIXED_CONDITIONS)) THEN
      IF(ALLOCATED(EQUATIONS_SET_FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS))  &
        & DEALLOCATE(EQUATIONS_SET_FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS)
      CALL DISTRIBUTED_VECTOR_DESTROY(EQUATIONS_SET_FIXED_CONDITIONS%BOUNDARY_CONDITIONS,ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS_SET_FIXED_CONDITIONS)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the fixed conditions for a problem.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the fixed conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%FIXED_CONDITIONS)) THEN
        CALL FLAG_ERROR("Fixed conditions is already associated for this problem.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(EQUATIONS_SET%FIXED_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate problem fixed conditions.",ERR,ERROR,*999)
        EQUATIONS_SET%FIXED_CONDITIONS%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED=.FALSE.
        NULLIFY(EQUATIONS_SET%FIXED_CONDITIONS%BOUNDARY_CONDITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets fixed conditions for the equations set on the specified dofs.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS(EQUATIONS_SET,DOF_INDICES,CONDITIONS,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the fixed conditions for.
    INTEGER(INTG), INTENT(IN) :: DOF_INDICES(:) !<DOF_INDICES(i). The dof index for the i'th dof to set the fixed conditions for
    INTEGER(INTG), INTENT(IN) :: CONDITIONS(:) !<CONDITIONS(i). The fixed condition for the i'th dof.
    REAL(DP), INTENT(IN) :: VALUES(:) !<VALUES(i). The value of the fixed condition for the i'th dof.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,local_ny,global_ny
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          CALL FLAG_ERROR("Fixed conditions have been finished for this problem.",ERR,ERROR,*999)
        ELSE
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            IF(SIZE(DOF_INDICES,1)==SIZE(CONDITIONS,1)) THEN
              IF(SIZE(DOF_INDICES,1)==SIZE(VALUES,1)) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
                  IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
                    DO i=1,SIZE(DOF_INDICES,1)
!!TODO: set by global dof???
                      local_ny=DOF_INDICES(i)
                      IF(local_ny>0.AND.local_ny<=DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                        global_ny=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                        IF(DEPENDENT_DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(global_ny)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN
                          SELECT CASE(CONDITIONS(i))
                          CASE(EQUATIONS_SET_NOT_FIXED)
                            FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED                     
                          CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
!!TODO: need to think how initial conditions and increments for non-linear equations set are set.
                            FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                            CALL DISTRIBUTED_VECTOR_VALUES_SET(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,local_ny, &
                              & EQUATIONS_SET_FIXED_BOUNDARY_CONDITION,ERR,ERROR,*999)
                            CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,local_ny, &
                              & VALUES(i),ERR,ERROR,*999)
                          CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                          CASE DEFAULT
                            LOCAL_ERROR="The condition for index number "//TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))// &
                              & " is "//TRIM(NUMBER_TO_VSTRING(CONDITIONS(i),"*",ERR,ERROR))//" which is invalid."
                            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                          END SELECT
                        ELSE
                          !?Error
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Invalid dof indices. The dof number for index number "// &
                          & TRIM(NUMBER_TO_VSTRING(i,"*",ERR,ERROR))//" is "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                          & ". The allowed range is 1 to "// &
                          & TRIM(NUMBER_TO_VSTRING(DEPENDENT_DOFS_MAPPING%NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !i
                  ELSE
                    CALL FLAG_ERROR("Dependent field dofs mapping is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Dependent field is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The size of the dof indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                  & ") does not match the size of the values array ("// &
                  & TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))//")."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The size of the dof indices array ("//TRIM(NUMBER_TO_VSTRING(SIZE(DOF_INDICES,1),"*",ERR,ERROR))// &
                & ") does not match the size of the fixed conditions array ("// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(CONDITIONS,1),"*",ERR,ERROR))//")."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The equation set dependent field has not been finished.",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set fixed conditions are not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOFS

  !
  !================================================================================================================================
  !

  !>Sets a fixed condition for the equation set on the specified dof.
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1(EQUATIONS_SET,DOF_INDEX,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the fixed condition for
    INTEGER(INTG), INTENT(IN) :: DOF_INDEX !<The dof index to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The fixed condition to set
    REAL(DP), INTENT(IN) :: VALUE !<The value of the fixed condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          CALL FLAG_ERROR("Fixed conditions have been finished for this equations set.",ERR,ERROR,*999)
        ELSE
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
              IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
!!TODO: set by global dof index???
                local_ny=DOF_INDEX
                IF(local_ny>0.AND.local_ny<=DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                  global_ny=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                  SELECT CASE(CONDITION)
                  CASE(EQUATIONS_SET_NOT_FIXED)
                    FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                  CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                    FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                    CALL DISTRIBUTED_VECTOR_VALUES_SET(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,local_ny, &
                      & EQUATIONS_SET_FIXED_BOUNDARY_CONDITION,ERR,ERROR,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_DOF(DEPENDENT_FIELD,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,local_ny, &
                      & VALUE,ERR,ERROR,*999)
                  CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                    CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    LOCAL_ERROR="The specified condition of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  END SELECT
                ELSE
                  LOCAL_ERROR="The specified dof of "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                    & "is invalid. The allowed range is 1 to "// &
                    & TRIM(NUMBER_TO_VSTRING(DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The dependent field dofs mapping is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The dependent field is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem dependent has not been finished",ERR,ERROR,*999)              
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set fixed conditions are not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_SET_DOF1

  !
  !================================================================================================================================
  !

!!TODO: Check we are only setting local nodes.
  
  !>Sets a fixed condition for the equation set on the specified node. TODO update global condition as well. 
  SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE(EQUATIONS_SET,SET_TYPE,DERIVATIVE_NUMBER,NODE_NUMBER,COMPONENT_NUMBER, &
    & VARIABLE_NUMBER,CONDITION,VALUE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the fixed condition for
    INTEGER(INTG), INTENT(IN) :: SET_TYPE !<The field parameter set type
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBER !<The node_number to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<The variable number to set the fixed condition at
    INTEGER(INTG), INTENT(IN) :: CONDITION !<The fixed condition to set
    REAL(DP), INTENT(IN) :: VALUE !<The value of the fixed condition to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_ny,global_ny,GLOBAL_NODE_NUMBER
    LOGICAL :: NODE_EXISTS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOFS_MAPPING
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(EQUATIONS_SET_FIXED_CONDITIONS_TYPE), POINTER :: FIXED_CONDITIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      FIXED_CONDITIONS=>EQUATIONS_SET%FIXED_CONDITIONS
      IF(ASSOCIATED(FIXED_CONDITIONS)) THEN
        IF(FIXED_CONDITIONS%FIXED_CONDITIONS_FINISHED) THEN
          CALL FLAG_ERROR("Fixed conditions have been finished for this equations set.",ERR,ERROR,*999)
        ELSE
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
            DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
              DEPENDENT_DOFS_MAPPING=>DEPENDENT_FIELD%MAPPINGS%DOMAIN_MAPPING
              IF(ASSOCIATED(DEPENDENT_DOFS_MAPPING)) THEN
                !!TODO: set by global dof index???
                CALL NODE_CHECK_EXISTS(NODE_NUMBER,DEPENDENT_FIELD%REGION,NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                IF(.NOT.NODE_EXISTS) THEN
                  CALL FLAG_ERROR("Invalid node number.",ERR,ERROR,*999)
                ELSE
                  local_ny=DEPENDENT_FIELD%VARIABLES(VARIABLE_NUMBER)%COMPONENTS(COMPONENT_NUMBER)%PARAM_TO_DOF_MAP% &
                      & NODE_PARAM2DOF_MAP(DERIVATIVE_NUMBER,NODE_NUMBER,0)
                  IF(local_ny>0.AND.local_ny<=DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL) THEN
                    global_ny=DEPENDENT_DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_ny)
                    SELECT CASE(CONDITION)
                    CASE(EQUATIONS_SET_NOT_FIXED)
                      FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_NOT_FIXED
                    CASE(EQUATIONS_SET_FIXED_BOUNDARY_CONDITION)
                      FIXED_CONDITIONS%GLOBAL_BOUNDARY_CONDITIONS(global_ny)=EQUATIONS_SET_FIXED_BOUNDARY_CONDITION
                      CALL DISTRIBUTED_VECTOR_VALUES_SET(FIXED_CONDITIONS%BOUNDARY_CONDITIONS,local_ny, &
                        & EQUATIONS_SET_FIXED_BOUNDARY_CONDITION,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,SET_TYPE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_NODE(DEPENDENT_FIELD,SET_TYPE,DERIVATIVE_NUMBER, &
                        & GLOBAL_NODE_NUMBER,COMPONENT_NUMBER,VARIABLE_NUMBER,VALUE,ERR,ERROR,*999)
                      CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,SET_TYPE,ERR,ERROR,*999)
                    CASE(EQUATIONS_SET_MIXED_BOUNDARY_CONDITION)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The specified condition of "//TRIM(NUMBER_TO_VSTRING(CONDITION,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT
                  ELSE
                    LOCAL_ERROR="The specified dof of "//TRIM(NUMBER_TO_VSTRING(local_ny,"*",ERR,ERROR))// &
                      & "is invalid. The allowed range is 1 to "// &
                      & TRIM(NUMBER_TO_VSTRING(DEPENDENT_DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ELSE
                CALL FLAG_ERROR("The dependent field dofs mapping is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The dependent field is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("The problem dependent has not been finished",ERR,ERROR,*999)              
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set fixed conditions are not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_FIXED_CONDITIONS_SET_NODE

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

  !>Finalise the equations set linear data and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE(EQUATIONS_LINEAR_DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_LINEAR_DATA_TYPE), POINTER :: EQUATIONS_LINEAR_DATA !<A pointer to the equations linear data to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_LINEAR_DATA)) THEN      
      DEALLOCATE(EQUATIONS_LINEAR_DATA)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the linear data information for an equations
  SUBROUTINE EQUATIONS_SET_EQUATIONS_LINEAR_DATA_INITIALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<The pointer to the equations to initialise the linear data for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%LINEAR_DATA)) THEN
        CALL FLAG_ERROR("Linear data is already associated for these equations.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS%LINEAR_DATA,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations linear data.",ERR,ERROR,*999)
        EQUATIONS%LINEAR_DATA%EQUATIONS=>EQUATIONS        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE(EQUATIONS%LINEAR_DATA,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_LINEAR_DATA_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_LINEAR_DATA_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the field component interpolation for a materials field of a problem.
  SUBROUTINE EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET,COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the materials component interpolation for
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component of the material field to set the interpolation for
    INTEGER(INTG), INTENT(IN) :: INTERPOLATION_TYPE !<The interpolation type to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Equations set materials has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
            & COMPONENT_NUMBER,INTERPOLATION_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_COMPONENT_INTERPOLATION_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field component mesh component for a materials field of a problem.
  SUBROUTINE EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET,COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the materials field component mesh component for
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The component number of the equations set materials field to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Equations set materials has been finished",ERR,ERROR,*999)
        ELSE
          CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,FIELD_STANDARD_VARIABLE_TYPE, &
            & COMPONENT_NUMBER,MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set materials is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_COMPONENT_MESH_COMPONENT_SET

  !
  !================================================================================================================================
  !

  !>Finish the creation of materials for an equations set.
  SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the materials field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_MATERIALS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Equations set materials has already been finished",ERR,ERROR,*999)
        ELSE
          !Finish equations set specific startup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_MATERIALS_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION, &
            & ERR,ERROR,*999)
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

  !>Start the creation of materials for a problem.
  SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_SET_MATERIALS_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL FLAG_ERROR("The equations set materials is already associated",ERR,ERROR,*998)        
      ELSE
        !Initialise the equations set materials
        CALL EQUATIONS_SET_MATERIALS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
        !Start equations set specific startup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_MATERIALS_TYPE,EQUATIONS_SET_SETUP_START_ACTION, &
          & ERR,ERROR,*999)
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

  !>Destroy the materials for an equations set.
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
      IF(ASSOCIATED(EQUATIONS_SET_MATERIALS%MATERIALS_FIELD))  &
        & CALL FIELD_DESTROY(EQUATIONS_SET_MATERIALS%MATERIALS_FIELD,ERR,ERROR,*999)
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
  !================================================================================================================================
  !

  !>Returns a pointer to the material field of the materials for a problem.
  SUBROUTINE EQUATIONS_SET_MATERIALS_MATERIAL_FIELD_GET(EQUATIONS_SET,MATERIAL_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to get the material field for
    TYPE(FIELD_TYPE), POINTER :: MATERIAL_FIELD !<On return, a pointer to the materials field for the equations set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_MATERIALS_MATERIAL_FIELD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(MATERIAL_FIELD)) THEN
        CALL FLAG_ERROR("Material field is already associated.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
          IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
            MATERIAL_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          ELSE
            CALL FLAG_ERROR("Materials has not been finished for this equations set.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Materials is not associated for this equations set.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_MATERIAL_FIELD_GET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_MATERIAL_FIELD_GET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_MATERIAL_FIELD_GET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_MATERIAL_FIELD_GET

  !
  !================================================================================================================================
  !

  !>Sets/changes the field scaling for a materials field of an equations set.
  SUBROUTINE EQUATIONS_SET_MATERIALS_SCALING_SET(EQUATIONS_SET,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the scaling on the material field
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE !<The scaling type to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_MATERIALS_SCALING_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FLAG_ERROR("Problem materials has been finished.",ERR,ERROR,*999)
        ELSE
          CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,SCALING_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set materials is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_MATERIALS_SCALING_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_MATERIALS_SCALING_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_MATERIALS_SCALING_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_SCALING_SET

  !
  !================================================================================================================================
  !

  !>Finalise the equations set nonlinear data and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE(EQUATIONS_NONLINEAR_DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_NONLINEAR_DATA_TYPE), POINTER :: EQUATIONS_NONLINEAR_DATA !<A pointer to the equations nonlinear data to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_NONLINEAR_DATA)) THEN      
      DEALLOCATE(EQUATIONS_NONLINEAR_DATA)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nonlinear data information for an equations
  SUBROUTINE EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_INITIALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<The pointer to the equations to initialise the nonlinear data for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%NONLINEAR_DATA)) THEN
        CALL FLAG_ERROR("Nonlinear data is already associated for these equations.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS%NONLINEAR_DATA,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations nonlinear data.",ERR,ERROR,*999)
        EQUATIONS%NONLINEAR_DATA%EQUATIONS=>EQUATIONS        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE(EQUATIONS%NONLINEAR_DATA,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the field component mesh component for a dependent field of an equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET,VARIABLE_NUMBER,COMPONENT_NUMBER, &
    & MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the dependent field component mesh component
    INTEGER(INTG), INTENT(IN) :: VARIABLE_NUMBER !<The dependent field variable number to set \todo this should be variable type???
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The dependent field component number to set
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Equations set dependent has been finished",ERR,ERROR,*999)
      ELSE
        CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,VARIABLE_NUMBER,COMPONENT_NUMBER, &
          & MESH_COMPONENT_NUMBER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_COMPONENT_MESH_COMPONENT_SET

  !
  !================================================================================================================================
  !

  !>Finish the creation of a dependent variables for an equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Equations set dependent has already been finished",ERR,ERROR,*999)
      ELSE
        !Finish equations set specific setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_DEPENDENT_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION, &
          & ERR,ERROR,*999)
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

  !>Start the creation of dependent variables for an equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a dependent field on
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_SET_DEPENDENT_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      !Start the equations set specfic solution setup
      CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_DEPENDENT_TYPE,EQUATIONS_SET_SETUP_START_ACTION, &
        & ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations_set is not associated",ERR,ERROR,*998)
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

  !>Returns a pointer to the dependent field of the dependent variables for an equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET(EQUATIONS_SET,DEPENDENT_FIELD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equation set to get the dependent field for
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD !<On return, a pointer to the dependent field for the equations set. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
        CALL FLAG_ERROR("Dependent field is already associated.",ERR,ERROR,*999)
      ELSE
        IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        ELSE
          CALL FLAG_ERROR("Dependent field has not been finished for this equations set.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_DEPENDENT_FIELD_GET

  !
  !================================================================================================================================
  !
  
  !>Destroy the dependent variables for an equations sety.
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

    IF(ASSOCIATED(EQUATIONS_SET_DEPENDENT%DEPENDENT_FIELD)) &
      & CALL FIELD_DESTROY(EQUATIONS_SET_DEPENDENT%DEPENDENT_FIELD,ERR,ERROR,*999)
    
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
      NULLIFY(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
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

  !>Sets/changes the field scaling for a dependent field of an equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_SCALING_SET(EQUATIONS_SET,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to the dependent field scaling for.
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE !<The scaling type to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_DEPENDENT_SCALING_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FLAG_ERROR("Equations set dependent has been finished",ERR,ERROR,*999)
      ELSE
        CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,SCALING_TYPE,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_DEPENDENT_SCALING_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_DEPENDENT_SCALING_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_DEPENDENT_SCALING_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_SCALING_SET

  !
  !================================================================================================================================
  !

  !>Sets up the specifices for an equation set.
  SUBROUTINE EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the setup on
    INTEGER(INTG), INTENT(IN) :: SETUP_TYPE !<The type of setup \see EQUATIONS_ROUTINES_SetupTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(IN) :: ACTION_TYPE !<The setup type action \see EQUATIONS_ROUTINES_SetupActionTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      SELECT CASE(EQUATIONS_SET%CLASS)
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_EQUATIONS_SET_SETUP(EQUATIONS_SET,SETUP_TYPE,ACTION_TYPE,ERR,ERROR,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%CLASS,"*",ERR,ERROR))//" is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
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

!!TODO: sort this out call of problem final setup???? 
  !>Finish the creation of equations for the equations set.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%EQUATIONS_FINISHED) THEN
          CALL FLAG_ERROR("Equations has already been finished.",ERR,ERROR,*999)
        ELSE
          !Finish the equations specific solution setup.
          CALL EQUATIONS_SET_SETUP(EQUATIONS%EQUATIONS_SET,EQUATIONS_SET_SETUP_EQUATIONS_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION, &
            & ERR,ERROR,*999)
          !Finish the problem solution creation
          EQUATIONS%EQUATIONS_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
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

  !>Start the creation of equations for the equation set. \todo Should this return a pointer to the equations???
  SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) THEN
        CALL FLAG_ERROR("The equations is already associated for the equations set.",ERR,ERROR,*999)        
      ELSE
        !Initialise the equations
        CALL EQUATIONS_SET_EQUATIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
        !Start the equations set specific solution setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_EQUATIONS_TYPE,EQUATIONS_SET_SETUP_START_ACTION, &
          & ERR,ERROR,*999)
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

  !>Finalise the equations and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_FINALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<A pointer to the equations to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS)) THEN
      CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
      CALL EQUATIONS_SET_EQUATIONS_LINEAR_DATA_FINALISE(EQUATIONS%LINEAR_DATA,ERR,ERROR,*999)
      CALL EQUATIONS_SET_EQUATIONS_NONLINEAR_DATA_FINALISE(EQUATIONS%NONLINEAR_DATA,ERR,ERROR,*999)
      CALL EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE(EQUATIONS%TIME_DATA,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS%EQUATIONS_MATRICES)) CALL EQUATIONS_MATRICES_DESTROY(EQUATIONS%EQUATIONS_MATRICES,ERR,ERROR,*999)
      DEALLOCATE(EQUATIONS)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the equations set.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS_SET,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The sparsity type to set \see EQUATIONS_ROUTINES_SparsityTypes,PROBLEM_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%EQUATIONS_FINISHED) THEN
          CALL FLAG_ERROR("Equations has already been finished.",ERR,ERROR,*999)
        ELSE
          SELECT CASE(SPARSITY_TYPE)
          CASE(EQUATIONS_SET_SPARSE_MATRICES)
            EQUATIONS%SPARSITY_TYPE=EQUATIONS_SET_SPARSE_MATRICES
          CASE(EQUATIONS_SET_FULL_MATRICES)
            EQUATIONS%SPARSITY_TYPE=EQUATIONS_SET_FULL_MATRICES
          CASE DEFAULT
            LOCAL_ERROR="The specified sparsity type of "//TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_SPARSITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Initialises the equations for an equations set.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_INITIALISE(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) THEN
        CALL FLAG_ERROR("Equations is already associated for this equations set.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(EQUATIONS_SET%EQUATIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations.",ERR,ERROR,*999)
        EQUATIONS_SET%EQUATIONS%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%EQUATIONS%OUTPUT_TYPE=EQUATIONS_SET_NO_OUTPUT
        EQUATIONS_SET%EQUATIONS%SPARSITY_TYPE=EQUATIONS_SET_SPARSE_MATRICES
        NULLIFY(EQUATIONS_SET%EQUATIONS%INTERPOLATION)
        NULLIFY(EQUATIONS_SET%EQUATIONS%LINEAR_DATA)
        NULLIFY(EQUATIONS_SET%EQUATIONS%NONLINEAR_DATA)
        NULLIFY(EQUATIONS_SET%EQUATIONS%TIME_DATA)
        NULLIFY(EQUATIONS_SET%EQUATIONS%EQUATIONS_MAPPING)
        NULLIFY(EQUATIONS_SET%EQUATIONS%EQUATIONS_MATRICES)
        EQUATIONS_SET%EQUATIONS%EQUATIONS_FINISHED=.FALSE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_INITIALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATION_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_INITIALISE
  
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
          SELECT CASE(EQUATIONS_SET%LINEARITY)
          CASE(EQUATIONS_SET_LINEAR)            
            CALL FLAG_ERROR("Can not evaluate a Jacobian for a linear equations set.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_NONLINEAR)
            SELECT CASE(EQUATIONS_SET%TIME_TYPE)
            CASE(EQUATIONS_SET_STATIC)
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
            CASE(EQUATIONS_SET_DYNAMIC)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_QUASISTATIC)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set time type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TIME_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations set linearity of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%LINEARITY,"*",ERR,ERROR))//" is invalid."
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
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
!!Do we need to transfer parameter sets???
            !Problem interpolation setup
            CALL EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS,ERR,ERROR,*999)
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL
              ne=ELEMENTS_MAPPING%INTERNAL_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finish the transfer of the solution values.
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
!!TODO: sort out combining boundary and ghost list
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY
              ne=ELEMENTS_MAPPING%BOUNDARY_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_GHOST
              ne=ELEMENTS_MAPPING%GHOST_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_JACOBIAN_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx          
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finalise the problem interpolation
            CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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

  !>Sets/changes the output type for the equations set.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS_SET,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the output type for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The output type to set \see EQUATIONS_ROUTINES_OutputTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(EQUATIONS%EQUATIONS_FINISHED) THEN
          CALL FLAG_ERROR("Equations has already been finished for this equations set.",ERR,ERROR,*999)
        ELSE
          SELECT CASE(OUTPUT_TYPE)
          CASE(EQUATIONS_SET_NO_OUTPUT)
            EQUATIONS%OUTPUT_TYPE=EQUATIONS_SET_NO_OUTPUT
          CASE(EQUATIONS_SET_TIMING_OUTPUT)
            EQUATIONS%OUTPUT_TYPE=EQUATIONS_SET_TIMING_OUTPUT
          CASE(EQUATIONS_SET_MATRIX_OUTPUT)
            EQUATIONS%OUTPUT_TYPE=EQUATIONS_SET_MATRIX_OUTPUT
          CASE(EQUATIONS_SET_ELEMENT_MATRIX_OUTPUT)
            EQUATIONS%OUTPUT_TYPE=EQUATIONS_SET_ELEMENT_MATRIX_OUTPUT
          CASE DEFAULT
            LOCAL_ERROR="The specified output type of "//TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid"
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_OUTPUT_TYPE_SET
  
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
          SELECT CASE(EQUATIONS_SET%LINEARITY)
          CASE(EQUATIONS_SET_LINEAR)            
            CALL FLAG_ERROR("Can not evaluate a residual for a linear equations set.",ERR,ERROR,*999)
          CASE(EQUATIONS_SET_NONLINEAR)
            SELECT CASE(EQUATIONS_SET%TIME_TYPE)
            CASE(EQUATIONS_SET_STATIC)
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
            CASE(EQUATIONS_SET_DYNAMIC)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(EQUATIONS_SET_QUASISTATIC)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set time type of "// &
                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%TIME_TYPE,"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_NONLINEAR_BCS)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equations set linearity of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%LINEARITY,"*",ERR,ERROR))//" is invalid."
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
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
              CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
              CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
            ENDIF
            !!Do we need to transfer parameter sets???
!Problem interpolation setup
            CALL EQUATIONS_INTERPOLATION_INITIALISE(EQUATIONS,ERR,ERROR,*999)
            !Initialise the matrices and rhs vector
            CALL EQUATIONS_MATRICES_VALUES_INITIALISE(EQUATIONS_MATRICES,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,ERR,ERROR,*999)
            !Assemble the elements
            !Allocate the element matrices 
            CALL EQUATIONS_MATRICES_ELEMENT_INITIALISE(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ELEMENTS_MAPPING=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
              & MAPPINGS%ELEMENTS
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL
              ne=ELEMENTS_MAPPING%INTERNAL_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx                  
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finish the transfer of the solution values.
            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
!!TODO: sort out combining boundary and ghost list
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY
              ne=ELEMENTS_MAPPING%BOUNDARY_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx
            DO element_idx=1,ELEMENTS_MAPPING%NUMBER_OF_GHOST
              ne=ELEMENTS_MAPPING%GHOST_LIST(element_idx)
              NUMBER_OF_TIMES=NUMBER_OF_TIMES+1
              CALL EQUATIONS_MATRICES_ELEMENT_CALCULATE(EQUATIONS_MATRICES,ne,ERR,ERROR,*999)
              CALL EQUATIONS_SET_FINITE_ELEMENT_RESIDUAL_EVALUATE(EQUATIONS_SET,ne,ERR,ERROR,*999)
              CALL EQUATIONS_MATRICES_ELEMENT_ADD(EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDDO !element_idx          
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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
            !Finalise the problem interpolation
            CALL EQUATIONS_INTERPOLATION_FINALISE(EQUATIONS%INTERPOLATION,ERR,ERROR,*999)
            !Output equations matrices and RHS vector if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_MATRIX_OUTPUT) THEN
              CALL EQUATIONS_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,EQUATIONS_MATRICES,ERR,ERROR,*999)
            ENDIF
            !Output timing information if required
            IF(EQUATIONS%OUTPUT_TYPE>=EQUATIONS_SET_TIMING_OUTPUT) THEN
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

  !>Finish the creation of a source for an equation set.
  SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a souce for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        IF(EQUATIONS_SET%SOURCE%SOURCE_FINISHED) THEN
          CALL FLAG_ERROR("Equations set source has already been finished.",ERR,ERROR,*999)
        ELSE
          !Finish the equation set specific source setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_SOURCE_TYPE,EQUATIONS_SET_SETUP_FINISH_ACTION, &
            & ERR,ERROR,*999)
          !Finish the source creation
          EQUATIONS_SET%SOURCE%SOURCE_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("The equations set source is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*999)
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

  !>Start the creation of a source for an equations set.
  SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_START(EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a source for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("EQUATIONS_SET_SOURCE_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL FLAG_ERROR("The equations set source is already associated.",ERR,ERROR,*998)        
      ELSE
        !Initialise the equations set source
        CALL EQUATIONS_SET_SOURCE_INITIALISE(EQUATIONS_SET,ERR,ERROR,*999)
        !Start the equation set specific source setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_SOURCE_TYPE,EQUATIONS_SET_SETUP_START_ACTION,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
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

  !>Destroy the source for an equations set.
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
      IF(ASSOCIATED(EQUATIONS_SET_SOURCE%SOURCE_FIELD)) CALL FIELD_DESTROY(EQUATIONS_SET_SOURCE%SOURCE_FIELD,ERR,ERROR,*999)
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
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated",ERR,ERROR,*998)
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

  !>Sets/changes the field scaling for a source field of an equations set.
  SUBROUTINE EQUATIONS_SET_SOURCE_SCALING_SET(EQUATIONS_SET,SCALING_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the scaling on the source field
    INTEGER(INTG), INTENT(IN) :: SCALING_TYPE !<The scaling type to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("EQUATIONS_SET_SOURCE_SCALING_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        IF(EQUATIONS_SET%SOURCE%SOURCE_FINISHED) THEN
          CALL FLAG_ERROR("Equations set source has been finished.",ERR,ERROR,*999)
        ELSE
          CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%SOURCE%SOURCE_FIELD,SCALING_TYPE,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Equations set source is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_SOURCE_SCALING_SET")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SOURCE_SCALING_SET",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SOURCE_SCALING_SET")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_SCALING_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation set specification i.e., equation set class, type and subtype for an equation set identified by a user number.
  SUBROUTINE EQUATIONS_SET_SPECIFICATION_SET_NUMBER(USER_NUMBER,REGION,EQUATIONS_SET_CLASS,EQUATIONS_SET_TYPE_, &
    & EQUATIONS_SET_SUBTYPE,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the equations set
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the equations set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_CLASS !<The equations set class to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_TYPE_ !<The equations set equation type to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equations set equation subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(EQUATIONS_SET)
    
    CALL ENTERS("EQUATIONS_SET_SPECIFICATION_SET_NUMBER",ERR,ERROR,*999)

!!TODO: Take in region number here and user FIND_REGION_NUMBER. This would require FIND_REGION_NUMBER to be moved from
!!REGION_ROUTINES otherwise there will be a circular module reference.
    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SET_USER_NUMBER_FIND(USER_NUMBER,REGION,EQUATIONS_SET,ERR,ERROR,*999)
      IF(ASSOCIATED(EQUATIONS_SET)) THEN
        CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_CLASS,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE, &
          & ERR,ERROR,*999)
      ELSE
        LOCAL_ERROR="Equation set user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
          & " is not defined on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
           
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_SET_NUMBER")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SPECIFICATION_SET_NUMBER",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_SET_NUMBER")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SPECIFICATION_SET_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the equations set specification i.e., equations set class, type and subtype for a equations set identified by a pointer.
  SUBROUTINE EQUATIONS_SET_SPECIFICATION_SET_PTR(EQUATIONS_SET,EQUATIONS_SET_CLASS,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_CLASS !<The equations set class to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_TYPE_ !<The equations set type to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SUBTYPE !<The equations set subtype to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("EQUATIONS_SET_SPECIFICATION_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        CALL FLAG_ERROR("Equations set has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(EQUATIONS_SET_CLASS)
        CASE(EQUATIONS_SET_ELASTICITY_CLASS)
          CALL ELASTICITY_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
          CALL CLASSICAL_FIELD_EQUATIONS_SET_CLASS_TYPE_SET(EQUATIONS_SET,EQUATIONS_SET_TYPE_,EQUATIONS_SET_SUBTYPE,ERR,ERROR,*999)
        CASE(EQUATIONS_SET_MODAL_CLASS)
          CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set class "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_CLASS,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_SET_PTR")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_SPECIFICATION_SET_PTR",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_SPECIFICATION_SET_PTR")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SPECIFICATION_SET_PTR
  
  !
  !================================================================================================================================
  !
  
  !>Finalises the time data information for an equations and deallocates all memory
  SUBROUTINE EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE(TIME_DATA,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TIME_DATA_TYPE), POINTER :: TIME_DATA !<A pointer to the equations set time data to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_DATA)) THEN
      DEALLOCATE(TIME_DATA)
    ENDIF
        
    CALL EXITS("EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE")
    RETURN
999 CALL ERRORS("EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the time data information for an equations
  SUBROUTINE EQUATIONS_SET_EQUATIONS_TIME_DATA_INITIALISE(EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS !<The pointer to the equations to initialise the time data for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("EQUATIONS_SET_EQUATIONS_TIME_DATA_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(EQUATIONS)) THEN
      IF(ASSOCIATED(EQUATIONS%TIME_DATA)) THEN
        CALL FLAG_ERROR("Time data is already associated for these equations.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(EQUATIONS%TIME_DATA,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate equations time data.",ERR,ERROR,*999)
        EQUATIONS%TIME_DATA%EQUATIONS=>EQUATIONS        
      ENDIF
    ELSE
      CALL FLAG_ERROR("Equations is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("EQUATIONS_SET_EQUATIONS_TIME_DATA_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_EQUATIONS_TIME_DATA_FINALISE(EQUATIONS%TIME_DATA,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EQUATIONS_SET_EQUATIONS_TIME_DATA_INITIALISE",ERR,ERROR)
    CALL EXITS("EQUATIONS_SET_EQUATIONS_TIME_DATA_INITIALISE")
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_TIME_DATA_INITIALISE

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
    INTEGER(INTG) ::USER_NUMBER

    CALL ENTERS("EQUATIONS_SETS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        DO WHILE(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>0)
          USER_NUMBER=REGION%EQUATIONS_SETS%EQUATIONS_SETS(1)%PTR%USER_NUMBER
          CALL EQUATIONS_SET_DESTROY(USER_NUMBER,REGION,ERR,ERROR,*999)
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
