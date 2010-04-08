!> \file
!> $Id: interface_conditions_routines.f90 690 2009-09-30 23:27:16Z chrispbradley $
!> \author Chris Bradley
!> \brief This module contains all interface conditions routines.
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

!>This module contains all interface conditions routines.
MODULE INTERFACE_CONDITIONS_ROUTINES

  USE BASE_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_EQUATIONS_ROUTINES
  USE INTERFACE_MATRICES_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC INTERFACE_CONDITION_CREATE_FINISH,INTERFACE_CONDITION_CREATE_START

  PUBLIC INTERFACE_CONDITION_DESTROY

  PUBLIC INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH,INTERFACE_CONDITION_EQUATIONS_CREATE_START

  PUBLIC INTERFACE_CONDITION_EQUATIONS_DESTROY

  PUBLIC INTERFACE_CONDITION_EQUATIONS_SET_ADD

  PUBLIC INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH,INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START

  PUBLIC INTERFACE_CONDITION_METHOD_GET,INTERFACE_CONDITION_METHOD_SET

  PUBLIC INTERFACE_CONDITION_OPERATOR_GET,INTERFACE_CONDITION_OPERATOR_SET

  PUBLIC INTERFACE_CONDITION_USER_NUMBER_FIND

  PUBLIC INTERFACE_CONDITIONS_FINALISE,INTERFACE_CONDITIONS_INITIALISE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an interface condition.
  SUBROUTINE INTERFACE_CONDITION_ASSEMBLE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
          SELECT CASE(INTERFACE_CONDITION%METHOD)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
            SELECT CASE(INTERFACE_CONDITION%SOLUTION_METHOD)
            CASE(INTERFACE_CONDITION_FEM_SOLUTION_METHOD)
              CALL INTERFACE_CONDITION_ASSEMBLE_FEM(INTERFACE_CONDITION,ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_FD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_FV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_GFEM_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_GFD_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_GFV_SOLUTION_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface condition solution method of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%SOLUTION_METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The interface condition method of "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FLAG_ERROR("Interface equations have not been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_ASSEMBLE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_ASSEMBLE

  !
  !================================================================================================================================
  !
  
  !>Assembles the interface matricesand rhs for using the finite element method.
  SUBROUTINE INTERFACE_CONDITION_ASSEMBLE_FEM(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,ne,NUMBER_OF_TIMES
    REAL(SP) :: ELEMENT_USER_ELAPSED,ELEMENT_SYSTEM_ELAPSED,USER_ELAPSED,USER_TIME1(1),USER_TIME2(1),USER_TIME3(1),USER_TIME4(1), &
      & USER_TIME5(1),USER_TIME6(1),SYSTEM_ELAPSED,SYSTEM_TIME1(1),SYSTEM_TIME2(1),SYSTEM_TIME3(1),SYSTEM_TIME4(1), &
      & SYSTEM_TIME5(1),SYSTEM_TIME6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD
    
!#ifdef TAUPROF
!    CHARACTER(28) :: CVAR
!    INTEGER :: PHASE(2) = (/ 0, 0 /)
!    SAVE PHASE
!#endif

    CALL ENTERS("INTERFACE_CONDITION_ASSEMBLE_FEM",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
        IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
            INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
            IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              !Initialise the matrices and rhs vector
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_VALUES_INITIALISE()")
#endif
              CALL INTERFACE_MATRICES_VALUES_INITIALISE(INTERFACE_MATRICES,0.0_DP,ERR,ERROR,*999)
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_VALUES_INITIALISE()")
#endif
              !Assemble the elements
              !Allocate the element matrices 
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_INITIALISE()")
#endif
              CALL INTERFACE_MATRICES_ELEMENT_INITIALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
              ELEMENTS_MAPPING=>LAGRANGE_FIELD%DECOMPOSITION%DOMAIN(LAGRANGE_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
                & MAPPINGS%ELEMENTS
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_INITIALISE()")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME2(1)-USER_TIME1(1)
                SYSTEM_ELAPSED=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for interface equations setup and initialisation = ", &
                  & USER_ELAPSED,ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for interface equations setup and initialisation = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
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
                CALL INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ne,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ne,ERR,ERROR,*999)
                CALL INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*999)
!#ifdef TAUPROF
!              CALL TAU_PHASE_STOP(PHASE)
!#endif
              ENDDO !element_idx
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("Internal Elements Loop")
#endif

              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
                CALL CPU_TIMER(USER_CPU,USER_TIME3,ERR,ERROR,*999)
                CALL CPU_TIMER(SYSTEM_CPU,SYSTEM_TIME3,ERR,ERROR,*999)
                USER_ELAPSED=USER_TIME3(1)-USER_TIME2(1)
                SYSTEM_ELAPSED=SYSTEM_TIME3(1)-SYSTEM_TIME2(1)
                ELEMENT_USER_ELAPSED=USER_ELAPSED
                ELEMENT_SYSTEM_ELAPSED=SYSTEM_ELAPSED
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"User time for internal interface equations assembly = ", &
                  & USER_ELAPSED, ERR,ERROR,*999)
                CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"System time for internal interface equations assembly = ", &
                  & SYSTEM_ELAPSED,ERR,ERROR,*999)
              ENDIF
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
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
                CALL INTERFACE_MATRICES_ELEMENT_CALCULATE(INTERFACE_MATRICES,ne,ERR,ERROR,*999)
                CALL INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ne,ERR,ERROR,*999)
                CALL INTERFACE_MATRICES_ELEMENT_ADD(INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDDO !element_idx
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("Boundary and Ghost Elements Loop")
#endif
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
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
              CALL TAU_STATIC_PHASE_START("INTERFACE_MATRICES_ELEMENT_FINALISE()")
#endif
              CALL INTERFACE_MATRICES_ELEMENT_FINALISE(INTERFACE_MATRICES,ERR,ERROR,*999)
#ifdef TAUPROF
              CALL TAU_STATIC_PHASE_STOP("INTERFACE_MATRICES_ELEMENT_FINALISE()")
#endif
              !Output equations matrices and vector if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_MATRIX_OUTPUT) THEN
                CALL INTERFACE_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,INTERFACE_MATRICES,ERR,ERROR,*999)
              ENDIF
              !Output timing information if required
              IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_TIMING_OUTPUT) THEN
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
              CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Lagrange field is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition Lagrange is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE_FEM")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_ASSEMBLE_FEM",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_ASSEMBLE_FEM")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_ASSEMBLE_FEM

  !
  !==================================================================================================================================
  !

  !>Finishes the process of creating an interface condition. \see OPENCMISS::CMISSInterfaceConditionCreateStart
  SUBROUTINE INTERFACE_CONDITION_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_CONDITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has already been finished.",ERR,ERROR,*999)
      ELSE            
        !Finish the interface condition creation
        INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("INTERFACE_CONDITION_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_CONDITION_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an interface condition on an interface. \see OPENCMISS::CMISSInterfaceConditionCreateStart
  SUBROUTINE INTERFACE_CONDITION_CREATE_START(USER_NUMBER,INTERFACE,INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the interface condition
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the interface condition on
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<On return, a pointer to the interface condition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,interface_conditions_idx
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: NEW_INTERFACE_CONDITION
    TYPE(INTERFACE_CONDITION_PTR_TYPE), POINTER :: NEW_INTERFACE_CONDITIONS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    NULLIFY(NEW_INTERFACE_CONDITION)
    NULLIFY(NEW_INTERFACE_CONDITIONS)

    CALL ENTERS("INTERFACE_CONDITION_CREATE_START",ERR,ERROR,*997)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) THEN
        CALL INTERFACE_CONDITION_USER_NUMBER_FIND(USER_NUMBER,INTERFACE,NEW_INTERFACE_CONDITION,ERR,ERROR,*997)
        IF(ASSOCIATED(NEW_INTERFACE_CONDITION)) THEN
          LOCAL_ERROR="Interface condition user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
        ELSE
          NULLIFY(NEW_INTERFACE_CONDITION)
          !Initialise the new interface condition
          CALL INTERFACE_CONDITION_INITIALISE(NEW_INTERFACE_CONDITION,ERR,ERROR,*999)
          !Set default interface condition values
          NEW_INTERFACE_CONDITION%USER_NUMBER=USER_NUMBER
          NEW_INTERFACE_CONDITION%GLOBAL_NUMBER=INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS+1
          NEW_INTERFACE_CONDITION%INTERFACE_CONDITIONS=>INTERFACE%INTERFACE_CONDITIONS
          NEW_INTERFACE_CONDITION%INTERFACE=>INTERFACE
          !Default attributes
          NEW_INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD
          NEW_INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_TEST_OPERATOR
          !Add new interface condition into list of interface conditions in the interface
          ALLOCATE(NEW_INTERFACE_CONDITIONS(INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS+1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface conditions.",ERR,ERROR,*999)
          DO interface_conditions_idx=1,INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS
            NEW_INTERFACE_CONDITIONS(interface_conditions_idx)%PTR=>INTERFACE%INTERFACE_CONDITIONS% &
              & INTERFACE_CONDITIONS(interface_conditions_idx)%PTR
          ENDDO !interface_conditions_idx
          NEW_INTERFACE_CONDITIONS(INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS+1)%PTR=>NEW_INTERFACE_CONDITION
          IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) DEALLOCATE(INTERFACE%INTERFACE_CONDITIONS% &
            & INTERFACE_CONDITIONS)
          INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS=>NEW_INTERFACE_CONDITIONS
          INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=INTERFACE%INTERFACE_CONDITIONS% &
            NUMBER_OF_INTERFACE_CONDITIONS+1
          !Return the pointer
          INTERFACE_CONDITION=>NEW_INTERFACE_CONDITION  
        ENDIF
      ELSE
        LOCAL_ERROR="The interface conditions on interface number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*997)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_INTERFACE_CONDITION)) CALL INTERFACE_CONDITION_FINALISE(NEW_INTERFACE_CONDITION,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_INTERFACE_CONDITIONS)) DEALLOCATE(NEW_INTERFACE_CONDITIONS)
997 CALL ERRORS("INTERFACE_CONDITION_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_CREATE_START")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroys an interface condition. \see OPENCMISS::CMISSInterfaceConditionDestroy
  SUBROUTINE INTERFACE_CONDITION_DESTROY(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_condition_idx,interface_condition_position
    TYPE(INTERFACE_CONDITION_PTR_TYPE), POINTER :: NEW_INTERFACE_CONDITIONS(:)
    TYPE(INTERFACE_CONDITIONS_TYPE), POINTER :: INTERFACE_CONDITIONS

    NULLIFY(NEW_INTERFACE_CONDITIONS)

    CALL ENTERS("INTERFACE_CONDITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE_CONDITIONS=>INTERFACE_CONDITION%INTERFACE_CONDITIONS
      IF(ASSOCIATED(INTERFACE_CONDITIONS)) THEN
        interface_condition_position=INTERFACE_CONDITION%GLOBAL_NUMBER

        !Destroy all the interface condition components
        CALL INTERFACE_CONDITION_FINALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
        
        !Remove the interface condition from the list of interface conditions
        IF(INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS>1) THEN
          ALLOCATE(NEW_INTERFACE_CONDITIONS(INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new interface conditions.",ERR,ERROR,*999)
          DO interface_condition_idx=1,INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS
            IF(interface_condition_idx<interface_condition_position) THEN
              NEW_INTERFACE_CONDITIONS(interface_condition_idx)%PTR=>INTERFACE_CONDITIONS% &
                & INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            ELSE IF(interface_condition_idx>interface_condition_position) THEN
              INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interface_condition_idx)%PTR%GLOBAL_NUMBER=INTERFACE_CONDITIONS% &
                & INTERFACE_CONDITIONS(interface_condition_idx)%PTR%GLOBAL_NUMBER-1
              NEW_INTERFACE_CONDITIONS(interface_condition_idx-1)%PTR=>INTERFACE_CONDITIONS% &
                & INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            ENDIF
          ENDDO !interface_conditions_idx
          IF(ASSOCIATED(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) DEALLOCATE(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
          INTERFACE_CONDITIONS%INTERFACE_CONDITIONS=>NEW_INTERFACE_CONDITIONS
          INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS-1
        ELSE
          DEALLOCATE(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
          INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=0
        ENDIF
        
      ELSE
        CALL FLAG_ERROR("Interface conditions interface conditions is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*998)
    ENDIF    

    CALL EXITS("INTERFACE_CONDITIONS_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_INTERFACE_CONDITIONS)) DEALLOCATE(NEW_INTERFACE_CONDITIONS)
998 CALL ERRORS("INTERFACE_CONDITION_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_DESTROY")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of interface equations for the interface condition. \see OPENCMISS::CMISSInterfaceConditionEquationsCreateFinish
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish the creation of the interface equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_CONDITIONS_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      !Initialise the setup
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of interface equations for the interface condition. \see CMISSInterfaceConditionEquationsCreateStart
  !>Default values set for the INTERFACE_EQUATIONS's attributes are:
  !>- OUTPUT_TYPE: 0 (INTERFACE_EQUATIONS_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (INTERFACE_EQUATIONS_SPARSE_MATRICES)
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_START(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create the interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<On exit, a pointer to the created interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        CALL FLAG_ERROR("Interface equations is already associated.",ERR,ERROR,*999)
      ELSE
        !Initialise the setup
        !Return the pointer
        INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_CREATE_START")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the interface equations for an interface condition. \see OPENCMISS::CMISSInterfaceConditionEquationsDestroy
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_DESTROY(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface conditions to destroy the interface equations for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) THEN        
        CALL INTERFACE_EQUATIONS_DESTROY(INTERFACE_CONDITION%INTERFACE_EQUATIONS,ERR,ERROR,*999)
      ELSE
        CALL FLAG_ERROR("Interface condition interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_DESTROY")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Adds an equations set to an interface condition. \see OPENCMISS::CMISSInterfaceConditionEquationsSetAdd
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_SET_ADD(INTERFACE_CONDITION,MESH_INDEX,EQUATIONS_SET,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to add the equations set to
    INTEGER(INTG), INTENT(IN) :: MESH_INDEX !<The mesh index in the interface conditions interface that the equations set corresponds to
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to add.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    INTEGER(INTG), POINTER :: NEW_EQUATIONS_SETS_MESH_INDEXES(:)
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(EQUATIONS_SET_TYPE), POINTER :: INTERFACE_EQUATIONS_SET
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: DEPENDENT_MESH,INTERFACE_MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_EQUATIONS_SET_ADD",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      INTERFACE=>INTERFACE_CONDITION%INTERFACE
      IF(ASSOCIATED(INTERFACE)) THEN
        IF(MESH_INDEX>0.AND.MESH_INDEX<=INTERFACE%NUMBER_OF_COUPLED_MESHES) THEN
          IF(ASSOCIATED(EQUATIONS_SET)) THEN
            !Check that the equations set hasn't already been added.
            equations_set_idx=1
            NULLIFY(INTERFACE_EQUATIONS_SET)
            DO WHILE(equations_set_idx<=INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS.AND..NOT.ASSOCIATED(INTERFACE_EQUATIONS_SET))
              IF(ASSOCIATED(EQUATIONS_SET,INTERFACE_CONDITION%EQUATIONS_SETS(equations_set_idx)%PTR)) THEN
                INTERFACE_EQUATIONS_SET=>INTERFACE_CONDITION%EQUATIONS_SETS(equations_set_idx)%PTR
              ELSE
                equations_set_idx=equations_set_idx+1
              ENDIF
            ENDDO
            IF(ASSOCIATED(INTERFACE_EQUATIONS_SET)) THEN
              LOCAL_ERROR="The equations set has already been added to the interface condition at position index "// &
                & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              !Check the equations set and the mesh index match.
              INTERFACE_MESH=>INTERFACE%COUPLED_MESHES(MESH_INDEX)%PTR
              IF(ASSOCIATED(INTERFACE_MESH)) THEN
                IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                  DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    DECOMPOSITION=>DEPENDENT_FIELD%DECOMPOSITION
                    IF(ASSOCIATED(DECOMPOSITION)) THEN
                      DEPENDENT_MESH=>DECOMPOSITION%MESH
                      IF(ASSOCIATED(DEPENDENT_MESH)) THEN
                        IF(ASSOCIATED(INTERFACE_MESH,DEPENDENT_MESH)) THEN
                          !The meshes match to add the equations set to the list
                          ALLOCATE(NEW_EQUATIONS_SETS(INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets.",ERR,ERROR,*999)
                          ALLOCATE(NEW_EQUATIONS_SETS_MESH_INDEXES(INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new equations sets mesh indexes.",ERR,ERROR,*999)
                          DO equations_set_idx=1,INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS
                            NEW_EQUATIONS_SETS(equations_set_idx)%PTR=>INTERFACE_CONDITION%EQUATIONS_SETS(equations_set_idx)%PTR
                            NEW_EQUATIONS_SETS_MESH_INDEXES(equations_set_idx)= &
                              & INTERFACE_CONDITION%EQUATIONS_SETS_MESH_INDEXES(equations_set_idx)
                          ENDDO !equations_set_idx
                          NEW_EQUATIONS_SETS(INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS+1)%PTR=>EQUATIONS_SET
                          NEW_EQUATIONS_SETS_MESH_INDEXES(INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS+1)=MESH_INDEX
                          IF(ASSOCIATED(INTERFACE_CONDITION%EQUATIONS_SETS)) DEALLOCATE(INTERFACE_CONDITION%EQUATIONS_SETS)
                          IF(ASSOCIATED(INTERFACE_CONDITION%EQUATIONS_SETS_MESH_INDEXES)) &
                            & DEALLOCATE(INTERFACE_CONDITION%EQUATIONS_SETS_MESH_INDEXES)
                          INTERFACE_CONDITION%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
                          INTERFACE_CONDITION%EQUATIONS_SETS_MESH_INDEXES=>NEW_EQUATIONS_SETS_MESH_INDEXES
                          INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS=INTERFACE_CONDITION%NUMBER_OF_EQUATIONS_SETS+1
                        ELSE
                          CALL FLAG_ERROR("The dependent field mesh does not match the interface mesh.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("The dependent field decomposition mesh is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("The dependent field decomposition is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("The equation set dependent field is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The equation set dependent field has not been finished.",ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The interface mesh for mesh index "//TRIM(NUMBER_TO_VSTRING(MESH_INDEX,"*",ERR,ERROR))// &
                  & " is not associated."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            CALL FLAG_ERROR("Equations set is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specificed mesh index of "//TRIM(NUMBER_TO_VSTRING(MESH_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The mesh index must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE%NUMBER_OF_COUPLED_MESHES,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface condition interface is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF    

    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_SET_ADD")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_SET_ADD",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_SET_ADD")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_SET_ADD
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_FINALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      CALL INTERFACE_CONDITION_LAGRANGE_FINALISE(INTERFACE_CONDITION%LAGRANGE,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) &
        & CALL INTERFACE_EQUATIONS_DESTROY(INTERFACE_CONDITION%INTERFACE_EQUATIONS,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_CONDITION)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_FINALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for a finite element interface equations.
  SUBROUTINE INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE(INTERFACE_CONDITION,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: matrix_idx
    TYPE(ELEMENT_MATRIX_TYPE), POINTER :: ELEMENT_MATRIX
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE()")
#endif

    CALL ENTERS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      SELECT CASE(INTERFACE_CONDITION%METHOD)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE(INTERFACE_CONDITION_PENALTY_METHOD)
        CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Interface condition method "//TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
          & " is not valid."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
      INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
      IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
        IF(INTERFACE_EQUATIONS%OUTPUT_TYPE>=INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
          INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
          IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",ERR,ERROR,*999)          
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element number = ",ELEMENT_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",INTERFACE_MATRICES% &
                & NUMBER_OF_INTERFACE_MATRICES,ERR,ERROR,*999)
            DO matrix_idx=1,INTERFACE_MATRICES%NUMBER_OF_INTERFACE_MATRICES
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrix_idx,ERR,ERROR,*999)
              CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Update matrix = ",INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR% &
                & UPDATE_MATRIX,ERR,ERROR,*999)
              IF(INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR%UPDATE_MATRIX) THEN
                ELEMENT_MATRIX=>INTERFACE_MATRICES%MATRICES(matrix_idx)%PTR%ELEMENT_MATRIX
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
        ELSE
          CALL FLAG_ERROR("Interface matrices is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF    

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE()")
#endif
       
    CALL EXITS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition.
  SUBROUTINE INTERFACE_CONDITION_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      CALL FLAG_ERROR("Interface condition is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(INTERFACE_CONDITION,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition.",ERR,ERROR,*999)
      INTERFACE_CONDITION%USER_NUMBER=0
      INTERFACE_CONDITION%GLOBAL_NUMBER=0
      INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED=.FALSE.
      NULLIFY(INTERFACE_CONDITION%INTERFACE_CONDITIONS)
      NULLIFY(INTERFACE_CONDITION%INTERFACE)
      INTERFACE_CONDITION%METHOD=0
      INTERFACE_CONDITION%OPERATOR=0
      NULLIFY(INTERFACE_CONDITION%LAGRANGE)
      NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_FINALISE(INTERFACE_CONDITION,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an interface condition. \see OPENCMISS::CMISSInterfaceConditionLagrangeConditionCreateFinish
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to finish creating the Lagrange field for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has already been finished.",ERR,ERROR,*999)
      ELSE            
        !Finish the interface condition creation
        INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating the Lagrange multiplyer field for interface condition. \see OPENCMISS::CMISSInterfaceConditionLagrangeFieldCreateStart
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START(INTERFACE_CONDITION,LAGRANGE_FIELD_USER_NUMBER,LAGRANGE_FIELD, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create the Lagrange field on
    INTEGER(INTG), INTENT(IN) :: LAGRANGE_FIELD_USER_NUMBER !<The user specified Lagrange field number
    TYPE(FIELD_TYPE), POINTER :: LAGRANGE_FIELD !<If associated on entry, a pointer to the user created Lagrange field which has the same user number as the specified Lagrange field user number. If not associated on entry, on exit, a pointer to the created Lagrange field for the interface condition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: FIELD
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(REGION_TYPE), POINTER :: INTERFACE_REGION,LAGRANGE_FIELD_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        CALL FLAG_ERROR("Interface condition Lagrange is already associated.",ERR,ERROR,*999)
      ELSE        
        INTERFACE=>INTERFACE_CONDITION%INTERFACE
        IF(ASSOCIATED(INTERFACE)) THEN
          INTERFACE_REGION=>INTERFACE%PARENT_REGION
          IF(ASSOCIATED(INTERFACE_REGION)) THEN
            IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
              !Check the Lagrange field has been finished
              IF(LAGRANGE_FIELD%FIELD_FINISHED) THEN
                !Check the user numbers match
                IF(LAGRANGE_FIELD_USER_NUMBER/=LAGRANGE_FIELD%USER_NUMBER) THEN
                  LOCAL_ERROR="The specified Lagrange field user number of "// &
                    & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                    & " does not match the user number of the specified Lagrange field of "// &
                    & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD%USER_NUMBER,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
                LAGRANGE_FIELD_REGION=>LAGRANGE_FIELD%REGION
                IF(ASSOCIATED(LAGRANGE_FIELD_REGION)) THEN                
                  !Check the field is defined on the same region as the interface
                  IF(LAGRANGE_FIELD_REGION%USER_NUMBER/=INTERFACE_REGION%USER_NUMBER) THEN
                    LOCAL_ERROR="Invalid region setup. The specified Lagrange field has been created on interface number "// &
                      & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" in parent region number "// &
                      & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD_REGION%USER_NUMBER,"*",ERR,ERROR))// &
                      & " and the specified interface has been created in parent region number "// &
                      & TRIM(NUMBER_TO_VSTRING(INTERFACE_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                  
!!TODO setup Lagrange field
                  
                ELSE
                  CALL FLAG_ERROR("The Lagrange field region is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The specified Lagrange field has not been finished.",ERR,ERROR,*999)
              ENDIF
            ELSE
              !Check the user number has not already been used for a field in this region.
              NULLIFY(FIELD)
              CALL FIELD_USER_NUMBER_FIND(LAGRANGE_FIELD_USER_NUMBER,INTERFACE,FIELD,ERR,ERROR,*999)
              IF(ASSOCIATED(FIELD)) THEN
                LOCAL_ERROR="The specified Lagrange field user number of "// &
                  & TRIM(NUMBER_TO_VSTRING(LAGRANGE_FIELD_USER_NUMBER,"*",ERR,ERROR))// &
                  & " has already been used to create a field on interface number "// &
                  & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDIF
            CALL INTERFACE_CONDITION_LAGRANGE_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
            IF(.NOT.ASSOCIATED(LAGRANGE_FIELD)) INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED=.TRUE.
            !Set pointers
            IF(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED) THEN            
              LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
            ELSE
              INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD=>LAGRANGE_FIELD
            ENDIF
          ELSE
            CALL FLAG_ERROR("The interface parent region is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("The interface interface conditions is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface conditions is not associated.",ERR,ERROR,*999)
    ENDIF
     
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START")
    RETURN 1   
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FIELD_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Finalise the interface condition Lagrange information and deallocate all memory.
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FINALISE(INTERFACE_LAGRANGE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: INTERFACE_LAGRANGE !<A pointer to the interface condition Lagrange information to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_LAGRANGE)) THEN
       DEALLOCATE(INTERFACE_LAGRANGE)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface condition Lagrange information.
  SUBROUTINE INTERFACE_CONDITION_LAGRANGE_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<The pointer to the interface condition to initialise to initialise the Lagrange information for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_CONDITION_LAGRANGE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%LAGRANGE)) THEN
        CALL FLAG_ERROR("Interface condition Lagrange is already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE_CONDITION%LAGRANGE,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface condition Lagrange.",ERR,ERROR,*999)
        INTERFACE_CONDITION%LAGRANGE%INTERFACE_CONDITION=>INTERFACE_CONDITION
        INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FINISHED=.FALSE.
        INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITION_LAGRANGE_FINALISE(INTERFACE_CONDITION%LAGRANGE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITION_LAGRANGE_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_LAGRANGE_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_LAGRANGE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the interface condition method \see OPENCMISS::CMISSInterfaceConditionMethodGet
  SUBROUTINE INTERFACE_CONDITION_METHOD_GET(INTERFACE_CONDITION,INTERFACE_CONDITION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to get the method for
    INTEGER(INTG), INTENT(OUT) :: INTERFACE_CONDITION_METHOD !<On return, the interface condition method. \see INTERFACE_CONDITIONS_Methods,INTERFACE_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_METHOD_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        INTERFACE_CONDITION_METHOD=INTERFACE_CONDITION%METHOD
      ELSE
        CALL FLAG_ERROR("Interface condition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_METHOD_GET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_METHOD_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_METHOD_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_METHOD_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition method \see OPENCMISS::CMISSInterfaceConditionMethodSet
  SUBROUTINE INTERFACE_CONDITION_METHOD_SET(INTERFACE_CONDITION,INTERFACE_CONDITION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to set the method for
    INTEGER(INTG), INTENT(IN) :: INTERFACE_CONDITION_METHOD !<The interface condition method to set. \see INTERFACE_CONDITIONS_Methods,INTERFACE_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_METHOD_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(INTERFACE_CONDITION_METHOD)
        CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_POINT_TO_POINT_METHOD
        CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD
        CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD
         CASE(INTERFACE_CONDITION_PENALTY_METHOD)
          INTERFACE_CONDITION%METHOD=INTERFACE_CONDITION_PENALTY_METHOD
       CASE DEFAULT
          LOCAL_ERROR="The specified interface condition method of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION_METHOD,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_METHOD_SET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_METHOD_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_METHOD_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_METHOD_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the interface condition operator \see OPENCMISS::CMISSInterfaceConditionOperatorGet
  SUBROUTINE INTERFACE_CONDITION_OPERATOR_GET(INTERFACE_CONDITION,INTERFACE_CONDITION_OPERATOR,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to get the operator for
    INTEGER(INTG), INTENT(OUT) :: INTERFACE_CONDITION_OPERATOR !<On return, the interface condition operator. \see INTERFACE_CONDITIONS_Operators,INTERFACE_CONDITIONS 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_CONDITION_OPERATOR_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        INTERFACE_CONDITION_OPERATOR=INTERFACE_CONDITION%OPERATOR
      ELSE
        CALL FLAG_ERROR("Interface condition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_GET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_OPERATOR_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_OPERATOR_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the interface condition operator \see OPENCMISS::CMISSInterfaceConditionOperatorSet
  SUBROUTINE INTERFACE_CONDITION_OPERATOR_SET(INTERFACE_CONDITION,INTERFACE_CONDITION_OPERATOR,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to set the operator for
    INTEGER(INTG), INTENT(IN) :: INTERFACE_CONDITION_OPERATOR !<The interface condition operator to set. \see INTERFACE_CONDITIONS_Operators,INTERFACE_CONDITIONS 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_OPERATOR_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        CALL FLAG_ERROR("Interface condition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(INTERFACE_CONDITION_OPERATOR)
        CASE(INTERFACE_CONDITION_TEST_OPERATOR)
          INTERFACE_CONDITION%OPERATOR=INTERFACE_CONDITION_TEST_OPERATOR
        CASE DEFAULT
          LOCAL_ERROR="The specified interface condition operator of "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION_OPERATOR,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_SET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_OPERATOR_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_OPERATOR_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_OPERATOR_SET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in INTERFACE_CONDITION a pointer to the interface condition identified by USER_NUMBER in the given INTERFACE. If no interface condition with that USER_NUMBER exists INTERFACE_CONDITION is left nullified.
  SUBROUTINE INTERFACE_CONDITION_USER_NUMBER_FIND(USER_NUMBER,INTERFACE,INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find.
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<The interface to find the interface condition in.
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<On return a pointer to the interface condition with the given user number. If no interface condition with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: interface_condition_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("INTERFACE_CONDITION_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
        CALL FLAG_ERROR("Interface condition is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(INTERFACE_CONDITION)
        IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) THEN
          interface_condition_idx=1
          DO WHILE(interface_condition_idx<=INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS.AND. &
            & .NOT.ASSOCIATED(INTERFACE_CONDITION))
            IF(INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interface_condition_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
              INTERFACE_CONDITION=>INTERFACE%INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(interface_condition_idx)%PTR
            ELSE
              interface_condition_idx=interface_condition_idx+1
            ENDIF
          ENDDO
        ELSE
          LOCAL_ERROR="The interface conditions on interface number "// &
            & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITION_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises an interface conditions and deallocates all memory.
  SUBROUTINE INTERFACE_CONDITIONS_FINALISE(INTERFACE_CONDITIONS,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_CONDITIONS_TYPE), POINTER :: INTERFACE_CONDITIONS !<A pointer to the interface conditions to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    
    CALL ENTERS("INTERFACE_CONDITIONS_FINALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(INTERFACE_CONDITIONS)) THEN
      DO WHILE(INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS>0)
        INTERFACE_CONDITION=>INTERFACE_CONDITIONS%INTERFACE_CONDITIONS(1)%PTR
        CALL INTERFACE_CONDITION_DESTROY(INTERFACE_CONDITION,ERR,ERROR,*999)
      ENDDO
      IF(ASSOCIATED(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)) DEALLOCATE(INTERFACE_CONDITIONS%INTERFACE_CONDITIONS)
      DEALLOCATE(INTERFACE_CONDITIONS)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITIONS_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITIONS_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITIONS_FINALISE

  !
  !================================================================================================================================
  !
  
  !>Initialises an interface conditions for an interface.
  SUBROUTINE INTERFACE_CONDITIONS_INITIALISE(INTERFACE,ERR,ERROR,*) 

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the conditions for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
     
    CALL ENTERS("INTERFACE_CONDITIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%INTERFACE_CONDITIONS)) THEN
        LOCAL_ERROR="Interface conditions is already associated for interface number "// &
          & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE%INTERFACE_CONDITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface interface conditions.",ERR,ERROR,*999)
        INTERFACE%INTERFACE_CONDITIONS%INTERFACE=>INTERFACE
        INTERFACE%INTERFACE_CONDITIONS%NUMBER_OF_INTERFACE_CONDITIONS=0          
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("INTERFACE_CONDITIONS_INITIALISE")
    RETURN
999 CALL INTERFACE_CONDITIONS_FINALISE(INTERFACE%INTERFACE_CONDITIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_CONDITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_CONDITIONS_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE INTERFACE_CONDITIONS_ROUTINES
