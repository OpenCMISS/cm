!> \file
!> $Id: interface_equations_routines.f90 818 2009-12-13 21:49:58Z chrispbradley $
!> \author Chris Bradley
!> \brief This module handles all interface equations routines.
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

!>This module handles all interface equations routines.
MODULE INTERFACE_EQUATIONS_ROUTINES

  USE BASE_ROUTINES
  USE EQUATIONS_ROUTINES
  USE INTERFACE_MAPPING_ROUTINES
  USE INTERFACE_MATRICES_ROUTINES
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup INTERFACE_EQUATIONS_ROUTINES_OutputTypes INTERFACE_EQUATIONS_ROUTINES::OutputTypes
  !> \brief The interface equations output types
  !> \see INTERFACE_EQUATIONS_ROUTINES,OPENCMISS_InterfaceEquationsConstants
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_NO_OUTPUT=0 !<No output. \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_TIMING_OUTPUT=1 !<Timing information output. \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output. \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output .\see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES 
  !>@}

  !> \addtogroup INTERFACE_EQUATIONS_ROUTINES_SparsityTypes INTERFACE_EQUATIONS_ROUTINES::SparsityTypes
  !> \brief Interface equations matrices sparsity types
  !> \see INTERFACE_EQUATIONS_ROUTINES,OPENCMISS_InterfaceEquationsSparsityTypes
  !>@{
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the interface equations. \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
  INTEGER(INTG), PARAMETER :: INTERFACE_EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the interface equations. \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
 !>@}

  !Module types

  !Module variables

  !Interfaces

  PUBLIC INTERFACE_EQUATIONS_NO_OUTPUT,INTERFACE_EQUATIONS_TIMING_OUTPUT,INTERFACE_EQUATIONS_MATRIX_OUTPUT, &
    & INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT

  PUBLIC INTERFACE_EQUATIONS_SPARSE_MATRICES,INTERFACE_EQUATIONS_FULL_MATRICES

  PUBLIC INTERFACE_EQUATIONS_CREATE_FINISH,INTERFACE_EQUATIONS_CREATE_START

  PUBLIC INTERFACE_EQUATIONS_DESTROY

  PUBLIC INTERFACE_EQUATIONS_OUTPUT_TYPE_GET,INTERFACE_EQUATIONS_OUTPUT_TYPE_SET
  
  PUBLIC INTERFACE_EQUATIONS_SPARSITY_TYPE_GET,INTERFACE_EQUATIONS_SPARSITY_TYPE_SET

  PUBLIC INTERFACE_CONDITION_EQUATIONS_GET
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_CREATE_FINISH(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_EQUATIONS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Interface equations have already been finished.",ERR,ERROR,*999)        
      ELSE
        !Set the finished flag
        INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_CREATE_FINISH")
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of interface equations for an interface condition.
  SUBROUTINE INTERFACE_EQUATIONS_CREATE_START(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to create interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<On exit, a pointer to the created interface equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_EQUATIONS_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) THEN
        CALL FLAG_ERROR("Interface equations are already associated for the interface condition.",ERR,ERROR,*999)        
      ELSE
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          CALL FLAG_ERROR("Interface equations is already associated.",ERR,ERROR,*999)
        ELSE
          !Initialise the equations
          CALL INTERFACE_EQUATIONS_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*999)
          !Return the pointer
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_CREATE_START",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_CREATE_START")
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_DESTROY(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      CALL INTERFACE_EQUATIONS_FINALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_DESTROY",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_DESTROY")
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the interface equations and deallocate all memory.
  SUBROUTINE INTERFACE_EQUATIONS_FINALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_EQUATIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MAPPING)) &
        & CALL INTERFACE_MAPPING_DESTROY(INTERFACE_EQUATIONS%INTERFACE_MAPPING,ERR,ERROR,*999)
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MATRICES)) &
        & CALL INTERFACE_MATRICES_DESTROY(INTERFACE_EQUATIONS%INTERFACE_MATRICES,ERR,ERROR,*999)
      DEALLOCATE(INTERFACE_EQUATIONS)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface equations for an interface condition
  SUBROUTINE INTERFACE_EQUATIONS_INITIALISE(INTERFACE_CONDITION,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface condition to initialise the interface equations for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    CALL ENTERS("INTERFACE_EQUATIONS_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(ASSOCIATED(INTERFACE_CONDITION%INTERFACE_EQUATIONS)) THEN
        CALL FLAG_ERROR("Interface equations is already associated for this interface condition.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_CONDITION%INTERFACE_EQUATIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface equations.",ERR,ERROR,*999)
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_CONDITION=>INTERFACE_CONDITION
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_NO_OUTPUT
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%SPARSITY_TYPE=INTERFACE_EQUATIONS_SPARSE_MATRICES
        NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_MAPPING)
        NULLIFY(INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_MATRICES)
        INTERFACE_CONDITION%INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED=.FALSE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface condition is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_INITIALISE")
    RETURN
999 CALL INTERFACE_EQUATIONS_FINALISE(INTERFACE_CONDITION%INTERFACE_EQUATIONS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_EQUATIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_INITIALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Gets the output type for interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_GET(INTERFACE_EQUATIONS,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to get the output type for
    INTEGER(INTG), INTENT(OUT) :: OUTPUT_TYPE !<On exit, the output type of the interface equations \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        OUTPUT_TYPE=INTERFACE_EQUATIONS%OUTPUT_TYPE
      ELSE
        CALL FLAG_ERROR("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for the interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_SET(INTERFACE_EQUATIONS,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to set the output type for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The output type to set \see INTERFACE_EQUATIONS_ROUTINES_OutputTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Interface equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(OUTPUT_TYPE)
        CASE(INTERFACE_EQUATIONS_NO_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_NO_OUTPUT
        CASE(INTERFACE_EQUATIONS_TIMING_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_TIMING_OUTPUT
        CASE(INTERFACE_EQUATIONS_MATRIX_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_MATRIX_OUTPUT
        CASE(INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT)
          INTERFACE_EQUATIONS%OUTPUT_TYPE=INTERFACE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
        CASE DEFAULT
          LOCAL_ERROR="The specified output type of "//TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_OUTPUT_TYPE_SET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_EQUATIONS_OUTPUT_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Gets the sparsity type for interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_GET(INTERFACE_EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to get the sparsity type for
    INTEGER(INTG), INTENT(OUT) :: SPARSITY_TYPE !<On exit, the sparsity type of the interface equations. \see INTERFACE_EQUATIONS_ROUTINES_SparsityTypes,INTERFACE_EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        SPARSITY_TYPE=INTERFACE_EQUATIONS%SPARSITY_TYPE
      ELSE
        CALL FLAG_ERROR("Interface equations has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_GET")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for the interface equations.
  SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_SET(INTERFACE_EQUATIONS,SPARSITY_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: SPARSITY_TYPE !<The sparsity type to set \see EQUATIONS_ROUTINES_SparsityTypes,EQUATIONS_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        CALL FLAG_ERROR("Interface equations has already been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(SPARSITY_TYPE)
        CASE(INTERFACE_EQUATIONS_SPARSE_MATRICES)
          INTERFACE_EQUATIONS%SPARSITY_TYPE=INTERFACE_EQUATIONS_SPARSE_MATRICES
        CASE(INTERFACE_EQUATIONS_FULL_MATRICES)
          INTERFACE_EQUATIONS%SPARSITY_TYPE=INTERFACE_EQUATIONS_FULL_MATRICES
        CASE DEFAULT
          LOCAL_ERROR="The specified sparsity type of "//TRIM(NUMBER_TO_VSTRING(SPARSITY_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN
999 CALL ERRORS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET",ERR,ERROR)
    CALL EXITS("INTERFACE_EQUATIONS_SPARSITY_TYPE_SET")
    RETURN 1
  END SUBROUTINE INTERFACE_EQUATIONS_SPARSITY_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the interface equations for an interface conditions.
  SUBROUTINE INTERFACE_CONDITION_EQUATIONS_GET(INTERFACE_CONDITION,INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION !<A pointer to the interface conditions to get the interface equations for
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<On exit, a pointer to the interface equations in the specified interface condition. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    CALL ENTERS("INTERFACE_CONDITION_EQUATIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
      IF(INTERFACE_CONDITION%INTERFACE_CONDITION_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          CALL FLAG_ERROR("Interface equations is already associated.",ERR,ERROR,*999)
        ELSE
          INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
          IF(.NOT.ASSOCIATED(INTERFACE_EQUATIONS)) &
            & CALL FLAG_ERROR("Interface equations set equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations set has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_GET")
    RETURN
999 CALL ERRORS("INTERFACE_CONDITION_EQUATIONS_GET",ERR,ERROR)
    CALL EXITS("INTERFACE_CONDITION_EQUATIONS_GET")
    RETURN 1
    
  END SUBROUTINE INTERFACE_CONDITION_EQUATIONS_GET

  !
  !================================================================================================================================
  !
  
END MODULE INTERFACE_EQUATIONS_ROUTINES
