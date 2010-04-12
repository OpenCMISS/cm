!> \file
!> $Id: interface_mapping_routines.f90 690 2009-09-30 23:27:16Z chrispbradley $
!> \author Chris Bradley
!> \brief This module contains all interface mapping routines.
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

!>This module contains all interface mapping routines.
MODULE INTERFACE_MAPPING_ROUTINES

  USE BASE_ROUTINES
  USE INPUT_OUTPUT
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  PUBLIC INTERFACE_MAPPING_CREATE_FINISH,INTERFACE_MAPPING_CREATE_START

  PUBLIC INTERFACE_MAPPING_DESTROY

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finishes the creation of interface mapping
  SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH(INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to finish the creation of.
    INTEGER(INTG), INTENT(OUT) ::       ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) ::    ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the process of creating the interface mapping for interface equations.
  SUBROUTINE INTERFACE_MAPPING_CREATE_START(INTERFACE_EQUATIONS,INTERFACE_MAPPING,ERR,ERROR,*)
    
    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to create the mapping for.
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<On exit, a pointer to the created interface mapping. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(INTERFACE_EQUATIONS%INTERFACE_EQUATIONS_FINISHED) THEN
        IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
          CALL FLAG_ERROR("Interface mapping is already associated.",ERR,ERROR,*999)
        ELSE
          NULLIFY(INTERFACE_MAPPING)
          CALL INTERFACE_MAPPING_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*999)
          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
        ENDIF
      ELSE
        CALL FLAG_ERROR("Interface equations have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface equations is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_START",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_CREATE_START")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_CREATE_START

  !
  !================================================================================================================================
  !

  !>Finalises an interface mapping create values cache and deallocates all memory
  SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE(CREATE_VALUES_CACHE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_CREATE_VALUES_CACHE_TYPE), POINTER :: CREATE_VALUES_CACHE !<A pointer to the create values cache to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(CREATE_VALUES_CACHE)) THEN
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COEFFICIENTS)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_INDEXES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_FIELD_INDEXES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COL_FIELD_INDEXES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COL_FIELD_INDEXES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_ROW_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_ROW_VARIABLE_TYPES)
      IF(ALLOCATED(CREATE_VALUES_CACHE%MATRIX_COL_VARIABLE_TYPES)) DEALLOCATE(CREATE_VALUES_CACHE%MATRIX_COL_VARIABLE_TYPES)
      DEALLOCATE(CREATE_VALUES_CACHE)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises an interface mapping create values cache 
  SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to create the create values cache for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,matrix_idx,VARIABLE_NUMBER
    TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      IF(ASSOCIATED(INTERFACE_MAPPING%CREATE_VALUES_CACHE)) THEN
        CALL FLAG_ERROR("Interface mapping create values cache is already associated.",ERR,ERROR,*998)
      ELSE
        INTERFACE_EQUATIONS=>INTERFACE_MAPPING%INTERFACE_EQUATIONS
        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
          INTERFACE_CONDITION=>INTERFACE_EQUATIONS%INTERFACE_CONDITION
          IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
            !Allocate and initialise the create values cache
            ALLOCATE(INTERFACE_MAPPING%CREATE_VALUES_CACHE,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface mapping create values cache.",ERR,ERROR,*999)
            INTERFACE_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES=0
            !Set the default interface mapping in the create values cache          
            !First calculate how many interface matrices we have and set the variable types
            SELECT CASE(INTERFACE_CONDITION%METHOD)
            CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
              INTERFACE_MAPPING%CREATE_VALUES_CACHE%NUMBER_OF_INTERFACE_MATRICES= &
                INTERFACE_CONDITION%NUMBER_OF_DEPENDENT_VARIABLES          
            CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)           
            CASE(INTERFACE_CONDITION_PENALTY_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The interface equations method of "// &
                & TRIM(NUMBER_TO_VSTRING(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE
            CALL FLAG_ERROR("Interface equations interface condition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Interface mapping interface equations is not associated.",ERR,ERROR,*998)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface mapping is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_CREATE_VALUES_CACHE_FINALISE(INTERFACE_MAPPING%CREATE_VALUES_CACHE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE")
    RETURN 1
  END SUBROUTINE INTERFACE_MAPPING_CREATE_VALUES_CACHE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Destroys an interface mapping.
  SUBROUTINE INTERFACE_MAPPING_DESTROY(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer the interface mapping to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      CALL INTERFACE_MAPPING_FINALISE(INTERFACE_MAPPING,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Equations mapping is not associated.",ERR,ERROR,*999)
    ENDIF
        
    CALL EXITS("INTERFACE_MAPPING_DESTROY")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_DESTROY",ERR,ERROR)    
    CALL EXITS("INTERFACE_MAPPING_DESTROY")
    RETURN 1
   
  END SUBROUTINE INTERFACE_MAPPING_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises the interface mapping and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_FINALISE(INTERFACE_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING !<A pointer to the interface mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("INTERFACE_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
      DEALLOCATE(INTERFACE_MAPPING)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("INTERFACE_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_FINALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the interface mapping and deallocates all memory.
  SUBROUTINE INTERFACE_MAPPING_INITIALISE(INTERFACE_EQUATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS !<A pointer to the interface equations to initialise the interface mapping for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("INTERFACE_MAPPING_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
      IF(ASSOCIATED(INTERFACE_EQUATIONS%INTERFACE_MAPPING)) THEN
        CALL FLAG_ERROR("Interface mapping is already associated.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(INTERFACE_EQUATIONS%INTERFACE_MAPPING,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface equations interface mapping.",ERR,ERROR,*999)
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_EQUATIONS=>INTERFACE_EQUATIONS
        INTERFACE_EQUATIONS%INTERFACE_MAPPING%INTERFACE_MAPPING_FINISHED=.FALSE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*998)
    ENDIF
       
    CALL EXITS("INTERFACE_MAPPING_INITIALISE")
    RETURN
999 CALL INTERFACE_MAPPING_FINALISE(INTERFACE_EQUATIONS%INTERFACE_MAPPING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("INTERFACE_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("INTERFACE_MAPPING_INITIALISE")
    RETURN 1
    
  END SUBROUTINE INTERFACE_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE INTERFACE_MAPPING_ROUTINES
