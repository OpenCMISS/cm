!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all region routines.
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

!> This module contains all region routines.
MODULE REGION_ROUTINES

  USE BASE_ROUTINES
  USE COORDINATE_ROUTINES
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  TYPE(REGION_TYPE), TARGET :: GLOBAL_REGION
  
  !Interfaces

  INTERFACE REGION_COORDINATE_SYSTEM_SET
    !!TODO: Allow for two more module procedurs to allow you to specify the coordinate system by number as well
    MODULE PROCEDURE REGION_COORDINATE_SYSTEM_SET_NUMBER
    MODULE PROCEDURE REGION_COORDINATE_SYSTEM_SET_PTR
  END INTERFACE !REGION_COORDINATE_SYSTEM_SET
  
  INTERFACE REGION_LABEL_SET
    MODULE PROCEDURE REGION_LABEL_SET_NUMBER
    MODULE PROCEDURE REGION_LABEL_SET_PTR
  END INTERFACE !REGION_LABEL_SET
  
  PUBLIC REGION_CREATE_START,REGION_CREATE_FINISH,REGION_SUB_REGION_CREATE_START,REGION_SUB_REGION_CREATE_FINISH,REGION_DESTROY, &
    & REGIONS_INITIALISE,REGIONS_FINALISE,REGION_USER_NUMBER_FIND,REGION_COORDINATE_SYSTEM_GET,REGION_LABEL_GET, &
    & REGION_COORDINATE_SYSTEM_SET,REGION_LABEL_SET

CONTAINS

  !
  !================================================================================================================================
  !

  FUNCTION REGION_COORDINATE_SYSTEM_GET(REGION,ERR,ERROR)

    !#### Function: REGION_COORDINATE_SYSTEM_GET
    !###  Type : TYPE(OORDINATE_SYSTEM_TYPE)
    !###  Description:
    !###    Returns the coordinate system of region REGION.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function result
    TYPE(COORDINATE_SYSTEM_TYPE) :: REGION_COORDINATE_SYSTEM_GET
    !Local Variables
    
    CALL ENTERS("REGION_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      REGION_COORDINATE_SYSTEM_GET=REGION%COORDINATE_SYSTEM
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_COORDINATE_SYSTEM_GET")
    RETURN
999 CALL ERRORS("REGION_COORDINATE_SYSTEM_GET",ERR,ERROR)
    CALL EXITS("REGION_COORDINATE_SYSTEM_GET")
    RETURN
  END FUNCTION REGION_COORDINATE_SYSTEM_GET
  
  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: REGION_COORDINATE_SYSTEM_SET
  !###  Description:
  !###    Sets the coordinate system of region.
  !###  Child-subroutines: REGION_COORDINATE_SYSTEM_SET_NUMBER,REGION_COORDINATE_SYSTEM_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_COORDINATE_SYSTEM_SET_NUMBER(USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*)

    !#### Subroutine: REGION_COORDINATE_SYSTEM_SET_NUMBER
    !###  Description:
    !###    Sets the coordinate system of region identified by USER_NUMBER to COORDINATE_SYSTEM. 

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: REGION

    CALL ENTERS("REGION_COORDINATE_SYSTEM_SET_NUMBER",ERR,ERROR,*999)

    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*999)
    CALL REGION_COORDINATE_SYSTEM_SET_PTR(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
    
    CALL EXITS("REGION_COORDINATE_SYSTEM_SET_NUMBER")
    RETURN
999 CALL ERRORS("REGION_COORDINATE_SYSTEM_SET_NUMBER",ERR,ERROR)
    CALL EXITS("REGION_COORDINATE_SYSTEM_SET_NUMBER")
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_SET_NUMBER
  
  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_COORDINATE_SYSTEM_SET_PTR(REGION,COORDINATE_SYSTEM,ERR,ERROR,*)

    !#### Subroutine: REGION_COORDINATE_SYSTEM_SET_PTR
    !###  Description:
    !###    Sets the coordinate system of region REGION identified by a pointer to COORDINATE_SYSTEM. 

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("REGION_COORDINATE_SYSTEM_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FLAG_ERROR("Region has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          REGION%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
        ELSE
          CALL FLAG_ERROR("Coordinate system is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_COORDINATE_SYSTEM_SET_PTR")
    RETURN
999 CALL ERRORS("REGION_COORDINATE_SYSTEM_SET_PTR",ERR,ERROR)
    CALL EXITS("REGION_COORDINATE_SYSTEM_SET_PTR")
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_SET_PTR
  
  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_CREATE_FINISH(REGION,ERR,ERROR,*)

    !#### Subroutine: REGION_CREATE_FINISH
    !###  Description:
    !###    Finishes the creation of a region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: sub_region_idx
    
    CALL ENTERS("REGION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      REGION%REGION_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of sub-regions in global region = ", &
        & GLOBAL_REGION%NUMBER_OF_SUB_REGIONS,ERR,ERROR,*999)    
      DO sub_region_idx=1,GLOBAL_REGION%NUMBER_OF_SUB_REGIONS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sub-region ",sub_region_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number = ",GLOBAL_REGION%SUB_REGIONS(sub_region_idx)%PTR%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Label = ",GLOBAL_REGION%SUB_REGIONS(sub_region_idx)%PTR%LABEL, &
          & ERR,ERROR,*999)
      ENDDO !sub_region_idx   
    ENDIF
    
    CALL EXITS("REGION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("REGION_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("REGION_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE REGION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_CREATE_START(USER_NUMBER,REGION,ERR,ERROR,*)

    !#### Subroutine: REGION_CREATE_START
    !###  Description:
    !###    Starts the creation of a new region number USER_NUMBER under the global region, initialises all variables and
    !###    inherits the global regions coordinate system etc.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: region_idx
    TYPE(REGION_TYPE), POINTER :: NEW_REGION
    TYPE(REGION_PTR_TYPE), POINTER :: NEW_SUB_REGIONS(:)
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_REGION)
    NULLIFY(NEW_SUB_REGIONS)
    
    CALL ENTERS("REGION_CREATE_START",ERR,ERROR,*999)

    NULLIFY(REGION)
    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN
      LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))//" has already been created"
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ELSE
      !Allocate the memory for the new region
      ALLOCATE(NEW_REGION,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new region",ERR,ERROR,*999)
      !Set default values
      NEW_REGION%USER_NUMBER=USER_NUMBER
      NEW_REGION%REGION_FINISHED=.FALSE.
      !CPB 21/02/07 The vstring operation crashes the AIX compiler so put a CHAR() around it.
      !NEW_REGION%LABEL="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
      LOCAL_STRING="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
      NEW_REGION%LABEL=CHAR(LOCAL_STRING)
      IF(ERR/=0) GOTO 999
      NULLIFY(NEW_REGION%COORDINATE_SYSTEM)
      NULLIFY(NEW_REGION%NODES)
      NULLIFY(NEW_REGION%MESHES)
      NULLIFY(NEW_REGION%FIELDS)
      NULLIFY(NEW_REGION%PARENT_REGION)
      NULLIFY(NEW_REGION%EQUATIONS_SETS)
      NEW_REGION%NUMBER_OF_SUB_REGIONS=0
      NULLIFY(NEW_REGION%SUB_REGIONS)
      CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(GLOBAL_REGION%COORDINATE_SYSTEM%USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
      IF(ASSOCIATED(GLOBAL_REGION%COORDINATE_SYSTEM)) THEN
        NEW_REGION%COORDINATE_SYSTEM=>GLOBAL_REGION%COORDINATE_SYSTEM
      ELSE
        CALL FLAG_ERROR("Global region does not have an associated coordinate system",ERR,ERROR,*999)
      ENDIF
      CALL NODES_INITIALISE(NEW_REGION,ERR,ERROR,*999)
      CALL MESHES_INITIALISE(NEW_REGION,ERR,ERROR,*999)
      CALL FIELDS_INITIALISE(NEW_REGION,ERR,ERROR,*999)
      CALL EQUATIONS_SETS_INITIALISE(NEW_REGION,ERR,ERROR,*999)
      NEW_REGION%NUMBER_OF_SUB_REGIONS=0
      NULLIFY(NEW_REGION%SUB_REGIONS)
      !Adjust the global region to include this new daughter region
      ALLOCATE(NEW_SUB_REGIONS(GLOBAL_REGION%NUMBER_OF_SUB_REGIONS+1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new_sub_regions",ERR,ERROR,*999)
      DO region_idx=1,GLOBAL_REGION%NUMBER_OF_SUB_REGIONS
        NEW_SUB_REGIONS(region_idx)%PTR=>GLOBAL_REGION%SUB_REGIONS(region_idx)%PTR
      ENDDO !region_idx
      GLOBAL_REGION%NUMBER_OF_SUB_REGIONS=GLOBAL_REGION%NUMBER_OF_SUB_REGIONS+1
      NEW_SUB_REGIONS(GLOBAL_REGION%NUMBER_OF_SUB_REGIONS)%PTR=>NEW_REGION
      IF(ASSOCIATED(GLOBAL_REGION%SUB_REGIONS)) DEALLOCATE(GLOBAL_REGION%SUB_REGIONS)
      GLOBAL_REGION%SUB_REGIONS=>NEW_SUB_REGIONS
      !Set the new regions parent region to the global region
      NEW_REGION%PARENT_REGION=>GLOBAL_REGION
      REGION=>NEW_REGION
    ENDIF
        
    CALL EXITS("REGION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_REGION)) DEALLOCATE(NEW_REGION)
    IF(ASSOCIATED(NEW_SUB_REGIONS)) DEALLOCATE(NEW_SUB_REGIONS)
    NULLIFY(REGION)
    CALL ERRORS("REGION_CREATE_START",ERR,ERROR)
    CALL EXITS("REGION_CREATE_START")
    RETURN 1
  END SUBROUTINE REGION_CREATE_START

  !
  !================================================================================================================================
  !

  RECURSIVE SUBROUTINE REGION_DESTROY(USER_NUMBER,ERR,ERROR,*)

    !#### Subroutine: REGION_DESTROY
    !###  Description:
    !###    Destroys a region given by USER_NUMBER and all sub-regions under it.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: count,nr
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(REGION_PTR_TYPE), POINTER :: NEW_SUB_REGIONS(:)

    CALL ENTERS("REGION_DESTROY",ERR,ERROR,*999)

    NULLIFY(REGION)
    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*999)
    IF(ASSOCIATED(REGION)) THEN

!!NOTE: We have to find a pointer to the region to destroy within this routine rather than passing in a pointer to a
!!DESTROY_REGION_PTR type routine because we need to change REGION%SUB_REGIONS of the PARENT region and this would violate section
!!12.4.1.6 of the Fortran standard if the dummy REGION pointer argument was associated with the SUB_REGIONS(x)%PTR actual
!!argument.
      
      IF(REGION%NUMBER_OF_SUB_REGIONS==0) THEN
        !No more daughter sub regions so delete this instance
        IF(ASSOCIATED(REGION%PARENT_REGION)) THEN
          NULLIFY(NEW_SUB_REGIONS)
          IF(REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS>1) THEN
            !If the parent region has more than one sub regions then remove this instance from its sub-regions list 
            ALLOCATE(NEW_SUB_REGIONS(REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-regions",ERR,ERROR,*999)
            count=0
            DO nr=1,REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS
              IF(REGION%PARENT_REGION%SUB_REGIONS(nr)%PTR%USER_NUMBER/=REGION%USER_NUMBER) THEN
                count=count+1
                NEW_SUB_REGIONS(count)%PTR=>REGION%PARENT_REGION%SUB_REGIONS(nr)%PTR
              ENDIF
            ENDDO !nr
            IF(ASSOCIATED(REGION%PARENT_REGION%SUB_REGIONS)) DEALLOCATE(REGION%PARENT_REGION%SUB_REGIONS)
          ENDIF
          REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS=REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS-1
          REGION%PARENT_REGION%SUB_REGIONS=>NEW_SUB_REGIONS
          !Destroy any data structures in this instance and deallocate any memory allocated.
          NULLIFY(REGION%COORDINATE_SYSTEM)
          CALL EQUATIONS_SETS_FINALISE(REGION,ERR,ERROR,*999)
          CALL FIELDS_FINALISE(REGION,ERR,ERROR,*999)
          CALL MESHES_FINALISE(REGION,ERR,ERROR,*999)
          CALL NODES_FINALISE(REGION,ERR,ERROR,*999)
          IF(ASSOCIATED(REGION%SUB_REGIONS)) DEALLOCATE(REGION%SUB_REGIONS)
          !Deallocate the current instance
          DEALLOCATE(REGION)
        ELSE
          CALL FLAG_ERROR("Parent region is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        !Recursively delete sub regions first
        DO WHILE(REGION%NUMBER_OF_SUB_REGIONS>0)
          CALL REGION_DESTROY(REGION%SUB_REGIONS(1)%PTR%USER_NUMBER,ERR,ERROR,*999)
        ENDDO
        !Now delete this instance
        CALL REGION_DESTROY(REGION%USER_NUMBER,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region number does not exist",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_DESTROY")
    RETURN
999 CALL ERRORS("REGION_DESTROY",ERR,ERROR)
    CALL EXITS("REGION_DESTROY")
    RETURN 1
  END SUBROUTINE REGION_DESTROY

  !
  !================================================================================================================================
  !

  FUNCTION REGION_LABEL_GET(REGION,ERR,ERROR)

    !#### Function: REGION_LABEL_GET
    !###  Type : TYPE(VARYING_STRING)
    !###  Description:
    !###    Returns the lable of region REGION.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Function result
    TYPE(VARYING_STRING) :: REGION_LABEL_GET
    !Local Variables

    CALL ENTERS("REGION_LABEL_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      !CPB 20/2/07 The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      REGION_LABEL_GET=VAR_STR(CHAR(REGION%LABEL))
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_LABEL_GET")
    RETURN
999 CALL ERRORS("REGION_LABEL_GET",ERR,ERROR)
    CALL EXITS("REGION_LABEL_GET")
    RETURN
  END FUNCTION REGION_LABEL_GET

  !
  !================================================================================================================================
  !
  
  !#### Generic-subroutine: REGION_LABEL_SET
  !###  Description:
  !###    Sets the label of region.
  !###  Child-subroutines: REGION_LABEL_SET_NUMBER,REGION_LABEL_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_LABEL_SET_NUMBER(USER_NUMBER,LABEL,ERR,ERROR,*)

    !#### Subroutine: REGION_LABEL_SET_NUMBER
    !###  Description:
    !###    Sets the lable of region identified by USER_NUMBER to LABEL 

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: REGION

    CALL ENTERS("REGION_LABEL_SET_NUMBER",ERR,ERROR,*999)

    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*999)
    CALL REGION_LABEL_SET_PTR(REGION,LABEL,ERR,ERROR,*999)
    
    CALL EXITS("REGION_LABEL_SET_NUMBER")
    RETURN
999 CALL ERRORS("REGION_LABEL_SET_NUMBER",ERR,ERROR)
    CALL EXITS("REGION_LABEL_SET_NUMBER")
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_LABEL_SET_PTR(REGION,LABEL,ERR,ERROR,*)

    !#### Subroutine: REGION_LABEL_SET_PTR
    !###  Description:
    !###    Sets the lable of region identified by the pointer REGION to LABEL 

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("REGION_LABEL_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FLAG_ERROR("Region has been finished",ERR,ERROR,*999)
      ELSE
        REGION%LABEL=LABEL
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_LABEL_SET_PTR")
    RETURN
999 CALL ERRORS("REGION_LABEL_SET_PTR",ERR,ERROR)
    CALL EXITS("REGION_LABEL_SET_PTR")
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_SUB_REGION_CREATE_START(USER_NUMBER,PARENT_REGION,SUB_REGION,ERR,ERROR,*)

    !#### Subroutine: REGION_SUB_REGION_CREATE_START
    !###  Description:
    !###    Starts the creation a new region number USER_NUMBER as a sub region to the given PARENT_REGION, initialises all
    !###    variables and inherits the PARENT_REGIONS coordinate system etc.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION
    TYPE(REGION_TYPE), POINTER :: SUB_REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: region_idx
    TYPE(REGION_TYPE), POINTER :: NEW_REGION
    TYPE(REGION_PTR_TYPE), POINTER :: NEW_SUB_REGIONS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_REGION)
    NULLIFY(NEW_SUB_REGIONS)
    
    CALL ENTERS("REGION_SUB_REGION_CREATE_START",ERR,ERROR,*999)

    NULLIFY(SUB_REGION)    
    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,SUB_REGION,ERR,ERROR,*999)
    IF(ASSOCIATED(SUB_REGION)) THEN
      LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has already been created"
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
    ELSE
      NULLIFY(NEW_REGION)
      IF(ASSOCIATED(PARENT_REGION)) THEN
        !Allocate the memory for the new region
        ALLOCATE(NEW_REGION,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new region",ERR,ERROR,*999)
        !Set default values
        NEW_REGION%USER_NUMBER=USER_NUMBER
        NEW_REGION%REGION_FINISHED=.FALSE.
        !CPB 21/02/07 The vstring operation crashes the AIX compiler so put a CHAR() etc. around it.
        !NEW_REGION%LABEL="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
        LOCAL_STRING="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
        NEW_REGION%LABEL=CHAR(LOCAL_STRING)
        IF(ERR/=0) GOTO 999
        NULLIFY(NEW_REGION%COORDINATE_SYSTEM)
        NULLIFY(NEW_REGION%NODES)
        NULLIFY(NEW_REGION%MESHES)
        NULLIFY(NEW_REGION%FIELDS)
        NULLIFY(NEW_REGION%PARENT_REGION)
        NEW_REGION%NUMBER_OF_SUB_REGIONS=0
        NULLIFY(NEW_REGION%SUB_REGIONS)
        IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN
          NEW_REGION%COORDINATE_SYSTEM=>PARENT_REGION%COORDINATE_SYSTEM
        ELSE
          CALL FLAG_ERROR("Parent region does not have an associated coordinate system",ERR,ERROR,*999)
        ENDIF
        CALL NODES_INITIALISE(NEW_REGION,ERR,ERROR,*999)
        CALL MESHES_INITIALISE(NEW_REGION,ERR,ERROR,*999)
        CALL FIELDS_INITIALISE(NEW_REGION,ERR,ERROR,*999)
        CALL EQUATIONS_SETS_INITIALISE(NEW_REGION,ERR,ERROR,*999)
        NEW_REGION%NUMBER_OF_SUB_REGIONS=0
        NULLIFY(NEW_REGION%SUB_REGIONS)
        !Adjust the parent region to include this new daughter
        ALLOCATE(NEW_SUB_REGIONS(PARENT_REGION%NUMBER_OF_SUB_REGIONS+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-regions",ERR,ERROR,*999)
        DO region_idx=1,PARENT_REGION%NUMBER_OF_SUB_REGIONS
          NEW_SUB_REGIONS(region_idx)%PTR=>PARENT_REGION%SUB_REGIONS(region_idx)%PTR
        ENDDO !region_no
        PARENT_REGION%NUMBER_OF_SUB_REGIONS=PARENT_REGION%NUMBER_OF_SUB_REGIONS+1
        NEW_SUB_REGIONS(PARENT_REGION%NUMBER_OF_SUB_REGIONS)%PTR=>NEW_REGION
        IF(ASSOCIATED(PARENT_REGION%SUB_REGIONS)) DEALLOCATE(PARENT_REGION%SUB_REGIONS)
        PARENT_REGION%SUB_REGIONS=>NEW_SUB_REGIONS
        !Set the new regions parent region to the parent region
        NEW_REGION%PARENT_REGION=>PARENT_REGION
        SUB_REGION=>NEW_REGION
      ELSE
        CALL FLAG_ERROR("Parent region is not associated",ERR,ERROR,*999)
      ENDIF
    ENDIF

    CALL EXITS("REGION_SUB_REGION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_REGION)) DEALLOCATE(NEW_REGION)
    IF(ASSOCIATED(NEW_SUB_REGIONS)) DEALLOCATE(NEW_SUB_REGIONS)
    NULLIFY(SUB_REGION)
    CALL ERRORS("REGION_SUB_REGION_CREATE_START",ERR,ERROR)
    CALL EXITS("REGION_SUB_REGION_CREATE_START")
    RETURN 1
  END SUBROUTINE REGION_SUB_REGION_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_SUB_REGION_CREATE_FINISH(REGION,ERR,ERROR,*)

    !#### Subroutine: REGION_SUB_REGION_CREATE_FINISH
    !###  Description:
    !###    Finishes the creation a new sub region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: sub_region_idx

    CALL ENTERS("REGION_SUB_REGION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      REGION%REGION_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Parent region number = ",REGION%PARENT_REGION%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of sub-regions in parent region = ", &
       & REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS,ERR,ERROR,*999)
      DO sub_region_idx=1,REGION%PARENT_REGION%NUMBER_OF_SUB_REGIONS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Sub-region ",sub_region_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number = ", &
          & REGION%PARENT_REGION%SUB_REGIONS(sub_region_idx)%PTR%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Label = ",REGION%PARENT_REGION%SUB_REGIONS(sub_region_idx)%PTR%LABEL, &
          & ERR,ERROR,*999)
      ENDDO !sub_region_idx
    ENDIF
    
    CALL EXITS("REGION_SUB_REGION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("REGION_SUB_REGION_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("REGION_SUB_REGION_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE REGION_SUB_REGION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  SUBROUTINE REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*)

    !#### Subroutine: REGION_USER_NUMBER_FIND
    !###  Description:
    !###    Finds and returns in REGION a pointer to the region with the number given in USER_NUMBER. If no region with that number
    !###    exits REGION is left nullified.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: nr

    CALL ENTERS("REGION_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(USER_NUMBER==0) THEN
      REGION=>GLOBAL_REGION
    ELSE
      NULLIFY(REGION)      
      nr=1
      DO WHILE(nr<=GLOBAL_REGION%NUMBER_OF_SUB_REGIONS.AND..NOT.ASSOCIATED(REGION))
        CALL REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,REGION,GLOBAL_REGION%SUB_REGIONS(nr)%PTR,ERR,ERROR,*999)
        IF(.NOT.ASSOCIATED(REGION)) nr=nr+1        
      END DO
    ENDIF
  
    CALL EXITS("REGION_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("REGION_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("REGION_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE REGION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  RECURSIVE SUBROUTINE REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,REGION,START_REGION,ERR,ERROR,*)

    !#### Subroutine: REGION_USER_NUMBER_FIND_PTR
    !###  Description:
    !###    Finds and returns in REGION a pointer to the region with the number given in USER_NUMBER starting from the
    !###    START_REGION and searching all sub-regions under the START_REGION. If no region with that number exit REGION is
    !###    left nullified.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION,START_REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: nr

    CALL ENTERS("REGION_USER_NUMBER_FIND_PTR",ERR,ERROR,*999)

    NULLIFY(REGION)
    IF(ASSOCIATED(START_REGION)) THEN
      IF(START_REGION%USER_NUMBER==USER_NUMBER) THEN
        REGION=>START_REGION
      ELSE
        nr=1
        DO WHILE(nr<=START_REGION%NUMBER_OF_SUB_REGIONS.AND..NOT.ASSOCIATED(REGION))
          CALL REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,REGION,START_REGION%SUB_REGIONS(nr)%PTR,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(REGION)) nr=nr+1
        END DO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Start region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_USER_NUMBER_FIND_PTR")
    RETURN
999 CALL ERRORS("REGION_USER_NUMBER_FIND_PTR",ERR,ERROR)
    CALL EXITS("REGION_USER_NUMBER_FIND_PTR")
    RETURN 1
  END SUBROUTINE REGION_USER_NUMBER_FIND_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE REGIONS_INITIALISE(ERR,ERROR,*)

    !#### Subroutine: REGIONS_INITIALISE
    !###  Description:
    !###    Initialises the regions and creates the global world region

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM

    CALL ENTERS("REGIONS_INITIALISE",ERR,ERROR,*999)

    GLOBAL_REGION%USER_NUMBER=0
    GLOBAL_REGION%LABEL="Global (World) Region"
    CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(0,COORDINATE_SYSTEM,ERR,ERROR,*999)
    IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
      GLOBAL_REGION%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
    ELSE
      CALL FLAG_ERROR("Could not find global coordinate system - Initialise coordinate systems first",ERR,ERROR,*999)
    ENDIF
    NULLIFY(GLOBAL_REGION%NODES)
    NULLIFY(GLOBAL_REGION%FIELDS)
    NULLIFY(GLOBAL_REGION%MESHES)
    GLOBAL_REGION%NUMBER_OF_SUB_REGIONS=0
    NULLIFY(GLOBAL_REGION%SUB_REGIONS)
    NULLIFY(GLOBAL_REGION%PARENT_REGION)
   
    CALL EXITS("REGIONS_INITIALISE")
    RETURN
999 CALL ERRORS("REGIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("REGIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE REGIONS_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE REGIONS_FINALISE(ERR,ERROR,*)

    !#### Subroutine: REGIONS_FINALISE
    !###  Description:
    !###    Finalises the regions and destroys any current regions.

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: nr

    CALL ENTERS("REGIONS_FINALISE",ERR,ERROR,*999)

    !Destroy any global region daughter regions first
    DO nr=1,GLOBAL_REGION%NUMBER_OF_SUB_REGIONS
      CALL REGION_DESTROY(GLOBAL_REGION%SUB_REGIONS(nr)%PTR%USER_NUMBER,ERR,ERROR,*999)
    ENDDO !region
    !Destroy global region and deallocated any memory allocated in the global region
    GLOBAL_REGION%LABEL=""
    NULLIFY(GLOBAL_REGION%COORDINATE_SYSTEM)
    GLOBAL_REGION%NUMBER_OF_SUB_REGIONS=0
    IF(ASSOCIATED(GLOBAL_REGION%SUB_REGIONS)) DEALLOCATE(GLOBAL_REGION%SUB_REGIONS)
   
    CALL EXITS("REGIONS_FINALISE")
    RETURN
999 CALL ERRORS("REGIONS_FINALISE",ERR,ERROR)
    CALL EXITS("REGIONS_FINALISE")
    RETURN 1
  END SUBROUTINE REGIONS_FINALISE
  
  !
  !================================================================================================================================
  !
  
END MODULE REGION_ROUTINES
