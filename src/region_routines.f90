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
  USE NODE_ROUTINES
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  TYPE(REGIONS_TYPE) :: REGIONS
  
  !Interfaces
  
  PUBLIC REGION_CREATE_START,REGION_CREATE_FINISH,REGION_DESTROY,REGIONS_INITIALISE,REGIONS_FINALISE,REGION_USER_NUMBER_FIND, &
    & REGION_COORDINATE_SYSTEM_GET,REGION_LABEL_GET,REGION_COORDINATE_SYSTEM_SET,REGION_LABEL_SET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns the coordinate system of region. \see OPENCMISS::CMISSRegionCoordinateSystemGet
  SUBROUTINE REGION_COORDINATE_SYSTEM_GET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<On exit, the coordinate system for the specified region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("REGION_COORDINATE_SYSTEM_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          CALL FLAG_ERROR("Coordinate system is already associated.",ERR,ERROR,*999)
        ELSE
          COORDINATE_SYSTEM=>REGION%COORDINATE_SYSTEM
        ENDIF
      ELSE
        CALL FLAG_ERROR("Region has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_COORDINATE_SYSTEM_GET")
    RETURN
999 CALL ERRORS("REGION_COORDINATE_SYSTEM_GET",ERR,ERROR)
    CALL EXITS("REGION_COORDINATE_SYSTEM_GET")
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_GET
  
  !
  !================================================================================================================================
  !

  !>Sets the coordinate system of region.  \see OPENCMISS::CMISSRegionCoordinateSystemSet
  SUBROUTINE REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to set the coordinate system for
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM !<The coordinate system to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("REGION_COORDINATE_SYSTEM_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FLAG_ERROR("Region has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(COORDINATE_SYSTEM)) THEN
          IF(COORDINATE_SYSTEM%COORDINATE_SYSTEM_FINISHED) THEN
            REGION%COORDINATE_SYSTEM=>COORDINATE_SYSTEM
          ELSE
            CALL FLAG_ERROR("Coordinate system has not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Coordinate system is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_COORDINATE_SYSTEM_SET")
    RETURN
999 CALL ERRORS("REGION_COORDINATE_SYSTEM_SET",ERR,ERROR)
    CALL EXITS("REGION_COORDINATE_SYSTEM_SET")
    RETURN 1
  END SUBROUTINE REGION_COORDINATE_SYSTEM_SET
  
  !
  !================================================================================================================================
  !

  !>Finishes the creation of a region. \see OPENCMISS::CMISSRegionCreateFinish
  SUBROUTINE REGION_CREATE_FINISH(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
     
    CALL ENTERS("REGION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FLAG_ERROR("Region has already been finished.",ERR,ERROR,*999)
      ELSE
        REGION%REGION_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Region : ",REGION%USER_NUMBER,ERR,ERROR,*999)      
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",REGION%LABEL,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%PARENT_REGION)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region user number = ",REGION%PARENT_REGION%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Parent region label = ",REGION%PARENT_REGION%LABEL, &
          & ERR,ERROR,*999)        
      ENDIF
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
  
  !>Starts the creation a new region number USER_NUMBER as a sub region to the given PARENT_REGION, initialises all
  !>variables and inherits the PARENT_REGIONS coordinate system. \see OPENCMISS::CMISSRegionCreateFinish
  !>Default values set for the REGION's attributes are:
  !>- COORDINATE_SYSTEM: parent coordinate system. See \ref COORDINATE_SYSTEM_TYPE
  !>- NODES: null
  !>- MESHES: 0 mesh
  !>- FIELDS: 0 field
  !>- EQUATIONS_SETS: 0 equation set
  !>- PARENT_REGION: global region
  !>- NUMBER_OF_SUB_REGIONS: 0
  !>- SUB_REGIONS: 0 region
  SUBROUTINE REGION_CREATE_START(USER_NUMBER,PARENT_REGION,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to create
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION !<A pointer to the parent region
    TYPE(REGION_TYPE), POINTER :: REGION !<On exit, a pointer to the created region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,region_idx
    TYPE(REGION_TYPE), POINTER :: NEW_REGION
    TYPE(REGION_PTR_TYPE), POINTER :: NEW_SUB_REGIONS(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR,LOCAL_STRING

    NULLIFY(NEW_REGION)
    NULLIFY(NEW_SUB_REGIONS)
    
    CALL ENTERS("REGION_CREATE_START",ERR,ERROR,*997)

    CALL REGION_USER_NUMBER_FIND(USER_NUMBER,NEW_REGION,ERR,ERROR,*997)
    IF(ASSOCIATED(NEW_REGION)) THEN
      LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
        & " has already been created."
      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*997)
    ELSE
      IF(ASSOCIATED(REGION)) THEN
        CALL FLAG_ERROR("Region is already associated.",ERR,ERROR,*997)
      ELSE
        NULLIFY(REGION)
        IF(ASSOCIATED(PARENT_REGION)) THEN
          IF(PARENT_REGION%REGION_FINISHED) THEN
            IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN
              !Initialise the region
              CALL REGION_INITIALISE(REGION,ERR,ERROR,*999)
              !Set the user number
              REGION%USER_NUMBER=USER_NUMBER
              !CPB 21/02/07 The vstring operation crashes the AIX compiler so put a CHAR() etc. around it.
              !REGION%LABEL="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
              LOCAL_STRING="Region "//NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR)
              REGION%LABEL=CHAR(LOCAL_STRING)
              IF(ERR/=0) GOTO 999
              REGION%COORDINATE_SYSTEM=>PARENT_REGION%COORDINATE_SYSTEM
              !Adjust the parent region to include this new daughter
              ALLOCATE(NEW_SUB_REGIONS(PARENT_REGION%NUMBER_OF_SUB_REGIONS+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-regions.",ERR,ERROR,*999)
              DO region_idx=1,PARENT_REGION%NUMBER_OF_SUB_REGIONS
                NEW_SUB_REGIONS(region_idx)%PTR=>PARENT_REGION%SUB_REGIONS(region_idx)%PTR
              ENDDO !region_no
              PARENT_REGION%NUMBER_OF_SUB_REGIONS=PARENT_REGION%NUMBER_OF_SUB_REGIONS+1
              NEW_SUB_REGIONS(PARENT_REGION%NUMBER_OF_SUB_REGIONS)%PTR=>REGION
              IF(ASSOCIATED(PARENT_REGION%SUB_REGIONS)) DEALLOCATE(PARENT_REGION%SUB_REGIONS)
              PARENT_REGION%SUB_REGIONS=>NEW_SUB_REGIONS
              !Set the new regions parent region to the parent region
              REGION%PARENT_REGION=>PARENT_REGION
            ELSE
              CALL FLAG_ERROR("Parent region does not have an associated coordinate system.",ERR,ERROR,*997)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Parent region has not been finished.",ERR,ERROR,*997)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Parent region is not associated.",ERR,ERROR,*997)
        ENDIF
      ENDIF
    ENDIF
    
    CALL EXITS("REGION_CREATE_START")
    RETURN
999 CALL REGION_FINALISE(REGION,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_SUB_REGIONS)) DEALLOCATE(NEW_SUB_REGIONS)
997 CALL ERRORS("REGION_CREATE_START",ERR,ERROR)
    CALL EXITS("REGION_CREATE_START")
    RETURN 1
  END SUBROUTINE REGION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a region given by USER_NUMBER and all sub-regions under it. \todo create destroy by pointer method. \see OPENCMISS::CMISSRegionDestroy
  RECURSIVE SUBROUTINE REGION_DESTROY(USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
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
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new sub-regions.",ERR,ERROR,*999)
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
          !Finalise the region
          CALL REGION_FINALISE(REGION,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Parent region is not associated.",ERR,ERROR,*999)
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
      CALL FLAG_ERROR("Region number does not exist.",ERR,ERROR,*999)
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

  !>Finalises a region and deallocates all memory
  SUBROUTINE REGION_FINALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
   
     CALL ENTERS("REGION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL EQUATIONS_SETS_FINALISE(REGION,ERR,ERROR,*999)
      CALL FIELDS_FINALISE(REGION,ERR,ERROR,*999)
      CALL MESHES_FINALISE(REGION,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%NODES)) CALL NODES_DESTROY(REGION%NODES,ERR,ERROR,*999)
      IF(ASSOCIATED(REGION%SUB_REGIONS)) DEALLOCATE(REGION%SUB_REGIONS)
      DEALLOCATE(REGION)
    ENDIF
    
    CALL EXITS("REGION_FINALISE")
    RETURN
999 CALL ERRORS("REGION_FINALISE",ERR,ERROR)
    CALL EXITS("REGION_FINALISE")
    RETURN 1
  END SUBROUTINE REGION_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a region.
  SUBROUTINE REGION_INITIALISE(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
   
    CALL ENTERS("REGION_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      CALL FLAG_ERROR("Region is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(REGION,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region.",ERR,ERROR,*999)
      REGION%USER_NUMBER=0
      REGION%REGION_FINISHED=.FALSE.
      NULLIFY(REGION%COORDINATE_SYSTEM)
      NULLIFY(REGION%NODES)
      NULLIFY(REGION%MESHES)
      NULLIFY(REGION%FIELDS)
      NULLIFY(REGION%EQUATIONS_SETS)
      NULLIFY(REGION%PARENT_REGION)
      REGION%NUMBER_OF_SUB_REGIONS=0
      NULLIFY(REGION%SUB_REGIONS)
      CALL MESHES_INITIALISE(REGION,ERR,ERROR,*999)
      CALL FIELDS_INITIALISE(REGION,ERR,ERROR,*999)
      CALL EQUATIONS_SETS_INITIALISE(REGION,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_INITIALISE")
    RETURN
999 CALL REGION_FINALISE(REGION,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("REGION_INITIALISE",ERR,ERROR)
    CALL EXITS("REGION_INITIALISE")
    RETURN 1
  END SUBROUTINE REGION_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the label of a region . \see OPENCMISS::CMISSRegionLabelGet
  FUNCTION REGION_LABEL_GET(REGION,ERR,ERROR)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to get the label for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function result
    TYPE(VARYING_STRING) :: REGION_LABEL_GET
    !Local Variables

    CALL ENTERS("REGION_LABEL_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      !CPB 20/2/07 The following line crashes the AIX compiler unless it has a VAR_STR(CHAR()) around it
      REGION_LABEL_GET=VAR_STR(CHAR(REGION%LABEL))
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
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

  !>Sets the label of a region.\todo need a varying string method \see OPENCMISS::CMISSRegionLabelSet
  SUBROUTINE REGION_LABEL_SET(REGION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("REGION_LABEL_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(REGION%REGION_FINISHED) THEN
        CALL FLAG_ERROR("Region has been finished.",ERR,ERROR,*999)
      ELSE
        REGION%LABEL=LABEL
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("REGION_LABEL_SET")
    RETURN
999 CALL ERRORS("REGION_LABEL_SET",ERR,ERROR)
    CALL EXITS("REGION_LABEL_SET")
    RETURN 1
  END SUBROUTINE REGION_LABEL_SET

  !
  !================================================================================================================================
  !

  !>Finds and returns in REGION a pointer to the region with the number given in USER_NUMBER. If no region with that number
  !>exits REGION is left nullified.
  SUBROUTINE REGION_USER_NUMBER_FIND(USER_NUMBER,REGION,ERR,ERROR,*)

     !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the region to find
    TYPE(REGION_TYPE), POINTER :: REGION !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nr
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION
    
    CALL ENTERS("REGION_USER_NUMBER_FIND",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL FLAG_ERROR("Region is already associated.",ERR,ERROR,*999)
    ELSE
      NULLIFY(REGION)
      WORLD_REGION=>REGIONS%WORLD_REGION
      IF(ASSOCIATED(WORLD_REGION)) THEN
        nr=1
        DO WHILE(nr<=WORLD_REGION%NUMBER_OF_SUB_REGIONS.AND..NOT.ASSOCIATED(REGION))
          CALL REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,WORLD_REGION%SUB_REGIONS(nr)%PTR,REGION,ERR,ERROR,*999)
          IF(.NOT.ASSOCIATED(REGION)) nr=nr+1        
        END DO
      ELSE
        CALL FLAG_ERROR("World region is not associated.",ERR,ERROR,*999)
      ENDIF
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

  !>Finds and returns in REGION a pointer to the region with the number given in USER_NUMBER starting from the
  !>START_REGION and searching all sub-regions under the START_REGION. If no region with that number exit REGION is
  !>left nullified.
  RECURSIVE SUBROUTINE REGION_USER_NUMBER_FIND_PTR(USER_NUMBER,START_REGION,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to find
    TYPE(REGION_TYPE), POINTER :: START_REGION !<A pointer to the region to start the search from
    TYPE(REGION_TYPE), POINTER :: REGION !<On exit, a pointer to the region with the specified user number if it exists. If no region exists with the specified user number a NULL pointer is returned.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
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

  !>Finalises the regions and destroys any current regions.
  SUBROUTINE REGIONS_FINALISE(ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nr

    CALL ENTERS("REGIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGIONS%WORLD_REGION)) THEN
      !Destroy any global region daughter regions first
      DO nr=1,REGIONS%WORLD_REGION%NUMBER_OF_SUB_REGIONS
        CALL REGION_DESTROY(REGIONS%WORLD_REGION%SUB_REGIONS(nr)%PTR%USER_NUMBER,ERR,ERROR,*999)
      ENDDO !region
      !Destroy global region and deallocated any memory allocated in the global region
      CALL REGION_FINALISE(REGIONS%WORLD_REGION,ERR,ERROR,*999)
      NULLIFY(REGIONS%WORLD_REGION)
    ENDIF
   
    CALL EXITS("REGIONS_FINALISE")
    RETURN
999 CALL ERRORS("REGIONS_FINALISE",ERR,ERROR)
    CALL EXITS("REGIONS_FINALISE")
    RETURN 1
  END SUBROUTINE REGIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the regions and creates the global world region.
  SUBROUTINE REGIONS_INITIALISE(WORLD_REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: WORLD_REGION !<On exit, a pointer to the world region. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: WORLD_COORDINATE_SYSTEM

    NULLIFY(WORLD_COORDINATE_SYSTEM)
    
    CALL ENTERS("REGIONS_INITIALISE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(WORLD_REGION)) THEN
      CALL FLAG_ERROR("World region is already associated.",ERR,ERROR,*999)
    ELSE
      CALL COORDINATE_SYSTEM_USER_NUMBER_FIND(0,WORLD_COORDINATE_SYSTEM,ERR,ERROR,*999)
      IF(ASSOCIATED(WORLD_COORDINATE_SYSTEM)) THEN        
        CALL REGION_INITIALISE(REGIONS%WORLD_REGION,ERR,ERROR,*999)
        REGIONS%WORLD_REGION%USER_NUMBER=0
        REGIONS%WORLD_REGION%LABEL="World Region"
        REGIONS%WORLD_REGION%COORDINATE_SYSTEM=>WORLD_COORDINATE_SYSTEM
        REGIONS%WORLD_REGION%REGION_FINISHED=.TRUE.
        !Return the pointer
        WORLD_REGION=>REGIONS%WORLD_REGION
      ELSE
        CALL FLAG_ERROR("World coordinate system has not been created.",ERR,ERROR,*999)
      ENDIF
    ENDIF
   
    CALL EXITS("REGIONS_INITIALISE")
    RETURN
999 CALL ERRORS("REGIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("REGIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE REGIONS_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE REGION_ROUTINES
