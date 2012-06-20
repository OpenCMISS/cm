!> \file
!> \author Tim Wu
!> \brief This module handles all data point routines.
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
!> Contributor(s): Chris Bradley
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

!> This module handles all data point routines.

MODULE DATA_POINT_ROUTINES

  USE BASE_ROUTINES
  USE COMP_ENVIRONMENT
  USE DATA_PROJECTION_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE STRINGS
  USE TREES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  !>Starts the process of creating data points for an interface or region
  INTERFACE DATA_POINTS_CREATE_START
    MODULE PROCEDURE DATA_POINTS_CREATE_START_REGION
    MODULE PROCEDURE DATA_POINTS_CREATE_START_INTERFACE
  END INTERFACE !DATA_POINTS_CREATE_START

  !>Initialises data points for an interface or region
  INTERFACE DATA_POINTS_INITIALISE
    MODULE PROCEDURE DATA_POINTS_INITIALISE_REGION
    MODULE PROCEDURE DATA_POINTS_INITIALISE_INTERFACE
  END INTERFACE !DATA_POINTS_INITIALIES

  !>Gets the label for a data point identified by a given global number.
  INTERFACE DATA_POINTS_LABEL_GET
    MODULE PROCEDURE DATA_POINTS_LABEL_GET_C
    MODULE PROCEDURE DATA_POINTS_LABEL_GET_VS
  END INTERFACE !DATA_POINTS_LABEL_SET

  !>Changes/sets the label for a data point identified by a given global number.
  INTERFACE DATA_POINTS_LABEL_SET
    MODULE PROCEDURE DATA_POINTS_LABEL_SET_C
    MODULE PROCEDURE DATA_POINTS_LABEL_SET_VS
  END INTERFACE !DATA_POINTS_LABEL_SET

  PUBLIC DATA_POINT_CHECK_EXISTS

  PUBLIC DATA_POINTS_CREATE_FINISH,DATA_POINTS_CREATE_START,DATA_POINTS_DESTROY
  
  PUBLIC DATA_POINTS_DATA_PROJECTION_GET

  PUBLIC DATA_POINTS_LABEL_GET,DATA_POINTS_LABEL_SET
  
  PUBLIC DATA_POINTS_VALUES_GET,DATA_POINTS_VALUES_SET

  PUBLIC DATA_POINTS_NUMBER_OF_DATA_POINTS_GET
  
  PUBLIC DATA_POINTS_USER_NUMBER_GET,DATA_POINTS_USER_NUMBER_SET
  
  PUBLIC DATA_POINTS_WEIGHTS_GET,DATA_POINTS_WEIGHTS_SET
  
  PUBLIC DATA_POINTS_PROJECTION_DISTANCE_GET,DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET
  
  PUBLIC DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET,DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET
  
  PUBLIC DATA_POINTS_PROJECTION_EXIT_TAG_GET,DATA_POINTS_PROJECTION_XI_GET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Checks that a user data point number is defined on the specified region.
  SUBROUTINE DATA_POINT_CHECK_EXISTS(DATA_POINTS,USER_NUMBER,DATA_POINT_EXISTS,GLOBAL_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to check
    INTEGER(INTG) :: USER_NUMBER !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: DATA_POINT_EXISTS !<On exit, is .TRUE. if the data point user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_NUMBER !<On exit, if the data point exists the global number corresponding to the user data point number. If the data point does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
   
    CALL ENTERS("DATA_POINT_CHECK_EXISTS",ERR,ERROR,*999)

    DATA_POINT_EXISTS=.FALSE.
    GLOBAL_NUMBER=0
    IF(ASSOCIATED(DATA_POINTS)) THEN
      NULLIFY(TREE_NODE)
      CALL TREE_SEARCH(DATA_POINTS%DATA_POINTS_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
      IF(ASSOCIATED(TREE_NODE)) THEN
        CALL TREE_NODE_VALUE_GET(DATA_POINTS%DATA_POINTS_TREE,TREE_NODE,GLOBAL_NUMBER,ERR,ERROR,*999)
        DATA_POINT_EXISTS=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_POINT_CHECK_EXISTS")
    RETURN
999 CALL ERRORS("DATA_POINT_CHECK_EXISTS",ERR,ERROR)
    CALL EXITS("DATA_POINT_CHECK_EXISTS")
    RETURN 1
   
  END SUBROUTINE DATA_POINT_CHECK_EXISTS

  !
  !================================================================================================================================
  !

  !>Finalises a data point and deallocates all memory
  SUBROUTINE DATA_POINT_FINALISE(DATA_POINT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DATA_POINT_TYPE),INTENT(OUT) :: DATA_POINT !<The data point to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DATA_POINT_FINALISE",ERR,ERROR,*999)

    DATA_POINT%GLOBAL_NUMBER=0
    DATA_POINT%USER_NUMBER=0
    IF(ALLOCATED(DATA_POINT%VALUES)) DEALLOCATE(DATA_POINT%VALUES)
    IF(ALLOCATED(DATA_POINT%WEIGHTS)) DEALLOCATE(DATA_POINT%WEIGHTS)
    IF(ALLOCATED(DATA_POINT%PROJECTION_XI)) DEALLOCATE(DATA_POINT%PROJECTION_XI)
    
    CALL EXITS("DATA_POINT_FINALISE")
    RETURN
999 CALL ERRORS("DATA_POINT_FINALISE",ERR,ERROR)
    CALL EXITS("DATA_POINT_FINALISE")
    RETURN 1  
 
  END SUBROUTINE DATA_POINT_FINALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating data points in the region.
  SUBROUTINE DATA_POINTS_CREATE_FINISH(DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to be finished
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: data_point_idx
    
    CALL ENTERS("DATA_POINTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        CALL FLAG_ERROR("Data points have already been finished.",ERR,ERROR,*999)
      ELSE
        DATA_POINTS%DATA_POINTS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN !<TODO Still Diagnostics 1??
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of data points = ",DATA_POINTS%NUMBER_OF_DATA_POINTS,ERR,ERROR,*999)
      DO data_point_idx=1,DATA_POINTS%NUMBER_OF_DATA_POINTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Data Points = ",data_point_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number    = ",DATA_POINTS%DATA_POINTS(data_point_idx)% &
          & GLOBAL_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number      = ",DATA_POINTS%DATA_POINTS(data_point_idx)% &
          & USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Label            = ",DATA_POINTS%DATA_POINTS(data_point_idx)%LABEL, &
          & ERR,ERROR,*999)
      ENDDO !data_point_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"User->Global number tree",ERR,ERROR,*999)
      CALL TREE_OUTPUT(DIAGNOSTIC_OUTPUT_TYPE,DATA_POINTS%DATA_POINTS_TREE,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_POINTS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("DATA_POINTS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("DATA_POINTS_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the process of creating generic data points
  SUBROUTINE DATA_POINTS_CREATE_START_GENERIC(DATA_POINTS,NUMBER_OF_DATA_POINTS,NUMBER_OF_DIMENSIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<On exit, a pointer to the created data points
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DATA_POINTS !<The number of data points to create
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions for data points values
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: INSERT_STATUS,data_point_idx,coord_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DATA_POINTS_CREATE_START_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(NUMBER_OF_DATA_POINTS>0) THEN
        ALLOCATE(DATA_POINTS%DATA_POINTS(NUMBER_OF_DATA_POINTS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data points.",ERR,ERROR,*999)
        DATA_POINTS%NUMBER_OF_DATA_POINTS=NUMBER_OF_DATA_POINTS
        CALL TREE_CREATE_START(DATA_POINTS%DATA_POINTS_TREE,ERR,ERROR,*999)
        CALL TREE_INSERT_TYPE_SET(DATA_POINTS%DATA_POINTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
        CALL TREE_CREATE_FINISH(DATA_POINTS%DATA_POINTS_TREE,ERR,ERROR,*999)
        !Set default data point numbers
        DO data_point_idx=1,DATA_POINTS%NUMBER_OF_DATA_POINTS
          DATA_POINTS%DATA_POINTS(data_point_idx)%GLOBAL_NUMBER=data_point_idx
          DATA_POINTS%DATA_POINTS(data_point_idx)%USER_NUMBER=data_point_idx
          DATA_POINTS%DATA_POINTS(data_point_idx)%LABEL=""
          ! initialise data points values to 0.0 and weights to 1.0
          ALLOCATE(DATA_POINTS%DATA_POINTS(data_point_idx)%VALUES(NUMBER_OF_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data points values("//TRIM(NUMBER_TO_VSTRING &
            & (data_point_idx,"*",ERR,ERROR))//").",ERR,ERROR,*999)
          ALLOCATE(DATA_POINTS%DATA_POINTS(data_point_idx)%WEIGHTS(NUMBER_OF_DIMENSIONS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data points weights("//TRIM(NUMBER_TO_VSTRING &
            & (data_point_idx,"*",ERR,ERROR))//").",ERR,ERROR,*999)              
          DO coord_idx=1,NUMBER_OF_DIMENSIONS
            DATA_POINTS%DATA_POINTS(data_point_idx)%VALUES(coord_idx)=0.0_DP
            DATA_POINTS%DATA_POINTS(data_point_idx)%WEIGHTS(coord_idx)=1.0_DP
          ENDDO
          CALL TREE_ITEM_INSERT(DATA_POINTS%DATA_POINTS_TREE,data_point_idx,data_point_idx,INSERT_STATUS,ERR,ERROR,*999)
        ENDDO !data_point_idx
        NULLIFY(DATA_POINTS%DATA_PROJECTION)
        DATA_POINTS%DATA_POINTS_PROJECTED=.FALSE.    
      ELSE
        LOCAL_ERROR="The specified number of data points of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))// &
          & " is invalid. The number of data points must be > 0."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_POINTS_CREATE_START_GENERIC")
    RETURN  
999 CALL ERRORS("DATA_POINTS_CREATE_START_GENERIC",ERR,ERROR)
    CALL EXITS("DATA_POINTS_CREATE_START_GENERIC")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating data points in an interface.
  SUBROUTINE DATA_POINTS_CREATE_START_INTERFACE(INTERFACE,NUMBER_OF_DATA_POINTS,DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface in which to create the data points
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DATA_POINTS !<The number of data points to create
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<On exit, a pointer to the created data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("DATA_POINTS_CREATE_START_INTERFACE",ERR,ERROR,*998)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(DATA_POINTS)) THEN
        CALL FLAG_ERROR("Data points is already associated.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(INTERFACE%DATA_POINTS)) THEN
          CALL FLAG_ERROR("Interface already has data points associated.",ERR,ERROR,*998)
        ELSE
          !Initialise the data points for the interface
          CALL DATA_POINTS_INITIALISE(INTERFACE,ERR,ERROR,*999)
          !Create the data points 
          CALL DATA_POINTS_CREATE_START_GENERIC(INTERFACE%DATA_POINTS,NUMBER_OF_DATA_POINTS,INTERFACE% &
            & COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
          !Return the pointer        
          DATA_POINTS=>INTERFACE%DATA_POINTS
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("DATA_POINTS_CREATE_START_INTERFACE")
    RETURN
999 CALL DATA_POINTS_FINALISE(INTERFACE%DATA_POINTS,DUMMY_ERR,DUMMY_ERROR,*998)    
998 CALL ERRORS("DATA_POINTS_CREATE_START_INTERFACE",ERR,ERROR)
    CALL EXITS("DATA_POINTS_CREATE_START_INTERFACE")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating data points in an region.
  SUBROUTINE DATA_POINTS_CREATE_START_REGION(REGION,NUMBER_OF_DATA_POINTS,DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region in which to create the data points
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DATA_POINTS !<The number of data points to create
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<On exit, a pointer to the created data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("DATA_POINTS_CREATE_START_REGION",ERR,ERROR,*998)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(DATA_POINTS)) THEN
        CALL FLAG_ERROR("Data points is already associated.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(REGION%DATA_POINTS)) THEN
          CALL FLAG_ERROR("Region already has data points associated.",ERR,ERROR,*998)
        ELSE
          !Initialise the data points for the region
          CALL DATA_POINTS_INITIALISE(REGION,ERR,ERROR,*999)
          !Create the data points 
          CALL DATA_POINTS_CREATE_START_GENERIC(REGION%DATA_POINTS,NUMBER_OF_DATA_POINTS,REGION%COORDINATE_SYSTEM% &
            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
          !Return the pointer        
          DATA_POINTS=>REGION%DATA_POINTS
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*998)
    ENDIF
    
    CALL EXITS("DATA_POINTS_CREATE_START_REGION")
    RETURN
999 CALL DATA_POINTS_FINALISE(REGION%DATA_POINTS,DUMMY_ERR,DUMMY_ERROR,*998)    
998 CALL ERRORS("DATA_POINTS_CREATE_START_REGION",ERR,ERROR)
    CALL EXITS("DATA_POINTS_CREATE_START_REGION")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_CREATE_START_REGION

     
  !
  !================================================================================================================================
  !

  !>Destroys data points. \see OPENCMISS::CMISSDataPointsDestroy
  SUBROUTINE DATA_POINTS_DESTROY(DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_POINTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(ASSOCIATED(DATA_POINTS%REGION)) THEN
        NULLIFY(DATA_POINTS%REGION%DATA_POINTS)
      ELSE
        IF(ASSOCIATED(DATA_POINTS%INTERFACE)) THEN
          NULLIFY(DATA_POINTS%INTERFACE%DATA_POINTS)
        ELSE
          CALL FLAG_ERROR("Data points region and interface are not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
      CALL DATA_POINTS_FINALISE(DATA_POINTS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_DESTROY")
    RETURN
999 CALL ERRORS("DATA_POINTS_DESTROY",ERR,ERROR)
    CALL EXITS("DATA_POINTS_DESTROY")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_DESTROY

  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. 
  SUBROUTINE DATA_POINTS_DATA_PROJECTION_GET(DATA_POINTS,DATA_PROJECTION,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to get the data projection for
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<On exit, a pointer to the data projection for the data points. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    CALL ENTERS("DATA_POINTS_DATA_PROJECTION_GET",ERR,ERROR,*999)
    
    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN 
        IF(ASSOCIATED(DATA_PROJECTION)) THEN
          CALL FLAG_ERROR("Data projection is already associated.",ERR,ERROR,*999)
        ELSE
          DATA_PROJECTION=>DATA_POINTS%DATA_PROJECTION
          IF(.NOT.ASSOCIATED(DATA_PROJECTION)) CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF    
    
    CALL EXITS("DATA_POINTS_DATA_PROJECTION_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_DATA_PROJECTION_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_DATA_PROJECTION_GET")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_DATA_PROJECTION_GET

  !
  !===============================================================================================================================
  !

  !>Finalises the data points and deallocates any memory. 
  SUBROUTINE DATA_POINTS_FINALISE(DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: data_point_idx

    CALL ENTERS("DATA_POINTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(ALLOCATED(DATA_POINTS%DATA_POINTS)) THEN
        DO data_point_idx=1,SIZE(DATA_POINTS%DATA_POINTS,1)
          CALL DATA_POINT_FINALISE(DATA_POINTS%DATA_POINTS(data_point_idx),ERR,ERROR,*999)
        ENDDO !data_point_idx
        DEALLOCATE(DATA_POINTS%DATA_POINTS)
      ENDIF
      IF(ASSOCIATED(DATA_POINTS%DATA_POINTS_TREE)) CALL TREE_DESTROY(DATA_POINTS%DATA_POINTS_TREE,ERR,ERROR,*999)
      DEALLOCATE(DATA_POINTS)
    ENDIF
    
    CALL EXITS("DATA_POINTS_FINALISE")
    RETURN
999 CALL ERRORS("DATA_POINTS_FINALISE",ERR,ERROR)
    CALL EXITS("DATA_POINTS_FINALISE")
    RETURN 1

  END SUBROUTINE DATA_POINTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the data points.
  SUBROUTINE DATA_POINTS_INITIALISE_GENERIC(DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DATA_POINTS_INITIALISE_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      NULLIFY(DATA_POINTS%REGION)
      NULLIFY(DATA_POINTS%INTERFACE)
      DATA_POINTS%DATA_POINTS_FINISHED=.FALSE.
      DATA_POINTS%NUMBER_OF_DATA_POINTS=0
      NULLIFY(DATA_POINTS%DATA_POINTS_TREE)
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_INITIALISE_GENERIC")
    RETURN
999 CALL ERRORS("DATA_POINTS_INITIALISE_GENERIC",ERR,ERROR)
    CALL EXITS("DATA_POINTS_INITIALISE_GENERIC")
    RETURN 1
  END SUBROUTINE DATA_POINTS_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the data points in a given interface.
  SUBROUTINE DATA_POINTS_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the data points for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DATA_POINTS_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%DATA_POINTS)) THEN
        CALL FLAG_ERROR("Interface already has associated data points.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(INTERFACE%DATA_POINTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate interface data points.",ERR,ERROR,*999)
        CALL DATA_POINTS_INITIALISE_GENERIC(INTERFACE%DATA_POINTS,ERR,ERROR,*999)
        INTERFACE%DATA_POINTS%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_INITIALISE_INTERFACE")
    RETURN
999 CALL ERRORS("DATA_POINTS_INITIALISE_INTERFACE",ERR,ERROR)
    CALL EXITS("DATA_POINTS_INITIALISE_INTERFACE")
    RETURN 1
    
  END SUBROUTINE DATA_POINTS_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the data points in a given region.
  SUBROUTINE DATA_POINTS_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the data points for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DATA_POINTS_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%DATA_POINTS)) THEN
        CALL FLAG_ERROR("Region has associated data points.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(REGION%DATA_POINTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate region data points.",ERR,ERROR,*999)
        CALL DATA_POINTS_INITIALISE_GENERIC(REGION%DATA_POINTS,ERR,ERROR,*999)
        REGION%DATA_POINTS%REGION=>REGION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_INITIALISE_REGION")
    RETURN
999 CALL ERRORS("DATA_POINTS_INITIALISE_REGION",ERR,ERROR)
    CALL EXITS("DATA_POINTS_INITIALISE_REGION")
    RETURN 1
  END SUBROUTINE DATA_POINTS_INITIALISE_REGION
       
  !
  !================================================================================================================================
  !

  !>Gets the character label for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_LABEL_GET_C(DATA_POINTS,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to get the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On exit, the label of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER :: C_LENGTH,VS_LENGTH
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          C_LENGTH=LEN(LABEL)
          VS_LENGTH=LEN_TRIM(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%LABEL)
          IF(C_LENGTH>VS_LENGTH) THEN
            LABEL=CHAR(LEN_TRIM(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%LABEL))
          ELSE
            LABEL=CHAR(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%LABEL,C_LENGTH)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_LABEL_GET_C")
    RETURN
999 CALL ERRORS("DATA_POINTS_LABEL_GET_C",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_LABEL_GET_C")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_LABEL_GET_C
        
  !
  !================================================================================================================================
  !

  !>Gets the varying string label for a data point identified by a given global number. 
  SUBROUTINE DATA_POINTS_LABEL_GET_VS(DATA_POINTS,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to get the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On exit, the label of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          LABEL=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%LABEL
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_LABEL_GET_VS")
    RETURN
999 CALL ERRORS("DATA_POINTS_LABEL_GET_VS",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_LABEL_GET_VS")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Changes/sets the character label for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_LABEL_SET_C(DATA_POINTS,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        CALL FLAG_ERROR("Data points have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%LABEL=LABEL
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_LABEL_SET_C")
    RETURN
999 CALL ERRORS("DATA_POINTS_LABEL_SET_C",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_LABEL_SET_C")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_LABEL_SET_C    
  
  !
  !================================================================================================================================
  !


  !>Changes/sets the varying string label for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_LABEL_SET_VS(DATA_POINTS,GLOBAL_NUMBER,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the label for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        CALL FLAG_ERROR("Data points have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%LABEL=LABEL
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_LABEL_SET_VS")
    RETURN
999 CALL ERRORS("DATA_POINTS_LABEL_SET_VS",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_LABEL_SET_VS")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_LABEL_SET_VS
        
  !
  !================================================================================================================================
  !
  
  !>Gets the values for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_VALUES_GET(DATA_POINTS,GLOBAL_NUMBER,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the values for
    REAL(DP), INTENT(OUT) :: VALUES(:) !<On exit, the values of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_VALUES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        !Check the data point global number exists
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          IF(SIZE(VALUES,1)==SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%VALUES,1)) THEN
            VALUES=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%VALUES
          ELSE
            CALL FLAG_ERROR("array values has size of "//TRIM(NUMBER_TO_VSTRING(SIZE(VALUES,1),"*",ERR,ERROR))// &
              & "but it needs to have size of "// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%VALUES,1),"*",ERR,ERROR))// &
              & "." ,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_VALUES_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_VALUES_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_VALUES_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_VALUES_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the values for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_VALUES_SET(DATA_POINTS,GLOBAL_NUMBER,VALUES,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the values for
    REAL(DP), INTENT(IN) :: VALUES(:) !<The global number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_VALUES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        CALL FLAG_ERROR("Data points have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          IF(SIZE(VALUES,1)==SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%VALUES,1)) THEN
            DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%VALUES=VALUES
          ELSE
            CALL FLAG_ERROR("The dimension of the input values does not match.",ERR,ERROR,*999)    
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_VALUES_SET")
    RETURN
999 CALL ERRORS("DATA_POINTS_VALUES_SET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_VALUES_SET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_VALUES_SET

  !
  !================================================================================================================================
  !

  !>Returns the number of data points. \see OPENCMISS::CMISSDataPointsNumberOfDataPointsGet
  SUBROUTINE DATA_POINTS_NUMBER_OF_DATA_POINTS_GET(DATA_POINTS,NUMBER_OF_DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to get the number of data points for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_DATA_POINTS !<On return, the number of data points
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_POINTS_NUMBER_OF_DATA_POINTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        NUMBER_OF_DATA_POINTS=DATA_POINTS%NUMBER_OF_DATA_POINTS
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_NUMBER_OF_DATA_POINTS_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_NUMBER_OF_DATA_POINTS_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_NUMBER_OF_DATA_POINTS_GET")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_NUMBER_OF_DATA_POINTS_GET

  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. 
  SUBROUTINE DATA_POINTS_USER_NUMBER_GET(DATA_POINTS,GLOBAL_NUMBER,USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to get the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the user number for
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<On exit, the user number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_USER_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          USER_NUMBER=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_USER_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_USER_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_USER_NUMBER_GET")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_USER_NUMBER_GET
        
  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_USER_NUMBER_SET(DATA_POINTS,GLOBAL_NUMBER,USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the user number for
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: INSERT_STATUS,OLD_GLOBAL_NUMBER
    LOGICAL :: DATA_POINT_EXISTS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_USER_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        CALL FLAG_ERROR("Data points have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          !Check the data point user number is not already used
          CALL DATA_POINT_CHECK_EXISTS(DATA_POINTS,USER_NUMBER,DATA_POINT_EXISTS,OLD_GLOBAL_NUMBER,ERR,ERROR,*999)
          IF(DATA_POINT_EXISTS) THEN
            IF(OLD_GLOBAL_NUMBER/=GLOBAL_NUMBER) THEN
              LOCAL_ERROR="The specified data point user number of "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                & " is already used by global data point number "//TRIM(NUMBER_TO_VSTRING(OLD_GLOBAL_NUMBER,"*",ERR,ERROR))// &
                & ". User data point numbers must be unique."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL TREE_ITEM_DELETE(DATA_POINTS%DATA_POINTS_TREE,DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
            CALL TREE_ITEM_INSERT(DATA_POINTS%DATA_POINTS_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
            IF(INSERT_STATUS/=TREE_NODE_INSERT_SUCESSFUL) CALL FLAG_ERROR("Unsucessful data points tree insert.",ERR,ERROR,*999)
            DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_USER_NUMBER_SET")
    RETURN
999 CALL ERRORS("DATA_POINTS_USER_NUMBER_SET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_USER_NUMBER_SET")
    RETURN 1
   
  END SUBROUTINE DATA_POINTS_USER_NUMBER_SET
  
  !
  !================================================================================================================================
  !
  
  !>Gets the weights for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_WEIGHTS_GET(DATA_POINTS,GLOBAL_NUMBER,WEIGHTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the weights for
    REAL(DP), INTENT(OUT) :: WEIGHTS(:) !<On exit, the weights of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_WEIGHTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        !Check the data point global number exists
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          IF(SIZE(WEIGHTS,1)==SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%WEIGHTS,1)) THEN
            WEIGHTS=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%WEIGHTS
          ELSE
            CALL FLAG_ERROR("array weights has size of "//TRIM(NUMBER_TO_VSTRING(SIZE(WEIGHTS,1),"*",ERR,ERROR))// &
              & "but it needs to have size of "// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%WEIGHTS,1),"*",ERR,ERROR))// &
              & "." ,ERR,ERROR,*999)
          ENDIF          
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF        
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_WEIGHTS_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_WEIGHTS_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_WEIGHTS_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_WEIGHTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the weights for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_WEIGHTS_SET(DATA_POINTS,GLOBAL_NUMBER,WEIGHTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to set the weights for
    REAL(DP), INTENT(IN) :: WEIGHTS(:) !<The global number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_WEIGHTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        CALL FLAG_ERROR("Data points have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
          IF(SIZE(WEIGHTS,1)==SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%WEIGHTS,1)) THEN
            DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%WEIGHTS=WEIGHTS
          ELSE
            CALL FLAG_ERROR("The dimension of the input weights does not match.",ERR,ERROR,*999)    
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global data point number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_WEIGHTS_SET")
    RETURN
999 CALL ERRORS("DATA_POINTS_WEIGHTS_SET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_WEIGHTS_SET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_WEIGHTS_SET  
        
  !
  !================================================================================================================================
  !
  
  !>Gets the projection distance for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_PROJECTION_DISTANCE_GET(DATA_POINTS,GLOBAL_NUMBER,PROJECTION_DISTANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the projection distance for
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE !<On exit, the projection distance of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_PROJECTION_DISTANCE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(DATA_POINTS%DATA_POINTS_PROJECTED) THEN
          !Check the data point global number exists
          IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
            PROJECTION_DISTANCE=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_DISTANCE
          ELSE
            LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The global data point number should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_PROJECTION_DISTANCE_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_PROJECTION_DISTANCE_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_PROJECTION_DISTANCE_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_PROJECTION_DISTANCE_GET


  !
  !================================================================================================================================
  !
  
  !>Gets the projection element number for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET(DATA_POINTS,GLOBAL_NUMBER,PROJECTION_ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the projection element number for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER !<On exit, the projection element number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(DATA_POINTS%DATA_POINTS_PROJECTED) THEN
          !Check the data point global number exists
          IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
            PROJECTION_ELEMENT_NUMBER=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_ELEMENT_NUMBER
          ELSE
            LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The global data point number should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_PROJECTION_ELEMENT_NUMBER_GET
  
  !
  !================================================================================================================================
  !
  
  !>Gets the projection element face number for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET(DATA_POINTS,GLOBAL_NUMBER,PROJECTION_ELEMENT_FACE_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the projection element face number for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_FACE_NUMBER !<On exit, the projection element face number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(DATA_POINTS%DATA_POINTS_PROJECTED) THEN
          ! Check if boundary faces projection type was set
          IF(DATA_POINTS%DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)THEN
            !Check the data point global number exists
            IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
              PROJECTION_ELEMENT_FACE_NUMBER=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_ELEMENT_FACE_NUMBER
            ELSE
              LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The global data point number should be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Data points data projection projection type is not set to boundary faces projection type.", &
              & ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_PROJECTION_ELEMENT_FACE_NUMBER_GET
  
  !
  !================================================================================================================================
  !
  
  !>Gets the projection element line number for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET(DATA_POINTS,GLOBAL_NUMBER,PROJECTION_ELEMENT_LINE_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the projection element line number for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_LINE_NUMBER !<On exit, the projection element line number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(DATA_POINTS%DATA_POINTS_PROJECTED) THEN
          ! Check if boundary lines projection type was set
          IF(DATA_POINTS%DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)THEN
            !Check the data point global number exists
            IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
              PROJECTION_ELEMENT_LINE_NUMBER=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_ELEMENT_LINE_NUMBER
            ELSE
              LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The global data point number should be between 1 and "// &
                & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Data points data projection projection type is not set to boundary lines projection type.", &
              & ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_PROJECTION_ELEMENT_LINE_NUMBER_GET

  !
  !================================================================================================================================
  !
  
  !>Gets the projection exit tag for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_PROJECTION_EXIT_TAG_GET(DATA_POINTS,GLOBAL_NUMBER,PROJECTION_EXIT_TAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the projection exit tag for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG !<On exit, the projection exit tag of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_PROJECTION_EXIT_TAG_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(DATA_POINTS%DATA_POINTS_PROJECTED) THEN
          !Check the data point global number exists
          IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
            PROJECTION_EXIT_TAG=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_EXIT_TAG
          ELSE
            LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The global data point number should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_PROJECTION_EXIT_TAG_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_PROJECTION_EXIT_TAG_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_PROJECTION_EXIT_TAG_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_PROJECTION_EXIT_TAG_GET

  !
  !================================================================================================================================
  !
  
  !>Gets the projection xi for a data point identified by a given global number.
  SUBROUTINE DATA_POINTS_PROJECTION_XI_GET(DATA_POINTS,GLOBAL_NUMBER,PROJECTION_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to set the number for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the projection xi for
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(:) !<On exit, the projection xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_POINTS_PROJECTION_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(DATA_POINTS%DATA_POINTS_PROJECTED) THEN
          !Check the data point global number exists
          IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DATA_POINTS%NUMBER_OF_DATA_POINTS) THEN
            IF(SIZE(PROJECTION_XI,1)==SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_XI,1)) THEN
              PROJECTION_XI=DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_XI
            ELSE
              CALL FLAG_ERROR("projection xi has size of "//TRIM(NUMBER_TO_VSTRING(SIZE(PROJECTION_XI,1),"*",ERR,ERROR))// &
                & "but it needs to have size of "// &
                & TRIM(NUMBER_TO_VSTRING(SIZE(DATA_POINTS%DATA_POINTS(GLOBAL_NUMBER)%PROJECTION_XI,1),"*",ERR,ERROR))// &
                & "." ,ERR,ERROR,*999)
            ENDIF
          ELSE
            LOCAL_ERROR="The specified global data point number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
              & " is invalid. The global data point number should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DATA_POINTS%NUMBER_OF_DATA_POINTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points have not been projected.",ERR,ERROR,*999)
        ENDIF   
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_POINTS_PROJECTION_XI_GET")
    RETURN
999 CALL ERRORS("DATA_POINTS_PROJECTION_XI_GET",ERR,ERROR)    
    CALL EXITS("DATA_POINTS_PROJECTION_XI_GET")
    RETURN 1

  END SUBROUTINE DATA_POINTS_PROJECTION_XI_GET
        
  !
  !================================================================================================================================
  !

END MODULE DATA_POINT_ROUTINES



