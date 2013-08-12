!> \file
!> \author Tim Wu
!> \brief This module handles all data projection routines
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
!> Contributor(s): Chris Bradley, Kumar Mithraratne, Prasad Babarenda Gamage
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

!> This module handles all data projection routines.

MODULE DATA_PROJECTION_ROUTINES

  USE BASE_ROUTINES
  USE CMISS_MPI  
  USE COMP_ENVIRONMENT
  USE DOMAIN_MAPPINGS
  USE FIELD_ROUTINES  
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE MESH_ROUTINES
  USE MPI
  USE SORTING
  USE STRINGS
  USE TREES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup DATA_POINT_PROJECTION_ROUTINES_DataProjectionTypes DATA_POINT_PROJECTION_ROUTINES::DataProjectionTypes
  !> \brief Datapoint projection definition type parameters
  !> \see DATA_POINT_PROJECTION_ROUTINES,OPENCMISS_DataProjectionTypes
  !>@{ 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE=1 !<The boundary line projection type for data projection, only projects to boundary lines of the mesh. \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE=2 !<The boundary face projection type for data projection, only projects to boundary faces of the mesh. \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE=3 !<The element projection type for data projection, projects to all elements in mesh. \see DATA_PROJECTION_ROUTINES 
  !>@}

  !> \addtogroup DATA_POINT_PROJECTION_ROUTINES_DataProjectionExitTags DATA_POINT_PROJECTION_ROUTINES::DataProjectionExitTags
  !> \brief Datapoint projection exit tags
  !> \see DATA_POINT_PROJECTION_ROUTINES,OPENCMISS_DataProjectionExitTags
  !>@{ 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_CONVERGED=1 !<Data projection exited due to it being converged \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_BOUNDS=2 !<Data projection exited due to it hitting the bound and continue to travel out of the element. \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_MAX_ITERATION=3 !<Data projection exited due to it attaining maximum number of iteration specified by user. \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_NO_ELEMENT=4 !<Data projection exited due to no local element found, this happens when none of the candidate elements are within this computational node, and before MPI communication with other nodes. \see DATA_PROJECTION_ROUTINES     
  !>@}

  !Module types

  !Module variables

  !Interfaces
  
  !>Gets the label for a data projection.
  INTERFACE DATA_PROJECTION_LABEL_GET
    MODULE PROCEDURE DATA_PROJECTION_LABEL_GET_C
    MODULE PROCEDURE DATA_PROJECTION_LABEL_GET_VS
  END INTERFACE !DATA_PROJECTION_LABEL_GET
  
  !>Sets/changes the label for a data projection.
  INTERFACE DATA_PROJECTION_LABEL_SET
    MODULE PROCEDURE DATA_PROJECTION_LABEL_SET_C
    MODULE PROCEDURE DATA_PROJECTION_LABEL_SET_VS
  END INTERFACE !DATA_PROJECTION_LABEL_SET  

  PUBLIC DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE,DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE, &
    & DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
  
  PUBLIC DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET,DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET

  PUBLIC DATA_PROJECTION_CREATE_FINISH,DATA_PROJECTION_CREATE_START_DATA_POINTS
  
  PUBLIC DATA_PROJECTION_DESTROY
  
  PUBLIC DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE
  
  PUBLIC DataProjection_DataPointsPositionEvaluate
  
  PUBLIC DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET,DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET
  
  PUBLIC DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET,DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET
  
  PUBLIC DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET,DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET
  
  PUBLIC DataProjection_ProjectionCandidatesSet
  
  PUBLIC DATA_PROJECTION_PROJECTION_TYPE_GET,DATA_PROJECTION_PROJECTION_TYPE_SET
  
  PUBLIC DATA_PROJECTION_RELATIVE_TOLERANCE_GET,DATA_PROJECTION_RELATIVE_TOLERANCE_SET
  
  PUBLIC DATA_PROJECTION_STARTING_XI_GET,DATA_PROJECTION_STARTING_XI_SET

  PUBLIC DATA_PROJECTION_RESULT_XI_GET, DATA_PROJECTION_RESULT_XI_SET

  PUBLIC DATA_PROJECTION_ELEMENT_SET
  
  PUBLIC DATA_PROJECTION_RESULT_DISTANCE_GET
  
  PUBLIC DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET,DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET
  
  PUBLIC DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET
  
  PUBLIC DATA_PROJECTION_RESULT_EXIT_TAG_GET
  
  PUBLIC DATA_PROJECTION_LABEL_GET,DATA_PROJECTION_LABEL_SET

  PUBLIC DataProjection_ElementGet,DataProjection_DistanceGet

CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Gets the absolute tolerance for a data projection.
  SUBROUTINE DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET(DATA_PROJECTION,ABSOLUTE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the absolute tolerance for
    REAL(DP), INTENT(OUT) :: ABSOLUTE_TOLERANCE !<On exit, the absolute tolerance of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE       
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_ABSOLUTE_TOLERANCE_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the absolute tolerance for a data projection.
  SUBROUTINE DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET(DATA_PROJECTION,ABSOLUTE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the absolute tolerance for
    REAL(DP), INTENT(IN) :: ABSOLUTE_TOLERANCE !<the absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE      
        IF(ABSOLUTE_TOLERANCE>=0) THEN
          DATA_PROJECTION%ABSOLUTE_TOLERANCE=ABSOLUTE_TOLERANCE
        ELSE
          CALL FLAG_ERROR("Data projection absolute tolerance must be a positive real number.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_ABSOLUTE_TOLERANCE_SET

  !
  !================================================================================================================================
  !
  
  !>Find the closest elements to a data point based on starting xi guess.
  SUBROUTINE DATA_PROJECTION_CLOSEST_ELEMENTS_FIND(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES, &
    & CANDIDATE_ELEMENTS,NUMBER_OF_CANDIDATES,CLOSEST_ELEMENTS,CLOSEST_DISTANCES,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_CANDIDATES
    INTEGER(INTG), INTENT(OUT) :: CLOSEST_ELEMENTS(:) !<List of shortest element number
    REAL(DP), INTENT(OUT) :: CLOSEST_DISTANCES(:) !<List of shortest distance (squared)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    INTEGER(INTG) :: REGION_DIMENSIONS !<Region coordinate dimension
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_CANDIDATES !<Number of closest elements to record
    REAL(DP) :: DISTANCE_VECTOR(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: DISTANCE2 !<Distance squared
    INTEGER(INTG) :: nce
    INTEGER(INTG) :: ELEMENT_NUMBER,insert_idx
      
    CALL ENTERS("DATA_PROJECTION_CLOSEST_ELEMENTS_FIND",ERR,ERROR,*999)
    
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN 
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        NUMBER_OF_CLOSEST_CANDIDATES=SIZE(CLOSEST_ELEMENTS,1)    
        !loop through the first few elements
        DO nce=1,NUMBER_OF_CLOSEST_CANDIDATES
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(nce)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999, &
            & FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS) = POINT_VALUES-INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          DO insert_idx=nce,1,-1
            IF(insert_idx>1) THEN
              IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
            ENDIF
            !insert the element number into the correct index
            IF(insert_idx<nce) THEN
              CLOSEST_DISTANCES(insert_idx+1:nce)=CLOSEST_DISTANCES(insert_idx:nce-1)
              CLOSEST_ELEMENTS(insert_idx+1:nce)=CLOSEST_ELEMENTS(insert_idx:nce-1) 
            ENDIF
            CLOSEST_DISTANCES(insert_idx)=DISTANCE2
            CLOSEST_ELEMENTS(insert_idx)=ELEMENT_NUMBER  
            EXIT                        
          ENDDO
        ENDDO  
        !Loop through the rest of the elements
        DO nce=NUMBER_OF_CLOSEST_CANDIDATES+1,NUMBER_OF_CANDIDATES
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(nce)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999, &
            & FIELD_GEOMETRIC_COMPONENTS_TYPE) 
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES - INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          IF(DISTANCE2<CLOSEST_DISTANCES(NUMBER_OF_CLOSEST_CANDIDATES))THEN
            DO insert_idx=NUMBER_OF_CLOSEST_CANDIDATES,1,-1
              IF(insert_idx>1) THEN
                IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
              ENDIF
              !insert the element into the correct index
              IF(insert_idx<NUMBER_OF_CLOSEST_CANDIDATES) THEN
                CLOSEST_DISTANCES(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_DISTANCES(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
                CLOSEST_ELEMENTS(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_ELEMENTS(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
              ENDIF
              CLOSEST_DISTANCES(insert_idx)=DISTANCE2
              CLOSEST_ELEMENTS(insert_idx)=ELEMENT_NUMBER  
              EXIT 
            ENDDO         
          ENDIF        
        ENDDO
        !CLOSEST_DISTANCES=SQRT(CLOSEST_DISTANCES) !return shortest distances
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF    
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_CLOSEST_ELEMENTS_FIND")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_CLOSEST_ELEMENTS_FIND",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_CLOSEST_ELEMENTS_FIND")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_CLOSEST_ELEMENTS_FIND
  
  !
  !================================================================================================================================
  !
  
  !>Find the closest faces to a data point base on starting xi guess.
  SUBROUTINE DATA_PROJECTION_CLOSEST_FACES_FIND(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & CANDIDATE_ELEMENT_FACES,NUMBER_OF_CANDIDATES,CLOSEST_ELEMENTS,CLOSEST_ELEMENT_FACES,CLOSEST_DISTANCES,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENT_FACES(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_CANDIDATES
    INTEGER(INTG), INTENT(OUT) :: CLOSEST_ELEMENTS(:) !<List of shortest element number
    INTEGER(INTG), INTENT(OUT) :: CLOSEST_ELEMENT_FACES(:) !<List of shortest element face number
    REAL(DP), INTENT(OUT) :: CLOSEST_DISTANCES(:) !<List of shortest distance (squared)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: REGION_DIMENSIONS !<Region coordinate dimension
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_CANDIDATES !<Number of closest elements to record
    REAL(DP) :: DISTANCE_VECTOR(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: DISTANCE2 !<Distance squared
    INTEGER(INTG) :: nce
    INTEGER(INTG) :: ELEMENT_NUMBER,ELEMENT_FACE_NUMBER,FACE_NUMBER,insert_idx
      
    CALL ENTERS("DATA_PROJECTION_CLOSEST_FACES_FIND",ERR,ERROR,*999)
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN 
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        NUMBER_OF_CLOSEST_CANDIDATES=SIZE(CLOSEST_ELEMENTS,1)    
        !loop through the first few faces
        DO nce=1,NUMBER_OF_CLOSEST_CANDIDATES
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(nce)
          ELEMENT_FACE_NUMBER=CANDIDATE_ELEMENT_FACES(nce)
          FACE_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS( &
            & ELEMENT_NUMBER)%ELEMENT_FACES(ELEMENT_FACE_NUMBER)
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999, &
            & FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS) = POINT_VALUES-INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          DO insert_idx=nce,1,-1
            IF(insert_idx>1) THEN
              IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
            ENDIF
            !insert the element number into the correct index
            IF(insert_idx<nce) THEN
              CLOSEST_DISTANCES(insert_idx+1:nce)=CLOSEST_DISTANCES(insert_idx:nce-1)
              CLOSEST_ELEMENTS(insert_idx+1:nce)=CLOSEST_ELEMENTS(insert_idx:nce-1)
              CLOSEST_ELEMENT_FACES(insert_idx+1:nce)=CLOSEST_ELEMENT_FACES(insert_idx:nce-1) 
            ENDIF
            CLOSEST_DISTANCES(insert_idx)=DISTANCE2
            CLOSEST_ELEMENTS(insert_idx)=ELEMENT_NUMBER
            CLOSEST_ELEMENT_FACES(insert_idx)=ELEMENT_FACE_NUMBER
            EXIT                        
          ENDDO
        ENDDO  
        !Loop through the rest of the faces
        DO nce=NUMBER_OF_CLOSEST_CANDIDATES+1,NUMBER_OF_CANDIDATES
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(nce)
          ELEMENT_FACE_NUMBER=CANDIDATE_ELEMENT_FACES(nce)
          FACE_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS( &
            & ELEMENT_NUMBER)%ELEMENT_FACES(ELEMENT_FACE_NUMBER)          
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999, &
            & FIELD_GEOMETRIC_COMPONENTS_TYPE) 
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES - INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          IF(DISTANCE2<CLOSEST_DISTANCES(NUMBER_OF_CLOSEST_CANDIDATES))THEN
            DO insert_idx=NUMBER_OF_CLOSEST_CANDIDATES,1,-1
              IF(insert_idx>1) THEN
                IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
              ENDIF
              !insert the element into the correct index
              IF(insert_idx<NUMBER_OF_CLOSEST_CANDIDATES) THEN
                CLOSEST_DISTANCES(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_DISTANCES(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
                CLOSEST_ELEMENTS(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_ELEMENTS(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
                CLOSEST_ELEMENT_FACES(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_ELEMENT_FACES(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
              ENDIF
              CLOSEST_DISTANCES(insert_idx)=DISTANCE2
              CLOSEST_ELEMENTS(insert_idx)=ELEMENT_NUMBER
              CLOSEST_ELEMENT_FACES(insert_idx)=ELEMENT_FACE_NUMBER
              EXIT 
            ENDDO         
          ENDIF        
        ENDDO
        !CLOSEST_DISTANCES=SQRT(CLOSEST_DISTANCES) !return shortest distances
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF    
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_CLOSEST_FACES_FIND")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_CLOSEST_FACES_FIND",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_CLOSEST_FACES_FIND")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_CLOSEST_FACES_FIND
  
  !
  !================================================================================================================================
  !
  
  !>Find the closest lines to a data point base on starting xi guess.
  SUBROUTINE DATA_PROJECTION_CLOSEST_LINES_FIND(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & CANDIDATE_ELEMENT_LINES,NUMBER_OF_CANDIDATES,CLOSEST_ELEMENTS,CLOSEST_ELEMENT_LINES,CLOSEST_DISTANCES,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENT_LINES(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_CANDIDATES
    INTEGER(INTG), INTENT(OUT) :: CLOSEST_ELEMENTS(:) !<List of shortest element number
    INTEGER(INTG), INTENT(OUT) :: CLOSEST_ELEMENT_LINES(:) !<List of shortest element line number
    REAL(DP), INTENT(OUT) :: CLOSEST_DISTANCES(:) !<List of shortest distance (squared)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: REGION_DIMENSIONS !<Region coordinate dimension
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_CANDIDATES !<Number of closest elements to record
    REAL(DP) :: DISTANCE_VECTOR(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: DISTANCE2 !<Distance squared
    INTEGER(INTG) :: nce
    INTEGER(INTG) :: ELEMENT_NUMBER,ELEMENT_LINE_NUMBER,LINE_NUMBER,insert_idx
      
    CALL ENTERS("DATA_PROJECTION_CLOSEST_LINES_FIND",ERR,ERROR,*999)
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        NUMBER_OF_CLOSEST_CANDIDATES=SIZE(CLOSEST_ELEMENTS,1)    
        !loop through the first few lines
        DO nce=1,NUMBER_OF_CLOSEST_CANDIDATES
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(nce)
          ELEMENT_LINE_NUMBER=CANDIDATE_ELEMENT_LINES(nce)
          LINE_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS( &
            & ELEMENT_NUMBER)%ELEMENT_LINES(ELEMENT_LINE_NUMBER)
          CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,LINE_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999, &
            & FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS) = POINT_VALUES-INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          DO insert_idx=nce,1,-1
            IF(insert_idx>1) THEN
              IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
            ENDIF
            !insert the element number into the correct index
            IF(insert_idx<nce) THEN
              CLOSEST_DISTANCES(insert_idx+1:nce)=CLOSEST_DISTANCES(insert_idx:nce-1)
              CLOSEST_ELEMENTS(insert_idx+1:nce)=CLOSEST_ELEMENTS(insert_idx:nce-1)
              CLOSEST_ELEMENT_LINES(insert_idx+1:nce)=CLOSEST_ELEMENT_LINES(insert_idx:nce-1) 
            ENDIF
            CLOSEST_DISTANCES(insert_idx)=DISTANCE2
            CLOSEST_ELEMENTS(insert_idx)=ELEMENT_NUMBER
            CLOSEST_ELEMENT_LINES(insert_idx)=ELEMENT_LINE_NUMBER
            EXIT                        
          ENDDO
        ENDDO  
        !Loop through the rest of the lines
        DO nce=NUMBER_OF_CLOSEST_CANDIDATES+1,NUMBER_OF_CANDIDATES
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(nce)
          ELEMENT_LINE_NUMBER=CANDIDATE_ELEMENT_LINES(nce)
          LINE_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS( &
            & ELEMENT_NUMBER)%ELEMENT_LINES(ELEMENT_LINE_NUMBER)          
          CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,LINE_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999, &
            & FIELD_GEOMETRIC_COMPONENTS_TYPE) 
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES - INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          IF(DISTANCE2<CLOSEST_DISTANCES(NUMBER_OF_CLOSEST_CANDIDATES))THEN
            DO insert_idx=NUMBER_OF_CLOSEST_CANDIDATES,1,-1
              IF(insert_idx>1) THEN
                IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
              ENDIF
              !insert the element into the correct index
              IF(insert_idx<NUMBER_OF_CLOSEST_CANDIDATES) THEN
                CLOSEST_DISTANCES(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_DISTANCES(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
                CLOSEST_ELEMENTS(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_ELEMENTS(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
                CLOSEST_ELEMENT_LINES(insert_idx+1:NUMBER_OF_CLOSEST_CANDIDATES)=CLOSEST_ELEMENT_LINES(insert_idx: &
                  & NUMBER_OF_CLOSEST_CANDIDATES-1)
              ENDIF
              CLOSEST_DISTANCES(insert_idx)=DISTANCE2
              CLOSEST_ELEMENTS(insert_idx)=ELEMENT_NUMBER
              CLOSEST_ELEMENT_LINES(insert_idx)=ELEMENT_LINE_NUMBER
              EXIT 
            ENDDO         
          ENDIF        
        ENDDO
        !CLOSEST_DISTANCES=SQRT(CLOSEST_DISTANCES) !return shortest distances
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF    
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_CLOSEST_LINES_FIND")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_CLOSEST_LINES_FIND",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_CLOSEST_LINES_FIND")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_CLOSEST_LINES_FIND
  
  !
  !================================================================================================================================
  !
  
  !>Finishes the process of creating data projection.
  SUBROUTINE DATA_PROJECTION_CREATE_FINISH(DATA_PROJECTION,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<On exit, a pointer to the created data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables   

    CALL ENTERS("DATA_PROJECTION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN !Has to be associated
      IF(ASSOCIATED(DATA_PROJECTION%DATA_POINTS)) THEN !Has to be associated
        IF(DATA_PROJECTION%DATA_POINTS%DATA_POINTS_FINISHED) THEN !Has to be finished
          IF(ASSOCIATED(DATA_PROJECTION%MESH)) THEN !Has to be associated
            IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN !Cannot be finished
              CALL FLAG_ERROR("Data projection have already been finished.",ERR,ERROR,*999)
            ELSE
              DATA_PROJECTION%DATA_PROJECTION_FINISHED=.TRUE.
              CALL DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE(DATA_PROJECTION,ERR,ERROR, &
                & *999) !<Initialise data projection part in data points
            ENDIF         
          ELSE
            CALL FLAG_ERROR("Data projection mesh is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection data points have not been finished.",ERR,ERROR,*999)
        ENDIF     
      ELSE
        CALL FLAG_ERROR("Data projection data points is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_CREATE_FINISH")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_CREATE_FINISH
  
  !
  !================================================================================================================================
  !  
  !>Starts the process of creating data projection.
  SUBROUTINE DATA_PROJECTION_CREATE_START_DATA_POINTS(DATA_PROJECTION_USER_NUMBER,DATA_POINTS,MESH,DATA_PROJECTION, &
    & ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: DATA_PROJECTION_USER_NUMBER !<the user number of the data projection
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points in which to create data projection
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh where the data points are projected
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<On exit, a pointer to the created data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DATA_PROJECTION_PTR_TYPE), ALLOCATABLE :: NEW_DATA_PROJECTIONS_PTR(:)
    INTEGER(INTG) :: DATA_POINTS_REGION_DIMENSIONS,MESH_REGION_DIMENSIONS,INSERT_STATUS
    INTEGER(INTG) :: xi_idx,data_projection_idx

    CALL ENTERS("DATA_PROJECTION_CREATE_START_DATA_POINTS",ERR,ERROR,*999)

    NULLIFY(DATA_PROJECTION)
    IF(ASSOCIATED(DATA_POINTS)) THEN !Has to be associated
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN !Has to be finished
        IF(ASSOCIATED(MESH)) THEN !Has to be associated
          IF(ASSOCIATED(DATA_PROJECTION)) THEN !Cannot be associated
            CALL FLAG_ERROR("Data projection is already associated.",ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(DATA_POINTS%REGION)) THEN
              DATA_POINTS_REGION_DIMENSIONS=DATA_POINTS%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            ELSEIF(ASSOCIATED(DATA_POINTS%INTERFACE)) THEN
              DATA_POINTS_REGION_DIMENSIONS=DATA_POINTS%INTERFACE%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            ELSE
              CALL FLAG_ERROR("Data points is not associated with a region or a interface.",ERR,ERROR,*999)
            ENDIF
            IF(ASSOCIATED(MESH%REGION)) THEN
              MESH_REGION_DIMENSIONS=MESH%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            ELSEIF(ASSOCIATED(MESH%INTERFACE)) THEN
              MESH_REGION_DIMENSIONS=MESH%INTERFACE%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
            ELSE
              CALL FLAG_ERROR("Mesh is not associated with a region or a interface.",ERR,ERROR,*999)
            ENDIF
            IF(DATA_POINTS_REGION_DIMENSIONS==MESH_REGION_DIMENSIONS) THEN !Dimension has to be equal
              ALLOCATE(DATA_PROJECTION,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data projection.",ERR,ERROR,*999)
              DATA_PROJECTION%USER_NUMBER=DATA_PROJECTION_USER_NUMBER
              DATA_PROJECTION%LABEL=""
              CALL TREE_ITEM_INSERT(DATA_POINTS%DATA_PROJECTIONS_TREE,DATA_PROJECTION_USER_NUMBER, &
                & DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS+1,INSERT_STATUS,ERR,ERROR,*999)
              DATA_PROJECTION%DATA_PROJECTION_FINISHED=.FALSE.
              DATA_PROJECTION%DATA_POINTS=>DATA_POINTS
              DATA_PROJECTION%MESH=>MESH
              DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS=DATA_POINTS_REGION_DIMENSIONS
              DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE=0.5_DP
              DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS=25                   
              !Default always project to boundaries faces/lines when mesh dimension is equal to region dimension. If mesh dimension is less, project to all elements            
              IF(MESH%NUMBER_OF_DIMENSIONS<DATA_POINTS_REGION_DIMENSIONS) THEN !mesh dimension < data dimension
                DATA_PROJECTION%NUMBER_OF_XI=MESH%NUMBER_OF_DIMENSIONS
                DATA_PROJECTION%PROJECTION_TYPE=DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
              ELSE
                SELECT CASE(MESH%NUMBER_OF_DIMENSIONS) !mesh dimension = data dimension
                CASE (2) 
                  DATA_PROJECTION%NUMBER_OF_XI=1
                  DATA_PROJECTION%PROJECTION_TYPE=DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE
                CASE (3)
                  DATA_PROJECTION%NUMBER_OF_XI=2
                  DATA_PROJECTION%PROJECTION_TYPE=DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE
                CASE DEFAULT
                  CALL FLAG_ERROR("Mesh dimensions out of bond [1,3].",ERR,ERROR,*999)
                END SELECT
              ENDIF
              SELECT CASE(DATA_PROJECTION%NUMBER_OF_XI) !mesh dimension = data dimension
              CASE (1)
                DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS=2
              CASE (2)
                DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS=4  
              CASE (3)
                DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS=8
              CASE DEFAULT
                CALL FLAG_ERROR("Mesh dimensions out of bond [1,3].",ERR,ERROR,*999)
              END SELECT 
              ALLOCATE(DATA_PROJECTION%STARTING_XI(DATA_PROJECTION%NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data projection starting xi.",ERR,ERROR,*999)
              DO xi_idx=1,DATA_PROJECTION%NUMBER_OF_XI
                DATA_PROJECTION%STARTING_XI(xi_idx)=0.5_DP !<initialised to 0.5 in each xi direction
              ENDDO !xi_idx              
              DATA_PROJECTION%ABSOLUTE_TOLERANCE=1.0E-8_DP
              DATA_PROJECTION%RELATIVE_TOLERANCE=1.0E-6_DP
              IF(DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS>0) THEN
                ALLOCATE(NEW_DATA_PROJECTIONS_PTR(DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS+1),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new data projections.",ERR,ERROR,*999)
                DO data_projection_idx=1,DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS
                  NEW_DATA_PROJECTIONS_PTR(data_projection_idx)%PTR=>DATA_POINTS%DATA_PROJECTIONS(data_projection_idx)%PTR
                ENDDO !xi_idx
              ELSE IF(DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS==0) THEN
                ALLOCATE(NEW_DATA_PROJECTIONS_PTR(DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS+1),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new data projections.",ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("The number of data projections is < 0.",ERR,ERROR,*999)
              ENDIF
              !Return the pointer
              NEW_DATA_PROJECTIONS_PTR(DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS+1)%PTR=>DATA_PROJECTION
              CALL MOVE_ALLOC(NEW_DATA_PROJECTIONS_PTR,DATA_POINTS%DATA_PROJECTIONS)
              DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS=DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS+1
              DATA_PROJECTION%GLOBAL_NUMBER=DATA_POINTS%NUMBER_OF_DATA_PROJECTIONS
            ELSE
              CALL FLAG_ERROR("Dimensions bewtween the mesh region/interface and data points region/interface does not match.", &
                & ERR,ERROR,*999)        
            ENDIF !DATA_POINTS_REGION_DIMENSIONS==MESH_REGION_DIMENSIONS
          ENDIF !ASSOCIATED(DATA_PROJECTION)
        ELSE
          CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
        ENDIF !ASSOCIATED(MESH)
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF !DATA_POINTS%DATA_POINTS_FINISHED
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF !ASSOCIATED(DATA_POINTS)
    
    CALL EXITS("DATA_PROJECTION_CREATE_START_DATA_POINTS")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_CREATE_START_DATA_POINTS",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_CREATE_START_DATA_POINTS")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_CREATE_START_DATA_POINTS 
  
  !
  !================================================================================================================================
  !

  !>Checks that a user data point number is defined for a specific data projection
  SUBROUTINE DataProjection_DataPointCheckExist(DataProjection,dataPointUserNumber,dataPointExist,dataPointGlobalNumber,err,error,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DataProjection !<A pointer to the data projection whose data points to check
    INTEGER(INTG) :: dataPointUserNumber !<The user data point number to check if it exists
    LOGICAL, INTENT(OUT) :: dataPointExist !<On exit, is .TRUE. if the data point user number exists in the region, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: dataPointGlobalNumber !<On exit, if the data point exists the global number corresponding to the user data point number. If the data point does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DATA_POINTS_TYPE), POINTER :: dataPoints
    TYPE(TREE_NODE_TYPE), POINTER :: treeNode
   
    CALL ENTERS("DataProjection_DataPointCheckExist",err,ERROR,*999)
    
    dataPointExist=.FALSE.
    dataPointGlobalNumber=0
    IF(ASSOCIATED(DataProjection)) THEN
      dataPoints=>DataProjection%DATA_POINTS
      IF(ASSOCIATED(dataPoints)) THEN
        NULLIFY(treeNode)
        CALL TREE_SEARCH(dataPoints%DATA_POINTS_TREE,dataPointUserNumber,treeNode,err,error,*999)
        IF(ASSOCIATED(treeNode)) THEN
          CALL TREE_NODE_VALUE_GET(dataPoints%DATA_POINTS_TREE,treeNode,dataPointGlobalNumber,err,error,*999)
          dataPointExist=.TRUE.
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",err,error,*999)
    ENDIF

    CALL EXITS("DataProjection_DataPointCheckExist")
    RETURN
999 CALL ERRORS("DataProjection_DataPointCheckExist",err,error)
    CALL EXITS("DataProjection_DataPointCheckExist")
    RETURN 1
   
  END SUBROUTINE DataProjection_DataPointCheckExist

  !
  !================================================================================================================================
  !
  
  !>Initialises the data projection part in a given data points. %%%%% THIS NEED TO BE CHANGED!!! TIM
  SUBROUTINE DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE(DATA_PROJECTION,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to initialise the data projection result for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_DATA_POINTS
    INTEGER(INTG) :: data_point_idx
    
    CALL ENTERS("DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE",ERR,ERROR,*999)

    !NULLIFY(PROJECTED_POINTS)
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(ASSOCIATED(DATA_PROJECTION%DATA_POINTS)) THEN
          IF(DATA_PROJECTION%DATA_POINTS%DATA_POINTS_FINISHED) THEN
            NUMBER_OF_DATA_POINTS = DATA_PROJECTION%DATA_POINTS%NUMBER_OF_DATA_POINTS
            ALLOCATE(DATA_PROJECTION%DATA_PROJECTION_RESULTS(NUMBER_OF_DATA_POINTS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data projection data projection results.",ERR,ERROR,*999)
            DO data_point_idx=1,NUMBER_OF_DATA_POINTS
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%USER_NUMBER=DATA_PROJECTION%DATA_POINTS%DATA_POINTS( &
                & data_point_idx)%USER_NUMBER
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%DISTANCE=0.0_DP
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_NUMBER=0
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_FACE_NUMBER=0
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_LINE_NUMBER=0
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
              ALLOCATE(DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI(DATA_PROJECTION%NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data projection data projection results "// &
                & "("//TRIM(NUMBER_TO_VSTRING (data_point_idx,"*",ERR,ERROR))//") xi.",ERR,ERROR,*999)
              DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI(1:DATA_PROJECTION%NUMBER_OF_XI)= &
                & DATA_PROJECTION%STARTING_XI(1:DATA_PROJECTION%NUMBER_OF_XI)
            ENDDO !data_point_idx
          ELSE
            CALL FLAG_ERROR("Data projection data points have not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection data points is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_DATA_PROJECTION_RESULT_INITIALISE 
  
  !
  !================================================================================================================================
  !

  !>Destroys a data projection.
  SUBROUTINE DATA_PROJECTION_DESTROY(DATA_PROJECTION,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to destroy
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
       
    CALL ENTERS("DATA_PROJECTION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      CALL DATA_PROJECTION_FINALISE(DATA_PROJECTION,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_DESTROY")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_DESTROY",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_DESTROY")
    RETURN 1
    
  END SUBROUTINE DATA_PROJECTION_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalise a data projection.
  SUBROUTINE DATA_PROJECTION_FINALISE(DATA_PROJECTION,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx
       
    CALL ENTERS("DATA_PROJECTION_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(ALLOCATED(DATA_PROJECTION%STARTING_XI)) THEN
        DEALLOCATE(DATA_PROJECTION%STARTING_XI)
      ENDIF
      IF(ALLOCATED(DATA_PROJECTION%DATA_PROJECTION_RESULTS)) THEN
        DO dataPointIdx=1,SIZE(DATA_PROJECTION%DATA_PROJECTION_RESULTS,1)
          DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%USER_NUMBER=0
          DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%DISTANCE=0
          DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%ELEMENT_NUMBER=0
          DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%ELEMENT_FACE_NUMBER=0
          DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%ELEMENT_LINE_NUMBER=0
          DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%EXIT_TAG=0
          IF(ALLOCATED(DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%XI)) THEN
            DEALLOCATE(DATA_PROJECTION%DATA_PROJECTION_RESULTS(dataPointIdx)%XI)
          ENDIF
        ENDDO !dataPointIdx
        DEALLOCATE(DATA_PROJECTION%DATA_PROJECTION_RESULTS)
      ENDIF
      DEALLOCATE(DATA_PROJECTION)
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_FINALISE")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_FINALISE",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_FINALISE")
    RETURN 1
    
  END SUBROUTINE DATA_PROJECTION_FINALISE
  
  !
  !================================================================================================================================
  !  

  !>Gets the user number for a data point identified by a given global number. 
  SUBROUTINE DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,USER_NUMBER,GLOBAL_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection whose data points to get the number for
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the data point to get the global number for
    INTEGER(INTG), INTENT(OUT) :: GLOBAL_NUMBER !<On exit, the global number of the specified user data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DATA_POINTS_TYPE), POINTER  :: DATA_POINTS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    CALL ENTERS("DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      DATA_POINTS=>DATA_PROJECTION%DATA_POINTS
      IF(ASSOCIATED(DATA_POINTS)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DATA_POINTS%DATA_POINTS_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DATA_POINTS%DATA_POINTS_TREE,TREE_NODE,GLOBAL_NUMBER,ERR,ERROR,*999)
        ELSE
          LOCAL_ERROR="Tree node is not associates (cannot find the user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR, &
            & ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_DATA_POINST_GLOBAL_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET")
    RETURN 1
   
  END SUBROUTINE DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET
  
  !
  !================================================================================================================================
  !


  !>Evaluates data projection.
  SUBROUTINE DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE(DATA_PROJECTION,PROJECTION_FIELD,ERR,ERROR,*) !optimising
    
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to evaluate
    TYPE(FIELD_TYPE), POINTER :: PROJECTION_FIELD !<A pointer to the projection field to evaluate, this would normally be geometric field or dependent field
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINTS(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    !TYPE(DATA_PROJECTION_PROJECTED_POINTS_TYPE), POINTER :: PROJECTED_POINTS
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE,NUMBER_COMPUTATIONAL_NODES !<computational node/rank of the current process    
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_DATA_POINTS
    INTEGER(INTG) :: NUMBER_OF_ELEMENTS,NUMBER_OF_FACES,NUMBER_OF_LINES,NUMBER_OF_CANDIDATES
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_CANDIDATES,TOTAL_NUMBER_OF_CLOSEST_CANDIDATES,REDUCED_NUMBER_OF_CLOSEST_CANDIDATES
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(:) 
    INTEGER(INTG), ALLOCATABLE :: CANDIDATE_ELEMENTS(:),CLOSEST_ELEMENTS(:,:),CANDIDATE_FACES(:),CLOSEST_FACES(:,:)
    !INTEGER(INTG) :: NUMBER_OF_CANDIDATE_LINES
    REAL(DP), ALLOCATABLE :: CLOSEST_DISTANCES(:,:),GLOBAL_CLOSEST_DISTANCES(:,:)
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES(:)
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_MPI_DISPLACEMENTS(:),SORTING_IND_1(:),SORTING_IND_2(:)
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_NUMBER_OF_PROJECTED_POINTS(:)
    INTEGER(INTG) :: MPI_CLOSEST_DISTANCES,DATA_PROJECTION_GLOBAL_NUMBER
    INTEGER(INTG) :: MPI_IERROR
    !LOGICAL :: ELEMENT_FOUND
    !LOGICAL, ALLOCATABLE :: PROJECTION_CONVERGED(:)
    INTEGER(INTG), ALLOCATABLE :: PROJECTED_ELEMENT(:),PROJECTED_FACE(:),PROJECTION_EXIT_TAG(:)
    REAL(DP), ALLOCATABLE :: PROJECTED_DISTANCE(:,:),PROJECTED_XI(:,:)
    
    INTEGER(INTG) :: ne,nse,ncn,ni,localElementNumber
    INTEGER(INTG) :: temp_number,start_idx,finish_idx,data_point_idx
    
    LOGICAL :: BOUNDARY_PROJECTION,elementExists,ghostElement    
    
    !INTEGER(INTG) :: NUMBER_OF_PROJECTED_POINTS,NUMBER_OF_PROJECTED_POINTS_2
  
    CALL ENTERS("DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE",ERR,ERROR,*999)
    
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATED_POINTS)
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(ASSOCIATED(PROJECTION_FIELD)) THEN
          IF(PROJECTION_FIELD%FIELD_FINISHED) THEN
            IF(ASSOCIATED(DATA_PROJECTION%MESH,PROJECTION_FIELD%DECOMPOSITION%MESH)) THEN
              DATA_PROJECTION%PROJECTION_FIELD=>PROJECTION_FIELD
              DATA_PROJECTION_GLOBAL_NUMBER=DATA_PROJECTION%GLOBAL_NUMBER
              DATA_POINTS=>DATA_PROJECTION%DATA_POINTS
              CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(PROJECTION_FIELD,INTERPOLATION_PARAMETERS,ERR,ERROR,*998, &
                & FIELD_GEOMETRIC_COMPONENTS_TYPE)
              CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINTS,ERR,ERROR,*998, &
                & FIELD_GEOMETRIC_COMPONENTS_TYPE)
              INTERPOLATED_POINT=>INTERPOLATED_POINTS(FIELD_U_VARIABLE_TYPE)%PTR
              DECOMPOSITION=>PROJECTION_FIELD%DECOMPOSITION  
              MESH_COMPONENT_NUMBER=PROJECTION_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
              DOMAIN=>PROJECTION_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
              NUMBER_OF_DATA_POINTS=DATA_POINTS%NUMBER_OF_DATA_POINTS
              MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
              NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
              BOUNDARY_PROJECTION=(DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE).OR.( &
                & DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)          
              !#####################################################################################################################
              !find elements/faces/lines inside current computational node, get the boundary faces/lines only if asked
              !the elements/faces/lines are required to perform projection of points in the current computational node
              !the are all pre-allocated to the maximum array length (i.e. NUMBER_OF_ELEMENTS), but only up to the NUMBER_OF_CANDIDATESth index are assigned
              NUMBER_OF_ELEMENTS=DOMAIN%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
              NUMBER_OF_FACES=DOMAIN%TOPOLOGY%FACES%NUMBER_OF_FACES
              NUMBER_OF_LINES=DOMAIN%TOPOLOGY%LINES%NUMBER_OF_LINES        
              NUMBER_OF_CANDIDATES=0
              IF(ALLOCATED(DATA_PROJECTION%candidateElementNumbers) .AND. ALLOCATED(DATA_PROJECTION%localFaceLineNumbers)) THEN
                SELECT CASE(DATA_PROJECTION%PROJECTION_TYPE)
                CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)  !identify all non-ghost boundary lines
                  ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_LINES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
                  ALLOCATE(CANDIDATE_FACES(NUMBER_OF_LINES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate lines.",ERR,ERROR,*999)
                  DO ne=1,SIZE(DATA_PROJECTION%candidateElementNumbers,1) !Loop through all candidate element defined by user number
                    CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(DECOMPOSITION%TOPOLOGY,DATA_PROJECTION% &
                      & candidateElementNumbers(ne),elementExists,localElementNumber,ghostElement,ERR,ERROR,*999) !Check if element exists on current domain, get local number
                    IF((elementExists) .AND. (.NOT.ghostElement)) THEN !Get non-ghost elements
                      NUMBER_OF_CANDIDATES=NUMBER_OF_CANDIDATES+1
                      CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATES)=localElementNumber
                      CANDIDATE_FACES(NUMBER_OF_CANDIDATES)=DATA_PROJECTION%localFaceLineNumbers(ne) !Store element line number for line projection type
                    ENDIF
                  ENDDO !ne
                CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) !identify all non-ghost boundary faces
                  IF(DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS>=2) THEN
                    ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_FACES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
                    ALLOCATE(CANDIDATE_FACES(NUMBER_OF_FACES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate faces.",ERR,ERROR,*999)
                    DO ne=1,SIZE(DATA_PROJECTION%candidateElementNumbers,1) !Loop through all candidate element defined by user number
                      CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(DECOMPOSITION%TOPOLOGY,DATA_PROJECTION% &
                        & candidateElementNumbers(ne),elementExists,localElementNumber,ghostElement,ERR,ERROR,*999) !Check if element exists on current domain, get local number
                      IF((elementExists) .AND. (.NOT.ghostElement)) THEN !Get non-ghost elements
                        NUMBER_OF_CANDIDATES=NUMBER_OF_CANDIDATES+1
                        CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATES)=localElementNumber
                        CANDIDATE_FACES(NUMBER_OF_CANDIDATES)=DATA_PROJECTION%localFaceLineNumbers(ne) !Store element face number for face projection type
                      ENDIF
                    ENDDO !ne
                  ELSE
                    CALL FLAG_ERROR("Decomposition mesh number of dimensions has to be 2 or greater.",ERR,ERROR,*999)        
                  ENDIF
                CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE) !identify all non-ghost elements
                  IF(DATA_PROJECTION%NUMBER_OF_XI==DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS) THEN
                    ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_ELEMENTS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
                    DO ne=1,SIZE(DATA_PROJECTION%candidateElementNumbers,1) !Loop through all candidate element defined by user number
                      CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(DECOMPOSITION%TOPOLOGY,DATA_PROJECTION% &
                        & candidateElementNumbers(ne),elementExists,localElementNumber,ghostElement,ERR,ERROR,*999) !Check if element exists on current domain, get local number
                      IF((elementExists) .AND. (.NOT.ghostElement)) THEN !Get non-ghost elements
                        NUMBER_OF_CANDIDATES=NUMBER_OF_CANDIDATES+1
                        CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATES)=localElementNumber
                      ENDIF
                    ENDDO !ne
                  ELSE
                    CALL FLAG_ERROR("Data projection number of xi has to equal to decomposition mesh number of dimensions",ERR, &
                      & ERROR,*999)
                  ENDIF
                CASE DEFAULT
                  CALL FLAG_ERROR("No match for data projection type found",ERR,ERROR,*999)
                END SELECT !DATA_PROJECTION%PROJECTION_TYPE           
              ELSE !If user didn't define candidate element number
                SELECT CASE(DATA_PROJECTION%PROJECTION_TYPE)
                  CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)  !identify all non-ghost boundary lines
                    ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_LINES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
                    ALLOCATE(CANDIDATE_FACES(NUMBER_OF_LINES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate lines.",ERR,ERROR,*999)
                    DO ne=1,DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
                      IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BOUNDARY_ELEMENT) THEN
                        DO nse=1,SIZE(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_LINES,1)
                          temp_number=DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_LINES(nse)
                          IF(DECOMPOSITION%TOPOLOGY%LINES%LINES(temp_number)%BOUNDARY_LINE) THEN
                            NUMBER_OF_CANDIDATES=NUMBER_OF_CANDIDATES+1
                            CANDIDATE_FACES(NUMBER_OF_CANDIDATES)=nse
                            CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATES)=ne
                          ENDIF !boundary line
                        ENDDO !nse
                      ENDIF !boundary element
                    ENDDO !ne
                  CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) !identify all non-ghost boundary faces
                    IF(DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS>=2) THEN
                      ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_FACES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
                      ALLOCATE(CANDIDATE_FACES(NUMBER_OF_FACES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate faces.",ERR,ERROR,*999)
                      DO ne=1,DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
                        IF(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BOUNDARY_ELEMENT) THEN
                          DO nse=1,SIZE(DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_FACES,1)
                            temp_number=DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_FACES(nse)
                            IF(DECOMPOSITION%TOPOLOGY%FACES%FACES(temp_number)%BOUNDARY_FACE) THEN
                              NUMBER_OF_CANDIDATES=NUMBER_OF_CANDIDATES+1
                              CANDIDATE_FACES(NUMBER_OF_CANDIDATES)=nse
                              CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATES)=ne
                            ENDIF !boundary face
                          ENDDO !nse
                        ENDIF !boundary element
                      ENDDO !ne
                    ELSE
                      CALL FLAG_ERROR("Decomposition mesh number of dimensions has to be 2 or greater.",ERR,ERROR,*999)        
                    ENDIF
                  CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE) !identify all non-ghost elements
                    IF(DATA_PROJECTION%NUMBER_OF_XI==DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS) THEN
                      ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_ELEMENTS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
                      DO ne=1,DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
                        NUMBER_OF_CANDIDATES=NUMBER_OF_CANDIDATES+1
                        CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATES)=ne
                      ENDDO
                    ELSE
                      CALL FLAG_ERROR("Data projection number of xi has to equal to decomposition mesh number of dimensions",ERR, &
                        & ERROR,*999)
                    ENDIF
                  CASE DEFAULT
                    CALL FLAG_ERROR("No match for data projection type found",ERR,ERROR,*999)
                END SELECT !DATA_PROJECTION%PROJECTION_TYPE
              ENDIF
              !#####################################################################################################################
              !find the clostest elements/faces/lines for each point in the current computational node base on starting xi
              !the clostest elements/faces/lines are required to shrink down on the list of possible projection candiates
              NUMBER_OF_CLOSEST_CANDIDATES=MIN(DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS,NUMBER_OF_CANDIDATES)
              ALLOCATE(CLOSEST_ELEMENTS(NUMBER_OF_DATA_POINTS,NUMBER_OF_CLOSEST_CANDIDATES),STAT=ERR)!the information for each data point has to be stored in the corresponding rows for them to be contiguous in memory for easy MPI access
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate closest elements.",ERR,ERROR,*999)
              IF(BOUNDARY_PROJECTION) THEN
                ALLOCATE(CLOSEST_FACES(NUMBER_OF_DATA_POINTS,NUMBER_OF_CLOSEST_CANDIDATES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate closest sub element.",ERR,ERROR,*999)          
              ENDIF
              ALLOCATE(CLOSEST_DISTANCES(NUMBER_OF_DATA_POINTS,NUMBER_OF_CLOSEST_CANDIDATES),STAT=ERR)!the information for each data point has to be stored in the corresponding rows for them to be contiguous in memory for easy MPI access
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate closest distances.",ERR,ERROR,*999) 
              SELECT CASE(DATA_PROJECTION%PROJECTION_TYPE)
                CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) !find closest candidate lines
                  DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                    CALL DATA_PROJECTION_CLOSEST_LINES_FIND(DATA_PROJECTION,INTERPOLATED_POINT, &
                      & DATA_POINTS%DATA_POINTS(data_point_idx)%position,CANDIDATE_ELEMENTS, &
                      & CANDIDATE_FACES,NUMBER_OF_CANDIDATES,CLOSEST_ELEMENTS(data_point_idx,:),CLOSEST_FACES( &
                      & data_point_idx,:),CLOSEST_DISTANCES(data_point_idx,:),ERR,ERROR,*999)
                  ENDDO !data_point_idx
                CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) !find closest candidate faces      
                  DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                    CALL DATA_PROJECTION_CLOSEST_FACES_FIND(DATA_PROJECTION,INTERPOLATED_POINT, &
                      & DATA_POINTS%DATA_POINTS(data_point_idx)%position,CANDIDATE_ELEMENTS, &
                      & CANDIDATE_FACES,NUMBER_OF_CANDIDATES,CLOSEST_ELEMENTS(data_point_idx,:),CLOSEST_FACES( &
                      & data_point_idx,:),CLOSEST_DISTANCES(data_point_idx,:),ERR,ERROR,*999)
                  ENDDO !data_point_idx
                CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE) !find closest candidate elements
                  DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                    CALL DATA_PROJECTION_CLOSEST_ELEMENTS_FIND(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS%DATA_POINTS( &
                      & data_point_idx)%position,CANDIDATE_ELEMENTS,NUMBER_OF_CANDIDATES,CLOSEST_ELEMENTS(data_point_idx,:), &
                      & CLOSEST_DISTANCES(data_point_idx,:),ERR,ERROR,*999)
                    !CLOSEST_ELEMENTS(data_point_idx,:)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(CLOSEST_ELEMENTS(data_point_idx,:)) !local to global element number mapping
                  ENDDO !data_point_idx
                CASE DEFAULT
                  CALL FLAG_ERROR("No match for data projection type found",ERR,ERROR,*999)
              END SELECT
              !#####################################################################################################################
              !Newton project data point to the list of closest elements, faces or lines
              !project the data points to each of the closest elements, use MPI if number of computational nodes is greater than 1
              IF(NUMBER_COMPUTATIONAL_NODES>1) THEN !use mpi
                !allocate arrays for mpi communication
                ALLOCATE(GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(NUMBER_OF_DATA_POINTS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global to local number of closest elements.",ERR,ERROR,*999)
                ALLOCATE(GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES(NUMBER_COMPUTATIONAL_NODES),STAT=ERR) 
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all number of closest candidates.",ERR,ERROR,*999)
                ALLOCATE(GLOBAL_MPI_DISPLACEMENTS(NUMBER_COMPUTATIONAL_NODES),STAT=ERR) 
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all displacements.",ERR,ERROR,*999)
                ALLOCATE(GLOBAL_NUMBER_OF_PROJECTED_POINTS(NUMBER_COMPUTATIONAL_NODES),STAT=ERR) 
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all number of projected points.",ERR,ERROR,*999)
                ALLOCATE(PROJECTION_EXIT_TAG(NUMBER_OF_DATA_POINTS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected.",ERR,ERROR,*999)
                ALLOCATE(PROJECTED_ELEMENT(NUMBER_OF_DATA_POINTS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected element.",ERR,ERROR,*999)
                IF(BOUNDARY_PROJECTION) THEN
                  ALLOCATE(PROJECTED_FACE(NUMBER_OF_DATA_POINTS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected sub element.",ERR,ERROR,*999)
                ENDIF
                ALLOCATE(PROJECTED_DISTANCE(2,NUMBER_OF_DATA_POINTS),STAT=ERR) !PROJECTED_DISTANCE(2,:) stores the compuational node number, the information for each data point has to be stored in the corresponding column for MPI_ALLREDUCE with location return to work
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected distance.",ERR,ERROR,*999)
                ALLOCATE(PROJECTED_XI(DATA_PROJECTION%NUMBER_OF_XI,NUMBER_OF_DATA_POINTS),STAT=ERR) !the information for each data point is stored in the corresponding column to be consistent with PROJECTED_DISTANCE
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected.",ERR,ERROR,*999)
                ALLOCATE(SORTING_IND_2(NUMBER_OF_DATA_POINTS),STAT=ERR) 
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sorting ind 2.",ERR,ERROR,*999)          
                !gather and distribute the number of closest elements from all computational nodes
                CALL MPI_ALLGATHER(NUMBER_OF_CLOSEST_CANDIDATES,1,MPI_INTEGER,GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES,1,MPI_INTEGER, &
                  & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
                TOTAL_NUMBER_OF_CLOSEST_CANDIDATES=SUM(GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES,1) !sum of all number of closest candidates from all computational nodes
                !allocate arrays to store information gathered from all computational node
                ALLOCATE(GLOBAL_CLOSEST_DISTANCES(NUMBER_OF_DATA_POINTS,TOTAL_NUMBER_OF_CLOSEST_CANDIDATES),STAT=ERR) !the information for each data point is stored in the corresponding row so they are contiguous in memory for easy MPI access
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all closest distances.",ERR,ERROR,*999)
                ALLOCATE(SORTING_IND_1(TOTAL_NUMBER_OF_CLOSEST_CANDIDATES),STAT=ERR) 
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sorting ind 1.",ERR,ERROR,*999)
                !MPI:create and commit MPI_TYPE_CONTIGUOUS      
                CALL MPI_TYPE_CONTIGUOUS(NUMBER_OF_DATA_POINTS,MPI_DOUBLE_PRECISION,MPI_CLOSEST_DISTANCES,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_TYPE_CONTIGUOUS",MPI_IERROR,ERR,ERROR,*999)        
                CALL MPI_TYPE_COMMIT(MPI_CLOSEST_DISTANCES,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPI_IERROR,ERR,ERROR,*999)
                !create displacement vectors for MPI_ALLGATHERV
                GLOBAL_MPI_DISPLACEMENTS(1)=0
                DO ncn=1,(NUMBER_COMPUTATIONAL_NODES-1)
                  GLOBAL_MPI_DISPLACEMENTS(ncn+1)=GLOBAL_MPI_DISPLACEMENTS(ncn)+GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES(ncn)
                ENDDO
                !MPI:shares closest element distances between all domains
                CALL MPI_ALLGATHERV(CLOSEST_DISTANCES(1,1),NUMBER_OF_CLOSEST_CANDIDATES,MPI_CLOSEST_DISTANCES, &
                  & GLOBAL_CLOSEST_DISTANCES,GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES,GLOBAL_MPI_DISPLACEMENTS,MPI_CLOSEST_DISTANCES, &
                  & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
                REDUCED_NUMBER_OF_CLOSEST_CANDIDATES=MIN(DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS, &
                  & TOTAL_NUMBER_OF_CLOSEST_CANDIDATES)
                PROJECTED_DISTANCE(2,:)=MY_COMPUTATIONAL_NODE
                ! find the globally closest distances in the current domain
                DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                  CALL BUBBLE_ISORT(GLOBAL_CLOSEST_DISTANCES(data_point_idx,:),SORTING_IND_1,ERR,ERROR,*999)
                  SORTING_IND_1(1:TOTAL_NUMBER_OF_CLOSEST_CANDIDATES)=SORTING_IND_1(1:TOTAL_NUMBER_OF_CLOSEST_CANDIDATES)- &
                    & GLOBAL_MPI_DISPLACEMENTS(MY_COMPUTATIONAL_NODE+1) !shift the index to current computational node
                  GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)=0
                  DO ne=1,REDUCED_NUMBER_OF_CLOSEST_CANDIDATES
                    !sorted index indicates it is in the current computational domain
                    IF((SORTING_IND_1(ne)>=1).AND.(SORTING_IND_1(ne)<=GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES(MY_COMPUTATIONAL_NODE &
                      & +1)))GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)= &
                      & GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)+1
                  ENDDO
                  PROJECTED_DISTANCE(1,data_point_idx)=GLOBAL_CLOSEST_DISTANCES(data_point_idx,TOTAL_NUMBER_OF_CLOSEST_CANDIDATES) !assign initial distance to something large                           
                ENDDO
                SELECT CASE(DATA_PROJECTION%PROJECTION_TYPE)
                  CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) !Newton project to closest lines, and find miminum projection
                    DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                      NUMBER_OF_CLOSEST_CANDIDATES=GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)
                      IF(NUMBER_OF_CLOSEST_CANDIDATES>0) THEN 
                        CALL DATA_PROJECTION_NEWTON_LINES_EVALUATE(DATA_PROJECTION,INTERPOLATED_POINT, &
                          & DATA_POINTS%DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS( &
                          & data_point_idx,1:NUMBER_OF_CLOSEST_CANDIDATES),CLOSEST_FACES(data_point_idx,1: &
                          & NUMBER_OF_CLOSEST_CANDIDATES),PROJECTION_EXIT_TAG(data_point_idx),PROJECTED_ELEMENT(data_point_idx),  &
                          & PROJECTED_FACE(data_point_idx),PROJECTED_DISTANCE(1,data_point_idx),PROJECTED_XI(:,data_point_idx), &
                          & ERR,ERROR,*999)
                        PROJECTED_ELEMENT(data_point_idx)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(PROJECTED_ELEMENT( &
                          & data_point_idx)) !map the element number to global number
                      ENDIF
                    ENDDO
                  CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) !find closest candidate faces
                    DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                      NUMBER_OF_CLOSEST_CANDIDATES=GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)
                      IF(NUMBER_OF_CLOSEST_CANDIDATES>0) THEN 
                        CALL DATA_PROJECTION_NEWTON_FACES_EVALUATE(DATA_PROJECTION,INTERPOLATED_POINT, &
                          & DATA_POINTS%DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS( &
                          & data_point_idx,1:NUMBER_OF_CLOSEST_CANDIDATES),CLOSEST_FACES(data_point_idx, &
                          & 1:NUMBER_OF_CLOSEST_CANDIDATES),PROJECTION_EXIT_TAG(data_point_idx),PROJECTED_ELEMENT(data_point_idx), &
                          & PROJECTED_FACE(data_point_idx),PROJECTED_DISTANCE(1,data_point_idx),PROJECTED_XI(:,data_point_idx), &
                          & ERR,ERROR,*999)
                        PROJECTED_ELEMENT(data_point_idx)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(PROJECTED_ELEMENT( &
                          & data_point_idx)) !map the element number to global number
                      ENDIF
                    ENDDO
                  CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE) !find closest candidate elements
                    SELECT CASE(DATA_PROJECTION%NUMBER_OF_XI)
                      CASE (1) !1D element
                        DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                          NUMBER_OF_CLOSEST_CANDIDATES=GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)
                          IF(NUMBER_OF_CLOSEST_CANDIDATES>0) THEN 
                            CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS% &
                              & DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx, &
                              & 1:NUMBER_OF_CLOSEST_CANDIDATES),PROJECTION_EXIT_TAG(data_point_idx), &
                              & PROJECTED_ELEMENT(data_point_idx),PROJECTED_DISTANCE(1,data_point_idx), &
                              & PROJECTED_XI(:,data_point_idx),ERR,ERROR,*999)
                            PROJECTED_ELEMENT(data_point_idx)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(PROJECTED_ELEMENT( &
                              & data_point_idx)) !map the element number to global number

                          ENDIF
                        ENDDO
                      CASE (2) !2D element
                        DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                          NUMBER_OF_CLOSEST_CANDIDATES=GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)
                          IF(NUMBER_OF_CLOSEST_CANDIDATES>0) THEN 
                            CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS% &
                              & DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx, &
                              & 1:NUMBER_OF_CLOSEST_CANDIDATES),PROJECTION_EXIT_TAG(data_point_idx), &
                              & PROJECTED_ELEMENT(data_point_idx),PROJECTED_DISTANCE(1,data_point_idx), &
                              & PROJECTED_XI(:,data_point_idx),ERR,ERROR,*999)                    
                            PROJECTED_ELEMENT(data_point_idx)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(PROJECTED_ELEMENT( &
                              & data_point_idx)) !map the element number to global number
                          ENDIF
                        ENDDO
                      CASE (3) !3D element
                        DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                          NUMBER_OF_CLOSEST_CANDIDATES=GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES(data_point_idx)
                          IF(NUMBER_OF_CLOSEST_CANDIDATES>0) THEN 
                            CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS% &
                              & DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx, &
                              & 1:NUMBER_OF_CLOSEST_CANDIDATES),PROJECTION_EXIT_TAG(data_point_idx), &
                              & PROJECTED_ELEMENT(data_point_idx),PROJECTED_DISTANCE(1,data_point_idx), &
                              & PROJECTED_XI(:,data_point_idx),ERR,ERROR,*999)
                            PROJECTED_ELEMENT(data_point_idx)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(PROJECTED_ELEMENT( &
                              & data_point_idx)) !map the element number to global number
                          ENDIF
                        ENDDO
                      CASE DEFAULT
                        CALL FLAG_ERROR("Data projection number of xi is invalid",ERR,ERROR,*999)
                    END SELECT
                  CASE DEFAULT
                    CALL FLAG_ERROR("No match for data projection type found",ERR,ERROR,*999)
                END SELECT          
                !MPI:find the shortest projected distance in all domains
                CALL MPI_ALLREDUCE(MPI_IN_PLACE,PROJECTED_DISTANCE,NUMBER_OF_DATA_POINTS,MPI_2DOUBLE_PRECISION,MPI_MINLOC, &
                  & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                !sort the computational node/rank from 0 to number of computational node
                CALL BUBBLE_ISORT(PROJECTED_DISTANCE(2,:),SORTING_IND_2,ERR,ERROR,*999)
                DO ncn=0,(NUMBER_COMPUTATIONAL_NODES-1)
                  GLOBAL_NUMBER_OF_PROJECTED_POINTS(ncn+1)=COUNT(ABS(PROJECTED_DISTANCE(2,:)-REAL(ncn))<ZERO_TOLERANCE)
                ENDDO !ncn
                start_idx=SUM(GLOBAL_NUMBER_OF_PROJECTED_POINTS(1:MY_COMPUTATIONAL_NODE))+1
                finish_idx=start_idx+GLOBAL_NUMBER_OF_PROJECTED_POINTS(MY_COMPUTATIONAL_NODE+1)-1
                !create displacement vectors for MPI_ALLGATHERV          
                DO ncn=1,(NUMBER_COMPUTATIONAL_NODES-1)
                  GLOBAL_MPI_DISPLACEMENTS(ncn+1)=GLOBAL_MPI_DISPLACEMENTS(ncn)+GLOBAL_NUMBER_OF_PROJECTED_POINTS(ncn)
                ENDDO !ncn  
                !MPI:shares minimum projection informationg bewteen all domains
                CALL MPI_ALLGATHERV(PROJECTED_ELEMENT(SORTING_IND_2(start_idx:finish_idx)),GLOBAL_NUMBER_OF_PROJECTED_POINTS( &
                  & MY_COMPUTATIONAL_NODE+1),MPI_INTEGER,PROJECTED_ELEMENT,GLOBAL_NUMBER_OF_PROJECTED_POINTS, &
                  & GLOBAL_MPI_DISPLACEMENTS,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR) !PROJECTED_ELEMENT
                CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
                IF(BOUNDARY_PROJECTION) THEN
                  CALL MPI_ALLGATHERV(PROJECTED_FACE(SORTING_IND_2(start_idx:finish_idx)),GLOBAL_NUMBER_OF_PROJECTED_POINTS( &
                    & MY_COMPUTATIONAL_NODE+1),MPI_INTEGER,PROJECTED_FACE,GLOBAL_NUMBER_OF_PROJECTED_POINTS, &
                    & GLOBAL_MPI_DISPLACEMENTS,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR) !PROJECTED_FACE
                  CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999) 
                ENDIF
                DO ni=1,DATA_PROJECTION%NUMBER_OF_XI
                  CALL MPI_ALLGATHERV(PROJECTED_XI(ni,SORTING_IND_2(start_idx:finish_idx)),GLOBAL_NUMBER_OF_PROJECTED_POINTS( &
                    & MY_COMPUTATIONAL_NODE+1),MPI_DOUBLE_PRECISION,PROJECTED_XI(ni,:),GLOBAL_NUMBER_OF_PROJECTED_POINTS, &
                    & GLOBAL_MPI_DISPLACEMENTS,MPI_DOUBLE_PRECISION,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR) !PROJECTED_XI
                  CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
                ENDDO                      
                CALL MPI_ALLGATHERV(PROJECTION_EXIT_TAG(SORTING_IND_2(start_idx:finish_idx)),GLOBAL_NUMBER_OF_PROJECTED_POINTS( &
                  & MY_COMPUTATIONAL_NODE+1),MPI_INTEGER,PROJECTION_EXIT_TAG,GLOBAL_NUMBER_OF_PROJECTED_POINTS, &
                  & GLOBAL_MPI_DISPLACEMENTS,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR) !PROJECTION_EXIT_TAG
                CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
                !assign projection information to projected points
                DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                  DATA_PROJECTION%DATA_PROJECTION_RESULTS(SORTING_IND_2(data_point_idx))%EXIT_TAG=PROJECTION_EXIT_TAG( &
                    & data_point_idx)
                  DATA_PROJECTION%DATA_PROJECTION_RESULTS(SORTING_IND_2(data_point_idx))%ELEMENT_NUMBER=PROJECTED_ELEMENT( &
                    & data_point_idx)
                  DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%DISTANCE=PROJECTED_DISTANCE(1,data_point_idx)
                  DATA_PROJECTION%DATA_PROJECTION_RESULTS(SORTING_IND_2(data_point_idx))%XI(1:DATA_PROJECTION%NUMBER_OF_XI)= &
                    & PROJECTED_XI(1:DATA_PROJECTION%NUMBER_OF_XI,data_point_idx)
                ENDDO !data_point_idx
                PROJECTED_XI(:,SORTING_IND_2)=PROJECTED_XI
                PROJECTED_ELEMENT(SORTING_IND_2)=PROJECTED_ELEMENT       
                IF(DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) THEN
                  DO data_point_idx=1,NUMBER_OF_DATA_POINTS          
                    DATA_PROJECTION%DATA_PROJECTION_RESULTS(SORTING_IND_2(data_point_idx))%ELEMENT_LINE_NUMBER=PROJECTED_FACE( &
                      & data_point_idx)
                  ENDDO !data_point_idx
                ELSEIF(DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) THEN
                  DO data_point_idx=1,NUMBER_OF_DATA_POINTS          
                    DATA_PROJECTION%DATA_PROJECTION_RESULTS(SORTING_IND_2(data_point_idx))%ELEMENT_FACE_NUMBER=PROJECTED_FACE( &
                      & data_point_idx)
                  ENDDO !data_point_idx
                ENDIF            
              ELSE !no need to use mpi
                SELECT CASE(DATA_PROJECTION%PROJECTION_TYPE)
                  CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) !Newton project to closest lines, and find miminum projection
                    DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                      CALL DATA_PROJECTION_NEWTON_LINES_EVALUATE(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS%DATA_POINTS( &
                        & data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx,:),CLOSEST_FACES(data_point_idx,:), &
                        & DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%EXIT_TAG,DATA_PROJECTION% &
                        & DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_NUMBER,DATA_PROJECTION% &
                        & DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_LINE_NUMBER,DATA_PROJECTION%DATA_PROJECTION_RESULTS( &
                        & data_point_idx)%DISTANCE,DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI,ERR,ERROR,*999)
                    ENDDO
                  CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) !find closest candidate faces
                    DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                      CALL DATA_PROJECTION_NEWTON_FACES_EVALUATE(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS%DATA_POINTS( &
                        & data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx,:),CLOSEST_FACES(data_point_idx,:), &
                        & DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%EXIT_TAG,DATA_PROJECTION% &
                        & DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_NUMBER,DATA_PROJECTION% &
                        & DATA_PROJECTION_RESULTS(data_point_idx)%ELEMENT_FACE_NUMBER,DATA_PROJECTION%DATA_PROJECTION_RESULTS( &
                        & data_point_idx)%DISTANCE,DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI,ERR,ERROR,*999)
                    ENDDO
                  CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE) !find closest candidate elements        
                    SELECT CASE(DATA_PROJECTION%NUMBER_OF_XI)
                      CASE (1) !1D mesh
                        DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                          CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS% &
                            & DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx,:),DATA_PROJECTION% &
                            & DATA_PROJECTION_RESULTS(data_point_idx)%EXIT_TAG,DATA_PROJECTION%DATA_PROJECTION_RESULTS( &
                            & data_point_idx)%ELEMENT_NUMBER, DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%DISTANCE, &
                            & DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI,ERR,ERROR,*999)
                        ENDDO
                      CASE (2) !2D mesh
                        DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                          CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS% &
                            & DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx,:),DATA_PROJECTION% &
                            & DATA_PROJECTION_RESULTS(data_point_idx)%EXIT_TAG,DATA_PROJECTION%DATA_PROJECTION_RESULTS( &
                            & data_point_idx)%ELEMENT_NUMBER,DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%DISTANCE, &
                            & DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI,ERR,ERROR,*999)
                        ENDDO
                      CASE (3) !3D mesh
                        DO data_point_idx=1,NUMBER_OF_DATA_POINTS
                          CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINTS% &
                            & DATA_POINTS(data_point_idx)%position,CLOSEST_ELEMENTS(data_point_idx,:),DATA_PROJECTION% &
                            & DATA_PROJECTION_RESULTS(data_point_idx)%EXIT_TAG,DATA_PROJECTION%DATA_PROJECTION_RESULTS( &
                            & data_point_idx)%ELEMENT_NUMBER,DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%DISTANCE, &
                            & DATA_PROJECTION%DATA_PROJECTION_RESULTS(data_point_idx)%XI,ERR,ERROR,*999)
                        ENDDO
                      CASE DEFAULT
                        CALL FLAG_ERROR("Data projection number of xi is invalid",ERR,ERROR,*999)
                    END SELECT !DATA_PROJECTION%NUMBER_OF_XI
                  CASE DEFAULT
                    CALL FLAG_ERROR("No match for data projection type found",ERR,ERROR,*999)
                END SELECT                  
              ENDIF !NUMBER_COMPUTATIONAL_NODES>1
              DATA_PROJECTION%DATA_PROJECTION_PROJECTED=.TRUE.
            ELSE
              CALL FLAG_ERROR("Data projection and projection field are not sharing the same mesh.",ERR,ERROR,*999)
            ENDIF !ASSOCIATED(DATA_PROJECTION%MESH,PROJECTION_FIELD%DECOMPOSITION%MESH)
          ELSE
            CALL FLAG_ERROR("Projection field have not been finished.",ERR,ERROR,*999)
          ENDIF !PROJECTION_FIELD%FIELD_FINISHED
        ELSE
          CALL FLAG_ERROR("Projection field is not associated.",ERR,ERROR,*999)
        ENDIF !ASSOCIATED(PROJECTION_FIELD)  
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF !DATA_PROJECTION%DATA_PROJECTION_FINISHED
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF !ASSOCIATED(DATA_PROJECTION)  

    ! Deallocate arrays used within this routine
    IF(ALLOCATED(CANDIDATE_ELEMENTS)) DEALLOCATE(CANDIDATE_ELEMENTS)
    IF(ALLOCATED(CANDIDATE_FACES)) DEALLOCATE(CANDIDATE_FACES)
    IF(ALLOCATED(CLOSEST_ELEMENTS)) DEALLOCATE(CLOSEST_ELEMENTS)
    IF(ALLOCATED(CLOSEST_FACES)) DEALLOCATE(CLOSEST_FACES)
    IF(ALLOCATED(CLOSEST_DISTANCES)) DEALLOCATE(CLOSEST_DISTANCES)
    IF(ALLOCATED(GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES)) DEALLOCATE(GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES)
    IF(ALLOCATED(GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES)) DEALLOCATE(GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES)
    IF(ALLOCATED(GLOBAL_MPI_DISPLACEMENTS)) DEALLOCATE(GLOBAL_MPI_DISPLACEMENTS)
    IF(ALLOCATED(GLOBAL_NUMBER_OF_PROJECTED_POINTS)) DEALLOCATE(GLOBAL_NUMBER_OF_PROJECTED_POINTS)
    IF(ALLOCATED(PROJECTION_EXIT_TAG)) DEALLOCATE(PROJECTION_EXIT_TAG)
    IF(ALLOCATED(PROJECTED_ELEMENT)) DEALLOCATE(PROJECTED_ELEMENT)
    IF(ALLOCATED(PROJECTED_FACE)) DEALLOCATE(PROJECTED_FACE)
    IF(ALLOCATED(PROJECTED_DISTANCE)) DEALLOCATE(PROJECTED_DISTANCE)
    IF(ALLOCATED(PROJECTED_XI)) DEALLOCATE(PROJECTED_XI)
    IF(ALLOCATED(SORTING_IND_2)) DEALLOCATE(SORTING_IND_2)
    IF(ALLOCATED(GLOBAL_CLOSEST_DISTANCES)) DEALLOCATE(GLOBAL_CLOSEST_DISTANCES)
    IF(ALLOCATED(SORTING_IND_1)) DEALLOCATE(SORTING_IND_1)

    CALL EXITS("DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE")
    RETURN
999 IF(ALLOCATED(CANDIDATE_ELEMENTS)) DEALLOCATE(CANDIDATE_ELEMENTS)
    IF(ALLOCATED(CANDIDATE_FACES)) DEALLOCATE(CANDIDATE_FACES)
    IF(ALLOCATED(CLOSEST_ELEMENTS)) DEALLOCATE(CLOSEST_ELEMENTS)
    IF(ALLOCATED(CLOSEST_FACES)) DEALLOCATE(CLOSEST_FACES)
    IF(ALLOCATED(CLOSEST_DISTANCES)) DEALLOCATE(CLOSEST_DISTANCES)
    IF(ALLOCATED(GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES)) DEALLOCATE(GLOBAL_TO_LOCAL_NUMBER_OF_CLOSEST_CANDIDATES)
    IF(ALLOCATED(GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES)) DEALLOCATE(GLOBAL_NUMBER_OF_CLOSEST_CANDIDATES)
    IF(ALLOCATED(GLOBAL_MPI_DISPLACEMENTS)) DEALLOCATE(GLOBAL_MPI_DISPLACEMENTS)
    IF(ALLOCATED(GLOBAL_NUMBER_OF_PROJECTED_POINTS)) DEALLOCATE(GLOBAL_NUMBER_OF_PROJECTED_POINTS)
    IF(ALLOCATED(PROJECTION_EXIT_TAG)) DEALLOCATE(PROJECTION_EXIT_TAG)
    IF(ALLOCATED(PROJECTED_ELEMENT)) DEALLOCATE(PROJECTED_ELEMENT)
    IF(ALLOCATED(PROJECTED_FACE)) DEALLOCATE(PROJECTED_FACE)
    IF(ALLOCATED(PROJECTED_DISTANCE)) DEALLOCATE(PROJECTED_DISTANCE)
    IF(ALLOCATED(PROJECTED_XI)) DEALLOCATE(PROJECTED_XI)
    IF(ALLOCATED(SORTING_IND_2)) DEALLOCATE(SORTING_IND_2)
    IF(ALLOCATED(GLOBAL_CLOSEST_DISTANCES)) DEALLOCATE(GLOBAL_CLOSEST_DISTANCES)
    IF(ALLOCATED(SORTING_IND_1)) DEALLOCATE(SORTING_IND_1)
998 CALL ERRORS("DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_DATA_POINTS_PROJECTION_EVALUATE
  
  !
  !================================================================================================================================
  !
  
  !>Evaluate the data points position in a field based on data projection
  SUBROUTINE DataProjection_DataPointsPositionEvaluate(dataProjection,field,fieldVariableType,err,error,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<Data projection to give the xi locations and element number for the data points
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to be interpolated
    INTEGER(INTG), INTENT(IN) :: fieldVariableType !<The field variable type to be interpolated
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dataPointIdx,elementNumber,coordIdx
    TYPE(DATA_POINTS_TYPE), POINTER :: dataPoints
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoints(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: interpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    
    CALL ENTERS("DataProjection_DataPointsPositionEvaluate",err,error,*999)
    
    IF(ASSOCIATED(field)) THEN 
      IF (FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.FIELD%TYPE==FIELD_GEOMETRIC_GENERAL_TYPE) THEN
        IF(ASSOCIATED(dataProjection)) THEN
          dataPoints=>dataProjection%DATA_POINTS
          IF(ASSOCIATED(dataPoints)) THEN
            NULLIFY(interpolatedPoints)
            NULLIFY(interpolationParameters)
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(field,interpolationParameters,err,error,*999, &
              & FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoints,err,error,*999, &
              & FIELD_GEOMETRIC_COMPONENTS_TYPE)
            interpolatedPoint=>interpolatedPoints(fieldVariableType)%PTR
            !Loop through data points 
             DO dataPointIdx=1,dataPoints%NUMBER_OF_DATA_POINTS
               elementNumber=dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)%ELEMENT_NUMBER
               CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber, &
                 & interpolationParameters(fieldVariableType)%PTR,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
               CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,dataProjection%DATA_PROJECTION_RESULTS(dataPointIdx)%XI, &
                 & interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
               DO coordIdx=1,SIZE(dataPoints%DATA_POINTS(dataPointIdx)%position)
                 dataPoints%DATA_POINTS(dataPointIdx)%position(coordIdx)=interpolatedPoint%VALUES(coordIdx,NO_PART_DERIV)
               ENDDO !coordIdx     
             ENDDO !dataPointIdx
           ELSE
             CALL FLAG_ERROR("Data points is not associated.",err,error,*999)
           ENDIF
         ELSE
           CALL FLAG_ERROR("Data projection is not associated.",err,error,*999)
         ENDIF
       ELSE
         CALL FLAG_ERROR("Cannot evaluate data points position on field other than geometric or geometric general type.", &
           & err,error,*999)
       ENDIF
    ELSE
      CALL FLAG_ERROR("Field is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("DataProjection_DataPointsPositionEvaluate")
    RETURN
999 CALL ERRORS("DataProjection_DataPointsPositionEvaluate",err,error)    
    CALL EXITS("DataProjection_DataPointsPositionEvaluate")
    RETURN 1

  END SUBROUTINE DataProjection_DataPointsPositionEvaluate
  
  !
  !================================================================================================================================
  !
  
  !>Gets the maximum iteration update for a data projection.
  SUBROUTINE DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET(DATA_PROJECTION,MAXIMUM_ITERATION_UPDATE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the maximum iteration update for
    REAL(DP), INTENT(OUT) :: MAXIMUM_ITERATION_UPDATE !<On exit, the maximum iteration update of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    CALL ENTERS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        MAXIMUM_ITERATION_UPDATE=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE       
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the maximum iteration update for a data projection.
  SUBROUTINE DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET(DATA_PROJECTION,MAXIMUM_ITERATION_UPDATE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the maximum iteration update for
    REAL(DP), INTENT(IN) :: MAXIMUM_ITERATION_UPDATE !<the maximum iteration update to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE
        IF((MAXIMUM_ITERATION_UPDATE>=0.1).AND.(MAXIMUM_ITERATION_UPDATE<=1)) THEN
          DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE=MAXIMUM_ITERATION_UPDATE
        ELSE
          CALL FLAG_ERROR("Data projection maximum iteration update must be between 0.1 and 1.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_MAXIMUM_ITERATION_UPDATE_SET
  

  !
  !================================================================================================================================
  !
  
  !>Gets the maximum number of iterations for a data projection.
  SUBROUTINE DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET(DATA_PROJECTION,MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the maximum number of iterations for
    INTEGER(INTG), INTENT(OUT) :: MAXIMUM_NUMBER_OF_ITERATIONS !<On exit, the maximum number of iterations of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN      
        MAXIMUM_NUMBER_OF_ITERATIONS=DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS       
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the maximum number of iterations for a data projection.
  SUBROUTINE DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET(DATA_PROJECTION,MAXIMUM_NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the maximum number of iterations for
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<the maximum number of iterations to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE      
        IF(MAXIMUM_NUMBER_OF_ITERATIONS>=1) THEN
          DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_NUMBER_OF_ITERATIONS
        ELSE
          CALL FLAG_ERROR("Data projection maximum number of iterations must be at least 1.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_MAXIMUM_NUMBER_OF_ITERATIONS_SET
  
  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 1D elements
  SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_DISTANCE,PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(1)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    LOGICAL :: INSIDE_REGION,CONVERGED
    INTEGER(INTG) :: ELEMENT_NUMBER !local element number in current computational domain
    INTEGER(INTG) :: REGION_DIMENSIONS
    INTEGER(INTG) :: BOUND,EXIT_TAG
    REAL(DP) :: XI(1),XI_NEW(1),XI_UPDATE(1),XI_UPDATE_NORM !<xi
    REAL(DP) :: RELATIVE_TOLERANCE,ABSOLUTE_TOLERANCE !<tolerances
    REAL(DP) :: DISTANCE_VECTOR(3),FUNCTION_VALUE,FUNCTION_VALUE_NEW
    REAL(DP) :: FUNCTION_GRADIENT,FUNCTION_HESSIAN
    REAL(DP) :: MAXIMUM_DELTA,MINIMUM_DELTA,DELTA !<trust region size
    REAL(DP) :: PREDICTED_REDUCTION,PREDICTION_ACCURACY
    
    INTEGER(INTG) :: ne,itr1,itr2
    
    CALL ENTERS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1",ERR,ERROR,*999)
              
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        PROJECTION_EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          CONVERGED=.FALSE.
          DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,INTERPOLATED_POINT% &
            & INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          XI=DATA_PROJECTION%STARTING_XI
          CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
          FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))       
          main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
            !Check for bounds [0,1]
            IF(ABS(XI(1))<ZERO_TOLERANCE) THEN
              BOUND=-1 !bound at negative direction             
            ELSEIF(ABS(XI(1)-1.0_DP)<ZERO_TOLERANCE) THEN
              BOUND=1 !bound at positive direction
            ELSE !inside the bounds
              BOUND=0
            ENDIF              
            !FUNCTION_GRADIENT 
            FUNCTION_GRADIENT=-2.0_DP* &
              & (DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,FIRST_PART_DERIV)))
            !FUNCTION_HESSIAN 
            FUNCTION_HESSIAN=-2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,SECOND_PART_DERIV))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,FIRST_PART_DERIV),INTERPOLATED_POINT%VALUES(:,FIRST_PART_DERIV)))
            !A model trust region approach, directly taken from CMISS CLOS22: V = -(H + EIGEN_SHIFT*I)g
            !The calculation of EIGEN_SHIFT are only approximated as oppose to the common trust region approach               
            DO itr2=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(inner loop: adjust region size) usually EXIT at 1 or 2 iterations
              INSIDE_REGION=.FALSE.
              IF(FUNCTION_HESSIAN>ABSOLUTE_TOLERANCE) THEN !positive: minimum exists
                XI_UPDATE(1)=-FUNCTION_GRADIENT/FUNCTION_HESSIAN
                XI_UPDATE_NORM=DABS(XI_UPDATE(1))
                INSIDE_REGION=XI_UPDATE_NORM<=DELTA
              ENDIF !positive                 
              IF(.NOT.INSIDE_REGION) THEN !minimum not in the region
                XI_UPDATE(1)=-DSIGN(DELTA,FUNCTION_GRADIENT)
                XI_UPDATE_NORM=DELTA
              ENDIF
              IF((BOUND/=0).AND.(BOUND>0.EQV.XI_UPDATE(1)>0.0_DP)) THEN !projection go out of element bound
                EXIT_TAG=DATA_PROJECTION_EXIT_TAG_BOUNDS
                EXIT main_loop
              ENDIF
              CONVERGED=XI_UPDATE_NORM<ABSOLUTE_TOLERANCE !first half of the convergence test (before collision detection)
              XI_NEW=XI+XI_UPDATE !update XI
              IF(XI_NEW(1)<0.0_DP) THEN !boundary collision check
                XI_NEW(1)=0.0_DP
              ELSEIF(XI_NEW(1)>1.0_DP) THEN
                XI_NEW(1)=1.0_DP  
              ENDIF
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              DISTANCE_VECTOR=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
              CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test
              IF(CONVERGED) EXIT !converged: exit inner loop first
              IF((FUNCTION_VALUE_NEW-FUNCTION_VALUE)>ABSOLUTE_TOLERANCE) THEN !bad model: reduce step size
                IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                  EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                  EXIT main_loop
                ENDIF
                DELTA=DMAX1(MINIMUM_DELTA,0.25_DP*DELTA)
              ELSE
                PREDICTED_REDUCTION=XI_UPDATE(1)*(FUNCTION_GRADIENT+0.5_DP*FUNCTION_HESSIAN*XI_UPDATE(1))
                PREDICTION_ACCURACY=(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/PREDICTED_REDUCTION
                IF(PREDICTION_ACCURACY<0.01_DP) THEN !bad model: reduce region size
                  IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                    EXIT main_loop
                  ENDIF
                  DELTA=DMAX1(MINIMUM_DELTA,0.5_DP*DELTA)
                ELSEIF(PREDICTION_ACCURACY>0.9_DP.AND.PREDICTION_ACCURACY<1.1_DP) THEN !good model: increase region size
                  DELTA=DMIN1(MAXIMUM_DELTA,2.0_DP*DELTA)
                  EXIT
                ELSE !ok model: keep the current region size
                  EXIT
                ENDIF
              ENDIF
            ENDDO !itr2 (inner loop: adjust region size)
            FUNCTION_VALUE=FUNCTION_VALUE_NEW
            XI=XI_NEW
            IF(CONVERGED) THEN
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_CONVERGED
              EXIT
            ENDIF
          ENDDO main_loop !itr1 (outer loop)
          IF(EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.itr1>=DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS) &
            & EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
          IF((PROJECTION_EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(FUNCTION_VALUE)<PROJECTION_DISTANCE)) THEN
            PROJECTION_EXIT_TAG=EXIT_TAG
            PROJECTION_ELEMENT_NUMBER=ELEMENT_NUMBER
            PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
            PROJECTION_XI=XI
          ENDIF
        ENDDO !ne
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF    
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1
  
  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 2D elements
  SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_DISTANCE,PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(2)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    LOGICAL :: FREE,CONVERGED,INSIDE_REGION
    INTEGER(INTG) :: ELEMENT_NUMBER
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,REGION_DIMENSIONS
    INTEGER(INTG) :: BOUND(2),EXIT_TAG
    REAL(DP) :: XI(2),XI_NEW(2),XI_UPDATE(2),XI_UPDATE_NORM !<xi
    REAL(DP) :: RELATIVE_TOLERANCE,ABSOLUTE_TOLERANCE !<tolerances
    REAL(DP) :: DISTANCE_VECTOR(3),FUNCTION_VALUE,FUNCTION_VALUE_NEW
    REAL(DP) :: FUNCTION_GRADIENT(2),FUNCTION_GRADIENT_NORM
    REAL(DP) :: FUNCTION_HESSIAN(2,2),HESSIAN_DIAGONAL(2)
    REAL(DP) :: TEMP1,TEMP2,DET,EIGEN_MIN,EIGEN_MAX,EIGEN_SHIFT
    REAL(DP) :: MAXIMUM_DELTA,MINIMUM_DELTA,DELTA !<trust region size
    REAL(DP) :: PREDICTED_REDUCTION,PREDICTION_ACCURACY
    
    
    INTEGER(INTG) :: ne,ni,nifix,itr1,itr2
    
    CALL ENTERS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2",ERR,ERROR,*999)
              
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        PROJECTION_EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        MESH_COMPONENT_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN_MAPPING=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)% &
          & PTR%MAPPINGS%ELEMENTS
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          CONVERGED=.FALSE.
          DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet            
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          XI=DATA_PROJECTION%STARTING_XI
          CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
          FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))   
          main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
            !Check for bounds [0,1]
            DO ni=1,2 
              IF(ABS(XI(ni))<ZERO_TOLERANCE) THEN
                BOUND(ni)=-1 !bound at negative direction             
              ELSEIF(ABS(XI(ni)-1.0_DP)<ZERO_TOLERANCE) THEN
                BOUND(ni)=1 !bound at positive direction
              ELSE !inside the bounds
                BOUND(ni)=0
              ENDIF
            ENDDO !ni              
            !FUNCTION_GRADIENT 
            FUNCTION_GRADIENT(1)= &
              & -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1)))
            FUNCTION_GRADIENT(2)= &
              & -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            !FUNCTION_HESSIAN 
            FUNCTION_HESSIAN(1,1)= -2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S1))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1))) 
            FUNCTION_HESSIAN(1,2)= -2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S2))- &         
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            FUNCTION_HESSIAN(2,2)= -2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2_S2))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            !A model trust region approach, a Newton step is taken if the minimum lies inside the trust region (DELTA), if not, shift the step towards the steepest descent
            TEMP1=0.5_DP*(FUNCTION_HESSIAN(1,1)+FUNCTION_HESSIAN(2,2))
            TEMP2=DSQRT((0.5_DP*(FUNCTION_HESSIAN(1,1)-FUNCTION_HESSIAN(2,2)))**2+FUNCTION_HESSIAN(1,2)**2)
            EIGEN_MIN=TEMP1-TEMP2
            EIGEN_MAX=TEMP1+TEMP2
            FUNCTION_GRADIENT_NORM=DSQRT(DOT_PRODUCT(FUNCTION_GRADIENT,FUNCTION_GRADIENT))
            DO itr2=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(inner loop: adjust region size) usually EXIT at 1 or 2 iterations
              TEMP1=FUNCTION_GRADIENT_NORM/DELTA
              INSIDE_REGION=(EIGEN_MIN>=TEMP1).AND.(EIGEN_MIN>ABSOLUTE_TOLERANCE) !estimate if the solution is inside the trust region without calculating a newton step, this also guarantees the hessian matrix is positive definite
              IF(INSIDE_REGION) THEN
                DET=EIGEN_MIN*EIGEN_MAX !det(H)
                HESSIAN_DIAGONAL(1)=FUNCTION_HESSIAN(1,1)
                HESSIAN_DIAGONAL(2)=FUNCTION_HESSIAN(2,2)
              ELSE
                EIGEN_SHIFT=MAX(TEMP1-EIGEN_MIN,ABSOLUTE_TOLERANCE) !shift towards steepest decent
                DET=TEMP1*(EIGEN_MAX+EIGEN_SHIFT) !det(H)
                HESSIAN_DIAGONAL(1)=FUNCTION_HESSIAN(1,1)+EIGEN_SHIFT
                HESSIAN_DIAGONAL(2)=FUNCTION_HESSIAN(2,2)+EIGEN_SHIFT
              ENDIF
              XI_UPDATE(1)=-(HESSIAN_DIAGONAL(2)*FUNCTION_GRADIENT(1)-FUNCTION_HESSIAN(1,2)*FUNCTION_GRADIENT(2))/DET
              XI_UPDATE(2)=(FUNCTION_HESSIAN(1,2)*FUNCTION_GRADIENT(1)-HESSIAN_DIAGONAL(1)*FUNCTION_GRADIENT(2))/DET
              XI_UPDATE_NORM=DSQRT(DOT_PRODUCT(XI_UPDATE,XI_UPDATE))
              FREE=.TRUE.
              DO ni=1,2
                IF((BOUND(ni)/=0).AND.(BOUND(ni)>0.EQV.XI_UPDATE(ni)>0.0_DP)) THEN !projection go out of element bound
                  IF(.NOT.FREE) THEN !both xi are fixed
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_BOUNDS
                    EXIT main_loop
                  ENDIF
                  FREE=.FALSE.
                  nifix=ni
                ENDIF
              ENDDO !ni
              IF(FREE) THEN !both xi are free
                IF(.NOT.INSIDE_REGION) THEN
                  IF(XI_UPDATE_NORM>0.0_DP) THEN
                    XI_UPDATE=DELTA/XI_UPDATE_NORM*XI_UPDATE !readjust XI_UPDATE to lie on the region bound                      
                  ENDIF
                ENDIF
              ELSE !xi are not free
                XI_UPDATE(nifix)=0.0_DP
                ni=3-nifix
                INSIDE_REGION=.FALSE.
                IF(FUNCTION_HESSIAN(ni,ni)>0.0_DP) THEN !positive: minimum exists in the unbounded direction                
                  XI_UPDATE(ni)=-FUNCTION_GRADIENT(ni)/FUNCTION_HESSIAN(ni,ni)
                  XI_UPDATE_NORM=DABS(XI_UPDATE(ni))
                  INSIDE_REGION=XI_UPDATE_NORM<=DELTA
                ENDIF
                IF(.NOT.INSIDE_REGION) THEN !minimum not in the region
                  XI_UPDATE(ni)=-DSIGN(DELTA,FUNCTION_GRADIENT(ni))
                  XI_UPDATE_NORM=DELTA
                ENDIF            
              ENDIF !if xi is free
              CONVERGED=XI_UPDATE_NORM<ABSOLUTE_TOLERANCE !first half of the convergence test
              XI_NEW=XI+XI_UPDATE !update XI
              DO ni=1,2
                IF(XI_NEW(ni)<0.0_DP) THEN !boundary collision check
                  XI_NEW(ni)=0.0_DP
                  XI_NEW(3-ni)=XI(3-ni)-XI_UPDATE(3-ni)*XI(ni)/XI_UPDATE(ni)
                ELSEIF(XI_NEW(ni)>1.0_DP) THEN
                  XI_NEW(ni)=1.0_DP  
                  XI_NEW(3-ni)=XI(3-ni)+XI_UPDATE(3-ni)*(1.0_DP-XI(ni))/XI_UPDATE(ni)
                ENDIF
              ENDDO
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
              CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test (before collision detection)
              IF(CONVERGED) EXIT !converged: exit inner loop first
              IF((FUNCTION_VALUE_NEW-FUNCTION_VALUE)>ABSOLUTE_TOLERANCE) THEN !bad model: reduce step size
                IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                  EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                  EXIT main_loop
                ENDIF
                DELTA=DMAX1(MINIMUM_DELTA,0.25_DP*DELTA)
              ELSE
                PREDICTED_REDUCTION=DOT_PRODUCT(FUNCTION_GRADIENT,XI_UPDATE)+ &
                  & 0.5_DP*(XI_UPDATE(1)*(XI_UPDATE(1)*FUNCTION_HESSIAN(1,1)+2.0_DP*XI_UPDATE(2)*FUNCTION_HESSIAN(1,2))+ &
                  & XI_UPDATE(2)**2*FUNCTION_HESSIAN(2,2))
                PREDICTION_ACCURACY=(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/PREDICTED_REDUCTION
                IF(PREDICTION_ACCURACY<0.01_DP) THEN !bad model: reduce region size
                  IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                    EXIT main_loop
                  ENDIF
                  DELTA=DMAX1(MINIMUM_DELTA,0.5_DP*DELTA)
                ELSEIF(PREDICTION_ACCURACY>0.9_DP.AND.PREDICTION_ACCURACY<1.1_DP) THEN !good model: increase region size
                  DELTA=DMIN1(MAXIMUM_DELTA,2.0_DP*DELTA)
                  EXIT
                ELSE !ok model: keep the current region size
                  EXIT
                ENDIF
              ENDIF
            ENDDO !itr2 (inner loop: adjust region size)
            FUNCTION_VALUE=FUNCTION_VALUE_NEW
            XI=XI_NEW
            IF(CONVERGED) THEN
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_CONVERGED
              EXIT
            ENDIF
          ENDDO main_loop !itr1 (outer loop)
          IF(EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.itr1>=DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS) &
            & EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
          IF((PROJECTION_EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(FUNCTION_VALUE)<PROJECTION_DISTANCE)) THEN
            !IF(.NOT.ELEMENT_FOUND) ELEMENT_FOUND=.TRUE.
            PROJECTION_EXIT_TAG=EXIT_TAG
            PROJECTION_ELEMENT_NUMBER=ELEMENT_NUMBER
            PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
            PROJECTION_XI=XI
          ENDIF
        ENDDO !ne
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2

  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 3D elements
  SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_DISTANCE,PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(3)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    LOGICAL :: FREE,CONVERGED,INSIDE_REGION
    INTEGER(INTG) :: ELEMENT_NUMBER
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: NBOUND,BOUND(3),EXIT_TAG
    REAL(DP) :: XI(3),XI_NEW(3),XI_UPDATE(3),XI_UPDATE_NORM !<xi
    REAL(DP) :: RELATIVE_TOLERANCE,ABSOLUTE_TOLERANCE !<tolerances
    REAL(DP) :: DISTANCE_VECTOR(3),FUNCTION_VALUE,FUNCTION_VALUE_NEW
    REAL(DP) :: FUNCTION_GRADIENT(3),FUNCTION_GRADIENT_NORM,FUNCTION_GRADIENT2(2)
    REAL(DP) :: FUNCTION_HESSIAN(3,3),HESSIAN_DIAGONAL(3),FUNCTION_HESSIAN2(2,2),HESSIAN_DIAGONAL2(2)
    REAL(DP) :: TEMP1,TEMP2,TEMP3,TEMP4,DET,TRACE,TRACE2,EIGEN_MIN,EIGEN_MAX,EIGEN_SHIFT    
    REAL(DP) :: MAXIMUM_DELTA,MINIMUM_DELTA,DELTA !<trust region size
    REAL(DP) :: PREDICTED_REDUCTION,PREDICTION_ACCURACY
    
    
    INTEGER(INTG) :: ne,ni,ni2(2),nb,nifix,nifix2(2),itr1,itr2,nbfix
    
    CALL ENTERS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3",ERR,ERROR,*999)
              
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        PROJECTION_EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        MESH_COMPONENT_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN_MAPPING=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)% &
          & PTR%MAPPINGS%ELEMENTS
        !REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          CONVERGED=.FALSE.
          DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet            
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          XI=DATA_PROJECTION%STARTING_XI
          CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
          FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR,DISTANCE_VECTOR)   
          main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
            !Check for bounds [0,1]
            DO ni=1,3
              IF(ABS(XI(ni))<ZERO_TOLERANCE) THEN
                BOUND(ni)=-1 !bound at negative direction             
              ELSEIF(ABS(XI(ni)-1.0_DP)<ZERO_TOLERANCE) THEN
                BOUND(ni)=1 !bound at positive direction
              ELSE !inside the bounds
                BOUND(ni)=0
              ENDIF
            ENDDO !ni              
            !FUNCTION_GRADIENT 
            FUNCTION_GRADIENT(1)=-2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1)))
            FUNCTION_GRADIENT(2)=-2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            FUNCTION_GRADIENT(3)=-2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S3)))
            !FUNCTION_HESSIAN 
            FUNCTION_HESSIAN(1,1)= -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S1))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1))) 
            FUNCTION_HESSIAN(1,2)= -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S2))- &         
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            FUNCTION_HESSIAN(1,3)= -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S3))- &         
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S3)))
            FUNCTION_HESSIAN(2,2)= -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2_S2))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            FUNCTION_HESSIAN(2,3)= -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2_S3))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S3)))
            FUNCTION_HESSIAN(3,3)= -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR,INTERPOLATED_POINT%VALUES(:,PART_DERIV_S3_S3))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S3),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S3)))    
            !A model trust region approach, a Newton step is taken if the solution lies inside the trust region (DELTA), if not, shift the step towards the steepest descent
            TRACE=FUNCTION_HESSIAN(1,1)+FUNCTION_HESSIAN(2,2)+FUNCTION_HESSIAN(3,3) !tr(H)
            TRACE2=FUNCTION_HESSIAN(1,1)*FUNCTION_HESSIAN(2,2)+FUNCTION_HESSIAN(1,1)*FUNCTION_HESSIAN(3,3)+ &
              & FUNCTION_HESSIAN(2,2)*FUNCTION_HESSIAN(3,3)-FUNCTION_HESSIAN(1,2)**2-FUNCTION_HESSIAN(1,3)**2- &
              & FUNCTION_HESSIAN(2,3)**2 !tr(H**2)-(tr(H))**2
            DET=FUNCTION_HESSIAN(1,1)*FUNCTION_HESSIAN(2,2)*FUNCTION_HESSIAN(3,3)- &
              & FUNCTION_HESSIAN(1,1)*FUNCTION_HESSIAN(2,3)**2-FUNCTION_HESSIAN(2,2)*FUNCTION_HESSIAN(1,3)**2- &
              & FUNCTION_HESSIAN(3,3)*FUNCTION_HESSIAN(1,2)**2+ &
              & 2.0_DP*FUNCTION_HESSIAN(1,2)*FUNCTION_HESSIAN(1,3)*FUNCTION_HESSIAN(2,3) !det(H)                 
            TEMP1=-TRACE/3.0_DP
            TEMP2=TRACE2/3.0_DP
            TEMP3=TEMP2-TEMP1**2 !<=0
            IF(TEMP3>-1.0E-5_DP) THEN !include some negatives for numerical errors
              EIGEN_MIN=-TEMP1 !all eigenvalues are the same                
            ELSE
              TEMP3=DSQRT(-TEMP3)
              TEMP4=(DET+3.0_DP*(TEMP1*TEMP2)-2.0_DP*TEMP1**3)/(2.0_DP*TEMP3**3)
              EIGEN_MIN=2.0_DP*TEMP3*DCOS((DACOS(TEMP4)+TWOPI)/3.0_DP)-TEMP1                
            ENDIF
            FUNCTION_GRADIENT_NORM=DSQRT(DOT_PRODUCT(FUNCTION_GRADIENT,FUNCTION_GRADIENT))
            DO itr2=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(inner loop: adjust region size) usually EXIT at 1 or 2 iterations
              TEMP1=FUNCTION_GRADIENT_NORM/DELTA
              INSIDE_REGION=(EIGEN_MIN>=TEMP1).AND.(EIGEN_MIN>ABSOLUTE_TOLERANCE) !estimate if the solution is inside the trust region without calculating a newton step, this also guarantees the hessian matrix is positive definite
              IF(INSIDE_REGION) THEN
                HESSIAN_DIAGONAL(1)=FUNCTION_HESSIAN(1,1)
                HESSIAN_DIAGONAL(2)=FUNCTION_HESSIAN(2,2)
                HESSIAN_DIAGONAL(3)=FUNCTION_HESSIAN(3,3)
              ELSE
                EIGEN_SHIFT=MAX(TEMP1-EIGEN_MIN,ABSOLUTE_TOLERANCE) !shift towards steepest decent
                DET=DET+EIGEN_SHIFT*(TRACE2+EIGEN_SHIFT*(TRACE+EIGEN_SHIFT)) !shift the determinant
                HESSIAN_DIAGONAL(1)=FUNCTION_HESSIAN(1,1)+EIGEN_SHIFT
                HESSIAN_DIAGONAL(2)=FUNCTION_HESSIAN(2,2)+EIGEN_SHIFT
                HESSIAN_DIAGONAL(3)=FUNCTION_HESSIAN(3,3)+EIGEN_SHIFT
              ENDIF   
              TEMP2=FUNCTION_HESSIAN(1,3)*FUNCTION_HESSIAN(2,3)-FUNCTION_HESSIAN(1,2)*HESSIAN_DIAGONAL(3)
              TEMP3=FUNCTION_HESSIAN(1,2)*FUNCTION_HESSIAN(2,3)-FUNCTION_HESSIAN(1,3)*HESSIAN_DIAGONAL(2)
              TEMP4=FUNCTION_HESSIAN(1,2)*FUNCTION_HESSIAN(1,3)-FUNCTION_HESSIAN(2,3)*HESSIAN_DIAGONAL(1)               
              XI_UPDATE(1)=((FUNCTION_HESSIAN(2,3)**2-HESSIAN_DIAGONAL(2)*HESSIAN_DIAGONAL(3))*FUNCTION_GRADIENT(1)- &
                & TEMP2*FUNCTION_GRADIENT(2)-TEMP3*FUNCTION_GRADIENT(3))/DET
              XI_UPDATE(2)=((FUNCTION_HESSIAN(1,3)**2-HESSIAN_DIAGONAL(1)*HESSIAN_DIAGONAL(3))*FUNCTION_GRADIENT(2)- &                    
                & TEMP2*FUNCTION_GRADIENT(1)-TEMP4*FUNCTION_GRADIENT(3))/DET
              XI_UPDATE(3)=((FUNCTION_HESSIAN(1,2)**2-HESSIAN_DIAGONAL(1)*HESSIAN_DIAGONAL(2))*FUNCTION_GRADIENT(3)- &                    
                & TEMP3*FUNCTION_GRADIENT(1)-TEMP4*FUNCTION_GRADIENT(2))/DET
              XI_UPDATE_NORM=DSQRT(DOT_PRODUCT(XI_UPDATE,XI_UPDATE))
              FREE=.TRUE.
              NBOUND=0
              DO ni=1,3
                IF((BOUND(ni)/=0).AND.(BOUND(ni)>0.EQV.XI_UPDATE(ni)>0.0_DP)) THEN !projection go out of element bound
                  NBOUND=NBOUND+1
                  FREE=.FALSE.
                  IF(NBOUND<=2) THEN
                    nifix2(NBOUND)=ni
                  ELSE !all xi are fixed
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_BOUNDS
                    EXIT main_loop
                  ENDIF
                ENDIF
              ENDDO !ni
              IF(FREE) THEN !all xi are free
                IF(.NOT.INSIDE_REGION) THEN
                  IF(XI_UPDATE_NORM>0.0_DP) THEN
                    XI_UPDATE=DELTA/XI_UPDATE_NORM*XI_UPDATE !readjust XI_UPDATE to lie on the region bound                      
                  ENDIF
                ENDIF
              ELSE !at least one of the xi are not free
                !try 2D projection
                FREE=.TRUE.
                nifix=nifix2(1)
                IF(NBOUND==2) THEN
                  IF(XI_UPDATE(nifix2(2))>XI_UPDATE(nifix2(1))) nifix=nifix2(2) !only fix the direction that is most strongly suggesting leaving the element
                ENDIF
                XI_UPDATE(nifix)=0.0_DP
                ni2(1)=1+MOD(nifix,3)
                ni2(2)=1+MOD(nifix+1,3)
                !FUNCTION_GRADIENT2
                FUNCTION_GRADIENT2(1)=FUNCTION_GRADIENT(ni2(1))
                FUNCTION_GRADIENT2(2)=FUNCTION_GRADIENT(ni2(2))
                !FUNCTION_HESSIAN2
                FUNCTION_HESSIAN2(1,1)=FUNCTION_HESSIAN(ni2(1),ni2(1))
                FUNCTION_HESSIAN2(1,2)=FUNCTION_HESSIAN(ni2(1),ni2(2))
                FUNCTION_HESSIAN2(2,2)=FUNCTION_HESSIAN(ni2(2),ni2(2))
                !re-estimate the trust solution in 2D
                TEMP1=0.5_DP*(FUNCTION_HESSIAN2(1,1)+FUNCTION_HESSIAN2(2,2))
                TEMP2=DSQRT((0.5_DP*(FUNCTION_HESSIAN2(1,1)-FUNCTION_HESSIAN2(2,2)))**2+FUNCTION_HESSIAN2(1,2)**2)
                EIGEN_MIN=TEMP1-TEMP2
                EIGEN_MAX=TEMP1+TEMP2
                TEMP3=DSQRT(DOT_PRODUCT(FUNCTION_GRADIENT2,FUNCTION_GRADIENT2))/DELTA
                INSIDE_REGION=(EIGEN_MIN>=TEMP3).AND.(EIGEN_MIN>ABSOLUTE_TOLERANCE) !estimate if the solution is inside the trust region without calculating a newton step, this also guarantees the hessian matrix is positive definite
                IF(INSIDE_REGION) THEN
                  DET=EIGEN_MIN*EIGEN_MAX !determinant of FUNCTION_HESSIAN
                  HESSIAN_DIAGONAL2(1)=FUNCTION_HESSIAN2(1,1)
                  HESSIAN_DIAGONAL2(2)=FUNCTION_HESSIAN2(2,2)
                ELSE
                  EIGEN_SHIFT=MAX(TEMP3-EIGEN_MIN,ABSOLUTE_TOLERANCE) !shift towards steepest decent
                  DET=TEMP3*(EIGEN_MAX+EIGEN_SHIFT) !determinant of shifted FUNCTION_HESSIAN
                  HESSIAN_DIAGONAL2(1)=FUNCTION_HESSIAN2(1,1)+EIGEN_SHIFT
                  HESSIAN_DIAGONAL2(2)=FUNCTION_HESSIAN2(2,2)+EIGEN_SHIFT
                ENDIF
                XI_UPDATE(ni2(1))=-(HESSIAN_DIAGONAL2(2)*FUNCTION_GRADIENT2(1)-FUNCTION_HESSIAN2(1,2)*FUNCTION_GRADIENT2(2))/DET
                XI_UPDATE(ni2(2))=(FUNCTION_HESSIAN2(1,2)*FUNCTION_GRADIENT2(1)-HESSIAN_DIAGONAL2(1)*FUNCTION_GRADIENT2(2))/DET
                XI_UPDATE_NORM=DSQRT(DOT_PRODUCT(XI_UPDATE,XI_UPDATE))
                !check again for bounds
                DO nb=1,2
                  IF((BOUND(ni2(nb))/=0).AND.(BOUND(ni2(nb))>0.EQV.XI_UPDATE(ni2(nb))>0.0_DP)) THEN !projection go out of element bound
                    IF(.NOT.FREE) THEN !both xi are fixed
                      EXIT_TAG=DATA_PROJECTION_EXIT_TAG_BOUNDS
                      EXIT main_loop
                    ENDIF
                    FREE=.FALSE.
                    nbfix=nb
                  ENDIF
                ENDDO !ni
                IF(FREE) THEN !both xi are free
                  IF(.NOT.INSIDE_REGION) THEN
                    IF(XI_UPDATE_NORM>0.0_DP) THEN
                      XI_UPDATE=DELTA/XI_UPDATE_NORM*XI_UPDATE !readjust XI_UPDATE to lie on the region bound                      
                    ENDIF
                  ENDIF
                ELSE !xi are not free
                  XI_UPDATE(ni2(nbfix))=0.0_DP
                  ni=ni2(3-nbfix)
                  INSIDE_REGION=.FALSE.
                  IF(FUNCTION_HESSIAN(ni,ni)>0.0_DP) THEN !positive: minimum exists in the unbounded direction                
                    XI_UPDATE(ni)=-FUNCTION_GRADIENT(ni)/FUNCTION_HESSIAN(ni,ni)
                    XI_UPDATE_NORM=DABS(XI_UPDATE(ni))
                    INSIDE_REGION=XI_UPDATE_NORM<=DELTA
                  ENDIF
                  IF(.NOT.INSIDE_REGION) THEN !minimum not in the region
                    XI_UPDATE(ni)=-DSIGN(DELTA,FUNCTION_GRADIENT(ni))
                    XI_UPDATE_NORM=DELTA
                  ENDIF            
                ENDIF !if xi are free (2D)
              ENDIF !if xi are free (3D)
              CONVERGED=XI_UPDATE_NORM<ABSOLUTE_TOLERANCE !first half of the convergence test
              XI_NEW=XI+XI_UPDATE !update XI
              DO ni=1,3
                IF(XI_NEW(ni)<0.0_DP) THEN !boundary collision check
                  XI_NEW(ni)=0.0_DP
                  XI_NEW=XI-XI_UPDATE*XI(ni)/XI_UPDATE(ni)
                ELSEIF(XI_NEW(ni)>1.0_DP) THEN
                  XI_NEW(ni)=1.0_DP  
                  XI_NEW=XI+XI_UPDATE*(1.0_DP-XI(ni))/XI_UPDATE(ni)
                ENDIF
              ENDDO !ni
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              DISTANCE_VECTOR=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR,DISTANCE_VECTOR)
              CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test (before collision detection)
              IF(CONVERGED) EXIT !converged: exit inner loop first
              IF((FUNCTION_VALUE_NEW-FUNCTION_VALUE)>ABSOLUTE_TOLERANCE) THEN !bad model: reduce step size
                IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                  EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                  EXIT main_loop
                ENDIF
                DELTA=DMAX1(MINIMUM_DELTA,0.25_DP*DELTA)
              ELSE
                PREDICTED_REDUCTION=DOT_PRODUCT(FUNCTION_GRADIENT,XI_UPDATE)+ &
                  & 0.5_DP*(XI_UPDATE(1)*(XI_UPDATE(1)*FUNCTION_HESSIAN(1,1)+2.0_DP*XI_UPDATE(2)*FUNCTION_HESSIAN(1,2)+ &
                  & 2.0_DP*XI_UPDATE(3)*FUNCTION_HESSIAN(1,3))+XI_UPDATE(2)*(XI_UPDATE(2)*FUNCTION_HESSIAN(2,2)+ &
                  & 2.0_DP*XI_UPDATE(3)*FUNCTION_HESSIAN(2,3))+XI_UPDATE(2)**2*FUNCTION_HESSIAN(2,2))
                PREDICTION_ACCURACY=(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/PREDICTED_REDUCTION
                IF(PREDICTION_ACCURACY<0.01_DP) THEN !bad model: reduce region size
                  IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                    EXIT main_loop
                  ENDIF
                  DELTA=DMAX1(MINIMUM_DELTA,0.5_DP*DELTA)
                ELSEIF(PREDICTION_ACCURACY>0.9_DP.AND.PREDICTION_ACCURACY<1.1_DP) THEN !good model: increase region size
                  DELTA=DMIN1(MAXIMUM_DELTA,2.0_DP*DELTA)
                  EXIT
                ELSE !ok model: keep the current region size
                  EXIT
                ENDIF
              ENDIF
            ENDDO !itr2 (inner loop: adjust region size)
            FUNCTION_VALUE=FUNCTION_VALUE_NEW
            XI=XI_NEW
            IF(CONVERGED) THEN
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_CONVERGED
              EXIT
            ENDIF
          ENDDO main_loop !itr1 (outer loop)
          IF(EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.itr1>=DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS) &
            & EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
          IF((PROJECTION_EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(FUNCTION_VALUE)<PROJECTION_DISTANCE)) THEN
            !IF(.NOT.ELEMENT_FOUND) ELEMENT_FOUND=.TRUE.
            PROJECTION_EXIT_TAG=EXIT_TAG
            PROJECTION_ELEMENT_NUMBER=ELEMENT_NUMBER
            PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
            PROJECTION_XI=XI
          ENDIF
        ENDDO !ne
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3
  
  !
  !================================================================================================================================
  !

  !>Find the projection of a data point onto element faces (slight difference to DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2)
  SUBROUTINE DATA_PROJECTION_NEWTON_FACES_EVALUATE(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & CANDIDATE_ELEMENT_FACES,PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_ELEMENT_FACE_NUMBER,PROJECTION_DISTANCE, &
    & PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENT_FACES(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_FACE_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(2)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    LOGICAL :: FREE,CONVERGED,INSIDE_REGION
    INTEGER(INTG) :: ELEMENT_NUMBER,ELEMENT_FACE_NUMBER,FACE_NUMBER
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,REGION_DIMENSIONS
    INTEGER(INTG) :: BOUND(2),EXIT_TAG
    REAL(DP) :: XI(2),XI_NEW(2),XI_UPDATE(2),XI_UPDATE_NORM !<xi
    REAL(DP) :: RELATIVE_TOLERANCE,ABSOLUTE_TOLERANCE !<tolerances
    REAL(DP) :: DISTANCE_VECTOR(3),FUNCTION_VALUE,FUNCTION_VALUE_NEW
    REAL(DP) :: FUNCTION_GRADIENT(2),FUNCTION_GRADIENT_NORM
    REAL(DP) :: FUNCTION_HESSIAN(2,2),HESSIAN_DIAGONAL(2)
    REAL(DP) :: TEMP1,TEMP2,DET,EIGEN_MIN,EIGEN_MAX,EIGEN_SHIFT
    REAL(DP) :: MAXIMUM_DELTA,MINIMUM_DELTA,DELTA !<trust region size
    REAL(DP) :: PREDICTED_REDUCTION,PREDICTION_ACCURACY
    
    
    INTEGER(INTG) :: ne,ni,nifix,itr1,itr2
    
    CALL ENTERS("DATA_PROJECTION_NEWTON_FACES_EVALUATE",ERR,ERROR,*999)
              
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        PROJECTION_EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        MESH_COMPONENT_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN_MAPPING=>INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)% &
          & PTR%MAPPINGS%ELEMENTS
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          ELEMENT_FACE_NUMBER=CANDIDATE_ELEMENT_FACES(ne)
          FACE_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS( &
            & ELEMENT_NUMBER)%ELEMENT_FACES(ELEMENT_FACE_NUMBER)     
          EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          CONVERGED=.FALSE.
          DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet            
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,FACE_NUMBER, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          XI=DATA_PROJECTION%STARTING_XI
          CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
          FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))   
          main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
            !Check for bounds [0,1]
            DO ni=1,2 
              IF(ABS(XI(ni))<ZERO_TOLERANCE) THEN
                BOUND(ni)=-1 !bound at negative direction             
              ELSEIF(ABS(XI(ni)-1.0_DP)<ZERO_TOLERANCE) THEN
                BOUND(ni)=1 !bound at positive direction
              ELSE !inside the bounds
                BOUND(ni)=0
              ENDIF
            ENDDO !ni              
            !FUNCTION_GRADIENT 
            FUNCTION_GRADIENT(1)= &
              & -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1)))
            FUNCTION_GRADIENT(2)= &
              & -2.0_DP*(DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            !FUNCTION_HESSIAN 
            FUNCTION_HESSIAN(1,1)= -2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S1))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1))) 
            FUNCTION_HESSIAN(1,2)= -2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1_S2))- &         
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S1),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            FUNCTION_HESSIAN(2,2)= -2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2_S2))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2),INTERPOLATED_POINT%VALUES(:,PART_DERIV_S2)))
            !A model trust region approach, a Newton step is taken if the minimum lies inside the trust region (DELTA), if not, shift the step towards the steepest descent
            TEMP1=0.5_DP*(FUNCTION_HESSIAN(1,1)+FUNCTION_HESSIAN(2,2))
            TEMP2=DSQRT((0.5_DP*(FUNCTION_HESSIAN(1,1)-FUNCTION_HESSIAN(2,2)))**2+FUNCTION_HESSIAN(1,2)**2)
            EIGEN_MIN=TEMP1-TEMP2
            EIGEN_MAX=TEMP1+TEMP2
            FUNCTION_GRADIENT_NORM=DSQRT(DOT_PRODUCT(FUNCTION_GRADIENT,FUNCTION_GRADIENT))
            DO itr2=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(inner loop: adjust region size) usually EXIT at 1 or 2 iterations
              TEMP1=FUNCTION_GRADIENT_NORM/DELTA
              INSIDE_REGION=(EIGEN_MIN>=TEMP1).AND.(EIGEN_MIN>ABSOLUTE_TOLERANCE) !estimate if the solution is inside the trust region without calculating a newton step, this also guarantees the hessian matrix is positive definite
              IF(INSIDE_REGION) THEN
                DET=EIGEN_MIN*EIGEN_MAX !det(H)
                HESSIAN_DIAGONAL(1)=FUNCTION_HESSIAN(1,1)
                HESSIAN_DIAGONAL(2)=FUNCTION_HESSIAN(2,2)
              ELSE
                EIGEN_SHIFT=MAX(TEMP1-EIGEN_MIN,ABSOLUTE_TOLERANCE) !shift towards steepest decent
                DET=TEMP1*(EIGEN_MAX+EIGEN_SHIFT) !det(H)
                HESSIAN_DIAGONAL(1)=FUNCTION_HESSIAN(1,1)+EIGEN_SHIFT
                HESSIAN_DIAGONAL(2)=FUNCTION_HESSIAN(2,2)+EIGEN_SHIFT
              ENDIF
              XI_UPDATE(1)=-(HESSIAN_DIAGONAL(2)*FUNCTION_GRADIENT(1)-FUNCTION_HESSIAN(1,2)*FUNCTION_GRADIENT(2))/DET
              XI_UPDATE(2)=(FUNCTION_HESSIAN(1,2)*FUNCTION_GRADIENT(1)-HESSIAN_DIAGONAL(1)*FUNCTION_GRADIENT(2))/DET
              XI_UPDATE_NORM=DSQRT(DOT_PRODUCT(XI_UPDATE,XI_UPDATE))
              FREE=.TRUE.
              DO ni=1,2
                IF((BOUND(ni)/=0).AND.(BOUND(ni)>0.EQV.XI_UPDATE(ni)>0.0_DP)) THEN !projection go out of element bound
                  IF(.NOT.FREE) THEN !both xi are fixed
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_BOUNDS
                    EXIT main_loop
                  ENDIF
                  FREE=.FALSE.
                  nifix=ni
                ENDIF
              ENDDO !ni
              IF(FREE) THEN !both xi are free
                IF(.NOT.INSIDE_REGION) THEN
                  IF(XI_UPDATE_NORM>0.0_DP) THEN
                    XI_UPDATE=DELTA/XI_UPDATE_NORM*XI_UPDATE !readjust XI_UPDATE to lie on the region bound                      
                  ENDIF
                ENDIF
              ELSE !xi are not free
                XI_UPDATE(nifix)=0.0_DP
                ni=3-nifix
                INSIDE_REGION=.FALSE.
                IF(FUNCTION_HESSIAN(ni,ni)>0.0_DP) THEN !positive: minimum exists in the unbounded direction                
                  XI_UPDATE(ni)=-FUNCTION_GRADIENT(ni)/FUNCTION_HESSIAN(ni,ni)
                  XI_UPDATE_NORM=DABS(XI_UPDATE(ni))
                  INSIDE_REGION=XI_UPDATE_NORM<=DELTA
                ENDIF
                IF(.NOT.INSIDE_REGION) THEN !minimum not in the region
                  XI_UPDATE(ni)=-DSIGN(DELTA,FUNCTION_GRADIENT(ni))
                  XI_UPDATE_NORM=DELTA
                ENDIF            
              ENDIF !if xi is free
              CONVERGED=XI_UPDATE_NORM<ABSOLUTE_TOLERANCE !first half of the convergence test
              XI_NEW=XI+XI_UPDATE !update XI
              DO ni=1,2
                IF(XI_NEW(ni)<0.0_DP) THEN !boundary collision check
                  XI_NEW(ni)=0.0_DP
                  XI_NEW(3-ni)=XI(3-ni)-XI_UPDATE(3-ni)*XI(ni)/XI_UPDATE(ni)
                ELSEIF(XI_NEW(ni)>1.0_DP) THEN
                  XI_NEW(ni)=1.0_DP  
                  XI_NEW(3-ni)=XI(3-ni)+XI_UPDATE(3-ni)*(1.0_DP-XI(ni))/XI_UPDATE(ni)
                ENDIF
              ENDDO
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              DISTANCE_VECTOR=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
              CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test (before collision detection)
              IF(CONVERGED) EXIT !converged: exit inner loop first
              IF((FUNCTION_VALUE_NEW-FUNCTION_VALUE)>ABSOLUTE_TOLERANCE) THEN !bad model: reduce step size
                IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                  EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                  EXIT main_loop
                ENDIF
                DELTA=DMAX1(MINIMUM_DELTA,0.25_DP*DELTA)
              ELSE
                PREDICTED_REDUCTION=DOT_PRODUCT(FUNCTION_GRADIENT,XI_UPDATE)+ &
                  & 0.5_DP*(XI_UPDATE(1)*(XI_UPDATE(1)*FUNCTION_HESSIAN(1,1)+2.0_DP*XI_UPDATE(2)*FUNCTION_HESSIAN(1,2))+ &
                  & XI_UPDATE(2)**2*FUNCTION_HESSIAN(2,2))
                PREDICTION_ACCURACY=(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/PREDICTED_REDUCTION
                IF(PREDICTION_ACCURACY<0.01_DP) THEN !bad model: reduce region size
                  IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                    EXIT main_loop
                  ENDIF
                  DELTA=DMAX1(MINIMUM_DELTA,0.5_DP*DELTA)
                ELSEIF(PREDICTION_ACCURACY>0.9_DP.AND.PREDICTION_ACCURACY<1.1_DP) THEN !good model: increase region size
                  DELTA=DMIN1(MAXIMUM_DELTA,2.0_DP*DELTA)
                  EXIT
                ELSE !ok model: keep the current region size
                  EXIT
                ENDIF
              ENDIF
            ENDDO !itr2 (inner loop: adjust region size)
            FUNCTION_VALUE=FUNCTION_VALUE_NEW
            XI=XI_NEW
            IF(CONVERGED) THEN
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_CONVERGED
              EXIT
            ENDIF
          ENDDO main_loop !itr1 (outer loop)
          IF(EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.itr1>=DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS) &
            & EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
          IF((PROJECTION_EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(FUNCTION_VALUE)<PROJECTION_DISTANCE)) THEN
            !IF(.NOT.ELEMENT_FOUND) ELEMENT_FOUND=.TRUE.
            PROJECTION_EXIT_TAG=EXIT_TAG
            PROJECTION_ELEMENT_NUMBER=ELEMENT_NUMBER
            PROJECTION_ELEMENT_FACE_NUMBER=ELEMENT_FACE_NUMBER            
            PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
            PROJECTION_XI=XI
          ENDIF
        ENDDO !ne
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NEWTON_FACES_EVALUATE")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NEWTON_FACES_EVALUATE",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NEWTON_FACES_EVALUATE")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NEWTON_FACES_EVALUATE

  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto element lines (slight difference to DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1)
  SUBROUTINE DATA_PROJECTION_NEWTON_LINES_EVALUATE(DATA_PROJECTION,INTERPOLATED_POINT,POINT_VALUES,CANDIDATE_ELEMENTS, &
    & CANDIDATE_ELEMENT_LINES,PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_ELEMENT_LINE_NUMBER,PROJECTION_DISTANCE, &
    & PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    REAL(DP), INTENT(IN) :: POINT_VALUES(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENT_LINES(:)  
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_LINE_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(1)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    LOGICAL :: INSIDE_REGION,CONVERGED
    INTEGER(INTG) :: ELEMENT_NUMBER,ELEMENT_LINE_NUMBER,LINE_NUMBER !local element number in current computational domain
    INTEGER(INTG) :: REGION_DIMENSIONS
    INTEGER(INTG) :: BOUND,EXIT_TAG
    REAL(DP) :: XI(1),XI_NEW(1),XI_UPDATE(1),XI_UPDATE_NORM !<xi
    REAL(DP) :: RELATIVE_TOLERANCE,ABSOLUTE_TOLERANCE !<tolerances
    REAL(DP) :: DISTANCE_VECTOR(3),FUNCTION_VALUE,FUNCTION_VALUE_NEW
    REAL(DP) :: FUNCTION_GRADIENT,FUNCTION_HESSIAN
    REAL(DP) :: MAXIMUM_DELTA,MINIMUM_DELTA,DELTA !<trust region size
    REAL(DP) :: PREDICTED_REDUCTION,PREDICTION_ACCURACY
    
    INTEGER(INTG) :: ne,itr1,itr2
    
    CALL ENTERS("DATA_PROJECTION_NEWTON_LINES_EVALUATE",ERR,ERROR,*999)
              
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        PROJECTION_EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          ELEMENT_LINE_NUMBER=CANDIDATE_ELEMENT_LINES(ne)
          LINE_NUMBER=INTERPOLATED_POINT%INTERPOLATION_PARAMETERS%FIELD%DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS( &
            & ELEMENT_NUMBER)%ELEMENT_LINES(ELEMENT_LINE_NUMBER)
          EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
          CONVERGED=.FALSE.
          DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet
          CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,LINE_NUMBER,INTERPOLATED_POINT% &
            & INTERPOLATION_PARAMETERS,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          XI=DATA_PROJECTION%STARTING_XI
          CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
          FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))       
          main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
            !Check for bounds [0,1]
            IF(ABS(XI(1))<ZERO_TOLERANCE) THEN
              BOUND=-1 !bound at negative direction             
            ELSEIF(ABS(XI(1)-1.0_DP)<ZERO_TOLERANCE) THEN
              BOUND=1 !bound at positive direction
            ELSE !inside the bounds
              BOUND=0
            ENDIF              
            !FUNCTION_GRADIENT 
            FUNCTION_GRADIENT=-2.0_DP* &
              & (DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,FIRST_PART_DERIV)))
            !FUNCTION_HESSIAN 
            FUNCTION_HESSIAN=-2.0_DP*(&
              & DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),INTERPOLATED_POINT%VALUES(:,SECOND_PART_DERIV))- &
              & DOT_PRODUCT(INTERPOLATED_POINT%VALUES(:,FIRST_PART_DERIV),INTERPOLATED_POINT%VALUES(:,FIRST_PART_DERIV)))
            !A model trust region approach, directly taken from CMISS CLOS22: V = -(H + EIGEN_SHIFT*I)g
            !The calculation of EIGEN_SHIFT are only approximated as oppose to the common trust region approach               
            DO itr2=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(inner loop: adjust region size) usually EXIT at 1 or 2 iterations
              INSIDE_REGION=.FALSE.
              IF(FUNCTION_HESSIAN>ABSOLUTE_TOLERANCE) THEN !positive: minimum exists
                XI_UPDATE(1)=-FUNCTION_GRADIENT/FUNCTION_HESSIAN
                XI_UPDATE_NORM=DABS(XI_UPDATE(1))
                INSIDE_REGION=XI_UPDATE_NORM<=DELTA
              ENDIF !positive                 
              IF(.NOT.INSIDE_REGION) THEN !minimum not in the region
                XI_UPDATE(1)=-DSIGN(DELTA,FUNCTION_GRADIENT)
                XI_UPDATE_NORM=DELTA
              ENDIF
              IF((BOUND/=0).AND.(BOUND>0.EQV.XI_UPDATE(1)>0.0_DP)) THEN !projection go out of element bound
                EXIT_TAG=DATA_PROJECTION_EXIT_TAG_BOUNDS
                EXIT main_loop
              ENDIF
              CONVERGED=XI_UPDATE_NORM<ABSOLUTE_TOLERANCE !first half of the convergence test (before collision detection)
              XI_NEW=XI+XI_UPDATE !update XI
              IF(XI_NEW(1)<0.0_DP) THEN !boundary collision check
                XI_NEW(1)=0.0_DP
              ELSEIF(XI_NEW(1)>1.0_DP) THEN
                XI_NEW(1)=1.0_DP  
              ENDIF
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              DISTANCE_VECTOR=POINT_VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
              CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test
              IF(CONVERGED) EXIT !converged: exit inner loop first
              IF((FUNCTION_VALUE_NEW-FUNCTION_VALUE)>ABSOLUTE_TOLERANCE) THEN !bad model: reduce step size
                IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                  EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                  EXIT main_loop
                ENDIF
                DELTA=DMAX1(MINIMUM_DELTA,0.25_DP*DELTA)
              ELSE
                PREDICTED_REDUCTION=XI_UPDATE(1)*(FUNCTION_GRADIENT+0.5_DP*FUNCTION_HESSIAN*XI_UPDATE(1))
                PREDICTION_ACCURACY=(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/PREDICTED_REDUCTION
                IF(PREDICTION_ACCURACY<0.01_DP) THEN !bad model: reduce region size
                  IF(DELTA<=MINIMUM_DELTA) THEN !something went wrong, MINIMUM_DELTA too large? not likely to happen if MINIMUM_DELTA is small
                    EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION ! it will get stucked!!
                    EXIT main_loop
                  ENDIF
                  DELTA=DMAX1(MINIMUM_DELTA,0.5_DP*DELTA)
                ELSEIF(PREDICTION_ACCURACY>0.9_DP.AND.PREDICTION_ACCURACY<1.1_DP) THEN !good model: increase region size
                  DELTA=DMIN1(MAXIMUM_DELTA,2.0_DP*DELTA)
                  EXIT
                ELSE !ok model: keep the current region size
                  EXIT
                ENDIF
              ENDIF
            ENDDO !itr2 (inner loop: adjust region size)
            FUNCTION_VALUE=FUNCTION_VALUE_NEW
            XI=XI_NEW
            IF(CONVERGED) THEN
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_CONVERGED
              EXIT
            ENDIF
          ENDDO main_loop !itr1 (outer loop)
          IF(EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT.AND.itr1>=DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS) &
            & EXIT_TAG=DATA_PROJECTION_EXIT_TAG_MAX_ITERATION
          IF((PROJECTION_EXIT_TAG==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT).OR.(DSQRT(FUNCTION_VALUE)<PROJECTION_DISTANCE)) THEN
            PROJECTION_EXIT_TAG=EXIT_TAG
            PROJECTION_ELEMENT_NUMBER=ELEMENT_NUMBER
            PROJECTION_ELEMENT_LINE_NUMBER=ELEMENT_LINE_NUMBER
            PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
            PROJECTION_XI=XI
          ENDIF
        ENDDO !ne
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF    
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NEWTON_LINES_EVALUATE")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NEWTON_LINES_EVALUATE",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NEWTON_LINES_EVALUATE")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NEWTON_LINES_EVALUATE
  
  !
  !================================================================================================================================
  !
  
  !>Gets the number of closest elements for a data projection.
  SUBROUTINE DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET(DATA_PROJECTION,NUMBER_OF_CLOSEST_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the number of closest elements for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_CLOSEST_ELEMENTS !<On exit, the number of closest elements of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        NUMBER_OF_CLOSEST_ELEMENTS=DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS       
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the number of closest elements for a data projection.
  SUBROUTINE DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET(DATA_PROJECTION,NUMBER_OF_CLOSEST_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the number of closest elements for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_CLOSEST_ELEMENTS !<the number of closest elements to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE
        IF(NUMBER_OF_CLOSEST_ELEMENTS>=1) THEN
          DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS=NUMBER_OF_CLOSEST_ELEMENTS
        ELSE
          CALL FLAG_ERROR("Data projection number of closest elements must be at least 1.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_NUMBER_OF_CLOSEST_ELEMENTS_SET
  
  !
  !================================================================================================================================
  !
  
  !>Sets the candidates element numbers and local line/face numbers for a data projection.
  SUBROUTINE DataProjection_ProjectionCandidatesSet(dataProjection,elementUserNumber,localFaceLineNumbers,err,error,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: elementUserNumber(:) !<the projection candidate user element numbers
    INTEGER(INTG), INTENT(IN) :: localFaceLineNumbers(:) !<the projection candidate element face/line numbers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,elementGlobalNumber
    INTEGER(INTG) :: meshComponentNumber=1 !<TODO:mesh component is harded coded to be 1, need to be removed once MeshComponentsElementsType is moved under MeshTopologyType
    LOGICAL :: elementExists
    
    CALL ENTERS("DataProjection_ProjectionCandidatesSet",err,error,*999)

    IF(ASSOCIATED(dataProjection)) THEN
      IF(SIZE(elementUserNumber,1)==SIZE(localFaceLineNumbers,1)) THEN
        ALLOCATE(dataProjection%candidateElementNumbers(SIZE(elementUserNumber,1)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidiate element numbers.",ERR,ERROR,*998)
        ALLOCATE(dataProjection%localFaceLineNumbers(SIZE(localFaceLineNumbers,1)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidiate local face/line numbers.",ERR,ERROR,*999)
        DO elementIdx=1,SIZE(elementUserNumber,1)
          CALL MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(dataProjection%MESH,meshComponentNumber,elementUserNumber(elementIdx), &
            & elementExists,elementGlobalNumber,err,error,*999)       
          IF(elementExists) THEN
            dataProjection%candidateElementNumbers(elementIdx)=elementUserNumber(elementIdx)
            dataProjection%localFaceLineNumbers(elementIdx)=localFaceLineNumbers(elementIdx)
          ELSE
            CALL FLAG_ERROR("Element with user number ("//TRIM(NUMBER_TO_VSTRING &
              & (elementUserNumber(elementIdx),"*",err,ERROR))//") does not exist.",err,error,*999)
          ENDIF
        ENDDO !elementIdx
      ELSE
        CALL FLAG_ERROR("Input user element numbers and face numbers sizes do not match.",err,error,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("DataProjection_ProjectionCandidatesSet")
    RETURN
999 IF(ALLOCATED(dataProjection%candidateElementNumbers)) THEN
      DEALLOCATE(dataProjection%candidateElementNumbers)
    END IF
    IF(ALLOCATED(dataProjection%localFaceLineNumbers)) THEN
      DEALLOCATE(dataProjection%localFaceLineNumbers)
    END IF
998 CALL ERRORS("DataProjection_ProjectionCandidatesSet",err,error)
    CALL EXITS("DataProjection_ProjectionCandidatesSet")
    RETURN 1

  END SUBROUTINE DataProjection_ProjectionCandidatesSet
  
  !
  !================================================================================================================================
  !
  
  !>Gets the projection type for a data projection.
  SUBROUTINE DATA_PROJECTION_PROJECTION_TYPE_GET(DATA_PROJECTION,PROJECTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the projection type for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_TYPE !<On exit, the projection type of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_PROJECTION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        PROJECTION_TYPE=DATA_PROJECTION%PROJECTION_TYPE       
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_PROJECTION_TYPE_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_PROJECTION_TYPE_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_PROJECTION_TYPE_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_PROJECTION_TYPE_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the projection type for a data projection.
  SUBROUTINE DATA_PROJECTION_PROJECTION_TYPE_SET(DATA_PROJECTION,PROJECTION_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the projection type for
    INTEGER(INTG), INTENT(IN) :: PROJECTION_TYPE !<the projection type to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP), ALLOCATABLE :: STARTING_XI(:)
    
    CALL ENTERS("DATA_PROJECTION_PROJECTION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE
        DATA_PROJECTION%PROJECTION_TYPE=PROJECTION_TYPE
        SELECT CASE(PROJECTION_TYPE)
          CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
            DATA_PROJECTION%NUMBER_OF_XI=1
          CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE)
            DATA_PROJECTION%NUMBER_OF_XI=2
          CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
            DATA_PROJECTION%NUMBER_OF_XI=DATA_PROJECTION%MESH%NUMBER_OF_DIMENSIONS
          CASE DEFAULT
            CALL FLAG_ERROR("Input projection type is undefined.",ERR,ERROR,*999)
        END SELECT
        IF(DATA_PROJECTION%NUMBER_OF_XI/=SIZE(DATA_PROJECTION%STARTING_XI,1)) THEN
          ALLOCATE(STARTING_XI(DATA_PROJECTION%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate starting xi.",ERR,ERROR,*999)
          IF(DATA_PROJECTION%NUMBER_OF_XI>SIZE(DATA_PROJECTION%STARTING_XI,1)) THEN
            STARTING_XI(1:SIZE(DATA_PROJECTION%STARTING_XI,1))=DATA_PROJECTION%STARTING_XI(1:SIZE(DATA_PROJECTION%STARTING_XI,1))
            STARTING_XI(SIZE(DATA_PROJECTION%STARTING_XI,1):DATA_PROJECTION%NUMBER_OF_XI)=0.5_DP
          ELSE
            STARTING_XI(1:SIZE(DATA_PROJECTION%STARTING_XI,1))=DATA_PROJECTION%STARTING_XI(1:DATA_PROJECTION%NUMBER_OF_XI)
          ENDIF
          DEALLOCATE(DATA_PROJECTION%STARTING_XI)
          ALLOCATE(DATA_PROJECTION%STARTING_XI(DATA_PROJECTION%NUMBER_OF_XI),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data projection starting xi.",ERR,ERROR,*999)
          DATA_PROJECTION%STARTING_XI(1:DATA_PROJECTION%NUMBER_OF_XI)=STARTING_XI(1:DATA_PROJECTION%NUMBER_OF_XI)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_PROJECTION_TYPE_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_PROJECTION_TYPE_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_PROJECTION_TYPE_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_PROJECTION_TYPE_SET
    
   !
  !================================================================================================================================
  !
  
  !>Gets the relative tolerance for a data projection.
  SUBROUTINE DATA_PROJECTION_RELATIVE_TOLERANCE_GET(DATA_PROJECTION,RELATIVE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the relative tolerance for
    REAL(DP), INTENT(OUT) :: RELATIVE_TOLERANCE !<On exit, the relative tolerance of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_RELATIVE_TOLERANCE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE       
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_RELATIVE_TOLERANCE_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RELATIVE_TOLERANCE_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RELATIVE_TOLERANCE_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RELATIVE_TOLERANCE_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the relative tolerance for a data projection.
  SUBROUTINE DATA_PROJECTION_RELATIVE_TOLERANCE_SET(DATA_PROJECTION,RELATIVE_TOLERANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the relative tolerance for
    REAL(DP), INTENT(IN) :: RELATIVE_TOLERANCE !<the relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    CALL ENTERS("DATA_PROJECTION_RELATIVE_TOLERANCE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE      
        IF(RELATIVE_TOLERANCE>=0) THEN
          DATA_PROJECTION%RELATIVE_TOLERANCE=RELATIVE_TOLERANCE
        ELSE
          CALL FLAG_ERROR("Data projection relative tolerance must be a positive real number.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_RELATIVE_TOLERANCE_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RELATIVE_TOLERANCE_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RELATIVE_TOLERANCE_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RELATIVE_TOLERANCE_SET
  

  !
  !================================================================================================================================
  !
  
  !>Gets the starting xi for a data projection.
  SUBROUTINE DATA_PROJECTION_STARTING_XI_GET(DATA_PROJECTION,STARTING_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the starting xi for
    REAL(DP), INTENT(OUT) :: STARTING_XI(:) !<On exit, the starting xi of the specified data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    CHARACTER(LEN=MAXSTRLEN) :: LOCAL_ERROR
    
    CALL ENTERS("DATA_PROJECTION_STARTING_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(SIZE(STARTING_XI,1)>=SIZE(DATA_PROJECTION%STARTING_XI,1)) THEN
          STARTING_XI(1:SIZE(DATA_PROJECTION%STARTING_XI,1))=DATA_PROJECTION%STARTING_XI(1:SIZE(DATA_PROJECTION%STARTING_XI,1))
        ELSE
          WRITE(LOCAL_ERROR,'("The size of the supplied starting xi  array of ",I2," is too small. The size must be >= ",I2,".")' &
            & )SIZE(STARTING_XI,1),SIZE(DATA_PROJECTION%STARTING_XI,1)
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_STARTING_XI_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_STARTING_XI_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_STARTING_XI_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_STARTING_XI_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the starting xi for a data projection.
  SUBROUTINE DATA_PROJECTION_STARTING_XI_SET(DATA_PROJECTION,STARTING_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the starting xi for
    REAL(DP), INTENT(IN) :: STARTING_XI(:) !<the starting xi to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ni
    LOGICAL :: VALID_INPUT
    
    CALL ENTERS("DATA_PROJECTION_STARTING_XI_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        CALL FLAG_ERROR("Data projection have been finished.",ERR,ERROR,*999)
      ELSE      
        IF(SIZE(STARTING_XI,1)==SIZE(DATA_PROJECTION%STARTING_XI,1)) THEN
          VALID_INPUT=.TRUE.
          DO ni=1,SIZE(STARTING_XI,1)
            IF((STARTING_XI(ni)<0).OR.(STARTING_XI(ni)>1)) THEN
              VALID_INPUT=.FALSE.
              EXIT !this do
            ENDIF
          ENDDO
          IF(VALID_INPUT) THEN
            DATA_PROJECTION%STARTING_XI(1:SIZE(STARTING_XI))=STARTING_XI(1:SIZE(STARTING_XI))
          ELSE
            CALL FLAG_ERROR("Data projection starting xi must be between 0 and 1.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection starting xi dimension mismatch.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_STARTING_XI_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_STARTING_XI_SET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_STARTING_XI_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_STARTING_XI_SET
  
  !
  !================================================================================================================================
  !
  
  !>Sets the element for a data projection.
  SUBROUTINE DATA_PROJECTION_ELEMENT_SET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,ELEMENT_USER_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the element for
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<data point user number
    INTEGER(INTG), INTENT(IN) :: ELEMENT_USER_NUMBER !<the user element number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    LOGICAL :: DATA_POINT_EXISTS
    
    CALL ENTERS("DATA_PROJECTION_ELEMENT_SET",err,error,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      CALL DataProjection_DataPointCheckExist(DATA_PROJECTION,DATA_POINT_USER_NUMBER,DATA_POINT_EXISTS, &
        & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
      IF(DATA_POINT_EXISTS) THEN
        DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%ELEMENT_NUMBER=ELEMENT_USER_NUMBER
      ELSE
        CALL FLAG_ERROR("Data point with user number ("//TRIM(NUMBER_TO_VSTRING &
            & (DATA_POINT_USER_NUMBER,"*",ERR,ERROR))//") does not exist.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_ELEMENT_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_ELEMENT_SET",err,error)
    CALL EXITS("DATA_PROJECTION_ELEMENT_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_ELEMENT_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the element a data point is projected on.
  SUBROUTINE DataProjection_ElementGet(dataProjection,dataPointNumber,elementNumber,err,error,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<A pointer to the data projection to set the element for
    INTEGER(INTG), INTENT(IN) :: dataPointNumber !<data point number
    INTEGER(INTG), INTENT(OUT) :: elementNumber !<the element number to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL ENTERS("DataProjection_ElementGet",err,error,*999)

    IF(ASSOCIATED(dataProjection)) THEN
      IF(dataProjection%data_projection_finished) THEN
        IF((dataPointNumber<=dataProjection%data_points%number_of_data_points) .AND. (dataPointNumber>0))  THEN
          elementNumber = dataProjection%DATA_PROJECTION_RESULTS(dataPointNumber)%element_number
        ELSE
          CALL FLAG_ERROR("Data point number out of range.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data projection has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("DataProjection_ElementGet")
    RETURN
999 CALL ERRORS("DataProjection_ElementGet",err,error)    
    CALL EXITS("DataProjection_ElementGet")
    RETURN 1

  END SUBROUTINE DataProjection_ElementGet

  !
  !================================================================================================================================
  !
  
  !>Gets the distance from a data point projection.
  SUBROUTINE DataProjection_DistanceGet(dataProjection,dataPointNumber,distance,err,error,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: dataProjection !<A pointer to the data projection to get the distance for
    INTEGER(INTG), INTENT(IN) :: dataPointNumber !<data point number
    REAL(DP), INTENT(OUT) :: distance !<the distance to return
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    CALL ENTERS("DataProjection_DistanceGet",err,error,*999)

    IF(ASSOCIATED(dataProjection)) THEN
      IF(dataProjection%data_projection_finished) THEN
        IF((dataPointNumber<=dataProjection%data_points%number_of_data_points) .AND. (dataPointNumber>0))  THEN
          distance = dataProjection%DATA_PROJECTION_RESULTS(dataPointNumber)%distance
        ELSE
          CALL FLAG_ERROR("Data point number out of range.",err,error,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data projection has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",err,error,*999)
    ENDIF
    
    CALL EXITS("DataProjection_DistanceGet")
    RETURN
999 CALL ERRORS("DataProjection_DistanceGet",err,error)    
    CALL EXITS("DataProjection_DistanceGet")
    RETURN 1

  END SUBROUTINE DataProjection_DistanceGet

  !
  !================================================================================================================================
  !
  
  !>Gets the projection distance for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_DISTANCE_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_DISTANCE,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to get the projection distance for
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE !<On exit, the projection distance of the specified data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    
    CALL ENTERS("DATA_PROJECTION_RESULT_DISTANCE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
            & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
          PROJECTION_DISTANCE=DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%DISTANCE
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_RESULT_DISTANCE_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_DISTANCE_GET",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_RESULT_DISTANCE_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_DISTANCE_GET

  !
  !================================================================================================================================
  !

  !>Gets the label for a data projection for varying string labels. \see OPENCMISS::CMISSDataProjectionLabelGet
  SUBROUTINE DATA_PROJECTION_LABEL_GET_VS(DATA_PROJECTION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<the label to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      LABEL=DATA_PROJECTION%LABEL
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_LABEL_GET_VS")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_LABEL_GET_VS",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_LABEL_GET_VS")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Gets the label for a data projection for character labels. \see OPENCMISS::CMISSDataProjectionLabelGet
  SUBROUTINE DATA_PROJECTION_LABEL_GET_C(DATA_PROJECTION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<the label to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH
    
    CALL ENTERS("DATA_PROJECTION_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(DATA_PROJECTION%LABEL)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(LEN_TRIM(DATA_PROJECTION%LABEL))
      ELSE
        LABEL=CHAR(DATA_PROJECTION%LABEL,C_LENGTH)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_LABEL_GET_C")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_LABEL_GET_C",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_LABEL_GET_C")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_LABEL_GET_C

  !
  !================================================================================================================================
  !
  
  !>Sets the label for a data projection for varying string labels. \see OPENCMISS::CMISSDataProjectionLabelSet
  SUBROUTINE DATA_PROJECTION_LABEL_SET_C(DATA_PROJECTION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the label for
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<the label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      DATA_PROJECTION%LABEL=LABEL
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_LABEL_SET_C")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_LABEL_SET_C",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_LABEL_SET_C")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label for a data projection for varying string labels. \see OPENCMISS::CMISSDataProjectionLabelSet
  SUBROUTINE DATA_PROJECTION_LABEL_SET_VS(DATA_PROJECTION,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection to set the label for
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<the label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DATA_PROJECTION_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      DATA_PROJECTION%LABEL=LABEL
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_LABEL_SET_VS")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_LABEL_SET_VS",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_LABEL_SET_VS")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Gets the projection element number for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_ELEMENT_NUMBER, &
    & ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to get the projection element number for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER !<On exit, the projection element number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
   
    CALL ENTERS("DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
            & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
          PROJECTION_ELEMENT_NUMBER=DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%ELEMENT_NUMBER
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_ELEMENT_NUMBER_GET
  
  !
  !================================================================================================================================
  !
  
  !>Gets the projection element face number for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_ELEMENT_FACE_NUMBER &
    & ,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to get the projection element face number for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_FACE_NUMBER !<On exit, the projection element face number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    
    CALL ENTERS("DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          ! Check if boundary faces projection type was set
          IF(DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) THEN
            CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
              & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
            PROJECTION_ELEMENT_FACE_NUMBER=DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%ELEMENT_FACE_NUMBER
          ELSE
            CALL FLAG_ERROR("Data projection data projection projection type is not set to boundary faces projection type.", &
              & ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_ELEMENT_FACE_NUMBER_GET
  
  !
  !================================================================================================================================
  !
  
  !>Gets the projection element line number for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_ELEMENT_LINE_NUMBER &
    & ,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to get the element line number distance for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_LINE_NUMBER !<On exit, the projection element line number of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    
    CALL ENTERS("DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          ! Check if boundary lines projection type was set
          IF(DATA_PROJECTION%PROJECTION_TYPE==DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE) THEN
            CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
              & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
            PROJECTION_ELEMENT_LINE_NUMBER=DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%ELEMENT_LINE_NUMBER
          ELSE
            CALL FLAG_ERROR("Data projection data projection projection type is not set to boundary lines projection type.", &
              & ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)  
        ENDIF    
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_ELEMENT_LINE_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Gets the projection exit tag for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_EXIT_TAG_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_EXIT_TAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to get the projection exit tag for
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG !<On exit, the projection exit tag of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    
    CALL ENTERS("DATA_PROJECTION_RESULT_EXIT_TAG_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
            & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
          PROJECTION_EXIT_TAG=DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%EXIT_TAG
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)  
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_RESULT_EXIT_TAG_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_EXIT_TAG_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RESULT_EXIT_TAG_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_EXIT_TAG_GET

  !
  !================================================================================================================================
  !
  
  !>Gets the projection xi for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_XI_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to get the projection xi for
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(:) !<On exit, the projection xi of the specified global data point
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    
    CALL ENTERS("DATA_PROJECTION_RESULT_XI_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
            & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
          IF(SIZE(PROJECTION_XI,1)==SIZE(DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%XI,1)) THEN
            PROJECTION_XI=DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%XI
          ELSE
            CALL FLAG_ERROR("projection xi has size of "//TRIM(NUMBER_TO_VSTRING(SIZE(PROJECTION_XI,1),"*",ERR,ERROR))// &
              & "but it needs to have size of "// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(DATA_PROJECTION%DATA_PROJECTION_RESULTS &
              & (DATA_POINT_GLOBAL_NUMBER)%XI,1),"*",ERR,ERROR))// &
              & "." ,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)
        ENDIF   
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_RESULT_XI_GET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_XI_GET",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_RESULT_XI_GET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_XI_GET

  !
  !================================================================================================================================
  !
  
  !>Sets the projection xi for a data point identified by a given global number.
  SUBROUTINE DATA_PROJECTION_RESULT_XI_SET(DATA_PROJECTION,DATA_POINT_USER_NUMBER,PROJECTION_XI,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<A pointer to the data projection for which projection result is stored
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_USER_NUMBER !<The Data projection user number to set the projection xi for
    REAL(DP), INTENT(IN) :: PROJECTION_XI(:) !<The projection xi of the specified global data point to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DATA_POINT_GLOBAL_NUMBER
    
    CALL ENTERS("DATA_PROJECTION_RESULT_XI_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        IF(DATA_PROJECTION%DATA_PROJECTION_PROJECTED) THEN
          CALL DATA_PROJECTION_DATA_POINTS_GLOBAL_NUMBER_GET(DATA_PROJECTION,DATA_POINT_USER_NUMBER, &
            & DATA_POINT_GLOBAL_NUMBER,ERR,ERROR,*999)
          IF(SIZE(PROJECTION_XI,1)==SIZE(DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%XI,1)) THEN
            DATA_PROJECTION%DATA_PROJECTION_RESULTS(DATA_POINT_GLOBAL_NUMBER)%XI(1:SIZE(PROJECTION_XI,1))= &
              & PROJECTION_XI(1:SIZE(PROJECTION_XI,1))
          ELSE
            CALL FLAG_ERROR("projection xi has size of "//TRIM(NUMBER_TO_VSTRING(SIZE(PROJECTION_XI,1),"*",ERR,ERROR))// &
              & "but it needs to have size of "// &
              & TRIM(NUMBER_TO_VSTRING(SIZE(DATA_PROJECTION%DATA_PROJECTION_RESULTS &
              & (DATA_POINT_GLOBAL_NUMBER)%XI,1),"*",ERR,ERROR))// &
              & "." ,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data projection have not been projected.",ERR,ERROR,*999)
        ENDIF   
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DATA_PROJECTION_RESULT_XI_SET")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_RESULT_XI_SET",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_RESULT_XI_SET")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_RESULT_XI_SET

  !
  !================================================================================================================================
  !

END MODULE DATA_PROJECTION_ROUTINES

