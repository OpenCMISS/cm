!> \file
!> $Id: data_projection_routines.f90 660 2009-09-17 04:05:21Z chrispbradley $
!> \author Chris Bradley
!> \brief This module handles all data projection routines.
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
  USE MPI
  USE SORTING
  USE STRINGS
  USE TREES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE=1 !<The boundary line projection type for data projection, only projects to boundary lines of the mesh \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE=2 !<The boundary face projection type for data projection, only projects to boundary faces of the mesh \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE=3 !<The element projection type for data projection, projects to all elements in mesh \see DATA_PROJECTION_ROUTINES 
  
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_CONVERGED=1 !<Data projection exited due to it being converged \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_BOUNDS=2 !<Data projection exited due to it hitting the bound and continue to travel out of the element \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_MAX_ITERATION=3 !<Data projection exited due to it attaining maximum number of iteration specified by user \see DATA_PROJECTION_ROUTINES 
  INTEGER(INTG), PARAMETER :: DATA_PROJECTION_EXIT_TAG_NO_ELEMENT=4 !<Data projection exited due to no local element found, this happens when none of the candidate elements are within this computational node, and before MPI communication with other nodes \see DATA_PROJECTION_ROUTINES     
    
  !Module types

  !Module variables

  !Interfaces
  
  !INTERFACE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE
  !  MODULE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1
  !  MODULE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2
  !END INTERFACE !DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE
  

  PUBLIC DATA_PROJECTION_CREATE_FINISH,DATA_PROJECTION_CREATE_START
  
  PUBLIC DATA_PROJECTION_DESTROY
  
  PUBLIC DATA_PROJECTION_EVALUATE
 
CONTAINS

  !
  !================================================================================================================================CHECKED
  !
  
  !>Find the closest elements to a data point base on starting xi guess.
  SUBROUTINE DATA_PROJECTION_CLOSEST_ELEMENTS_FIND(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINT_NUMBER, &
    & CANDIDATE_ELEMENTS,NUMBER_OF_CANDIDATE_ELEMENTS,CLOSEST_ELEMENTS,CLOSEST_DISTANCES,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_NUMBER
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_CANDIDATE_ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: CLOSEST_ELEMENTS(:) !<List of shortest element number
    REAL(DP), INTENT(OUT) :: CLOSEST_DISTANCES(:) !<List of shortest distance (squared)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DATA_POINT_TYPE) :: DATA_POINT
    INTEGER(INTG) :: REGION_DIMENSIONS !<Region coordinate dimension
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_ELEMENTS !<Number of closest elements to record
    REAL(DP) :: DISTANCE_VECTOR(3) !<distance vector between data point and interpolated point, maximum dimension is 3
    REAL(DP) :: DISTANCE2 !<Distance squared
    INTEGER(INTG) :: nce
    INTEGER(INTG) :: element_number,insert_idx
      
    CALL ENTERS("DATA_PROJECTION_CLOSEST_ELEMENTS_FIND",ERR,ERROR,*999)
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN 
        DATA_POINT=DATA_PROJECTION%DATA_POINTS%DATA_POINTS(DATA_POINT_NUMBER)
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        NUMBER_OF_CLOSEST_ELEMENTS=SIZE(CLOSEST_ELEMENTS,1)    
        !loop through the first few elements
        DO nce=1,NUMBER_OF_CLOSEST_ELEMENTS
          element_number=CANDIDATE_ELEMENTS(nce)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_number, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999)
          DISTANCE_VECTOR(1:REGION_DIMENSIONS) = DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,1)
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
            CLOSEST_ELEMENTS(insert_idx)=element_number  
            EXIT                        
          ENDDO
        ENDDO  
        !Loop through the rest of the faces
        DO nce=NUMBER_OF_CLOSEST_ELEMENTS+1,NUMBER_OF_CANDIDATE_ELEMENTS
          element_number=CANDIDATE_ELEMENTS(nce)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_number, &
            & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,DATA_PROJECTION%STARTING_XI,INTERPOLATED_POINT,ERR,ERROR,*999) 
          DISTANCE_VECTOR(1:REGION_DIMENSIONS)=DATA_POINT%VALUES - INTERPOLATED_POINT%VALUES(:,1)
          DISTANCE2 = DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
          IF(DISTANCE2<CLOSEST_DISTANCES(NUMBER_OF_CLOSEST_ELEMENTS))THEN
            DO insert_idx=NUMBER_OF_CLOSEST_ELEMENTS,1,-1
              IF(insert_idx>1) THEN
                IF(DISTANCE2<CLOSEST_DISTANCES(insert_idx-1)) CYCLE
              ENDIF
              !insert the element into the correct index
              IF(insert_idx<NUMBER_OF_CLOSEST_ELEMENTS) THEN
                CLOSEST_DISTANCES(insert_idx+1:NUMBER_OF_CLOSEST_ELEMENTS)=CLOSEST_DISTANCES(insert_idx: &
                  & NUMBER_OF_CLOSEST_ELEMENTS-1)
                CLOSEST_ELEMENTS(insert_idx+1:NUMBER_OF_CLOSEST_ELEMENTS)=CLOSEST_ELEMENTS(insert_idx:NUMBER_OF_CLOSEST_ELEMENTS-1)
              ENDIF
              CLOSEST_DISTANCES(insert_idx)=DISTANCE2
              CLOSEST_ELEMENTS(insert_idx)=element_number  
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
  !================================================================================================================================CHECKED2
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
          IF(ASSOCIATED(DATA_PROJECTION%GEOMETRIC_FIELD)) THEN !Has to be associated
            IF(DATA_PROJECTION%GEOMETRIC_FIELD%FIELD_FINISHED) THEN !Has to be finished
              IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN !Cannot be finished
                CALL FLAG_ERROR("Data projection have already been finished.",ERR,ERROR,*999)
              ELSE
                DATA_PROJECTION%DATA_PROJECTION_FINISHED=.TRUE.
                CALL DATA_PROJECTION_DATA_POINTS_INITIALISE(DATA_PROJECTION%DATA_POINTS,ERR,ERROR,*999) !<Initialise data projection part in data points
              ENDIF    
            ELSE
              CALL FLAG_ERROR("Data projection geometric field have not been finished.",ERR,ERROR,*999)
            ENDIF      
          ELSE
            CALL FLAG_ERROR("Data projection geometric field is not associated.",ERR,ERROR,*999)
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
  !================================================================================================================================CHECKED2
  !  
  !>Starts the process of creating data projection.
  SUBROUTINE DATA_PROJECTION_CREATE_START(DATA_POINTS,GEOMETRIC_FIELD,DATA_PROJECTION,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points in which to create data projection    
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD !<A pointer to the geometric field to be interpolated
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<On exit, a pointer to the created data projection
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    !TYPE(DATA_PROJECTION_TYPE), POINTER :: NEW_DATA_PROJECTION
    !TYPE(DATA_PROJECTION_PTR_TYPE), POINTER :: NEW_DATA_PROJECTION_PTRS(:)
    !TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: INTERPOLATION_PARAMETERS    
    !TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    !INTEGER(INTG) :: TOTAL_NUMBER_OF_PROJECTION_ELEMENTS
    INTEGER(INTG) :: DATA_POINTS_REGION_DIMENSIONS,GEOMETRIC_FIELD_REGION_DIMENSIONS
    INTEGER(INTG) :: ni

    CALL ENTERS("DATA_PROJECTION_CREATE_START",ERR,ERROR,*999)

    NULLIFY(DATA_PROJECTION)
    IF(ASSOCIATED(DATA_POINTS)) THEN !Has to be associated
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN !Has to be finished
        IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN !Has to be associated
          IF(GEOMETRIC_FIELD%FIELD_FINISHED) THEN !Has to be finished
            IF(ASSOCIATED(DATA_PROJECTION)) THEN !Cannot be associated
              CALL FLAG_ERROR("Data projection is already associated.",ERR,ERROR,*999)
            ELSE
              DATA_POINTS_REGION_DIMENSIONS=DATA_POINTS%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              GEOMETRIC_FIELD_REGION_DIMENSIONS=GEOMETRIC_FIELD%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
              IF(DATA_POINTS_REGION_DIMENSIONS==GEOMETRIC_FIELD_REGION_DIMENSIONS) THEN !Dimension has to be equal
                IF(ASSOCIATED(DATA_POINTS%DATA_PROJECTION)) THEN
                  CALL FLAG_ERROR("Data points already has data projection associated.",ERR,ERROR,*999)
                ELSE
                  IF(ASSOCIATED(DATA_PROJECTION)) THEN
                    CALL FLAG_ERROR("Data projection is already associated.",ERR,ERROR,*999)
                  ELSE
                    ALLOCATE(DATA_POINTS%DATA_PROJECTION,STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data projection.",ERR,ERROR,*999)
                    DATA_POINTS%DATA_PROJECTION%DATA_PROJECTION_FINISHED=.FALSE.
                    DATA_POINTS%DATA_PROJECTION%DATA_POINTS=>DATA_POINTS
                    DATA_POINTS%DATA_PROJECTION%GEOMETRIC_FIELD=>GEOMETRIC_FIELD
                    DATA_POINTS%DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS=DATA_POINTS_REGION_DIMENSIONS
                    DATA_POINTS%DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE=0.5_DP
                    DATA_POINTS%DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS=25
                    DATA_POINTS%DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS=4
                    !Default always project to boundaries faces/lines when mesh dimension is equal to region dimension. If mesh dimension is less, project to all elements                
                    IF(GEOMETRIC_FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS<DATA_POINTS_REGION_DIMENSIONS) THEN
                      DATA_POINTS%DATA_PROJECTION%NUMBER_OF_XI=GEOMETRIC_FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS
                      DATA_POINTS%DATA_PROJECTION%PROJECTION_TYPE=DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE
                    ELSE
                      SELECT CASE(GEOMETRIC_FIELD%DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS)
                        CASE (2) 
                          DATA_POINTS%DATA_PROJECTION%NUMBER_OF_XI=1
                          DATA_POINTS%DATA_PROJECTION%PROJECTION_TYPE=DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE
                        CASE (3)
                          DATA_POINTS%DATA_PROJECTION%NUMBER_OF_XI=2
                          DATA_POINTS%DATA_PROJECTION%PROJECTION_TYPE=DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE
                        CASE DEFAULT
                          CALL FLAG_ERROR("Mesh dimensions out of bond [1,3].",ERR,ERROR,*999)
                        END SELECT
                    ENDIF
                    ALLOCATE(DATA_POINTS%DATA_PROJECTION%STARTING_XI(DATA_POINTS%DATA_PROJECTION%NUMBER_OF_XI),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data projection starting xi.",ERR,ERROR,*999)
                    DO ni=1,DATA_POINTS%DATA_PROJECTION%NUMBER_OF_XI
                      DATA_POINTS%DATA_PROJECTION%STARTING_XI(ni)=0.5_DP !<initialised to 0.5 in each xi direction
                    ENDDO !ni              
                    DATA_POINTS%DATA_PROJECTION%ABSOLUTE_TOLERANCE=1.0E-8_DP
                    DATA_POINTS%DATA_PROJECTION%RELATIVE_TOLERANCE=1.0E-6_DP
                    !Return the pointer        
                    DATA_PROJECTION=>DATA_POINTS%DATA_PROJECTION                
                  ENDIF !ASSOCIATED(DATA_PROJECTION)
                ENDIF !(ASSOCIATED(DATA_POINTS%DATA_PROJECTION))
              ELSE
                CALL FLAG_ERROR("Dimensions bewtween the geomtric field region and data points region does not match.", &
                  & ERR,ERROR,*999)        
              ENDIF !DATA_POINTS_REGION_DIMENSIONS==GEOMETRIC_FIELD_REGION_DIMENSIONS
            ENDIF !ASSOCIATED(DATA_PROJECTION)
          ELSE
            CALL FLAG_ERROR("Geometric field have not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Geometric field is not associated.",ERR,ERROR,*999)
        ENDIF !ASSOCIATED(GEOMETRIC_FIELD)
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF !ASSOCIATED(DATA_POINTS)
    
    CALL EXITS("DATA_PROJECTION_CREATE_START")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_CREATE_START",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_CREATE_START")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_CREATE_START 

  !
  !================================================================================================================================CHECKED2
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
      IF(ASSOCIATED(DATA_PROJECTION%DATA_POINTS)) NULLIFY(DATA_PROJECTION%DATA_POINTS)
      IF(ASSOCIATED(DATA_PROJECTION%GEOMETRIC_FIELD)) NULLIFY(DATA_PROJECTION%GEOMETRIC_FIELD)
      IF(ALLOCATED(DATA_PROJECTION%STARTING_XI)) DEALLOCATE(DATA_PROJECTION%STARTING_XI)
      DEALLOCATE(DATA_PROJECTION)
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

  !>Evaluates data projection.
  SUBROUTINE DATA_PROJECTION_EVALUATE(DATA_PROJECTION,ERR,ERROR,*) !optimising
    
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINTS(:)
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS
    !TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    !TYPE(DATA_PROJECTION_PROJECTED_POINTS_TYPE), POINTER :: PROJECTED_POINTS
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE,NUMBER_COMPUTATIONAL_NODES !<computational node/rank of the current process    
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,NUMBER_OF_DATA_POINTS,NUMBER_OF_ELEMENTS,NUMBER_OF_CANDIDATE_ELEMENTS
    INTEGER(INTG) :: NUMBER_OF_CLOSEST_ELEMENTS,TOTAL_NUMBER_OF_CLOSEST_ELEMENTS,REDUCED_NUMBER_OF_CLOSEST_ELEMENTS
    INTEGER(INTG), ALLOCATABLE :: CANDIDATE_ELEMENTS(:),CLOSEST_ELEMENTS(:,:),ALL_CLOSEST_ELEMENTS(:,:)
    REAL(DP), ALLOCATABLE :: CLOSEST_DISTANCES(:,:),ALL_CLOSEST_DISTANCES(:,:)
    INTEGER(INTG), ALLOCATABLE :: ALL_NUMBER_OF_CLOSEST_ELEMENTS(:)
    INTEGER(INTG), ALLOCATABLE :: ALL_MPI_DISPLACEMENTS(:),SORTING_IND_1(:),SORTING_IND_2(:)
    INTEGER(INTG), ALLOCATABLE :: ALL_NUMBER_OF_PROJECTED_POINTS(:)
    INTEGER(INTG) :: MPI_CLOSEST_ELEMENTS,MPI_CLOSEST_DISTANCES
    INTEGER(INTG) :: MPI_IERROR
    !LOGICAL :: ELEMENT_FOUND
    !LOGICAL, ALLOCATABLE :: PROJECTION_CONVERGED(:)
    INTEGER(INTG), ALLOCATABLE :: PROJECTED_ELEMENT(:),PROJECTION_EXIT_TAG(:)
    REAL(DP), ALLOCATABLE :: PROJECTED_DISTANCE(:,:),PROJECTED_XI(:,:)
    
    INTEGER(INTG) :: ne,ndp,ncn
    INTEGER(INTG) :: start_idx,finish_idx
    
    
    !INTEGER(INTG) :: NUMBER_OF_PROJECTED_POINTS,NUMBER_OF_PROJECTED_POINTS_2
  
    CALL ENTERS("DATA_PROJECTION_EVALUATE",ERR,ERROR,*999)
    
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(INTERPOLATED_POINTS)
    IF(ASSOCIATED(DATA_PROJECTION)) THEN
      IF(DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
        DATA_POINTS=>DATA_PROJECTION%DATA_POINTS
        !CALL DATA_PROJECTION_DATA_POINTS_INITIALISE(DATA_PROJECTION,PROJECTED_POINTS,ERR,ERROR,*999) !<Initialise data projection points
        !CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DATA_PROJECTION%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
        !  & INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(DATA_PROJECTION%GEOMETRIC_FIELD,INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
        CALL FIELD_INTERPOLATED_POINTS_INITIALISE(INTERPOLATION_PARAMETERS,INTERPOLATED_POINTS,ERR,ERROR,*999)
        INTERPOLATED_POINT=>INTERPOLATED_POINTS(FIELD_U_VARIABLE_TYPE)%PTR
        !DECOMPOSITION=>DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION  
        MESH_COMPONENT_NUMBER=DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN=>DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
        NUMBER_OF_DATA_POINTS=DATA_POINTS%NUMBER_OF_DATA_POINTS
        MY_COMPUTATIONAL_NODE=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
        NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
         !extract the boundary elements(faces/lines) if required
        NUMBER_OF_CANDIDATE_ELEMENTS=0
        SELECT CASE(DATA_PROJECTION%PROJECTION_TYPE)
          CASE (DATA_PROJECTION_BOUNDARY_LINES_PROJECTION_TYPE)
            !TODO extract boundary lines
        !    NUMBER_OF_DECOMPOSITION_ELEMENTS=DECOMPOSITION%TOPOLOGY%LINES%NUMBER_OF_LINES
        !    ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_DECOMPOSITION_ELEMENTS),STAT=ERR)
        !    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
        !    NUMBER_OF_CANDIDATE_ELEMENTS=0
        !    DO nde=1,NUMBER_OF_DECOMPOSITION_ELEMENTS
        !      IF(DECOMPOSITION%TOPOLOGY%LINES%LINES(nde)%BOUNDARY_LINE) THEN
        !        NUMBER_OF_CANDIDATE_ELEMENTS=NUMBER_OF_CANDIDATE_ELEMENTS+1
        !        CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATE_ELEMENTS)=nde
        !      ENDIF
        !    ENDDO
          CASE (DATA_PROJECTION_BOUNDARY_FACES_PROJECTION_TYPE) !identify the boundary faces
            !TODO extract boundary faces
            !NUMBER_OF_CANDIDATE_ELEMENTS=DECOMPOSITION%TOPOLOGY%FACES%NUMBER_OF_FACES
          CASE (DATA_PROJECTION_ALL_ELEMENTS_PROJECTION_TYPE)
            NUMBER_OF_ELEMENTS=DOMAIN%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS !non-ghost elements
            ALLOCATE(CANDIDATE_ELEMENTS(NUMBER_OF_ELEMENTS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate candidate elements.",ERR,ERROR,*999)
            DO ne=DOMAIN%MAPPINGS%ELEMENTS%INTERNAL_START,DOMAIN%MAPPINGS%ELEMENTS%INTERNAL_FINISH
              NUMBER_OF_CANDIDATE_ELEMENTS=NUMBER_OF_CANDIDATE_ELEMENTS+1
              CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATE_ELEMENTS)=ne
            ENDDO
            DO ne=DOMAIN%MAPPINGS%ELEMENTS%BOUNDARY_START,DOMAIN%MAPPINGS%ELEMENTS%BOUNDARY_FINISH
              NUMBER_OF_CANDIDATE_ELEMENTS=NUMBER_OF_CANDIDATE_ELEMENTS+1
              CANDIDATE_ELEMENTS(NUMBER_OF_CANDIDATE_ELEMENTS)=ne
            ENDDO
          CASE DEFAULT
            CALL FLAG_ERROR("No match for data projection type found",ERR,ERROR,*999)
        END SELECT !DATA_PROJECTION%PROJECTION_TYPE
        !find the clostest elements for each data point in the current computational node base on starting xi
        NUMBER_OF_CLOSEST_ELEMENTS=MIN(DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS,NUMBER_OF_CANDIDATE_ELEMENTS)
        ALLOCATE(CLOSEST_ELEMENTS(NUMBER_OF_DATA_POINTS,NUMBER_OF_CLOSEST_ELEMENTS),STAT=ERR)!the information for each data point has to be stored in the corresponding rows for them to be contiguous in memory for easy MPI access
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate closest elements.",ERR,ERROR,*999)
        ALLOCATE(CLOSEST_DISTANCES(NUMBER_OF_DATA_POINTS,NUMBER_OF_CLOSEST_ELEMENTS),STAT=ERR)!the information for each data point has to be stored in the corresponding rows for them to be contiguous in memory for easy MPI access
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate closest distances.",ERR,ERROR,*999)   
        DO ndp=1,NUMBER_OF_DATA_POINTS 
          CALL DATA_PROJECTION_CLOSEST_ELEMENTS_FIND(DATA_PROJECTION,INTERPOLATED_POINT,ndp,CANDIDATE_ELEMENTS, &
            &  NUMBER_OF_CANDIDATE_ELEMENTS,CLOSEST_ELEMENTS(ndp,:),CLOSEST_DISTANCES(ndp,:), &
            & ERR,ERROR,*999)
          CLOSEST_ELEMENTS(ndp,:)=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(CLOSEST_ELEMENTS(ndp,:)) !local to global element number mapping
        ENDDO !ndp
        !project the data points to each of the closest elements, use MPI if number of computational nodes is greater than 1
        IF(NUMBER_COMPUTATIONAL_NODES>1) THEN !use mpi
          !allocate arrays for mpi communication
          ALLOCATE(ALL_NUMBER_OF_CLOSEST_ELEMENTS(NUMBER_COMPUTATIONAL_NODES),STAT=ERR) 
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all number of closest elements.",ERR,ERROR,*999)
          ALLOCATE(ALL_MPI_DISPLACEMENTS(NUMBER_COMPUTATIONAL_NODES),STAT=ERR) 
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all displacements.",ERR,ERROR,*999)
          ALLOCATE(ALL_NUMBER_OF_PROJECTED_POINTS(NUMBER_COMPUTATIONAL_NODES),STAT=ERR) 
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all number of projected points.",ERR,ERROR,*999)
          ALLOCATE(PROJECTION_EXIT_TAG(NUMBER_OF_DATA_POINTS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected.",ERR,ERROR,*999)
          ALLOCATE(PROJECTED_ELEMENT(NUMBER_OF_DATA_POINTS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected element.",ERR,ERROR,*999)
          ALLOCATE(PROJECTED_DISTANCE(2,NUMBER_OF_DATA_POINTS),STAT=ERR) !PROJECTED_DISTANCE(2,:) stores the compuational node number, the information for each data point has to be stored in the corresponding column for MPI_ALLREDUCE with location return to work
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected distance.",ERR,ERROR,*999)
          ALLOCATE(PROJECTED_XI(DATA_PROJECTION%NUMBER_OF_XI,NUMBER_OF_DATA_POINTS),STAT=ERR) !the information for each data point is stored in the corresponding column to be consistent with PROJECTED_DISTANCE
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate projected.",ERR,ERROR,*999)
          ALLOCATE(SORTING_IND_2(NUMBER_OF_DATA_POINTS),STAT=ERR) 
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sorting ind 2.",ERR,ERROR,*999)          
          !gather and distribute the number of closest elements from all computational nodes
          CALL MPI_ALLGATHER(NUMBER_OF_CLOSEST_ELEMENTS,1,MPI_INTEGER,ALL_NUMBER_OF_CLOSEST_ELEMENTS,1,MPI_INTEGER, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
          TOTAL_NUMBER_OF_CLOSEST_ELEMENTS=SUM(ALL_NUMBER_OF_CLOSEST_ELEMENTS,1)
          !allocate arrays to store information gathered from all computational node
          ALLOCATE(ALL_CLOSEST_ELEMENTS(NUMBER_OF_DATA_POINTS,TOTAL_NUMBER_OF_CLOSEST_ELEMENTS),STAT=ERR) !the information for each data point is stored in the corresponding row so they are contiguous in memory for easy MPI access
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all closest elements.",ERR,ERROR,*999) 
          ALLOCATE(ALL_CLOSEST_DISTANCES(NUMBER_OF_DATA_POINTS,TOTAL_NUMBER_OF_CLOSEST_ELEMENTS),STAT=ERR) !the information for each data point is stored in the corresponding row so they are contiguous in memory for easy MPI access
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate all closest distances.",ERR,ERROR,*999)
          ALLOCATE(SORTING_IND_1(TOTAL_NUMBER_OF_CLOSEST_ELEMENTS),STAT=ERR) 
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate sorting ind 1.",ERR,ERROR,*999)
          !create and commit MPI_TYPE_CONTIGUOUS
          CALL MPI_TYPE_CONTIGUOUS(NUMBER_OF_DATA_POINTS,MPI_INTEGER,MPI_CLOSEST_ELEMENTS,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_TYPE_CONTIGUOUS",MPI_IERROR,ERR,ERROR,*999)
          CALL MPI_TYPE_COMMIT(MPI_CLOSEST_ELEMENTS,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPI_IERROR,ERR,ERROR,*999)        
          CALL MPI_TYPE_CONTIGUOUS(NUMBER_OF_DATA_POINTS,MPI_DOUBLE_PRECISION,MPI_CLOSEST_DISTANCES,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_TYPE_CONTIGUOUS",MPI_IERROR,ERR,ERROR,*999)        
          CALL MPI_TYPE_COMMIT(MPI_CLOSEST_DISTANCES,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPI_IERROR,ERR,ERROR,*999)
          !create displacement vectors for MPI_ALLGATHERV
          ALL_MPI_DISPLACEMENTS(1)=0
          DO ncn=1,(NUMBER_COMPUTATIONAL_NODES-1)
            ALL_MPI_DISPLACEMENTS(ncn+1)=ALL_MPI_DISPLACEMENTS(ncn)+ALL_NUMBER_OF_CLOSEST_ELEMENTS(ncn)
          ENDDO
          !gather and distribute the information of closest elements from all computational nodes
          CALL MPI_ALLGATHERV(CLOSEST_ELEMENTS(1,1),NUMBER_OF_CLOSEST_ELEMENTS,MPI_CLOSEST_ELEMENTS, &
            & ALL_CLOSEST_ELEMENTS,ALL_NUMBER_OF_CLOSEST_ELEMENTS,ALL_MPI_DISPLACEMENTS,MPI_CLOSEST_ELEMENTS, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
          CALL MPI_ALLGATHERV(CLOSEST_DISTANCES(1,1),NUMBER_OF_CLOSEST_ELEMENTS,MPI_CLOSEST_DISTANCES, &
            & ALL_CLOSEST_DISTANCES,ALL_NUMBER_OF_CLOSEST_ELEMENTS,ALL_MPI_DISPLACEMENTS,MPI_CLOSEST_DISTANCES, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
          REDUCED_NUMBER_OF_CLOSEST_ELEMENTS=MIN(DATA_PROJECTION%NUMBER_OF_CLOSEST_ELEMENTS,TOTAL_NUMBER_OF_CLOSEST_ELEMENTS) 
          PROJECTED_DISTANCE(2,:)=MY_COMPUTATIONAL_NODE
          SELECT CASE(DATA_PROJECTION%NUMBER_OF_XI)
            CASE (1) !1D mesh
              DO ndp=1,NUMBER_OF_DATA_POINTS
                CALL BUBBLE_ISORT(ALL_CLOSEST_DISTANCES(ndp,:),SORTING_IND_1,ERR,ERROR,*999)
                CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1(DATA_PROJECTION,INTERPOLATED_POINT,ndp, &
                  & ALL_CLOSEST_ELEMENTS(ndp,SORTING_IND_1(1:REDUCED_NUMBER_OF_CLOSEST_ELEMENTS)), &
                  & PROJECTION_EXIT_TAG(ndp),PROJECTED_ELEMENT(ndp),PROJECTED_DISTANCE(1,ndp),PROJECTED_XI(:,ndp),ERR,ERROR,*999)
                IF(PROJECTION_EXIT_TAG(ndp)==DATA_PROJECTION_EXIT_TAG_NO_ELEMENT) &
                  & PROJECTED_DISTANCE(1,ndp)=ALL_CLOSEST_DISTANCES(ndp,TOTAL_NUMBER_OF_CLOSEST_ELEMENTS) !assign a large value if no element is found in this computational node
              ENDDO
            CASE (2) !2D mesh
              !TODO project to 2D mesh
            CASE (3) !3D mesh
              !TODO project to 3D mesh
            CASE DEFAULT
              CALL FLAG_ERROR("Data projection number of xi is invalid",ERR,ERROR,*999)
          END SELECT !DATA_PROJECTION%NUMBER_OF_XI
          !find the shortest projection         
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,PROJECTED_DISTANCE,NUMBER_OF_DATA_POINTS,MPI_2DOUBLE_PRECISION,MPI_MINLOC, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
          !sort the computational node/rank from 0 to number of computational node
          CALL BUBBLE_ISORT(PROJECTED_DISTANCE(2,:),SORTING_IND_2,ERR,ERROR,*999)
          DO ncn=0,(NUMBER_COMPUTATIONAL_NODES-1)
            ALL_NUMBER_OF_PROJECTED_POINTS(ncn+1)=COUNT(PROJECTED_DISTANCE(2,:).EQ.ncn)
          ENDDO !ncn
          start_idx=SUM(ALL_NUMBER_OF_PROJECTED_POINTS(1:MY_COMPUTATIONAL_NODE))+1
          finish_idx=start_idx+ALL_NUMBER_OF_PROJECTED_POINTS(MY_COMPUTATIONAL_NODE+1)-1
          !create displacement vectors for MPI_ALLGATHERV          
          DO ncn=1,(NUMBER_COMPUTATIONAL_NODES-1)
            ALL_MPI_DISPLACEMENTS(ncn+1)=ALL_MPI_DISPLACEMENTS(ncn)+ALL_NUMBER_OF_PROJECTED_POINTS(ncn)
          ENDDO !ncn  
          !gather and distribute the shortest projection information
          CALL MPI_ALLGATHERV(PROJECTED_ELEMENT(SORTING_IND_2(start_idx:finish_idx)), &
            & ALL_NUMBER_OF_PROJECTED_POINTS(MY_COMPUTATIONAL_NODE+1),MPI_INTEGER,PROJECTED_ELEMENT, &
            & ALL_NUMBER_OF_PROJECTED_POINTS,ALL_MPI_DISPLACEMENTS,MPI_INTEGER, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)         
          CALL MPI_ALLGATHERV(PROJECTED_XI(:,SORTING_IND_2(start_idx:finish_idx)), &
            & ALL_NUMBER_OF_PROJECTED_POINTS(MY_COMPUTATIONAL_NODE+1),MPI_DOUBLE_PRECISION,PROJECTED_XI, &
            & ALL_NUMBER_OF_PROJECTED_POINTS,ALL_MPI_DISPLACEMENTS,MPI_DOUBLE_PRECISION, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)                           
          CALL MPI_ALLGATHERV(PROJECTION_EXIT_TAG(SORTING_IND_2(start_idx:finish_idx)), &
            & ALL_NUMBER_OF_PROJECTED_POINTS(MY_COMPUTATIONAL_NODE+1),MPI_INTEGER,PROJECTION_EXIT_TAG, &
            & ALL_NUMBER_OF_PROJECTED_POINTS,ALL_MPI_DISPLACEMENTS,MPI_INTEGER, &
            & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
          !assign projection information to projected points
          DO ndp=1,NUMBER_OF_DATA_POINTS          
            DATA_POINTS%DATA_POINTS(SORTING_IND_2(ndp))%PROJECTION_EXIT_TAG=PROJECTION_EXIT_TAG(ndp)
            DATA_POINTS%DATA_POINTS(SORTING_IND_2(ndp))%PROJECTION_ELEMENT_NUMBER=PROJECTED_ELEMENT(ndp)
            DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_DISTANCE=PROJECTED_DISTANCE(1,ndp)
            DATA_POINTS%DATA_POINTS(SORTING_IND_2(ndp))%PROJECTION_XI=PROJECTED_XI(:,ndp)
          ENDDO !ndp
        ELSE !no need to use mpi
          SELECT CASE(DATA_PROJECTION%NUMBER_OF_XI)
            CASE (1) !1D mesh
              DO ndp=1,NUMBER_OF_DATA_POINTS
                CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1(DATA_PROJECTION,INTERPOLATED_POINT,ndp,CLOSEST_ELEMENTS(ndp,:), &
                  & DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_EXIT_TAG,DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_ELEMENT_NUMBER, &
                  & DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_DISTANCE,DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_XI,ERR,ERROR,*999)
              ENDDO
            CASE (2) !2D mesh
              DO ndp=1,NUMBER_OF_DATA_POINTS
                CALL DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2(DATA_PROJECTION,INTERPOLATED_POINT,ndp,CLOSEST_ELEMENTS(ndp,:), &
                  & DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_EXIT_TAG,DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_ELEMENT_NUMBER, &
                  & DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_DISTANCE,DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_XI,ERR,ERROR,*999)
              ENDDO
            CASE (3) !3D mesh
              !TODO project to 3D mesh
            CASE DEFAULT
              CALL FLAG_ERROR("Data projection number of xi is invalid",ERR,ERROR,*999)
          END SELECT !DATA_PROJECTION%NUMBER_OF_XI
        ENDIF !NUMBER_COMPUTATIONAL_NODES>1
        DATA_PROJECTION%DATA_POINTS%DATA_POINTS_PROJECTED=.TRUE.               
      ELSE
        CALL FLAG_ERROR("Data projection have not been finished.",ERR,ERROR,*999)
      ENDIF !DATA_PROJECTION%DATA_PROJECTION_FINISHED
    ELSE
      CALL FLAG_ERROR("Data projection is not associated.",ERR,ERROR,*999)
    ENDIF !ASSOCIATED(DATA_PROJECTION) 
    !Test
    

    CALL EXITS("DATA_PROJECTION_EVALUATE")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_EVALUATE",ERR,ERROR)    
    CALL EXITS("DATA_PROJECTION_EVALUATE")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_EVALUATE
  
  !
  !================================================================================================================================
  !
  
  !>Find the projection of a data point onto 1D elements
  SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_1(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINT_NUMBER, CANDIDATE_ELEMENTS, &
    & PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_DISTANCE,PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_NUMBER
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(1)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(DATA_POINT_TYPE) :: DATA_POINT
    LOGICAL :: ELEMENT_EXISTS,INSIDE_REGION,CONVERGED
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,LOCAL_ELEMENT_NUMBER
    INTEGER(INTG) :: MESH_COMPONENT_NUMBER,REGION_DIMENSIONS
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
        MESH_COMPONENT_NUMBER=DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN_MAPPING=>DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DATA_POINT=DATA_PROJECTION%DATA_POINTS%DATA_POINTS(DATA_POINT_NUMBER)
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          GLOBAL_ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(DOMAIN_MAPPING,GLOBAL_ELEMENT_NUMBER,ELEMENT_EXISTS,LOCAL_ELEMENT_NUMBER, &
            & ERR,ERROR,*999)
          IF(ELEMENT_EXISTS) THEN !element exists
            IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_ELEMENT_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN !not a ghost element
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
              CONVERGED=.FALSE.
              DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,LOCAL_ELEMENT_NUMBER, &
                & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
              XI=DATA_PROJECTION%STARTING_XI
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
              DISTANCE_VECTOR=DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))       
              main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
                !Check for bounds [0,1]
                IF(XI(1)==0.0_DP) THEN
                  BOUND=-1 !bound at negative direction             
                ELSEIF(XI(1)==1.0_DP) THEN
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
                  CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999)
                  DISTANCE_VECTOR=DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
                  FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
                  CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test
                  IF(CONVERGED) EXIT !converged: exit inner loop first
                  IF(FUNCTION_VALUE_NEW>FUNCTION_VALUE) THEN !bad model: reduce step size
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
                !IF(.NOT.ELEMENT_FOUND) ELEMENT_FOUND=.TRUE.
                PROJECTION_EXIT_TAG=EXIT_TAG
                PROJECTION_ELEMENT_NUMBER=GLOBAL_ELEMENT_NUMBER
                PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
                PROJECTION_XI=XI
              ENDIF         
            ENDIF !not a ghost element
          ENDIF !local element exists
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
  SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_2(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINT_NUMBER, CANDIDATE_ELEMENTS, &
    & PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_DISTANCE,PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_NUMBER
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(2)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(DATA_POINT_TYPE) :: DATA_POINT
    LOGICAL :: ELEMENT_EXISTS,FREE,CONVERGED,INSIDE_REGION
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,LOCAL_ELEMENT_NUMBER
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
        MESH_COMPONENT_NUMBER=DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN_MAPPING=>DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
        REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DATA_POINT=DATA_PROJECTION%DATA_POINTS%DATA_POINTS(DATA_POINT_NUMBER)
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          GLOBAL_ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(DOMAIN_MAPPING,GLOBAL_ELEMENT_NUMBER,ELEMENT_EXISTS,LOCAL_ELEMENT_NUMBER, &
            & ERR,ERROR,*999)
          IF(ELEMENT_EXISTS) THEN !element exists
            IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_ELEMENT_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN !not a ghost element
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
              CONVERGED=.FALSE.
              DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet            
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,LOCAL_ELEMENT_NUMBER, &
                & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
              XI=DATA_PROJECTION%STARTING_XI
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
              DISTANCE_VECTOR=DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))   
              main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
                !Check for bounds [0,1]
                DO ni=1,2 
                  IF(XI(ni)==0.0_DP) THEN
                    BOUND(ni)=-1 !bound at negative direction             
                  ELSEIF(XI(ni)==1.0_DP) THEN
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
                  CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999)
                  DISTANCE_VECTOR=DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
                  FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR(1:REGION_DIMENSIONS),DISTANCE_VECTOR(1:REGION_DIMENSIONS))
                  CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test (before collision detection)
                  IF(CONVERGED) EXIT !converged: exit inner loop first
                  IF(FUNCTION_VALUE_NEW>FUNCTION_VALUE) THEN !bad model: reduce step size
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
                PROJECTION_ELEMENT_NUMBER=GLOBAL_ELEMENT_NUMBER
                PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
                PROJECTION_XI=XI
              ENDIF
            ENDIF !not a ghost element
          ENDIF !local element exists
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
  !================================================================================================================================CHECKED
  !
  
  !>Find the projection of a data point onto 3D elements
  SUBROUTINE DATA_PROJECTION_NEWTON_ELEMENTS_EVALUATE_3(DATA_PROJECTION,INTERPOLATED_POINT,DATA_POINT_NUMBER, CANDIDATE_ELEMENTS, &
    & PROJECTION_EXIT_TAG,PROJECTION_ELEMENT_NUMBER,PROJECTION_DISTANCE,PROJECTION_XI,ERR,ERROR,*)
    !Argument variables
    TYPE(DATA_PROJECTION_TYPE), POINTER :: DATA_PROJECTION !<Data projection problem to evaluate
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: INTERPOLATED_POINT    
    INTEGER(INTG), INTENT(IN) :: DATA_POINT_NUMBER
    INTEGER(INTG), INTENT(IN) :: CANDIDATE_ELEMENTS(:)
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_EXIT_TAG
    INTEGER(INTG), INTENT(OUT) :: PROJECTION_ELEMENT_NUMBER
    REAL(DP), INTENT(OUT) :: PROJECTION_DISTANCE
    REAL(DP), INTENT(OUT) :: PROJECTION_XI(3)
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
   
    !Local Variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(DATA_POINT_TYPE) :: DATA_POINT
    LOGICAL :: ELEMENT_EXISTS,FREE,CONVERGED,INSIDE_REGION
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,LOCAL_ELEMENT_NUMBER
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
        MESH_COMPONENT_NUMBER=DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER
        DOMAIN_MAPPING=>DATA_PROJECTION%GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
        !REGION_DIMENSIONS=DATA_PROJECTION%COORDINATE_SYSTEM_DIMENSIONS
        RELATIVE_TOLERANCE=DATA_PROJECTION%RELATIVE_TOLERANCE
        ABSOLUTE_TOLERANCE=DATA_PROJECTION%ABSOLUTE_TOLERANCE
        MAXIMUM_DELTA=DATA_PROJECTION%MAXIMUM_ITERATION_UPDATE
        MINIMUM_DELTA=0.025_DP*MAXIMUM_DELTA !need to set a minimum, in case if it gets too small      
        DATA_POINT=DATA_PROJECTION%DATA_POINTS%DATA_POINTS(DATA_POINT_NUMBER)
        DO ne=1,SIZE(CANDIDATE_ELEMENTS,1) !project on each candidate elements
          GLOBAL_ELEMENT_NUMBER=CANDIDATE_ELEMENTS(ne)
          CALL DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(DOMAIN_MAPPING,GLOBAL_ELEMENT_NUMBER,ELEMENT_EXISTS,LOCAL_ELEMENT_NUMBER, &
            & ERR,ERROR,*999)
          IF(ELEMENT_EXISTS) THEN !element exists
            IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_ELEMENT_NUMBER)%LOCAL_TYPE(1)/=DOMAIN_LOCAL_GHOST) THEN !not a ghost element
              EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
              CONVERGED=.FALSE.
              DELTA=0.5_DP*MAXIMUM_DELTA !start at half the MAXIMUM_DELTA as we do not know if quadratic model is a good approximation yet            
              CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,LOCAL_ELEMENT_NUMBER, &
                & INTERPOLATED_POINT%INTERPOLATION_PARAMETERS,ERR,ERROR,*999)
              XI=DATA_PROJECTION%STARTING_XI
              CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI,INTERPOLATED_POINT,ERR,ERROR,*999)
              DISTANCE_VECTOR=DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
              FUNCTION_VALUE=DOT_PRODUCT(DISTANCE_VECTOR,DISTANCE_VECTOR)   
              main_loop: DO itr1=1,DATA_PROJECTION%MAXIMUM_NUMBER_OF_ITERATIONS !(outer loop)
                !Check for bounds [0,1]
                DO ni=1,3
                  IF(XI(ni)==0.0_DP) THEN
                    BOUND(ni)=-1 !bound at negative direction             
                  ELSEIF(XI(ni)==1.0_DP) THEN
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
                IF(TEMP3>=0.0_DP) THEN !include>0 for numerical error
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
                  XI_UPDATE(1)=((HESSIAN_DIAGONAL(2)*HESSIAN_DIAGONAL(3)-FUNCTION_HESSIAN(2,3)**2)*FUNCTION_GRADIENT(1)+ &
                    & TEMP2*FUNCTION_GRADIENT(2)+TEMP3*FUNCTION_GRADIENT(3))/DET
                  XI_UPDATE(2)=((HESSIAN_DIAGONAL(1)*HESSIAN_DIAGONAL(3)-FUNCTION_HESSIAN(1,3)**2)*FUNCTION_GRADIENT(2)+ &                    
                    & TEMP2*FUNCTION_GRADIENT(1)+TEMP4*FUNCTION_GRADIENT(3))/DET
                  XI_UPDATE(3)=((HESSIAN_DIAGONAL(1)*HESSIAN_DIAGONAL(2)-FUNCTION_HESSIAN(1,2)**2)*FUNCTION_GRADIENT(3)+ &                    
                    & TEMP3*FUNCTION_GRADIENT(1)+TEMP4*FUNCTION_GRADIENT(2))/DET
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
                    IF((NBOUND==2).AND.(XI_UPDATE(nifix2(2))>XI_UPDATE(nifix2(1)))) nifix=nifix2(2) !only fix the direction that is most strongly suggesting leaving the element
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
                  CALL FIELD_INTERPOLATE_XI(SECOND_PART_DERIV,XI_NEW,INTERPOLATED_POINT,ERR,ERROR,*999)
                  DISTANCE_VECTOR=DATA_POINT%VALUES-INTERPOLATED_POINT%VALUES(:,NO_PART_DERIV)
                  FUNCTION_VALUE_NEW=DOT_PRODUCT(DISTANCE_VECTOR,DISTANCE_VECTOR)
                  CONVERGED=CONVERGED.AND.(DABS(FUNCTION_VALUE_NEW-FUNCTION_VALUE)/(1.0_DP+FUNCTION_VALUE)<RELATIVE_TOLERANCE) !second half of the convergence test (before collision detection)
                  IF(CONVERGED) EXIT !converged: exit inner loop first
                  IF(FUNCTION_VALUE_NEW>FUNCTION_VALUE) THEN !bad model: reduce step size
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
                PROJECTION_ELEMENT_NUMBER=GLOBAL_ELEMENT_NUMBER
                PROJECTION_DISTANCE=DSQRT(FUNCTION_VALUE)
                PROJECTION_XI=XI
              ENDIF
            ENDIF !not a ghost element
          ENDIF !local element exists
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
  !================================================================================================================================CHECKED2
  !
  
  !>Initialises the data projection part in a given data points.
  SUBROUTINE DATA_PROJECTION_DATA_POINTS_INITIALISE(DATA_POINTS,ERR,ERROR,*)

    !Argument variables
    TYPE(DATA_POINTS_TYPE), POINTER :: DATA_POINTS !<A pointer to the data points to initialise the data projection part for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ndp
    
    CALL ENTERS("DATA_PROJECTION_DATA_POINTS_INITIALISE",ERR,ERROR,*999)

    !NULLIFY(PROJECTED_POINTS)
    IF(ASSOCIATED(DATA_POINTS)) THEN
      IF(DATA_POINTS%DATA_POINTS_FINISHED) THEN
        IF(ASSOCIATED(DATA_POINTS%DATA_PROJECTION)) THEN
          IF(DATA_POINTS%DATA_PROJECTION%DATA_PROJECTION_FINISHED) THEN
            DO ndp=1,DATA_POINTS%NUMBER_OF_DATA_POINTS
              DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_COMPUTATIONAL_NODE_NUMBER=0
              DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_DISTANCE=0.0_DP
              DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_ELEMENT_NUMBER=0
              DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_EXIT_TAG=DATA_PROJECTION_EXIT_TAG_NO_ELEMENT
              ALLOCATE(DATA_POINTS%DATA_POINTS(ndp)%PROJECTION_XI(DATA_POINTS%DATA_PROJECTION%NUMBER_OF_XI),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate data points data points "// &
                & "("//TRIM(NUMBER_TO_VSTRING (ndp,"*",ERR,ERROR))//").",ERR,ERROR,*999)
            ENDDO !ndp
          ELSE
            CALL FLAG_ERROR("Data points data projection have not been finished.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Data points data projection is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Data points have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Data points is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DATA_PROJECTION_DATA_POINTS_INITIALISE")
    RETURN
999 CALL ERRORS("DATA_PROJECTION_DATA_POINTS_INITIALISE",ERR,ERROR)
    CALL EXITS("DATA_PROJECTION_DATA_POINTS_INITIALISE")
    RETURN 1

  END SUBROUTINE DATA_PROJECTION_DATA_POINTS_INITIALISE 
  
END MODULE DATA_PROJECTION_ROUTINES

