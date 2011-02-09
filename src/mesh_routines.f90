!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all mesh (node and element) routines.
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

!> This module handles all mesh (node and element) routines.
MODULE MESH_ROUTINES

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS_MPI
  USE CMISS_PARMETIS
  USE COMP_ENVIRONMENT
  USE COORDINATE_ROUTINES
  USE DOMAIN_MAPPINGS
  USE KINDS
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE LISTS
  USE MPI
  USE NODE_ROUTINES
  USE STRINGS
  USE TREES
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup MESH_ROUTINES_DecompositionTypes MESH_ROUTINES::DecompositionTypes
  !> \brief The Decomposition types parameters
  !> \see MESH_ROUTINES
  !>@{
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_ALL_TYPE=1 !<The decomposition contains all elements. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_CALCULATED_TYPE=2 !<The element decomposition is calculated by graph partitioning. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  INTEGER(INTG), PARAMETER :: DECOMPOSITION_USER_DEFINED_TYPE=3 !<The user will set the element decomposition. \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  !>Starts the process of creating a mesh
  INTERFACE MESH_CREATE_START
    MODULE PROCEDURE MESH_CREATE_START_INTERFACE
    MODULE PROCEDURE MESH_CREATE_START_REGION
  END INTERFACE !MESH_CREATE_START

  !>Initialises the meshes for a region or interface.
  INTERFACE MESHES_INITIALISE
    MODULE PROCEDURE MESHES_INITIALISE_INTERFACE
    MODULE PROCEDURE MESHES_INITIALISE_REGION
  END INTERFACE !MESHES_INITIALISE

  INTERFACE MESH_USER_NUMBER_FIND
    MODULE PROCEDURE MESH_USER_NUMBER_FIND_INTERFACE
    MODULE PROCEDURE MESH_USER_NUMBER_FIND_REGION
  END INTERFACE !MESH_USER_NUMBER_FIND

  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE

  PUBLIC DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH

  PUBLIC DECOMPOSITION_DESTROY

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_GET,DECOMPOSITION_ELEMENT_DOMAIN_SET

  PUBLIC DECOMPOSITION_MESH_COMPONENT_NUMBER_GET,DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET
  
  PUBLIC DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS
  
  PUBLIC DECOMPOSITION_TYPE_GET,DECOMPOSITION_TYPE_SET
  
  PUBLIC DECOMPOSITION_USER_NUMBER_FIND
  
  PUBLIC DECOMPOSITION_NODE_DOMAIN_GET

  PUBLIC DECOMPOSITION_CALCULATE_LINES_SET,DECOMPOSITION_CALCULATE_FACES_SET

  PUBLIC DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS
  
  PUBLIC MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS,MESH_TOPOLOGY_NODE_CHECK_EXISTS
  
  PUBLIC MESH_CREATE_START,MESH_CREATE_FINISH

  PUBLIC MESH_DESTROY
  
  PUBLIC MESH_NUMBER_OF_COMPONENTS_GET,MESH_NUMBER_OF_COMPONENTS_SET

  PUBLIC MESH_NUMBER_OF_ELEMENTS_GET,MESH_NUMBER_OF_ELEMENTS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_CREATE_START,MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH

  PUBLIC MESH_TOPOLOGY_ELEMENTS_DESTROY

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET

  PUBLIC MESH_USER_NUMBER_FIND
  
  PUBLIC MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  PUBLIC MESH_EMBEDDING_CREATE,MESH_EMBEDDING_SET_CHILD_NODE_POSITION
 
  PUBLIC MESH_EMBEDDING_SET_GAUSS_POINT_DATA

  PUBLIC MESHES_INITIALISE,MESHES_FINALISE

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finialises a decomposition adjacent element information and deallocates all memory
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(DECOMPOSITION_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    CALL EXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    CALL EXITS("DECOMPOSITION_ADJACENT_ELEMENT_FINALISE")
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises a decomposition adjacent element information.
  SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(DECOMPOSITION_ADJACENT_ELEMENT_TYPE) :: DECOMPOSITION_ADJACENT_ELEMENT !<The decomposition adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    DECOMPOSITION_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    CALL EXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    CALL EXITS("DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the creation of a domain decomposition on a given mesh. \see OPENCMISS::CMISSDecompositionCreateFinish
  SUBROUTINE DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(MESH_TYPE), POINTER :: MESH

    CALL ENTERS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      !Calculate which elements belong to which domain
      CALL DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the topology information for this decomposition
      CALL DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Initialise the domain for this computational node
      CALL DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
      !Calculate the decomposition topology
      CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
      DECOMPOSITION%DECOMPOSITION_FINISHED=.TRUE.
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      MESH=>DECOMPOSITION%MESH
      IF(ASSOCIATED(MESH)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Mesh = ",MESH%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of decompositions = ", &
          & MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS,ERR,ERROR,*999)
        DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Decomposition number = ",decomposition_no,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%GLOBAL_NUMBER,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
            & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR%USER_NUMBER,ERR,ERROR,*999)
        ENDDO !decomposition_no
      ELSE
        CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    CALL EXITS("DECOMPOSITION_CREATE_FINISH")
    RETURN
999 CALL ERRORS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_CREATE_FINISH")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Starts the creation of a domain decomposition for a given mesh. \see OPENCMISS::CMISSDecompositionCreateStart
  SUBROUTINE DECOMPOSITION_CREATE_START(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition 
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to decompose
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the created decomposition. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: NEW_DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITION)
    NULLIFY(NEW_DECOMPOSITIONS)

    CALL ENTERS("DECOMPOSITION_CREATE_START",ERR,ERROR,*999)

    NULLIFY(DECOMPOSITION)
    
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
          IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
            CALL DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
            IF(ASSOCIATED(DECOMPOSITION)) THEN
              LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
                & " has already been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              ALLOCATE(NEW_DECOMPOSITION,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decomposition.",ERR,ERROR,*999)
              !Set default decomposition properties
              NEW_DECOMPOSITION%GLOBAL_NUMBER=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1
              NEW_DECOMPOSITION%USER_NUMBER=USER_NUMBER
              NEW_DECOMPOSITION%DECOMPOSITION_FINISHED=.FALSE.
              NEW_DECOMPOSITION%CALCULATE_LINES=.TRUE. !Default
              NEW_DECOMPOSITION%CALCULATE_FACES=.FALSE. !Default
              NEW_DECOMPOSITION%DECOMPOSITIONS=>MESH%DECOMPOSITIONS
              NEW_DECOMPOSITION%MESH=>MESH
              NEW_DECOMPOSITION%MESH_COMPONENT_NUMBER=1
              !Default decomposition is all the mesh with one domain.
              NEW_DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
              NEW_DECOMPOSITION%NUMBER_OF_DOMAINS=1          
              ALLOCATE(NEW_DECOMPOSITION%ELEMENT_DOMAIN(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decomposition element domain.",ERR,ERROR,*999)
              NEW_DECOMPOSITION%ELEMENT_DOMAIN=0          
              !Nullify the domain
              NULLIFY(NEW_DECOMPOSITION%DOMAIN)
              !Nullify the topology
              NULLIFY(NEW_DECOMPOSITION%TOPOLOGY)
              !Add new decomposition into list of decompositions on the mesh
              ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decompositions.",ERR,ERROR,*999)
              DO decomposition_no=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
                NEW_DECOMPOSITIONS(decomposition_no)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_no)%PTR
              ENDDO !decomposition_no
              NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1)%PTR=>NEW_DECOMPOSITION
              IF(ASSOCIATED(MESH%DECOMPOSITIONS%DECOMPOSITIONS)) DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
              MESH%DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
              MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1        
              DECOMPOSITION=>NEW_DECOMPOSITION
            ENDIF
          ELSE
            LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
              & " are not associated."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITION)) THEN
      IF(ALLOCATED(NEW_DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(NEW_DECOMPOSITION%ELEMENT_DOMAIN)
      DEALLOCATE(NEW_DECOMPOSITION)
    ENDIF
    IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    NULLIFY(DECOMPOSITION)
    CALL ERRORS("DECOMPOSITION_CREATE_START",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_CREATE_START")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a user number and deallocates all memory. \see OPENCMISS::CMISSDecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY_NUMBER(USER_NUMBER,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to destroy.
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    CALL ENTERS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN

        !Find the decomposition identified by the user number
        FOUND=.FALSE.
        decomposition_position=0
        DO WHILE(decomposition_position<MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.FOUND)
          decomposition_position=decomposition_position+1
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          DECOMPOSITION=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_position)%PTR

          !Destroy all the decomposition components          
          IF(ALLOCATED(DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(DECOMPOSITION%ELEMENT_DOMAIN)
          CALL DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
          CALL DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
          
          DEALLOCATE(DECOMPOSITION)

          !Remove the decomposition from the list of decompositions
          IF(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
            ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decompositions.",ERR,ERROR,*999)
            DO decomposition_idx=1,MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
              IF(decomposition_idx<decomposition_position) THEN
                NEW_DECOMPOSITIONS(decomposition_idx)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ELSE IF(decomposition_idx>decomposition_position) THEN
                MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                  & MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
                NEW_DECOMPOSITIONS(decomposition_idx-1)%PTR=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
              ENDIF
            ENDDO !decomposition_idx
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
          ELSE
            DEALLOCATE(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
            MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Decomposition number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    CALL ERRORS("DECOMPOSITION_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_DESTROY_NUMBER")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys a domain decomposition identified by a pointer and deallocates all memory. \see OPENCMISS::CMISSDecompositionDestroy
  SUBROUTINE DECOMPOSITION_DESTROY(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    TYPE(DECOMPOSITIONS_TYPE), POINTER :: DECOMPOSITIONS
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    CALL ENTERS("DECOMPOSITION_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      DECOMPOSITIONS=>DECOMPOSITION%DECOMPOSITIONS
      IF(ASSOCIATED(DECOMPOSITIONS)) THEN
        decomposition_position=DECOMPOSITION%GLOBAL_NUMBER

        !Destroy all the decomposition components          
        IF(ALLOCATED(DECOMPOSITION%ELEMENT_DOMAIN)) DEALLOCATE(DECOMPOSITION%ELEMENT_DOMAIN)
        CALL DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
        CALL DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*999)
        
        DEALLOCATE(DECOMPOSITION)
        
        !Remove the decomposition from the list of decompositions
        IF(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>1) THEN
          ALLOCATE(NEW_DECOMPOSITIONS(DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decompositions.",ERR,ERROR,*999)
          DO decomposition_idx=1,DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS
            IF(decomposition_idx<decomposition_position) THEN
              NEW_DECOMPOSITIONS(decomposition_idx)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
            ELSE IF(decomposition_idx>decomposition_position) THEN
              DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER= &
                & DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%GLOBAL_NUMBER-1
              NEW_DECOMPOSITIONS(decomposition_idx-1)%PTR=>DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
            ENDIF
          ENDDO !decomposition_idx
          DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
          DECOMPOSITIONS%DECOMPOSITIONS=>NEW_DECOMPOSITIONS
          DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS-1
        ELSE
          DEALLOCATE(DECOMPOSITIONS%DECOMPOSITIONS)
          DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition decompositions is not associated.",ERR,ERROR,*999)
       ENDIF
    ELSE
      CALL FLAG_ERROR("Decompositions is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_DECOMPOSITIONS)) DEALLOCATE(NEW_DECOMPOSITIONS)
    CALL ERRORS("DECOMPOSITION_DESTROY",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_DESTROY")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_DESTROY

  !
  !================================================================================================================================
  !

  !>Calculates the element domains for a decomposition of a mesh. \see OPENCMISS::CMISSDecompositionElementDomainCalculate
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to calculate the element domains for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_elem_indicies,elem_index,elem_count,ne,nn,my_computational_node_number,number_computational_nodes, &
      & no_computational_node,ELEMENT_START,ELEMENT_STOP,MY_ELEMENT_START,MY_ELEMENT_STOP,NUMBER_OF_ELEMENTS, &
      & MY_NUMBER_OF_ELEMENTS,MPI_IERROR,MAX_NUMBER_ELEMENTS_PER_NODE,component_idx
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_COUNT(:),ELEMENT_PTR(:),ELEMENT_INDICIES(:),ELEMENT_DISTANCE(:),DISPLACEMENTS(:), &
      & RECEIVE_COUNTS(:)
    INTEGER(INTG) :: ELEMENT_WEIGHT(1),WEIGHT_FLAG,NUMBER_FLAG,NUMBER_OF_CONSTRAINTS, &
      & NUMBER_OF_COMMON_NODES,PARMETIS_OPTIONS(0:2)
    REAL(SP) :: UBVEC(1)
    REAL(SP), ALLOCATABLE :: TPWGTS(:)
    REAL(DP) :: NUMBER_ELEMENTS_PER_NODE
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN

          component_idx=DECOMPOSITION%MESH_COMPONENT_NUMBER
          
          number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
          IF(ERR/=0) GOTO 999
          my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
          IF(ERR/=0) GOTO 999
          
          SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
          CASE(DECOMPOSITION_ALL_TYPE)
            !Do nothing. Decomposition checked below.
           CASE(DECOMPOSITION_CALCULATED_TYPE)
            !Calculate the general decomposition

            IF(DECOMPOSITION%NUMBER_OF_DOMAINS==1) THEN
              DECOMPOSITION%ELEMENT_DOMAIN=0
            ELSE
              number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              
              NUMBER_ELEMENTS_PER_NODE=REAL(MESH%NUMBER_OF_ELEMENTS,DP)/REAL(number_computational_nodes,DP)
              ELEMENT_START=1
              ELEMENT_STOP=0
              MAX_NUMBER_ELEMENTS_PER_NODE=-1
              ALLOCATE(RECEIVE_COUNTS(0:number_computational_nodes-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate recieve counts.",ERR,ERROR,*999)
              ALLOCATE(DISPLACEMENTS(0:number_computational_nodes-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate displacements.",ERR,ERROR,*999)
              ALLOCATE(ELEMENT_DISTANCE(0:number_computational_nodes),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element distance.",ERR,ERROR,*999)
              ELEMENT_DISTANCE(0)=0
              DO no_computational_node=0,number_computational_nodes-1
                ELEMENT_START=ELEMENT_STOP+1
                IF(no_computational_node==number_computational_nodes-1) THEN
                  ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS
                ELSE
                  ELEMENT_STOP=ELEMENT_START+NINT(NUMBER_ELEMENTS_PER_NODE,INTG)-1
                ENDIF
                IF((number_computational_nodes-1-no_computational_node)>(MESH%NUMBER_OF_ELEMENTS-ELEMENT_STOP)) &
                  & ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS-(number_computational_nodes-1-no_computational_node)
                IF(ELEMENT_START>MESH%NUMBER_OF_ELEMENTS) ELEMENT_START=MESH%NUMBER_OF_ELEMENTS
                IF(ELEMENT_STOP>MESH%NUMBER_OF_ELEMENTS) ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS
                DISPLACEMENTS(no_computational_node)=ELEMENT_START-1
                ELEMENT_DISTANCE(no_computational_node+1)=ELEMENT_STOP !C numbering
                NUMBER_OF_ELEMENTS=ELEMENT_STOP-ELEMENT_START+1
                RECEIVE_COUNTS(no_computational_node)=NUMBER_OF_ELEMENTS
                IF(NUMBER_OF_ELEMENTS>MAX_NUMBER_ELEMENTS_PER_NODE) MAX_NUMBER_ELEMENTS_PER_NODE=NUMBER_OF_ELEMENTS
                IF(no_computational_node==my_computational_node_number) THEN
                  MY_ELEMENT_START=ELEMENT_START
                  MY_ELEMENT_STOP=ELEMENT_STOP
                  MY_NUMBER_OF_ELEMENTS=ELEMENT_STOP-ELEMENT_START+1
                  number_elem_indicies=0
                  DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
                    BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                    number_elem_indicies=number_elem_indicies+BASIS%NUMBER_OF_NODES
                  ENDDO !ne
                ENDIF
              ENDDO !no_computational_node
              
              ALLOCATE(ELEMENT_PTR(0:MY_NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element pointer list.",ERR,ERROR,*999)
              ALLOCATE(ELEMENT_INDICIES(0:number_elem_indicies-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element indicies list.",ERR,ERROR,*999)
              ALLOCATE(TPWGTS(1:DECOMPOSITION%NUMBER_OF_DOMAINS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate tpwgts.",ERR,ERROR,*999)
              elem_index=0
              elem_count=0
              ELEMENT_PTR(0)=0
              DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
                elem_count=elem_count+1
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                DO nn=1,BASIS%NUMBER_OF_NODES
                  ELEMENT_INDICIES(elem_index)=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)% &
                    & MESH_ELEMENT_NODES(nn)-1 !C numbering
                  elem_index=elem_index+1
                ENDDO !nn
                ELEMENT_PTR(elem_count)=elem_index !C numbering
              ENDDO !ne
              
              !Set up ParMETIS variables
              WEIGHT_FLAG=0 !No weights
              ELEMENT_WEIGHT(1)=1 !Isn't used due to weight flag
              NUMBER_FLAG=0 !C Numbering as there is a bug with Fortran numbering
              NUMBER_OF_CONSTRAINTS=1
              NUMBER_OF_COMMON_NODES=2
              TPWGTS=1.0_SP/REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,SP)
              UBVEC=1.05_SP
              PARMETIS_OPTIONS(0)=1 !If zero, defaults are used, otherwise next two values are used
              PARMETIS_OPTIONS(1)=7 !Level of information to output
              PARMETIS_OPTIONS(2)=CMISS_RANDOM_SEEDS(1) !Seed for random number generator
              
              !Call ParMETIS to calculate the partitioning of the mesh graph.
              CALL PARMETIS_PARTMESHKWAY(ELEMENT_DISTANCE,ELEMENT_PTR,ELEMENT_INDICIES,ELEMENT_WEIGHT,WEIGHT_FLAG,NUMBER_FLAG, &
                & NUMBER_OF_CONSTRAINTS,NUMBER_OF_COMMON_NODES,DECOMPOSITION%NUMBER_OF_DOMAINS,TPWGTS,UBVEC,PARMETIS_OPTIONS, &
                & DECOMPOSITION%NUMBER_OF_EDGES_CUT,DECOMPOSITION%ELEMENT_DOMAIN(DISPLACEMENTS(my_computational_node_number)+1:), &
                & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,ERR,ERROR,*999)
              
              !Transfer all the element domain information to the other computational nodes so that each rank has all the info
              IF(number_computational_nodes>1) THEN
                !This should work on a single processor but doesn't for mpich2 under windows. Maybe a bug? Avoid for now.
                CALL MPI_ALLGATHERV(MPI_IN_PLACE,MAX_NUMBER_ELEMENTS_PER_NODE,MPI_INTEGER,DECOMPOSITION%ELEMENT_DOMAIN, &
                  & RECEIVE_COUNTS,DISPLACEMENTS,MPI_INTEGER,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
                CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
              ENDIF
              
              DEALLOCATE(DISPLACEMENTS)
              DEALLOCATE(RECEIVE_COUNTS)
              DEALLOCATE(ELEMENT_DISTANCE)
              DEALLOCATE(ELEMENT_PTR)
              DEALLOCATE(ELEMENT_INDICIES)
              DEALLOCATE(TPWGTS)

            ENDIF
            
          CASE(DECOMPOSITION_USER_DEFINED_TYPE)
            !Do nothing. Decomposition checked below.          
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid domain decomposition type.",ERR,ERROR,*999)            
          END SELECT

          !Check decomposition and check that each domain has an element in it.
          ALLOCATE(ELEMENT_COUNT(0:number_computational_nodes-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element count.",ERR,ERROR,*999)
          ELEMENT_COUNT=0
          DO elem_index=1,MESH%NUMBER_OF_ELEMENTS
            no_computational_node=DECOMPOSITION%ELEMENT_DOMAIN(elem_index)
            IF(no_computational_node>=0.AND.no_computational_node<number_computational_nodes) THEN
              ELEMENT_COUNT(no_computational_node)=ELEMENT_COUNT(no_computational_node)+1
            ELSE
              LOCAL_ERROR="The computational node number of "//TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
                & " for element number "//TRIM(NUMBER_TO_VSTRING(elem_index,"*",ERR,ERROR))// &
                & " is invalid. The computational node number must be between 0 and "// &
                & TRIM(NUMBER_TO_VSTRING(number_computational_nodes-1,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !elem_index
          DO no_computational_node=0,number_computational_nodes-1
            IF(ELEMENT_COUNT(no_computational_node)==0) THEN
              LOCAL_ERROR="Invalid decomposition. There are no elements in computational node "// &
                & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !no_computational_node
          DEALLOCATE(ELEMENT_COUNT)
          
        ELSE
          CALL FLAG_ERROR("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition for mesh number ",DECOMPOSITION%MESH%USER_NUMBER, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains = ", DECOMPOSITION%NUMBER_OF_DOMAINS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Element domains:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Decomposition type = ", DECOMPOSITION%DECOMPOSITION_TYPE, &
        & ERR,ERROR,*999)
      IF(DECOMPOSITION%DECOMPOSITION_TYPE==DECOMPOSITION_CALCULATED_TYPE) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of edges cut = ",DECOMPOSITION%NUMBER_OF_EDGES_CUT, &
          & ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of elements = ",DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS, &
        & ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Domain = ",DECOMPOSITION%ELEMENT_DOMAIN(ne), &
          & ERR,ERROR,*999)
      ENDDO !ne
    ENDIF
    
    CALL EXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE")
    RETURN
999 IF(ALLOCATED(RECEIVE_COUNTS)) DEALLOCATE(RECEIVE_COUNTS)
    IF(ALLOCATED(DISPLACEMENTS)) DEALLOCATE(DISPLACEMENTS)
    IF(ALLOCATED(ELEMENT_DISTANCE)) DEALLOCATE(ELEMENT_DISTANCE)
    IF(ALLOCATED(ELEMENT_PTR)) DEALLOCATE(ELEMENT_PTR)
    IF(ALLOCATED(ELEMENT_INDICIES)) DEALLOCATE(ELEMENT_INDICIES)
    IF(ALLOCATED(TPWGTS)) DEALLOCATE(TPWGTS)
    CALL ERRORS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Gets the domain for a given element in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OPENCMISS::CMISSDecompositionElementDomainGet
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET(DECOMPOSITION,USER_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables`
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS


    CALL ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR,*999)

    GLOBAL_ELEMENT_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
            IF(ASSOCIATED(MESH_ELEMENTS)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
                  DOMAIN_NUMBER=DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Decomposition mesh element corresponding to user number not found.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Decomposition mesh elements are not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_ELEMENT_DOMAIN_GET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_ELEMENT_DOMAIN_GET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_GET
  
  !
  !================================================================================================================================
  !

  !>Sets the domain for a given element in a decomposition of a mesh. \todo move to user number, should be able to specify lists of elements. \see OPENCMISS::CMISSDecompositionElementDomainSet 
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET(DECOMPOSITION,GLOBAL_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: GLOBAL_ELEMENT_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: DOMAIN_NUMBER !<The domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: number_computational_nodes
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            IF(GLOBAL_ELEMENT_NUMBER>0.AND.GLOBAL_ELEMENT_NUMBER<=MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
              number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              IF(DOMAIN_NUMBER>=0.AND.DOMAIN_NUMBER<number_computational_nodes) THEN
                DECOMPOSITION%ELEMENT_DOMAIN(GLOBAL_ELEMENT_NUMBER)=DOMAIN_NUMBER
              ELSE
                LOCAL_ERROR="Domain number "//TRIM(NUMBER_TO_VSTRING(DOMAIN_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid. The limits are 0 to "//TRIM(NUMBER_TO_VSTRING(number_computational_nodes,"*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The limits are 1 to "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Decomposition mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_ELEMENT_DOMAIN_SET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::CMISSDecompositionMeshComponentGet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the mesh component for
    INTEGER(INTG), INTENT(OUT) :: MESH_COMPONENT_NUMBER !<On return, the mesh component number to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          MESH_COMPONENT_NUMBER=DECOMPOSITION%MESH_COMPONENT_NUMBER
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_GET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_GET
  
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the mesh component number which will be used for the decomposition of a mesh. \see OPENCMISS::CMISSDecompositionMeshComponentSet
  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the mesh component for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS) THEN
            DECOMPOSITION%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
          ELSE
            LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
              & "is invalid. The component number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  !
  !================================================================================================================================
  !

  !!MERGE: ditto
  
  !>Gets the number of domains for a decomposition. \see OPENCMISS::CMISSDecompositionNumberOfDomainsGet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the number of domains for.
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_DOMAINS !<On return, the number of domains to get.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_DOMAINS=DECOMPOSITION%NUMBER_OF_DOMAINS
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_GET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_GET
  

  !
  !================================================================================================================================
  !

  !>Sets/changes the number of domains for a decomposition. \see OPENCMISS::CMISSDecompositionNumberOfDomainsSet
  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the number of domains for.
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(NUMBER_OF_DOMAINS==1) THEN
            DECOMPOSITION%NUMBER_OF_DOMAINS=1
          ELSE
            CALL FLAG_ERROR("Can only have one domain for all decomposition type.",ERR,ERROR,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE)
          IF(NUMBER_OF_DOMAINS>=1) THEN
            !wolfye???<=?
            IF(NUMBER_OF_DOMAINS<=DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS) THEN
              !Get the number of computational nodes
              NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999
              !!TODO: relax this later
              !IF(NUMBER_OF_DOMAINS==NUMBER_COMPUTATIONAL_NODES) THEN
                DECOMPOSITION%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS             
              !ELSE
              !  LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINSS,"*",ERR,ERROR))// &
              !    & ") is not equal to the number of computational nodes ("// &
              !    & TRIM(NUMBER_TO_VSTRING(NUMBER_COMPUTATIONAL_NODES,"*",ERR,ERROR))//")"
              !  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              !ENDIF
            ELSE
              LOCAL_ERROR="The number of domains ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
                & ") must be <= the number of global elements ("// &
                & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//") in the mesh."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
             CALL FLAG_ERROR("Number of domains must be >= 1.",ERR,ERROR,*999)
           ENDIF
         CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%DECOMPOSITION_TYPE,"*",ERR,ERROR))// &
            & " is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET

  !
  !================================================================================================================================
  !

  !>Calculates the topology for a decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finish creating.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      IF(DECOMPOSITION%CALCULATE_LINES)THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      !Calculate the face topology
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_CALCULATE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_CALCULATE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a decomposition. 
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(DECOMPOSITION_TOPOLOGY,USER_ELEMENT_NUMBER,ELEMENT_EXISTS, &
    & DECOMPOSITION_LOCAL_ELEMENT_NUMBER,GHOST_ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: DECOMPOSITION_TOPOLOGY !<A pointer to the decomposition topology to check the element exists on
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NUMBER !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: ELEMENT_EXISTS !<On exit, is .TRUE. if the element user number exists in the decomposition topology, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DECOMPOSITION_LOCAL_ELEMENT_NUMBER !<On exit, if the element exists the local number corresponding to the user element number. If the element does not exist then local number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_ELEMENT !<On exit, is .TRUE. if the local element (if it exists) is a ghost element, .FALSE. if not.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR,*999)

    ELEMENT_EXISTS=.FALSE.
    DECOMPOSITION_LOCAL_ELEMENT_NUMBER=0
    GHOST_ELEMENT=.FALSE.
    IF(ASSOCIATED(DECOMPOSITION_TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>DECOMPOSITION_TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,TREE_NODE,DECOMPOSITION_LOCAL_ELEMENT_NUMBER,ERR,ERROR,*999)
          ELEMENT_EXISTS=.TRUE.
          GHOST_ELEMENT=DECOMPOSITION_LOCAL_ELEMENT_NUMBER>DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS")
    RETURN 1
    
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic !\todo add comment

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL DECOMPOSITION_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
    IF(ALLOCATED(ELEMENT%ELEMENT_LINES)) DEALLOCATE(ELEMENT%ELEMENT_LINES)
    IF(ALLOCATED(ELEMENT%ELEMENT_FACES)) DEALLOCATE(ELEMENT%ELEMENT_FACES)
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)
 
    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%LOCAL_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
  
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element numbers adjacent to an element in a decomposition topology.
  SUBROUTINE DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the adjacent elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: candidate_idx,face_node_idx,node_idx,surrounding_el_idx,candidate_el,idx
    INTEGER(INTG) :: SURROUNDING_ELEMENTS(100) !Fixed size array... should be large enough
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4)
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    CALL ENTERS("DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      IF(ASSOCIATED(DECOMPOSITION)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                  !Loop over the elements in the decomposition
                  DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                    !Create a list for every xi direction (plus and minus)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
                      CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
                    ENDDO !nic
                    NUMBER_OF_NODES_XIC=1
                    NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)= &
                      & BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
                    !Place the current element in the surrounding list
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER, &
                      & ERR,ERROR,*999)
                    SELECT CASE(BASIS%TYPE)
                    CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
!!TODO: Calculate this and set it as part of the basis type
                      !Determine the collapsed "faces" if any
                      NODE_POSITION_INDEX=1
                      !Loop over the face normals of the element
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the face xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Reset the node_position_index in this xi direction
                        NODE_POSITION_INDEX(ni)=1
                        !Loop over the two faces with this normal
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni
                          FACE_COLLAPSED(xi_direction)=.FALSE.
                          DO j=1,2
                            xi_dir_check=FACE_XI(j)
                            IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                              xi_dir_search=FACE_XI(3-j)
                              NODE_POSITION_INDEX(xi_dir_search)=1
                              XI_COLLAPSED=.TRUE.
                              DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                                !Get the first local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=1
                                nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                !Get the second local node along the xi check direction
                                NODE_POSITION_INDEX(xi_dir_check)=2
                                nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                  & NODE_POSITION_INDEX(3),1)
                                IF(nn1/=0.AND.nn2/=0) THEN
                                  IF(DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn1)/= &
                                    & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn2)) XI_COLLAPSED=.FALSE.
                                ENDIF
                                NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                              ENDDO !xi_dir_search
                              IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                            ENDIF
                          ENDDO !j
                          NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                        ENDDO !direction_index
                      ENDDO !ni
                      !Loop over the xi directions and calculate the surrounding elements
                      DO ni=1,BASIS%NUMBER_OF_XI
                        !Determine the xi directions that lie in this xi direction
                        FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                        FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                        !Loop over the two faces
                        DO direction_index=-1,1,2
                          xi_direction=direction_index*ni                  
                          !Find nodes in the element on the appropriate face/line/point
                          NULLIFY(NODE_MATCH_LIST)
                          CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                          IF(direction_index==-1) THEN
                            NODE_POSITION_INDEX(ni)=1
                          ELSE
                            NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                          ENDIF
                          !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
                          !also collapsed. This may indicate that we have a funny element in non-rc coordinates that goes around the
                          !central axis back to itself
                          IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                            !Do nothing - the match lists are already empty
                          ELSE
                            !Find the nodes to match and add them to the node match list
                            SELECT CASE(BASIS%NUMBER_OF_XI)
                            CASE(1)
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),1,1,1)
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            CASE(2)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),1,1)
                                IF(nn/=0) THEN
                                  np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                  CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                ENDIF
                              ENDDO !nn1
                            CASE(3)
                              DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1)),NUMBER_OF_NODES_XIC(FACE_XI(1))-1
                                NODE_POSITION_INDEX(FACE_XI(1))=nn1
                                DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2)),NUMBER_OF_NODES_XIC(FACE_XI(2))-1
                                  NODE_POSITION_INDEX(FACE_XI(2))=nn2
                                  nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                    & NODE_POSITION_INDEX(3),1)
                                  IF(nn/=0) THEN
                                    np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                    CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                                  ENDIF
                                ENDDO !nn2
                              ENDDO !nn1
                            CASE DEFAULT
                              LOCAL_ERROR="The number of xi directions in the basis of "// &
                                & TRIM(NUMBER_TO_VSTRING(BASIS%NUMBER_OF_XI,"*",ERR,ERROR))//" is invalid."
                              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                            END SELECT
                          ENDIF
                          CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                          CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=0
                          IF(NUMBER_NODE_MATCHES>0) THEN
                            !NODE_MATCHES now contain the list of corner nodes in the current face with normal_xi of ni.
                            !Look at the surrounding elements of each of these nodes, if there is a repeated element that
                            !is not the current element ne, it's an adjacent element.
                            candidate_idx=0
                            DO face_node_idx=1,NUMBER_NODE_MATCHES
                              !Dump all the surrounding elements into an array, see if any are repeated
                              node_idx=NODE_MATCHES(face_node_idx)
                              DO surrounding_el_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                                candidate_el=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(surrounding_el_idx)
                                IF(candidate_el/=ne) THEN
                                  candidate_idx=candidate_idx+1
                                  SURROUNDING_ELEMENTS(candidate_idx)=candidate_el
                                ENDIF
                              ENDDO
                            ENDDO !face_node_idx
                            DO idx=1,candidate_idx
                              ne1=SURROUNDING_ELEMENTS(idx)
                              !In a 2D mesh, match 2 nodes, in a 3D mesh, match 3 nodes (line/face)
                              IF(COUNT(SURROUNDING_ELEMENTS(1:candidate_idx)==ne1)>=BASIS%NUMBER_OF_XI) THEN
                                !Found it, just exit
                                CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                                NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                                EXIT
                              ENDIF
                            ENDDO
                          ENDIF
                          IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                        ENDDO !direction_index
                      ENDDO !ni
                    CASE(BASIS_SIMPLEX_TYPE)
                      !Loop over the xi coordinates and calculate the surrounding elements
                      DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                        !Find the other coordinates of the face/line/point
                        FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                        FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                        FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                        !Find nodes in the element on the appropriate face/line/point
                        NULLIFY(NODE_MATCH_LIST)
                        CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                        CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                        CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                        NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                        !Find the nodes to match and add them to the node match list
                        DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                          NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                          DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                            NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                            DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                              NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                & NODE_POSITION_INDEX(3),NODE_POSITION_INDEX(4))
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !nn3
                          ENDDO !nn2
                        ENDDO !nn1
                        CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                        IF(NUMBER_NODE_MATCHES>0) THEN
                          !Find list of elements surrounding those nodes
                          DO node_idx=1,NUMBER_NODE_MATCHES
                            np1=NODE_MATCHES(node_idx)
                            DO nep1=1,DOMAIN_NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                              ne1=DOMAIN_NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                              IF(ne1/=ne) THEN !Don't want the current element
                                ! grab the nodes list for current and this surrouding elements
                                ! current face : NODE_MATCHES
                                ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                                ! if all of current face belongs to the candidate element, we will have found the neighbour
                                CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),DOMAIN_ELEMENTS%ELEMENTS(ne1)% &
                                  & ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                                IF(SUBSET) THEN
                                  CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                            ENDDO !nep1
                          ENDDO !node_idx
                        ENDIF
                        IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                      ENDDO !nic
                    CASE(BASIS_SERENDIPITY_TYPE)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_AUXILLIARY_TYPE)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_B_SPLINE_TP_TYPE)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    CASE DEFAULT
                      LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                        & " is invalid."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)                     
                    END SELECT
                    !Set the surrounding elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
                      BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent elements.",ERR,ERROR,*999)
                    DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      CALL DECOMPOSITION_ADJACENT_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic), &
                        & ERR,ERROR,*999)
                      CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
                      ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS( &
                        DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element adjacent elements.",ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1: &
                        & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)= &
                        ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
                      IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
                    ENDDO !nic
                  ENDDO !ne           
                ELSE
                  CALL FLAG_ERROR("Domain topology elements is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology decomposition is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not allocated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",DECOMPOSITION_ELEMENTS% &
        & TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF

    CALL EXITS("DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 CALL ERRORS("DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE")
    RETURN 1   
  END SUBROUTINE DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Calculates the decomposition element topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: global_element,INSERT_STATUS,local_element
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
        DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(DOMAIN)) THEN
            DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
            IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
              DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
              IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                DOMAIN_MAPPINGS=>DOMAIN%MAPPINGS
                IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
                  DOMAIN_ELEMENTS_MAPPING=>DOMAIN_MAPPINGS%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS_MAPPING)) THEN
                    MESH=>DECOMPOSITION%MESH
                    IF(ASSOCIATED(MESH)) THEN
                      MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
                      IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
                        MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(MESH_ELEMENTS)) THEN
                          !Allocate the element topology arrays
                          ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition elements elements.",ERR,ERROR,*999)
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS
                          CALL TREE_CREATE_START(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                          CALL TREE_INSERT_TYPE_SET(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                          CALL TREE_CREATE_FINISH(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                          DO local_element=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(local_element), &
                              & ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(local_element)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%USER_NUMBER=MESH_ELEMENTS%ELEMENTS(global_element)% &
                              & USER_NUMBER
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%LOCAL_NUMBER=local_element
                            CALL TREE_ITEM_INSERT(DECOMPOSITION_ELEMENTS%ELEMENTS_TREE,DECOMPOSITION_ELEMENTS% &
                              & ELEMENTS(local_element)%USER_NUMBER,local_element,INSERT_STATUS,ERR,ERROR,*999)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%GLOBAL_NUMBER=global_element
                            DECOMPOSITION_ELEMENTS%ELEMENTS(local_element)%BOUNDARY_ELEMENT=MESH_ELEMENTS% &
                              & ELEMENTS(global_element)%BOUNDARY_ELEMENT
                          ENDDO !local_element
                          !Calculate the elements surrounding the elements in the decomposition topology
                          CALL DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
                        ELSE
                          CALL FLAG_ERROR("Mesh elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain mappings elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain topology elements is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Topology decomposition domain is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Topology decomposition is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given decomposition topology. \todo Pass in the decomposition elements pointer.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)) CALL TREE_DESTROY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FLAG_ERROR("Decomposition already has topology elements associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology elements.",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given decomposition. \todo Pass in a pointer to the decomposition topology
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      IF(DECOMPOSITION%CALCULATE_LINES) THEN
        CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      IF(DECOMPOSITION%CALCULATE_FACES) THEN
        CALL DECOMPOSITION_TOPOLOGY_FACES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
      DEALLOCATE(DECOMPOSITION%TOPOLOGY)
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FINALISE")
    RETURN 1
   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given decomposition.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
        CALL FLAG_ERROR("Decomposition already has topology associated.",ERR,ERROR,*999)
      ELSE
        !Allocate decomposition topology
        ALLOCATE(DECOMPOSITION%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Decomposition topology could not be allocated.",ERR,ERROR,*999)
        DECOMPOSITION%TOPOLOGY%DECOMPOSITION=>DECOMPOSITION
        NULLIFY(DECOMPOSITION%TOPOLOGY%ELEMENTS)
        NULLIFY(DECOMPOSITION%TOPOLOGY%LINES)
        NULLIFY(DECOMPOSITION%TOPOLOGY%FACES)
        !Initialise the topology components
        CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        IF(DECOMPOSITION%CALCULATE_LINES) THEN !Default is currently true
          CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF
        IF(DECOMPOSITION%CALCULATE_FACES) THEN !Default is currently false
          CALL DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(LINE%SURROUNDING_ELEMENTS)) DEALLOCATE(LINE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(LINE%ELEMENT_LINES)) DEALLOCATE(LINE%ELEMENT_LINES)
    LINE%ADJACENT_LINES=0
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINE_FINALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a decomposition topology line.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE !<The decomposition line to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    LINE%ADJACENT_LINES=0
    LINE%BOUNDARY_LINE=.FALSE.    

    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the lines in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,element_idx,surrounding_element_idx,basis_local_line_idx, &
      & surrounding_element_basis_local_line_idx,element_local_node_idx,basis_local_line_node_idx,derivative_idx,version_idx, &
      & local_line_idx,surrounding_element_local_line_idx,node_idx,local_node_idx,elem_idx,line_end_node_idx,basis_node_idx, &
      & NODES_IN_LINE(4),NUMBER_OF_LINES,MAX_NUMBER_OF_LINES,NEW_MAX_NUMBER_OF_LINES,LINE_NUMBER,COUNT
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_LINES(:)
    INTEGER(INTG), POINTER :: TEMP_LINES(:,:),NEW_TEMP_LINES(:,:)
    REAL(DP) :: APPROX_DIMENSION
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_LINE_TYPE), POINTER :: DECOMPOSITION_LINE,DECOMPOSITION_LINE2
    TYPE(DECOMPOSITION_LINES_TYPE), POINTER :: DECOMPOSITION_LINES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_LINE_TYPE), POINTER :: DOMAIN_LINE,DOMAIN_LINE2
    TYPE(DOMAIN_LINES_TYPE), POINTER :: DOMAIN_LINES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_LINES)
    NULLIFY(NEW_TEMP_LINES)
    
    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_LINES=>TOPOLOGY%LINES
      IF(ASSOCIATED(DECOMPOSITION_LINES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish line
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Guestimate the number of lines
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      MAX_NUMBER_OF_LINES=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                    CASE(2)
                      APPROX_DIMENSION=SQRT(REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP))
                      !This should give the maximum and will over estimate the number of lines for a "square mesh" by approx 33%
                      MAX_NUMBER_OF_LINES=NINT(3.0_DP*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE(3)
                      !This should give the maximum and will over estimate the number of lines for a "cube mesh" by approx 73%
                      APPROX_DIMENSION=REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP)**(1.0_DP/3.0_DP)
                      MAX_NUMBER_OF_LINES=NINT(11.0_DP*APPROX_DIMENSION*APPROX_DIMENSION*(APPROX_DIMENSION+1),INTG)
                    CASE DEFAULT
                      CALL FLAG_ERROR("Invalid number of dimensions for a topology domain.",ERR,ERROR,*999)
                    END SELECT
                    DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES
                    IF(ASSOCIATED(DOMAIN_LINES)) THEN
                      ALLOCATE(TEMP_LINES(4,MAX_NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary lines array.",ERR,ERROR,*999)
                      ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                      NODES_NUMBER_OF_LINES=0
                      NUMBER_OF_LINES=0
                      TEMP_LINES=0
                      !Loop over the elements in the topology
                      DO element_idx=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_LINES(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element element lines.",ERR,ERROR,*999)
                        !Loop over the local lines of the element
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          !Calculate the topology node numbers that make up the line
                          NODES_IN_LINE=0
                          DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                            NODES_IN_LINE(basis_local_line_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                              & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                          ENDDO !basis_local_line_node_idx
                          !Try and find a previously created line that matches in the adjacent elements
                          FOUND=.FALSE.
                          node_idx=NODES_IN_LINE(1)
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            IF(surrounding_element_idx/=element_idx) THEN
                              IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_LINES)) THEN
                                BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                DO surrounding_element_basis_local_line_idx=1,BASIS2%NUMBER_OF_LOCAL_LINES
                                  local_line_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)% &
                                    & ELEMENT_LINES(surrounding_element_basis_local_line_idx)
                                  IF(ALL(NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx))== &
                                    & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx),local_line_idx))) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                  ENDIF
                                ENDDO !surrounding_element_basis_local_line_idx
                                IF(FOUND) EXIT
                              ENDIF
                            ENDIF
                          ENDDO !elem_idx
                          IF(FOUND) THEN
                            !Line has already been created
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=local_line_idx
                          ELSE
                            !Line has not been created
                            IF(NUMBER_OF_LINES==MAX_NUMBER_OF_LINES) THEN
                              !We are at maximum. Reallocate the LINES array to be 20% bigger and try again.
                              NEW_MAX_NUMBER_OF_LINES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_LINES,DP),INTG)
                              ALLOCATE(NEW_TEMP_LINES(4,NEW_MAX_NUMBER_OF_LINES),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new number of lines.",ERR,ERROR,*999)
                              NEW_TEMP_LINES(:,1:NUMBER_OF_LINES)=TEMP_LINES(:,1:NUMBER_OF_LINES)
                              NEW_TEMP_LINES(:,NUMBER_OF_LINES+1:NEW_MAX_NUMBER_OF_LINES)=0
                              DEALLOCATE(TEMP_LINES)
                              TEMP_LINES=>NEW_TEMP_LINES
                              NULLIFY(NEW_TEMP_LINES)
                              MAX_NUMBER_OF_LINES=NEW_MAX_NUMBER_OF_LINES
                            ENDIF
                            NUMBER_OF_LINES=NUMBER_OF_LINES+1
                            TEMP_LINES(:,NUMBER_OF_LINES)=NODES_IN_LINE
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)=NUMBER_OF_LINES
                            DO basis_local_line_node_idx=1,SIZE(NODES_IN_LINE,1)
                              IF(NODES_IN_LINE(basis_local_line_node_idx)/=0) &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))= &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(basis_local_line_node_idx))+1
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      !Allocate the line arrays and set them from the LINES and NODE_LINES arrays
                      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                        ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node lines array.",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                      ENDDO !node_idx
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition topology lines.",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      ALLOCATE(DOMAIN_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain topology lines.",ERR,ERROR,*999)
                      DOMAIN_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                        CALL DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(DECOMPOSITION_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                        DO basis_local_line_node_idx=1,SIZE(TEMP_LINES,1)
                          IF(TEMP_LINES(basis_local_line_node_idx,local_line_idx)/=0) THEN
                            node_idx=TEMP_LINES(basis_local_line_node_idx,local_line_idx)
                            DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES+1
                            DOMAIN_NODES%NODES(node_idx)%NODE_LINES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES)= &
                              & local_line_idx
                          ENDIF
                        ENDDO !basis_local_line_node_idx  
                      ENDDO !local_line_idx
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          IF(.NOT.ASSOCIATED(DOMAIN_LINE%BASIS)) THEN
                            DECOMPOSITION_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%ELEMENT_NUMBER=element_idx !Needs checking
                            DECOMPOSITION_LINE%XI_DIRECTION=BASIS%LOCAL_LINE_XI_DIRECTION(basis_local_line_idx)
                            DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                            ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line nodes in line.",ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
                            ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                              & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line derivatives in line.",ERR,ERROR,*999)
                            DOMAIN_LINE%NODES_IN_LINE=TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx), &
                              & LINE_NUMBER)
                            DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                              !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                              !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                              version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(2,1, &
                                & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                              IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                  & 1,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                  & 2,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
      !#######################################################VERSIONS END##########################################################
                              ENDIF
                            ENDDO !basis_local_line_node_idx
                          ENDIF
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                      DEALLOCATE(TEMP_LINES)
                      !Calculate adjacent lines and the surrounding elements for each line
                      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                        BASIS=>DOMAIN_LINE%BASIS
                        IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                          DECOMPOSITION_LINE%BOUNDARY_LINE=.TRUE.
                          DOMAIN_LINE%BOUNDARY_LINE=.TRUE.
                        ENDIF
                        !Allocate the elements surrounding the line
                        ALLOCATE(DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line surrounding elements.",ERR,ERROR,*999)
                        ALLOCATE(DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line element lines.",ERR,ERROR,*999)
                        DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
                        DECOMPOSITION_LINE%ADJACENT_LINES=0
                        !Loop over the nodes at each end of the line
                        DO line_end_node_idx=0,1
                          FOUND=.FALSE.
                          node_idx=DOMAIN_LINE%NODES_IN_LINE(line_end_node_idx*(BASIS%NUMBER_OF_NODES-1)+1)
                          !Loop over the elements surrounding the node.
                          DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                            element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                            DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                            DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                            !Loop over the local lines of the element
                            DO basis_local_line_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_LINES
                              surrounding_element_local_line_idx=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                              IF(surrounding_element_local_line_idx/=local_line_idx) THEN
                                DECOMPOSITION_LINE2=>DECOMPOSITION_LINES%LINES(surrounding_element_local_line_idx)
                                DOMAIN_LINE2=>DOMAIN_LINES%LINES(surrounding_element_local_line_idx)
                                IF(DECOMPOSITION_LINE2%XI_DIRECTION==DECOMPOSITION_LINE%XI_DIRECTION) THEN
                                  !Lines run in the same direction.
                                  BASIS2=>DOMAIN_LINE2%BASIS
                                  IF(line_end_node_idx==0) THEN
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES)
                                  ELSE
                                    local_node_idx=DOMAIN_LINE2%NODES_IN_LINE(1)
                                  ENDIF
                                  IF(local_node_idx==node_idx) THEN
                                    !The node at the 'other' end of this line matches the node at the current end of the line.
                                    !Check it is not a coexistant line running the other way
                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
                                      COUNT=0
                                      DO basis_node_idx=1,BASIS%NUMBER_OF_NODES
                                        IF(DOMAIN_LINE2%NODES_IN_LINE(basis_node_idx)== &
                                          & DOMAIN_LINE%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES-basis_node_idx+1)) &
                                          & COUNT=COUNT+1
                                      ENDDO !basis_node_idx
                                      IF(COUNT<BASIS%NUMBER_OF_NODES) THEN
                                        FOUND=.TRUE.
                                        EXIT
                                      ENDIF
                                    ELSE
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDIF
                                ENDIF
                              ENDIF
                            ENDDO !basis_local_line_idx
                            IF(FOUND) EXIT
                          ENDDO !element_idx
                          IF(FOUND) DECOMPOSITION_LINE%ADJACENT_LINES(line_end_node_idx)=surrounding_element_local_line_idx
                        ENDDO !line_end_node_idx
                      ENDDO !local_line_idx
                      !Set the surrounding elements
                      DO element_idx=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(element_idx)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO basis_local_line_idx=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(basis_local_line_idx)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=element_idx
                          DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=basis_local_line_idx
                        ENDDO !basis_local_line_idx
                      ENDDO !element_idx
                    ELSE
                      CALL FLAG_ERROR("Domain topology lines is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain topology elements is not associated.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Topology decomposition domain topology is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Topology decomposition domain is not associated.",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain lines
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES                      
                          IF(ASSOCIATED(DOMAIN_LINES)) THEN
                            ALLOCATE(DOMAIN_LINES%LINES(DECOMPOSITION_LINES%NUMBER_OF_LINES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain lines lines.",ERR,ERROR,*999)
                            DOMAIN_LINES%NUMBER_OF_LINES=DECOMPOSITION_LINES%NUMBER_OF_LINES
                            ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes number of lines array.",ERR,ERROR,*999)
                            NODES_NUMBER_OF_LINES=0
                            !Loop over the lines in the topology
                            DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                              DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                element_idx=DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(1)
                                basis_local_line_idx=DECOMPOSITION_LINE%ELEMENT_LINES(1)
                                CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(local_line_idx),ERR,ERROR,*999)
                                DOMAIN_LINE%NUMBER=local_line_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_LINE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                                DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                                ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes in line.",ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
                                ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(2,DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivatives in line.",ERR,ERROR,*999)
                                DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(basis_local_line_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx, &
                                    & basis_local_line_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain line
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(1,1,basis_local_line_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain line
                                  version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(2,1, &
                                    & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(2,1,basis_local_line_node_idx)=version_idx
                                  IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                      & 1,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(1,2,basis_local_line_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                      & 2,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx), &
                                      & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(basis_local_line_node_idx,basis_local_line_idx))
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(2,2,basis_local_line_node_idx)=version_idx
                                  ENDIF
                                  NODES_NUMBER_OF_LINES(node_idx)=NODES_NUMBER_OF_LINES(node_idx)+1
      !#######################################################VERSIONS END##########################################################
                                ENDDO !basis_local_line_node_idx
                              ELSE
                                CALL FLAG_ERROR("Line is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF
                            ENDDO !local_line_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_LINES(NODES_NUMBER_OF_LINES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node lines.",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_LINES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO local_line_idx=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO basis_local_line_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_LINE%NODES_IN_LINE(basis_local_line_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_LINES=DOMAIN_NODE%NUMBER_OF_NODE_LINES+1
                                DOMAIN_NODE%NODE_LINES(DOMAIN_NODE%NUMBER_OF_NODE_LINES)=local_line_idx
                              ENDDO !basis_local_line_node_idx
                            ENDDO !local_line_idx
                          ELSE
                            CALL FLAG_ERROR("Domain lines is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Domain elements is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain nodes is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FLAG_ERROR("Topology decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Topology decomposition elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology lines is not associated.",ERR,ERROR,*999)

      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology lines:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of lines = ",DECOMPOSITION_LINES%NUMBER_OF_LINES,ERR,ERROR,*999)
      DO local_line_idx=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(local_line_idx)
        DOMAIN_LINE=>DOMAIN_LINES%LINES(local_line_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Line number = ",DECOMPOSITION_LINE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",DECOMPOSITION_LINE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_LINE%ELEMENT_LINES,'("      Element lines        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_LINE%ADJACENT_LINES, &
          & '("      Adjacent lines       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary line = ",DECOMPOSITION_LINE%BOUNDARY_LINE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_LINE=>DOMAIN%TOPOLOGY%LINES%LINES(local_line_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_LINE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_LINE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_LINE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in lines = ",DOMAIN_LINE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_LINE%NODES_IN_LINE, &
            & '("        Nodes in line        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_line_node_idx=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_line_node_idx,ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(1,:,basis_local_line_node_idx),'("            Derivatives in line  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_line_node_idx),4,4, &
              & DOMAIN_LINE%DERIVATIVES_IN_LINE(2,:,basis_local_line_node_idx), &
              & '("            Derivatives Versions in line  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
      !#######################################################VERSIONS END##########################################################
          ENDDO !basis_local_line_node_idx
        ENDDO !component_idx
      ENDDO !local_line_idx
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_LINES)) DEALLOCATE(TEMP_LINES)
    IF(ASSOCIATED(NEW_TEMP_LINES)) DEALLOCATE(NEW_TEMP_LINES)
    IF(ALLOCATED(NODES_NUMBER_OF_LINES)) DEALLOCATE(NODES_NUMBER_OF_LINES)
    CALL ERRORS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINES_CALCULATE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given decomposition topology. \todo Pass in the topology lines
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DECOMPOSITION_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINES_FINALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FLAG_ERROR("Decomposition already has topology lines associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology lines.",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !
 
  !>Finalises a face in the given decomposition topology and deallocates all memory.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    IF(ALLOCATED(FACE%SURROUNDING_ELEMENTS)) DEALLOCATE(FACE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(FACE%ELEMENT_FACES)) DEALLOCATE(FACE%ELEMENT_FACES)
!    FACE%ADJACENT_FACES=0
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACE_FINALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a decomposition topology face.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_FACE_TYPE) :: FACE !<The decomposition face to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    FACE%XI_DIRECTION=0
    FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!    FACE%ADJACENT_FACES=0
    FACE%BOUNDARY_FACE=.FALSE.
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACE_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,ne,surrounding_element_idx,basis_local_face_idx,surrounding_element_basis_local_face_idx, &
      & element_local_node_idx,basis_local_face_node_idx,derivative_idx,version_idx,face_idx,node_idx,elem_idx,NODES_IN_FACE(16), &
      & NUMBER_OF_FACES,MAX_NUMBER_OF_FACES,NEW_MAX_NUMBER_OF_FACES,FACE_NUMBER!,NODE_COUNT,node_idx1,node_idx2,node_idx3,node_idx4,nf2,np2
    INTEGER(INTG), ALLOCATABLE :: NODES_NUMBER_OF_FACES(:)
    INTEGER(INTG), POINTER :: TEMP_FACES(:,:),NEW_TEMP_FACES(:,:)
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS,BASIS2
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: DECOMPOSITION_ELEMENT
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: DECOMPOSITION_FACE!,DECOMPOSITION_FACE2
    TYPE(DECOMPOSITION_FACES_TYPE), POINTER :: DECOMPOSITION_FACES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENT_TYPE), POINTER :: DOMAIN_ELEMENT
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_FACE_TYPE), POINTER :: DOMAIN_FACE!,DOMAIN_FACE2
    TYPE(DOMAIN_FACES_TYPE), POINTER :: DOMAIN_FACES
    TYPE(DOMAIN_NODE_TYPE), POINTER :: DOMAIN_NODE
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
    TYPE(MESH_TYPE), POINTER :: MESH

    NULLIFY(TEMP_FACES)
    NULLIFY(NEW_TEMP_FACES)
    
    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      DECOMPOSITION_FACES=>TOPOLOGY%FACES
      IF(ASSOCIATED(DECOMPOSITION_FACES)) THEN
        DECOMPOSITION_ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(DECOMPOSITION_ELEMENTS)) THEN
          DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
          IF(ASSOCIATED(DECOMPOSITION)) THEN
            !Process the mesh component number (component number the decomposition was calculated from) first to establish face
            !topology then process the other mesh components.
            DOMAIN=>DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(DOMAIN)) THEN
              DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
              IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                IF(ASSOCIATED(DOMAIN_NODES)) THEN
                  DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                  IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                    !Guestimate the number of faces
                    SELECT CASE(DOMAIN%NUMBER_OF_DIMENSIONS)
                    CASE(1)
                      ! Faces not calculated in 1D 
                    CASE(2)
                      ! Faces not calculated in 2D
                    CASE(3)
                      !This should give the maximum and will over estimate the number of faces for a "cube mesh" by approx 33%
                      MAX_NUMBER_OF_FACES=NINT(((REAL(DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,DP)*5.0_DP)+1.0_DP)&
                                                                                                 & *(4.0_DP/3.0_DP),INTG)

                      DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES
                      IF(ASSOCIATED(DOMAIN_FACES)) THEN
                        ALLOCATE(TEMP_FACES(16,MAX_NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary faces array",ERR,ERROR,*999)
                        ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                        NODES_NUMBER_OF_FACES=0
                        NUMBER_OF_FACES=0
                        TEMP_FACES=0
                        !Loop over the elements in the topology and fill temp_faces with node numbers for each element
                        DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_FACES(BASIS%NUMBER_OF_LOCAL_FACES),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element faces of element",ERR,ERROR,*999)
                          !Loop over the local faces of the element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            !Calculate the topology node numbers that make up the face
                            NODES_IN_FACE=0
                            !Check whether face has already been read out
                            DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                              !Read out node numbers of local face from ELEMENT_NODES
                              NODES_IN_FACE(basis_local_face_node_idx)=DOMAIN_ELEMENT%ELEMENT_NODES( &
                                & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                            ENDDO !basis_local_face_node_idx
                            !Try and find a previously created face that matches in the adjacent elements
                            FOUND=.FALSE.
                            node_idx=NODES_IN_FACE(1)
                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                              surrounding_element_idx=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
                              IF(surrounding_element_idx/=ne) THEN
                                IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES)) THEN
                                  BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(surrounding_element_idx)%BASIS
                                  DO surrounding_element_basis_local_face_idx=1,BASIS2%NUMBER_OF_LOCAL_FACES
                                    face_idx=DECOMPOSITION_ELEMENTS%ELEMENTS(surrounding_element_idx)%ELEMENT_FACES( &
                                      & surrounding_element_basis_local_face_idx)
                                    IF(ALL(NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx))== &
                                      & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx),face_idx))) THEN
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !surrounding_element_basis_local_face_idx
                                  IF(FOUND) EXIT
                                ENDIF
                              ENDIF
                            ENDDO !elem_idx
                            IF(FOUND) THEN
                              !Face has already been created
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=face_idx
                            ELSE
                              !Face has not been created
                              IF(NUMBER_OF_FACES==MAX_NUMBER_OF_FACES) THEN
                                !We are at maximum. Reallocate the FACES array to be 20% bigger and try again.
                                NEW_MAX_NUMBER_OF_FACES=NINT(1.20_DP*REAL(MAX_NUMBER_OF_FACES,DP),INTG)
!HERE: Change 16 to a variable and above for NODES_IN_FACE
                                ALLOCATE(NEW_TEMP_FACES(16,NEW_MAX_NUMBER_OF_FACES),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new number of faces",ERR,ERROR,*999)
                                NEW_TEMP_FACES(:,1:NUMBER_OF_FACES)=TEMP_FACES(:,1:NUMBER_OF_FACES)
                                NEW_TEMP_FACES(:,NUMBER_OF_FACES+1:NEW_MAX_NUMBER_OF_FACES)=0
                                DEALLOCATE(TEMP_FACES)
                                TEMP_FACES=>NEW_TEMP_FACES
                                NULLIFY(NEW_TEMP_FACES)
                                MAX_NUMBER_OF_FACES=NEW_MAX_NUMBER_OF_FACES
                              ENDIF
                              NUMBER_OF_FACES=NUMBER_OF_FACES+1
                              TEMP_FACES(:,NUMBER_OF_FACES)=NODES_IN_FACE(:)
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)=NUMBER_OF_FACES
                              DO basis_local_face_node_idx=1,SIZE(NODES_IN_FACE,1)
                                IF(NODES_IN_FACE(basis_local_face_node_idx)/=0) &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))= &
                                  & NODES_NUMBER_OF_FACES(NODES_IN_FACE(basis_local_face_node_idx))+1
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        !Allocate the face arrays and set them from the FACES and NODE_FACES arrays
                        DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                          ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node faces array",ERR,ERROR,*999)
                          DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                        ENDDO !node_idx
                        DEALLOCATE(NODES_NUMBER_OF_FACES)
                        ALLOCATE(DECOMPOSITION_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition topology faces",ERR,ERROR,*999)
                        DECOMPOSITION_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        ALLOCATE(DOMAIN_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain topology faces",ERR,ERROR,*999)
                        DOMAIN_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                          CALL DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(DECOMPOSITION_FACES%FACES(face_idx),ERR,ERROR,*999)
                          CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                          DO basis_local_face_node_idx=1,SIZE(TEMP_FACES,1)
                            IF(TEMP_FACES(basis_local_face_node_idx,face_idx)/=0) THEN
                              node_idx=TEMP_FACES(basis_local_face_node_idx,face_idx)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES+1
                              DOMAIN_NODES%NODES(node_idx)%NODE_FACES(DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES)=face_idx
                            ENDIF
                          ENDDO !basis_local_face_node_idx  
                        ENDDO !face_idx

                        !Set nodes in face and derivatives of nodes in face for domain faces
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          !Loop over local faces of element
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                            DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                            IF(.NOT.ASSOCIATED(DOMAIN_FACE%BASIS)) THEN
                              DECOMPOSITION_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%ELEMENT_NUMBER=ne !! Needs checking
!                              DECOMPOSITION_FACE%ELEMENT_NUMBER=DECOMPOSITION_ELEMENT%NUMBER
!                              DOMAIN_FACE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                              DECOMPOSITION_FACE%XI_DIRECTION=BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx)
                              DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(ABS(DECOMPOSITION_FACE%XI_DIRECTION))%PTR
                              ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                & STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face nodes in face",ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
                              ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face derivatives in face",ERR,ERROR,*999)
                              !Set nodes in face based upon face number
                              DOMAIN_FACE%NODES_IN_FACE=TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx), &
                                & FACE_NUMBER)
                              !Set derivatives of nodes in domain face from derivatives of nodes in element
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(2,1, &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                  derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                    & 1,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx), &
                                    & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(1,2,basis_local_face_node_idx)=derivative_idx
                                  version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                    & 2,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx), &
                                    & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(2,2,basis_local_face_node_idx)=version_idx
      !#######################################################VERSIONS END##########################################################
                                ENDIF
                              ENDDO !basis_local_face_node_idx
                            ENDIF
                          ENDDO !basis_local_face_idx
                        ENDDO !ne

                        DEALLOCATE(TEMP_FACES)
! Note: Adjacency will be left out of faces calculation for the time being
                        !Calculate adjacent faces and the surrounding elements for each face
                        DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                          DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                          DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                          BASIS=>DOMAIN_FACE%BASIS
                          IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS==1) THEN
                            DECOMPOSITION_FACE%BOUNDARY_FACE=.TRUE.
                            DOMAIN_FACE%BOUNDARY_FACE=.TRUE.
                          ENDIF
                          !Allocate the elements surrounding the face
                          ALLOCATE(DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face surrounding elements",ERR,ERROR,*999)

                          ALLOCATE(DECOMPOSITION_FACE%ELEMENT_FACES(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                            & STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face element faces",ERR,ERROR,*999)
!                          DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=0
!                          DECOMPOSITION_FACE%ADJACENT_FACES=0

                           !Loop over the nodes at each end of the face
!                          DO node_idx1=0,1
!                           DO node_idx2=0,1
!                            FOUND=.FALSE.
!                            node_idx=DOMAIN_FACE%NODES_IN_FACE((node_idx2*BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION*(BASIS%NUMBER_OF_FACES-1))&
!                                                                             &+(node_idx1*(BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
                             !Loop over the elements surrounding the node.
!                            DO elem_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
!                              ne=DOMAIN_NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
!                              DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
!                              DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                               !Loop over the local faces of the element
!                              DO basis_local_face_idx=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_FACES
!                                nf2=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
!                                IF(nf2/=face_idx) THEN
!                                  DECOMPOSITION_FACE2=>DECOMPOSITION_FACES%FACES(nf2)
!                                  DOMAIN_FACE2=>DOMAIN_FACES%FACES(nf2)
                                   !Check whether XI of face have same direction
!                                  IF ((OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),2,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),2,1)).OR.&
!                                     &(OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),3,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(basis_local_face_idx),3,1))) THEN
                                     !Loop over nodes in face of surrounding element
!                                    BASIS2=>DOMAIN_FACE2%BASIS
!                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
!                                      NODE_COUNT=0
!                                      DO node_idx3=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                        DO node_idx4=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                          np2=DOMAIN_FACE2%NODES_IN_FACE((node_idx4*(BASIS2%NUMBER_OF_FACES-1))&
!                                                                      &+(node_idx3*(BASIS2%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
!                                          IF(np2==node_idx) NODE_COUNT=NODE_COUNT+1
!                                        ENDDO !node_idx4
!                                      ENDDO !node_idx3
!                                      IF(NODE_COUNT<BASIS%NUMBER_OF_NODES) THEN
!                                        FOUND=.TRUE.
!                                        EXIT
!                                      ENDIF
!                                    ENDIF
!                                  ENDIF
!                                ENDIF
!                              ENDDO !basis_local_face_idx
!                                IF(FOUND) EXIT
!                            ENDDO !elem_idx
!                            IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx2)=nf2
!                           ENDDO !node_idx2
!                           IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx1)=nf2
!                          ENDDO !node_idx1
                        ENDDO !face_idx

                        !Set the surrounding elements
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          DO basis_local_face_idx=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(basis_local_face_idx)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DO face_idx=1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS
                              DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(face_idx)=ne
                              DECOMPOSITION_FACE%ELEMENT_FACES(face_idx)=basis_local_face_idx
                            ENDDO
                          ENDDO !basis_local_face_idx
                        ENDDO !ne
                      ELSE
                        CALL FLAG_ERROR("Domain topology faces is not associated",ERR,ERROR,*999)
                      ENDIF
                    CASE DEFAULT
                      CALL FLAG_ERROR("Invalid number of dimensions for a topology domain",ERR,ERROR,*999)
                    END SELECT
                 ELSE
                    CALL FLAG_ERROR("Domain topology elements is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain topology nodes is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Topology decomposition domain topology is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Topology decomposition domain is not associated",ERR,ERROR,*999)
            ENDIF
            !Now loop over the other mesh components in the decomposition and calculate the domain faces
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=>DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_FACES=>DOMAIN_TOPOLOGY%FACES                      
                          IF(ASSOCIATED(DOMAIN_FACES)) THEN
                            ALLOCATE(DOMAIN_FACES%FACES(DECOMPOSITION_FACES%NUMBER_OF_FACES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain faces faces",ERR,ERROR,*999)
                            DOMAIN_FACES%NUMBER_OF_FACES=DECOMPOSITION_FACES%NUMBER_OF_FACES
                            ALLOCATE(NODES_NUMBER_OF_FACES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes number of faces array",ERR,ERROR,*999)
                            NODES_NUMBER_OF_FACES=0
                            !Loop over the faces in the topology
                            DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                              DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                ne=DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(1)
                                basis_local_face_idx=DECOMPOSITION_FACE%ELEMENT_FACES(1)
                                CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(face_idx),ERR,ERROR,*999)
                                DOMAIN_FACE%NUMBER=face_idx
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(ABS(DECOMPOSITION_FACE%XI_DIRECTION))%PTR
                                ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)), &
                                  & STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes in face",ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
                                ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(2,DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivatives in face",ERR,ERROR,*999)
                                !Set derivatives of nodes in domain face from derivatives of nodes in element
                                DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(basis_local_face_idx)
                                  element_local_node_idx=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx, &
                                    & basis_local_face_idx)
                                  node_idx=DOMAIN_ELEMENT%ELEMENT_NODES(element_local_node_idx)
                                  DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)=node_idx
                                  !Set derivative number of u (NO_GLOBAL_DERIV) for the domain face
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(1,1,basis_local_face_node_idx)=NO_GLOBAL_DERIV
                                  !Set version number of u (NO_GLOBAL_DERIV) for the domain face
                                  version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(2,1, &
                                    & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(2,1,basis_local_face_node_idx)=version_idx
                                  IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    derivative_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                      & 1,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx), &
                                      & element_local_node_idx)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(1,2,basis_local_face_node_idx)=derivative_idx
                                    version_idx=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES( &
                                      & 2,BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx), &
                                      & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(basis_local_face_node_idx,basis_local_face_idx))
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(2,2,basis_local_face_node_idx)=version_idx
                                  ENDIF
                                  NODES_NUMBER_OF_FACES(node_idx)=NODES_NUMBER_OF_FACES(node_idx)+1
      !#######################################################VERSIONS END##########################################################
                                ENDDO !basis_local_face_node_idx
                              ELSE
                                CALL FLAG_ERROR("Face is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF                              
                            ENDDO !face_idx
                            DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(node_idx)%NODE_FACES(NODES_NUMBER_OF_FACES(node_idx)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node faces",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_NODE_FACES=0
                            ENDDO !node_idx
                            DEALLOCATE(NODES_NUMBER_OF_FACES)
                            DO face_idx=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
                              BASIS=>DOMAIN_FACE%BASIS
                              DO basis_local_face_node_idx=1,BASIS%NUMBER_OF_NODES
                                node_idx=DOMAIN_FACE%NODES_IN_FACE(basis_local_face_node_idx)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(node_idx)
                                DOMAIN_NODE%NUMBER_OF_NODE_FACES=DOMAIN_NODE%NUMBER_OF_NODE_FACES+1
                                !Set the face numbers a node is on
                                DOMAIN_NODE%NODE_FACES(DOMAIN_NODE%NUMBER_OF_NODE_FACES)=face_idx
                              ENDDO !basis_local_face_node_idx
                            ENDDO !face_idx
                          ELSE
                            CALL FLAG_ERROR("Domain faces is not associated",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FLAG_ERROR("Domain elements is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Domain nodes is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Domain topology is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
                  ENDIF
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
            ENDIF                        
          ELSE
            CALL FLAG_ERROR("Topology decomposition is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Topology decomposition elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology faces is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Decomposition topology faces:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of mesh components = ",MESH%NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of faces = ",DECOMPOSITION_FACES%NUMBER_OF_FACES,ERR,ERROR,*999)
      DO face_idx=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
        DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(face_idx)
        DOMAIN_FACE=>DOMAIN_FACES%FACES(face_idx)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Face number = ",DECOMPOSITION_FACE%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction (Normal to Face) = &
                                                                         &",DECOMPOSITION_FACE%XI_DIRECTION,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of surrounding elements = ", &
          & DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%SURROUNDING_ELEMENTS,'("      Surrounding elements :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS,4,4, &
          & DECOMPOSITION_FACE%ELEMENT_FACES,'("      Element faces        :",4(X,I8))','(28X,4(X,I8))',ERR,ERROR,*999)
!        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2,DECOMPOSITION_FACE%ADJACENT_FACES, &
!          & '("      Adjacent faces       :",2(X,I8))','(28X,2(X,I8))',ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary face = ",DECOMPOSITION_FACE%BOUNDARY_FACE,ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Mesh component : ",component_idx,ERR,ERROR,*999)
          DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
          DOMAIN_FACE=>DOMAIN%TOPOLOGY%FACES%FACES(face_idx)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis user number = ",DOMAIN_FACE%BASIS%USER_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis family number = ",DOMAIN_FACE%BASIS%FAMILY_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation type = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_TYPE(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Basis interpolation order = ",DOMAIN_FACE%BASIS% &
            & INTERPOLATION_ORDER(1),ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of nodes in faces = ",DOMAIN_FACE%BASIS%NUMBER_OF_NODES, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES,4,4,DOMAIN_FACE%NODES_IN_FACE, &
            & '("        Nodes in face        :",4(X,I8))','(30X,4(X,I8))',ERR,ERROR,*999)
          DO basis_local_face_node_idx=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",basis_local_face_node_idx,ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
            !/TODO::Loop over local_derivative index so this output makes more sense !<DERIVATIVES_IN_LINE(i,local_derivative_idx,local_node_idx)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(1,:,basis_local_face_node_idx),'("            Derivatives in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
              & DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(basis_local_face_node_idx),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(2,:,basis_local_face_node_idx),'("            Derivatives Versions in face  :",4(X,I8))', &
              & '(34X,4(X,I8))',ERR,ERROR,*999)
      !#######################################################VERSIONS END##########################################################
          ENDDO !basis_local_face_node_idx
        ENDDO !component_idx
      ENDDO !face_idx
    ENDIF

    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(TEMP_FACES)) DEALLOCATE(TEMP_FACES)
    IF(ASSOCIATED(NEW_TEMP_FACES)) DEALLOCATE(NEW_TEMP_FACES)
    IF(ALLOCATED(NODES_NUMBER_OF_FACES)) DEALLOCATE(NODES_NUMBER_OF_FACES)
    CALL ERRORS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACES_CALCULATE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_CALCULATE

  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DECOMPOSITION_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACES_FINALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a decomposition topology.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FLAG_ERROR("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_FACES_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_FACES_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the decomposition type for a decomposition. \see OPENCMISS::CMISSDecompositionTypeGet
  SUBROUTINE DECOMPOSITION_TYPE_GET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to get the type for
    INTEGER(INTG), INTENT(OUT) :: TYPE !<On return, the decomposition type for the specified decomposition \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        TYPE=DECOMPOSITION%DECOMPOSITION_TYPE
      ELSE
        CALL FLAG_ERROR("Decomposition has not finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TYPE_GET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TYPE_GET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TYPE_GET")
    RETURN
  END SUBROUTINE DECOMPOSITION_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the decomposition type for a decomposition.  \see OPENCMISS::CMISSDecompositionTypeSet
  SUBROUTINE DECOMPOSITION_TYPE_SET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the type for
    INTEGER(INTG), INTENT(IN) :: TYPE !<The decomposition type to set \see MESH_ROUTINES_DecompositionTypes,MESH_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          !heye: three types for decomposition--decompostion_all_type means no decomposition
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_CALCULATED_TYPE
        CASE(DECOMPOSITION_USER_DEFINED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_USER_DEFINED_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is not valid."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TYPE_SET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TYPE_SET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TYPE_SET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TYPE_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether lines should be calculated in the the decomposition. \see OPENCMISS::CMISSDecompositionCalculateLinesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET(DECOMPOSITION,CALCULATE_LINES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_LINES_FLAG !<The boolean flag to determine whether the lines should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_LINES=CALCULATE_LINES_FLAG
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_CALCULATE_LINES_SET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_CALCULATE_LINES_SET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_CALCULATE_LINES_SET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_LINES_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes whether faces should be calculated in the the decomposition. \see OPENCMISS::CMISSDecompositionCalculateFacesSet
  SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET(DECOMPOSITION,CALCULATE_FACES_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition
    LOGICAL, INTENT(IN) :: CALCULATE_FACES_FLAG !<The boolean flag to determine whether the faces should be calculated or not
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished.",ERR,ERROR,*999)
      ELSE
        DECOMPOSITION%CALCULATE_FACES=CALCULATE_FACES_FLAG
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_CALCULATE_FACES_SET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_CALCULATE_FACES_SET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_CALCULATE_FACES_SET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_CALCULATE_FACES_SET
  
  !
  !================================================================================================================================
  !

  !>Finds and returns in DECOMPOSITION a pointer to the decomposition identified by USER_NUMBER in the given MESH. If no decomposition with that USER_NUMBER exists DECOMPOSITION is left nullified.
  SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the decomposition to find
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh containing the decomposition to find
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<On return a pointer to the decomposition with the specified user number. If no decomposition with that user number exists then DECOMPOSITION is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: decomposition_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_USER_NUMBER_FIND",ERR,ERROR,*999)

    NULLIFY(DECOMPOSITION)
    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        decomposition_idx=1
        DO WHILE(decomposition_idx<=MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS.AND..NOT.ASSOCIATED(DECOMPOSITION))
          IF(MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            DECOMPOSITION=>MESH%DECOMPOSITIONS%DECOMPOSITIONS(decomposition_idx)%PTR
          ELSE
            decomposition_idx=decomposition_idx+1
          ENDIF
        ENDDO
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITION_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("DECOMPOSITION_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  !>Finalises the domain decompositions for a given mesh. \todo Pass in a pointer to the decompositions.
  SUBROUTINE DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise the decomposition for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        DO WHILE(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>0)
          CALL DECOMPOSITION_DESTROY(MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR,ERR,ERROR,*999)
        ENDDO !no_decomposition
       DEALLOCATE(MESH%DECOMPOSITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITIONS_FINALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITIONS_FINALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITIONS_FINALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain decompositions for a given mesh.
  SUBROUTINE DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise the decompositions for

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        CALL FLAG_ERROR("Mesh already has decompositions associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(MESH%DECOMPOSITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Mesh decompositions could not be allocated.",ERR,ERROR,*999)
        MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        NULLIFY(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
        MESH%DECOMPOSITIONS%MESH=>MESH
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DECOMPOSITIONS_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITIONS_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITIONS_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITIONS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the domain for a given decomposition and deallocates all memory. \todo Pass in a pointer to the domain
  SUBROUTINE DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to finalise the domain for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    
    CALL ENTERS("DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS
            IF(ALLOCATED(DECOMPOSITION%DOMAIN(component_idx)%PTR%NODE_DOMAIN))  &
              & DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR%NODE_DOMAIN)
            CALL DOMAIN_MAPPINGS_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)        
            CALL DOMAIN_TOPOLOGY_FINALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            DEALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR)
          ENDDO !component_idx
          DEALLOCATE(DECOMPOSITION%DOMAIN)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_FINALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the domain for a given decomposition.
  SUBROUTINE DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to initialise the domain for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    CALL ENTERS("DOMAIN_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          CALL FLAG_ERROR("Decomposition already has a domain associated.",ERR,ERROR,*999)
        ELSE
          ALLOCATE(DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Decomposition domain could not be allocated.",ERR,ERROR,*999)
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS
            ALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Decomposition domain component could not be allocated.",ERR,ERROR,*999)
            DECOMPOSITION%DOMAIN(component_idx)%PTR%DECOMPOSITION=>DECOMPOSITION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH=>DECOMPOSITION%MESH
            DECOMPOSITION%DOMAIN(component_idx)%PTR%MESH_COMPONENT_NUMBER=component_idx
            DECOMPOSITION%DOMAIN(component_idx)%PTR%REGION=>DECOMPOSITION%MESH%REGION
            DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_DIMENSIONS=DECOMPOSITION%MESH%NUMBER_OF_DIMENSIONS
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_ELEMENTS=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_FACES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_LINES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_NODES=0
            !DECOMPOSITION%DOMAIN(component_idx)%PTR%NUMBER_OF_MESH_DOFS=0
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%MAPPINGS)
            NULLIFY(DECOMPOSITION%DOMAIN(component_idx)%PTR%TOPOLOGY)
            CALL DOMAIN_MAPPINGS_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
            CALL DOMAIN_TOPOLOGY_INITIALISE(DECOMPOSITION%DOMAIN(component_idx)%PTR,ERR,ERROR,*999)
          ENDDO !component_idx
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_INITIALISE

  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs mapping in the given domain mappings. \todo Pass in the domain mappings dofs
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%DOFS,ERR,ERROR,*999)
        DEALLOCATE(DOMAIN_MAPPINGS%DOFS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_DOFS_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_DOFS_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises the dofs mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL FLAG_ERROR("Domain dofs mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings dofs.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%DOFS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_DOFS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the local/global element mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the element mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,adjacent_element,domain_no,domain_idx,ne,nn,np,NUMBER_OF_DOMAINS, &
      & NUMBER_OF_ADJACENT_ELEMENTS,my_computational_node_number,component_idx
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_ELEMENTS(:),DOMAINS(:),LOCAL_ELEMENT_NUMBERS(:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
          ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
          IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
            DECOMPOSITION=>DOMAIN%DECOMPOSITION
            IF(ASSOCIATED(DOMAIN%MESH)) THEN
              MESH=>DOMAIN%MESH
              component_idx=DOMAIN%MESH_COMPONENT_NUMBER
              my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
              IF(ERR/=0) GOTO 999        
              
              !Calculate the local and global numbers and set up the mappings
              ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element mapping global to local map.",ERR,ERROR,*999)
              ELEMENTS_MAPPING%NUMBER_OF_GLOBAL=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS
              !Loop over the global elements and calculate local numbers
              ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local element numbers.",ERR,ERROR,*999)
              LOCAL_ELEMENT_NUMBERS=0
              ALLOCATE(ADJACENT_ELEMENTS_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent elements list.",ERR,ERROR,*999)
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                NULLIFY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR)
                CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,MAX(INT(MESH%NUMBER_OF_ELEMENTS/2),1), &
                  & ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
              ENDDO !domain_idx
            
              DO ne=1,MESH%NUMBER_OF_ELEMENTS
                !Calculate the local numbers
                domain_no=DECOMPOSITION%ELEMENT_DOMAIN(ne)
                LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                !Calculate the adjacent elements to the computational domains and the adjacent domain numbers themselves
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                NULLIFY(ADJACENT_DOMAINS_LIST)
                CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                DO nn=1,BASIS%NUMBER_OF_NODES
                  np=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                  DO no_adjacent_element=1,MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                    adjacent_element=MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%SURROUNDING_ELEMENTS(no_adjacent_element)
                    IF(DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element)/=domain_no) THEN
                      CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(domain_no)%PTR,adjacent_element,ERR,ERROR,*999)
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element),ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_adjacent_element
                ENDDO !nn
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                DEALLOCATE(DOMAINS)
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne),ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element global to local map local number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element global to local map domain number.",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element global to local map local type.",ERR,ERROR,*999)
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS=1
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(1)=LOCAL_ELEMENT_NUMBERS(domain_no)
                ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(1)=DECOMPOSITION%ELEMENT_DOMAIN(ne) 
                IF(NUMBER_OF_DOMAINS==1) THEN
                  !Element is an internal element
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                ELSE
                  !Element is on the boundary of computational domains
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                ENDIF
              ENDDO !ne
              
              !Compute ghost element mappings
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,NUMBER_OF_ADJACENT_ELEMENTS, &
                  & ADJACENT_ELEMENTS,ERR,ERROR,*999)
                DO no_adjacent_element=1,NUMBER_OF_ADJACENT_ELEMENTS
                  adjacent_element=ADJACENT_ELEMENTS(no_adjacent_element)
                  LOCAL_ELEMENT_NUMBERS(domain_idx)=LOCAL_ELEMENT_NUMBERS(domain_idx)+1
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS= &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS+1
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=LOCAL_ELEMENT_NUMBERS(domain_idx)
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%DOMAIN_NUMBER( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)=domain_idx
                  ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%LOCAL_TYPE( &
                    & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS)= &
                    & DOMAIN_LOCAL_GHOST
                ENDDO !no_adjacent_element
                IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
              ENDDO !domain_idx
              
              DEALLOCATE(ADJACENT_ELEMENTS_LIST)
              DEALLOCATE(LOCAL_ELEMENT_NUMBERS)
              
              !Calculate element local to global maps from global to local map
              CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(ELEMENTS_MAPPING,ERR,ERROR,*999)
                            
            ELSE
              CALL FLAG_ERROR("Domain mesh is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Domain decomposition is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain mappings elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Element mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ne=1,MESH%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER, &
          & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)% &
          & NUMBER_OF_DOMAINS,8,8,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ne
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ne=1,ELEMENTS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global element = ", &
          & ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ne),ERR,ERROR,*999)
      ENDDO !ne
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%INTERNAL_START:ELEMENTS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%BOUNDARY_START:ELEMENTS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
          & ELEMENTs_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & ELEMENTS_MAPPING%DOMAIN_LIST(ELEMENTS_MAPPING%GHOST_START:ELEMENTS_MAPPING%GHOST_FINISH), &
          & '("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS_PTR( &
        & ELEMENTS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,ELEMENTS_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,ELEMENTS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,ELEMENTS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)    
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the mappings in the given domain. \todo pass in the domain mappings
  SUBROUTINE DOMAIN_MAPPINGS_FINALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%MAPPINGS)
    ELSE
      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the element mapping in the given domain mapping. \todo pass in the domain mappings elements
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to finalise the elements for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
        DEALLOCATE(DOMAIN_MAPPINGS%ELEMENTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element mapping in the given domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL FLAG_ERROR("Domain elements mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings elements.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%ELEMENTS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the mappings for a domain decomposition. \todo finalise on error.
  SUBROUTINE DOMAIN_MAPPINGS_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the mappings for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        CALL FLAG_ERROR("Domain already has mappings associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN%MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings.",ERR,ERROR,*999)
        DOMAIN%MAPPINGS%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%MAPPINGS%ELEMENTS)
        NULLIFY(DOMAIN%MAPPINGS%NODES)
        NULLIFY(DOMAIN%MAPPINGS%DOFS)
        !Calculate the node and element mappings
        CALL DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the local/global node and dof mappings for a domain decomposition.
  SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to calculate the node dofs for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,no_adjacent_element,no_computational_node,no_ghost_node,adjacent_element,ghost_node, &
      & NUMBER_OF_NODES_PER_DOMAIN,domain_idx,domain_idx2,domain_no,node_idx,derivative_idx,version_idx,ny,NUMBER_OF_DOMAINS, &
      & MAX_NUMBER_DOMAINS,NUMBER_OF_GHOST_NODES,my_computational_node_number,number_computational_nodes,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:),NODE_COUNT(:),NUMBER_INTERNAL_NODES(:), &
      & NUMBER_BOUNDARY_NODES(:)
    INTEGER(INTG), ALLOCATABLE :: DOMAINS(:),ALL_DOMAINS(:),GHOST_NODES(:)
    LOGICAL :: BOUNDARY_DOMAIN
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST,ALL_ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_NODES_LIST(:)
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: NODES_MAPPING
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOFS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    CALL ENTERS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) THEN
          NODES_MAPPING=>DOMAIN%MAPPINGS%NODES
          IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) THEN
            DOFS_MAPPING=>DOMAIN%MAPPINGS%DOFS
            IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) THEN
              ELEMENTS_MAPPING=>DOMAIN%MAPPINGS%ELEMENTS
              IF(ASSOCIATED(DOMAIN%DECOMPOSITION)) THEN
                DECOMPOSITION=>DOMAIN%DECOMPOSITION
                IF(ASSOCIATED(DOMAIN%MESH)) THEN
                  MESH=>DOMAIN%MESH
                  component_idx=DOMAIN%MESH_COMPONENT_NUMBER
                  MESH_TOPOLOGY=>MESH%TOPOLOGY(component_idx)%PTR
                  
                  number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  
                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node mapping global to local map.",ERR,ERROR,*999)
                  NODES_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                  ALLOCATE(LOCAL_NODE_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local node numbers.",ERR,ERROR,*999)
                  LOCAL_NODE_NUMBERS=0
                  ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%DOFS%NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dofs mapping global to local map.",ERR,ERROR,*999)
                  DOFS_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%DOFS%NUMBER_OF_DOFS
                  ALLOCATE(LOCAL_DOF_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local dof numbers.",ERR,ERROR,*999)
                  LOCAL_DOF_NUMBERS=0
                  ALLOCATE(GHOST_NODES_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ghost nodes list.",ERR,ERROR,*999)
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    NULLIFY(GHOST_NODES_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(GHOST_NODES_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(GHOST_NODES_LIST(domain_idx)%PTR,INT(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES/2), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx
                  ALLOCATE(NUMBER_INTERNAL_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of internal nodes.",ERR,ERROR,*999)
                  NUMBER_INTERNAL_NODES=0
                  ALLOCATE(NUMBER_BOUNDARY_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of boundary nodes.",ERR,ERROR,*999)
                  NUMBER_BOUNDARY_NODES=0

                  !For the first pass just determine the internal and boundary nodes
                  DO node_idx=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                    NULLIFY(ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    NULLIFY(ALL_ADJACENT_DOMAINS_LIST)
                    CALL LIST_CREATE_START(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ALL_ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ALL_ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    DO no_adjacent_element=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
                      adjacent_element=MESH_TOPOLOGY%NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(no_adjacent_element)
                      domain_no=DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element)
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                      DO domain_idx=1,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)%NUMBER_OF_DOMAINS
                        CALL LIST_ITEM_ADD(ALL_ADJACENT_DOMAINS_LIST,ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(adjacent_element)% &
                          & DOMAIN_NUMBER(domain_idx),ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDDO !no_adjacent_element
                    CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                    CALL LIST_REMOVE_DUPLICATES(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(ALL_ADJACENT_DOMAINS_LIST,MAX_NUMBER_DOMAINS,ALL_DOMAINS,ERR,ERROR,*999)
                    IF(NUMBER_OF_DOMAINS/=MAX_NUMBER_DOMAINS) THEN !Ghost node
                      DO domain_idx=1,MAX_NUMBER_DOMAINS
                        domain_no=ALL_DOMAINS(domain_idx)
                        BOUNDARY_DOMAIN=.FALSE.
                        DO domain_idx2=1,NUMBER_OF_DOMAINS
                          IF(domain_no==DOMAINS(domain_idx2)) THEN
                            BOUNDARY_DOMAIN=.TRUE.
                            EXIT
                          ENDIF
                        ENDDO !domain_idx2
                        IF(.NOT.BOUNDARY_DOMAIN) CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDIF
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map local type.",ERR,ERROR,*999)
                    DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                      ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)
!                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
!                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
!                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
!                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
!                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
!                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                      DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                        ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)




                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
                        ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                      ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                    ENDDO !derivative_idx
                    IF(NUMBER_OF_DOMAINS==1) THEN
                      !Node is an internal node
                      domain_no=DOMAINS(1)
                      NUMBER_INTERNAL_NODES(domain_no)=NUMBER_INTERNAL_NODES(domain_no)+1
                      !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=1
                      !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(DOMAINS(1))
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=-1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)=DOMAINS(1) 
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                        ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)
!                        !LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
!                        !DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=-1
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=-1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                        ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                      ENDDO !derivative_idx
                    ELSE
                      !Node is on the boundary of computational domains
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                        ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                        ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                      ENDDO !derivative_idx
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=DOMAINS(domain_idx)
                        !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                        !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=LOCAL_NODE_NUMBERS(domain_no)
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(domain_idx)=-1
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(domain_idx)=domain_no
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                        DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)
!                          !LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
!                          !DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=LOCAL_DOF_NUMBERS(domain_no)
!                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=-1
!                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)=domain_no
!                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                          DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                            ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=-1
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)=domain_no
                            DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                          ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                        ENDDO !derivative_idx
                      ENDDO !domain_idx
                    ENDIF
                    DEALLOCATE(DOMAINS) 
                    DEALLOCATE(ALL_DOMAINS)
                  ENDDO !node_idx

                  !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
                  NUMBER_OF_NODES_PER_DOMAIN=FLOOR(REAL(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES,DP)/ &
                    & REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP))
                  ALLOCATE(DOMAIN%NODE_DOMAIN(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node domain",ERR,ERROR,*999)
                  DOMAIN%NODE_DOMAIN=-1
                  DO node_idx=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                    IF(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS==1) THEN !Internal node
                      domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)
                      DOMAIN%NODE_DOMAIN(node_idx)=domain_no
                      LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                        ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)
!                        LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                          LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                        ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                      ENDDO !derivative_idx
                    ELSE !Boundary node
                      NUMBER_OF_DOMAINS=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(domain_idx)
                        IF(DOMAIN%NODE_DOMAIN(node_idx)<0) THEN
                          IF((NUMBER_INTERNAL_NODES(domain_no)+NUMBER_BOUNDARY_NODES(domain_no)<NUMBER_OF_NODES_PER_DOMAIN).OR. &
                            & (domain_idx==NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS)) THEN
                            !Allocate the node to this domain
                            DOMAIN%NODE_DOMAIN(node_idx)=domain_no
                            NUMBER_BOUNDARY_NODES(domain_no)=NUMBER_BOUNDARY_NODES(domain_no)+1
                            LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                            !Reset the boundary information to be in the first domain index. The remaining domain indicies will
                            !be overwritten when the ghost nodes are calculated below. 
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS=1
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER(1)=domain_no
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                            DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                              ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)
!                              LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
!                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
!                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
!                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
!                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                              DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                                ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                                LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                                DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                              ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                            ENDDO !derivative_idx
                          ELSE
                            !The node as already been assigned to a domain so it must be a ghost node in this domain
                            CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          !The node as already been assigned to a domain so it must be a ghost node in this domain
                          CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,node_idx,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !domain_idx
                    ENDIF
                  ENDDO !node_idx
                  DEALLOCATE(NUMBER_INTERNAL_NODES)
                  
                  !Calculate ghost node and dof mappings
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    CALL LIST_REMOVE_DUPLICATES(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DETACH_AND_DESTROY(GHOST_NODES_LIST(domain_idx)%PTR,NUMBER_OF_GHOST_NODES,GHOST_NODES,ERR,ERROR,*999)
                    DO no_ghost_node=1,NUMBER_OF_GHOST_NODES
                      ghost_node=GHOST_NODES(no_ghost_node)
                      LOCAL_NODE_NUMBERS(domain_idx)=LOCAL_NODE_NUMBERS(domain_idx)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS= &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%LOCAL_NUMBER( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)= &
                        & LOCAL_NODE_NUMBERS(domain_idx)
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%DOMAIN_NUMBER( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)=domain_idx
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%LOCAL_TYPE( &
                        & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(ghost_node)%NUMBER_OF_DOMAINS)= &
                        & DOMAIN_LOCAL_GHOST
                      DO derivative_idx=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                         !\todo DELETE THIS
!                        ny=MESH_TOPOLOGY%NODES%NODES(ghost_node)%DOF_INDEX(derivative_idx)
!                        LOCAL_DOF_NUMBERS(domain_idx)=LOCAL_DOF_NUMBERS(domain_idx)+1
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS= &
!                          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS+1
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER( &
!                          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
!                          & LOCAL_DOF_NUMBERS(domain_idx)
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER( &
!                          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)=domain_idx
!                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE( &
!                          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
!                          & DOMAIN_LOCAL_GHOST
                        DO version_idx=1,MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                          ny=MESH_TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)
                          LOCAL_DOF_NUMBERS(domain_idx)=LOCAL_DOF_NUMBERS(domain_idx)+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS= &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS+1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & LOCAL_DOF_NUMBERS(domain_idx)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)=domain_idx
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE( &
                            & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS)= &
                            & DOMAIN_LOCAL_GHOST
                        ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                      ENDDO !derivative_idx
                    ENDDO !no_ghost_node
                    DEALLOCATE(GHOST_NODES)
                  ENDDO !domain_idx
                  
                  !Check decomposition and check that each domain has a node in it.
                  ALLOCATE(NODE_COUNT(0:number_computational_nodes-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node count.",ERR,ERROR,*999)
                  NODE_COUNT=0
                  DO node_idx=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                    no_computational_node=DOMAIN%NODE_DOMAIN(node_idx)
                    IF(no_computational_node>=0.AND.no_computational_node<number_computational_nodes) THEN
                      NODE_COUNT(no_computational_node)=NODE_COUNT(no_computational_node)+1
                    ELSE
                      LOCAL_ERROR="The computational node number of "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
                        & " for node number "//TRIM(NUMBER_TO_VSTRING(node_idx,"*",ERR,ERROR))// &
                        & " is invalid. The computational node number must be between 0 and "// &
                        & TRIM(NUMBER_TO_VSTRING(number_computational_nodes-1,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !node_idx
                  DO no_computational_node=0,number_computational_nodes-1
                    IF(NODE_COUNT(no_computational_node)==0) THEN
                      LOCAL_ERROR="Invalid decomposition. There are no nodes in computational node "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_computational_node
                  DEALLOCATE(NODE_COUNT)
          
                  DEALLOCATE(GHOST_NODES_LIST)
                  DEALLOCATE(LOCAL_NODE_NUMBERS)
                  
                  !Calculate node and dof local to global maps from global to local map
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(NODES_MAPPING,ERR,ERROR,*999)
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOFS_MAPPING,ERR,ERROR,*999)
                  
                ELSE
                  CALL FLAG_ERROR("Domain mesh is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain decomposition is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Domain mappings elements is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Domain mappings dofs is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain mappings nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain mappings is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node decomposition :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",DOMAIN%NODE_DOMAIN(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO node_idx=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
            & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(node_idx)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !node_idx     
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO node_idx=1,NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(node_idx),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal nodes = ", &
          & NODES_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%INTERNAL_START:NODES_MAPPING%INTERNAL_FINISH), &
          & '("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & NODES_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%BOUNDARY_START:NODES_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & NODES_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_GHOST,8,8, &
          & NODES_MAPPING%DOMAIN_LIST(NODES_MAPPING%GHOST_START:NODES_MAPPING%GHOST_FINISH), &
          & '("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & NODES_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS_PTR( &
        & NODES_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,NODES_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,NODES_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,NODES_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Dofs mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO ny=1,MESH_TOPOLOGY%DOFS%NUMBER_OF_DOFS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
            & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)% &
          & NUMBER_OF_DOMAINS,8,8,DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !ny   
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO ny=1,DOFS_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local dof = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global dof = ", &
          & DOFS_MAPPING%LOCAL_TO_GLOBAL_MAP(ny),ERR,ERROR,*999)
      ENDDO !node_idx
      IF(DIAGNOSTICS2) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Internal dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_INTERNAL,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%INTERNAL_START:DOFS_MAPPING%INTERNAL_FINISH), &
          & '("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%BOUNDARY_START:DOFS_MAPPING%BOUNDARY_FINISH), &
          & '("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & DOFS_MAPPING%DOMAIN_LIST(DOFS_MAPPING%GHOST_START:DOFS_MAPPING%GHOST_FINISH), &
          & '("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
      ENDIF
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
        & DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
        & DOFS_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS_PTR( &
        & DOFS_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,DOFS_MAPPING%ADJACENT_DOMAINS_LIST, &
        '("    Adjacent domains list :",8(X,I7))','(27X,8(X,I7))',ERR,ERROR,*999)
      DO domain_idx=1,DOFS_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_SEND_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
          & '("      Local send ghost indicies       :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
          & DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
          & NUMBER_OF_RECEIVE_GHOSTS,6,6,DOFS_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
          & '("      Local receive ghost indicies    :",6(X,I7))','(39X,6(X,I7))',ERR,ERROR,*999)              
      ENDDO !domain_idx
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE")
    RETURN
999 IF(ALLOCATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ALLOCATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ALLOCATED(GHOST_NODES)) DEALLOCATE(GHOST_NODES)
    IF(ALLOCATED(NUMBER_INTERNAL_NODES)) DEALLOCATE(NUMBER_INTERNAL_NODES)
    IF(ALLOCATED(NUMBER_BOUNDARY_NODES)) DEALLOCATE(NUMBER_BOUNDARY_NODES)
    IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,DUMMY_ERR,DUMMY_ERROR,*997)
997 CALL ERRORS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the node mapping in the given domain mappings. \todo pass in the nodes mapping
  SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mapping to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
        DEALLOCATE(DOMAIN_MAPPINGS%NODES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_NODES_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_NODES_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the node mapping in the given domain mapping. \todo finalise on error
  SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS !<A pointer to the domain mappings to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL FLAG_ERROR("Domain nodes mappings are already associated.",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%NODES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings nodes.",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%NODES,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_MAPPINGS_NODES_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_NODES_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the domain topology.
  SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne,np
    TYPE(BASIS_TYPE), POINTER :: BASIS

    CALL ENTERS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Find maximum number of element parameters for all elements in the domain toplogy
      TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=-1
      DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        IF(ASSOCIATED(BASIS)) THEN
          IF(BASIS%NUMBER_OF_ELEMENT_PARAMETERS>TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS) &
            & TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=BASIS%NUMBER_OF_ELEMENT_PARAMETERS
        ELSE
          CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
        ENDIF
      ENDDO !ne
      !Find maximum number of derivatives for all nodes in the domain toplogy
      TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=-1
      DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
        IF(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES>TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES) &
            & TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
      ENDDO !np
      !Calculate the elements surrounding the nodes in the domain topology
      CALL DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_CALCULATE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_CALCULATE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Initialises the local domain topology from the mesh topology.
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to initialise the domain topology from the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: local_element,global_element,local_node,global_node,version_idx,derivative_idx,node_idx,dof_idx,component_idx
    INTEGER(INTG) :: ne,nn,nkk,INSERT_STATUS
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(MESH_NODES_TYPE), POINTER :: MESH_NODES
    TYPE(DOMAIN_DOFS_TYPE), POINTER :: DOMAIN_DOFS

    CALL ENTERS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR,*999)
    
    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
          IF(ASSOCIATED(DOMAIN%MESH)) THEN
            MESH=>DOMAIN%MESH
            component_idx=DOMAIN%MESH_COMPONENT_NUMBER
            IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
              MESH_ELEMENTS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS
              DOMAIN_ELEMENTS=>DOMAIN%TOPOLOGY%ELEMENTS
              MESH_NODES=>MESH%TOPOLOGY(component_idx)%PTR%NODES
              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
              DOMAIN_DOFS=>DOMAIN%TOPOLOGY%DOFS
              ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain elements elements.",ERR,ERROR,*999)
              DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_GLOBAL
              ALLOCATE(DOMAIN_NODES%NODES(DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain nodes nodes.",ERR,ERROR,*999)
              DOMAIN_NODES%NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_LOCAL
              DOMAIN_NODES%TOTAL_NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_NODES%NUMBER_OF_GLOBAL_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_GLOBAL
              ALLOCATE(DOMAIN_DOFS%DOF_INDEX(3,DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain dofs dof index.",ERR,ERROR,*999)
  !#######################################################VERSIONS END##########################################################
              DOMAIN_DOFS%NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL
              DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL
              DOMAIN_DOFS%NUMBER_OF_GLOBAL_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_GLOBAL
              !Loop over the domain nodes and calculate the parameters from the mesh nodes
              CALL TREE_CREATE_START(DOMAIN_NODES%NODES_TREE,ERR,ERROR,*999)
              CALL TREE_INSERT_TYPE_SET(DOMAIN_NODES%NODES_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
              CALL TREE_CREATE_FINISH(DOMAIN_NODES%NODES_TREE,ERR,ERROR,*999)               
              dof_idx=0
              DO local_node=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                CALL DOMAIN_TOPOLOGY_NODE_INITIALISE(DOMAIN_NODES%NODES(local_node),ERR,ERROR,*999)
                global_node=DOMAIN%MAPPINGS%NODES%LOCAL_TO_GLOBAL_MAP(local_node)
                DOMAIN_NODES%NODES(local_node)%LOCAL_NUMBER=local_node
                DOMAIN_NODES%NODES(local_node)%MESH_NUMBER=global_node
                DOMAIN_NODES%NODES(local_node)%GLOBAL_NUMBER=MESH_NODES%NODES(global_node)%GLOBAL_NUMBER
                DOMAIN_NODES%NODES(local_node)%USER_NUMBER=MESH_NODES%NODES(global_node)%USER_NUMBER
                CALL TREE_ITEM_INSERT(DOMAIN_NODES%NODES_TREE,DOMAIN_NODES%NODES(local_node)%USER_NUMBER,local_node, &
                  & INSERT_STATUS,ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_SURROUNDING_ELEMENTS=0
                NULLIFY(DOMAIN_NODES%NODES(local_node)%SURROUNDING_ELEMENTS)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES=MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES
                !ALLOCATE(DOMAIN_NODES%NODES(local_node)%GLOBAL_DERIVATIVE_INDEX( &
                !  & MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR) !TODO: VERSIONS Delete line
                !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global derivative index.",ERR,ERROR,*999)
                !DOMAIN_NODES%NODES(local_node)%GLOBAL_DERIVATIVE_INDEX=MESH_NODES%NODES(global_node)%GLOBAL_DERIVATIVE_INDEX !TODO: VERSIONS Delete line
                !ALLOCATE(DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX( &
                !  & MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR) !TODO: VERSIONS Delete line
                !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node partial derivative index.",ERR,ERROR,*999)
                !DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX=MESH_NODES%NODES(global_node)%PARTIAL_DERIVATIVE_INDEX !TODO: VERSIONS Delete line
                !ALLOCATE(DOMAIN_NODES%NODES(local_node)%DOF_INDEX(MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node dof index.",ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain node derivatives.",ERR,ERROR,*999)
  !#######################################################VERSIONS END##########################################################
                DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                  CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE( &
                    & DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX= &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX
                  DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS= & 
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                  ALLOCATE(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%DOF_INDEX( &
                    & MESH_NODES%NODES(global_node)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node dervative versions dof index.",ERR,ERROR,*999)
                  DO version_idx=1,DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                    dof_idx=dof_idx+1
                    DOMAIN_NODES%NODES(local_node)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)=dof_idx
                    DOMAIN_DOFS%DOF_INDEX(1,dof_idx)=version_idx
                    DOMAIN_DOFS%DOF_INDEX(2,dof_idx)=derivative_idx
                    DOMAIN_DOFS%DOF_INDEX(3,dof_idx)=local_node
                  ENDDO !version_idx
  !#######################################################VERSIONS END##########################################################
                ENDDO !derivative_idx
                DOMAIN_NODES%NODES(local_node)%BOUNDARY_NODE=MESH_NODES%NODES(global_node)%BOUNDARY_NODE
              ENDDO !local_node
              !Loop over the domain elements and renumber from the mesh elements
              DO local_element=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                CALL DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(DOMAIN_ELEMENTS%ELEMENTS(local_element),ERR,ERROR,*999)
                global_element=DOMAIN%MAPPINGS%ELEMENTS%LOCAL_TO_GLOBAL_MAP(local_element)               
                BASIS=>MESH_ELEMENTS%ELEMENTS(global_element)%BASIS
                DOMAIN_ELEMENTS%ELEMENTS(local_element)%BASIS=>BASIS
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain elements element nodes",ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(2,BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                  & BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain elements element derivatives",ERR,ERROR,*999)
  !#######################################################VERSIONS END##########################################################
                DO nn=1,BASIS%NUMBER_OF_NODES
                  global_node=MESH_ELEMENTS%ELEMENTS(global_element)%MESH_ELEMENT_NODES(nn)
                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(global_node)%LOCAL_NUMBER(1)
                  DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(nn)=local_node
                  DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                    !Find equivalent node derivative by matching partial derivative index
  !#######################################################VERSIONS START########################################################
                    !/todo Take a look at this later- is it needed?
                    FOUND=.FALSE.
                    DO nkk=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                      IF(DOMAIN_NODES%NODES(local_node)%DERIVATIVES(nkk)%PARTIAL_DERIVATIVE_INDEX == &
                        & BASIS%PARTIAL_DERIVATIVE_INDEX(derivative_idx,nn)) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nkk
                    IF(FOUND) THEN
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(1,derivative_idx,nn)=nkk
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(2,derivative_idx,nn) = & 
                        & MESH_ELEMENTS%ELEMENTS(global_element)%USER_ELEMENT_NODE_VERSIONS(derivative_idx,nn)
  !#######################################################VERSIONS END##########################################################
                    ELSE
                      CALL FLAG_ERROR("Could not find equivalent node derivative",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !derivative_idx
                ENDDO !nn
              ENDDO !local_element
            ELSE
              CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)

          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Initialised domain topology :",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain nodes = ",DOMAIN_NODES%TOTAL_NUMBER_OF_NODES, &
        & ERR,ERROR,*999)
      DO node_idx=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",DOMAIN_NODES%NODES(node_idx)%LOCAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",DOMAIN_NODES%NODES(node_idx)%MESH_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",DOMAIN_NODES%NODES(node_idx)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",DOMAIN_NODES%NODES(node_idx)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES,ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
        DO derivative_idx=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node local derivative number = ",derivative_idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Global derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Partial derivative index = ", &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS,4,4, &
            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX, &
            & '("        Degree-of-freedom index(version_idx)  :",4(X,I9))','(36X,4(X,I9))',ERR,ERROR,*999)
        ENDDO !derivative_idx
  !#######################################################VERSIONS END##########################################################
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary node = ", &
          & DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE,ERR,ERROR,*999)
     ENDDO !node_idx
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total Number of domain dofs = ",DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS, &
        & ERR,ERROR,*999)
      DO dof_idx=1,DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dof number = ",dof_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,3,3, &
          & DOMAIN_DOFS%DOF_INDEX(:,dof_idx),'("    Degree-of-freedom index :",3(X,I9))','(29X,3(X,I9))', &
          & ERR,ERROR,*999)
      ENDDO !dof_idx
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of domain elements = ", &
        & DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element number = ",DOMAIN_ELEMENTS%ELEMENTS(ne)%NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis user number = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of local nodes = ", &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
          & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES,'("    Element nodes(nn) :",8(X,I9))','(23X,8(X,I9))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Element derivatives :",ERR,ERROR,*999)
        DO nn=1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Local node number = ",nn,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(1,:,nn), &
            & '("        Element derivatives(derivative_idx) :",8(X,I2))','(33X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
  !#######################################################VERSIONS START########################################################
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"    Element versions :",ERR,ERROR,*999)
        DO nn=1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Local node number = ",nn,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_DERIVATIVES(nn),8,8, &
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(2,:,nn), &
            & '("        Element derivatives(version_idx) :",8(X,I2))','(33X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
  !#######################################################VERSIONS END##########################################################
      ENDDO !ne
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH")
    RETURN 1   
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH
  
  !
  !================================================================================================================================
  !

  !>Finalises the dofs in the given domain topology. \todo pass in the dofs topolgy
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        IF(ALLOCATED(TOPOLOGY%DOFS%DOF_INDEX)) DEALLOCATE(TOPOLOGY%DOFS%DOF_INDEX)
        DEALLOCATE(TOPOLOGY%DOFS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_DOFS_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_DOFS_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the dofs data structures for a domain topology. \todo finalise on exit
  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        CALL FLAG_ERROR("Decomposition already has topology dofs associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology dofs",ERR,ERROR,*999)
        TOPOLOGY%DOFS%NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%TOTAL_NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%NUMBER_OF_GLOBAL_DOFS=0
        TOPOLOGY%DOFS%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_DOFS_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%ELEMENT_NODES)) DEALLOCATE(ELEMENT%ELEMENT_NODES)
 
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENT_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology element.
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT !<The domain element to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%NUMBER=0
    NULLIFY(ELEMENT%BASIS)
  
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENT_INITALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the elements in the given domain topology. \todo pass in the domain elements
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    CALL ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
          CALL DOMAIN_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the element data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FLAG_ERROR("Decomposition already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%NUMBER_OF_GLOBAL_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        TOPOLOGY%ELEMENTS%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS=0
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given domain. \todo pass in domain topology
  SUBROUTINE DOMAIN_TOPOLOGY_FINALISE(DOMAIN,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !<A pointer to the domain to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_TOPOLOGY_NODES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_DOFS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_LINES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_FACES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%TOPOLOGY)
    ELSE
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_TOPOLOGY_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given domain. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE(DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN !A pointer to the domain to initialise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
        CALL FLAG_ERROR("Domain already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate domain topology
        ALLOCATE(DOMAIN%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Domain topology could not be allocated",ERR,ERROR,*999)
        DOMAIN%TOPOLOGY%DOMAIN=>DOMAIN
        NULLIFY(DOMAIN%TOPOLOGY%ELEMENTS)
        NULLIFY(DOMAIN%TOPOLOGY%NODES)
        NULLIFY(DOMAIN%TOPOLOGY%DOFS)
        NULLIFY(DOMAIN%TOPOLOGY%LINES)
        NULLIFY(DOMAIN%TOPOLOGY%FACES)
        !Initialise the topology components
        CALL DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_NODES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_DOFS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_LINES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_FACES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        !Initialise the domain topology from the domain mappings and the mesh it came from
        CALL DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*999)
        !Calculate the topological information.
        CALL DOMAIN_TOPOLOGY_CALCULATE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises a line in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    IF(ALLOCATED(LINE%NODES_IN_LINE)) DEALLOCATE(LINE%NODES_IN_LINE)
    IF(ALLOCATED(LINE%DERIVATIVES_IN_LINE)) DEALLOCATE(LINE%DERIVATIVES_IN_LINE)
 
    CALL EXITS("DOMAIN_TOPOLOGY_LINE_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_LINE_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_LINE_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structure for a domain topology line.
  SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE !<The domain line to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    NULLIFY(LINE%BASIS)
    
    CALL EXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_LINE_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the lines in the given domain topology. \todo pass in domain lines
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nl
    
    CALL ENTERS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        DO nl=1,TOPOLOGY%LINES%NUMBER_OF_LINES
          CALL DOMAIN_TOPOLOGY_LINE_FINALISE(TOPOLOGY%LINES%LINES(nl),ERR,ERROR,*999)
        ENDDO !nl
        IF(ALLOCATED(TOPOLOGY%LINES%LINES)) DEALLOCATE(TOPOLOGY%LINES%LINES)
        DEALLOCATE(TOPOLOGY%LINES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_TOPOLOGY_LINES_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_LINES_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_LINES_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the line data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the lines for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FLAG_ERROR("Decomposition already has topology lines associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology lines",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_LINES_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_LINES_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE
  
  !
  !================================================================================================================================
  !
  !>Finalises a face in the given domain topology and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    IF(ALLOCATED(FACE%NODES_IN_FACE)) DEALLOCATE(FACE%NODES_IN_FACE)
    IF(ALLOCATED(FACE%DERIVATIVES_IN_FACE)) DEALLOCATE(FACE%DERIVATIVES_IN_FACE)
 
    CALL EXITS("DOMAIN_TOPOLOGY_FACE_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_FACE_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_FACE_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structure for a domain topology face.
  SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE(FACE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_FACE_TYPE) :: FACE !<The domain face to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR,*999)

    FACE%NUMBER=0
    NULLIFY(FACE%BASIS)
    
    CALL EXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_FACE_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_FACE_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACE_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finalises the faces in the given domain topology. \todo pass in domain faces
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to finalise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nf
    
    CALL ENTERS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        DO nf=1,TOPOLOGY%FACES%NUMBER_OF_FACES
          CALL DOMAIN_TOPOLOGY_FACE_FINALISE(TOPOLOGY%FACES%FACES(nf),ERR,ERROR,*999)
        ENDDO !nf
        IF(ALLOCATED(TOPOLOGY%FACES%FACES)) DEALLOCATE(TOPOLOGY%FACES%FACES)
        DEALLOCATE(TOPOLOGY%FACES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_TOPOLOGY_FACES_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_FACES_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_FACES_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the face data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the faces for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%FACES)) THEN
        CALL FLAG_ERROR("Decomposition already has topology faces associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%FACES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology faces",ERR,ERROR,*999)
        TOPOLOGY%FACES%NUMBER_OF_FACES=0
        TOPOLOGY%FACES%DOMAIN=>TOPOLOGY%DOMAIN
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_FACES_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_FACES_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_FACES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a domain nodes topolgoy. 
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(DOMAIN_TOPOLOGY,USER_NODE_NUMBER,NODE_EXISTS,DOMAIN_LOCAL_NODE_NUMBER, &
    & GHOST_NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY !<A pointer to the domain topology to check the node exists on
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the domain nodes topolgoy, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_LOCAL_NODE_NUMBER !<On exit, if the node exists the local number corresponding to the user node number. If the node does not exist then global number will be 0.
    LOGICAL, INTENT(OUT) :: GHOST_NODE !<On exit, is .TRUE. if the local node (if it exists) is a ghost node, .FALSE. if not. 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    
    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    DOMAIN_LOCAL_NODE_NUMBER=0
    GHOST_NODE=.FALSE.
    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
      IF(ASSOCIATED(DOMAIN_NODES)) THEN
        NULLIFY(TREE_NODE)
        CALL TREE_SEARCH(DOMAIN_NODES%NODES_TREE,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
        IF(ASSOCIATED(TREE_NODE)) THEN
          CALL TREE_NODE_VALUE_GET(DOMAIN_NODES%NODES_TREE,TREE_NODE,DOMAIN_LOCAL_NODE_NUMBER,ERR,ERROR,*999)
          NODE_EXISTS=.TRUE.
          GHOST_NODE=DOMAIN_LOCAL_NODE_NUMBER>DOMAIN_NODES%NUMBER_OF_NODES
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain topology is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN 1
    
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !
  !#######################################################VERSIONS START########################################################

  !>Finalises the given domain topology node derivative and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE_DERIVATIVE%USER_VERSION_NUMBERS)) DEALLOCATE(NODE_DERIVATIVE%USER_VERSION_NUMBERS)
    IF(ALLOCATED(NODE_DERIVATIVE%LOCAL_VERSION_NUMBERS)) DEALLOCATE(NODE_DERIVATIVE%LOCAL_VERSION_NUMBERS)
    IF(ALLOCATED(NODE_DERIVATIVE%DOF_INDEX)) DEALLOCATE(NODE_DERIVATIVE%DOF_INDEX)

    CALL EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR,*999)

    NODE_DERIVATIVE%NUMBER_OF_VERSIONS=0
    NODE_DERIVATIVE%GLOBAL_DERIVATIVE_INDEX=0
    NODE_DERIVATIVE%PARTIAL_DERIVATIVE_INDEX=0

    CALL EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_DERIVATIVE_INITIALISE

  !#######################################################VERSIONS END##########################################################

  !
  !================================================================================================================================
  !

  !>Finalises the given domain topology node and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

!  !#######################################################VERSIONS START########################################################
    DO derivative_idx=1,NODE%NUMBER_OF_DERIVATIVES
       CALL DOMAIN_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
    ENDDO !derivative_idx
    !IF(ALLOCATED(NODE%GLOBAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%GLOBAL_DERIVATIVE_INDEX) !TODO: VERSIONS Delete line
    !IF(ALLOCATED(NODE%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%PARTIAL_DERIVATIVE_INDEX) !TODO: VERSIONS Delete line
!  !#######################################################VERSIONS END##########################################################
    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(NODE%NODE_LINES)) DEALLOCATE(NODE%NODE_LINES)
 
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given domain topology node.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%LOCAL_NUMBER=0
    NODE%MESH_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%USER_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NODE%NUMBER_OF_NODE_LINES=0
    NODE%BOUNDARY_NODE=.FALSE.
    
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the nodees in the given domain topology. \todo pass in domain nodes
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    CALL ENTERS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
          CALL DOMAIN_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) DEALLOCATE(TOPOLOGY%NODES%NODES)
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES_TREE)) CALL TREE_DESTROY(TOPOLOGY%NODES%NODES_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%NODES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("DOMAIN_TOPOLOGY_NODES_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODES_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODES_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes data structures for a domain topology. \todo finalise on error
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        CALL FLAG_ERROR("Decomposition already has topology nodes associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%NODES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology nodes",ERR,ERROR,*999)
        TOPOLOGY%NODES%NUMBER_OF_NODES=0
        TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES=0
        TOPOLOGY%NODES%NUMBER_OF_GLOBAL_NODES=0
        TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=0
        TOPOLOGY%NODES%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%NODES%NODES)
        NULLIFY(TOPOLOGY%NODES%NODES_TREE)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODES_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODES_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Calculates the element numbers surrounding a node for a domain.
  SUBROUTINE DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the domain topology to calculate the elements surrounding the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_no,insert_position,ne,nn,np,surrounding_elem_no
    INTEGER(INTG), POINTER :: NEW_SURROUNDING_ELEMENTS(:)
    LOGICAL :: FOUND_ELEMENT
    TYPE(BASIS_TYPE), POINTER :: BASIS

    NULLIFY(NEW_SURROUNDING_ELEMENTS)

    CALL ENTERS("DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) THEN
            DO np=1,TOPOLOGY%NODES%TOTAL_NUMBER_OF_NODES
              TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS=0
              IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) &
                & DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
            ENDDO !np
            DO ne=1,TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
              BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
              DO nn=1,BASIS%NUMBER_OF_NODES
                np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                FOUND_ELEMENT=.FALSE.
                element_no=1
                insert_position=1
                DO WHILE(element_no<=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.FOUND_ELEMENT)
                  surrounding_elem_no=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(element_no)
                  IF(surrounding_elem_no==ne) THEN
                    FOUND_ELEMENT=.TRUE.
                  ENDIF
                  element_no=element_no+1
                  IF(ne>=surrounding_elem_no) THEN
                    insert_position=element_no
                  ENDIF
                ENDDO
                IF(.NOT.FOUND_ELEMENT) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(NEW_SURROUNDING_ELEMENTS(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new surrounding elements",ERR,ERROR,*999)
                  IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) THEN
                    NEW_SURROUNDING_ELEMENTS(1:insert_position-1)=TOPOLOGY%NODES%NODES(np)% &
                      & SURROUNDING_ELEMENTS(1:insert_position-1)
                    NEW_SURROUNDING_ELEMENTS(insert_position)=ne
                    NEW_SURROUNDING_ELEMENTS(insert_position+1:TOPOLOGY%NODES%NODES(np)% &
                      & NUMBER_OF_SURROUNDING_ELEMENTS+1)=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(insert_position: &
                      & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS)
                    DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
                  ELSE
                    NEW_SURROUNDING_ELEMENTS(1)=ne
                  ENDIF
                  TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS=>NEW_SURROUNDING_ELEMENTS
                  TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS= &
                    & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1
                ENDIF
              ENDDO !nn
            ENDDO !ne
          ELSE
            CALL FLAG_ERROR("Domain topology nodes nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain topology nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain topology is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE")
    RETURN
999 IF(ASSOCIATED(NEW_SURROUNDING_ELEMENTS)) DEALLOCATE(NEW_SURROUNDING_ELEMENTS)
    CALL ERRORS("DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE")
    RETURN 1   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finialises the mesh adjacent elements information and deallocates all memory
  SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ADJACENT_ELEMENT_TYPE) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
    IF(ALLOCATED(MESH_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(MESH_ADJACENT_ELEMENT%ADJACENT_ELEMENTS)
       
    CALL EXITS("MESH_ADJACENT_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("MESH_ADJACENT_ELEMENT_FINALISE",ERR,ERROR)    
    CALL EXITS("MESH_ADJACENT_ELEMENT_FINALISE")
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !
  !>Initalises the mesh adjacent elements information.
  SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE(MESH_ADJACENT_ELEMENT,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESH_ADJACENT_ELEMENT_TYPE) :: MESH_ADJACENT_ELEMENT !<The mesh adjacent element to initialise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR,*999)

    MESH_ADJACENT_ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS=0
       
    CALL EXITS("MESH_ADJACENT_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_ADJACENT_ELEMENT_INITIALISE",ERR,ERROR)    
    CALL EXITS("MESH_ADJACENT_ELEMENT_INITIALISE")
    RETURN 1
   
  END SUBROUTINE MESH_ADJACENT_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a mesh. \see OPENCMISS::CMISSMeshCreateFinish
  SUBROUTINE MESH_CREATE_FINISH(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,mesh_idx
    LOGICAL :: FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
        !Check that the mesh component elements have been finished
        FINISHED=.TRUE.
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
            IF(.NOT.MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS_FINISHED) THEN
              LOCAL_ERROR="The elements for mesh component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
                & " have not been finished"
              FINISHED=.FALSE.
              EXIT
            ENDIF
          ELSE
            LOCAL_ERROR="The elements for mesh topology component "//TRIM(NUMBER_TO_VSTRING(component_idx,"*",ERR,ERROR))// &
              & " are not associated"
            FINISHED=.FALSE.
            EXIT
          ENDIF
        ENDDO !component_idx
        IF(.NOT.FINISHED) CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        MESH%MESH_FINISHED=.TRUE.
        !Calulcate the mesh topology
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          CALL MESH_TOPOLOGY_CALCULATE(MESH%TOPOLOGY(component_idx)%PTR,MESH%SURROUNDING_ELEMENTS_CALCULATE,ERR,ERROR,*999)
        ENDDO !component_idx
      ELSE
        CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh number = ",mesh_idx,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ",MESH%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",MESH%GLOBAL_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ",MESH%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_CREATE_FINISH")
    RETURN
999 CALL ERRORS("MESH_CREATE_FINISH",ERR,ERROR)    
    CALL EXITS("MESH_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh 
  SUBROUTINE MESH_CREATE_START_GENERIC(MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The pointer to the meshes
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,mesh_idx
    TYPE(MESH_TYPE), POINTER :: NEW_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(NEW_MESH)
    NULLIFY(NEW_MESHES)

    CALL ENTERS("MESH_CREATE_START_GENERIC",ERR,ERROR,*997)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FLAG_ERROR("Mesh is already associated.",ERR,ERROR,*997)
      ELSE
        CALL MESH_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Set default mesh values
        NEW_MESH%USER_NUMBER=USER_NUMBER
        NEW_MESH%GLOBAL_NUMBER=MESHES%NUMBER_OF_MESHES+1
        NEW_MESH%MESHES=>MESHES
        NEW_MESH%NUMBER_OF_DIMENSIONS=NUMBER_OF_DIMENSIONS
        NEW_MESH%NUMBER_OF_COMPONENTS=1
        NEW_MESH%SURROUNDING_ELEMENTS_CALCULATE=.true. !default true
        !Initialise mesh topology and decompositions
        CALL MESH_TOPOLOGY_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        CALL DECOMPOSITIONS_INITIALISE(NEW_MESH,ERR,ERROR,*999)
        !Add new mesh into list of meshes 
        ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES+1),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new meshes",ERR,ERROR,*999)
        DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
          NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
        ENDDO !mesh_idx
        NEW_MESHES(MESHES%NUMBER_OF_MESHES+1)%PTR=>NEW_MESH
        IF(ASSOCIATED(MESHES%MESHES)) DEALLOCATE(MESHES%MESHES)
        MESHES%MESHES=>NEW_MESHES
        MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES+1
        MESH=>NEW_MESH
      ENDIF
    ELSE
      CALL FLAG_ERROR("Meshes is not associated.",ERR,ERROR,*997)
    ENDIF
      
    CALL EXITS("MESH_CREATE_START_GENERIC")
    RETURN
999 CALL MESH_FINALISE(NEW_MESH,DUMMY_ERR,DUMMY_ERROR,*998)
998 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    NULLIFY(MESH)    
997 CALL ERRORS("MESH_CREATE_START_GENERIC",ERR,ERROR)    
    CALL EXITS("MESH_CREATE_START_GENERIC")
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_GENERIC

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in an interface. \see OPENCMISS::CMISSMeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_INTERFACE(USER_NUMBER,INTERFACE,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_CREATE_START_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FLAG_ERROR("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(INTERFACE%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on interface number "// &
              & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(INTERFACE%INTERFACES)) THEN
              PARENT_REGION=>INTERFACE%INTERFACES%PARENT_REGION
              IF(ASSOCIATED(PARENT_REGION)) THEN
                IF(ASSOCIATED(PARENT_REGION%COORDINATE_SYSTEM)) THEN                  
                  IF(NUMBER_OF_DIMENSIONS>0) THEN
                    IF(NUMBER_OF_DIMENSIONS<=PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                      CALL MESH_CREATE_START_GENERIC(INTERFACE%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                      MESH%INTERFACE=>INTERFACE
                    ELSE
                      LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                        & ") must be <= number of parent region dimensions ("// &
                        & TRIM(NUMBER_TO_VSTRING(PARENT_REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Parent region coordinate system is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Interfaces parent region is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Interface interfaces is not associated.",ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_CREATE_START_INTERFACE")
    RETURN
999 CALL ERRORS("MESH_CREATE_START_INTERFACE",ERR,ERROR)    
    CALL EXITS("MESH_CREATE_START_INTERFACE")
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_INTERFACE

  !
  !================================================================================================================================
  !

  !>Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in the region identified by REGION. \see OPENCMISS::CMISSMeshCreateStart
  !>Default values set for the MESH's attributes are:
  !>- NUMBER_OF_COMPONENTS: 1
  SUBROUTINE MESH_CREATE_START_REGION(USER_NUMBER,REGION,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to create
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the mesh on
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of dimensions in the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On exit, a pointer to the created mesh. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_CREATE_START_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FLAG_ERROR("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        IF(ASSOCIATED(REGION%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
              IF(NUMBER_OF_DIMENSIONS>0) THEN
                IF(NUMBER_OF_DIMENSIONS<=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
                  CALL MESH_CREATE_START_GENERIC(REGION%MESHES,USER_NUMBER,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*999)
                  MESH%REGION=>REGION
                ELSE
                  LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                    & ") must be <= number of region dimensions ("// &
                    & TRIM(NUMBER_TO_VSTRING(REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Number of mesh dimensions must be > 0.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="The coordinate system on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
                & " are not associated."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MESH_CREATE_START_REGION")
    RETURN
999 CALL ERRORS("MESH_CREATE_START_REGION",ERR,ERROR)    
    CALL EXITS("MESH_CREATE_START_REGION")
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START_REGION

  !
  !================================================================================================================================
  !

  !>Destroys the mesh identified by a user number on the given region and deallocates all memory. \see OPENCMISS::CMISSMeshDestroy
  SUBROUTINE MESH_DESTROY_NUMBER(USER_NUMBER,REGION,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region containing the mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    CALL ENTERS("MESH_DESTROY_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN

!!TODO: have a mesh_destory_ptr and mesh_destroy_number
        
        !Find the problem identified by the user number
        FOUND=.FALSE.
        mesh_position=0
        DO WHILE(mesh_position<REGION%MESHES%NUMBER_OF_MESHES.AND..NOT.FOUND)
          mesh_position=mesh_position+1
          IF(REGION%MESHES%MESHES(mesh_position)%PTR%USER_NUMBER==USER_NUMBER) FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          MESH=>REGION%MESHES%MESHES(mesh_position)%PTR

          CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

          !Remove the mesh from the list of meshes
          IF(REGION%MESHES%NUMBER_OF_MESHES>1) THEN
            ALLOCATE(NEW_MESHES(REGION%MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new meshes",ERR,ERROR,*999)
            DO mesh_idx=1,REGION%MESHES%NUMBER_OF_MESHES
              IF(mesh_idx<mesh_position) THEN
                NEW_MESHES(mesh_idx)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ELSE IF(mesh_idx>mesh_position) THEN
                REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
                NEW_MESHES(mesh_idx-1)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ENDIF
            ENDDO !mesh_idx
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%MESHES=>NEW_MESHES
            REGION%MESHES%NUMBER_OF_MESHES=REGION%MESHES%NUMBER_OF_MESHES-1
          ELSE
            DEALLOCATE(REGION%MESHES%MESHES)
            REGION%MESHES%NUMBER_OF_MESHES=0
          ENDIF
          
        ELSE
          LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
            & " has not been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MESH_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    CALL ERRORS("MESH_DESTROY_NUMBER",ERR,ERROR)
    CALL EXITS("MESH_DESTROY_NUMBER")
    RETURN 1   
  END SUBROUTINE MESH_DESTROY_NUMBER

  !
  !================================================================================================================================
  !

  !>Destroys the mesh and deallocates all memory. \see OPENCMISS::CMISSMeshDestroy
  SUBROUTINE MESH_DESTROY(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to destroy.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    TYPE(MESHES_TYPE), POINTER :: MESHES
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    CALL ENTERS("MESH_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      MESHES=>MESH%MESHES
      IF(ASSOCIATED(MESHES)) THEN
        mesh_position=MESH%GLOBAL_NUMBER
          
        CALL MESH_FINALISE(MESH,ERR,ERROR,*999)

        !Remove the mesh from the list of meshes
        IF(MESHES%NUMBER_OF_MESHES>1) THEN
          ALLOCATE(NEW_MESHES(MESHES%NUMBER_OF_MESHES-1),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new meshes.",ERR,ERROR,*999)
          DO mesh_idx=1,MESHES%NUMBER_OF_MESHES
            IF(mesh_idx<mesh_position) THEN
              NEW_MESHES(mesh_idx)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ELSE IF(mesh_idx>mesh_position) THEN
              MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER=MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER-1
              NEW_MESHES(mesh_idx-1)%PTR=>MESHES%MESHES(mesh_idx)%PTR
            ENDIF
          ENDDO !mesh_idx
          DEALLOCATE(MESHES%MESHES)
          MESHES%MESHES=>NEW_MESHES
          MESHES%NUMBER_OF_MESHES=MESHES%NUMBER_OF_MESHES-1
        ELSE
          DEALLOCATE(MESHES%MESHES)
          MESHES%NUMBER_OF_MESHES=0
        ENDIF
      ELSE
        CALL FLAG_ERROR("The mesh meshes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MESH_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    CALL ERRORS("MESH_DESTROY",ERR,ERROR)
    CALL EXITS("MESH_DESTROY")
    RETURN 1   
  END SUBROUTINE MESH_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalises a mesh and deallocates all memory.
  SUBROUTINE MESH_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL MESH_TOPOLOGY_FINALISE(MESH,ERR,ERROR,*999)
      CALL DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*999)
!      IF(ASSOCIATED(MESH%INTF)) CALL INTERFACE_MESH_FINALISE(MESH,ERR,ERROR,*999)  ! <<??>>
      DEALLOCATE(MESH)
    ENDIF
 
    CALL EXITS("MESH_FINALISE")
    RETURN
999 CALL ERRORS("MESH_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_FINALISE")
    RETURN 1
   
  END SUBROUTINE MESH_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a mesh.
  SUBROUTINE MESH_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL FLAG_ERROR("Mesh is already associated.",ERR,ERROR,*999)
    ELSE
      ALLOCATE(MESH,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new mesh.",ERR,ERROR,*999)
      MESH%USER_NUMBER=0
      MESH%GLOBAL_NUMBER=0
      MESH%MESH_FINISHED=.FALSE.
      NULLIFY(MESH%MESHES)
      NULLIFY(MESH%REGION)
      NULLIFY(MESH%INTERFACE)
      NULLIFY(MESH%GENERATED_MESH)
      MESH%NUMBER_OF_DIMENSIONS=0
      MESH%NUMBER_OF_COMPONENTS=0
      MESH%MESH_EMBEDDED=.FALSE.
      NULLIFY(MESH%EMBEDDING_MESH)
      MESH%NUMBER_OF_EMBEDDED_MESHES=0
      NULLIFY(MESH%EMBEDDED_MESHES)
      MESH%NUMBER_OF_ELEMENTS=0
      NULLIFY(MESH%TOPOLOGY)
      NULLIFY(MESH%DECOMPOSITIONS)
    ENDIF
    
    CALL EXITS("MESH_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the number of mesh components for a mesh identified by a pointer. \see OPENCMISS::CMISSMeshNumberOfComponentsGet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of components for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_COMPONENTS !<On return, the number of components in the specified mesh.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_COMPONENTS=MESH%NUMBER_OF_COMPONENTS
      ELSE
        CALL FLAG_ERROR("Mesh has not finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_GET")
    RETURN
999 CALL ERRORS("MESH_NUMBER_OF_COMPONENTS_GET",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_GET")
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of mesh components for a mesh. \see OPENCMISS::CMISSMeshNumberOfComponentsSet
  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of components for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components to set.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TOPOLOGY_PTR_TYPE), POINTER :: NEW_TOPOLOGY(:)

    NULLIFY(NEW_TOPOLOGY)
    
    CALL ENTERS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_COMPONENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FLAG_ERROR("Mesh has been finished",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_COMPONENTS/=MESH%NUMBER_OF_COMPONENTS) THEN
            ALLOCATE(NEW_TOPOLOGY(NUMBER_OF_COMPONENTS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new topology",ERR,ERROR,*999)
            IF(NUMBER_OF_COMPONENTS<MESH%NUMBER_OF_COMPONENTS) THEN
              DO component_idx=1,NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
            ELSE !NUMBER_OF_COMPONENTS>MESH%NUMBER_OF_COMPONENTS
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                NEW_TOPOLOGY(component_idx)%PTR=>MESH%TOPOLOGY(component_idx)%PTR
              ENDDO !component_idx
!!TODO \todo sort out mesh_topology initialise/finalise so that they allocate and deal with this below then call that routine
              DO component_idx=MESH%NUMBER_OF_COMPONENTS+1,NUMBER_OF_COMPONENTS
                ALLOCATE(NEW_TOPOLOGY(component_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new topology component",ERR,ERROR,*999)
                NEW_TOPOLOGY(component_idx)%PTR%MESH=>MESH
                NEW_TOPOLOGY(component_idx)%PTR%MESH_COMPONENT_NUMBER=component_idx
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%ELEMENTS)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%NODES)
                NULLIFY(NEW_TOPOLOGY(component_idx)%PTR%DOFS)
                !Initialise the topology components
                CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_NODES_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_DOFS_INITIALISE(NEW_TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
            IF(ASSOCIATED(MESH%TOPOLOGY)) DEALLOCATE(MESH%TOPOLOGY)
            MESH%TOPOLOGY=>NEW_TOPOLOGY
            MESH%NUMBER_OF_COMPONENTS=NUMBER_OF_COMPONENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of mesh components ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_COMPONENTS,"*",ERR,ERROR))// &
          & ") is illegal. You must have >0 mesh components"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_SET")
    RETURN
!!TODO: tidy up memory deallocation on error
999 CALL ERRORS("MESH_NUMBER_OF_COMPONENTS_SET",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_SET")
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET

  !
  !================================================================================================================================
  !
  
  !>Gets the number of elements for a mesh identified by a pointer. \see OPENCMISS::CMISSMeshNumberOfElementsGet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the number of elements for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ELEMENTS !<On return, the number of elements in the specified mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS
      ELSE
        CALL FLAG_ERROR("Mesh has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
999 CALL ERRORS("MESH_NUMBER_OF_ELEMENTS_GET",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_GET")
    RETURN
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the number of elements for a mesh. \see OPENCMISS::CMISSMeshNumberOfElementsSet
  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the number of elements for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS !<The number of elements to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_ELEMENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FLAG_ERROR("Mesh has been finished.",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_ELEMENTS/=MESH%NUMBER_OF_ELEMENTS) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
                  IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
                    IF(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS>0) THEN
!!TODO: Reallocate the elements and copy information. 
                      CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Mesh topology component pointer is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*999)
            ENDIF
            MESH%NUMBER_OF_ELEMENTS=NUMBER_OF_ELEMENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of elements ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
          & ") is invalid. You must have > 0 elements."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN
999 CALL ERRORS("MESH_NUMBER_OF_ELEMENTS_SET",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_SET")
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET

  !
  !================================================================================================================================
  !

  !>Changes/sets the surrounding elements calculate flag. \see OPENCMISS::CMISSMeshSurroundingElementsCalculateSet
  SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET(MESH,SURROUNDING_ELEMENTS_CALCULATE_FLAG,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to set the surrounding elements calculate flag for
    LOGICAL, INTENT(IN) :: SURROUNDING_ELEMENTS_CALCULATE_FLAG !<The surrounding elements calculate flag
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        CALL FLAG_ERROR("Mesh has been finished.",ERR,ERROR,*999)
      ELSE
        MESH%SURROUNDING_ELEMENTS_CALCULATE=SURROUNDING_ELEMENTS_CALCULATE_FLAG
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET")
    RETURN
999 CALL ERRORS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET",ERR,ERROR)    
    CALL EXITS("MESH_SURROUNDING_ELEMENTS_CALCULATE_SET")
    RETURN 1
   
  END SUBROUTINE MESH_SURROUNDING_ELEMENTS_CALCULATE_SET

  !
  !================================================================================================================================
  !

  !>Calculates the mesh topology.
  SUBROUTINE MESH_TOPOLOGY_CALCULATE(TOPOLOGY,SURROUNDING_ELEMENTS_CALCULATE,ERR,ERROR,*)    

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate
    LOGICAL, INTENT(IN) :: SURROUNDING_ELEMENTS_CALCULATE !<Flag to determine whether surrounding elements should be calculated
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Calculate the nodes used in the mesh
      CALL MESH_TOPOLOGY_NODES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the elements surrounding the nodes in a mesh
      !IF(SURROUNDING_ELEMENTS_CALCULATE) THEN
        CALL MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the number of derivatives at each node in a mesh
      CALL MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the number of versions for each derivative at each node in a mesh
      CALL MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the elements surrounding the elements in the mesh
      !IF(SURROUNDING_ELEMENTS_CALCULATE) THEN
        CALL MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
        !Calculate the boundary nodes and elements in the mesh
        CALL MESH_TOPOLOGY_BOUNDARY_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !ENDIF
      !Calculate the elements surrounding the elements in the mesh
      CALL MESH_TOPOLOGY_DOFS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_CALCULATE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_CALCULATE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_CALCULATE
  
  !
  !===============================================================================================================================
  !

  !>Calculates the boundary nodes and elements for a mesh topology. 
  SUBROUTINE MESH_TOPOLOGY_BOUNDARY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the boundary for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_idx,MATCH_INDEX,nic,local_node_idx,node_idx,XI_DIRECTION
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MESH_TOPOLOGY_BOUNDARY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN

      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN        
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          DO element_idx=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BASIS
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                IF(nic/=0) THEN
                  IF(TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                    TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BOUNDARY_ELEMENT=.TRUE.
                    IF(nic<0) THEN
                      XI_DIRECTION=-nic
                      MATCH_INDEX=1
                    ELSE
                      XI_DIRECTION=nic
                      MATCH_INDEX=BASIS%NUMBER_OF_NODES_XIC(nic)
                    ENDIF                    
                    DO local_node_idx=1,BASIS%NUMBER_OF_NODES
                      IF(BASIS%NODE_POSITION_INDEX(local_node_idx,XI_DIRECTION)==MATCH_INDEX) THEN
                        node_idx=TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%MESH_ELEMENT_NODES(local_node_idx)
                        TOPOLOGY%NODES%NODES(node_idx)%BOUNDARY_NODE=.TRUE.
                      ENDIF
                    ENDDO !nn
                  ENDIF
                ENDIF
              ENDDO !nic            
            CASE(BASIS_SIMPLEX_TYPE)
              TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BOUNDARY_ELEMENT=.FALSE.
              DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BOUNDARY_ELEMENT=TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)% &
                  & BOUNDARY_ELEMENT.OR.TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%ADJACENT_ELEMENTS(nic)% &
                  & NUMBER_OF_ADJACENT_ELEMENTS==0
                IF(TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS==0) THEN
                  DO local_node_idx=1,BASIS%NUMBER_OF_NODES
                    IF(BASIS%NODE_POSITION_INDEX(local_node_idx,nic)==1) THEN
                      node_idx=TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%MESH_ELEMENT_NODES(local_node_idx)
                      TOPOLOGY%NODES%NODES(node_idx)%BOUNDARY_NODE=.TRUE.
                    ENDIF
                  ENDDO !local_node_idx
                  
                ENDIF
              ENDDO !nic
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDDO !element_idx
        ELSE
          CALL FLAG_ERROR("Topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated.",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Boundary elements:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS, &
        ERR,ERROR,*999)
      DO element_idx=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Element : ",element_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary element = ",TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)% &
          & BOUNDARY_ELEMENT,ERR,ERROR,*999)        
      ENDDO !element_idx
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Boundary nodes:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of nodes = ",TOPOLOGY%NODES%NUMBER_OF_NODES, &
        ERR,ERROR,*999)
      DO node_idx=1,TOPOLOGY%NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Node : ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary node = ",TOPOLOGY%NODES%NODES(node_idx)% &
          & BOUNDARY_NODE,ERR,ERROR,*999)        
      ENDDO !element_idx            
    ENDIF
 
    CALL EXITS("MESH_TOPOLOGY_BOUNDARY_CALCULATE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_BOUNDARY_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_BOUNDARY_CALCULATE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_BOUNDARY_CALCULATE

  !
  !===============================================================================================================================
  !

  !>Calculates the degrees-of-freedom for a mesh topology. 
  SUBROUTINE MESH_TOPOLOGY_DOFS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: node_idx,derivative_idx,version_idx,NUMBER_OF_DOFS

    CALL ENTERS("MESH_TOPOLOGY_DOFS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
          NUMBER_OF_DOFS=0
          DO node_idx=1,TOPOLOGY%NODES%NUMBER_OF_NODES
!  !#######################################################VERSIONS START########################################################
             !\todo DELETE THIS
!            ALLOCATE(TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES),STAT=ERR)
!            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh topology node dof index",ERR,ERROR,*999)
            DO derivative_idx=1,TOPOLOGY%NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
               !\todo DELETE THIS
!              NUMBER_OF_DOFS=NUMBER_OF_DOFS+1
!              TOPOLOGY%NODES%NODES(node_idx)%DOF_INDEX(derivative_idx)=NUMBER_OF_DOFS
              ALLOCATE(TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX( &
                & TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh topology node derivative version dof index",ERR,ERROR,*999)
              DO version_idx=1,TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS
                NUMBER_OF_DOFS=NUMBER_OF_DOFS+1
                TOPOLOGY%NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%DOF_INDEX(version_idx)=NUMBER_OF_DOFS
              ENDDO !version_idx
!  !#######################################################VERSIONS END##########################################################
            ENDDO !derivative_idx
          ENDDO !node_idx
          TOPOLOGY%DOFS%NUMBER_OF_DOFS=NUMBER_OF_DOFS
        ELSE
          CALL FLAG_ERROR("Topology dofs is not assocaited",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology nodes is not assocaited",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("MESH_TOPOLOGY_DOFS_CALCULATE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_DOFS_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_DOFS_CALCULATE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DOFS_CALCULATE

  !
  !===============================================================================================================================
  !

  !>Finalises the dof data structures for a mesh topology and deallocates any memory. \todo pass in dofs
  SUBROUTINE MESH_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to finalise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        DEALLOCATE(TOPOLOGY%DOFS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("MESH_TOPOLOGY_DOFS_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_DOFS_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_DOFS_FINALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DOFS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the dofs in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the dofs for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
        CALL FLAG_ERROR("Mesh already has topology dofs associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology dofs",ERR,ERROR,*999)
        TOPOLOGY%DOFS%NUMBER_OF_DOFS=0
        TOPOLOGY%DOFS%MESH=>TOPOLOGY%MESH
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_DOFS_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_DOFS_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_DOFS_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_DOFS_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating elements for a specified mesh component in a mesh topology. \see OPENCMISS::CMISSMeshElementsCreateFinish
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the mesh elements to finish creating
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(MESH_TYPE), POINTER :: MESH

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Mesh elements have already been finished.",ERR,ERROR,*999)
      ELSE        
        ELEMENTS%ELEMENTS_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      MESH=>ELEMENTS%MESH
      IF(ASSOCIATED(MESH)) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",MESH%NUMBER_OF_ELEMENTS, &
          & ERR,ERROR,*999)
        DO ne=1,MESH%NUMBER_OF_ELEMENTS
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",ne,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ",ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER, &
            & ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ",ELEMENTS%ELEMENTS(ne)%USER_NUMBER, &
            & ERR,ERROR,*999)
          IF(ASSOCIATED(ELEMENTS%ELEMENTS(ne)%BASIS)) THEN
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ",ELEMENTS%ELEMENTS(ne)%BASIS%USER_NUMBER, &
              & ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Basis is not associated.",ERR,ERROR,*999)
          ENDIF
          IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES)) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)% BASIS%NUMBER_OF_NODES,8,8, &
              & ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES,'("    User element nodes   =",8(X,I6))','(26X,8(X,I6))',ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("User element nodes are not associated.",ERR,ERROR,*999)
          ENDIF
          IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES)) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
              & ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES,'("    Global element nodes =",8(X,I6))','(26X,8(X,I6))',ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Global element nodes are not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !ne
      ELSE
        CALL FLAG_ERROR("Mesh elements mesh is not associated.",ERR,ERROR,*999)
      ENDIF
    ENDIF

    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH
    
  !
  !================================================================================================================================
  !

  !>Starts the process of creating elements in the mesh component identified by MESH and component_idx. The elements will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure. \see OPENCMISS::CMISSMeshElementsCreateStart
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMPONENT_NUMBER,BASIS,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to start creating the elements on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the default basis to use
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<On return, a pointer to the created mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,INSERT_STATUS,ne
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR
 
    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR,*999)
    
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FLAG_ERROR("Elements is already associated.",ERR,ERROR,*999)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
              IF(ASSOCIATED(ELEMENTS%ELEMENTS)) THEN
                CALL FLAG_ERROR("Mesh topology already has elements associated",ERR,ERROR,*998)
              ELSE
                IF(ASSOCIATED(BASIS)) THEN
                  MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
                  ALLOCATE(ELEMENTS%ELEMENTS(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate individual elements",ERR,ERROR,*999)
                  ELEMENTS%NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS !Psuedo inheritance of the number of elements
                  CALL TREE_CREATE_START(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  CALL TREE_INSERT_TYPE_SET(ELEMENTS%ELEMENTS_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                  CALL TREE_CREATE_FINISH(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
                  ELEMENTS%ELEMENTS_FINISHED=.FALSE.
                  !Set up the default values and allocate element structures
                  DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                    CALL MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=ne
                    ELEMENTS%ELEMENTS(ne)%USER_NUMBER=ne
                    CALL TREE_ITEM_INSERT(ELEMENTS%ELEMENTS_TREE,ne,ne,INSERT_STATUS,ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%BASIS=>BASIS
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate user element nodes",ERR,ERROR,*999)
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global element nodes",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES=1
                    ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES=1
  !#######################################################VERSIONS START########################################################
                    ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                      & BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global element nodes versions",ERR,ERROR,*999)
                    ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS = 1
  !#######################################################VERSIONS END##########################################################
                  ENDDO !ne
                ELSE
                  CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*998)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN
999 CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,DUMMY_ERR,DUMMY_ERROR,*998)
998 NULLIFY(ELEMENTS)
    CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroys the elements in a mesh topology. \todo as this is a user routine it should take a mesh pointer like create start and finish? Split this into destroy and finalise?
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the mesh elements to destroy 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_DESTROY")
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nic

    CALL ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODES)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%GLOBAL_ELEMENT_NODES)) DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%MESH_ELEMENT_NODES)) DEALLOCATE(ELEMENT%MESH_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) THEN
      DO nic=LBOUND(ELEMENT%ADJACENT_ELEMENTS,1),UBOUND(ELEMENT%ADJACENT_ELEMENTS,1)
        CALL MESH_ADJACENT_ELEMENT_FINALISE(ELEMENT%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
      ENDDO !nic
      DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    ENDIF
  
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_FINALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the mesh elements for a given mesh component.
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_GET(MESH,MESH_COMPONENT_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to get the elements for
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to get the elements for
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<On return, a pointer to the mesh elements
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_GET",ERR,ERROR,*998)
    
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(ELEMENTS)) THEN
          CALL FLAG_ERROR("Elements is already associated.",ERR,ERROR,*998)
        ELSE
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
              ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
            ELSE
              CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*998)
    ENDIF

    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_GET")
    RETURN
999 NULLIFY(ELEMENTS)
998 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_GET",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_GET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_GET

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology element.
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT !<The mesh element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    NULLIFY(ELEMENT%BASIS)
    ELEMENT%BOUNDARY_ELEMENT=.FALSE.
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: Take user number
  
  !>Gets the basis for a mesh element identified by a given global number. \todo should take user number \see OPENCMISS::CMISSMeshElementsBasisGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to get the basis for
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to get the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<On return, a pointer to the basis to get
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
          BASIS=>ELEMENT%BASIS
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the basis for a mesh element identified by a given global number. \todo should take user number
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the basis for
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to set the basis for \todo before number?
    TYPE(BASIS_TYPE), POINTER :: BASIS !<A pointer to the basis to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG), ALLOCATABLE :: OLD_USER_ELEMENT_NODES(:),OLD_GLOBAL_ELEMENT_NODES(:)
    TYPE(MESH_ELEMENT_TYPE), POINTER :: ELEMENT
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ASSOCIATED(BASIS)) THEN
            ELEMENT=>ELEMENTS%ELEMENTS(GLOBAL_NUMBER)
            IF(ELEMENT%BASIS%NUMBER_OF_NODES/=BASIS%NUMBER_OF_NODES) THEN
              !Reallocate the user and global element nodes
              ALLOCATE(OLD_USER_ELEMENT_NODES(ELEMENT%BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old user element nodes",ERR,ERROR,*999)
              ALLOCATE(OLD_GLOBAL_ELEMENT_NODES(ELEMENT%BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate old user element nodes",ERR,ERROR,*999)
              OLD_USER_ELEMENT_NODES=ELEMENT%USER_ELEMENT_NODES
              OLD_GLOBAL_ELEMENT_NODES=ELEMENT%GLOBAL_ELEMENT_NODES
              DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
              DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
              ALLOCATE(ELEMENT%USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate user element nodes",ERR,ERROR,*999)
              ALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global element nodes",ERR,ERROR,*999)
              IF(ELEMENT%BASIS%NUMBER_OF_NODES<BASIS%NUMBER_OF_NODES) THEN
                ELEMENT%USER_ELEMENT_NODES(1:ELEMENT%BASIS%NUMBER_OF_NODES)=OLD_USER_ELEMENT_NODES(1:ELEMENT%BASIS%NUMBER_OF_NODES)
                ELEMENT%GLOBAL_ELEMENT_NODES(1:ELEMENT%BASIS%NUMBER_OF_NODES)= &
                  & OLD_GLOBAL_ELEMENT_NODES(1:ELEMENT%BASIS%NUMBER_OF_NODES)
                ELEMENT%USER_ELEMENT_NODES(ELEMENT%BASIS%NUMBER_OF_NODES+1:BASIS%NUMBER_OF_NODES)=OLD_USER_ELEMENT_NODES(1)
                ELEMENT%GLOBAL_ELEMENT_NODES(ELEMENT%BASIS%NUMBER_OF_NODES+1:BASIS%NUMBER_OF_NODES)=OLD_GLOBAL_ELEMENT_NODES(1)
              ELSE !ELEMENT%BASIS%NUMBER_OF_NODES>BASIS%NUMBER_OF_NODES
                ELEMENT%USER_ELEMENT_NODES(1:BASIS%NUMBER_OF_NODES)=OLD_USER_ELEMENT_NODES(1:BASIS%NUMBER_OF_NODES)
                ELEMENT%GLOBAL_ELEMENT_NODES(1:BASIS%NUMBER_OF_NODES)=OLD_GLOBAL_ELEMENT_NODES(1:BASIS%NUMBER_OF_NODES)
              ENDIF
              DEALLOCATE(OLD_USER_ELEMENT_NODES)
              DEALLOCATE(OLD_GLOBAL_ELEMENT_NODES)
            ENDIF            
            ELEMENT%BASIS=>BASIS
          ELSE
            CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET
  
  !
  !================================================================================================================================
  !

!!MERGE: user number. Dont use a pointer or allocate.
  
  !>Gets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(OUT) :: USER_ELEMENT_NODES(:) !<On return, USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(.NOT.ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have not been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)>=SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1)) THEN
            USER_ELEMENT_NODES=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES
          ELSE
            LOCAL_ERROR="The size of USER_ELEMENT_NODES is too small. The supplied size is "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(USER_ELEMENT_NODES,1),"*",ERR,ERROR))//" and it needs to be >= "// &
            & TRIM(NUMBER_TO_VSTRING(SIZE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES,1),"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the element nodes for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODES(:) !<USER_ELEMENT_NODES(i). USER_ELEMENT_NODES(i) is the i'th user node number for the element
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: nn,NUMBER_OF_BAD_NODES,GLOBAL_NODE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:),BAD_NODES(:)
    LOGICAL :: ELEMENT_NODES_OK,NODE_EXISTS
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)==ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            IF(ASSOCIATED(ELEMENTS%MESH)) THEN              
              REGION=>ELEMENTS%MESH%REGION
              IF(ASSOCIATED(REGION)) THEN
                NODES=>REGION%NODES
              ELSE
                INTERFACE=>ELEMENTS%MESH%INTERFACE
                IF(ASSOCIATED(INTERFACE)) THEN
                  NODES=>INTERFACE%NODES
                  PARENT_REGION=>INTERFACE%PARENT_REGION
                  IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FLAG_ERROR("Mesh interface has no parent region.",ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                ENDIF    
              ENDIF
              IF(ASSOCIATED(NODES)) THEN
                ELEMENT_NODES_OK=.TRUE.
                ALLOCATE(GLOBAL_ELEMENT_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global element nodes.",ERR,ERROR,*999)
                ALLOCATE(BAD_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bad nodes.",ERR,ERROR,*999)
                NUMBER_OF_BAD_NODES=0
                DO nn=1,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES
                  CALL NODE_CHECK_EXISTS(NODES,USER_ELEMENT_NODES(nn),NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                  IF(NODE_EXISTS) THEN
                    GLOBAL_ELEMENT_NODES(nn)=GLOBAL_NODE_NUMBER
                  ELSE
                    NUMBER_OF_BAD_NODES=NUMBER_OF_BAD_NODES+1
                    BAD_NODES(NUMBER_OF_BAD_NODES)=USER_ELEMENT_NODES(nn)
                    ELEMENT_NODES_OK=.FALSE.
                  ENDIF
                ENDDO !nn
                IF(ELEMENT_NODES_OK) THEN
                  ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES=USER_ELEMENT_NODES
                  ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES=GLOBAL_ELEMENT_NODES
                ELSE
                  IF(NUMBER_OF_BAD_NODES==1) THEN
                    IF(ASSOCIATED(REGION)) THEN
                      LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                        & " is not defined in region "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                    ELSE
                      LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                        & " is not defined in interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
                        & " of parent region number "//TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                    ENDIF
                  ELSE
                    LOCAL_ERROR="The element user node number of "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))
                    DO nn=2,NUMBER_OF_BAD_NODES-1
                      LOCAL_ERROR=LOCAL_ERROR//","//TRIM(NUMBER_TO_VSTRING(BAD_NODES(nn),"*",ERR,ERROR))
                    ENDDO !nn
                    IF(ASSOCIATED(REGION)) THEN
                      LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                        & " are not defined in region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                    ELSE
                      LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                        & " are not defined in interface number "// &
                        & TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))//" of parent region number "// &
                        &  TRIM(NUMBER_TO_VSTRING(PARENT_REGION%USER_NUMBER,"*",ERR,ERROR))//"."
                    ENDIF
                  ENDIF
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                IF(ASSOCIATED(REGION)) THEN

                  CALL FLAG_ERROR("The elements mesh region does not have any associated nodes.",ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("The elements mesh interface does not have any associated nodes.",ERR,ERROR,*999) 
                ENDIF
              ENDIF
            ELSE
              CALL FLAG_ERROR("The elements do not have an associated mesh.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Number of element nodes does not match number of basis nodes for this element.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  !
  !================================================================================================================================
  !
  !#######################################################VERSIONS START########################################################

  !>Changes/sets an element node's version for a mesh element identified by a given global number. \todo specify by user number not global number \see OPENCMISS::CMISSMeshElementsNodesSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET(GLOBAL_NUMBER,ELEMENTS,VERSION_NUMBER,DERIVATIVE_NUMBER, &
      & USER_ELEMENT_NODE_INDEX,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the element to set the nodes for
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to set \todo before number?
    INTEGER(INTG), INTENT(IN) :: VERSION_NUMBER !<The version number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: DERIVATIVE_NUMBER !<The derivative number of the specified element node to set.
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODE_INDEX !< The node index of the specified element node to set a version for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: PARENT_REGION,REGION
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(USER_ELEMENT_NODE_INDEX>=1.AND.USER_ELEMENT_NODE_INDEX<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            IF(ASSOCIATED(ELEMENTS%MESH)) THEN
              REGION=>ELEMENTS%MESH%REGION
              IF(ASSOCIATED(REGION)) THEN
                NODES=>REGION%NODES
              ELSE
                INTERFACE=>ELEMENTS%MESH%INTERFACE
                IF(ASSOCIATED(INTERFACE)) THEN
                  NODES=>INTERFACE%NODES
                  PARENT_REGION=>INTERFACE%PARENT_REGION
                  IF(.NOT.ASSOCIATED(PARENT_REGION)) CALL FLAG_ERROR("Mesh interface has no parent region.",ERR,ERROR,*999)
                ELSE
                  CALL FLAG_ERROR("Elements mesh has no associated region or interface.",ERR,ERROR,*999)
                ENDIF    
              ENDIF
              IF(DERIVATIVE_NUMBER>=1.AND.DERIVATIVE_NUMBER<=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS% &
                & NUMBER_OF_DERIVATIVES(USER_ELEMENT_NODE_INDEX)) THEN !Check if the specified derivative exists
                IF(VERSION_NUMBER>=1) THEN !Check if the specified version is greater than 1
                  ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODE_VERSIONS(DERIVATIVE_NUMBER,USER_ELEMENT_NODE_INDEX) & 
                    & = VERSION_NUMBER
                  !TODO: There is redunancy in USER_ELEMENT_NODE_VERSIONS since it was allocated in MESH_TOPOLOGY_ELEMENTS_CREATE_START based on MAXIMUM_NUMBER_OF_DERIVATIVES for that elements basis:ALLOCATE(ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES,BASIS%NUMBER_OF_NODES),STAT=ERR)
                ELSE
                  LOCAL_ERROR="The specified node version number of "//TRIM(NUMBER_TO_VSTRING(VERSION_NUMBER,"*", & 
                    & ERR,ERROR))//" is invalid. The element node index should be greater than 1."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                LOCAL_ERROR="The specified node derivative number of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVE_NUMBER,"*", & 
                  & ERR,ERROR))//" is invalid. The element node derivative index should be between 1 and "// &
                  & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_DERIVATIVES(USER_ELEMENT_NODE_INDEX), &
                  & "*",ERR,ERROR))//"."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDIF
          ELSE
            LOCAL_ERROR="The specified element node index of "//TRIM(NUMBER_TO_VSTRING(USER_ELEMENT_NODE_INDEX,"*",ERR,ERROR))// &
              & " is invalid. The element node index should be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified global element number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The global element number should be between 1 and "// &
            & TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODE_VERSION_SET

  !
  !================================================================================================================================
  !

  !#######################################################VERSIONS END##########################################################

 !>Calculates the element numbers surrounding an element in a mesh topology.
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the elements adjacent to elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,nep1,ni,nic,nn,nn1,nn2,nn3,node_idx,np,np1,DUMMY_ERR,FACE_XI(2),FACE_XIC(3),NODE_POSITION_INDEX(4)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: NUMBER_SURROUNDING,NUMBER_OF_NODES_XIC(4)
    INTEGER(INTG), ALLOCATABLE :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:)
    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(NODE_MATCH_LIST)
    DO nic=-4,4
      NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
    ENDDO !nic
    
    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          !Loop over the global elements in the mesh         
          DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
          !%%%% first we initialize lists that are required to find the adjacent elements list
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
            DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              NULLIFY(ADJACENT_ELEMENTS_LIST(nic)%PTR)
              CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(nic)%PTR,5,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
            ENDDO !ni
            NUMBER_OF_NODES_XIC=1
            NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)=BASIS%NUMBER_OF_NODES_XIC(1:BASIS%NUMBER_OF_XI_COORDINATES)
            !Place the current element in the surrounding list
            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER,ERR,ERROR,*999)
            SELECT CASE(BASIS%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              !Determine the collapsed "faces" if any
              NODE_POSITION_INDEX=1
              !Loop over the face normals of the element
              DO ni=1,BASIS%NUMBER_OF_XI
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Reset the node_position_index in this xi direction
                NODE_POSITION_INDEX(ni)=1
                !Loop over the two faces with this normal
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni
                  FACE_COLLAPSED(xi_direction)=.FALSE.
                  DO j=1,2
                    xi_dir_check=FACE_XI(j)
                    IF(xi_dir_check<=BASIS%NUMBER_OF_XI) THEN
                      xi_dir_search=FACE_XI(3-j)
                      NODE_POSITION_INDEX(xi_dir_search)=1
                      XI_COLLAPSED=.TRUE.
                      DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XIC(xi_dir_search).AND.XI_COLLAPSED)
                        !Get the first local node along the xi check direction
                        NODE_POSITION_INDEX(xi_dir_check)=1
                        nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        !Get the second local node along the xi check direction
                        NODE_POSITION_INDEX(xi_dir_check)=2
                        nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        IF(nn1/=0.AND.nn2/=0) THEN
                          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn1)/= &
                            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn2)) XI_COLLAPSED=.TRUE.
                        ENDIF
                        NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                      ENDDO !xi_dir_search
                      IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                    ENDIF
                  ENDDO !j
                  NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                ENDDO !direction_index
              ENDDO !ni
              !Loop over the xi directions and calculate the surrounding elements
              DO ni=1,BASIS%NUMBER_OF_XI
                !Determine the xi directions that lie in this xi direction
                FACE_XI(1)=OTHER_XI_DIRECTIONS3(ni,2,1)
                FACE_XI(2)=OTHER_XI_DIRECTIONS3(ni,3,1)
                !Loop over the two faces                
                DO direction_index=-1,1,2
                  xi_direction=direction_index*ni                  
                  !Find nodes in the element on the appropriate face/line/point
                  NULLIFY(NODE_MATCH_LIST)
                  CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)

                  CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                  IF(direction_index==-1) THEN
                    NODE_POSITION_INDEX(ni)=1
                  ELSE
                    NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XIC(ni)
                  ENDIF
                  !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is also
                  !collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the central
                  !axis back to itself
                  IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                    !Do nothing - the match lists are already empty
                  ELSE
                    !Find the nodes to match and add them to the node match list
                    DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XI(1))
                      NODE_POSITION_INDEX(FACE_XI(1))=nn1
                      DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XI(2))
                        NODE_POSITION_INDEX(FACE_XI(2))=nn2
                        nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                        IF(nn/=0) THEN
                          np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                          CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !nn2
                    ENDDO !nn1
                  ENDIF
                  CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                  CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                  NUMBER_SURROUNDING=0
                  IF(NUMBER_NODE_MATCHES>0) THEN
                    !Find list of elements surrounding those nodes
                    np1=NODE_MATCHES(1)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                      ne1=TOPOLOGY%NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES ! should this be GLOBAL_ELEMENT_NODES?
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & MESH_ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                          NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDIF
                  IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                ENDDO !direction_index
              ENDDO !ni
            CASE(BASIS_SIMPLEX_TYPE)
              !Loop over the xi coordinates and calculate the surrounding elements
              DO nic=1,BASIS%NUMBER_OF_XI_COORDINATES
                !Find the other coordinates of the face/line/point
                FACE_XIC(1)=OTHER_XI_DIRECTIONS4(nic,1)
                FACE_XIC(2)=OTHER_XI_DIRECTIONS4(nic,2)
                FACE_XIC(3)=OTHER_XI_DIRECTIONS4(nic,3)
                !Find nodes in the element on the appropriate face/line/point
                NULLIFY(NODE_MATCH_LIST)
                CALL LIST_CREATE_START(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(NODE_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(NODE_MATCH_LIST,16,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(NODE_MATCH_LIST,ERR,ERROR,*999)
                NODE_POSITION_INDEX(nic)=1 !Furtherest away from node with the nic'th coordinate
                !Find the nodes to match and add them to the node match list
                DO nn1=1,NUMBER_OF_NODES_XIC(FACE_XIC(1))
                  NODE_POSITION_INDEX(FACE_XIC(1))=nn1
                  DO nn2=1,NUMBER_OF_NODES_XIC(FACE_XIC(2))
                    NODE_POSITION_INDEX(FACE_XIC(2))=nn2
                    DO nn3=1,NUMBER_OF_NODES_XIC(FACE_XIC(3))
                      NODE_POSITION_INDEX(FACE_XIC(3))=nn3
                      nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3), &
                        NODE_POSITION_INDEX(4))
                      IF(nn/=0) THEN
                        np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                        CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !nn3
                  ENDDO !nn2
                ENDDO !nn1
                CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH_AND_DESTROY(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                IF(NUMBER_NODE_MATCHES>0) THEN
                  !Find list of elements surrounding those nodes
                  DO node_idx=1,NUMBER_NODE_MATCHES
                    np1=NODE_MATCHES(node_idx)
                    DO nep1=1,TOPOLOGY%NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                      ne1=TOPOLOGY%NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                      IF(ne1/=ne) THEN !Don't want the current element
                        ! grab the nodes list for current and this surrouding elements
                        ! current face : NODE_MATCHES
                        ! candidate elem : TOPOLOGY%ELEMENTS%ELEMENTS(ne1)%MESH_ELEMENT_NODES 
                        ! if all of current face belongs to the candidate element, we will have found the neighbour
                        CALL LIST_SUBSET_OF(NODE_MATCHES(1:NUMBER_NODE_MATCHES),TOPOLOGY%ELEMENTS%ELEMENTS(ne1)% &
                          & MESH_ELEMENT_NODES,SUBSET,ERR,ERROR,*999)
                        IF(SUBSET) THEN
                          CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(nic)%PTR,ne1,ERR,ERROR,*999)
                        ENDIF
                      ENDIF
                    ENDDO !nep1
                  ENDDO !node_idx
                ENDIF
                IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
              ENDDO !nic
            CASE(BASIS_SERENDIPITY_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_AUXILLIARY_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_B_SPLINE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_FOURIER_LAGRANGE_HERMITE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE(BASIS_EXTENDED_LAGRANGE_TP_TYPE)
              CALL FLAG_ERROR("Not implemented.",ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The basis type of "//TRIM(NUMBER_TO_VSTRING(BASIS%TYPE,"*",ERR,ERROR))// &
                & " is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
            !Set the surrounding elements for this element
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
              & BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent elements.",ERR,ERROR,*999)
            DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              CALL MESH_ADJACENT_ELEMENT_INITIALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic),ERR,ERROR,*999)
              CALL LIST_REMOVE_DUPLICATES(ADJACENT_ELEMENTS_LIST(nic)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS,ERR,ERROR,*999)
              ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element adjacent elements.",ERR,ERROR,*999)
              TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS) = ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS)
              IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
            ENDDO !nic
          ENDDO !ne           
        ELSE
          CALL FLAG_ERROR("Mesh topology elements is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology nodes is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not allocated.",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global element number : ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi coordinates = ",BASIS%NUMBER_OF_XI_COORDINATES, &
          & ERR,ERROR,*999)
        DO nic=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi coordinate : ",nic,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,ERR,ERROR,*999)
          IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS>0) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
              & ADJACENT_ELEMENTS(nic)%NUMBER_OF_ADJACENT_ELEMENTS,8,8,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(nic)% &
              & ADJACENT_ELEMENTS,'("        Adjacent elements :",8(X,I8))','(30x,8(X,I8))',ERR,ERROR,*999)
          ENDIF
        ENDDO !nic
      ENDDO !ne
    ENDIF
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE")
    RETURN
999 IF(ALLOCATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ALLOCATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO nic=-4,4
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(nic)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(nic)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE")
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the elements data structures for a mesh topology and deallocates any memory. \todo pass in elements
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE(ELEMENTS,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the mesh topology to finalise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ne

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
        CALL MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
      ENDDO !ne
      DEALLOCATE(ELEMENTS%ELEMENTS)
      IF(ASSOCIATED(ELEMENTS%ELEMENTS_TREE)) CALL TREE_DESTROY(ELEMENTS%ELEMENTS_TREE,ERR,ERROR,*999)
      DEALLOCATE(ELEMENTS)
    ENDIF
 
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_FINALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the elements in a given mesh topology. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FLAG_ERROR("Mesh already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%MESH=>TOPOLOGY%MESH
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS_TREE)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE
  
  !
  !================================================================================================================================
  !

!!MERGE: ditto.
  
  !>Gets the user number for a global element identified by a given global number. \todo Check that the user number doesn't already exist. \see OPENCMISS::CMISSMeshElementsUserNumberGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<On return, a pointer to the elements to get the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_GET")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Returns the user number for a global element identified by a given global number. \see OPENCMISS::CMISSMeshElementsUserNumberGet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to get.
    INTEGER(INTG), INTENT(OUT) :: USER_NUMBER !<The user number of the element to get
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN          
          USER_NUMBER=ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Elements have not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET")
    RETURN 1

   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET

  !
  !================================================================================================================================
  !

  !>Changes/sets the user number for a global element identified by a given global number. \see OPENCMISS::CMISSMeshElementsUserNumberSet
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number of the elements to set.
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the element to set
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS !<A pointer to the elements to set the user number for \todo This should be the first parameter.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code

    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GLOBAL_ELEMENT_NUMBER,INSERT_STATUS
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished.",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          NULLIFY(TREE_NODE)
          CALL TREE_SEARCH(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,TREE_NODE,ERR,ERROR,*999)
          IF(ASSOCIATED(TREE_NODE)) THEN
            CALL TREE_NODE_VALUE_GET(ELEMENTS%ELEMENTS_TREE,TREE_NODE,GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
            LOCAL_ERROR="Element user number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " is already used by global element number "// &
              & TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))//"."
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            CALL TREE_ITEM_DELETE(ELEMENTS%ELEMENTS_TREE,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER,ERR,ERROR,*999)
            CALL TREE_ITEM_INSERT(ELEMENTS%ELEMENTS_TREE,USER_NUMBER,GLOBAL_NUMBER,INSERT_STATUS,ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
          ENDIF
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET

  !
  !================================================================================================================================
  !

  !>Finalises the topology in the given mesh. \todo pass in the mesh topology
  SUBROUTINE MESH_TOPOLOGY_FINALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to finalise the topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx

    CALL ENTERS("MESH_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
        CALL MESH_TOPOLOGY_NODES_FINALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS,ERR,ERROR,*999)
        CALL MESH_TOPOLOGY_DOFS_FINALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        DEALLOCATE(MESH%TOPOLOGY(component_idx)%PTR)
      ENDDO !component_idx
      DEALLOCATE(MESH%TOPOLOGY)
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("MESH_TOPOLOGY_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_FINALISE")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the topology for a given mesh. \todo finalise on error
  SUBROUTINE MESH_TOPOLOGY_INITIALISE(MESH,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to initialise the mesh topology for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx
    
    CALL ENTERS("MESH_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
        CALL FLAG_ERROR("Mesh already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate mesh topology
        ALLOCATE(MESH%TOPOLOGY(MESH%NUMBER_OF_COMPONENTS),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Mesh topology could not be allocated",ERR,ERROR,*999)
        DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
          ALLOCATE(MESH%TOPOLOGY(component_idx)%PTR,STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Mesh topology component could not be allocated",ERR,ERROR,*999)
          MESH%TOPOLOGY(component_idx)%PTR%MESH=>MESH
          NULLIFY(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)
          NULLIFY(MESH%TOPOLOGY(component_idx)%PTR%NODES)
          NULLIFY(MESH%TOPOLOGY(component_idx)%PTR%DOFS)
          !Initialise the topology components
          CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
          CALL MESH_TOPOLOGY_NODES_INITIALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
          CALL MESH_TOPOLOGY_DOFS_INITIALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        ENDDO !component_idx
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Checks that a user element number exists in a mesh component. 
  SUBROUTINE MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS(MESH,MESH_COMPONENT_NUMBER,USER_ELEMENT_NUMBER,ELEMENT_EXISTS, &
    & MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to check the element exists on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component to check the element exits on
    INTEGER(INTG), INTENT(IN) :: USER_eLEMENT_NUMBER !<The user element number to check if it exists
    LOGICAL, INTENT(OUT) :: ELEMENT_EXISTS !<On exit, is .TRUE. if the element user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: MESH_GLOBAL_ELEMENT_NUMBER !<On exit, if the element exists the global number corresponding to the user element number. If the element does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR,*999)

    ELEMENT_EXISTS=.FALSE.
    MESH_GLOBAL_ELEMENT_NUMBER=0
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        IF(MESH_COMPONENT_NUMBER>=1.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
          IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
            MESH_TOPOLOGY=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
              MESH_ELEMENTS=>MESH_TOPOLOGY%ELEMENTS
              IF(ASSOCIATED(MESH_ELEMENTS)) THEN
                NULLIFY(TREE_NODE)
                CALL TREE_SEARCH(MESH_ELEMENTS%ELEMENTS_TREE,USER_ELEMENT_NUMBER,TREE_NODE,ERR,ERROR,*999)
                IF(ASSOCIATED(TREE_NODE)) THEN
                  CALL TREE_NODE_VALUE_GET(MESH_ELEMENTS%ELEMENTS_TREE,TREE_NODE,MESH_GLOBAL_ELEMENT_NUMBER,ERR,ERROR,*999)
                  ELEMENT_EXISTS=.TRUE.
                ENDIF
              ELSE
                CALL FLAG_ERROR("Mesh topology elements is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Mesh topology is not associated for component number "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified mesh component of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//" has "// &
            & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Checks that a user node number exists in a mesh component. 
  SUBROUTINE MESH_TOPOLOGY_NODE_CHECK_EXISTS(MESH,MESH_COMPONENT_NUMBER,USER_NODE_NUMBER,NODE_EXISTS, &
    & MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH !<A pointer to the mesh to check the node exists on
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component to check the node exits on
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The user node number to check if it exists
    LOGICAL, INTENT(OUT) :: NODE_EXISTS !<On exit, is .TRUE. if the node user number exists in the mesh component, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: MESH_GLOBAL_NODE_NUMBER !<On exit, if the node exists the global number corresponding to the user node number. If the node does not exist then global number will be 0.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_NODES_TYPE), POINTER :: MESH_NODES
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MESH_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR,*999)

    NODE_EXISTS=.FALSE.
    MESH_GLOBAL_NODE_NUMBER=0
    IF(ASSOCIATED(MESH)) THEN
      IF(MESH%MESH_FINISHED) THEN
        IF(MESH_COMPONENT_NUMBER>=1.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
          IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
            MESH_TOPOLOGY=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR
            IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
              MESH_NODES=>MESH_TOPOLOGY%NODES
              IF(ASSOCIATED(MESH_NODES)) THEN
                NULLIFY(TREE_NODE)
                CALL TREE_SEARCH(MESH_NODES%NODES_TREE,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
                IF(ASSOCIATED(TREE_NODE)) THEN
                  CALL TREE_NODE_VALUE_GET(MESH_NODES%NODES_TREE,TREE_NODE,MESH_GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                  NODE_EXISTS=.TRUE.
                ENDIF
              ELSE
                CALL FLAG_ERROR("Mesh topology nodes is not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Mesh topology is not associated for component number "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))//"."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          LOCAL_ERROR="The specified mesh component of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))//" has "// &
            & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))//"."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODE_CHECK_EXISTS",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODE_CHECK_EXISTS")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_NODE_CHECK_EXISTS
  
  !
  !================================================================================================================================
  !

  !>Finalises the given mesh topology node. 
  SUBROUTINE MESH_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_NODE_TYPE) :: NODE !<The mesh node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: derivative_idx

    CALL ENTERS("MESH_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

!  !#######################################################VERSIONS START########################################################
    DO derivative_idx=1,NODE%NUMBER_OF_DERIVATIVES
       CALL MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
    ENDDO !derivative_idx
!  !#######################################################VERSIONS END##########################################################
    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    !IF(ALLOCATED(NODE%GLOBAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%GLOBAL_DERIVATIVE_INDEX) !TODO: VERSIONS Delete line
    !IF(ALLOCATED(NODE%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%PARTIAL_DERIVATIVE_INDEX) !TODO: VERSIONS Delete line
  
    CALL EXITS("MESH_TOPOLOGY_NODE_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODE_FINALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MESH_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_NODE_TYPE) :: NODE !<The mesh node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%USER_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NULLIFY(NODE%SURROUNDING_ELEMENTS)
    NODE%NUMBER_OF_DERIVATIVES=0
    NODE%BOUNDARY_NODE=.FALSE.
    
    CALL EXITS("MESH_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODE_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Calculates the nodes used the mesh identified by a given mesh topology.
  SUBROUTINE MESH_TOPOLOGY_NODES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR,INSERT_STATUS,MESH_NUMBER,NUMBER_OF_MESH_NODES,ne,nn,np  
    INTEGER(INTG), POINTER :: MESH_NODES(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(MESH_NODES_TYPE), POINTER :: MESHNODES
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(TREE_TYPE), POINTER :: MESH_NODES_TREE
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(VARYING_STRING) :: DUMMY_ERROR,LOCAL_ERROR

    NULLIFY(MESH_NODES)
    NULLIFY(MESH_NODES_TREE)
    
    CALL ENTERS("MESH_TOPOLOGY_NODES_CALCULATE",ERR,ERROR,*998)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      ELEMENTS=>TOPOLOGY%ELEMENTS
      IF(ASSOCIATED(ELEMENTS)) THEN
        MESHNODES=>TOPOLOGY%NODES
        IF(ASSOCIATED(MESHNODES)) THEN
          MESH=>TOPOLOGY%MESH
          IF(ASSOCIATED(MESH)) THEN
            NULLIFY(INTERFACE)
            REGION=>MESH%REGION
            IF(ASSOCIATED(REGION)) THEN
              NODES=>REGION%NODES
            ELSE
              INTERFACE=>MESH%INTERFACE
              IF(ASSOCIATED(INTERFACE)) THEN
                NODES=>INTERFACE%NODES
              ELSE
                LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not have an associated region or interface."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDIF
            IF(ASSOCIATED(NODES)) THEN
              IF(ASSOCIATED(MESHNODES%NODES)) THEN
                LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                  & " already has associated mesh topology nodes."
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
              ELSE
                !Work out what nodes are in the mesh
                CALL TREE_CREATE_START(MESH_NODES_TREE,ERR,ERROR,*999)
                CALL TREE_INSERT_TYPE_SET(MESH_NODES_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                CALL TREE_CREATE_FINISH(MESH_NODES_TREE,ERR,ERROR,*999)
                DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                  BASIS=>ELEMENTS%ELEMENTS(ne)%BASIS
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    np=ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)
                    CALL TREE_ITEM_INSERT(MESH_NODES_TREE,np,np,INSERT_STATUS,ERR,ERROR,*999)
                  ENDDO !nn
                ENDDO !ne
                CALL TREE_DETACH_AND_DESTROY(MESH_NODES_TREE,NUMBER_OF_MESH_NODES,MESH_NODES,ERR,ERROR,*999)
                !Set up the mesh nodes.
                ALLOCATE(MESHNODES%NODES(NUMBER_OF_MESH_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh topology nodes nodes.",ERR,ERROR,*999)
                CALL TREE_CREATE_START(MESHNODES%NODES_TREE,ERR,ERROR,*999)
                CALL TREE_INSERT_TYPE_SET(MESHNODES%NODES_TREE,TREE_NO_DUPLICATES_ALLOWED,ERR,ERROR,*999)
                CALL TREE_CREATE_FINISH(MESHNODES%NODES_TREE,ERR,ERROR,*999) 
                DO np=1,NUMBER_OF_MESH_NODES
                  CALL MESH_TOPOLOGY_NODE_INITIALISE(MESHNODES%NODES(np),ERR,ERROR,*999)
                  MESHNODES%NODES(np)%MESH_NUMBER=np
                  MESHNODES%NODES(np)%GLOBAL_NUMBER=MESH_NODES(np)
                  MESHNODES%NODES(np)%USER_NUMBER=NODES%NODES(MESH_NODES(np))%USER_NUMBER
                  CALL TREE_ITEM_INSERT(MESHNODES%NODES_TREE,MESH_NODES(np),np,INSERT_STATUS,ERR,ERROR,*999)
                ENDDO !np
                MESHNODES%NUMBER_OF_NODES=NUMBER_OF_MESH_NODES
                IF(ASSOCIATED(MESH_NODES)) DEALLOCATE(MESH_NODES)
                !Now recalculate the mesh element nodes
                DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                  BASIS=>ELEMENTS%ELEMENTS(ne)%BASIS
                  ALLOCATE(ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh topology elements mesh element nodes.",ERR,ERROR,*999)
                  DO nn=1,BASIS%NUMBER_OF_NODES
                    np=ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)
                    NULLIFY(TREE_NODE)
                    CALL TREE_SEARCH(MESHNODES%NODES_TREE,np,TREE_NODE,ERR,ERROR,*999)
                    IF(ASSOCIATED(TREE_NODE)) THEN
                      CALL TREE_NODE_VALUE_GET(MESHNODES%NODES_TREE,TREE_NODE,MESH_NUMBER,ERR,ERROR,*999)
                      ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)=MESH_NUMBER
                    ELSE
                      LOCAL_ERROR="Could not find global node "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))//" (user node "// &
                        & TRIM(NUMBER_TO_VSTRING(NODES%NODES(np)%USER_NUMBER,"*",ERR,ERROR))//") in the mesh nodes."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !nn
                ENDDO !ne
              ENDIF
            ELSE
              IF(ASSOCIATED(REGION)) THEN
                LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not have any nodes associated with the mesh region."
              ELSE
                LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                  & " does not have any nodes associated with the mesh interface."
              ENDIF
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*998)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh topology mesh is not associated.",ERR,ERROR,*998)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh topology nodes is not associated.",ERR,ERROR,*998)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated.",ERR,ERROR,*998)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*998)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",MESHNODES%NUMBER_OF_NODES,ERR,ERROR,*999)
!Using "NODES%NUMBER_OF_NODES" seems to be a bug, since it exceeds the size of the array "MESHNODES%NODES(np)"
!       DO np=1,NODES%NUMBER_OF_NODES
      DO np=1,MESHNODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node number = ",MESHNODES%NODES(np)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)        
      ENDDO !np
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(MESH_NODES)) DEALLOCATE(MESH_NODES)
    IF(ASSOCIATED(MESH_NODES_TREE)) CALL TREE_DESTROY(MESH_NODES_TREE,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("MESH_TOPOLOGY_NODES_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_CALCULATE")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_NODES_CALCULATE

  !
  !================================================================================================================================
  !

  !#######################################################VERSIONS START########################################################

  !>Finalises the given mesh topology node. 
  SUBROUTINE MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The mesh node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE_DERIVATIVE%USER_VERSION_NUMBERS)) DEALLOCATE(NODE_DERIVATIVE%USER_VERSION_NUMBERS)
    IF(ALLOCATED(NODE_DERIVATIVE%DOF_INDEX)) DEALLOCATE(NODE_DERIVATIVE%DOF_INDEX)

    CALL EXITS("MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODE_DERIVATIVE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the given mesh topology node.
  SUBROUTINE MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODE_DERIVATIVE,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_NODE_DERIVATIVE_TYPE) :: NODE_DERIVATIVE !<The mesh node to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR,*999)

    NODE_DERIVATIVE%NUMBER_OF_VERSIONS=0
    NODE_DERIVATIVE%GLOBAL_DERIVATIVE_INDEX=0
    NODE_DERIVATIVE%PARTIAL_DERIVATIVE_INDEX=0

    CALL EXITS("MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE

  !#######################################################VERSIONS END##########################################################

  !
  !================================================================================================================================
  !

  !>Calculates the number of derivatives at each node in a topology.
  SUBROUTINE MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: elem_idx,global_deriv,MAX_NUMBER_OF_DERIVATIVES,ne,derivative_idx,nn,node_idx,NUMBER_OF_DERIVATIVES
    INTEGER(INTG), ALLOCATABLE :: DERIVATIVES(:)
    LOGICAL :: FOUND
    TYPE(LIST_TYPE), POINTER :: NODE_DERIVATIVE_LIST
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(MESH_NODES_TYPE), POINTER :: NODES
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE",ERR,ERROR,*999)

     IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          NODES=>TOPOLOGY%NODES
          !Loop over the mesh nodes
          DO node_idx=1,NODES%NUMBER_OF_NODES
            !VERSIONS - revert comment back to original
            !Calculate the number of derivatives and versions at each node. This needs to be calculated by looking at the mesh elements
            !as we may have an adjacent element in another domain with a higher order basis also with versions.
            NULLIFY(NODE_DERIVATIVE_LIST)
            CALL LIST_CREATE_START(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            CALL LIST_DATA_TYPE_SET(NODE_DERIVATIVE_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
            CALL LIST_INITIAL_SIZE_SET(NODE_DERIVATIVE_LIST,8,ERR,ERROR,*999)
            CALL LIST_CREATE_FINISH(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            MAX_NUMBER_OF_DERIVATIVES=-1
            DO elem_idx=1,NODES%NODES(node_idx)%NUMBER_OF_SURROUNDING_ELEMENTS
              ne=NODES%NODES(node_idx)%SURROUNDING_ELEMENTS(elem_idx)
              BASIS=>ELEMENTS%ELEMENTS(ne)%BASIS
              !Find the local node corresponding to this node
              FOUND=.FALSE.
              DO nn=1,BASIS%NUMBER_OF_NODES
                IF(ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)==node_idx) THEN
                  FOUND=.TRUE.
                  EXIT
                ENDIF
              ENDDO !nn
              IF(FOUND) THEN
                DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                  CALL LIST_ITEM_ADD(NODE_DERIVATIVE_LIST,BASIS%PARTIAL_DERIVATIVE_INDEX(derivative_idx,nn),ERR,ERROR,*999)
                ENDDO !derivative_idx
                IF(BASIS%NUMBER_OF_DERIVATIVES(nn)>MAX_NUMBER_OF_DERIVATIVES)MAX_NUMBER_OF_DERIVATIVES= &
                  & BASIS%NUMBER_OF_DERIVATIVES(nn)
              ELSE
                CALL FLAG_ERROR("Could not find local node.",ERR,ERROR,*999)
              ENDIF
            ENDDO !elem_idx
            CALL LIST_REMOVE_DUPLICATES(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            CALL LIST_DETACH_AND_DESTROY(NODE_DERIVATIVE_LIST,NUMBER_OF_DERIVATIVES,DERIVATIVES,ERR,ERROR,*999)
            IF(NUMBER_OF_DERIVATIVES==MAX_NUMBER_OF_DERIVATIVES) THEN
  !#######################################################VERSIONS START########################################################
              !Set up the node derivatives.
              ALLOCATE(NODES%NODES(node_idx)%DERIVATIVES(MAX_NUMBER_OF_DERIVATIVES),STAT=ERR)
  !#######################################################VERSIONS END##########################################################
              NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES=MAX_NUMBER_OF_DERIVATIVES
              !ALLOCATE(NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(MAX_NUMBER_OF_DERIVATIVES),STAT=ERR) !TODO: VERSIONS Delete line
              !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global derivative index.",ERR,ERROR,*999) !TODO: VERSIONS Delete line
              !ALLOCATE(NODES%NODES(node_idx)%PARTIAL_DERIVATIVE_INDEX(MAX_NUMBER_OF_DERIVATIVES),STAT=ERR) !TODO: VERSIONS Delete line
              !IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node partial derivative index.",ERR,ERROR,*999) !TODO: VERSIONS Delete line
              DO derivative_idx=1,NUMBER_OF_DERIVATIVES
  !#######################################################VERSIONS START########################################################
                CALL MESH_TOPOLOGY_NODE_DERIVATIVE_INITIALISE(NODES%NODES(node_idx)%DERIVATIVES(derivative_idx),ERR,ERROR,*999)
                NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX = DERIVATIVES(derivative_idx)
  !#######################################################VERSIONS END##########################################################
                !NODES%NODES(node_idx)%PARTIAL_DERIVATIVE_INDEX(derivative_idx)=DERIVATIVES(derivative_idx) !TODO: VERSIONS Delete line
                global_deriv=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(DERIVATIVES(derivative_idx))
                IF(global_deriv/=0) THEN
                  !NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX(derivative_idx)=global_deriv !TODO: VERSIONS Delete line
  !#######################################################VERSIONS START########################################################
                  NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX=global_deriv
  !#######################################################VERSIONS END##########################################################
                ELSE
                  LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVES(derivative_idx),"*", &
                    & ERR,ERROR))//" for derivative number "//TRIM(NUMBER_TO_VSTRING(derivative_idx,"*",ERR,ERROR))// &
                    & " does not have a corresponding global derivative."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !derivative_idx
              DEALLOCATE(DERIVATIVES)
            ELSE
              LOCAL_ERROR="Invalid mesh configuration. User node "// &
                & TRIM(NUMBER_TO_VSTRING(NODES%NODES(node_idx)%USER_NUMBER,"*",ERR,ERROR))// &
                & " has inconsistent derivative directions."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ENDDO !node_idx
        ELSE
          CALL FLAG_ERROR("Mesh topology nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*999)
    ENDIF
   
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",NODES%NUMBER_OF_NODES,ERR,ERROR,*999)
      DO node_idx=1,NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES, &
          & ERR,ERROR,*999)
  !#######################################################VERSIONS START########################################################
        DO derivative_idx=1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
          !TODO: change output string below so that it writes out derivative_idx index as well
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivative_idx) = ", &
            & NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivative_idx) = ", &
            & NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
        ENDDO !derivative_idx
  !#######################################################VERSIONS END##########################################################
        !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES,8,8, & !TODO: VERSIONS Delete line
        !  & NODES%NODES(node_idx)%GLOBAL_DERIVATIVE_INDEX,'("    Global derivative index(derivative_idx) :",8(X,I2))', &
        !  & '(36X,8(X,I2))',ERR,ERROR,*999) !TODO: VERSIONS Delete line
        !CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES,8,8, & !TODO: VERSIONS Delete line
        !  & NODES%NODES(node_idx)%PARTIAL_DERIVATIVE_INDEX,'("    Partial derivative index(derivative_idx) :",8(X,I2))', &
        !  & '(36X,8(X,I2))',ERR,ERROR,*999) !TODO: VERSIONS Delete line
      ENDDO !node_idx
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE")
    RETURN
999 IF(ALLOCATED(DERIVATIVES)) DEALLOCATE(DERIVATIVES)
    IF(ASSOCIATED(NODE_DERIVATIVE_LIST)) CALL LIST_DESTROY(NODE_DERIVATIVE_LIST,ERR,ERROR,*998)
998 CALL ERRORS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE

  !#######################################################VERSIONS END##########################################################

  !
  !================================================================================================================================
  !


  !>Calculates the number of versions at each node in a topology.
  SUBROUTINE MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the versions at each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_VERSIONS,ne,nn,derivative_idx,node_idx
    INTEGER(INTG), ALLOCATABLE :: VERSIONS(:)
    TYPE(LIST_PTR_TYPE), POINTER :: NODE_VERSION_LIST(:,:)
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(MESH_NODES_TYPE), POINTER :: NODES
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE",ERR,ERROR,*999)

     IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          NODES=>TOPOLOGY%NODES
          !Loop over the mesh elements
          !Calculate the number of versions at each node. This needs to be calculated by looking at all the mesh elements
          !as we may have an adjacent elements in another domain with a higher order basis along with different versions being assigned to its derivatives.
          !TODO: See if there are any constraints that can be applied to restrict the amount of memory being allocated here
          ALLOCATE(NODE_VERSION_LIST(MAXIMUM_GLOBAL_DERIV_NUMBER,NODES%NUMBER_OF_NODES),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node version list.",ERR,ERROR,*999)
          DO node_idx=1,NODES%NUMBER_OF_NODES
            DO derivative_idx=1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
              NULLIFY(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR)
              CALL LIST_CREATE_START(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,8,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,ERR,ERROR,*999)
            ENDDO!derivative_idx
          ENDDO!node_idx
          DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
            BASIS=>ELEMENTS%ELEMENTS(ne)%BASIS
            DO nn=1,BASIS%NUMBER_OF_NODES
              DO derivative_idx=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                CALL LIST_ITEM_ADD(NODE_VERSION_LIST(derivative_idx,ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn))%PTR, & 
                  & ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODE_VERSIONS(derivative_idx,nn),ERR,ERROR,*999)
              ENDDO!derivative_idx
            ENDDO!nn
          ENDDO!ne
          DO node_idx=1,NODES%NUMBER_OF_NODES
            DO derivative_idx=1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
              CALL LIST_REMOVE_DUPLICATES(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,ERR,ERROR,*999)
              CALL LIST_DETACH_AND_DESTROY(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,NUMBER_OF_VERSIONS,VERSIONS, &
                & ERR,ERROR,*999)
              NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS = NUMBER_OF_VERSIONS
              ALLOCATE(NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%USER_VERSION_NUMBERS(NUMBER_OF_VERSIONS),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global derivative index.",ERR,ERROR,*999)
              !Possible bug in LIST_REMOVE_DUPLICATES, when it comes time to LIST_DETACH_AND_DESTROY the deleted list values are still present. Or is this how lists are supposed to work - you shift the deleted one to the end of the list and just change the number of list values to -1 
              NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%USER_VERSION_NUMBERS = VERSIONS(1:NUMBER_OF_VERSIONS)
              DEALLOCATE(VERSIONS)
            ENDDO!derivative_idx
          ENDDO!node_idx
        ELSE
          CALL FLAG_ERROR("Mesh topology nodes is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not associated.",ERR,ERROR,*999)
    ENDIF
   
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",NODES%NUMBER_OF_NODES,ERR,ERROR,*999)
      DO node_idx=1,NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",node_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ", &
          & NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES,ERR,ERROR,*999)
        DO derivative_idx=1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
          !TODO: change output string below so that it writes out derivative_idx index as well
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global derivative index(derivative_idx) = ", &
            & NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%GLOBAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Partial derivative index(derivative_idx) = ", &
            & NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%PARTIAL_DERIVATIVE_INDEX,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1, &
            & NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%NUMBER_OF_VERSIONS,8,8, &
            & NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%USER_VERSION_NUMBERS, & 
            & '("    User Version index(derivative_idx,:) :",8(X,I2))','(36X,8(X,I2))',ERR,ERROR,*999)
        ENDDO!derivative_idx
      ENDDO !node_idx
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE")
    RETURN
999 IF(ALLOCATED(VERSIONS)) DEALLOCATE(VERSIONS)
    DO node_idx=1,NODES%NUMBER_OF_NODES
      DO derivative_idx=1,NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
        IF(ASSOCIATED(NODE_VERSION_LIST)) THEN
          CALL LIST_DESTROY(NODE_VERSION_LIST(derivative_idx,node_idx)%PTR,ERR,ERROR,*998)
        ENDIF
      ENDDO !derivative_idx
    ENDDO !node_idx
998 CALL ERRORS("MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_NODES_VERSIONS_CALCULATE

  !#######################################################VERSIONS END##########################################################

  !
  !================================================================================================================================
  !


  !>Calculates the element numbers surrounding a node for a mesh.
  SUBROUTINE MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the elements surrounding each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: element_no,insert_position,ne,nn,np,surrounding_elem_no
    INTEGER(INTG), POINTER :: NEW_SURROUNDING_ELEMENTS(:)
    LOGICAL :: FOUND_ELEMENT
    TYPE(BASIS_TYPE), POINTER :: BASIS

    NULLIFY(NEW_SURROUNDING_ELEMENTS)

    CALL ENTERS("MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) THEN
            DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
              BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
              DO nn=1,BASIS%NUMBER_OF_NODES
                np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)
                FOUND_ELEMENT=.FALSE.
                element_no=1
                insert_position=1
                DO WHILE(element_no<=TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.FOUND_ELEMENT)
                  surrounding_elem_no=TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(element_no)
                  IF(surrounding_elem_no==ne) THEN
                    FOUND_ELEMENT=.TRUE.
                  ENDIF
                  element_no=element_no+1
                  IF(ne>=surrounding_elem_no) THEN
                    insert_position=element_no
                  ENDIF
                ENDDO
                IF(.NOT.FOUND_ELEMENT) THEN
                  !Insert element into surrounding elements
                  ALLOCATE(NEW_SURROUNDING_ELEMENTS(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new surrounding elements",ERR,ERROR,*999)
                  IF(ASSOCIATED(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)) THEN
                    NEW_SURROUNDING_ELEMENTS(1:insert_position-1)=TOPOLOGY%NODES%NODES(np)% &
                      & SURROUNDING_ELEMENTS(1:insert_position-1)
                    NEW_SURROUNDING_ELEMENTS(insert_position)=ne
                    NEW_SURROUNDING_ELEMENTS(insert_position+1:TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1)= &
                      & TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(insert_position:TOPOLOGY%NODES%NODES(np)% &
                      & NUMBER_OF_SURROUNDING_ELEMENTS)
                    DEALLOCATE(TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS)
                  ELSE
                    NEW_SURROUNDING_ELEMENTS(1)=ne
                  ENDIF
                  TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS=>NEW_SURROUNDING_ELEMENTS
                  TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS= &
                    & TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS+1
                ENDIF
              ENDDO !nn
            ENDDO !ne
          ELSE
            CALL FLAG_ERROR("Mesh topology nodes nodes are not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh topology nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE")
    RETURN
999 IF(ASSOCIATED(NEW_SURROUNDING_ELEMENTS)) DEALLOCATE(NEW_SURROUNDING_ELEMENTS)
    CALL ERRORS("MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE")
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE
  
  !
  !===============================================================================================================================
  !

  !>Finalises the nodes data structures for a mesh topology and deallocates any memory. \todo pass in nodes
  SUBROUTINE MESH_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to finalise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: np

    CALL ENTERS("MESH_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%NUMBER_OF_NODES
          CALL MESH_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        DEALLOCATE(TOPOLOGY%NODES%NODES)
        IF(ASSOCIATED(TOPOLOGY%NODES%NODES_TREE)) CALL TREE_DESTROY(TOPOLOGY%NODES%NODES_TREE,ERR,ERROR,*999)
        DEALLOCATE(TOPOLOGY%NODES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("MESH_TOPOLOGY_NODES_FINALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODES_FINALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_FINALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the nodes in a given mesh topology. \todo finalise on errors
  SUBROUTINE MESH_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to initialise the nodes for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        CALL FLAG_ERROR("Mesh already has topology nodes associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%NODES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology nodes",ERR,ERROR,*999)
        TOPOLOGY%NODES%NUMBER_OF_NODES=0
        TOPOLOGY%NODES%MESH=>TOPOLOGY%MESH
        NULLIFY(TOPOLOGY%NODES%NODES)
        NULLIFY(TOPOLOGY%NODES%NODES_TREE)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODES_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODES_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODES_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given list of MESHES. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,MESHES,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(MESHES_TYPE), POINTER :: MESHES !<The list of meshes containing the mesh.
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: mesh_idx

    CALL ENTERS("MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      IF(ASSOCIATED(MESH)) THEN
        CALL FLAG_ERROR("Mesh is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(MESH)
        mesh_idx=1
        DO WHILE(mesh_idx<=MESHES%NUMBER_OF_MESHES.AND..NOT.ASSOCIATED(MESH))
          IF(MESHES%MESHES(mesh_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            MESH=>MESHES%MESHES(mesh_idx)%PTR
          ELSE
            mesh_idx=mesh_idx+1
          ENDIF
        ENDDO
      ENDIF
    ELSE
      CALL FLAG_ERROR("Meshes is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_USER_NUMBER_FIND_GENERIC")
    RETURN
999 CALL ERRORS("MESH_USER_NUMBER_FIND_GENERIC",ERR,ERROR)
    CALL EXITS("MESH_USER_NUMBER_FIND_GENERIC")
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND_GENERIC

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given INTERFACE. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_INTERFACE(USER_NUMBER,INTERFACE,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
  
    CALL ENTERS("MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,INTERFACE%MESHES,MESH,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN
999 CALL ERRORS("MESH_USER_NUMBER_FIND_INTERFACE",ERR,ERROR)
    CALL EXITS("MESH_USER_NUMBER_FIND_INTERFACE")
    RETURN 1
    
  END SUBROUTINE MESH_USER_NUMBER_FIND_INTERFACE

  !
  !================================================================================================================================
  !

  !>Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given REGION. If no mesh with that number exits MESH is left nullified.
  SUBROUTINE MESH_USER_NUMBER_FIND_REGION(USER_NUMBER,REGION,MESH,ERR,ERROR,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the mesh to find
    TYPE(REGION_TYPE), POINTER :: REGION !<The region containing the mesh
    TYPE(MESH_TYPE), POINTER :: MESH !<On return, a pointer to the mesh of the specified user number. In no mesh with the specified user number exists the pointer is returned NULL.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    
    CALL ENTERS("MESH_USER_NUMBER_FIND_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      CALL MESH_USER_NUMBER_FIND_GENERIC(USER_NUMBER,REGION%MESHES,MESH,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_USER_NUMBER_FIND_REGION")
    RETURN
999 CALL ERRORS("MESH_USER_NUMBER_FIND_REGION",ERR,ERROR)
    CALL EXITS("MESH_USER_NUMBER_FIND_REGION")
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND_REGION

  !
  !================================================================================================================================
  !

  !>Finalises the meshes and deallocates all memory
  SUBROUTINE MESHES_FINALISE(MESHES,ERR,ERROR,*)

   !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to finalise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH
 
    CALL ENTERS("MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESHES)) THEN
      DO WHILE(MESHES%NUMBER_OF_MESHES>0)
        MESH=>MESHES%MESHES(1)%PTR
        CALL MESH_DESTROY(MESH,ERR,ERROR,*999)
      ENDDO !mesh_idx
      DEALLOCATE(MESHES)
    ELSE
      CALL FLAG_ERROR("Meshes is not associated.",ERR,ERROR,*999)
    ENDIF
 
    CALL EXITS("MESHES_FINALISE")
    RETURN
999 CALL ERRORS("MESHES_FINALISE",ERR,ERROR)
    CALL EXITS("MESHES_FINALISE")
    RETURN 1
   
  END SUBROUTINE MESHES_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the generic meshes.
  SUBROUTINE MESHES_INITIALISE_GENERIC(MESHES,ERR,ERROR,*)

    !Argument variables
    TYPE(MESHES_TYPE), POINTER :: MESHES !<A pointer to the meshes to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("MESHES_INITIALISE_GENERIC",ERR,ERROR,*998)

    IF(ASSOCIATED(MESHES)) THEN
      CALL FLAG_ERROR("Meshes is already associated.",ERR,ERROR,*998)
    ELSE
      ALLOCATE(MESHES,STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Meshes could not be allocated",ERR,ERROR,*999)
      NULLIFY(MESHES%REGION)
      NULLIFY(MESHES%INTERFACE)
      MESHES%NUMBER_OF_MESHES=0
      NULLIFY(MESHES%MESHES)
    ENDIF
    
    CALL EXITS("MESHES_INITIALISE_GENERIC")
    RETURN
999 CALL MESHES_FINALISE(MESHES,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("MESHES_INITIALISE_GENERIC",ERR,ERROR)
    CALL EXITS("MESHES_INITIALISE_GENERIC")
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_GENERIC

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given interface.
  SUBROUTINE MESHES_INITIALISE_INTERFACE(INTERFACE,ERR,ERROR,*)

    !Argument variables
    TYPE(INTERFACE_TYPE), POINTER :: INTERFACE !<A pointer to the interface to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESHES_INITIALISE_INTERFACE",ERR,ERROR,*999)

    IF(ASSOCIATED(INTERFACE)) THEN
      IF(ASSOCIATED(INTERFACE%MESHES)) THEN
        LOCAL_ERROR="Interface number "//TRIM(NUMBER_TO_VSTRING(INTERFACE%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(INTERFACE%MESHES,ERR,ERROR,*999)
        INTERFACE%MESHES%INTERFACE=>INTERFACE
      ENDIF
    ELSE
      CALL FLAG_ERROR("Interface is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESHES_INITIALISE_INTERFACE")
    RETURN
999 CALL ERRORS("MESHES_INITIALISE_INTERFACE",ERR,ERROR)
    CALL EXITS("MESHES_INITIALISE_INTERFACE")
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_INTERFACE

  !
  !================================================================================================================================
  !

  !>Initialises the meshes for the given region.
  SUBROUTINE MESHES_INITIALISE_REGION(REGION,ERR,ERROR,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the meshes for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESHES_INITIALISE_REGION",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        CALL MESHES_INITIALISE_GENERIC(REGION%MESHES,ERR,ERROR,*999)
        REGION%MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESHES_INITIALISE_REGION")
    RETURN
999 CALL ERRORS("MESHES_INITIALISE_REGION",ERR,ERROR)
    CALL EXITS("MESHES_INITIALISE_REGION")
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE_REGION

  !
  !================================================================================================================================
  !

    !>Gets the domain for a given node in a decomposition of a mesh. \todo should be able to specify lists of elements. \see OPENCMISS::CMISSDecompositionNodeDomainGet
  SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET(DECOMPOSITION,USER_NODE_NUMBER,MESH_COMPONENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION !<A pointer to the decomposition to set the element domain for
    INTEGER(INTG), INTENT(IN) :: USER_NODE_NUMBER !<The global element number to set the domain for.
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER !<The mesh component number to set the domain for.
    INTEGER(INTG), INTENT(OUT) :: DOMAIN_NUMBER !<On return, the domain of the global element.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables`
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: GLOBAL_NODE_NUMBER
    TYPE(TREE_NODE_TYPE), POINTER :: TREE_NODE
    TYPE(MESH_NODES_TYPE), POINTER :: MESH_NODES
    TYPE(DOMAIN_TYPE), POINTER :: MESH_DOMAIN

    CALL ENTERS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    GLOBAL_NODE_NUMBER=0
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH)) THEN
          MESH_TOPOLOGY=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR
          IF(ASSOCIATED(MESH_TOPOLOGY)) THEN
            MESH_NODES=>MESH_TOPOLOGY%NODES
            IF(ASSOCIATED(MESH_NODES)) THEN
              NULLIFY(TREE_NODE)
              CALL TREE_SEARCH(MESH_NODES%NODES_TREE,USER_NODE_NUMBER,TREE_NODE,ERR,ERROR,*999)
              IF(ASSOCIATED(TREE_NODE)) THEN
                CALL TREE_NODE_VALUE_GET(MESH_NODES%NODES_TREE,TREE_NODE,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                IF(GLOBAL_NODE_NUMBER>0.AND.GLOBAL_NODE_NUMBER<=MESH_TOPOLOGY%NODES%NUMBER_OF_NODES) THEN
                  IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
                    MESH_DOMAIN=>DECOMPOSITION%DOMAIN(MESH_COMPONENT_NUMBER)%PTR
                    IF(ASSOCIATED(MESH_DOMAIN)) THEN
                      DOMAIN_NUMBER=MESH_DOMAIN%NODE_DOMAIN(GLOBAL_NODE_NUMBER)
                    ELSE
                      CALL FLAG_ERROR("Decomposition domain is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    LOCAL_ERROR="Mesh Component number "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
                      & " is invalid. The limits are 1 to "// &
                      & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  LOCAL_ERROR="Global element number found "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NODE_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid. The limits are 1 to "// &
                    & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES,"*",ERR,ERROR))
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Decomposition mesh node corresponding to user number not found",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Decomposition mesh nodes are not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Decomposition mesh topology is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition has not been finished",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("DECOMPOSITION_NODE_DOMAIN_GET")
    RETURN
999 CALL ERRORS("DECOMPOSITION_NODE_DOMAIN_GET",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_NODE_DOMAIN_GET")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_NODE_DOMAIN_GET

!
  !================================================================================================================================
  !

  !>Initialises the embedded meshes.
  SUBROUTINE EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*)

    !Argument variables
    !TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to initialise
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("EMBEDDED_MESH_INITIALISE",ERR,ERROR,*998)
    
    ALLOCATE(MESH_EMBEDDING,STAT=ERR)
    NULLIFY(MESH_EMBEDDING%PARENT_MESH)
    NULLIFY(MESH_EMBEDDING%CHILD_MESH)
    
    CALL EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN
!999 CALL EMBEDDED_MESH_FINALISE(MESH_EMBEDDING,DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("EMBEDDED_MESH_INITIALISE",ERR,ERROR)
    CALL EXITS("EMBEDDED_MESH_INITIALISE")
    RETURN 1
  END SUBROUTINE EMBEDDED_MESH_INITIALISE

  !
  !
  !================================================================================================================================
  !

  !>Creates an embedding of one mesh in another
  SUBROUTINE MESH_EMBEDDING_CREATE(MESH_EMBEDDING, PARENT_MESH, CHILD_MESH,ERR,ERROR,*)
!    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_EMBEDDING_TYPE), POINTER :: MESH_EMBEDDING !<Mesh embedding to create.
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: PARENT_MESH !<The parent mesh
    TYPE(MESH_TYPE), POINTER, INTENT(IN) :: CHILD_MESH  !<The child mesh
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: NGP = 0, ne

    CALL ENTERS("MESH_EMBEDDING_CREATE",ERR,ERROR,*999) 
    
    WRITE(*,*) 'parent mesh', PARENT_MESH%NUMBER_OF_ELEMENTS
    WRITE(*,*) 'child mesh', child_MESH%NUMBER_OF_ELEMENTS
    CALL EMBEDDED_MESH_INITIALISE(MESH_EMBEDDING,ERR,ERROR,*999)

    DO ne=1,PARENT_MESH%NUMBER_OF_ELEMENTS
      NGP = MAX(NGP,PARENT_MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS%QUADRATURE%&
        & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR%NUMBER_OF_GAUSS)
    ENDDO !ne

    MESH_EMBEDDING%PARENT_MESH => PARENT_MESH
    MESH_EMBEDDING%CHILD_MESH  => CHILD_MESH
    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate child node positions.",ERR,ERROR,*999)
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(NGP,PARENT_MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate gauss point positions.",ERR,ERROR,*999)

999 CALL ERRORS("MESH_EMBEDDING_CREATE",ERR,ERROR)
    CALL EXITS("MESH_EMBEDDING_CREATE")
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_CREATE

  !
  !================================================================================================================================
  !

  !>Sets the positions of nodes in the child mesh for one element in the parent mesh
  SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION(MESH_EMBEDDING, ELEMENT_NUMBER, NODE_NUMBERS, XI_COORDS,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER  !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: NODE_NUMBERS(:) !<NODE_NUMBERS(node_idx) Node numbers in child mesh for the node_idx'th embedded node in the ELEMENT_NUMBER'th element of the parent mesh
    REAL(DP), INTENT(IN)      :: XI_COORDS(:,:)  !<XI_COORDS(:,node_idx) Xi coordinates of the node_idx'th embedded node in the ELEMENT_NUMBER'th

    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    CALL ENTERS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR,*999)

    IF(ELEMENT_NUMBER<1 .OR. ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FLAG_ERROR("Element number out of range",ERR,ERROR,*999)
    ENDIF

    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NUMBER_OF_NODES = SIZE(NODE_NUMBERS)

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS(SIZE(NODE_NUMBERS)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%NODE_NUMBERS = NODE_NUMBERS

    ALLOCATE(MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS(SIZE(XI_COORDS,1),SIZE(XI_COORDS,2)))
    MESH_EMBEDDING%CHILD_NODE_XI_POSITION(ELEMENT_NUMBER)%XI_COORDS = XI_COORDS
    RETURN
999 CALL ERRORS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION",ERR,ERROR)
    CALL EXITS("MESH_EMBEDDING_SET_CHILD_NODE_POSITION")
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_CHILD_NODE_POSITION

  !
  !================================================================================================================================
  !

  !>Sets the positions of a Gauss point of the parent mesh in terms of element/xi coordinate in the child mesh
  SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA(MESH_EMBEDDING, PARENT_ELEMENT_NUMBER, GAUSSPT_NUMBER,&
    & PARENT_XI_COORD, CHILD_ELEMENT_NUMBER, CHILD_XI_COORD,ERR,ERROR,*)
    TYPE(MESH_EMBEDDING_TYPE), INTENT(INOUT) :: MESH_EMBEDDING   !<The mesh embedding object
    INTEGER(INTG), INTENT(IN) :: PARENT_ELEMENT_NUMBER           !<Element number in the parent mesh
    INTEGER(INTG), INTENT(IN) :: GAUSSPT_NUMBER                  !<Gauss point number in this element
    REAL(DP), INTENT(IN) :: PARENT_XI_COORD(:)              !<Xi coordinate in parent element

    INTEGER(INTG), INTENT(IN) :: CHILD_ELEMENT_NUMBER !<Element number in the child mesh
    REAL(DP), INTENT(IN) :: CHILD_XI_COORD(:)    !<Xi coordinate in child element

    INTEGER(INTG), INTENT(OUT) :: ERR           !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR  !<The error string

    CALL ENTERS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR,*999)

    IF(PARENT_ELEMENT_NUMBER<1 .OR. PARENT_ELEMENT_NUMBER > MESH_EMBEDDING%PARENT_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FLAG_ERROR("Parent element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(CHILD_ELEMENT_NUMBER<1 .OR. CHILD_ELEMENT_NUMBER > MESH_EMBEDDING%CHILD_MESH%NUMBER_OF_ELEMENTS) THEN
      CALL FLAG_ERROR("Child element number out of range",ERR,ERROR,*999)
    ENDIF
    IF(GAUSSPT_NUMBER<1 .OR. GAUSSPT_NUMBER > SIZE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION,1)) THEN
      CALL FLAG_ERROR("Gauss point number out of range",ERR,ERROR,*999)
    ENDIF

    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %PARENT_XI_COORD(SIZE(PARENT_XI_COORD)))
    ALLOCATE(MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)&
     & %CHILD_XI_COORD(SIZE(CHILD_XI_COORD)))


    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%PARENT_XI_COORD = PARENT_XI_COORD
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%CHILD_XI_COORD = CHILD_XI_COORD
    MESH_EMBEDDING%GAUSS_POINT_XI_POSITION(GAUSSPT_NUMBER,PARENT_ELEMENT_NUMBER)%ELEMENT_NUMBER = CHILD_ELEMENT_NUMBER

    RETURN
999 CALL ERRORS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA",ERR,ERROR)
    CALL EXITS("MESH_EMBEDDING_SET_GAUSS_POINT_DATA")
    RETURN 1
  END SUBROUTINE MESH_EMBEDDING_SET_GAUSS_POINT_DATA


!   !
!   !================================================================================================================================
!   !

END MODULE MESH_ROUTINES

