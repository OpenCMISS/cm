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

  PUBLIC MESH_TOPOLOGY_ELEMENT_CHECK_EXISTS,MESH_TOPOLOGY_NODE_CHECK_EXISTS
  
  PUBLIC MESHES_INITIALISE,MESHES_FINALISE

  PUBLIC MESH_CREATE_START,MESH_CREATE_FINISH

  PUBLIC MESH_DESTROY
  
  PUBLIC MESH_NUMBER_OF_COMPONENTS_GET,MESH_NUMBER_OF_COMPONENTS_SET

  PUBLIC MESH_NUMBER_OF_ELEMENTS_GET,MESH_NUMBER_OF_ELEMENTS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_CREATE_START,MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH

  PUBLIC MESH_TOPOLOGY_ELEMENTS_DESTROY

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_GET

  PUBLIC MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_GET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_USER_NUMBER_SET

  PUBLIC MESH_USER_NUMBER_FIND
  
  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE

  PUBLIC DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH

  PUBLIC DECOMPOSITION_DESTROY

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE

  PUBLIC DECOMPOSITION_ELEMENT_DOMAIN_GET,DECOMPOSITION_ELEMENT_DOMAIN_SET

  PUBLIC DECOMPOSITION_MESH_COMPONENT_NUMBER_GET,DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
  
  PUBLIC DECOMPOSITION_NUMBER_OF_DOMAINS_GET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET
  
  PUBLIC DECOMPOSITION_TYPE_GET,DECOMPOSITION_TYPE_SET
  
  PUBLIC DECOMPOSITION_USER_NUMBER_FIND
  
  PUBLIC DECOMPOSITION_NODE_DOMAIN_GET
  
CONTAINS

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
      CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
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
              NUMBER_FLAG=0 !C Numbering as there is a bug with Fortran numbering
              NUMBER_OF_CONSTRAINTS=1
              NUMBER_OF_COMMON_NODES=2
              TPWGTS=1.0_SP/REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,SP)
              UBVEC=1.05_SP
              PARMETIS_OPTIONS(0)=1
              PARMETIS_OPTIONS(1)=7
              
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
        CASE(DECOMPOSITION_CALCULATED_TYPE)
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
  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the decomposition topology to calculate.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the face topology
!!      CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
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

  !>Finalises the given decomposition topology element.
  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT !<The decomposition element to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS)) DEALLOCATE(ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS)
    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
    IF(ALLOCATED(ELEMENT%ELEMENT_LINES)) DEALLOCATE(ELEMENT%ELEMENT_LINES)
 
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
    INTEGER(INTG) :: i,j,ne1,ne2,nep1,nep2,ni1,nn,nn1,nn2,nnl,nng,nnpl,np,np1,np2,DUMMY_ERR,FACE_XI(2)    
    INTEGER(INTG) :: direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    
    INTEGER(INTG) :: ne,ni,is,COUNTER1,COUNTER2,COUNTER3,NUMBER_OF_ITEM,GLOBAL_NUMBER_OF_NODE,XI_DIRECTION
    INTEGER(INTG) :: NODE_NUMBER_IN_LOCAL_INTERFACE,NUMBER_OF_LOCAL_INTERFACE,NUMBER_OF_NODES_IN_LOCAL_INTERFACE
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS
    INTEGER(INTG) :: NUMBER_SURROUNDING,MAX_NUMBER_SURROUNDING,NUMBER_OF_NODES_XI(3),NODE_POSITION_INDEX(3)
    INTEGER(INTG), POINTER :: ADJACENT_ELEMENTS(:),SURROUNDING_ELEMENTS(:),ELEMENTS_MATCH(:)
    LOGICAL :: FOUND,XI_COLLAPSED,FACE_COLLAPSED(-3:3)
    TYPE(VARYING_STRING) :: DUMMY_ERROR
!    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(LIST_TYPE), POINTER :: ELEMENTS_MATCH_LIST ! JUST CONTAINS LIST OF SURROUNDING ELEMENTS AND NO MORE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: SURROUNDING_ELEMENTS_LIST(:) 
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY

    NULLIFY(SURROUNDING_ELEMENTS)
    NULLIFY(ADJACENT_ELEMENTS)
    NULLIFY(ELEMENTS_MATCH)
    NULLIFY(BASIS)
    
    CALL ENTERS("DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE",ERR,ERROR,*999)

    !%%% Here is the algorithm

!1  DO i=1,total_number_of_elements
!2    DO node=1,number_of_nodes_per_this_element
!3      ADD (surrounding element number to ELEMENTS_MATCH_LIST) and (local node number to SURROUNDING_ELEMENTS_LIST(i)%ptr)
!4    CHECK data in SURROUNDING_ELEMENTS_LIST(i)%ptr with faces of element to find the an ELEMENTS_MATCH_LIST(i) that is adjacent in face. 
!5    FIND the corresponding xi direction for this face
!6    ADD ELEMENTS_MATCH_LIST(i) to ADJACENT_ELEMENT and the corresponding xi data in TOPOLOGY%ELEMENTS%ELEMENTS dataset

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
                  ALLOCATE(SURROUNDING_ELEMENTS_LIST(DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS),STAT=ERR)
                  !Loop over the elements in the decomposition
                  DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                  !%%%% first we initialize lists that are required to find the adjacent elements list
                    BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
                    DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      NULLIFY(ADJACENT_ELEMENTS_LIST(ni)%PTR)
                      CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,5,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
                    ENDDO !ni
                                                            
                    NULLIFY(ELEMENTS_MATCH_LIST) ! Contains element number of the surrounding elements 
                    CALL LIST_CREATE_START(ELEMENTS_MATCH_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ELEMENTS_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ELEMENTS_MATCH_LIST,16,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ELEMENTS_MATCH_LIST,ERR,ERROR,*999)
                    DO COUNTER1=1,ELEMENTS_MATCH_LIST%INITIAL_SIZE ! SET INITIAL VALUES TO ZERO
                      ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1)=0
                    ENDDO              
                    
                    NUMBER_OF_SURROUNDING_ELEMENTS=0
                    !%%% now loop over the local nodes                    
                    DO COUNTER1=1,DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES
                      GLOBAL_NUMBER_OF_NODE=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(COUNTER1)  
                      !%%% loop over the surrounding elements to the selected node
                      DO COUNTER2=1,DOMAIN_NODES%NODES(GLOBAL_NUMBER_OF_NODE)%NUMBER_OF_SURROUNDING_ELEMENTS
                        CALL LIST_SEARCH(ELEMENTS_MATCH_LIST%LIST_INTG,DOMAIN_NODES%NODES(GLOBAL_NUMBER_OF_NODE) &
                         & %SURROUNDING_ELEMENTS(COUNTER2),NUMBER_OF_ITEM,ERR,ERROR,*999)
                        IF (NUMBER_OF_ITEM .EQ. 0) THEN
                          NUMBER_OF_SURROUNDING_ELEMENTS=NUMBER_OF_SURROUNDING_ELEMENTS+1
                          
                          NULLIFY(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR)
                          CALL LIST_CREATE_START(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,ERR,ERROR,*999)
                          CALL LIST_DATA_TYPE_SET(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,LIST_INTG_TYPE, &
                           & ERR,ERROR,*999)
                          CALL LIST_INITIAL_SIZE_SET(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,5,ERR,ERROR,*999)
                          CALL LIST_CREATE_FINISH(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,ERR,ERROR,*999)

                          DO COUNTER3=1,SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR%INITIAL_SIZE ! SET INITIAL VALUES ZERO
                            SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR%LIST_INTG(COUNTER3)=0
                          ENDDO                          

                          CALL LIST_ITEM_ADD(ELEMENTS_MATCH_LIST,DOMAIN_NODES%NODES(GLOBAL_NUMBER_OF_NODE) &
                           & %SURROUNDING_ELEMENTS(COUNTER2),ERR,ERROR,*999)
                          CALL LIST_ITEM_ADD(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,COUNTER1,ERR,ERROR,*999)
                        ELSE
                          CALL LIST_ITEM_ADD(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_ITEM)%PTR,COUNTER1,ERR,ERROR,*999)
                        ENDIF
                      ENDDO
                    ENDDO
                    
                    !%%% now surrounding elements and shared nodes are stored, we want to find the adjacent elements and correspounding xi direction
                    ! loading ADJACENT_ELEMENTS_LIST on xi=0  
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER, &
                     & ERR,ERROR,*999)
                                            
                    CALL LIST_NUMBER_OF_ITEMS_GET(ELEMENTS_MATCH_LIST,NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
                    DO COUNTER1=1,NUMBER_OF_SURROUNDING_ELEMENTS
                      CALL LIST_NUMBER_OF_ITEMS_GET(SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR,NUMBER_OF_ITEM,ERR,ERROR,*999)
                      IF (NUMBER_OF_ITEM .LT. BASIS%NUMBER_OF_NODES .AND. NUMBER_OF_ITEM .GT. 1) THEN ! if /= element ITSELF and not attached in 1 node
                      
                        ! first lets find the interface number
                        IF (BASIS%NUMBER_OF_XI .EQ. 2) NUMBER_OF_LOCAL_INTERFACE = BASIS%NUMBER_OF_LOCAL_LINES
                        IF (BASIS%NUMBER_OF_XI .EQ. 3) NUMBER_OF_LOCAL_INTERFACE = BASIS%NUMBER_OF_LOCAL_FACES
                        
                        DO COUNTER2=1,NUMBER_OF_LOCAL_INTERFACE
                          is = 1
                          IF (BASIS%NUMBER_OF_XI .EQ. 2) NUMBER_OF_NODES_IN_LOCAL_INTERFACE = &
                           & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(COUNTER2)
                          IF (BASIS%NUMBER_OF_XI .EQ. 3) NUMBER_OF_NODES_IN_LOCAL_INTERFACE = &
                           & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(COUNTER2)
                          DO COUNTER3=1,NUMBER_OF_NODES_IN_LOCAL_INTERFACE
                            IF (BASIS%NUMBER_OF_XI .EQ. 2) NODE_NUMBER_IN_LOCAL_INTERFACE = &
                             & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(COUNTER3,COUNTER2)
                            IF (BASIS%NUMBER_OF_XI .EQ. 3) NODE_NUMBER_IN_LOCAL_INTERFACE = & 
                             & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(COUNTER3,COUNTER2)
                            CALL LIST_SEARCH(SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR%LIST_INTG,NODE_NUMBER_IN_LOCAL_INTERFACE, &
                            & NUMBER_OF_ITEM,ERR,ERROR,*999)
                            IF (NUMBER_OF_ITEM .EQ. 0) is=0
                          ENDDO !COUNTER3
                          IF (is .EQ. 1) THEN ! (COUNTER2) is the interface number of the FACE/LINE that is shared between element DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER and ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1)
                            IF (BASIS%TYPE .EQ. 2) XI_DIRECTION = COUNTER2
                            IF (BASIS%TYPE .EQ. 1) THEN
                              IF (BASIS%NUMBER_OF_XI .EQ. 2) XI_DIRECTION = BASIS%LOCAL_LINE_XI_DIRECTION(COUNTER2)
                              IF (BASIS%NUMBER_OF_XI .EQ. 3) XI_DIRECTION = BASIS%LOCAL_FACE_XI_DIRECTION(COUNTER2)
                              XI_DIRECTION = XI_DIRECTION * (-1)**(COUNTER2+1)! FOR LAGRANGE TYPE SHOULD DETECT NEGATIVE DIRECTION
                            ENDIF
                            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(XI_DIRECTION)%PTR,ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1), &
                             & ERR,ERROR,*999)
                          ENDIF
                        ENDDO !COUNTER2
                                                    
                      ENDIF
                    ENDDO !COUNTER1
                    
                    !%%% set maximum number of adjacent elements
                    MAX_NUMBER_SURROUNDING=1
                    DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      CALL LIST_NUMBER_OF_ITEMS_GET(ADJACENT_ELEMENTS_LIST(ni)%PTR,i,ERR,ERROR,*999)
                      IF(i>MAX_NUMBER_SURROUNDING) MAX_NUMBER_SURROUNDING=i
                    ENDDO
                    
                    !%%% set the adjacent elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS( &
                      & -BASIS%NUMBER_OF_XI_COORDINATES:BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of surrounding elements.",ERR,ERROR,*999)
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(MAX_NUMBER_SURROUNDING, &
                      & -BASIS%NUMBER_OF_XI_COORDINATES:BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate surrounding elements.",ERR,ERROR,*999)
                    DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
                      CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & NUMBER_OF_ADJACENT_ELEMENTS(ni),ADJACENT_ELEMENTS,ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & NUMBER_OF_ADJACENT_ELEMENTS(ni),ni)=ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS% &
                        & ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni))
                      IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
                    ENDDO !ni
                    
                    ! destroy lists information
                    DO COUNTER1=1,NUMBER_OF_SURROUNDING_ELEMENTS
                      ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1)=0
                    ENDDO              
!                    CALL LIST_DETACH_AND_DESTROY(ELEMENTS_MATCH_LIST,NUMBER_OF_SURROUNDING_ELEMENTS,ELEMENTS_MATCH,ERR,ERROR,*999)
!                    DEALLOCATE(ELEMENTS_MATCH)

                    DO COUNTER1=1,NUMBER_OF_SURROUNDING_ELEMENTS
                      CALL LIST_NUMBER_OF_ITEMS_GET(SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR,NUMBER_OF_ITEM,ERR,ERROR,*999)
                      DO COUNTER2=1,NUMBER_OF_ITEM
                        SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR%LIST_INTG(COUNTER2)=0
                      ENDDO
                    ENDDO                          
!                    DO ni=1,NUMBER_OF_SURROUNDING_ELEMENTS
!                      CALL LIST_NUMBER_OF_ITEMS_GET(SURROUNDING_ELEMENTS_LIST(ni)%PTR,NUMBER_OF_ITEM,ERR,ERROR,*999)
!                      CALL LIST_DETACH_AND_DESTROY(SURROUNDING_ELEMENTS_LIST(ni)%PTR,NUMBER_OF_ITEM,SURROUNDING_ELEMENTS,ERR,ERROR,*999)
!                      DEALLOCATE(SURROUNDING_ELEMENTS)
!                    ENDDO !ni

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
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi directions = ", &
         & BASIS%NUMBER_OF_XI_COORDINATES,ERR,ERROR,*999)
        DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",ni,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni),ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
            & NUMBER_OF_ADJACENT_ELEMENTS(ni),8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(:,ni), &
            & '("        Adjacent elements =",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
        ENDDO !ni
      ENDDO !ne
    ENDIF

    CALL EXITS("DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE")
    RETURN
999 IF(ASSOCIATED(SURROUNDING_ELEMENTS)) DEALLOCATE(SURROUNDING_ELEMENTS)
    IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(ELEMENTS_MATCH_LIST)) CALL LIST_DESTROY(ELEMENTS_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*997)
!998 DO ni=-4,4
!      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(ni)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
!    ENDDO !ni
!    DO ni=1,100
!      IF(ASSOCIATED(SURROUNDING_ELEMENTS_LIST(ni)%PTR)) CALL LIST_DESTROY(SURROUNDING_ELEMENTS_LIST(ni)%PTR, &
!       & DUMMY_ERR,DUMMY_ERROR,*997)
!    ENDDO !ni
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
    INTEGER(INTG) :: ne,global_element
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
                          DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ne)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%USER_NUMBER=MESH_ELEMENTS%ELEMENTS(global_element)%USER_NUMBER
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER=ne
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=global_element
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%BOUNDARY_ELEMENT=MESH_ELEMENTS%ELEMENTS(global_element)% &
                              & BOUNDARY_ELEMENT
                          ENDDO !ne
                          !Calculate the elements surrounding the elements in the decomposition topology
                          !CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
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
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
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
      CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      CALL DECOMPOSITION_TOPOLOGY_FACES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
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
        CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        CALL DECOMPOSITION_TOPOLOGY_FACES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
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
    INTEGER(INTG) :: component_idx,ne,ne2,nae,nae2,nn,nnl,nk,nl,nl2,np,np2,elem_idx,node_idx,node_idx2,NODES_IN_LINE(4), &
      & NUMBER_OF_LINES,MAX_NUMBER_OF_LINES,NEW_MAX_NUMBER_OF_LINES,LINE_NUMBER,COUNT
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
                      DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_LINES(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element element lines.",ERR,ERROR,*999)
                        !Loop over the local lines of the element
                        DO nae=1,BASIS%NUMBER_OF_LOCAL_LINES
                          !Calculate the topology node numbers that make up the line
                          NODES_IN_LINE=0
                          DO nnl=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)
                            NODES_IN_LINE(nnl)=DOMAIN_ELEMENT%ELEMENT_NODES(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,nae))
                          ENDDO !nnl
                          !Try and find a previously created line that matches in the adjacent elements
                          FOUND=.FALSE.
                          np=NODES_IN_LINE(1)
                          DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                            ne2=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
                            IF(ne2/=ne) THEN
                              IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(ne2)%ELEMENT_LINES)) THEN
                                BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(ne2)%BASIS
                                DO nae2=1,BASIS2%NUMBER_OF_LOCAL_LINES
                                  nl=DECOMPOSITION_ELEMENTS%ELEMENTS(ne2)%ELEMENT_LINES(nae2)
                                  IF(ALL(NODES_IN_LINE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae))== &
                                    & TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae),nl))) THEN
                                    FOUND=.TRUE.
                                    EXIT
                                  ENDIF
                                ENDDO !nae2
                                IF(FOUND) EXIT
                              ENDIF
                            ENDIF
                          ENDDO !elem_idx
                          IF(FOUND) THEN
                            !Line has already been created
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(nae)=nl
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
                            DECOMPOSITION_ELEMENT%ELEMENT_LINES(nae)=NUMBER_OF_LINES
                            DO nnl=1,SIZE(NODES_IN_LINE,1)
                              IF(NODES_IN_LINE(nnl)/=0) &
                                & NODES_NUMBER_OF_LINES(NODES_IN_LINE(nnl))=NODES_NUMBER_OF_LINES(NODES_IN_LINE(nnl))+1
                            ENDDO !nnl
                          ENDIF
                        ENDDO !nae
                      ENDDO !ne
                      !Allocate the line arrays and set them from the LINES and NODE_LINES arrays
                      DO np=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                        ALLOCATE(DOMAIN_NODES%NODES(np)%NODE_LINES(NODES_NUMBER_OF_LINES(np)),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node lines array.",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES=0
                      ENDDO !np
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition topology lines.",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      ALLOCATE(DOMAIN_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain topology lines.",ERR,ERROR,*999)
                      DOMAIN_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      DO nl=1,DOMAIN_LINES%NUMBER_OF_LINES
                        CALL DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(DECOMPOSITION_LINES%LINES(nl),ERR,ERROR,*999)
                        CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(nl),ERR,ERROR,*999)
                        DO nnl=1,SIZE(TEMP_LINES,1)
                          IF(TEMP_LINES(nnl,nl)/=0) THEN
                            np=TEMP_LINES(nnl,nl)
                            DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES=DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES+1
                            DOMAIN_NODES%NODES(np)%NODE_LINES(DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES)=nl
                          ENDIF
                        ENDDO !nnl  
                      ENDDO !nl
                      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO nae=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(nae)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DOMAIN_LINE=>DOMAIN_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          IF(.NOT.ASSOCIATED(DOMAIN_LINE%BASIS)) THEN
                            DECOMPOSITION_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%NUMBER=LINE_NUMBER
                            DOMAIN_LINE%ELEMENT_NUMBER=ne !Needs checking
                            DECOMPOSITION_LINE%XI_DIRECTION=BASIS%LOCAL_LINE_XI_DIRECTION(nae)
                            DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                            ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line nodes in line.",ERR,ERROR,*999)
                            ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                              & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line derivatives in line.",ERR,ERROR,*999)
                            DOMAIN_LINE%NODES_IN_LINE=TEMP_LINES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae),LINE_NUMBER)
                            DO nnl=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)
                              DOMAIN_LINE%DERIVATIVES_IN_LINE(1,nnl)=NO_GLOBAL_DERIV
                              IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                nk=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nnl,nae), &
                                  & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,nae))
                                DOMAIN_LINE%DERIVATIVES_IN_LINE(2,nnl)=nk
                              ENDIF
                            ENDDO !nn
                          ENDIF
                        ENDDO !nae
                      ENDDO !ne
                      DEALLOCATE(TEMP_LINES)
                      !Calculate adjacent lines and the surrounding elements for each line
                      DO nl=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(nl)
                        DOMAIN_LINE=>DOMAIN_LINES%LINES(nl)
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
                        DO node_idx=0,1
                          FOUND=.FALSE.
                          np=DOMAIN_LINE%NODES_IN_LINE(node_idx*(BASIS%NUMBER_OF_NODES-1)+1)
                          !Loop over the elements surrounding the node.
                          DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                            ne=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
                            DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                            DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                            !Loop over the local lines of the element
                            DO nae=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_LINES
                              nl2=DECOMPOSITION_ELEMENT%ELEMENT_LINES(nae)
                              IF(nl2/=nl) THEN
                                DECOMPOSITION_LINE2=>DECOMPOSITION_LINES%LINES(nl2)
                                DOMAIN_LINE2=>DOMAIN_LINES%LINES(nl2)
                                IF(DECOMPOSITION_LINE2%XI_DIRECTION==DECOMPOSITION_LINE%XI_DIRECTION) THEN
                                  !Lines run in the same direction.
                                  BASIS2=>DOMAIN_LINE2%BASIS
                                  IF(node_idx==0) THEN
                                    np2=DOMAIN_LINE2%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES)
                                  ELSE
                                    np2=DOMAIN_LINE2%NODES_IN_LINE(1)
                                  ENDIF
                                  IF(np2==np) THEN
                                    !The node at the 'other' end of this line matches the node at the current end of the line.
                                    !Check it is not a coexistant line running the other way
                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
                                      COUNT=0
                                      DO node_idx2=1,BASIS%NUMBER_OF_NODES
                                        IF(DOMAIN_LINE2%NODES_IN_LINE(node_idx2)== &
                                          & DOMAIN_LINE%NODES_IN_LINE(BASIS2%NUMBER_OF_NODES-node_idx2+1)) &
                                          & COUNT=COUNT+1
                                      ENDDO !node_idx
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
                            ENDDO !nae
                            IF(FOUND) EXIT
                          ENDDO !element_idx
                          IF(FOUND) DECOMPOSITION_LINE%ADJACENT_LINES(node_idx)=nl2
                        ENDDO !node_idx
                      ENDDO !nl
                      !Set the surrounding elements
                      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        DO nae=1,BASIS%NUMBER_OF_LOCAL_LINES
                          LINE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_LINES(nae)
                          DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(LINE_NUMBER)
                          DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                          DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=ne
                          DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS)=nae
                        ENDDO !nae
                      ENDDO !ne
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
                            DO nl=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
                              DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(nl)
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(nl)
                              IF(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                ne=DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(1)
                                nae=DECOMPOSITION_LINE%ELEMENT_LINES(1)
                                CALL DOMAIN_TOPOLOGY_LINE_INITIALISE(DOMAIN_LINES%LINES(nl),ERR,ERROR,*999)
                                DOMAIN_LINE%NUMBER=nl
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_LINE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                                DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                                ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes in line.",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivatives in line.",ERR,ERROR,*999)
                                DO nnl=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)
                                  nn=BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nnl,nae)
                                  np=DOMAIN_ELEMENT%ELEMENT_NODES(nn)
                                  DOMAIN_LINE%NODES_IN_LINE(nnl)=np
                                  DOMAIN_LINE%DERIVATIVES_IN_LINE(1,nnl)=NO_GLOBAL_DERIV
                                  IF(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    nk=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nnl,nae),nn)
                                    DOMAIN_LINE%DERIVATIVES_IN_LINE(2,nnl)=nk
                                  ENDIF
                                  NODES_NUMBER_OF_LINES(np)=NODES_NUMBER_OF_LINES(np)+1
                                ENDDO !nnl
                              ELSE
                                CALL FLAG_ERROR("Line is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF                              
                            ENDDO !nl
                            DO np=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(np)%NODE_LINES(NODES_NUMBER_OF_LINES(np)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node lines.",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES=0
                            ENDDO !np
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO nl=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(nl)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO nnl=1,BASIS%NUMBER_OF_NODES
                                np=DOMAIN_LINE%NODES_IN_LINE(nnl)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(np)
                                DOMAIN_NODE%NUMBER_OF_NODE_LINES=DOMAIN_NODE%NUMBER_OF_NODE_LINES+1
                                DOMAIN_NODE%NODE_LINES(DOMAIN_NODE%NUMBER_OF_NODE_LINES)=nl
                              ENDDO !nnl
                            ENDDO !nl
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
      DO nl=1,DECOMPOSITION_LINES%NUMBER_OF_LINES
        DECOMPOSITION_LINE=>DECOMPOSITION_LINES%LINES(nl)
        DOMAIN_LINE=>DOMAIN_LINES%LINES(nl)
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
          DOMAIN_LINE=>DOMAIN%TOPOLOGY%LINES%LINES(nl)
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
          DO nnl=1,DOMAIN_LINE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",nnl,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_LINE%BASIS%NUMBER_OF_DERIVATIVES(nnl),4,4,DOMAIN_LINE% &
              & DERIVATIVES_IN_LINE(:,nnl),'("            Derivatives in line  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !nnl
        ENDDO !component_idx
      ENDDO !nl
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
    INTEGER(INTG) :: component_idx,ne,ne2,nae,nae2,nn,nnf,nk,nf,np,elem_idx,NODES_IN_FACE(16),NUMBER_OF_FACES, &
      & MAX_NUMBER_OF_FACES,NEW_MAX_NUMBER_OF_FACES,FACE_NUMBER!,NODE_COUNT,node_idx1,node_idx2,node_idx3,node_idx4,nf2,np2
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
                          DO nae=1,BASIS%NUMBER_OF_LOCAL_FACES
                            !Calculate the topology node numbers that make up the face
                            NODES_IN_FACE=0
                            !Check whether face has already been read out
                            DO nnf=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)
                              !Read out node numbers of local face from ELEMENT_NODES
                              NODES_IN_FACE(nnf)=DOMAIN_ELEMENT%ELEMENT_NODES(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(nnf,nae))
                            ENDDO !nnf
                            !Try and find a previously created face that matches in the adjacent elements
                            FOUND=.FALSE.
                            np=NODES_IN_FACE(1)
                            DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                              ne2=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
                              IF(ne2/=ne) THEN
                                IF(ALLOCATED(DECOMPOSITION_ELEMENTS%ELEMENTS(ne2)%ELEMENT_FACES)) THEN
                                  BASIS2=>DOMAIN_ELEMENTS%ELEMENTS(ne2)%BASIS
                                  DO nae2=1,BASIS2%NUMBER_OF_LOCAL_FACES
                                    nf=DECOMPOSITION_ELEMENTS%ELEMENTS(ne2)%ELEMENT_FACES(nae2)
                                    IF(ALL(NODES_IN_FACE(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae))== &
                                      & TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae),nf))) THEN
                                      FOUND=.TRUE.
                                      EXIT
                                    ENDIF
                                  ENDDO !nae2
                                  IF(FOUND) EXIT
                                ENDIF
                              ENDIF
                            ENDDO !elem_idx
                            IF(FOUND) THEN
                              !Face has already been created
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(nae)=nf
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
                              DECOMPOSITION_ELEMENT%ELEMENT_FACES(nae)=NUMBER_OF_FACES
                              DO nnf=1,SIZE(NODES_IN_FACE,1)
                                IF(NODES_IN_FACE(nnf)/=0) &
                                 & NODES_NUMBER_OF_FACES(NODES_IN_FACE(nnf))=NODES_NUMBER_OF_FACES(NODES_IN_FACE(nnf))+1
                              ENDDO !nnf
                            ENDIF
                          ENDDO !nae
                        ENDDO !ne

                        !Allocate the face arrays and set them from the FACES and NODE_FACES arrays
                        DO np=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                          ALLOCATE(DOMAIN_NODES%NODES(np)%NODE_FACES(NODES_NUMBER_OF_FACES(np)),STAT=ERR)
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node faces array",ERR,ERROR,*999)
                          DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_FACES=0
                        ENDDO !np
                        DEALLOCATE(NODES_NUMBER_OF_FACES)
                        ALLOCATE(DECOMPOSITION_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition topology faces",ERR,ERROR,*999)
                        DECOMPOSITION_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        ALLOCATE(DOMAIN_FACES%FACES(NUMBER_OF_FACES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain topology faces",ERR,ERROR,*999)
                        DOMAIN_FACES%NUMBER_OF_FACES=NUMBER_OF_FACES
                        DO nf=1,DOMAIN_FACES%NUMBER_OF_FACES
                          CALL DECOMPOSITION_TOPOLOGY_FACE_INITIALISE(DECOMPOSITION_FACES%FACES(nf),ERR,ERROR,*999)
                          CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(nf),ERR,ERROR,*999)
                          DO nnf=1,SIZE(TEMP_FACES,1)
                            IF(TEMP_FACES(nnf,nf)/=0) THEN
                              np=TEMP_FACES(nnf,nf)
                              DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_FACES=DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_FACES+1
                              DOMAIN_NODES%NODES(np)%NODE_FACES(DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_FACES)=nf
                            ENDIF
                          ENDDO !nnf  
                        ENDDO !nf

                        !Set nodes in face and derivatives of nodes in face for domain faces
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          !Loop over local faces of element
                          DO nae=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(nae)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DOMAIN_FACE=>DOMAIN_FACES%FACES(FACE_NUMBER)
                            DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS=DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS+1
                            IF(.NOT.ASSOCIATED(DOMAIN_FACE%BASIS)) THEN
                              DECOMPOSITION_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%NUMBER=FACE_NUMBER
                              DOMAIN_FACE%ELEMENT_NUMBER=ne !! Needs checking
!                              DECOMPOSITION_FACE%ELEMENT_NUMBER=DECOMPOSITION_ELEMENT%NUMBER
!                              DOMAIN_FACE%ELEMENT_NUMBER=DOMAIN_ELEMENT%NUMBER
                              DECOMPOSITION_FACE%XI_DIRECTION=BASIS%LOCAL_FACE_XI_DIRECTION(nae)
                              DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(DECOMPOSITION_FACE%XI_DIRECTION)%PTR
                              ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face nodes in face",ERR,ERROR,*999)
                              ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate face derivatives in face",ERR,ERROR,*999)
                              !Set nodes in face based upon face number
                              DOMAIN_FACE%NODES_IN_FACE=TEMP_FACES(1:BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae),FACE_NUMBER)
                              !Set derivatives of nodes in domain face from derivatives of nodes in element
                              DO nnf=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)
                                DOMAIN_FACE%DERIVATIVES_IN_FACE(1,nnf)=NO_GLOBAL_DERIV
                                IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                  nk=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(nnf,nae), &
                                    & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(nnf,nae))
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(2,nnf)=nk
                                ENDIF
                              ENDDO !nn
                            ENDIF
                          ENDDO !nae
                        ENDDO !ne

                        DEALLOCATE(TEMP_FACES)
! Note: Adjacency will be left out of faces calculation for the time being
                        !Calculate adjacent faces and the surrounding elements for each face
                        DO nf=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                          DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(nf)
                          DOMAIN_FACE=>DOMAIN_FACES%FACES(nf)
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
!                            np=DOMAIN_FACE%NODES_IN_FACE((node_idx2*BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION*(BASIS%NUMBER_OF_FACES-1))&
!                                                                             &+(node_idx1*(BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
                             !Loop over the elements surrounding the node.
!                            DO elem_idx=1,DOMAIN_NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
!                              ne=DOMAIN_NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
!                              DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
!                              DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                               !Loop over the local faces of the element
!                              DO nae=1,DOMAIN_ELEMENT%BASIS%NUMBER_OF_LOCAL_FACES
!                                nf2=DECOMPOSITION_ELEMENT%ELEMENT_FACES(nae)
!                                IF(nf2/=nf) THEN
!                                  DECOMPOSITION_FACE2=>DECOMPOSITION_FACES%FACES(nf2)
!                                  DOMAIN_FACE2=>DOMAIN_FACES%FACES(nf2)
                                   !Check whether XI of face have same direction
!                                  IF ((OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(nae),2,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(nae),2,1)).OR.&
!                                     &(OTHER_XI_DIRECTIONS3(BASIS%LOCAL_FACE_XI_DIRECTION(nae),3,1)==&
!                                     &OTHER_XI_DIRECTIONS3(BASIS2%LOCAL_FACE_XI_DIRECTION(nae),3,1))) THEN
                                     !Loop over nodes in face of surrounding element
!                                    BASIS2=>DOMAIN_FACE2%BASIS
!                                    IF(BASIS2%INTERPOLATION_ORDER(1)==BASIS%INTERPOLATION_ORDER(1)) THEN
!                                      NODE_COUNT=0
!                                      DO node_idx3=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                        DO node_idx4=1,BASIS%NUMBER_OF_NODES_IN_XI_DIRECTION
!                                          np2=DOMAIN_FACE2%NODES_IN_FACE((node_idx4*(BASIS2%NUMBER_OF_FACES-1))&
!                                                                      &+(node_idx3*(BASIS2%NUMBER_OF_NODES_IN_XI_DIRECTION-1))+1)
!                                          IF(np2==np) NODE_COUNT=NODE_COUNT+1
!                                        ENDDO !node_idx4
!                                      ENDDO !node_idx3
!                                      IF(NODE_COUNT<BASIS%NUMBER_OF_NODES) THEN
!                                        FOUND=.TRUE.
!                                        EXIT
!                                      ENDIF
!                                    ENDIF
!                                  ENDIF
!                                ENDIF
!                              ENDDO !nae
!                                IF(FOUND) EXIT
!                            ENDDO !elem_idx
!                            IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx2)=nf2
!                           ENDDO !node_idx2
!                           IF(FOUND) DECOMPOSITION_FACE%ADJACENT_FACES(node_idx1)=nf2
!                          ENDDO !node_idx1
                        ENDDO !nf

                        !Set the surrounding elements
                        DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                          DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                          BASIS=>DOMAIN_ELEMENT%BASIS
                          DO nae=1,BASIS%NUMBER_OF_LOCAL_FACES
                            FACE_NUMBER=DECOMPOSITION_ELEMENT%ELEMENT_FACES(nae)
                            DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(FACE_NUMBER)
                            DO nf=1,DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS
                              DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(nf)=ne
                              DECOMPOSITION_FACE%ELEMENT_FACES(nf)=nae
                            ENDDO
                          ENDDO !nae
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
                            DO nf=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
                              DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(nf)
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(nf)
                              IF(DECOMPOSITION_FACE%NUMBER_OF_SURROUNDING_ELEMENTS>0) THEN
                                ne=DECOMPOSITION_FACE%SURROUNDING_ELEMENTS(1)
                                nae=DECOMPOSITION_FACE%ELEMENT_FACES(1)
                                CALL DOMAIN_TOPOLOGY_FACE_INITIALISE(DOMAIN_FACES%FACES(nf),ERR,ERROR,*999)
                                DOMAIN_FACE%NUMBER=nf
                                DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                                BASIS=>DOMAIN_ELEMENT%BASIS
                                DOMAIN_FACE%BASIS=>BASIS%FACE_BASES(DECOMPOSITION_FACE%XI_DIRECTION)%PTR
                                ALLOCATE(DOMAIN_FACE%NODES_IN_FACE(BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes in face",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_FACE%DERIVATIVES_IN_FACE(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivatives in face",ERR,ERROR,*999)
                                !Set derivatives of nodes in domain face from derivatives of nodes in element
                                DO nnf=1,BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(nae)
                                  nn=BASIS%NODE_NUMBERS_IN_LOCAL_FACE(nnf,nae)
                                  np=DOMAIN_ELEMENT%ELEMENT_NODES(nn)
                                  DOMAIN_FACE%NODES_IN_FACE(nnf)=np
                                  DOMAIN_FACE%DERIVATIVES_IN_FACE(1,nnf)=NO_GLOBAL_DERIV
                                  IF(DOMAIN_FACE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES>1) THEN
                                    nk=DOMAIN_ELEMENT%ELEMENT_DERIVATIVES(BASIS%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(nnf,nae),nn)
                                    DOMAIN_FACE%DERIVATIVES_IN_FACE(2,nnf)=nk
                                  ENDIF
                                  NODES_NUMBER_OF_FACES(np)=NODES_NUMBER_OF_FACES(np)+1
                                ENDDO !nnf
                              ELSE
                                CALL FLAG_ERROR("Face is not surrounded by any elements?",ERR,ERROR,*999)
                              ENDIF                              
                            ENDDO !nf
                            DO np=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(np)%NODE_FACES(NODES_NUMBER_OF_FACES(np)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node faces",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_FACES=0
                            ENDDO !np
                            DEALLOCATE(NODES_NUMBER_OF_FACES)
                            DO nf=1,DOMAIN_FACES%NUMBER_OF_FACES
                              DOMAIN_FACE=>DOMAIN_FACES%FACES(nf)
                              BASIS=>DOMAIN_FACE%BASIS
                              DO nnf=1,BASIS%NUMBER_OF_NODES
                                np=DOMAIN_FACE%NODES_IN_FACE(nnf)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(np)
                                DOMAIN_NODE%NUMBER_OF_NODE_FACES=DOMAIN_NODE%NUMBER_OF_NODE_FACES+1
                                !Set the face numbers a node is on
                                DOMAIN_NODE%NODE_FACES(DOMAIN_NODE%NUMBER_OF_NODE_FACES)=nf
                              ENDDO !nnf
                            ENDDO !nf
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
      DO nf=1,DECOMPOSITION_FACES%NUMBER_OF_FACES
        DECOMPOSITION_FACE=>DECOMPOSITION_FACES%FACES(nf)
        DOMAIN_FACE=>DOMAIN_FACES%FACES(nf)
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
          DOMAIN_FACE=>DOMAIN%TOPOLOGY%FACES%FACES(nf)
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
          DO nnf=1,DOMAIN_FACE%BASIS%NUMBER_OF_NODES
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"          Node : ",nnf,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_FACE%BASIS%NUMBER_OF_DERIVATIVES(nnf),4,4,DOMAIN_FACE% &
              & DERIVATIVES_IN_FACE(:,nnf),'("            Derivatives in face  :",4(X,I8))','(34X,4(X,I8))',ERR,ERROR,*999)
          ENDDO !nnf
        ENDDO !component_idx
      ENDDO !nf
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
    INTEGER(INTG), ALLOCATABLE :: LOCAL_ELEMENT_NUMBERS(:)
    INTEGER(INTG), POINTER :: DOMAINS(:),ADJACENT_ELEMENTS(:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    NULLIFY(DOMAINS)
    NULLIFY(ADJACENT_ELEMENTS)
    
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
                IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
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
999 IF(ASSOCIATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)    
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
      & NUMBER_OF_NODES_PER_DOMAIN,domain_idx,domain_idx2,domain_no,nk,np,ny,NUMBER_OF_DOMAINS,MAX_NUMBER_DOMAINS, &
      & NUMBER_OF_GHOST_NODES,my_computational_node_number,number_computational_nodes,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:),NODE_COUNT(:),NUMBER_INTERNAL_NODES(:), &
      & NUMBER_BOUNDARY_NODES(:)
    INTEGER(INTG), POINTER :: DOMAINS(:),ALL_DOMAINS(:),GHOST_NODES(:)
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

    NULLIFY(DOMAINS)
    NULLIFY(ALL_DOMAINS)
    NULLIFY(GHOST_NODES)
    
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
                  DO np=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
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
                    DO no_adjacent_element=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                      adjacent_element=MESH_TOPOLOGY%NODES%NODES(np)%SURROUNDING_ELEMENTS(no_adjacent_element)
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
                        IF(.NOT.BOUNDARY_DOMAIN) CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,np,ERR,ERROR,*999)
                      ENDDO !domain_idx
                    ENDIF
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map local number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map domain number.",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map local type.",ERR,ERROR,*999)
                    DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                      ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local number.",ERR,ERROR,*999)
                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map domain number.",ERR,ERROR,*999)
                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local type.",ERR,ERROR,*999)
                    ENDDO !nk
                    IF(NUMBER_OF_DOMAINS==1) THEN
                      !Node is an internal node
                      domain_no=DOMAINS(1)
                      NUMBER_INTERNAL_NODES(domain_no)=NUMBER_INTERNAL_NODES(domain_no)+1
                      !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS=1
                      !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(DOMAINS(1))
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(1)=-1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(1)=DOMAINS(1) 
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                      DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                        ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                        !LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                        !DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=-1
                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_INTERNAL
                      ENDDO !nk
                    ELSE
                      !Node is on the boundary of computational domains
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                      DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                        ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
                      ENDDO !nk
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=DOMAINS(domain_idx)
                        !LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                        !NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(domain_idx)=LOCAL_NODE_NUMBERS(domain_no)
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(domain_idx)=-1
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(domain_idx)=domain_no
                        NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                        DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                          ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                          !LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                          !DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=LOCAL_DOF_NUMBERS(domain_no)
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(domain_idx)=-1
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(domain_idx)=domain_no
                          DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(domain_idx)=DOMAIN_LOCAL_BOUNDARY
                        ENDDO !nk
                      ENDDO !domain_idx
                    ENDIF
                    DEALLOCATE(DOMAINS) 
                    DEALLOCATE(ALL_DOMAINS)
                  ENDDO !np

                  !For the second pass assign boundary nodes to one domain on the boundary and set local node numbers.
                  NUMBER_OF_NODES_PER_DOMAIN=FLOOR(REAL(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES,DP)/ &
                    & REAL(DECOMPOSITION%NUMBER_OF_DOMAINS,DP))
                  ALLOCATE(DOMAIN%NODE_DOMAIN(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node domain",ERR,ERROR,*999)
                  DOMAIN%NODE_DOMAIN=-1
                  DO np=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                    IF(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS==1) THEN !Internal node
                      domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(1)
                      DOMAIN%NODE_DOMAIN(np)=domain_no
                      LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                      NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                      DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                        ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                        LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                        DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                      ENDDO !nk
                    ELSE !Boundary node
                      NUMBER_OF_DOMAINS=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS
                      DO domain_idx=1,NUMBER_OF_DOMAINS
                        domain_no=NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(domain_idx)
                        IF(DOMAIN%NODE_DOMAIN(np)<0) THEN
                          IF((NUMBER_INTERNAL_NODES(domain_no)+NUMBER_BOUNDARY_NODES(domain_no)<NUMBER_OF_NODES_PER_DOMAIN).OR. &
                            & (domain_idx==NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS)) THEN
                            !Allocate the node to this domain
                            DOMAIN%NODE_DOMAIN(np)=domain_no
                            NUMBER_BOUNDARY_NODES(domain_no)=NUMBER_BOUNDARY_NODES(domain_no)+1
                            LOCAL_NODE_NUMBERS(domain_no)=LOCAL_NODE_NUMBERS(domain_no)+1
                            !Reset the boundary information to be in the first domain index. The remaining domain indicies will
                            !be overwritten when the ghost nodes are calculated below. 
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS=1                            
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER(1)=LOCAL_NODE_NUMBERS(domain_no)
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(1)=domain_no
                            NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                            DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                              ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                              LOCAL_DOF_NUMBERS(domain_no)=LOCAL_DOF_NUMBERS(domain_no)+1
                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%NUMBER_OF_DOMAINS=1
                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(1)=LOCAL_DOF_NUMBERS(domain_no)
                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(1)=domain_no
                              DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(1)=DOMAIN_LOCAL_BOUNDARY
                            ENDDO !nk
                          ELSE
                            !The node as already been assigned to a domain so it must be a ghost node in this domain
                            CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,np,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          !The node as already been assigned to a domain so it must be a ghost node in this domain
                          CALL LIST_ITEM_ADD(GHOST_NODES_LIST(domain_no)%PTR,np,ERR,ERROR,*999)
                        ENDIF
                      ENDDO !domain_idx
                    ENDIF
                  ENDDO !np
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
                      DO nk=1,MESH_TOPOLOGY%NODES%NODES(ghost_node)%NUMBER_OF_DERIVATIVES
                        ny=MESH_TOPOLOGY%NODES%NODES(ghost_node)%DOF_INDEX(nk)
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
                      ENDDO !nk
                    ENDDO !no_ghost_node
                    DEALLOCATE(GHOST_NODES)
                  ENDDO !domain_idx
                  
                  !Check decomposition and check that each domain has a node in it.
                  ALLOCATE(NODE_COUNT(0:number_computational_nodes-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node count.",ERR,ERROR,*999)
                  NODE_COUNT=0
                  DO np=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                    no_computational_node=DOMAIN%NODE_DOMAIN(np)
                    IF(no_computational_node>=0.AND.no_computational_node<number_computational_nodes) THEN
                      NODE_COUNT(no_computational_node)=NODE_COUNT(no_computational_node)+1
                    ELSE
                      LOCAL_ERROR="The computational node number of "// &
                        & TRIM(NUMBER_TO_VSTRING(no_computational_node,"*",ERR,ERROR))// &
                        & " for node number "//TRIM(NUMBER_TO_VSTRING(np,"*",ERR,ERROR))// &
                        & " is invalid. The computational node number must be between 0 and "// &
                        & TRIM(NUMBER_TO_VSTRING(number_computational_nodes-1,"*",ERR,ERROR))//"."
                      CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !np
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
      DO np=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Node = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Domain = ",DOMAIN%NODE_DOMAIN(np),ERR,ERROR,*999)
      ENDDO !np
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Node mappings :",ERR,ERROR,*999)
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map :",ERR,ERROR,*999)
      DO np=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
          & NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_NUMBER, &
          & '("      Local number :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)% &
            & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER, &
            & '("      Domain number:",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)% &
          & NUMBER_OF_DOMAINS,8,8,NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_TYPE, &
          & '("      Local type   :",8(X,I7))','(20X,8(X,I7))',ERR,ERROR,*999)      
      ENDDO !np     
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map :",ERR,ERROR,*999)
      DO np=1,NODES_MAPPING%TOTAL_NUMBER_OF_LOCAL
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local node = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global node = ", &
          & NODES_MAPPING%LOCAL_TO_GLOBAL_MAP(np),ERR,ERROR,*999)
      ENDDO !np
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
      ENDDO !np
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
999 IF(ASSOCIATED(DOMAINS)) DEALLOCATE(DOMAINS)
    IF(ASSOCIATED(ALL_DOMAINS)) DEALLOCATE(ALL_DOMAINS)
    IF(ASSOCIATED(GHOST_NODES)) DEALLOCATE(GHOST_NODES)
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
    INTEGER(INTG) :: local_element,global_element,local_node,global_node,ne,nn,nk,nkk,np,ny,component_idx
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
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain elements elements",ERR,ERROR,*999)
              DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%NUMBER_OF_LOCAL
              DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN%MAPPINGS%ELEMENTS%TOTAL_NUMBER_OF_LOCAL
              ALLOCATE(DOMAIN_NODES%NODES(DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain nodes nodes",ERR,ERROR,*999)
              DOMAIN_NODES%NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%NUMBER_OF_LOCAL
              DOMAIN_NODES%TOTAL_NUMBER_OF_NODES=DOMAIN%MAPPINGS%NODES%TOTAL_NUMBER_OF_LOCAL
              ALLOCATE(DOMAIN_DOFS%DOF_INDEX(2,DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain dofs dof index",ERR,ERROR,*999)
              DOMAIN_DOFS%NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%NUMBER_OF_LOCAL
              DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS=DOMAIN%MAPPINGS%DOFS%TOTAL_NUMBER_OF_LOCAL
              !Loop over the domain nodes and calculate the parameters from the mesh nodes
              ny=0
              DO local_node=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                CALL DOMAIN_TOPOLOGY_NODE_INITIALISE(DOMAIN_NODES%NODES(local_node),ERR,ERROR,*999)
                global_node=DOMAIN%MAPPINGS%NODES%LOCAL_TO_GLOBAL_MAP(local_node)
                DOMAIN_NODES%NODES(local_node)%LOCAL_NUMBER=local_node
                DOMAIN_NODES%NODES(local_node)%MESH_NUMBER=global_node
                DOMAIN_NODES%NODES(local_node)%GLOBAL_NUMBER=MESH_NODES%NODES(global_node)%GLOBAL_NUMBER
                DOMAIN_NODES%NODES(local_node)%USER_NUMBER=MESH_NODES%NODES(global_node)%USER_NUMBER
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_SURROUNDING_ELEMENTS=0
                NULLIFY(DOMAIN_NODES%NODES(local_node)%SURROUNDING_ELEMENTS)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES=MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%GLOBAL_DERIVATIVE_INDEX( &
                  & MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global derivative index.",ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%GLOBAL_DERIVATIVE_INDEX=MESH_NODES%NODES(global_node)%GLOBAL_DERIVATIVE_INDEX
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX( &
                  & MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node partial derivative index.",ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX=MESH_NODES%NODES(global_node)%PARTIAL_DERIVATIVE_INDEX
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%DOF_INDEX(MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node dof index.",ERR,ERROR,*999)
                DO nk=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                  ny=ny+1
                  DOMAIN_NODES%NODES(local_node)%DOF_INDEX(nk)=ny
                  DOMAIN_DOFS%DOF_INDEX(1,ny)=nk
                  DOMAIN_DOFS%DOF_INDEX(2,ny)=local_node
                ENDDO !nk
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
                ALLOCATE(DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                  & BASIS%NUMBER_OF_NODES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain elements element derivatives",ERR,ERROR,*999)
                DO nn=1,BASIS%NUMBER_OF_NODES
                  global_node=MESH_ELEMENTS%ELEMENTS(global_element)%MESH_ELEMENT_NODES(nn)
                  local_node=DOMAIN%MAPPINGS%NODES%GLOBAL_TO_LOCAL_MAP(global_node)%LOCAL_NUMBER(1)
                  DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_NODES(nn)=local_node
                  DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                    !Find equivalent node derivative by matching partial derivative index
                    FOUND=.FALSE.
                    DO nkk=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                      IF(DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX(nkk)==BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn)) THEN
                        FOUND=.TRUE.
                        EXIT
                      ENDIF
                    ENDDO !nkk
                    IF(FOUND) THEN
                      DOMAIN_ELEMENTS%ELEMENTS(local_element)%ELEMENT_DERIVATIVES(nk,nn)=nkk
                    ELSE
                      CALL FLAG_ERROR("Could not find equivalent node derivative",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !nk
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
      DO np=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node number = ",DOMAIN_NODES%NODES(np)%LOCAL_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node mesh number = ",DOMAIN_NODES%NODES(np)%MESH_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node global number = ",DOMAIN_NODES%NODES(np)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",DOMAIN_NODES%NODES(np)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,8,8, &
          & DOMAIN_NODES%NODES(np)%GLOBAL_DERIVATIVE_INDEX,'("      Global derivative index(nk) :",8(X,I2))','(36X,8(X,I2))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,8,8, &
          & DOMAIN_NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX,'("      Partial derivative index(nk) :",8(X,I2))','(36X,8(X,I2))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,4,4, &
          & DOMAIN_NODES%NODES(np)%DOF_INDEX,'("      Degree-of-freedom index(nk)  :",4(X,I9))','(36X,4(X,I9))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Boundary node = ",DOMAIN_NODES%NODES(np)%BOUNDARY_NODE,ERR,ERROR,*999)
     ENDDO !np
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total Number of domain dofs = ",DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS, &
        & ERR,ERROR,*999)
      DO ny=1,DOMAIN_DOFS%TOTAL_NUMBER_OF_DOFS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Dof number = ",ny,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,2,2,2, &
          & DOMAIN_DOFS%DOF_INDEX(:,ny),'("    Degree-of-freedom index :",2(X,I9))','(29X,2(X,I9))', &
          & ERR,ERROR,*999)
      ENDDO !ny
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
            & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_DERIVATIVES(:,nn), &
            & '("        Element derivatives(nk) :",8(X,I2))','(33X,8(X,I2))',ERR,ERROR,*999)
        ENDDO !nn
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

  !>Finalises the given domain topology node and deallocates all memory.
  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE !<The domain node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%GLOBAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%GLOBAL_DERIVATIVE_INDEX)
    IF(ALLOCATED(NODE%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%PARTIAL_DERIVATIVE_INDEX)
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
        TOPOLOGY%NODES%MAXIMUM_NUMBER_OF_DERIVATIVES=0
        TOPOLOGY%NODES%DOMAIN=>TOPOLOGY%DOMAIN
        NULLIFY(TOPOLOGY%NODES%NODES)
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
          CALL MESH_TOPOLOGY_CALCULATE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
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
      MESH%NUMBER_OF_DIMENSIONS=0
      MESH%NUMBER_OF_COMPONENTS=0
      MESH%MESH_EMBEDDED=.FALSE.
      NULLIFY(MESH%EMBEDDING_MESH)
      MESH%NUMBER_OF_EMBEDDED_MESHES=0
      NULLIFY(MESH%EMBEDDED_MESHES)
      MESH%NUMBER_OF_ELEMENTS=0
      MESH%NUMBER_OF_FACES=0
      MESH%NUMBER_OF_LINES=0
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

  !>Calculates the mesh topology.
  SUBROUTINE MESH_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)    

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Calculate the nodes used in the mesh
      CALL MESH_TOPOLOGY_NODES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the elements surrounding the nodes in a mesh
      CALL MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the number of derivatives at each nodes in a mesh
      CALL MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the elements surrounding the elements in the mesh
      CALL MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the boundary nodes and elements in the mesh
      CALL MESH_TOPOLOGY_BOUNDARY_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
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
    INTEGER(INTG) :: element_idx,MATCH_INDEX,ni,nn,node_idx,XI_DIRECTION,n0
    INTEGER(INTG) :: LOCAL_FACE_NUMBER,LOCAL_INTERFACE_NUMBER,NUMBER_OF_NODES_IN_LOCAL_INTERFACE
    TYPE(BASIS_TYPE), POINTER :: BASIS
    
    CALL ENTERS("MESH_TOPOLOGY_BOUNDARY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN

      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN        
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          DO element_idx=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BASIS
            n0=-BASIS%NUMBER_OF_XI_COORDINATES !
            IF(BASIS%TYPE .EQ. 2) n0=1 ! simplex starts 1
            DO ni=n0,BASIS%NUMBER_OF_XI_COORDINATES ! coordinates should be correct
              IF(ni/=0) THEN
                IF(TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%NUMBER_OF_ADJACENT_ELEMENTS(ni)==0) THEN
                  TOPOLOGY%ELEMENTS%ELEMENTS(element_idx)%BOUNDARY_ELEMENT=.TRUE.
                  
                  ! lets find local face number, we find local interface from its xi data
                  IF (BASIS%TYPE .EQ. 1) LOCAL_INTERFACE_NUMBER = 2*ABS(ni)-(1+ABS(ni)/ni)/2
                  IF (BASIS%TYPE .EQ. 2) LOCAL_INTERFACE_NUMBER = ni
                  
                  ! lets find number of nodes on that interface 
                  IF (BASIS%NUMBER_OF_XI .EQ. 2) NUMBER_OF_NODES_IN_LOCAL_INTERFACE = & 
                   & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(LOCAL_INTERFACE_NUMBER)
                  IF (BASIS%NUMBER_OF_XI .EQ. 3) NUMBER_OF_NODES_IN_LOCAL_INTERFACE = &
                   & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(LOCAL_INTERFACE_NUMBER)
                  
                  ! now set them as boundary nodes 
                  DO nn=1,NUMBER_OF_NODES_IN_LOCAL_INTERFACE
                    IF (BASIS%NUMBER_OF_XI .EQ. 2) node_idx=TOPOLOGY%ELEMENTS%ELEMENTS(element_idx) & 
                     & %MESH_ELEMENT_NODES(BASIS%NODE_NUMBERS_IN_LOCAL_LINE(nn,LOCAL_INTERFACE_NUMBER))
                    IF (BASIS%NUMBER_OF_XI .EQ. 3) node_idx=TOPOLOGY%ELEMENTS%ELEMENTS(element_idx) & 
                     & %MESH_ELEMENT_NODES(BASIS%NODE_NUMBERS_IN_LOCAL_FACE(nn,LOCAL_INTERFACE_NUMBER))
                    TOPOLOGY%NODES%NODES(node_idx)%BOUNDARY_NODE=.TRUE.
                  ENDDO !nn
                  
                ENDIF
              ENDIF
            ENDDO !ni            
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
    INTEGER(INTG) :: np,nk,NUMBER_OF_DOFS

    CALL ENTERS("MESH_TOPOLOGY_DOFS_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN        
        IF(ASSOCIATED(TOPOLOGY%DOFS)) THEN
          NUMBER_OF_DOFS=0
          DO np=1,TOPOLOGY%NODES%NUMBER_OF_NODES
            ALLOCATE(TOPOLOGY%NODES%NODES(np)%DOF_INDEX(TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh topology node dof index",ERR,ERROR,*999)
            DO nk=1,TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
              NUMBER_OF_DOFS=NUMBER_OF_DOFS+1
              TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)=NUMBER_OF_DOFS
            ENDDO !nk
          ENDDO !np
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
          IF(ALLOCATED(ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES)) THEN
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES,8,8, &
              & ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES,'("    Mesh element nodes   =",8(X,I6))','(26X,8(X,I6))',ERR,ERROR,*999)
          ELSE
            CALL FLAG_ERROR("Mesh element nodes are not associated.",ERR,ERROR,*999)
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
        LOCAL_ERROR="The speficied mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & "is invalid. The component number must be between 1 and "// &
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

    CALL ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODES)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%GLOBAL_ELEMENT_NODES)) DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%MESH_ELEMENT_NODES)) DEALLOCATE(ELEMENT%MESH_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS)) DEALLOCATE(ELEMENT%NUMBER_OF_ADJACENT_ELEMENTS)
    IF(ALLOCATED(ELEMENT%ADJACENT_ELEMENTS)) DEALLOCATE(ELEMENT%ADJACENT_ELEMENTS)
  
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
        LOCAL_ERROR="The speficied mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & "is invalid. The component number must be between 1 and "// &
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
 !>Calculates the element numbers surrounding an element in a mesh topology.
  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the elements adjacent to elements for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: j,ne1,nep1,nn,nn1,nn2,np,np1,DUMMY_ERR,FACE_XI(2),NODE_POSITION_INDEX(3)
    INTEGER(INTG) :: direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: i,ne,ni,is,COUNTER1,COUNTER2,COUNTER3,NUMBER_OF_ITEM,XI_DIRECTION
    INTEGER(INTG) :: NODE_NUMBER_IN_LOCAL_INTERFACE,NUMBER_OF_LOCAL_INTERFACE,NUMBER_OF_NODES_IN_LOCAL_INTERFACE
    INTEGER(INTG) :: NUMBER_OF_SURROUNDING_ELEMENTS
    INTEGER(INTG) :: NUMBER_SURROUNDING,MAX_NUMBER_SURROUNDING,NUMBER_OF_NODES_XI(3)
    INTEGER(INTG), POINTER :: ADJACENT_ELEMENTS(:),SURROUNDING_ELEMENTS(:),ELEMENTS_MATCH(:)
!    LOGICAL :: XI_COLLAPSED,FACE_COLLAPSED(-3:3),SUBSET
!    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    TYPE(LIST_PTR_TYPE)  :: ADJACENT_ELEMENTS_LIST(-4:4)
    TYPE(LIST_TYPE), POINTER :: ELEMENTS_MATCH_LIST ! JUST CONTAINS LIST OF SURROUNDING ELEMENTS AND NO MORE
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: SURROUNDING_ELEMENTS_LIST(:) 
    TYPE(BASIS_TYPE), POINTER :: BASIS

    NULLIFY(SURROUNDING_ELEMENTS)
    NULLIFY(ADJACENT_ELEMENTS)
    NULLIFY(ELEMENTS_MATCH)
    NULLIFY(BASIS)

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR,*999)
    !%%% Here is the algorithm

!1  DO i=1,total_number_of_elements
!2    DO node=1,number_of_nodes_per_this_element
!3      ADD (surrounding element number to ELEMENTS_MATCH_LIST) and (local node number to SURROUNDING_ELEMENTS_LIST(i)%ptr)
!4    CHECK data in SURROUNDING_ELEMENTS_LIST(i)%ptr with faces of element to find the an ELEMENTS_MATCH_LIST(i) that is adjacent in face. 
!5    FIND the corresponding xi direction for this face
!6    ADD ELEMENTS_MATCH_LIST(i) to ADJACENT_ELEMENT and the corresponding xi data in TOPOLOGY%ELEMENTS%ELEMENTS dataset
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          ALLOCATE(SURROUNDING_ELEMENTS_LIST(TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS),STAT=ERR)
          !Loop over the global elements in the mesh         
          DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
          !%%%% first we initialize lists that are required to find the adjacent elements list
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
            DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              NULLIFY(ADJACENT_ELEMENTS_LIST(ni)%PTR)
              CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,5,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
            ENDDO !ni
                                     
            NULLIFY(ELEMENTS_MATCH_LIST) ! Contains element number of the surrounding elements 
            CALL LIST_CREATE_START(ELEMENTS_MATCH_LIST,ERR,ERROR,*999)
            CALL LIST_DATA_TYPE_SET(ELEMENTS_MATCH_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
            CALL LIST_INITIAL_SIZE_SET(ELEMENTS_MATCH_LIST,16,ERR,ERROR,*999)
            CALL LIST_CREATE_FINISH(ELEMENTS_MATCH_LIST,ERR,ERROR,*999)
            DO COUNTER1=1,ELEMENTS_MATCH_LIST%INITIAL_SIZE ! SET INITIAL VALUES TO ZERO
              ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1)=0
            ENDDO              
            
            NUMBER_OF_SURROUNDING_ELEMENTS=0
            !%%% now loop over the local nodes                    
            DO COUNTER1=1,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS%NUMBER_OF_NODES
              nn=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(COUNTER1)  
              !%%% loop over the surrounding elements to the selected node
              DO COUNTER2=1,TOPOLOGY%NODES%NODES(nn)%NUMBER_OF_SURROUNDING_ELEMENTS
                CALL LIST_SEARCH(ELEMENTS_MATCH_LIST%LIST_INTG,TOPOLOGY%NODES%NODES(nn)&
                 &%SURROUNDING_ELEMENTS(COUNTER2),NUMBER_OF_ITEM,ERR,ERROR,*999)
                IF (NUMBER_OF_ITEM .EQ. 0) THEN
                  NUMBER_OF_SURROUNDING_ELEMENTS=NUMBER_OF_SURROUNDING_ELEMENTS+1
                  NULLIFY(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR)
                  CALL LIST_CREATE_START(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,ERR,ERROR,*999)
                  CALL LIST_DATA_TYPE_SET(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,LIST_INTG_TYPE, &
                   & ERR,ERROR,*999)
                  CALL LIST_INITIAL_SIZE_SET(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,5,ERR,ERROR,*999)
                  CALL LIST_CREATE_FINISH(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,ERR,ERROR,*999)
     
                  DO COUNTER3=1,SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR%INITIAL_SIZE ! SET INITIAL VALUES ZERO
                    SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR%LIST_INTG(COUNTER3)=0
                  ENDDO                          
                  
                  CALL LIST_ITEM_ADD(ELEMENTS_MATCH_LIST,TOPOLOGY%NODES%NODES(nn)&
                   &%SURROUNDING_ELEMENTS(COUNTER2),ERR,ERROR,*999)
                  CALL LIST_ITEM_ADD(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_SURROUNDING_ELEMENTS)%PTR,COUNTER1,ERR,ERROR,*999)
                ELSE
                  CALL LIST_ITEM_ADD(SURROUNDING_ELEMENTS_LIST(NUMBER_OF_ITEM)%PTR,COUNTER1,ERR,ERROR,*999)
                ENDIF
              ENDDO ! COUNTER2
            ENDDO ! COUNTER1
            
            !%%% now surrounding elements and shared nodes are stored, we want to find the adjacent elements and correspounding xi direction
            ! loading ADJACENT_ELEMENTS_LIST on xi=0  
            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER, &
             & ERR,ERROR,*999)
             
            CALL LIST_NUMBER_OF_ITEMS_GET(ELEMENTS_MATCH_LIST,NUMBER_OF_SURROUNDING_ELEMENTS,ERR,ERROR,*999)
            DO COUNTER1=1,NUMBER_OF_SURROUNDING_ELEMENTS
              CALL LIST_NUMBER_OF_ITEMS_GET(SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR,NUMBER_OF_ITEM,ERR,ERROR,*999)
              IF (NUMBER_OF_ITEM .LT. BASIS%NUMBER_OF_NODES .AND. NUMBER_OF_ITEM .GT. 1) THEN ! if /= element ITSELF and not attached in 1 node
                ! first lets find the interface number
                IF (BASIS%NUMBER_OF_XI .EQ. 2) NUMBER_OF_LOCAL_INTERFACE = BASIS%NUMBER_OF_LOCAL_LINES
                IF (BASIS%NUMBER_OF_XI .EQ. 3) NUMBER_OF_LOCAL_INTERFACE = BASIS%NUMBER_OF_LOCAL_FACES
                
                DO COUNTER2=1,NUMBER_OF_LOCAL_INTERFACE
                  is = 1
                  IF (BASIS%NUMBER_OF_XI .EQ. 2) NUMBER_OF_NODES_IN_LOCAL_INTERFACE = &
                   & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(COUNTER2)
                  IF (BASIS%NUMBER_OF_XI .EQ. 3) NUMBER_OF_NODES_IN_LOCAL_INTERFACE = &
                   & BASIS%NUMBER_OF_NODES_IN_LOCAL_FACE(COUNTER2)
                  DO COUNTER3=1,NUMBER_OF_NODES_IN_LOCAL_INTERFACE
                    IF (BASIS%NUMBER_OF_XI .EQ. 2) NODE_NUMBER_IN_LOCAL_INTERFACE = &
                     & BASIS%NODE_NUMBERS_IN_LOCAL_LINE(COUNTER3,COUNTER2)
                    IF (BASIS%NUMBER_OF_XI .EQ. 3) NODE_NUMBER_IN_LOCAL_INTERFACE = & 
                     & BASIS%NODE_NUMBERS_IN_LOCAL_FACE(COUNTER3,COUNTER2)
                    CALL LIST_SEARCH(SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR%LIST_INTG,NODE_NUMBER_IN_LOCAL_INTERFACE, &
                     & NUMBER_OF_ITEM,ERR,ERROR,*999)
                    IF (NUMBER_OF_ITEM .EQ. 0) is=0
                  ENDDO !COUNTER3
                  IF (is .EQ. 1) THEN ! (COUNTER2) is the interface number of the FACE/LINE that is shared between element DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER and ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1)
                    IF (BASIS%TYPE .EQ. 2) XI_DIRECTION = COUNTER2
                    IF (BASIS%TYPE .EQ. 1) THEN
                      IF (BASIS%NUMBER_OF_XI .EQ. 2) XI_DIRECTION = BASIS%LOCAL_LINE_XI_DIRECTION(COUNTER2)
                      IF (BASIS%NUMBER_OF_XI .EQ. 3) XI_DIRECTION = BASIS%LOCAL_FACE_XI_DIRECTION(COUNTER2)
                      XI_DIRECTION = XI_DIRECTION * (-1)**(COUNTER2+1)! FOR LAGRANGE TYPE SHOULD DETECT NEGATIVE DIRECTION
                    ENDIF
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(XI_DIRECTION)%PTR,ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1), &
                     & ERR,ERROR,*999)
                  ENDIF
                ENDDO !COUNTER2
                 
              ENDIF
            ENDDO !COUNTER1
             
            !%%% set maximum number of adjacent elements
            MAX_NUMBER_SURROUNDING=1
            DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              CALL LIST_NUMBER_OF_ITEMS_GET(ADJACENT_ELEMENTS_LIST(ni)%PTR,i,ERR,ERROR,*999)
              IF(i>MAX_NUMBER_SURROUNDING) MAX_NUMBER_SURROUNDING=i
            ENDDO ! ni
            
            !Set the surrounding elements for this element
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI_COORDINATES: &
              & BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of surrounding elements",ERR,ERROR,*999)
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(MAX_NUMBER_SURROUNDING, &
              & -BASIS%NUMBER_OF_XI_COORDINATES:BASIS%NUMBER_OF_XI_COORDINATES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate surrounding elements",ERR,ERROR,*999)
            DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
              CALL LIST_DETACH_AND_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & NUMBER_OF_ADJACENT_ELEMENTS(ni),ADJACENT_ELEMENTS,ERR,ERROR,*999)
              TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & NUMBER_OF_ADJACENT_ELEMENTS(ni),ni) = ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni))
              IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
            ENDDO !ni
            
            ! destroy lists information
            DO COUNTER1=1,NUMBER_OF_SURROUNDING_ELEMENTS
              ELEMENTS_MATCH_LIST%LIST_INTG(COUNTER1)=0
            ENDDO              
!            CALL LIST_DETACH_AND_DESTROY(ELEMENTS_MATCH_LIST,NUMBER_OF_SURROUNDING_ELEMENTS,ELEMENTS_MATCH,ERR,ERROR,*999)
!            DEALLOCATE(ELEMENTS_MATCH)

            DO COUNTER1=1,NUMBER_OF_SURROUNDING_ELEMENTS
              CALL LIST_NUMBER_OF_ITEMS_GET(SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR,NUMBER_OF_ITEM,ERR,ERROR,*999)
              DO COUNTER2=1,NUMBER_OF_ITEM
                SURROUNDING_ELEMENTS_LIST(COUNTER1)%PTR%LIST_INTG(COUNTER2)=0
              ENDDO
            ENDDO                          
!            DO ni=1,NUMBER_OF_SURROUNDING_ELEMENTS
!              CALL LIST_NUMBER_OF_ITEMS_GET(SURROUNDING_ELEMENTS_LIST(ni)%PTR,NUMBER_OF_ITEM,ERR,ERROR,*999)
!              CALL LIST_DETACH_AND_DESTROY(SURROUNDING_ELEMENTS_LIST(ni)%PTR,NUMBER_OF_ITEM,SURROUNDING_ELEMENTS,ERR,ERROR,*999)
!              DEALLOCATE(SURROUNDING_ELEMENTS)
!            ENDDO !ni
          ENDDO !ne
        ELSE
          CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology nodes is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not allocated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of elements = ",TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
        BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Global element number = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi directions = ", &
         & BASIS%NUMBER_OF_XI_COORDINATES,ERR,ERROR,*999)
        DO ni=-BASIS%NUMBER_OF_XI_COORDINATES,BASIS%NUMBER_OF_XI_COORDINATES
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",ni,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni),ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
            & NUMBER_OF_ADJACENT_ELEMENTS(ni),8,8,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(:,ni), &
            & '("        Adjacent elements =",8(X,I8))','(30x,8(X,I8))',ERR,ERROR,*999)
        ENDDO !ni
      ENDDO !ne
    ENDIF
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE")
    RETURN
    
999 IF(ASSOCIATED(SURROUNDING_ELEMENTS)) DEALLOCATE(SURROUNDING_ELEMENTS)
    IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(ELEMENTS_MATCH_LIST)) CALL LIST_DESTROY(ELEMENTS_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*997)
!998 DO ni=-4,4
!      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(ni)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
!    ENDDO !ni
!    DO ni=1,100
!      IF(ASSOCIATED(SURROUNDING_ELEMENTS_LIST(ni)%PTR)) &
!       & CALL LIST_DESTROY(SURROUNDING_ELEMENTS_LIST(ni)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
!    ENDDO !ni
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

    CALL ENTERS("MESH_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(NODE%GLOBAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%GLOBAL_DERIVATIVE_INDEX)
    IF(ALLOCATED(NODE%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%PARTIAL_DERIVATIVE_INDEX)
    IF(ALLOCATED(NODE%DOF_INDEX)) DEALLOCATE(NODE%DOF_INDEX)
  
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
      DO np=1,NODES%NUMBER_OF_NODES
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

  !>Calculates the number of derivatives at each node in a topology.
  SUBROUTINE MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY !<A pointer to the mesh topology to calculate the derivates at each node for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: elem_idx,global_deriv,MAX_NUMBER_OF_DERIVATIVES,ne,nk,nn,np,NUMBER_OF_DERIVATIVES
    INTEGER(INTG), POINTER :: DERIVATIVES(:)
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
          DO np=1,NODES%NUMBER_OF_NODES
            !Calculate the number of derivatives at each node. This needs to be calculated by looking at the mesh elements
            !as we may have an adjacent element in another domain with a higher order basis.
            NULLIFY(NODE_DERIVATIVE_LIST)
            CALL LIST_CREATE_START(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            CALL LIST_DATA_TYPE_SET(NODE_DERIVATIVE_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
            CALL LIST_INITIAL_SIZE_SET(NODE_DERIVATIVE_LIST,8,ERR,ERROR,*999)
            CALL LIST_CREATE_FINISH(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            MAX_NUMBER_OF_DERIVATIVES=-1
            DO elem_idx=1,NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
              ne=NODES%NODES(np)%SURROUNDING_ELEMENTS(elem_idx)
              BASIS=>ELEMENTS%ELEMENTS(ne)%BASIS
              !Find the local node corresponding to this node
              FOUND=.FALSE.
              DO nn=1,BASIS%NUMBER_OF_NODES
                IF(ELEMENTS%ELEMENTS(ne)%MESH_ELEMENT_NODES(nn)==np) THEN
                  FOUND=.TRUE.
                  EXIT
                ENDIF
              ENDDO !nn
              IF(FOUND) THEN
                DO nk=1,BASIS%NUMBER_OF_DERIVATIVES(nn)
                  CALL LIST_ITEM_ADD(NODE_DERIVATIVE_LIST,BASIS%PARTIAL_DERIVATIVE_INDEX(nk,nn),ERR,ERROR,*999)
                ENDDO !nk
                IF(BASIS%NUMBER_OF_DERIVATIVES(nn)>MAX_NUMBER_OF_DERIVATIVES) &
                  & MAX_NUMBER_OF_DERIVATIVES=BASIS%NUMBER_OF_DERIVATIVES(nn)
              ELSE
                CALL FLAG_ERROR("Could not find local node.",ERR,ERROR,*999)
              ENDIF
            ENDDO !elem_idx
            CALL LIST_REMOVE_DUPLICATES(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            CALL LIST_DETACH_AND_DESTROY(NODE_DERIVATIVE_LIST,NUMBER_OF_DERIVATIVES,DERIVATIVES,ERR,ERROR,*999)
            IF(NUMBER_OF_DERIVATIVES==MAX_NUMBER_OF_DERIVATIVES) THEN
              NODES%NODES(np)%NUMBER_OF_DERIVATIVES=MAX_NUMBER_OF_DERIVATIVES
              ALLOCATE(NODES%NODES(np)%GLOBAL_DERIVATIVE_INDEX(MAX_NUMBER_OF_DERIVATIVES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global derivative index.",ERR,ERROR,*999)
              ALLOCATE(NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(MAX_NUMBER_OF_DERIVATIVES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node partial derivative index.",ERR,ERROR,*999)
              DO nk=1,NUMBER_OF_DERIVATIVES                
                NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(nk)=DERIVATIVES(nk)
                global_deriv=PARTIAL_DERIVATIVE_GLOBAL_DERIVATIVE_MAP(DERIVATIVES(nk))
                IF(global_deriv/=0) THEN
                   NODES%NODES(np)%GLOBAL_DERIVATIVE_INDEX(nk)=global_deriv
                ELSE
                  LOCAL_ERROR="The partial derivative index of "//TRIM(NUMBER_TO_VSTRING(DERIVATIVES(nk),"*",ERR,ERROR))// &
                    & " for derivative number "//TRIM(NUMBER_TO_VSTRING(nk,"*",ERR,ERROR))// &
                    & " does not have a corresponding global derivative."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDDO !nk
              DEALLOCATE(DERIVATIVES)
            ELSE
              LOCAL_ERROR="Invalid mesh configuration. User node "// &
                & TRIM(NUMBER_TO_VSTRING(NODES%NODES(np)%USER_NUMBER,"*",ERR,ERROR))// &
                & " has inconsistent derivative directions."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF            
          ENDDO !np
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
      DO np=1,NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",NODES%NODES(np)%NUMBER_OF_DERIVATIVES, &
          & ERR,ERROR,*999)        
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES%NODES(np)%NUMBER_OF_DERIVATIVES,8,8, &
          & NODES%NODES(np)%GLOBAL_DERIVATIVE_INDEX,'("    Global derivative index(nk) :",8(X,I2))','(36X,8(X,I2))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES%NODES(np)%NUMBER_OF_DERIVATIVES,8,8, &
          & NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX,'("    Partial derivative index(nk) :",8(X,I2))','(36X,8(X,I2))', &
          & ERR,ERROR,*999)
      ENDDO !np
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE")
    RETURN
999 IF(ASSOCIATED(DERIVATIVES)) DEALLOCATE(DERIVATIVES)
    IF(ASSOCIATED(NODE_DERIVATIVE_LIST)) CALL LIST_DESTROY(NODE_DERIVATIVE_LIST,ERR,ERROR,*998)
998 CALL ERRORS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE

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
          MESH_TOPOLOGY=>MESH%TOPOLOGY(DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
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

END MODULE MESH_ROUTINES
