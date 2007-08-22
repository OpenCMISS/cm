!> \file
!> $Id: mesh_routines.f90 28 2007-07-27 08:35:14Z cpb $
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

  INTERFACE MESH_NUMBER_OF_COMPONENTS_SET
    MODULE PROCEDURE MESH_NUMBER_OF_COMPONENTS_SET_NUMBER
    MODULE PROCEDURE MESH_NUMBER_OF_COMPONENTS_SET_PTR
  END INTERFACE !MESH_NUMBER_OF_COMPONENTS_SET
  
  INTERFACE MESH_NUMBER_OF_ELEMENTS_SET
    MODULE PROCEDURE MESH_NUMBER_OF_ELEMENTS_SET_NUMBER
    MODULE PROCEDURE MESH_NUMBER_OF_ELEMENTS_SET_PTR
  END INTERFACE !MESH_NUMBER_OF_ELEMENTS_SET
  
  PUBLIC DECOMPOSITION_ALL_TYPE,DECOMPOSITION_CALCULATED_TYPE,DECOMPOSITION_USER_DEFINED_TYPE

  PUBLIC MESHES_INITIALISE,MESHES_FINALISE,MESH_CREATE_START,MESH_CREATE_FINISH,MESH_DESTROY,MESH_CREATE_REGULAR, &
    & MESH_NUMBER_OF_COMPONENTS_SET,MESH_NUMBER_OF_ELEMENTS_SET,MESH_TOPOLOGY_ELEMENTS_CREATE_START, &
    & MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH,MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET,MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  PUBLIC NODES_INITIALISE,NODES_FINALISE,NODES_CREATE_START,NODES_CREATE_FINISH

  PUBLIC DECOMPOSITIONS_INITIALISE,DECOMPOSITIONS_FINALISE,DECOMPOSITION_CREATE_START,DECOMPOSITION_CREATE_FINISH, &
    & DECOMPOSITION_DESTROY,DECOMPOSITION_USER_NUMBER_FIND,DECOMPOSITION_TYPE_SET,DECOMPOSITION_NUMBER_OF_DOMAINS_SET
  
CONTAINS

  !
  !================================================================================================================================
  !

  SUBROUTINE DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_CREATE_FINISH
    !###  Description:
    !###    Finishes the creation of a domain decomposition on a given mesh.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: decomposition_no
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        IF(ASSOCIATED(DECOMPOSITION)) THEN
          IF(DECOMPOSITION%DECOMPOSITIONS%MESH%USER_NUMBER==MESH%USER_NUMBER) THEN
            DECOMPOSITION%DECOMPOSITION_FINISHED=.TRUE.
            !Calculate which elements belong to which domain
            CALL DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*999)
            !Initialise the topology information for this decomposition
            CALL DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
            !Initialise the domain for this computational node            
            CALL DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*999)
            !Calculate the decomposition topology
            CALL DECOMPOSITION_TOPOLOGY_CALCULATE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
          ELSE
            LOCAL_ERROR="The specified decomposition was created on mesh number "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%DECOMPOSITIONS%MESH%USER_NUMBER,"*",ERR,ERROR))// &
              & " which is different from the specified mesh number of "// &
              & TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
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

  SUBROUTINE DECOMPOSITION_CREATE_START(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_CREATE_START
    !###  Description:
    !###    Starts the creation of a domain decomposition for a given mesh.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
              & " has already been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            ALLOCATE(NEW_DECOMPOSITION,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decomposition",ERR,ERROR,*999)
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
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decomposition element domain",ERR,ERROR,*999)
            NEW_DECOMPOSITION%ELEMENT_DOMAIN=0          
            !Nullify the domain
            NULLIFY(NEW_DECOMPOSITION%DOMAIN)
            !Add new decomposition into list of decompositions on the mesh
            ALLOCATE(NEW_DECOMPOSITIONS(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS+1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decompositions",ERR,ERROR,*999)
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
            & " are not associated"
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

  SUBROUTINE DECOMPOSITION_DESTROY(USER_NUMBER,MESH,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_DESTROY
    !###  Description:
    !###    Destroys a domain decomposition identified by a user number and deallocates all memory.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: decomposition_idx,decomposition_position
    LOGICAL :: FOUND    
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_PTR_TYPE), POINTER :: NEW_DECOMPOSITIONS(:)

    NULLIFY(NEW_DECOMPOSITIONS)

    CALL ENTERS("DECOMPOSITION_DESTROY",ERR,ERROR,*999)

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
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new decompositions",ERR,ERROR,*999)
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
            & " has not been created on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The decompositions on mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
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
  
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE(DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE
    !###    Calculates a domains for a decomposition of a mesh.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: number_elem_indicies,elem_index,elem_count,ne,nn,my_computational_node_number,number_computational_nodes, &
      & no_computational_node,ELEMENT_START,ELEMENT_STOP,MY_ELEMENT_START,MY_ELEMENT_STOP,NUMBER_OF_ELEMENTS, &
      & MY_NUMBER_OF_ELEMENTS,MPI_IERROR,MAX_NUMBER_ELEMENTS_PER_NODE,component_idx
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_PTR(:),ELEMENT_INDICIES(:),ELEMENT_DISTANCE(:),DISPLACEMENTS(:),RECEIVE_COUNTS(:)
    INTEGER(INTG) :: ELEMENT_WEIGHT(1),WEIGHT_FLAG,NUMBER_FLAG,NUMBER_OF_CONSTRAINTS, &
      & NUMBER_OF_COMMON_NODES,PARMETIS_OPTIONS(0:2)
    REAL(SP) :: UBVEC(1)
    REAL(SP), ALLOCATABLE :: TPWGTS(:)
    REAL(DP) :: NUMBER_ELEMENTS_PER_NODE
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH

    CALL ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        MESH=>DECOMPOSITION%MESH
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN

          component_idx=DECOMPOSITION%MESH_COMPONENT_NUMBER
          
          my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
          IF(ERR/=0) GOTO 999
          
          SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)          
          CASE(DECOMPOSITION_ALL_TYPE)
            !Do nothing
            
          CASE(DECOMPOSITION_CALCULATED_TYPE)
            !Calculate the general decomposition
            
            number_computational_nodes=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
            IF(ERR/=0) GOTO 999
            
            NUMBER_ELEMENTS_PER_NODE=REAL(MESH%NUMBER_OF_ELEMENTS,DP)/REAL(number_computational_nodes,DP)
            ELEMENT_START=1
            ELEMENT_STOP=0
            MAX_NUMBER_ELEMENTS_PER_NODE=-1
            ALLOCATE(RECEIVE_COUNTS(0:number_computational_nodes-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate recieve counts",ERR,ERROR,*999)
            ALLOCATE(DISPLACEMENTS(0:number_computational_nodes-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate displacements",ERR,ERROR,*999)
            ALLOCATE(ELEMENT_DISTANCE(0:number_computational_nodes),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element distance",ERR,ERROR,*999)
            ELEMENT_DISTANCE(0)=0
            DO no_computational_node=0,number_computational_nodes-1
              ELEMENT_START=ELEMENT_STOP+1
              IF(no_computational_node==number_computational_nodes-1) THEN
                ELEMENT_STOP=MESH%NUMBER_OF_ELEMENTS
              ELSE
                ELEMENT_STOP=ELEMENT_START+NINT(NUMBER_ELEMENTS_PER_NODE,INTG)-1
              ENDIF
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
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element pointer list",ERR,ERROR,*999)
            ALLOCATE(ELEMENT_INDICIES(0:number_elem_indicies-1),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element indicies list",ERR,ERROR,*999)
            ALLOCATE(TPWGTS(1:DECOMPOSITION%NUMBER_OF_DOMAINS),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate tpwgts",ERR,ERROR,*999)
            elem_index=0
            elem_count=0
            ELEMENT_PTR(0)=0
            DO ne=MY_ELEMENT_START,MY_ELEMENT_STOP
              elem_count=elem_count+1
              BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
              DO nn=1,BASIS%NUMBER_OF_NODES
                ELEMENT_INDICIES(elem_index)=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)% &
                  & GLOBAL_ELEMENT_NODES(nn)-1 !C numbering
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
              & MPI_COMM_WORLD,ERR,ERROR,*999)

            !Transfer all the element domain information to the other computational nodes so that each rank has all the information
            IF(number_computational_nodes>1) THEN
              !This should work on a single processor but doesn't for mpich2 under windows. Maybe a bug? Avoid for now.
              CALL MPI_ALLGATHERV(MPI_IN_PLACE,MAX_NUMBER_ELEMENTS_PER_NODE,MPI_INTEGER,DECOMPOSITION%ELEMENT_DOMAIN, &
                & RECEIVE_COUNTS,DISPLACEMENTS,MPI_INTEGER,MPI_COMM_WORLD,MPI_IERROR)
              CALL MPI_ERROR_CHECK("MPI_ALLGATHERV",MPI_IERROR,ERR,ERROR,*999)
            ENDIF
            
            DEALLOCATE(DISPLACEMENTS)
            DEALLOCATE(RECEIVE_COUNTS)
            DEALLOCATE(ELEMENT_DISTANCE)
            DEALLOCATE(ELEMENT_PTR)
            DEALLOCATE(ELEMENT_INDICIES)
            DEALLOCATE(TPWGTS)
            
          CASE(DECOMPOSITION_USER_DEFINED_TYPE)
!!TODO: Check decomposition setup.
          
          CASE DEFAULT
            CALL FLAG_ERROR("Invalid domain decomposition type",ERR,ERROR,*999)
            
          END SELECT
        ELSE
          CALL FLAG_ERROR("Decomposition mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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
  
  SUBROUTINE DECOMPOSITION_ELEMENT_DOMAIN_SET(DECOMPOSITION,GLOBAL_ELEMENT_NUMBER,DOMAIN_NUMBER,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_ELEMENT_DOMAIN_SET
    !###    Sets the domain for a given element in a decomposition of a mesh.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(IN) :: GLOBAL_ELEMENT_NUMBER
    INTEGER(INTG), INTENT(IN) :: DOMAIN_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: number_computational_nodes
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: MESH_TOPOLOGY
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_ELEMENT_DOMAIN_SET",ERR,ERROR,*999)

!!TODO: interface should specify user element number ???
    
    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished",ERR,ERROR,*999)
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
                  & " is invalid. The limits are 0 to "//TRIM(NUMBER_TO_VSTRING(number_computational_nodes,"*",ERR,ERROR))
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_ELEMENT_NUMBER,"*",ERR,ERROR))// &
                & " is invalid. The limits are 1 to "// &
                & TRIM(NUMBER_TO_VSTRING(MESH_TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Decomposition mesh topology is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_MESH_COMPONENT_NUMBER_SET(DECOMPOSITION,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_MESH_COMPONENT_NUMBER_SET
    !###  Description:
    !###    Sets/changes the mesh component number which will be used for the decomposition of a mesh.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_MESH_COMPONENT_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN     
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished",ERR,ERROR,*999)
      ELSE
        IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
          IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS) THEN
            DECOMPOSITION%MESH_COMPONENT_NUMBER=MESH_COMPONENT_NUMBER
          ELSE
            LOCAL_ERROR="The speficied mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
              & "is invalid. The component number must be between 1 and "// &
              & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_NUMBER_OF_DOMAINS_SET
    !###  Description:
    !###    Sets/changes the number of domains for a decomposition.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_NUMBER_OF_DOMAINS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(DECOMPOSITION%DECOMPOSITION_TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          IF(NUMBER_OF_DOMAINS==1) THEN
            DECOMPOSITION%NUMBER_OF_DOMAINS=1
          ELSE
            CALL FLAG_ERROR("Can only have one domain for all decomposition type",ERR,ERROR,*999)
          ENDIF
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          IF(NUMBER_OF_DOMAINS>=1) THEN
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
                & TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%MESH%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))//") in the mesh"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
             CALL FLAG_ERROR("Number of domains must be >= 1",ERR,ERROR,*999)
           ENDIF
         CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(DECOMPOSITION%DECOMPOSITION_TYPE,"*",ERR,ERROR))// &
            & " is not valid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_CALCULATE
    !###  Description:
    !###    Calculates the decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      !Calculate the elements topology
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the line topology
      CALL DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
      !Calculate the face topology
      !CALL DECOMPOSITION_TOPOLOGY_FACES_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_ELEMENT_FINALISE
    !###  Description:
    !###    Finalises the given decomposition topology element.

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE
    !###  Description:
    !###    Initialises the given decomposition topology element.

    !Argument variables
    TYPE(DECOMPOSITION_ELEMENT_TYPE) :: ELEMENT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%LOCAL_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
  
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENT_INITALISE")
    RETURN 1   
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)
  SUBROUTINE DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE
    !###  Description:
    !###    Calculates the element numbers adjacent to an element in a decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,ne2,nep1,nep2,ni,nn,nn1,nn2,np,np1,np2,DUMMY_ERR,FACE_XI(2),NODE_POSITION_INDEX(3)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: NUMBER_SURROUNDING,MAX_NUMBER_SURROUNDING,NUMBER_OF_NODES_XI(3)
    INTEGER(INTG), POINTER :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:)
    LOGICAL :: FOUND,XI_COLLAPSED,FACE_COLLAPSED(-3:3)
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-3:3)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DECOMPOSITION_ELEMENTS_TYPE), POINTER :: DECOMPOSITION_ELEMENTS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY

    NULLIFY(NODE_MATCHES)
    NULLIFY(ADJACENT_ELEMENTS)

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR,*999)
    
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
                    DO ni=-BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
                      CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
                      CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                      CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,5,ERR,ERROR,*999)
                      CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
                    ENDDO !ni
                    NUMBER_OF_NODES_XI=1
                    DO ni=1,BASIS%NUMBER_OF_XI
                      NUMBER_OF_NODES_XI(ni)=BASIS%NUMBER_OF_NODES_XI(ni)
                    ENDDO !ni
                    !Place the current element in the surrounding list
                    CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER, &
                      & ERR,ERROR,*999)
                    MAX_NUMBER_SURROUNDING=1
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
                            DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XI(xi_dir_search).AND.XI_COLLAPSED)
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
                                  & DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn2)) XI_COLLAPSED=.TRUE.
                              ENDIF
                              NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                            ENDDO !xi_dir_search
                            IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                          ENDIF
                        ENDDO !j
                        NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XI(ni)
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
                          NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XI(ni)
                        ENDIF
                        !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is
                        !also collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the
                        !central axis back to itself
                        IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                          !Do nothing - the match lists are already empty
                        ELSE
                          !Find the nodes to match and add them to the node match list
                          DO nn1=1,NUMBER_OF_NODES_XI(FACE_XI(1))
                            NODE_POSITION_INDEX(FACE_XI(1))=nn1
                            DO nn2=1,NUMBER_OF_NODES_XI(FACE_XI(2))
                              NODE_POSITION_INDEX(FACE_XI(2))=nn2
                              nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2), &
                                & NODE_POSITION_INDEX(3),1)
                              IF(nn/=0) THEN
                                np=DOMAIN_ELEMENTS%ELEMENTS(ne)%ELEMENT_NODES(nn)
                                CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                              ENDIF
                            ENDDO !nn2
                          ENDDO !nn1
                        ENDIF
                        CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                        CALL LIST_DETACH(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                        NUMBER_SURROUNDING=0
                        IF(NUMBER_NODE_MATCHES>0) THEN
                          !Find list of elements surrounding those nodes
                          np1=NODE_MATCHES(1)
                          DO nep1=1,DOMAIN_NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                            ne1=DOMAIN_NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                            IF(ne1/=ne) THEN !Don't want the current element
                              FOUND=.FALSE.
                              nn2=2
                              DO WHILE(nn2<=NUMBER_NODE_MATCHES.AND..NOT.FOUND)
                                np2=NODE_MATCHES(nn2)
                                nep2=1
                                DO WHILE(nep2<=DOMAIN_NODES%NODES(np2)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.FOUND)
                                  ne2=DOMAIN_NODES%NODES(np2)%SURROUNDING_ELEMENTS(nep2)
                                  IF(ne1==ne2) THEN
                                    FOUND=.TRUE.
                                  ELSE
                                    nep2=nep2+1
                                  ENDIF
                                ENDDO !nep2
                                nn2=nn2+1
                              ENDDO !nn2
                              IF(FOUND) THEN
                                CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                                NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                              ENDIF
                            ENDIF
                          ENDDO !nep1
                        ENDIF
                        IF(ASSOCIATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                        CALL LIST_DESTROY(NODE_MATCH_LIST,ERR,ERROR,*999)
                        IF(NUMBER_SURROUNDING>MAX_NUMBER_SURROUNDING) MAX_NUMBER_SURROUNDING=NUMBER_SURROUNDING
                      ENDDO !direction_index
                    ENDDO !ni
                    !Set the surrounding elements for this element
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS( &
                      & -BASIS%NUMBER_OF_XI:BASIS%NUMBER_OF_XI),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of surrounding elements",ERR,ERROR,*999)
                    ALLOCATE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(MAX_NUMBER_SURROUNDING, &
                      & -BASIS%NUMBER_OF_XI:BASIS%NUMBER_OF_XI),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate surrounding elements",ERR,ERROR,*999)
                    DO ni=-BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
                      CALL LIST_DETACH(ADJACENT_ELEMENTS_LIST(ni)%PTR,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & NUMBER_OF_ADJACENT_ELEMENTS(ni),ADJACENT_ELEMENTS,ERR,ERROR,*999)
                      DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
                        & NUMBER_OF_ADJACENT_ELEMENTS(ni),ni)=ADJACENT_ELEMENTS(1:DECOMPOSITION_ELEMENTS% &
                        & ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni))
                      IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
                      CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
                    ENDDO !ni            
                  ENDDO !ne           
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
        ELSE
          CALL FLAG_ERROR("Topology elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology decomposition is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not allocated",ERR,ERROR,*999)
    ENDIF

    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Total number of elements = ",DECOMPOSITION_ELEMENTS% &
        & TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
        BASIS=>DOMAIN_ELEMENTS%ELEMENTS(ne)%BASIS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Local element number = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi directions = ",BASIS%NUMBER_OF_XI,ERR,ERROR,*999)
        DO ni=-BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Xi direction = ",ni,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"        Number of adjacent elements = ", &
            & DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni),ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)% &
            & NUMBER_OF_ADJACENT_ELEMENTS(ni),8,8,DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(:,ni), &
            & '("        Adjacent elements =",8(X,I6))','(30x,8(X,I6))',ERR,ERROR,*999)
        ENDDO !ni
      ENDDO !ne
    ENDIF
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE")
    RETURN
999 IF(ASSOCIATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO ni=-3,3
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(ni)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 CALL ERRORS("DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE")
    RETURN 1   
  !END SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE
  END SUBROUTINE DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE
  
  !
  !================================================================================================================================
  !

  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_ELEMENTS_CALCULATE
    !###  Description:
    !###    Calculates the decomposition element topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
                          IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition elements elements",ERR,ERROR,*999)
                          DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                          DO ne=1,DECOMPOSITION_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                            CALL DECOMPOSITION_TOPOLOGY_ELEMENT_INITIALISE(DECOMPOSITION_ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                            global_element=DOMAIN_ELEMENTS_MAPPING%LOCAL_TO_GLOBAL_MAP(ne)
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%USER_NUMBER=MESH_ELEMENTS%ELEMENTS(global_element)%USER_NUMBER
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%LOCAL_NUMBER=ne
                            DECOMPOSITION_ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=global_element
                          ENDDO !ne
                          !Calculate the elements surrounding the elements in the decomposition topology
                          !CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
                          CALL DECOMP_TOPOLOGY_ELEM_ADJACENT_ELEM_CALCULATE(TOPOLOGY,ERR,ERROR,*999)
                        ELSE
                          CALL FLAG_ERROR("Mesh elements is not associated",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FLAG_ERROR("Domain mappings elements is not associated",ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Domain mappings is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain topology elements is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Topology decomposition domain topology is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Topology decomposition domain is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Topology decomposition is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE
    !###  Description:
    !###    Finalises the elements in the given decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE
    !###  Description:
    !###    Initialises the element data structures for a decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        CALL FLAG_ERROR("Decomposition already has topology elements associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology elements",ERR,ERROR,*999)
        TOPOLOGY%ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS=0
        TOPOLOGY%ELEMENTS%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
        NULLIFY(TOPOLOGY%ELEMENTS%ELEMENTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_FINALISE
    !###  Description:
    !###    Finalises the topology in the given decomposition.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      CALL DECOMPOSITION_TOPOLOGY_LINES_FINALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      DEALLOCATE(DECOMPOSITION%TOPOLOGY)
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_INITIALISE
    !###  Description:
    !###    Initialises the topology for a given decomposition.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%TOPOLOGY)) THEN
        CALL FLAG_ERROR("Decomposition already has topology associated",ERR,ERROR,*999)
      ELSE
        !Allocate decomposition topology
        ALLOCATE(DECOMPOSITION%TOPOLOGY,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Decomposition topology could not be allocated",ERR,ERROR,*999)
        DECOMPOSITION%TOPOLOGY%DECOMPOSITION=>DECOMPOSITION
        NULLIFY(DECOMPOSITION%TOPOLOGY%ELEMENTS)
        NULLIFY(DECOMPOSITION%TOPOLOGY%LINES)
        !Initialise the topology components
        CALL DECOMPOSITION_TOPOLOGY_ELEMENTS_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
        CALL DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(DECOMPOSITION%TOPOLOGY,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_LINE_FINALISE
    !###  Description:
    !###    Finalises a line in the given decomposition topology and deallocates all memory.

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
    !###  Description:
    !###    Initialises the line data structure for a decomposition topology line.

    !Argument variables
    TYPE(DECOMPOSITION_LINE_TYPE) :: LINE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR,*999)

    LINE%NUMBER=0
    LINE%XI_DIRECTION=0
    LINE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    LINE%ADJACENT_LINES=0
    
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN
999 CALL ERRORS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE",ERR,ERROR)
    CALL EXITS("DECOMPOSITION_TOPOLOGY_LINE_INITIALISE")
    RETURN 1
  END SUBROUTINE DECOMPOSITION_TOPOLOGY_LINE_INITIALISE
  
  !
  !================================================================================================================================
  !

  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_LINES_CALCULATE
    !###  Description:
    !###    Calculates the lines in the given decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
                      CALL FLAG_ERROR("Invalid number of dimensions for a topology domain",ERR,ERROR,*999)
                    END SELECT
                    DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES
                    IF(ASSOCIATED(DOMAIN_LINES)) THEN
                      ALLOCATE(TEMP_LINES(4,MAX_NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary lines array",ERR,ERROR,*999)
                      ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes number of lines array",ERR,ERROR,*999)
                      NODES_NUMBER_OF_LINES=0
                      NUMBER_OF_LINES=0
                      TEMP_LINES=0
                      !Loop over the elements in the topology
                      DO ne=1,DOMAIN_ELEMENTS%TOTAL_NUMBER_OF_ELEMENTS
                        DOMAIN_ELEMENT=>DOMAIN_ELEMENTS%ELEMENTS(ne)
                        DECOMPOSITION_ELEMENT=>DECOMPOSITION_ELEMENTS%ELEMENTS(ne)
                        BASIS=>DOMAIN_ELEMENT%BASIS
                        ALLOCATE(DECOMPOSITION_ELEMENT%ELEMENT_LINES(BASIS%NUMBER_OF_LOCAL_LINES),STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element element lines",ERR,ERROR,*999)
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
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new number of lines",ERR,ERROR,*999)
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
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node lines array",ERR,ERROR,*999)
                        DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES=0
                      ENDDO !np
                      DEALLOCATE(NODES_NUMBER_OF_LINES)
                      ALLOCATE(DECOMPOSITION_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate decomposition topology lines",ERR,ERROR,*999)
                      DECOMPOSITION_LINES%NUMBER_OF_LINES=NUMBER_OF_LINES
                      ALLOCATE(DOMAIN_LINES%LINES(NUMBER_OF_LINES),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain topology lines",ERR,ERROR,*999)
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
                            DECOMPOSITION_LINE%XI_DIRECTION=BASIS%LOCAL_LINE_XI_DIRECTION(nae)
                            DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                            ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line nodes in line",ERR,ERROR,*999)
                            ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                              & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line derivatives in line",ERR,ERROR,*999)
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
                        !Allocate the elements surrounding the line
                        ALLOCATE(DECOMPOSITION_LINE%SURROUNDING_ELEMENTS(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line surrounding elements",ERR,ERROR,*999)
                        ALLOCATE(DECOMPOSITION_LINE%ELEMENT_LINES(DECOMPOSITION_LINE%NUMBER_OF_SURROUNDING_ELEMENTS), &
                          & STAT=ERR)
                        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate line element lines",ERR,ERROR,*999)
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
                      CALL FLAG_ERROR("Domain topology lines is not associated",ERR,ERROR,*999)
                    ENDIF
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
            !Now loop over the other mesh components in the decomposition and calculate the domain lines
            MESH=>DECOMPOSITION%MESH
            IF(ASSOCIATED(MESH)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(component_idx/=DECOMPOSITION%MESH_COMPONENT_NUMBER) THEN
                  DOMAIN=>DECOMPOSITION%DOMAIN(component_idx)%PTR
                  IF(ASSOCIATED(DOMAIN)) THEN
                    DOMAIN_TOPOLOGY=DOMAIN%TOPOLOGY
                    IF(ASSOCIATED(DOMAIN_TOPOLOGY)) THEN
                      DOMAIN_NODES=>DOMAIN_TOPOLOGY%NODES
                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                        DOMAIN_ELEMENTS=>DOMAIN_TOPOLOGY%ELEMENTS
                        IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                          DOMAIN_LINES=>DOMAIN_TOPOLOGY%LINES                      
                          IF(ASSOCIATED(DOMAIN_LINES)) THEN
                            ALLOCATE(DOMAIN_LINES%LINES(DECOMPOSITION_LINES%NUMBER_OF_LINES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain lines lines",ERR,ERROR,*999)
                            DOMAIN_LINES%NUMBER_OF_LINES=DECOMPOSITION_LINES%NUMBER_OF_LINES
                            ALLOCATE(NODES_NUMBER_OF_LINES(DOMAIN_NODES%TOTAL_NUMBER_OF_NODES),STAT=ERR)
                            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes number of lines array",ERR,ERROR,*999)
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
                                DOMAIN_LINE%BASIS=>BASIS%LINE_BASES(DECOMPOSITION_LINE%XI_DIRECTION)%PTR
                                ALLOCATE(DOMAIN_LINE%NODES_IN_LINE(BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate nodes in line",ERR,ERROR,*999)
                                ALLOCATE(DOMAIN_LINE%DERIVATIVES_IN_LINE(DOMAIN_LINE%BASIS%MAXIMUM_NUMBER_OF_DERIVATIVES, &
                                  & BASIS%NUMBER_OF_NODES_IN_LOCAL_LINE(nae)),STAT=ERR)
                                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate derivatives in line",ERR,ERROR,*999)
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
                                CALL FLAG_ERROR("Line is not surrounded by any elements",ERR,ERROR,*999)
                              ENDIF                              
                            ENDDO !nl
                            DO np=1,DOMAIN_NODES%TOTAL_NUMBER_OF_NODES
                              ALLOCATE(DOMAIN_NODES%NODES(np)%NODE_LINES(NODES_NUMBER_OF_LINES(np)),STAT=ERR)
                              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node lines",ERR,ERROR,*999)
                              DOMAIN_NODES%NODES(np)%NUMBER_OF_NODE_LINES=0
                            ENDDO !np
                            DEALLOCATE(NODES_NUMBER_OF_LINES)
                            DO nl=1,DOMAIN_LINES%NUMBER_OF_LINES
                              DOMAIN_LINE=>DOMAIN_LINES%LINES(nl)
                              BASIS=>DOMAIN_LINE%BASIS
                              DO nnl=1,BASIS%NUMBER_OF_NODES
                                np=DOMAIN_LINE%NODES_IN_LINE(nl)
                                DOMAIN_NODE=>DOMAIN_NODES%NODES(np)
                                DOMAIN_NODE%NUMBER_OF_NODE_LINES=DOMAIN_NODE%NUMBER_OF_NODE_LINES+1
                                DOMAIN_NODE%NODE_LINES(DOMAIN_NODE%NUMBER_OF_NODE_LINES)=nl
                              ENDDO !nnl
                            ENDDO !nl
                          ELSE
                            CALL FLAG_ERROR("Domain lines is not associated",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Topology lines is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_LINES_FINALISE
    !###  Description:
    !###    Finalises the lines in the given decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TOPOLOGY_LINES_INITIALISE
    !###  Description:
    !###    Initialises the line data structures for a decomposition topology.

    !Argument variables
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITION_TOPOLOGY_LINES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%LINES)) THEN
        CALL FLAG_ERROR("Decomposition already has topology lines associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(TOPOLOGY%LINES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate topology lines",ERR,ERROR,*999)
        TOPOLOGY%LINES%NUMBER_OF_LINES=0
        TOPOLOGY%LINES%DECOMPOSITION=>TOPOLOGY%DECOMPOSITION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_TYPE_SET(DECOMPOSITION,TYPE,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_TYPE_SET
    !###  Description:
    !###    Sets/changes the decomposition type for a decomposition.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(IN) :: TYPE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DECOMPOSITION_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(DECOMPOSITION%DECOMPOSITION_FINISHED) THEN
        CALL FLAG_ERROR("Decomposition has been finished",ERR,ERROR,*999)
      ELSE
        SELECT CASE(TYPE)
        CASE(DECOMPOSITION_ALL_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_ALL_TYPE
        CASE(DECOMPOSITION_CALCULATED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_CALCULATED_TYPE
        CASE(DECOMPOSITION_USER_DEFINED_TYPE)
          DECOMPOSITION%DECOMPOSITION_TYPE=DECOMPOSITION_USER_DEFINED_TYPE
        CASE DEFAULT
          LOCAL_ERROR="Decomposition type "//TRIM(NUMBER_TO_VSTRING(TYPE,"*",ERR,ERROR))//" is not valid"
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITION_USER_NUMBER_FIND(USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITION_USER_NUMBER_FIND
    !###  Description:
    !###    Finds and returns in DECOMPOSITION a pointer to the decomposition identified by USER_NUMBER in the given MESH.
    !###    If no decomposition with that USER_NUMBER exists DECOMPOSITION is left nullified.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITIONS_FINALISE
    !###  Description:
    !###    Finalises the domain decompositions for a given mesh.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: USER_NUMBER

    CALL ENTERS("DECOMPOSITIONS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        DO WHILE(MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS>0)
          USER_NUMBER=MESH%DECOMPOSITIONS%DECOMPOSITIONS(1)%PTR%USER_NUMBER
          CALL DECOMPOSITION_DESTROY(USER_NUMBER,MESH,ERR,ERROR,*999)
        ENDDO !no_decomposition
       DEALLOCATE(MESH%DECOMPOSITIONS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DECOMPOSITIONS_INITIALISE(MESH,ERR,ERROR,*)

    !#### Subroutine: DECOMPOSITIONS_INITIALISE
    !###  Description:
    !###    Initialises the domain decompositions for a given mesh.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DECOMPOSITIONS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(ASSOCIATED(MESH%DECOMPOSITIONS)) THEN
        CALL FLAG_ERROR("Mesh already has decompositions associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(MESH%DECOMPOSITIONS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Mesh decompositions could not be allocated",ERR,ERROR,*999)
        MESH%DECOMPOSITIONS%NUMBER_OF_DECOMPOSITIONS=0
        NULLIFY(MESH%DECOMPOSITIONS%DECOMPOSITIONS)
        MESH%DECOMPOSITIONS%MESH=>MESH
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_FINALISE(DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_FINALISE
    !###  Description:
    !###    Finalises the domain for a given decomposition and deallocates all memory.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_INITIALISE(DECOMPOSITION,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_INITIALISE
    !###  Description:
    !###    Initialises the domain for a given decomposition.

    !Argument variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: component_idx

    CALL ENTERS("DOMAIN_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DECOMPOSITION)) THEN
      IF(ASSOCIATED(DECOMPOSITION%MESH)) THEN
        IF(ASSOCIATED(DECOMPOSITION%DOMAIN)) THEN
          CALL FLAG_ERROR("Decomposition already has a domain associated",ERR,ERROR,*999)
        ELSE
          ALLOCATE(DECOMPOSITION%DOMAIN(DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS),STAT=ERR)
          IF(ERR/=0) CALL FLAG_ERROR("Decomposition domain could not be allocated",ERR,ERROR,*999)
          DO component_idx=1,DECOMPOSITION%MESH%NUMBER_OF_COMPONENTS
            ALLOCATE(DECOMPOSITION%DOMAIN(component_idx)%PTR,STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Decomposition domain component could not be allocated",ERR,ERROR,*999)
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
        CALL FLAG_ERROR("Decomposition mesh is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Decomposition is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_DOFS_FINALISE
    !###  Description:
    !###    Finalises the dofs mapping in the given domain mappings.

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_DOFS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%DOFS,ERR,ERROR,*999)
        DEALLOCATE(DOMAIN_MAPPINGS%DOFS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_MAPPINGS_DOFS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_DOFS_INITIALISE
    !###  Description:
    !###    Intialises the dofs mapping in the given domain mapping.

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_DOFS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%DOFS)) THEN
        CALL FLAG_ERROR("Domain dofs mappings are already associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%DOFS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings dofs",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%DOFS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
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
  
  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
    !###    Calculates the local/global element mappings for a domain decomposition.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: no_adjacent_element,adjacent_element,domain_no,domain_idx,ne,nn,np,NUMBER_OF_DOMAINS, &
      & NUMBER_OF_ADJACENT_ELEMENTS,my_computational_node_number,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_ELEMENT_NUMBERS(:)
    INTEGER(INTG), POINTER :: DOMAINS(:),ADJACENT_ELEMENTS(:)
    TYPE(LIST_TYPE), POINTER :: ADJACENT_DOMAINS_LIST
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: ADJACENT_ELEMENTS_LIST(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: ELEMENTS_MAPPING

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
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element mapping global to local map",ERR,ERROR,*999)
              ELEMENTS_MAPPING%NUMBER_OF_GLOBAL=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS
              !Loop over the global elements and calculate local numbers
              ALLOCATE(LOCAL_ELEMENT_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local element numbers",ERR,ERROR,*999)
              LOCAL_ELEMENT_NUMBERS=0
              ALLOCATE(ADJACENT_ELEMENTS_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent elements list",ERR,ERROR,*999)
              DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                NULLIFY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR)
                CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,INT(MESH%NUMBER_OF_ELEMENTS/2),ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
              ENDDO !domain_idx
            
              DO ne=1,MESH%NUMBER_OF_ELEMENTS
                !Calculate the local numbers
                domain_no=DECOMPOSITION%ELEMENT_DOMAIN(ne)
                LOCAL_ELEMENT_NUMBERS(domain_no)=LOCAL_ELEMENT_NUMBERS(domain_no)+1
                !Calculate the adjacent elements to the computational domains and the adjacent domain numbers themselves
                BASIS=>MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS
                CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,domain_no,ERR,ERROR,*999)
                DO nn=1,BASIS%NUMBER_OF_NODES
                  np=MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)
                  DO no_adjacent_element=1,MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%NUMBER_OF_SURROUNDING_ELEMENTS
                    adjacent_element=MESH%TOPOLOGY(component_idx)%PTR%NODES%NODES(np)%SURROUNDING_ELEMENTS(no_adjacent_element)
                    IF(DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element)/=domain_no) THEN
                      CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(domain_no)%PTR,adjacent_element,ERR,ERROR,*999)
                      CALL LIST_ITEM_ADD(ADJACENT_DOMAINS_LIST,DECOMPOSITION%ELEMENT_DOMAIN(adjacent_element),ERR,ERROR,*999)
                    ENDIF
                  ENDDO !no_adjacent_element
                ENDDO !nn
                CALL LIST_REMOVE_DUPLICATES(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                CALL LIST_DESTROY(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                DEALLOCATE(DOMAINS)
                CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne),ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element global to local map local number",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%DOMAIN_NUMBER(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element global to local map domain number",ERR,ERROR,*999)
                ALLOCATE(ELEMENTS_MAPPING%GLOBAL_TO_LOCAL_MAP(ne)%LOCAL_TYPE(NUMBER_OF_DOMAINS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element global to local map local type",ERR,ERROR,*999)
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
                CALL LIST_DETACH(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,NUMBER_OF_ADJACENT_ELEMENTS,ADJACENT_ELEMENTS, &
                  & ERR,ERROR,*999)
                CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(domain_idx)%PTR,ERR,ERROR,*999)
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
              CALL FLAG_ERROR("Domain mesh is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Domain decomposition is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain mappings elements is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain mappings is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*998)
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
          & ELEMENTS_MAPPING%INTERNAL_LIST,'("    Internal elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary elements = ", &
          & ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & ELEMENTs_MAPPING%BOUNDARY_LIST,'("    Boundary elements:",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost elements :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost elements = ", &
          & ELEMENTs_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,ELEMENTS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & ELEMENTS_MAPPING%GHOST_LIST,'("    Ghost elements   :",8(X,I7))','(22X,8(X,I7))',ERR,ERROR,*999)
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
    IF(ASSOCIATED(DOMAIN%MAPPINGS%ELEMENTS)) CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*998)
998 CALL ERRORS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_ELEMENTS_CALCULATE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  SUBROUTINE DOMAIN_MAPPINGS_FINALISE(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: 
    !###  Description:
    !###    Finalises the mappings in the given domain.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*999)
      DEALLOCATE(DOMAIN%MAPPINGS)
    ELSE
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_ELEMENTS_FINALISE
    !###  Description:
    !###    Finalises the element mapping in the given domain mapping.

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%ELEMENTS,ERR,ERROR,*999)
        DEALLOCATE(DOMAIN_MAPPINGS%ELEMENTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_MAPPINGS_ELEMENTS_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_ELEMENTS_INITIALISE
    !###  Description:
    !###    Intialises the element mapping in the given domain mapping.

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_ELEMENTS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%ELEMENTS)) THEN
        CALL FLAG_ERROR("Domain elements mappings are already associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%ELEMENTS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings elements",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%ELEMENTS,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_MAPPINGS_INITIALISE(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_INITIALISE
    !###  Description:
    !###    Initialises the mappings for a domain decomposition.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      IF(ASSOCIATED(DOMAIN%MAPPINGS)) THEN
        CALL FLAG_ERROR("Domain already has mappings associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN%MAPPINGS,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings",ERR,ERROR,*999)
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
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*999)
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
  
  SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE
    !###    Calculates the local/global node and dof mappings for a domain decomposition.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: no_adjacent_element,no_ghost_node,adjacent_element,ghost_node,NUMBER_OF_NODES_PER_DOMAIN, &
      & domain_idx,domain_idx2,domain_no,nk,np,ny,NUMBER_OF_DOMAINS,MAX_NUMBER_DOMAINS,NUMBER_OF_GHOST_NODES, &
      & my_computational_node_number,component_idx
    INTEGER(INTG), ALLOCATABLE :: LOCAL_NODE_NUMBERS(:),LOCAL_DOF_NUMBERS(:),NUMBER_INTERNAL_NODES(:),NUMBER_BOUNDARY_NODES(:)
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
                  
                  my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
                  IF(ERR/=0) GOTO 999
                  
                  !Calculate the local and global numbers and set up the mappings
                  ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node mapping global to local map",ERR,ERROR,*999)
                  NODES_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                  ALLOCATE(LOCAL_NODE_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local node numbers",ERR,ERROR,*999)
                  LOCAL_NODE_NUMBERS=0
                  ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(MESH_TOPOLOGY%DOFS%NUMBER_OF_DOFS),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dofs mapping global to local map",ERR,ERROR,*999)
                  DOFS_MAPPING%NUMBER_OF_GLOBAL=MESH_TOPOLOGY%DOFS%NUMBER_OF_DOFS
                  ALLOCATE(LOCAL_DOF_NUMBERS(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local dof numbers",ERR,ERROR,*999)
                  LOCAL_DOF_NUMBERS=0
                  ALLOCATE(GHOST_NODES_LIST(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ghost nodes list",ERR,ERROR,*999)
                  DO domain_idx=0,DECOMPOSITION%NUMBER_OF_DOMAINS-1
                    NULLIFY(GHOST_NODES_LIST(domain_idx)%PTR)
                    CALL LIST_CREATE_START(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(GHOST_NODES_LIST(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(GHOST_NODES_LIST(domain_idx)%PTR,INT(MESH_TOPOLOGY%NODES%NUMBER_OF_NODES/2), &
                      & ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
                  ENDDO !domain_idx
                  ALLOCATE(NUMBER_INTERNAL_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of internal nodes",ERR,ERROR,*999)
                  NUMBER_INTERNAL_NODES=0
                  ALLOCATE(NUMBER_BOUNDARY_NODES(0:DECOMPOSITION%NUMBER_OF_DOMAINS-1),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of boundary nodes",ERR,ERROR,*999)
                  NUMBER_BOUNDARY_NODES=0

                  !For the first pass just determine the internal and boundary nodes
                  DO np=1,MESH_TOPOLOGY%NODES%NUMBER_OF_NODES
                    CALL LIST_CREATE_START(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DATA_TYPE_SET(ADJACENT_DOMAINS_LIST,LIST_INTG_TYPE,ERR,ERROR,*999)
                    CALL LIST_INITIAL_SIZE_SET(ADJACENT_DOMAINS_LIST,DECOMPOSITION%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_CREATE_FINISH(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
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
                    CALL LIST_DETACH(ADJACENT_DOMAINS_LIST,NUMBER_OF_DOMAINS,DOMAINS,ERR,ERROR,*999)
                    CALL LIST_DESTROY(ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_REMOVE_DUPLICATES(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
                    CALL LIST_DETACH(ALL_ADJACENT_DOMAINS_LIST,MAX_NUMBER_DOMAINS,ALL_DOMAINS,ERR,ERROR,*999)
                    CALL LIST_DESTROY(ALL_ADJACENT_DOMAINS_LIST,ERR,ERROR,*999)
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
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map local number",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map domain number",ERR,ERROR,*999)
                    ALLOCATE(NODES_MAPPING%GLOBAL_TO_LOCAL_MAP(np)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node global to local map local type",ERR,ERROR,*999)
                    DO nk=1,MESH_TOPOLOGY%NODES%NODES(np)%NUMBER_OF_DERIVATIVES
                      ny=MESH_TOPOLOGY%NODES%NODES(np)%DOF_INDEX(nk)
                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local number",ERR,ERROR,*999)
                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%DOMAIN_NUMBER(MAX_NUMBER_DOMAINS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map domain number",ERR,ERROR,*999)
                      ALLOCATE(DOFS_MAPPING%GLOBAL_TO_LOCAL_MAP(ny)%LOCAL_TYPE(MAX_NUMBER_DOMAINS),STAT=ERR)
                      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate dof global to local map local type",ERR,ERROR,*999)
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
                    CALL LIST_DETACH(GHOST_NODES_LIST(domain_idx)%PTR,NUMBER_OF_GHOST_NODES,GHOST_NODES,ERR,ERROR,*999)
                    CALL LIST_DESTROY(GHOST_NODES_LIST(domain_idx)%PTR,ERR,ERROR,*999)
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
                  
                  DEALLOCATE(GHOST_NODES_LIST)
                  DEALLOCATE(LOCAL_NODE_NUMBERS)
                  
                  !Calculate node and dof local to global maps from global to local map
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(NODES_MAPPING,ERR,ERROR,*999)
                  CALL DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOFS_MAPPING,ERR,ERROR,*999)
                  
                ELSE
                  CALL FLAG_ERROR("Domain mesh is not associated",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("Domain decomposition is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Domain mappings elements is not associated",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Domain mappings dofs is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Domain mappings nodes is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Domain mappings is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain is not associated",ERR,ERROR,*998)
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
          & NODES_MAPPING%INTERNAL_LIST,'("    Internal nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary nodes = ", &
          & NODES_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & NODES_MAPPING%BOUNDARY_LIST,'("    Boundary nodes:",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost nodes :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost nodes = ", &
          & NODES_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NODES_MAPPING%NUMBER_OF_GHOST,8,8, &
          & NODES_MAPPING%GHOST_LIST,'("    Ghost nodes   :",8(X,I7))','(19X,8(X,I7))',ERR,ERROR,*999)
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
          & DOFS_MAPPING%INTERNAL_LIST,'("    Internal dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Boundary dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_BOUNDARY,8,8, &
          & DOFS_MAPPING%BOUNDARY_LIST,'("    Boundary dofs:",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Ghost dofs :",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost dofs = ", &
          & DOFS_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOFS_MAPPING%NUMBER_OF_GHOST,8,8, &
          & DOFS_MAPPING%GHOST_LIST,'("    Ghost dofs   :",8(X,I7))','(18X,8(X,I7))',ERR,ERROR,*999)
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
    IF(ASSOCIATED(DOMAIN%MAPPINGS%NODES)) CALL DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*998)
998 IF(ASSOCIATED(DOMAIN%MAPPINGS%DOFS)) CALL DOMAIN_MAPPINGS_DOFS_FINALISE(DOMAIN%MAPPINGS,ERR,ERROR,*997)
997 CALL ERRORS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_NODES_DOFS_CALCULATE
  
  !
  !================================================================================================================================
  !

  SUBROUTINE DOMAIN_MAPPINGS_NODES_FINALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_NODES_FINALISE
    !###  Description:
    !###    Finalises the node mapping in the given domain mappings.

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPINGS%NODES,ERR,ERROR,*999)
        DEALLOCATE(DOMAIN_MAPPINGS%NODES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_MAPPINGS_NODES_INITIALISE(DOMAIN_MAPPINGS,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_MAPPINGS_NODES_INITIALISE
    !###  Description:
    !###    Intialises the node mapping in the given domain mapping.

    !Argument variables
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: DOMAIN_MAPPINGS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    
    CALL ENTERS("DOMAIN_MAPPINGS_NODES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPINGS)) THEN
      IF(ASSOCIATED(DOMAIN_MAPPINGS%NODES)) THEN
        CALL FLAG_ERROR("Domain nodes mappings are already associated",ERR,ERROR,*999)
      ELSE
        ALLOCATE(DOMAIN_MAPPINGS%NODES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain mappings nodes",ERR,ERROR,*999)
        CALL DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPINGS%NODES,DOMAIN_MAPPINGS%DOMAIN%DECOMPOSITION%NUMBER_OF_DOMAINS, &
          & ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_CALCULATE
    !###  Description:
    !###    Calculates the domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
          CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
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
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_INITIALISE_FROM_MESH
    !###  Description:
    !###    Initialises the local domain topology from the mesh topology.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
                DOMAIN_NODES%NODES(local_node)%USER_NUMBER=MESH_NODES%NODES(global_node)%USER_NUMBER
                DOMAIN_NODES%NODES(local_node)%LOCAL_NUMBER=local_node
                DOMAIN_NODES%NODES(local_node)%GLOBAL_NUMBER=global_node
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_SURROUNDING_ELEMENTS=0
                NULLIFY(DOMAIN_NODES%NODES(local_node)%SURROUNDING_ELEMENTS)
                DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES=MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX( &
                  & MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node partial derivative index",ERR,ERROR,*999)
                DOMAIN_NODES%NODES(local_node)%PARTIAL_DERIVATIVE_INDEX=MESH_NODES%NODES(global_node)%PARTIAL_DERIVATIVE_INDEX
                ALLOCATE(DOMAIN_NODES%NODES(local_node)%DOF_INDEX(MESH_NODES%NODES(global_node)%NUMBER_OF_DERIVATIVES),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node dof index",ERR,ERROR,*999)
                DO nk=1,DOMAIN_NODES%NODES(local_node)%NUMBER_OF_DERIVATIVES
                  ny=ny+1
                  DOMAIN_NODES%NODES(local_node)%DOF_INDEX(nk)=ny
                  DOMAIN_DOFS%DOF_INDEX(1,ny)=nk
                  DOMAIN_DOFS%DOF_INDEX(2,ny)=local_node
                ENDDO !nk
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
                  global_node=MESH_ELEMENTS%ELEMENTS(global_element)%GLOBAL_ELEMENT_NODES(nn)
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
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Node user number = ",DOMAIN_NODES%NODES(np)%USER_NUMBER, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of derivatives = ", &
          & DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,4,4, &
          & DOMAIN_NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX,'("      Partial derivative index(nk) :",4(X,I9))','(36X,4(X,I9))', &
          & ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_NODES%NODES(np)%NUMBER_OF_DERIVATIVES,4,4, &
          & DOMAIN_NODES%NODES(np)%DOF_INDEX,'("      Degree-of-freedom index(nk)  :",4(X,I9))','(36X,4(X,I9))', &
          & ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_DOFS_FINALISE
    !###  Description:
    !###    Finalises the dofs in the given domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_DOFS_INITIALISE
    !###  Description:
    !###    Initialises the dofs data structures for a domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_ELEMENT_FINALISE
    !###  Description:
    !###    Finalises the given domain topology element.

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_ELEMENT_INITIALISE
    !###  Description:
    !###    Initialises the given domain topology element.

    !Argument variables
    TYPE(DOMAIN_ELEMENT_TYPE) :: ELEMENT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_ELEMENTS_FINALISE
    !###  Description:
    !###    Finalises the elements in the given domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE
    !###  Description:
    !###    Initialises the element data structures for a domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_FINALISE(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_FINALISE
    !###  Description:
    !###    Finalises the topology in the given domain.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN)) THEN
      CALL DOMAIN_TOPOLOGY_NODES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_DOFS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_ELEMENTS_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
      CALL DOMAIN_TOPOLOGY_LINES_FINALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_TOPOLOGY_INITIALISE(DOMAIN,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_INITIALISE
    !###  Description:
    !###    Initialises the topology for a given domain.

    !Argument variables
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
        !Initialise the topology components
        CALL DOMAIN_TOPOLOGY_ELEMENTS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_NODES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_DOFS_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
        CALL DOMAIN_TOPOLOGY_LINES_INITIALISE(DOMAIN%TOPOLOGY,ERR,ERROR,*999)
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

  SUBROUTINE DOMAIN_TOPOLOGY_LINE_FINALISE(LINE,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_LINE_FINALISE
    !###  Description:
    !###    Finalises a line in the given domain topology and deallocates all memory.

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_LINE_INITIALISE(LINE,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_LINE_INITIALISE
    !###  Description:
    !###    Initialises the line data structure for a domain topology line.

    !Argument variables
    TYPE(DOMAIN_LINE_TYPE) :: LINE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_LINES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_LINES_FINALISE
    !###  Description:
    !###    Finalises the lines in the given domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_LINES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_LINES_INITIALISE
    !###  Description:
    !###    Initialises the line data structures for a domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_NODE_FINALISE
    !###  Description:
    !###    Finalises the given domain topology node and deallocates all memory.

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(NODE%PARTIAL_DERIVATIVE_INDEX)) DEALLOCATE(NODE%PARTIAL_DERIVATIVE_INDEX)
    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
    IF(ALLOCATED(NODE%NODE_Lines)) DEALLOCATE(NODE%NODE_LINES)
 
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_FINALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_NODE_INITIALISE
    !###  Description:
    !###    Initialises the given domain topology node.

    !Argument variables
    TYPE(DOMAIN_NODE_TYPE) :: NODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%USER_NUMBER=0
    NODE%LOCAL_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NODE%NUMBER_OF_NODE_LINES=0
    
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_TOPOLOGY_NODE_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE DOMAIN_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_NODES_FINALISE
    !###  Description:
    !###    Finalises the nodees in the given domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_NODES_INITIALISE
    !###  Description:
    !###    Initialises the nodes data structures for a domain topology.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: DOMAIN_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE
    !###  Description:
    !###    Calculates the element numbers surrounding a node for a domain.

    !Argument variables
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_CREATE_FINISH
    !###  Description:
    !###    Finishes the process of creating a mesh on a region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: component_idx,mesh_idx
    LOGICAL :: FINISHED
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(MESH)) THEN
        IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
          IF(ASSOCIATED(MESH%REGION)) THEN
            IF(MESH%REGION%USER_NUMBER==REGION%USER_NUMBER) THEN            
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
!!TODO: Really create an all decomposition????
              !Create the default domain decomposition
              !CALL DECOMPOSITION_CREATE_START(0,MESH,DECOMPOSITION,ERR,ERROR,*999)
              !CALL DECOMPOSITION_CREATE_FINISH(MESH,DECOMPOSITION,ERR,ERROR,*999)            
            ELSE
              LOCAL_ERROR="The region user number of the specified mesh ("// &
                & TRIM(NUMBER_TO_VSTRING(MESH%REGION%USER_NUMBER,"*",ERR,ERROR))// &
                & ") does not match the user number of the specified region ("// &
                & TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))//")"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh region is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Region = ",REGION%USER_NUMBER,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of meshes = ",REGION%MESHES%NUMBER_OF_MESHES,ERR,ERROR,*999)
      DO mesh_idx=1,REGION%MESHES%NUMBER_OF_MESHES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh number = ",mesh_idx,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
          & REGION%MESHES%MESHES(mesh_idx)%PTR%USER_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
          & REGION%MESHES%MESHES(mesh_idx)%PTR%GLOBAL_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of dimensions = ", &
          & REGION%MESHES%MESHES(mesh_idx)%PTR%NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
      ENDDO !problem_idx    
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

  SUBROUTINE MESH_CREATE_REGULAR(USER_NUMBER,REGION,ORIGIN,MAXIMUM_EXTENT,NUMBER_ELEMENTS_XI,BASIS,MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_CREATE_REGULAR
    !###  Description:
    !###    Creates the regular mesh with the given USER_NUMBER in the specifed REGION. The mesh starts at the ORIGIN(:) and has
    !###    a maximum extent position of MAXIMUM_EXTENT(:) with the NUMBER_OF_ELEMENTS(:) in each direction. Each element is of
    !###    the specified BASIS type. A pointer to the finished mesh is returned in MESH.  

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    REAL(DP), INTENT(IN) :: ORIGIN(:),MAXIMUM_EXTENT(:)
    INTEGER(INTG), INTENT(IN) :: NUMBER_ELEMENTS_XI(:)
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: ni,ne,ne1,ne2,ne3,NN,nn1,nn2,nn3,np,np1,np2,np3,TOTAL_NUMBER_OF_NODES_XI(3),TOTAL_NUMBER_ELEMENTS_XI(3), &
      & TOTAL_NUMBER_OF_NODES,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_DIMENSIONS
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_NODES(:)
    REAL(DP) :: INITIAL_POSITION(3),DELTA_COORD(3),MY_ORIGIN(3),MY_EXTENT(3),MESH_SIZE(3)
    LOGICAL :: BASIS_OK,ELEMENTS_OK
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: MESH_ELEMENTS
    TYPE(NODES_TYPE), POINTER :: NODES

    CALL ENTERS("MESH_CREATE_REGULAR",ERR,ERROR,*999)

    NULLIFY(MESH)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
        !Determine the coordinate system and create the regular mesh for that system
        SELECT CASE(REGION%COORDINATE_SYSTEM%TYPE)
        CASE(COORDINATE_RECTANGULAR_CARTESIAN_TYPE)
          NUMBER_OF_DIMENSIONS=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
          IF(SIZE(ORIGIN)==NUMBER_OF_DIMENSIONS) THEN
            IF(SIZE(MAXIMUM_EXTENT)==NUMBER_OF_DIMENSIONS) THEN
              IF(ASSOCIATED(BASIS)) THEN
                SELECT CASE(BASIS%TYPE)
                CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                  IF(BASIS%NUMBER_OF_XI==SIZE(NUMBER_ELEMENTS_XI,1)) THEN
                    ELEMENTS_OK=ALL(NUMBER_ELEMENTS_XI>0)
                    IF(.NOT.ELEMENTS_OK) CALL FLAG_ERROR("Must have 1 or more elements in all directions",ERR,ERROR,*999)
                    BASIS_OK=ALL(BASIS%COLLAPSED_XI==BASIS_NOT_COLLAPSED)
                    IF(.NOT.BASIS_OK) CALL FLAG_ERROR("Degenerate (collapsed) basis not implemented",ERR,ERROR,*999)
                    !Calculate sizes
                    TOTAL_NUMBER_OF_NODES=1
                    TOTAL_NUMBER_OF_ELEMENTS=1
                    TOTAL_NUMBER_OF_NODES_XI=1
                    TOTAL_NUMBER_ELEMENTS_XI=0
                    DELTA_COORD=0.0_DP
                    DO ni=1,BASIS%NUMBER_OF_XI
                      TOTAL_NUMBER_OF_NODES_XI(ni)=(BASIS%NUMBER_OF_NODES_XI(ni)-2)*NUMBER_ELEMENTS_XI(ni)+ &
                        & NUMBER_ELEMENTS_XI(ni)+1
                      TOTAL_NUMBER_ELEMENTS_XI(ni)=NUMBER_ELEMENTS_XI(ni)
                      TOTAL_NUMBER_OF_NODES=TOTAL_NUMBER_OF_NODES*TOTAL_NUMBER_OF_NODES_XI(ni)
                      TOTAL_NUMBER_OF_ELEMENTS=TOTAL_NUMBER_OF_ELEMENTS*TOTAL_NUMBER_ELEMENTS_XI(ni)
                    ENDDO !ni
                    !Create the default node set
                    CALL NODES_CREATE_START(TOTAL_NUMBER_OF_NODES,REGION,NODES,ERR,ERROR,*999)
                    MY_ORIGIN=0.0_DP
                    MY_EXTENT=0.0_DP
                    MY_ORIGIN(1:NUMBER_OF_DIMENSIONS)=ORIGIN
                    MY_EXTENT(1:NUMBER_OF_DIMENSIONS)=MAXIMUM_EXTENT
                    MESH_SIZE=MY_EXTENT-MY_ORIGIN
                    DO ni=1,BASIS%NUMBER_OF_XI
                      !This assumes that the xi directions are aligned with the coordinate directions
                      DELTA_COORD(ni)=MESH_SIZE(ni)/REAL(NUMBER_ELEMENTS_XI(ni),DP)
                    ENDDO !ni
                    DO np3=1,TOTAL_NUMBER_OF_NODES_XI(3)
                      DO np2=1,TOTAL_NUMBER_OF_NODES_XI(2)
                        DO np1=1,TOTAL_NUMBER_OF_NODES_XI(1)
                          np=np1+(np2-1)*TOTAL_NUMBER_OF_NODES_XI(1)+(np3-1)*TOTAL_NUMBER_OF_NODES_XI(2)
                          INITIAL_POSITION(1)=MY_ORIGIN(1)+REAL(np1-1,DP)*DELTA_COORD(1)
                          INITIAL_POSITION(2)=MY_ORIGIN(2)+REAL(np2-1,DP)*DELTA_COORD(2)
                          INITIAL_POSITION(3)=MY_ORIGIN(3)+REAL(np3-1,DP)*DELTA_COORD(3)
                          CALL NODE_INITIAL_POSITION_SET(np,INITIAL_POSITION(1:NUMBER_OF_DIMENSIONS),NODES,ERR,ERROR,*999)
                        ENDDO !np1
                      ENDDO !np2
                    ENDDO !np3
                    CALL NODES_CREATE_FINISH(REGION,ERR,ERROR,*999)
                    !Create the mesh
                    CALL MESH_CREATE_START(USER_NUMBER,REGION,SIZE(NUMBER_ELEMENTS_XI,1),MESH,ERR,ERROR,*999)
                    !Create the elements
                    CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,TOTAL_NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
                    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,1,BASIS,MESH_ELEMENTS,ERR,ERROR,*999)
                    !Set the elements for the regular mesh
                    ALLOCATE(ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
                    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate element nodes",ERR,ERROR,*999)                
                    !Step in the xi(3) direction
                    DO ne3=1,TOTAL_NUMBER_ELEMENTS_XI(3)+1
                      DO ne2=1,TOTAL_NUMBER_ELEMENTS_XI(2)+1
                        DO ne1=1,TOTAL_NUMBER_ELEMENTS_XI(1)+1
                          IF(BASIS%NUMBER_OF_XI<3.OR.ne3<=TOTAL_NUMBER_ELEMENTS_XI(3)) THEN
                            IF(BASIS%NUMBER_OF_XI<2.OR.ne2<=TOTAL_NUMBER_ELEMENTS_XI(2)) THEN
                              IF(ne1<=TOTAL_NUMBER_ELEMENTS_XI(1)) THEN
                                ne=ne1
                                np=1+(ne1-1)*(BASIS%NUMBER_OF_NODES_XI(1)-1)
                                IF(BASIS%NUMBER_OF_XI>1) THEN
                                  ne=ne+(ne2-1)*TOTAL_NUMBER_ELEMENTS_XI(1)
                                  np=np+(ne2-1)*TOTAL_NUMBER_OF_NODES_XI(1)*(BASIS%NUMBER_OF_NODES_XI(2)-1)
                                  IF(BASIS%NUMBER_OF_XI>2) THEN
                                    ne=ne+(ne3-1)*TOTAL_NUMBER_ELEMENTS_XI(1)*TOTAl_NUMBER_ELEMENTS_XI(2)
                                    np=np+(ne3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)* &
                                      & (BASIS%NUMBER_OF_NODES_XI(3)-1)
                                  ENDIF
                                ENDIF
                                nn=0
                                DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                  nn=nn+1
                                  ELEMENT_NODES(nn)=np+(nn1-1)                              
                                ENDDO !nn1
                                IF(BASIS%NUMBER_OF_XI>1) THEN
                                  DO nn2=2,BASIS%NUMBER_OF_NODES_XI(2)
                                    DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                      nn=nn+1
                                      ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)
                                    ENDDO !nn1
                                  ENDDO !nn2
                                  IF(BASIS%NUMBER_OF_XI>2) THEN
                                    DO nn3=2,BASIS%NUMBER_OF_NODES_XI(3)
                                      DO nn2=1,BASIS%NUMBER_OF_NODES_XI(2)
                                        DO nn1=1,BASIS%NUMBER_OF_NODES_XI(1)
                                          nn=nn+1
                                          ELEMENT_NODES(nn)=np+(nn1-1)+(nn2-1)*TOTAL_NUMBER_OF_NODES_XI(1)+ &
                                            & (nn3-1)*TOTAL_NUMBER_OF_NODES_XI(1)*TOTAL_NUMBER_OF_NODES_XI(2)
                                        ENDDO !nn1
                                      ENDDO !nn2
                                    ENDDO !nn3
                                  ENDIF
                                ENDIF
                                CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,MESH_ELEMENTS,ELEMENT_NODES, &
                                  & ERR,ERROR,*999)
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDDO !ne1
                      ENDDO !ne2
                    ENDDO !ne3
                    CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,1,ERR,ERROR,*999)
                    !Finish the mesh
                    CALL MESH_CREATE_FINISH(REGION,MESH,ERR,ERROR,*999)
                  ELSE
                    CALL FLAG_ERROR("The number of xi directions of the given basis does not match the size of &
                      &the number of elements for the mesh",ERR,ERROR,*999)
                  ENDIF
                CASE(BASIS_SIMPLEX_TYPE)                  
                  CALL FLAG_ERROR("Regular meshes with simplex basis types is not implemented",ERR,ERROR,*999)
                CASE DEFAULT
                  CALL FLAG_ERROR("Basis type is either invalid or not implemented",ERR,ERROR,*999)
                END SELECT
              ELSE
                CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("Size of the mesh extent does not match the region number of dimensions",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Size of origin does not match the region number of dimensions",ERR,ERROR,*999)
          ENDIF
        CASE DEFAULT
          CALL FLAG_ERROR("Non rectangular cartesian coordinate systems are not implemented",ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FLAG_ERROR("Region coordinate system is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_CREATE_REGULAR")
    RETURN
999 IF(ALLOCATED(ELEMENT_NODES)) DEALLOCATE(ELEMENT_NODES)
    CALL ERRORS("MESH_CREATE_REGULAR",ERR,ERROR)
    CALL EXITS("MESH_CREATE_REGULAR")
    RETURN 1   
  END SUBROUTINE MESH_CREATE_REGULAR
  
  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_CREATE_START(USER_NUMBER,REGION,NUMBER_OF_DIMENSIONS,MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_CREATE_START
    !###  Description:
    !###    Starts the process of creating a mesh defined by a user number with the specified NUMBER_OF_DIMENSIONS in the region
    !###    identified by REGION.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: mesh_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER :: NEW_MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESH)
    NULLIFY(NEW_MESHES)

    CALL ENTERS("",ERR,ERROR,*999)

    NULLIFY(MESH)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%COORDINATE_SYSTEM)) THEN
        IF(ASSOCIATED(REGION%MESHES)) THEN
          CALL MESH_USER_NUMBER_FIND(USER_NUMBER,REGION,MESH,ERR,ERROR,*999)
          IF(ASSOCIATED(MESH)) THEN
            LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(USER_NUMBER,"*",ERR,ERROR))// &
              & " has already been created on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
            CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(NUMBER_OF_DIMENSIONS>0.AND.NUMBER_OF_DIMENSIONS<=REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS) THEN
              ALLOCATE(NEW_MESH,STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new mesh",ERR,ERROR,*999)
              CALL MESH_INITIALISE(NEW_MESH,ERR,ERROR,*999)
              !Set default mesh values
              NEW_MESH%USER_NUMBER=USER_NUMBER
              NEW_MESH%GLOBAL_NUMBER=REGION%MESHES%NUMBER_OF_MESHES+1
              NEW_MESH%MESHES=>REGION%MESHES
              NEW_MESH%REGION=>REGION
              NEW_MESH%NUMBER_OF_DIMENSIONS=NUMBER_OF_DIMENSIONS
              NEW_MESH%NUMBER_OF_COMPONENTS=1
              !Initialise mesh topology and decompositions
              CALL MESH_TOPOLOGY_INITIALISE(NEW_MESH,ERR,ERROR,*999)
              CALL DECOMPOSITIONS_INITIALISE(NEW_MESH,ERR,ERROR,*999)
              !Add new mesh into list of meshes in the region
              ALLOCATE(NEW_MESHES(REGION%MESHES%NUMBER_OF_MESHES+1),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new meshes",ERR,ERROR,*999)
              DO mesh_idx=1,REGION%MESHES%NUMBER_OF_MESHES
                NEW_MESHES(mesh_idx)%PTR=>REGION%MESHES%MESHES(mesh_idx)%PTR
              ENDDO !mesh_idx
              NEW_MESHES(REGION%MESHES%NUMBER_OF_MESHES+1)%PTR=>NEW_MESH
              IF(ASSOCIATED(REGION%MESHES%MESHES)) DEALLOCATE(REGION%MESHES%MESHES)
              REGION%MESHES%MESHES=>NEW_MESHES
              REGION%MESHES%NUMBER_OF_MESHES=REGION%MESHES%NUMBER_OF_MESHES+1
              MESH=>NEW_MESH
            ELSE
              IF(NUMBER_OF_DIMENSIONS>0) THEN
                LOCAL_ERROR="Number of mesh dimensions ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
                  & ") must be <= number of region dimensions ("// &
                  & TRIM(NUMBER_TO_VSTRING(REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))//")"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ELSE
                CALL FLAG_ERROR("Number of mesh dimensions must be > 0",ERR,ERROR,*999)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
            & " are not associated. Initialise meshes first."
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The coordinate system on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF

    CALL EXITS("MESH_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_MESH)) DEALLOCATE(NEW_MESH)
    IF(ASSOCIATED(NEW_MESHES)) DEALLOCATE(NEW_MESHES)
    NULLIFY(MESH)
    CALL ERRORS("MESH_CREATE_START",ERR,ERROR)    
    CALL EXITS("MESH_CREATE_START")
    RETURN 1
   
  END SUBROUTINE MESH_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_DESTROY(USER_NUMBER,REGION,ERR,ERROR,*)

    !#### Subroutine: MESH_DESTROY
    !###  Description:
    !###    Destroys the mesh identified by a user number on the given region and deallocates all memory.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: mesh_idx,mesh_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(MESH_PTR_TYPE), POINTER :: NEW_MESHES(:)

    NULLIFY(NEW_MESHES)

    CALL ENTERS("MESH_DESTROY",ERR,ERROR,*999)

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

  SUBROUTINE MESH_FINALISE(MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_FINALISE
    !###  Description:
    !###    Finalises a mesh and deallocates all memory

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("MESH_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      CALL MESH_TOPOLOGY_FINALISE(MESH,ERR,ERROR,*999)
      CALL DECOMPOSITIONS_FINALISE(MESH,ERR,ERROR,*999)
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

  SUBROUTINE MESH_INITIALISE(MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_INITIALISE
    !###  Description:
    !###    Initialises a mesh.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("MESH_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
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
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
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
  
  !#### Generic-Subroutine: MESH_NUMBER_OF_COMPONENTS_SET
  !###  Description:
  !###    Sets/changes the number of mesh components for a mesh.
  !###  Child-subroutines: MESH_NUMBER_OF_COMPONENTS_SET_NUMBER,MESH_NUMBER_OF_COMPONENTS_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET_NUMBER(USER_NUMBER,REGION,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !#### Subroutine: MESH_NUMBER_OF_COMPONENTS_SET_NUMBER
    !###  Description:
    !###    Changes/sets the number of mesh components for a mesh identified by a given user number on a region.
    !###  Parent-subroutine: MESH_NUMBER_OF_COMPONENTS_SET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH

    CALL ENTERS("MESH_NUMBER_OF_COMPONENTS_SET_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN

!!TODO: Take in region number here and user FIND_REGION_NUMBER. This would require FIND_REGION_NUMBER to be moved from
!!REGION_ROUTINES otherwise there will be a circular module reference.

      CALL MESH_USER_NUMBER_FIND(USER_NUMBER,REGION,MESH,ERR,ERROR,*999)
      CALL MESH_NUMBER_OF_COMPONENTS_SET_PTR(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
      
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_SET_NUMBER")
    RETURN
999 CALL ERRORS("MESH_NUMBER_OF_COMPONENTS_SET_NUMBER",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_SET_NUMBER")
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET_PTR(MESH,NUMBER_OF_COMPONENTS,ERR,ERROR,*)

    !#### Subroutine: MESH_NUMBER_OF_COMPONENTS_SET_PTR
    !###  Description:
    !###    Changes/sets the number of mesh components for a mesh identified by a pointer.
    !###  Parent-subroutine: MESH_NUMBER_OF_COMPONENTS_SET

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(MESH_TOPOLOGY_PTR_TYPE), POINTER :: NEW_TOPOLOGY(:)

    NULLIFY(NEW_TOPOLOGY)
    
    CALL ENTERS("MESH_NUMBER_OF_COMPONENTS_SET_PTR",ERR,ERROR,*999)

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
              DO component_idx=MESH%NUMBER_OF_COMPONENTS+1,NUMBER_OF_COMPONENTS
                ALLOCATE(NEW_TOPOLOGY(component_idx)%PTR,STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate new topology component",ERR,ERROR,*999)
                NULLIFY(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)
                NULLIFY(MESH%TOPOLOGY(component_idx)%PTR%NODES)
                !Initialise the topology components
                CALL MESH_TOPOLOGY_ELEMENTS_INITIALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
                CALL MESH_TOPOLOGY_NODES_INITIALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
              ENDDO !component_idx
            ENDIF
            DEALLOCATE(MESH%TOPOLOGY)
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
    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_SET_PTR")
    RETURN
!!TODO: tidy up memory deallocation on error
999 CALL ERRORS("MESH_NUMBER_OF_COMPONENTS_SET_PTR",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_COMPONENTS_SET_PTR")
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_COMPONENTS_SET_PTR

  !
  !================================================================================================================================
  !
  
  !#### Generic-Subroutine: MESH_NUMBER_OF_ELEMENTS_SET
  !###  Description:
  !###    Sets/changes the number of elements for a mesh.
  !###  Child-subroutines: MESH_NUMBER_OF_ELEMENTS_SET_NUMBER,MESH_NUMBER_OF_ELEMENTS_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET_NUMBER(USER_NUMBER,REGION,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !#### Subroutine: MESH_NUMBER_OF_ELEMENTS_SET_NUMBER
    !###  Description:
    !###    Changes/sets the number of elements for a mesh identified by a given user number on a region.
    !###  Parent-subroutine: MESH_NUMBER_OF_ELEMENTS_SET

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(MESH_TYPE), POINTER :: MESH

    CALL ENTERS("MESH_NUMBER_OF_ELEMENTS_SET_NUMBER",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      
!!TODO: Take in region number here and user FIND_REGION_NUMBER. This would require FIND_REGION_NUMBER to be moved from
!!REGION_ROUTINES otherwise there will be a circular module reference.

      CALL MESH_USER_NUMBER_FIND(USER_NUMBER,REGION,MESH,ERR,ERROR,*999)
      CALL MESH_NUMBER_OF_ELEMENTS_SET_PTR(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
      
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_SET_NUMBER")
    RETURN
999 CALL ERRORS("MESH_NUMBER_OF_ELEMENTS_SET_NUMBER",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_SET_NUMBER")
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET_NUMBER

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET_PTR(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*)

    !#### Subroutine: MESH_NUMBER_OF_ELEMENTS_SET_PTR
    !###  Description:
    !###    Changes/sets the number of elements for a mesh identified by a pointer.
    !###  Parent-subroutine: MESH_NUMBER_OF_ELEMENTS_SET

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: component_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_NUMBER_OF_ELEMENTS_SET_PTR",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(NUMBER_OF_ELEMENTS>0) THEN
        IF(MESH%MESH_FINISHED) THEN
          CALL FLAG_ERROR("Mesh has been finished",ERR,ERROR,*999)
        ELSE
          IF(NUMBER_OF_ELEMENTS/=MESH%NUMBER_OF_ELEMENTS) THEN
            IF(ASSOCIATED(MESH%TOPOLOGY)) THEN
              DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
                IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR)) THEN
                  IF(ASSOCIATED(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS)) THEN
                    IF(MESH%TOPOLOGY(component_idx)%PTR%ELEMENTS%NUMBER_OF_ELEMENTS>0) THEN
!!TODO: Reallocate the elements and copy information. 
                      CALL FLAG_ERROR("Not implemented",ERR,ERROR,*999)
                    ENDIF
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("Mesh topology component pointer is not associated",ERR,ERROR,*999)
                ENDIF
              ENDDO !component_idx
            ELSE
              CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
            ENDIF
            MESH%NUMBER_OF_ELEMENTS=NUMBER_OF_ELEMENTS
          ENDIF
        ENDIF
      ELSE
        LOCAL_ERROR="The specified number of elements ("//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ELEMENTS,"*",ERR,ERROR))// &
          & ") is illegal. You must have >0 elements"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF      
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_SET_PTR")
    RETURN
999 CALL ERRORS("MESH_NUMBER_OF_ELEMENTS_SET_PTR",ERR,ERROR)    
    CALL EXITS("MESH_NUMBER_OF_ELEMENTS_SET_PTR")
    RETURN 1
   
  END SUBROUTINE MESH_NUMBER_OF_ELEMENTS_SET_PTR

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_CALCULATE
    !###  Description:
    !###    Calculates the mesh topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_DOFS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_DOFS_CALCULATE
    !###  Description:
    !###    Calculates the degrees-of-freedom for a mesh topology. 

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_DOFS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_DOFS_FINALISE
    !###  Description:
    !###    Finalises the dof data structures for a mesh topology and deallocates any memory. 

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_DOFS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_DOFS_INITIALISE
    !###  Description:
    !###    Initialises the dofs in a given mesh topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_BASIS_SET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_BASIS_SET
    !###  Description:
    !###    Changes/sets the basis for a global element identified by a given global number.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_BASIS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>0.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(ASSOCIATED(BASIS)) THEN
            IF(ALLOCATED(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES))  &
              & DEALLOCATE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES)
            ALLOCATE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate user element nodes",ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES=1
            IF(ALLOCATED(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES))  &
              & DEALLOCATE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES)
            ALLOCATE(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES(BASIS%NUMBER_OF_NODES),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global element nodes",ERR,ERROR,*999)
            ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES=1
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
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_BASIS_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_BASIS_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_BASIS_SET")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_BASIS_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(MESH,MESH_COMPONENT_NUMBER,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH
    !###  Description:
    !###    Finishes the process of creating elements for a specified mesh component in a mesh topology.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
            MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS_FINISHED=.TRUE.
          ELSE
            CALL FLAG_ERROR("Mesh elements is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        LOCAL_ERROR="The speficied mesh component number of "//TRIM(NUMBER_TO_VSTRING(MESH_COMPONENT_NUMBER,"*",ERR,ERROR))// &
          & "is invalid. The component number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(MESH%NUMBER_OF_COMPONENTS,"*",ERR,ERROR))
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of global elements = ",MESH%NUMBER_OF_ELEMENTS, &
        & ERR,ERROR,*999)
      DO ne=1,MESH%NUMBER_OF_ELEMENTS
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Element = ",ne,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global number        = ", &
          & MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    User number          = ", &
          & MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%USER_NUMBER,ERR,ERROR,*999)
        IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS)) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Basis number         = ", &
            & MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%BASIS%USER_NUMBER,ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Basis is not associated",ERR,ERROR,*999)
        ENDIF
        IF(ALLOCATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES)) THEN
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)% &
            & BASIS%NUMBER_OF_NODES,8,8,MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES, &
            & '("    User element nodes   =",8(X,I6))','(26X,8(X,I6))',ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("User element nodes are not associated",ERR,ERROR,*999)
        ENDIF
        IF(ALLOCATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES)) THEN
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)% &
            & BASIS%NUMBER_OF_NODES,8,8,MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES, &
            & '("    Global element nodes =",8(X,I6))','(26X,8(X,I6))',ERR,ERROR,*999)
        ELSE
          CALL FLAG_ERROR("Global element nodes are not associated",ERR,ERROR,*999)
        ENDIF
      ENDDO !ne
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMPONENT_NUMBER,BASIS,ELEMENTS,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_CREATE_START
    !###  Description:
    !###    Starts the process of creating elements in the mesh component identified by MESH and component_idx. The elements will
    !###    will be created with a default basis of BASIS. ELEMENTS is the returned pointer to the MESH_ELEMENTS data structure.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(IN) :: MESH_COMPONENT_NUMBER
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: ne
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR,*999)

    NULLIFY(ELEMENTS)
    IF(ASSOCIATED(MESH)) THEN     
      IF(MESH_COMPONENT_NUMBER>0.AND.MESH_COMPONENT_NUMBER<=MESH%NUMBER_OF_COMPONENTS) THEN
        IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR)) THEN
          IF(ASSOCIATED(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS)) THEN
            ELEMENTS=>MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR%ELEMENTS
            IF(ASSOCIATED(ELEMENTS%ELEMENTS)) THEN
              CALL FLAG_ERROR("Mesh topology already has elements associated",ERR,ERROR,*998)
            ELSE
              IF(ASSOCIATED(BASIS)) THEN
                ALLOCATE(ELEMENTS%ELEMENTS(MESH%NUMBER_OF_ELEMENTS),STAT=ERR)
                IF(ERR/=0) CALL FLAG_ERROR("Could not allocate individual elements",ERR,ERROR,*999)
                ELEMENTS%NUMBER_OF_ELEMENTS=MESH%NUMBER_OF_ELEMENTS !Psuedo inheritance of the number of elements
                ELEMENTS%ELEMENTS_FINISHED=.FALSE.
                !Set up the default values and allocate element structures
                DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                  CALL MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
                  ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER=ne
                  ELEMENTS%ELEMENTS(ne)%USER_NUMBER=ne
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
999 CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(MESH%TOPOLOGY(MESH_COMPONENT_NUMBER)%PTR,ERR,ERROR,*998)
998 NULLIFY(ELEMENTS)
 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_CREATE_START",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_CREATE_START")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_CREATE_START

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_DESTROY(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_DESTROY
    !###  Description:
    !###    Destroys the elements in a mesh topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: ne

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        IF(TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS>0) THEN
          IF(ASSOCIATED(TOPOLOGY%ELEMENTS%ELEMENTS)) THEN
            DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
              DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%USER_ELEMENT_NODES)
              DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES)
            ENDDO !ne
            DEALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS)
            DEALLOCATE(TOPOLOGY%ELEMENTS)
          ELSE
            CALL FLAG_ERROR("Mesh topology elements elements is not associated",ERR,ERROR,*999)
          ENDIF
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*999)
      ENDIF
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENT_FINALISE(ELEMENT,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENT_FINALISE
    !###  Description:
    !###    Finalises the given mesh topology element.

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_ELEMENT_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ELEMENT%USER_ELEMENT_NODES)) DEALLOCATE(ELEMENT%USER_ELEMENT_NODES)
    IF(ALLOCATED(ELEMENT%GLOBAL_ELEMENT_NODES)) DEALLOCATE(ELEMENT%GLOBAL_ELEMENT_NODES)
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE(ELEMENT,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENT_INITIALISE
    !###  Description:
    !###    Initialises the given mesh topology element.

    !Argument variables
    TYPE(MESH_ELEMENT_TYPE) :: ELEMENT
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR,*999)

    ELEMENT%USER_NUMBER=0
    ELEMENT%GLOBAL_NUMBER=0
    NULLIFY(ELEMENT%BASIS)
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_ELEMENT_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET(GLOBAL_NUMBER,ELEMENTS,BASIS,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_ELEMENT_BASIS_SET
    !###  Description:
    !###    Changes/sets the basis for a mesh element identified by a given global number.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(BASIS_TYPE), POINTER :: BASIS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(GLOBAL_NUMBER,ELEMENTS,USER_ELEMENT_NODES,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET
    !###  Description:
    !###    Changes/sets the element nodes for a mesh element identified by a given global number.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    INTEGER(INTG), INTENT(IN) :: USER_ELEMENT_NODES(:)
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: nn,NUMBER_OF_BAD_NODES,GLOBAL_NODE_NUMBER
    INTEGER(INTG), ALLOCATABLE :: GLOBAL_ELEMENT_NODES(:),BAD_NODES(:)
    LOGICAL :: ELEMENT_NODES_OK,NODE_EXISTS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(NODES_TYPE), POINTER :: NODES
    TYPE(REGION_TYPE), POINTER :: REGION

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          IF(SIZE(USER_ELEMENT_NODES,1)==ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES) THEN
            IF(ASSOCIATED(ELEMENTS%MESH)) THEN
              IF(ASSOCIATED(ELEMENTS%MESH%REGION)) THEN
                REGION=>ELEMENTS%MESH%REGION
                IF(ASSOCIATED(REGION%NODES)) THEN
                  NODES=>REGION%NODES
                  ELEMENT_NODES_OK=.TRUE.
                  ALLOCATE(GLOBAL_ELEMENT_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate global element nodes",ERR,ERROR,*999)
                  ALLOCATE(BAD_NODES(ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate bad nodes",ERR,ERROR,*999)
                  NUMBER_OF_BAD_NODES=0
                  DO nn=1,ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%BASIS%NUMBER_OF_NODES
                    CALL NODE_CHECK_EXISTS(USER_ELEMENT_NODES(nn),REGION,NODE_EXISTS,GLOBAL_NODE_NUMBER,ERR,ERROR,*999)
                    IF(.NOT.NODE_EXISTS) THEN
                      NUMBER_OF_BAD_NODES=NUMBER_OF_BAD_NODES+1
                      BAD_NODES(NUMBER_OF_BAD_NODES)=USER_ELEMENT_NODES(nn)
                      ELEMENT_NODES_OK=.FALSE.
                    ELSE
                      GLOBAL_ELEMENT_NODES(nn)=GLOBAL_NODE_NUMBER
                    ENDIF
                  ENDDO !nn
                  IF(ELEMENT_NODES_OK) THEN
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_ELEMENT_NODES=USER_ELEMENT_NODES
                    ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%GLOBAL_ELEMENT_NODES=GLOBAL_ELEMENT_NODES
                  ELSE
                    IF(NUMBER_OF_BAD_NODES==1) THEN
                      LOCAL_ERROR="The element node "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))// &
                        & " is not defined on region "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
                    ELSE
                      LOCAL_ERROR="The element nodes "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(1),"*",ERR,ERROR))
                      DO nn=2,NUMBER_OF_BAD_NODES-1
                        LOCAL_ERROR=LOCAL_ERROR//","//TRIM(NUMBER_TO_VSTRING(BAD_NODES(nn),"*",ERR,ERROR))
                      ENDDO !nn
                      LOCAL_ERROR=LOCAL_ERROR//" & "//TRIM(NUMBER_TO_VSTRING(BAD_NODES(NUMBER_OF_BAD_NODES),"*",ERR,ERROR))// &
                        & " are not defined on region "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))
                    ENDIF
                    CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                  ENDIF
                ELSE
                  CALL FLAG_ERROR("The elements mesh region does not have any associated nodes",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FLAG_ERROR("The elements mesh region is not associated",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FLAG_ERROR("The elements do not have an associated mesh",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Number of element nodes does not match number of basis nodes for this element",ERR,ERROR,*999)
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
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET")
    RETURN 1
    
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE
    !###  Description:
    !###    Calculates the element numbers surrounding an element in a mesh topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: j,ne,ne1,ne2,nep1,nep2,ni,nn,nn1,nn2,np,np1,np2,DUMMY_ERR,FACE_XI(2),NODE_POSITION_INDEX(3)
    INTEGER(INTG) :: xi_direction,direction_index,xi_dir_check,xi_dir_search,NUMBER_NODE_MATCHES
    INTEGER(INTG) :: NUMBER_SURROUNDING,MAX_NUMBER_SURROUNDING,NUMBER_OF_NODES_XI(3)
    INTEGER(INTG), POINTER :: NODE_MATCHES(:),ADJACENT_ELEMENTS(:)
    LOGICAL :: FOUND,XI_COLLAPSED,FACE_COLLAPSED(-3:3)
    TYPE(VARYING_STRING) :: DUMMY_ERROR
    TYPE(LIST_TYPE), POINTER :: NODE_MATCH_LIST
    TYPE(LIST_PTR_TYPE) :: ADJACENT_ELEMENTS_LIST(-3:3)
    TYPE(BASIS_TYPE), POINTER :: BASIS

    NULLIFY(NODE_MATCHES)
    NULLIFY(BASIS)

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR,*999)
    
    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
          !Loop over the global elements in the mesh
          DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS            
            BASIS=>TOPOLOGY%ELEMENTS%ELEMENTS(ne)%BASIS
            DO ni=-BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
              NULLIFY(ADJACENT_ELEMENTS_LIST(ni)%PTR)
              CALL LIST_CREATE_START(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
              CALL LIST_DATA_TYPE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
              CALL LIST_INITIAL_SIZE_SET(ADJACENT_ELEMENTS_LIST(ni)%PTR,5,ERR,ERROR,*999)
              CALL LIST_CREATE_FINISH(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
            ENDDO !ni
            NUMBER_OF_NODES_XI=1
            NUMBER_OF_NODES_XI(1:BASIS%NUMBER_OF_XI)=BASIS%NUMBER_OF_NODES_XI(1:BASIS%NUMBER_OF_XI)
            !Place the current element in the surrounding list
            CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(0)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_NUMBER,ERR,ERROR,*999)
            MAX_NUMBER_SURROUNDING=1
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
                    DO WHILE(NODE_POSITION_INDEX(xi_dir_search)<=NUMBER_OF_NODES_XI(xi_dir_search).AND.XI_COLLAPSED)
                      !Get the first local node along the xi check direction
                      NODE_POSITION_INDEX(xi_dir_check)=1
                      nn1=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                      !Get the second local node along the xi check direction
                      NODE_POSITION_INDEX(xi_dir_check)=2
                      nn2=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                      IF(nn1/=0.AND.nn2/=0) THEN
                        IF(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn1)/= &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn2)) XI_COLLAPSED=.TRUE.
                      ENDIF
                      NODE_POSITION_INDEX(xi_dir_search)=NODE_POSITION_INDEX(xi_dir_search)+1
                    ENDDO !xi_dir_search
                    IF(XI_COLLAPSED) FACE_COLLAPSED(xi_direction)=.TRUE.
                  ENDIF
                ENDDO !j
                NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XI(ni)
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
                  NODE_POSITION_INDEX(ni)=NUMBER_OF_NODES_XI(ni)
                ENDIF
                !If the face is collapsed then don't look in this xi direction. The exception is if the opposite face is also
                !collpased. This may indicate that we have a funny element in non-rc coordinates that goes around the central
                !axis back to itself
                IF(FACE_COLLAPSED(xi_direction).AND..NOT.FACE_COLLAPSED(-xi_direction)) THEN
                  !Do nothing - the match lists are already empty
                ELSE
                  !Find the nodes to match and add them to the node match list
                  DO nn1=1,NUMBER_OF_NODES_XI(FACE_XI(1))
                    NODE_POSITION_INDEX(FACE_XI(1))=nn1
                    DO nn2=1,NUMBER_OF_NODES_XI(FACE_XI(2))
                      NODE_POSITION_INDEX(FACE_XI(2))=nn2
                      nn=BASIS%NODE_POSITION_INDEX_INV(NODE_POSITION_INDEX(1),NODE_POSITION_INDEX(2),NODE_POSITION_INDEX(3),1)
                      IF(nn/=0) THEN
                        np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)
                        CALL LIST_ITEM_ADD(NODE_MATCH_LIST,np,ERR,ERROR,*999)
                      ENDIF
                    ENDDO !nn2
                  ENDDO !nn1
                ENDIF
                CALL LIST_REMOVE_DUPLICATES(NODE_MATCH_LIST,ERR,ERROR,*999)
                CALL LIST_DETACH(NODE_MATCH_LIST,NUMBER_NODE_MATCHES,NODE_MATCHES,ERR,ERROR,*999)
                NUMBER_SURROUNDING=0
                IF(NUMBER_NODE_MATCHES>0) THEN
                  !Find list of elements surrounding those nodes
                  np1=NODE_MATCHES(1)
                  DO nep1=1,TOPOLOGY%NODES%NODES(np1)%NUMBER_OF_SURROUNDING_ELEMENTS
                    ne1=TOPOLOGY%NODES%NODES(np1)%SURROUNDING_ELEMENTS(nep1)
                    IF(ne1/=ne) THEN !Don't want the current element
                      FOUND=.FALSE.
                      nn2=2
                      DO WHILE(nn2<=NUMBER_NODE_MATCHES.AND..NOT.FOUND)
                        np2=NODE_MATCHES(nn2)
                        nep2=1
                        DO WHILE(nep2<=TOPOLOGY%NODES%NODES(np2)%NUMBER_OF_SURROUNDING_ELEMENTS.AND..NOT.FOUND)
                          ne2=TOPOLOGY%NODES%NODES(np2)%SURROUNDING_ELEMENTS(nep2)
                          IF(ne1==ne2) THEN
                            FOUND=.TRUE.
                          ELSE
                            nep2=nep2+1
                          ENDIF
                        ENDDO !nep2
                        nn2=nn2+1
                      ENDDO !nn2
                      IF(FOUND) THEN
                        CALL LIST_ITEM_ADD(ADJACENT_ELEMENTS_LIST(xi_direction)%PTR,ne1,ERR,ERROR,*999)
                        NUMBER_SURROUNDING=NUMBER_SURROUNDING+1
                      ENDIF
                    ENDIF
                  ENDDO !nep1
                ENDIF
                IF(ASSOCIATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
                CALL LIST_DESTROY(NODE_MATCH_LIST,ERR,ERROR,*999)
                IF(NUMBER_SURROUNDING>MAX_NUMBER_SURROUNDING) MAX_NUMBER_SURROUNDING=NUMBER_SURROUNDING
              ENDDO !direction_index
            ENDDO !ni
            !Set the surrounding elements for this element
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(-BASIS%NUMBER_OF_XI: &
              & BASIS%NUMBER_OF_XI),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of surrounding elements",ERR,ERROR,*999)
            ALLOCATE(TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(MAX_NUMBER_SURROUNDING, &
              & -BASIS%NUMBER_OF_XI:BASIS%NUMBER_OF_XI),STAT=ERR)
            IF(ERR/=0) CALL FLAG_ERROR("Could not allocate surrounding elements",ERR,ERROR,*999)
            DO ni=-BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
              CALL LIST_DETACH(ADJACENT_ELEMENTS_LIST(ni)%PTR,TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & NUMBER_OF_ADJACENT_ELEMENTS(ni),ADJACENT_ELEMENTS,ERR,ERROR,*999)
              TOPOLOGY%ELEMENTS%ELEMENTS(ne)%ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS%ELEMENTS(ne)% &
                & NUMBER_OF_ADJACENT_ELEMENTS(ni),ni) = ADJACENT_ELEMENTS(1:TOPOLOGY%ELEMENTS% &
                & ELEMENTS(ne)%NUMBER_OF_ADJACENT_ELEMENTS(ni))
              IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
              CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,ERR,ERROR,*999)
            ENDDO !ni            
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
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of xi directions = ",BASIS%NUMBER_OF_XI,ERR,ERROR,*999)
        DO ni=-BASIS%NUMBER_OF_XI,BASIS%NUMBER_OF_XI
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
999 IF(ASSOCIATED(NODE_MATCHES)) DEALLOCATE(NODE_MATCHES)
    IF(ASSOCIATED(ADJACENT_ELEMENTS)) DEALLOCATE(ADJACENT_ELEMENTS)
    IF(ASSOCIATED(NODE_MATCH_LIST)) CALL LIST_DESTROY(NODE_MATCH_LIST,DUMMY_ERR,DUMMY_ERROR,*998)
998 DO ni=-3,3
      IF(ASSOCIATED(ADJACENT_ELEMENTS_LIST(ni)%PTR)) CALL LIST_DESTROY(ADJACENT_ELEMENTS_LIST(ni)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
    ENDDO !ni
997 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE")
    RETURN 1   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_ADJACENT_ELEMENTS_CALCULATE
  
  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_FINALISE
    !###  Description:
    !###    Finalises the elements data structures for a mesh topology and deallocates any memory. 

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: ne

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        DO ne=1,TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS
          CALL MESH_TOPOLOGY_ELEMENT_FINALISE(TOPOLOGY%ELEMENTS%ELEMENTS(ne),ERR,ERROR,*999)
        ENDDO !ne
        DEALLOCATE(TOPOLOGY%ELEMENTS)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Topology is not associated",ERR,ERROR,*999)
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_INITIALISE
    !###  Description:
    !###    Initialises the elements in a given mesh topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_SET(GLOBAL_NUMBER,USER_NUMBER,ELEMENTS,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_ELEMENTS_NUMBER_SET
    !###  Description:
    !###    Changes/sets the user number for a global element identified by a given global number.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER,USER_NUMBER
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(ELEMENTS)) THEN
      IF(ELEMENTS%ELEMENTS_FINISHED) THEN
        CALL FLAG_ERROR("Elements have been finished",ERR,ERROR,*999)
      ELSE
        IF(GLOBAL_NUMBER>1.AND.GLOBAL_NUMBER<ELEMENTS%NUMBER_OF_ELEMENTS) THEN
          ELEMENTS%ELEMENTS(GLOBAL_NUMBER)%USER_NUMBER=USER_NUMBER
        ELSE
          LOCAL_ERROR="Global element number "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
            & " is invalid. The limits are 1 to "//TRIM(NUMBER_TO_VSTRING(ELEMENTS%NUMBER_OF_ELEMENTS,"*",ERR,ERROR))
          CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FLAG_ERROR("Elements is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET",ERR,ERROR)    
    CALL EXITS("MESH_TOPOLOGY_ELEMENTS_NUMBER_SET")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_ELEMENTS_NUMBER_SET

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_FINALISE(MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_FINALISE
    !###  Description:
    !###    Finalises the topology in the given mesh.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: component_idx

    CALL ENTERS("MESH_TOPOLOGY_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(MESH)) THEN
      DO component_idx=1,MESH%NUMBER_OF_COMPONENTS
        CALL MESH_TOPOLOGY_NODES_FINALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
        CALL MESH_TOPOLOGY_ELEMENTS_FINALISE(MESH%TOPOLOGY(component_idx)%PTR,ERR,ERROR,*999)
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

  SUBROUTINE MESH_TOPOLOGY_INITIALISE(MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_INITIALISE
    !###  Description:
    !###    Initialises the topology for a given mesh.

    !Argument variables
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_TOPOLOGY_NODE_FINALISE(NODE,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODE_FINALISE
    !###  Description:
    !###    Finalises the given mesh topology node.

    !Argument variables
    TYPE(MESH_NODE_TYPE) :: NODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_NODE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(NODE%SURROUNDING_ELEMENTS)) DEALLOCATE(NODE%SURROUNDING_ELEMENTS)
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

  SUBROUTINE MESH_TOPOLOGY_NODE_INITIALISE(NODE,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODE_INITIALISE
    !###  Description:
    !###    Initialises the given mesh topology node.

    !Argument variables
    TYPE(MESH_NODE_TYPE) :: NODE
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables

    CALL ENTERS("MESH_TOPOLOGY_NODE_INITIALISE",ERR,ERROR,*999)

    NODE%USER_NUMBER=0
    NODE%GLOBAL_NUMBER=0
    NODE%NUMBER_OF_SURROUNDING_ELEMENTS=0
    NULLIFY(NODE%SURROUNDING_ELEMENTS)
    NODE%NUMBER_OF_DERIVATIVES=0
    
    CALL EXITS("MESH_TOPOLOGY_NODE_INITIALISE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODE_INITIALISE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODE_INITIALISE")
    RETURN 1
  END SUBROUTINE MESH_TOPOLOGY_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_NODES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODES_CALCULATE
    !###  Description:
    !###    Calculates the nodes used the mesh identified by MESH. NODES is the pointer to the MESH_NODES data structure.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: NUMBER_OF_MESH_NODES,ne,nn,np,node_no
    INTEGER(INTG), ALLOCATABLE :: MESH_NODES(:)
    LOGICAL :: FOUND
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(MESH_TYPE), POINTER :: MESH
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(MESH_ELEMENTS_TYPE), POINTER :: ELEMENTS
    TYPE(MESH_NODES_TYPE), POINTER :: NODES
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    CALL ENTERS("MESH_TOPOLOGY_NODES_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%ELEMENTS)) THEN
        ELEMENTS=>TOPOLOGY%ELEMENTS
        IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
          NODES=>TOPOLOGY%NODES
          IF(ASSOCIATED(TOPOLOGY%MESH)) THEN
            MESH=>TOPOLOGY%MESH
            IF(ASSOCIATED(MESH%REGION)) THEN
              REGION=>MESH%REGION
              IF(ASSOCIATED(REGION%NODES)) THEN
                IF(ASSOCIATED(TOPOLOGY%NODES%NODES)) THEN
                  LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                    & " already has associated mesh topology nodes"
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ELSE
                  ALLOCATE(MESH_NODES(REGION%NODES%NUMBER_OF_NODES), STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate temporary mesh nodes",ERR,ERROR,*999)
                  NUMBER_OF_MESH_NODES=0
                  DO ne=1,ELEMENTS%NUMBER_OF_ELEMENTS
                    BASIS=>ELEMENTS%ELEMENTS(ne)%BASIS
                    DO nn=1,BASIS%NUMBER_OF_NODES
                      np=ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)
!!TODO: need to be more efficient here
                      FOUND=.FALSE.
                      node_no=1
                      DO WHILE(node_no<=NUMBER_OF_MESH_NODES.AND..NOT.FOUND)
                        IF(MESH_NODES(node_no)==np) THEN
                          FOUND=.TRUE.
                        ELSE
                          node_no=node_no+1
                        ENDIF
                      ENDDO
                      IF(.NOT.FOUND) THEN
                        NUMBER_OF_MESH_NODES=NUMBER_OF_MESH_NODES+1
                        MESH_NODES(NUMBER_OF_MESH_NODES)=np
                      ENDIF
                    ENDDO !nn
                  ENDDO !ne
                  CALL LIST_SORT(MESH_NODES(1:NUMBER_OF_MESH_NODES),ERR,ERROR,*999)
                  ALLOCATE(NODES%NODES(NUMBER_OF_MESH_NODES),STAT=ERR)
                  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate mesh topology nodes nodes",ERR,ERROR,*999)
                  DO np=1,NUMBER_OF_MESH_NODES
                    CALL MESH_TOPOLOGY_NODE_INITIALISE(NODES%NODES(np),ERR,ERROR,*999)
                    NODES%NODES(np)%GLOBAL_NUMBER=MESH_NODES(np)
                    NODES%NODES(np)%USER_NUMBER=REGION%NODES%NODES(MESH_NODES(np))%USER_NUMBER
                  ENDDO !np
                  NODES%NUMBER_OF_NODES=NUMBER_OF_MESH_NODES
                ENDIF
              ELSE
                LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                & " does not have any nodes associated with the mesh region"
                CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              LOCAL_ERROR="Mesh number "//TRIM(NUMBER_TO_VSTRING(MESH%USER_NUMBER,"*",ERR,ERROR))// &
                & " does not have an associated region"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FLAG_ERROR("Mesh topology mesh is not associated",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FLAG_ERROR("Mesh topology nodes is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF
    
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",NODES%NUMBER_OF_NODES,ERR,ERROR,*999)
      DO np=1,NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global node number = ",NODES%NODES(np)%GLOBAL_NUMBER, &
          & ERR,ERROR,*999)        
      ENDDO !np
    ENDIF
    
    CALL EXITS("MESH_TOPOLOGY_NODES_CALCULATE")
    RETURN
999 CALL ERRORS("MESH_TOPOLOGY_NODES_CALCULATE",ERR,ERROR)
    CALL EXITS("MESH_TOPOLOGY_NODES_CALCULATE")
    RETURN 1
   
  END SUBROUTINE MESH_TOPOLOGY_NODES_CALCULATE

  !
  !================================================================================================================================
  !

  SUBROUTINE MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODES_DERIVATIVES_CALCULATE
    !###  Description:
    !###    Calculates the number of derivatives at each node in a topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: elem_idx,MAX_NUMBER_OF_DERIVATIVES,ne,nk,nn,np,NUMBER_OF_DERIVATIVES
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
                IF(ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)==np) THEN
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
                CALL FLAG_ERROR("Could not find local node",ERR,ERROR,*999)
              ENDIF
            ENDDO !elem_idx
            CALL LIST_REMOVE_DUPLICATES(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            CALL LIST_DETACH(NODE_DERIVATIVE_LIST,NUMBER_OF_DERIVATIVES,DERIVATIVES,ERR,ERROR,*999)
            CALL LIST_DESTROY(NODE_DERIVATIVE_LIST,ERR,ERROR,*999)
            IF(NUMBER_OF_DERIVATIVES==MAX_NUMBER_OF_DERIVATIVES) THEN
              NODES%NODES(np)%NUMBER_OF_DERIVATIVES=MAX_NUMBER_OF_DERIVATIVES
              ALLOCATE(NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX(MAX_NUMBER_OF_DERIVATIVES),STAT=ERR)
              IF(ERR/=0) CALL FLAG_ERROR("Could not allocate node partial derivative index",ERR,ERROR,*999)
              NODES%NODES(np)%PARTIAL_DERIVATIVE_INDEX=DERIVATIVES(1:NUMBER_OF_DERIVATIVES)
              DEALLOCATE(DERIVATIVES)
            ELSE
              LOCAL_ERROR="Invalid mesh configuration. User node "// &
                & TRIM(NUMBER_TO_VSTRING(NODES%NODES(np)%USER_NUMBER,"*",ERR,ERROR))// &
                & " has inconsistent derivative directions"
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            ENDIF            
          ENDDO !np
        ELSE
          CALL FLAG_ERROR("Mesh topology nodes is not associated",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FLAG_ERROR("Mesh topology elements is not associated",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Mesh topology is not associated",ERR,ERROR,*999)
    ENDIF
   
    IF(DIAGNOSTICS1) THEN
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Number of mesh global nodes = ",NODES%NUMBER_OF_NODES,ERR,ERROR,*999)
      DO np=1,NODES%NUMBER_OF_NODES
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Mesh global node number = ",np,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of derivatives = ",NODES%NODES(np)%NUMBER_OF_DERIVATIVES, &
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

  SUBROUTINE MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODES_SURROUNDING_ELEMENTS_CALCULATE
    !###  Description:
    !###    Calculates the element numbers surrounding a node for a mesh.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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
                np=TOPOLOGY%ELEMENTS%ELEMENTS(ne)%GLOBAL_ELEMENT_NODES(nn)
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

  SUBROUTINE MESH_TOPOLOGY_NODES_FINALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODES_FINALISE
    !###  Description:
    !###    Finalises the nodes data structures for a mesh topology and deallocates any memory. 

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: np

    CALL ENTERS("MESH_TOPOLOGY_NODES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TOPOLOGY)) THEN
      IF(ASSOCIATED(TOPOLOGY%NODES)) THEN
        DO np=1,TOPOLOGY%NODES%NUMBER_OF_NODES
          CALL MESH_TOPOLOGY_NODE_FINALISE(TOPOLOGY%NODES%NODES(np),ERR,ERROR,*999)
        ENDDO !np
        DEALLOCATE(TOPOLOGY%NODES%NODES)
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

  SUBROUTINE MESH_TOPOLOGY_NODES_INITIALISE(TOPOLOGY,ERR,ERROR,*)

    !#### Subroutine: MESH_TOPOLOGY_NODES_INITIALISE
    !###  Description:
    !###    Initialises the nodes in a given mesh topology.

    !Argument variables
    TYPE(MESH_TOPOLOGY_TYPE), POINTER :: TOPOLOGY
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
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

  SUBROUTINE MESH_USER_NUMBER_FIND(USER_NUMBER,REGION,MESH,ERR,ERROR,*)

    !#### Subroutine: MESH_USER_NUMBER_FIND
    !###  Description:
    !###    Finds and returns in MESH a pointer to the mesh identified by USER_NUMBER in the given REGION. If no mesh with that
    !###    number exits MESH is left nullified.

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER
    TYPE(REGION_TYPE), POINTER :: REGION
    TYPE(MESH_TYPE), POINTER :: MESH
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: mesh_idx
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESH_USER_NUMBER_FIND",ERR,ERROR,*999)

    NULLIFY(MESH)
    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        mesh_idx=1
        DO WHILE(mesh_idx<=REGION%MESHES%NUMBER_OF_MESHES.AND..NOT.ASSOCIATED(MESH))
          IF(REGION%MESHES%MESHES(mesh_idx)%PTR%USER_NUMBER==USER_NUMBER) THEN
            MESH=>REGION%MESHES%MESHES(mesh_idx)%PTR
          ELSE
            mesh_idx=mesh_idx+1
          ENDIF
        ENDDO
      ELSE
        LOCAL_ERROR="The meshes on region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " are not associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESH_USER_NUMBER_FIND")
    RETURN
999 CALL ERRORS("MESH_USER_NUMBER_FIND",ERR,ERROR)
    CALL EXITS("MESH_USER_NUMBER_FIND")
    RETURN 1
  END SUBROUTINE MESH_USER_NUMBER_FIND

  !
  !================================================================================================================================
  !

  SUBROUTINE MESHES_FINALISE(REGION,ERR,ERROR,*)

    !#### Subroutine: MESHES_FINALISE
    !###  Description:
    !###    Finalises the meshes in the given region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    INTEGER(INTG) :: USER_NUMBER

    CALL ENTERS("MESHES_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        DO WHILE(REGION%MESHES%NUMBER_OF_MESHES>0)
          USER_NUMBER=REGION%MESHES%MESHES(1)%PTR%USER_NUMBER
          CALL MESH_DESTROY(USER_NUMBER,REGION,ERR,ERROR,*999)
        ENDDO !mesh_idx
        DEALLOCATE(REGION%MESHES)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
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

  SUBROUTINE MESHES_INITIALISE(REGION,ERR,ERROR,*)

    !#### Subroutine: MESHES_INITIALISE
    !###  Description:
    !###    Initialises the meshes for the given region.

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("MESHES_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%MESHES)) THEN
        LOCAL_ERROR="Region number "//TRIM(NUMBER_TO_VSTRING(REGION%USER_NUMBER,"*",ERR,ERROR))// &
          & " already has a mesh associated"
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ELSE
        ALLOCATE(REGION%MESHES,STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Region meshes could not be allocated",ERR,ERROR,*999)
        REGION%MESHES%NUMBER_OF_MESHES=0
        NULLIFY(REGION%MESHES%MESHES)
        REGION%MESHES%REGION=>REGION
      ENDIF
    ELSE
      CALL FLAG_ERROR("Region is not associated",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("MESHES_INITIALISE")
    RETURN
999 CALL ERRORS("MESHES_INITIALISE",ERR,ERROR)
    CALL EXITS("MESHES_INITIALISE")
    RETURN 1
  END SUBROUTINE MESHES_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE MESH_ROUTINES
