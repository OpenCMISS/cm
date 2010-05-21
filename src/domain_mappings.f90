!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module handles all domain mappings routines.
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

!> This module handles all domain mappings routines.
MODULE DOMAIN_MAPPINGS

  USE BASE_ROUTINES
  USE COMP_ENVIRONMENT
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE STRINGS
  USE TYPES

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup DOMAIN_MAPPINGS_DomainType DOMAIN_MAPPINGS::DomainType
  !> \see DOMAIN_MAPPINGS
  !>@{
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_INTERNAL=1 !<The domain item is internal to the domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_BOUNDARY=2 !<The domain item is on the boundary of the domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  INTEGER(INTG), PARAMETER :: DOMAIN_LOCAL_GHOST=3 !<The domain item is ghosted from another domain \see DOMAIN_MAPPINGS_DomainType,DOMAIN_MAPPINGS
  !>@}
  
  !Module types

  !Module variables

  !Interfaces

  PUBLIC DOMAIN_LOCAL_INTERNAL,DOMAIN_LOCAL_BOUNDARY,DOMAIN_LOCAL_GHOST
  
  PUBLIC DOMAIN_MAPPINGS_MAPPING_FINALISE,DOMAIN_MAPPINGS_MAPPING_INITIALISE,DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE, &
    & DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET,DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE

CONTAINS
  
  !
  !================================================================================================================================
  !

  !>Finalises the adjacent domain and deallocates all memory for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(ADJACENT_DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE) :: ADJACENT_DOMAIN !<The adjacent domain to finalise.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(ADJACENT_DOMAIN%LOCAL_GHOST_SEND_INDICES)) DEALLOCATE(ADJACENT_DOMAIN%LOCAL_GHOST_SEND_INDICES)
    IF(ALLOCATED(ADJACENT_DOMAIN%LOCAL_GHOST_RECEIVE_INDICES)) DEALLOCATE(ADJACENT_DOMAIN%LOCAL_GHOST_RECEIVE_INDICES)
    ADJACENT_DOMAIN%NUMBER_OF_SEND_GHOSTS=0
    ADJACENT_DOMAIN%NUMBER_OF_RECEIVE_GHOSTS=0
    ADJACENT_DOMAIN%DOMAIN_NUMBER=0
    
    CALL EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialise the adjacent domain for a domain mapping.
  SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(ADJACENT_DOMAIN,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_ADJACENT_DOMAIN_TYPE) :: ADJACENT_DOMAIN !<The adjacent domain to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",ERR,ERROR,*999)

    ADJACENT_DOMAIN%NUMBER_OF_SEND_GHOSTS=0
    ADJACENT_DOMAIN%NUMBER_OF_RECEIVE_GHOSTS=0
    ADJACENT_DOMAIN%DOMAIN_NUMBER=0
    
    CALL EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Returns the local number, if it exists on the rank, for the specifed global number
  SUBROUTINE DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET(DOMAIN_MAPPING,GLOBAL_NUMBER,LOCAL_EXISTS,LOCAL_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to get the local number from
    INTEGER(INTG), INTENT(IN) :: GLOBAL_NUMBER !<The global number to get the local number for
    LOGICAL, INTENT(OUT) :: LOCAL_EXISTS !<On exit, is .TRUE. if the specifed global number exists on the local rank, .FALSE. if not
    INTEGER(INTG), INTENT(OUT) :: LOCAL_NUMBER !<On exit, the local number corresponding to the global number if it exists. If it doesn't exist then 0.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET",ERR,ERROR,*999)

    LOCAL_EXISTS=.FALSE.
    LOCAL_NUMBER=0
    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(GLOBAL_NUMBER>=1.AND.GLOBAL_NUMBER<=DOMAIN_MAPPING%NUMBER_OF_GLOBAL) THEN
        IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_NUMBER)%DOMAIN_NUMBER(1)== &
          & COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER) THEN
          LOCAL_NUMBER=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(GLOBAL_NUMBER)%LOCAL_NUMBER(1)
          LOCAL_EXISTS=.TRUE.
        ENDIF
      ELSE
        LOCAL_ERROR="The specified global number of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_NUMBER,"*",ERR,ERROR))// &
          & " is invalid. The number must be between 1 and "// &
          & TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%NUMBER_OF_GLOBAL,"*",ERR,ERROR))//"."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_GLOBAL_TO_LOCAL_GET
  
  !
  !================================================================================================================================
  !

  !>Calculates the domain mappings local map from a domain mappings global map.
  SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE(DOMAIN_MAPPING,ERR,ERROR,*)

   !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<The domain mapping to calculate the local mappings
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: domain_idx,domain_idx2,domain_no,domain_no2,global_number,idx,local_number,local_number2,NUMBER_INTERNAL, &
      & NUMBER_BOUNDARY,NUMBER_GHOST,my_computational_node_number,MY_DOMAIN_INDEX,TEMP,NUMBER_OF_ADJACENT_DOMAINS, &
      & RECEIVE_FROM_DOMAIN,DUMMY_ERR,NUMBER_OF_GHOST_RECEIVE,NUMBER_OF_GHOST_SEND,local_type,COUNT, &
      & TOTAL_NUMBER_OF_ADJACENT_DOMAINS
    INTEGER(INTG), ALLOCATABLE :: ADJACENT_DOMAIN_MAP(:),ADJACENT_DOMAINS(:,:)
    INTEGER(INTG), POINTER :: SEND_LIST(:),RECEIVE_LIST(:)
    LOGICAL :: OWNED_BY_ALL,SEND_GLOBAL
    TYPE(LIST_PTR_TYPE), ALLOCATABLE :: GHOST_SEND_LISTS(:),GHOST_RECEIVE_LISTS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR,DUMMY_ERROR

    NULLIFY(SEND_LIST)
    NULLIFY(RECEIVE_LIST)

    CALL ENTERS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      my_computational_node_number=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      IF(ERR/=0) GOTO 999        
      !Calculate local to global maps from global to local map
      ALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of domain local.",ERR,ERROR,*999)
      DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL=0
      ALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate number of domain ghost.",ERR,ERROR,*999)
      DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST=0
      NUMBER_INTERNAL=0
      NUMBER_BOUNDARY=0
      NUMBER_GHOST=0
      ALLOCATE(ADJACENT_DOMAINS(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1,0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent domains.",ERR,ERROR,*999)
      ADJACENT_DOMAINS=0
      DO global_number=1,DOMAIN_MAPPING%NUMBER_OF_GLOBAL
        !If necessary, reset global domain index so that my computational node is in the first index position
        MY_DOMAIN_INDEX=1
        DO domain_idx=2,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          IF(domain_no==my_computational_node_number) THEN
            MY_DOMAIN_INDEX=domain_idx
            EXIT
          ENDIF
        ENDDO !domain_idx
        IF(MY_DOMAIN_INDEX/=1) THEN !Swap domain index in the global to local map
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(MY_DOMAIN_INDEX) = TEMP
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(MY_DOMAIN_INDEX) = TEMP
          TEMP=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(1)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(1) = &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(MY_DOMAIN_INDEX)
          DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(MY_DOMAIN_INDEX) = TEMP
        ENDIF
        DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          DO domain_idx2=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
            domain_no2=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx2)
            ADJACENT_DOMAINS(domain_no,domain_no2)=1
          ENDDO !domain_idx2
        ENDDO !domain_idx
        DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
          IF(local_type==DOMAIN_LOCAL_GHOST) THEN
            DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST(domain_no)=DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST(domain_no)+1
          ELSE
            DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(domain_no)=DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL(domain_no)+1
          ENDIF
          IF(domain_no==my_computational_node_number) THEN
            SELECT CASE(local_type)
            CASE(DOMAIN_LOCAL_INTERNAL)
              NUMBER_INTERNAL=NUMBER_INTERNAL+1
            CASE(DOMAIN_LOCAL_BOUNDARY)
              NUMBER_BOUNDARY=NUMBER_BOUNDARY+1
            CASE(DOMAIN_LOCAL_GHOST)
              NUMBER_GHOST=NUMBER_GHOST+1
            CASE DEFAULT
              LOCAL_ERROR="The domain local type of "//TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                & global_number)%LOCAL_TYPE(domain_idx),"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ENDIF
        ENDDO !domain_idx
      ENDDO !global_number

      !!TODO: move adjacent domains calculation back to where the global to local array is set up????
      NUMBER_OF_ADJACENT_DOMAINS=0
      TOTAL_NUMBER_OF_ADJACENT_DOMAINS=0
      DO domain_no=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
        DO domain_no2=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
          IF(domain_no/=domain_no2) THEN 
            IF(ADJACENT_DOMAINS(domain_no,domain_no2)>0) THEN
              TOTAL_NUMBER_OF_ADJACENT_DOMAINS=TOTAL_NUMBER_OF_ADJACENT_DOMAINS+1
              IF(domain_no==my_computational_node_number) NUMBER_OF_ADJACENT_DOMAINS=NUMBER_OF_ADJACENT_DOMAINS+1
            ENDIF
          ENDIF
        ENDDO !domain_no2
      ENDDO !domain_no
      ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent domains ptr.",ERR,ERROR,*999)
      ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST(TOTAL_NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent domains list.",ERR,ERROR,*999)
      COUNT=1
      DO domain_no=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
        DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(domain_no)=COUNT
        DO domain_no2=0,DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1
          IF(domain_no/=domain_no2) THEN
            IF(ADJACENT_DOMAINS(domain_no,domain_no2)>0) THEN
              DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST(COUNT)=domain_no2
              COUNT=COUNT+1
            ENDIF
          ENDIF
        ENDDO !domain_no2
      ENDDO !domain_no
      DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(DOMAIN_MAPPING%NUMBER_OF_DOMAINS)=COUNT
      DEALLOCATE(ADJACENT_DOMAINS)
      
      ALLOCATE(DOMAIN_MAPPING%DOMAIN_LIST(NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain map domain list.",ERR,ERROR,*999)
      ALLOCATE(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate domain map local to global list.",ERR,ERROR,*999)
      DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL=NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST
      DOMAIN_MAPPING%NUMBER_OF_LOCAL=NUMBER_INTERNAL+NUMBER_BOUNDARY
      DOMAIN_MAPPING%NUMBER_OF_INTERNAL=NUMBER_INTERNAL
      DOMAIN_MAPPING%NUMBER_OF_BOUNDARY=NUMBER_BOUNDARY
      DOMAIN_MAPPING%NUMBER_OF_GHOST=NUMBER_GHOST
      DOMAIN_MAPPING%INTERNAL_START=1
      DOMAIN_MAPPING%INTERNAL_FINISH=NUMBER_INTERNAL
      DOMAIN_MAPPING%BOUNDARY_START=NUMBER_INTERNAL+1
      DOMAIN_MAPPING%BOUNDARY_FINISH=NUMBER_INTERNAL+NUMBER_BOUNDARY
      DOMAIN_MAPPING%GHOST_START=NUMBER_INTERNAL+NUMBER_BOUNDARY+1
      DOMAIN_MAPPING%GHOST_FINISH=NUMBER_INTERNAL+NUMBER_BOUNDARY+NUMBER_GHOST
      NUMBER_INTERNAL=0
      NUMBER_BOUNDARY=0
      NUMBER_GHOST=0
      ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS(NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent domains.",ERR,ERROR,*999)
      DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=NUMBER_OF_ADJACENT_DOMAINS
      ALLOCATE(ADJACENT_DOMAIN_MAP(0:DOMAIN_MAPPING%NUMBER_OF_DOMAINS-1),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate adjacent domain map.",ERR,ERROR,*999)
      ALLOCATE(GHOST_SEND_LISTS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ghost send list.",ERR,ERROR,*999)
      ALLOCATE(GHOST_RECEIVE_LISTS(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS),STAT=ERR)
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ghost recieve list.",ERR,ERROR,*999)
      DO domain_idx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_INITIALISE(DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx),ERR,ERROR,*999)
        domain_no= &
          & DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR(my_computational_node_number)+domain_idx-1)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER=domain_no
        ADJACENT_DOMAIN_MAP(domain_no)=domain_idx
        NULLIFY(GHOST_SEND_LISTS(domain_idx)%PTR)
        CALL LIST_CREATE_START(GHOST_SEND_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(GHOST_SEND_LISTS(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(GHOST_SEND_LISTS(domain_idx)%PTR,MAX(DOMAIN_MAPPING%NUMBER_OF_GHOST,1),ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(GHOST_SEND_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        NULLIFY(GHOST_RECEIVE_LISTS(domain_idx)%PTR)
        CALL LIST_CREATE_START(GHOST_RECEIVE_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DATA_TYPE_SET(GHOST_RECEIVE_LISTS(domain_idx)%PTR,LIST_INTG_TYPE,ERR,ERROR,*999)
        CALL LIST_INITIAL_SIZE_SET(GHOST_RECEIVE_LISTS(domain_idx)%PTR,MAX(DOMAIN_MAPPING%NUMBER_OF_GHOST,1),ERR,ERROR,*999)
        CALL LIST_CREATE_FINISH(GHOST_RECEIVE_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
      ENDDO !domain_idx
      DO global_number=1,DOMAIN_MAPPING%NUMBER_OF_GLOBAL
        SEND_GLOBAL=.FALSE.
        IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS>1) THEN
          !Check if we have a special case where the global number is owned by all domains e.g., as in a constant field
          IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS==DOMAIN_MAPPING%NUMBER_OF_DOMAINS) THEN
            OWNED_BY_ALL=.TRUE.
            DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
              local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
              OWNED_BY_ALL=OWNED_BY_ALL.AND.local_type==DOMAIN_LOCAL_INTERNAL
            ENDDO !domain_idx
          ELSE
            OWNED_BY_ALL=.FALSE.            
          ENDIF
          IF(.NOT.OWNED_BY_ALL) THEN
            RECEIVE_FROM_DOMAIN=-1
            DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
              domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
              local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
              IF(local_type/=DOMAIN_LOCAL_GHOST) THEN
                IF(domain_no==my_computational_node_number) SEND_GLOBAL=.TRUE.
                IF(RECEIVE_FROM_DOMAIN==-1) THEN
                  RECEIVE_FROM_DOMAIN=domain_no
                ELSE
                  LOCAL_ERROR="Invalid domain mapping. Global number "//TRIM(NUMBER_TO_VSTRING(global_number,"*",ERR,ERROR))// &
                    & " is owned by domain number "//TRIM(NUMBER_TO_VSTRING(RECEIVE_FROM_DOMAIN,"*",ERR,ERROR))// &
                    & " as well as domain number "//TRIM(NUMBER_TO_VSTRING(domain_no,"*",ERR,ERROR))//"."
                  CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDDO !domain_idx
            IF(RECEIVE_FROM_DOMAIN==-1) THEN
              LOCAL_ERROR="Invalid domain mapping. Global number "//TRIM(NUMBER_TO_VSTRING(global_number,"*",ERR,ERROR))// &
                & " is not owned by any domain."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)          
            ENDIF
          ENDIF
        ENDIF
        DO domain_idx=1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%NUMBER_OF_DOMAINS
          domain_no=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%DOMAIN_NUMBER(domain_idx)
          local_number=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(domain_idx)
          local_type=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx)
          IF(domain_no==my_computational_node_number) THEN
            DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(local_number)=global_number
            SELECT CASE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_TYPE(domain_idx))
            CASE(DOMAIN_LOCAL_INTERNAL)
              NUMBER_INTERNAL=NUMBER_INTERNAL+1
              DOMAIN_MAPPING%DOMAIN_LIST(NUMBER_INTERNAL)=local_number
            CASE(DOMAIN_LOCAL_BOUNDARY)
              NUMBER_BOUNDARY=NUMBER_BOUNDARY+1
              DOMAIN_MAPPING%DOMAIN_LIST(DOMAIN_MAPPING%INTERNAL_FINISH+NUMBER_BOUNDARY)=local_number
            CASE(DOMAIN_LOCAL_GHOST)
              NUMBER_GHOST=NUMBER_GHOST+1
              DOMAIN_MAPPING%DOMAIN_LIST(DOMAIN_MAPPING%BOUNDARY_FINISH+NUMBER_GHOST)=local_number
              CALL LIST_ITEM_ADD(GHOST_RECEIVE_LISTS(ADJACENT_DOMAIN_MAP(RECEIVE_FROM_DOMAIN))%PTR,local_number,ERR,ERROR,*999)
            CASE DEFAULT
              LOCAL_ERROR="The domain local type of "//TRIM(NUMBER_TO_VSTRING(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP( &
                & global_number)%LOCAL_TYPE(domain_idx),"*",ERR,ERROR))//" is invalid."
              CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          ELSE IF(SEND_GLOBAL.AND.local_type==DOMAIN_LOCAL_GHOST) THEN
            local_number2=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(global_number)%LOCAL_NUMBER(1) !The local number for this node
            CALL LIST_ITEM_ADD(GHOST_SEND_LISTS(ADJACENT_DOMAIN_MAP(domain_no))%PTR,local_number2,ERR,ERROR,*999)
          ENDIF
        ENDDO !domain_idx
      ENDDO !global_number
      
      DO domain_idx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
        CALL LIST_REMOVE_DUPLICATES(GHOST_SEND_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(GHOST_SEND_LISTS(domain_idx)%PTR,NUMBER_OF_GHOST_SEND,SEND_LIST,ERR,ERROR,*999)
        ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES(NUMBER_OF_GHOST_SEND),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local ghost send inidices.",ERR,ERROR,*999)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES=SEND_LIST(1:NUMBER_OF_GHOST_SEND)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS=NUMBER_OF_GHOST_SEND
        DEALLOCATE(SEND_LIST)
        CALL LIST_REMOVE_DUPLICATES(GHOST_RECEIVE_LISTS(domain_idx)%PTR,ERR,ERROR,*999)
        CALL LIST_DETACH_AND_DESTROY(GHOST_RECEIVE_LISTS(domain_idx)%PTR,NUMBER_OF_GHOST_RECEIVE,RECEIVE_LIST,ERR,ERROR,*999)
        ALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES(NUMBER_OF_GHOST_RECEIVE),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate local ghost receive inidices.",ERR,ERROR,*999)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES=RECEIVE_LIST(1:NUMBER_OF_GHOST_RECEIVE)
        DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS=NUMBER_OF_GHOST_RECEIVE
        DEALLOCATE(RECEIVE_LIST)        
      ENDDO !domain_idx

      DEALLOCATE(ADJACENT_DOMAIN_MAP)
      DEALLOCATE(GHOST_SEND_LISTS)
      DEALLOCATE(GHOST_RECEIVE_LISTS)
      
      IF(DIAGNOSTICS1) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Domain mappings:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of domains  = ",DOMAIN_MAPPING%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of global = ",DOMAIN_MAPPING%NUMBER_OF_GLOBAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of local = ",DOMAIN_MAPPING%NUMBER_OF_LOCAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Total number of local = ",DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL, &
          & ERR,ERROR,*999)
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Domain numbers:",ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%NUMBER_OF_DOMAINS,8,8,DOMAIN_MAPPING% &
          & NUMBER_OF_DOMAIN_LOCAL,'("    Number of domain local :",8(X,I10))','(26X,8(X,I10))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%NUMBER_OF_DOMAINS,8,8,DOMAIN_MAPPING% &
          & NUMBER_OF_DOMAIN_GHOST,'("    Number of domain ghost :",8(X,I10))','(26X,8(X,I10))',ERR,ERROR,*999)      
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Domain list:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of internal = ",DOMAIN_MAPPING%NUMBER_OF_INTERNAL,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of boundary = ",DOMAIN_MAPPING%NUMBER_OF_BOUNDARY,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of ghost = ",DOMAIN_MAPPING%NUMBER_OF_GHOST,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Internal start = ",DOMAIN_MAPPING%INTERNAL_START,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Internal finish = ",DOMAIN_MAPPING%INTERNAL_FINISH,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Boundary start = ",DOMAIN_MAPPING%BOUNDARY_START,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Boundary finish = ",DOMAIN_MAPPING%BOUNDARY_FINISH,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Ghost start = ",DOMAIN_MAPPING%GHOST_START,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Ghost finish = ",DOMAIN_MAPPING%GHOST_FINISH,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,DOMAIN_MAPPING%INTERNAL_START,1,DOMAIN_MAPPING%INTERNAL_FINISH,8,8, &
          & DOMAIN_MAPPING%DOMAIN_LIST,'("    Internal list :",8(X,I10))','(19X,8(X,I10))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,DOMAIN_MAPPING%BOUNDARY_START,1,DOMAIN_MAPPING%BOUNDARY_FINISH,8,8, &
          & DOMAIN_MAPPING%DOMAIN_LIST,'("    Boundary list :",8(X,I10))','(19X,8(X,I10))',ERR,ERROR,*999)      
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,DOMAIN_MAPPING%GHOST_START,1,DOMAIN_MAPPING%GHOST_FINISH,8,8, &
          & DOMAIN_MAPPING%DOMAIN_LIST,'("    Ghost list    :",8(X,I10))','(19X,8(X,I10))',ERR,ERROR,*999)      
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Local to global map:",ERR,ERROR,*999)
        DO idx=1,DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Local index : ",idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Global index = ",DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(idx), &
            & ERR,ERROR,*999)
        ENDDO !idx
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Global to local map:",ERR,ERROR,*999)
        DO idx=1,DOMAIN_MAPPING%NUMBER_OF_GLOBAL
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Global idx : ",idx,ERR,ERROR,*999)
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of domains  = ", &
            & DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx)%NUMBER_OF_DOMAINS,ERR,ERROR,*999)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx)% &
            & NUMBER_OF_DOMAINS,8,8,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx)%LOCAL_NUMBER, &
            & '("      Local number  :",8(X,I10))','(21X,8(X,I10))',ERR,ERROR,*999)      
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx)% &
            & NUMBER_OF_DOMAINS,8,8,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx)%DOMAIN_NUMBER, &
            & '("      Domain number :",8(X,I10))','(21X,8(X,I10))',ERR,ERROR,*999)      
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx)% &
            & NUMBER_OF_DOMAINS,8,8,DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(IDX)%LOCAL_TYPE, &
            & '("      Local type    :",8(X,I10))','(21X,8(X,I10))',ERR,ERROR,*999)      
        ENDDO !ne
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Adjacent domains:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of adjacent domains = ", &
          & DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS,ERR,ERROR,*999)
        CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%NUMBER_OF_DOMAINS+1,8,8, &
          & DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR,'("    Adjacent domains ptr  :",8(X,I5))','(27X,8(X,I5))',ERR,ERROR,*999)
        IF(DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS>0) THEN
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR( &
            & DOMAIN_MAPPING%NUMBER_OF_DOMAINS)-1,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST, &
            '("    Adjacent domains list :",8(X,I5))','(27X,8(X,I5))',ERR,ERROR,*999)
          DO domain_idx=1,DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Adjacent domain idx : ",domain_idx,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Domain number = ", &
              & DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%DOMAIN_NUMBER,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of send ghosts    = ", &
              & DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_SEND_GHOSTS,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
              & NUMBER_OF_SEND_GHOSTS,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_SEND_INDICES, &
              & '("      Local send ghost indicies       :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)      
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"      Number of recieve ghosts = ", &
              & DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%NUMBER_OF_RECEIVE_GHOSTS,ERR,ERROR,*999)
            CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)% &
              & NUMBER_OF_RECEIVE_GHOSTS,8,8,DOMAIN_MAPPING%ADJACENT_DOMAINS(domain_idx)%LOCAL_GHOST_RECEIVE_INDICES, &
              & '("      Local receive ghost indicies    :",8(X,I10))','(39X,8(X,I10))',ERR,ERROR,*999)              
          ENDDO !domain_idx
        ENDIF
      ENDIF
      
    ELSE
      CALL FLAG_ERROR("Domain mapping is not associated.",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE")
    RETURN
999 IF(ASSOCIATED(SEND_LIST)) DEALLOCATE(SEND_LIST)
    IF(ASSOCIATED(RECEIVE_LIST)) DEALLOCATE(RECEIVE_LIST)
    IF(ALLOCATED(ADJACENT_DOMAIN_MAP)) DEALLOCATE(ADJACENT_DOMAIN_MAP)
    IF(ALLOCATED(ADJACENT_DOMAINS)) DEALLOCATE(ADJACENT_DOMAINS)
    IF(ALLOCATED(GHOST_SEND_LISTS)) THEN
      DO domain_idx=1,SIZE(GHOST_SEND_LISTS)
        IF(ASSOCIATED(GHOST_SEND_LISTS(domain_idx)%PTR)) &
          & CALL LIST_DESTROY(GHOST_SEND_LISTS(domain_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*998)
      ENDDO ! domain_idx
998   DEALLOCATE(GHOST_SEND_LISTS)
    ENDIF
    IF(ALLOCATED(GHOST_RECEIVE_LISTS)) THEN
      DO domain_idx=1,SIZE(GHOST_RECEIVE_LISTS)
        IF(ASSOCIATED(GHOST_RECEIVE_LISTS(domain_idx)%PTR)) &
          & CALL LIST_DESTROY(GHOST_RECEIVE_LISTS(domain_idx)%PTR,DUMMY_ERR,DUMMY_ERROR,*997)
      ENDDO ! domain_idx
997   DEALLOCATE(GHOST_RECEIVE_LISTS)
    ENDIF
    CALL ERRORS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE")
    RETURN 1
    
  END SUBROUTINE DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE
  
  !
  !================================================================================================================================
  !

  !>Finalises the mapping for a domain mappings mapping and deallocates all memory.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE(DOMAIN_MAPPING,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: idx

    CALL ENTERS("DOMAIN_MAPPINGS_MAPPING_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(ALLOCATED(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL)) DEALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_LOCAL)
      IF(ALLOCATED(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST)) DEALLOCATE(DOMAIN_MAPPING%NUMBER_OF_DOMAIN_GHOST)
      IF(ALLOCATED(DOMAIN_MAPPING%DOMAIN_LIST)) DEALLOCATE(DOMAIN_MAPPING%DOMAIN_LIST)
      IF(ALLOCATED(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP)) DEALLOCATE(DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP)
      IF(ALLOCATED(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP,1)
          CALL DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(idx),ERR,ERROR,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP)
      ENDIF
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR)) DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_PTR)
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST)) DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS_LIST)
      IF(ALLOCATED(DOMAIN_MAPPING%ADJACENT_DOMAINS)) THEN
        DO idx=1,SIZE(DOMAIN_MAPPING%ADJACENT_DOMAINS,1)
          CALL DOMAIN_MAPPINGS_ADJACENT_DOMAIN_FINALISE(DOMAIN_MAPPING%ADJACENT_DOMAINS(idx),ERR,ERROR,*999)
        ENDDO !idx
        DEALLOCATE(DOMAIN_MAPPING%ADJACENT_DOMAINS)
      ENDIF
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_MAPPING_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_FINALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_FINALISE
  
  !
  !================================================================================================================================
  !

  !> Finalises the global mapping in the given domain mappings.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE(MAPPING_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: MAPPING_GLOBAL_MAP !<The domain global mapping to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<Th error string
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE",ERR,ERROR,*999)

    IF(ALLOCATED(MAPPING_GLOBAL_MAP%LOCAL_NUMBER)) DEALLOCATE(MAPPING_GLOBAL_MAP%LOCAL_NUMBER)
    IF(ALLOCATED(MAPPING_GLOBAL_MAP%DOMAIN_NUMBER)) DEALLOCATE(MAPPING_GLOBAL_MAP%DOMAIN_NUMBER)
    IF(ALLOCATED(MAPPING_GLOBAL_MAP%LOCAL_TYPE)) DEALLOCATE(MAPPING_GLOBAL_MAP%LOCAL_TYPE)
 
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_FINALISE

  !
  !================================================================================================================================
  !

  !>Finalises the global mapping in the given domain mappings.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE(MAPPING_GLOBAL_MAP,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_GLOBAL_MAPPING_TYPE) :: MAPPING_GLOBAL_MAP !<The domain global mapping to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE",ERR,ERROR,*999)

    MAPPING_GLOBAL_MAP%NUMBER_OF_DOMAINS=0
    
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE")
    RETURN 1
   
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_GLOBAL_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises the mapping for a domain mappings mapping.
  SUBROUTINE DOMAIN_MAPPINGS_MAPPING_INITIALISE(DOMAIN_MAPPING,NUMBER_OF_DOMAINS,ERR,ERROR,*)

    !Argument variables
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING !<A pointer to the domain mapping to initialise the mappings for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DOMAINS !<The number of domains 
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    CALL ENTERS("DOMAIN_MAPPINGS_MAPPING_INITIALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(DOMAIN_MAPPING)) THEN
      IF(NUMBER_OF_DOMAINS>0) THEN
        DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL=0
        DOMAIN_MAPPING%NUMBER_OF_LOCAL=0
        DOMAIN_MAPPING%NUMBER_OF_GLOBAL=0
        DOMAIN_MAPPING%NUMBER_OF_DOMAINS=NUMBER_OF_DOMAINS
        DOMAIN_MAPPING%NUMBER_OF_INTERNAL=0
        DOMAIN_MAPPING%NUMBER_OF_BOUNDARY=0
        DOMAIN_MAPPING%NUMBER_OF_GHOST=0
        DOMAIN_MAPPING%INTERNAL_START=0
        DOMAIN_MAPPING%INTERNAL_FINISH=0
        DOMAIN_MAPPING%BOUNDARY_START=0
        DOMAIN_MAPPING%BOUNDARY_FINISH=0
        DOMAIN_MAPPING%GHOST_START=0
        DOMAIN_MAPPING%GHOST_FINISH=0
        DOMAIN_MAPPING%NUMBER_OF_ADJACENT_DOMAINS=0
      ELSE
        LOCAL_ERROR="The specified number of domains of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DOMAINS,"*",ERR,ERROR))// &
          & " is invalid. The number of domains must be > 0."
        CALL FLAG_ERROR(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ENDIF
    
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_INITIALISE")
    RETURN
999 CALL ERRORS("DOMAIN_MAPPINGS_MAPPING_INITIALISE",ERR,ERROR)
    CALL EXITS("DOMAIN_MAPPINGS_MAPPING_INITIALISE")
    RETURN 1
  END SUBROUTINE DOMAIN_MAPPINGS_MAPPING_INITIALISE

  !
  !================================================================================================================================
  !
  
END MODULE DOMAIN_MAPPINGS
