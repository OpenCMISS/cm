!> \file
!> $Id$
!> \author Chris Bradley
!> \brief This module contains all computational environment variables.
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

!> This module contains all computational environment variables.
MODULE COMP_ENVIRONMENT
  
  USE BASE_ROUTINES
  USE CMISS_MPI
  USE CMISS_PETSC
  USE CONSTANTS
  USE KINDS
  USE MPI
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  
  IMPLICIT NONE
 
  PRIVATE

  !Module parameters

  !Module types

  !>Contains information on a cache heirarchy
  TYPE CACHE_TYPE
    INTEGER(INTG) :: NUMBER_LEVELS !<The number of levels in the cache hierarchy
    INTEGER(INTG),ALLOCATABLE :: SIZE(:) !<SIZE(level_idx). The size of the level_idx'th cache level.
  END TYPE CACHE_TYPE

  !>Contains information on a computational node containing a number of processors
  TYPE COMPUTATIONAL_NODE_TYPE
    INTEGER(INTG) :: NUMBER_PROCESSORS !<The number of processors for this computational node
    INTEGER(INTG) :: RANK !<The MPI rank of this computational node
!    TYPE(CACHE_TYPE) :: CACHE 
    INTEGER(INTG) :: NODE_NAME_LENGTH !<The length of the name of the computational node
    CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: NODE_NAME !<The name of the computational node
  END TYPE COMPUTATIONAL_NODE_TYPE

  !>Contains information on the MPI type to transfer information about a computational node
  TYPE MPI_COMPUTATIONAL_NODE_TYPE
    INTEGER(INTG) :: MPI_TYPE !<The MPI data type
    INTEGER(INTG) :: NUM_BLOCKS !<The number of blocks in the MPI data type. This will be equal to 4.
    INTEGER(INTG) :: BLOCK_LENGTHS(4) !<The length of each block.
    INTEGER(INTG) :: TYPES(4) !<The data types of each block.
    INTEGER(MPI_ADDRESS_KIND) :: DISPLACEMENTS(4) !<The address displacements to each block.
  END TYPE MPI_COMPUTATIONAL_NODE_TYPE

  !>Contains information on the computational environment the program is running in.
  TYPE COMPUTATIONAL_ENVIRONMENT_TYPE
    INTEGER(INTG) :: MPI_COMM !<The MPI communicator for cmiss
    INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES !<The number of computational nodes
    INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER !<The index of the running process
    TYPE(COMPUTATIONAL_NODE_TYPE), ALLOCATABLE :: COMPUTATIONAL_NODES(:) !<COMPUTATIONAL_NODES(node_idx). Contains information on the node_idx'th computational node. 
  END TYPE COMPUTATIONAL_ENVIRONMENT_TYPE

  !Module variables

  TYPE(COMPUTATIONAL_ENVIRONMENT_TYPE) :: COMPUTATIONAL_ENVIRONMENT !<The computational environment the program is running in.
  TYPE(MPI_COMPUTATIONAL_NODE_TYPE) :: MPI_COMPUTATIONAL_NODE_TYPE_DATA !<The MPI data on the computational nodes.

  !Interfaces

  PUBLIC COMPUTATIONAL_ENVIRONMENT_TYPE

  PUBLIC COMPUTATIONAL_ENVIRONMENT

  PUBLIC COMPUTATIONAL_ENVIRONMENT_INITIALISE,COMPUTATIONAL_ENVIRONMENT_FINALISE,COMPUTATIONAL_NODES_NUMBER_GET, &
    & COMPUTATIONAL_NODE_NUMBER_GET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finalises the computational node data structures and deallocates all memory.
  SUBROUTINE COMPUTATIONAL_NODE_FINALISE(COMPUTATIONAL_NODE,ERR,ERROR,*)
  
    !Argument Variables
    TYPE(COMPUTATIONAL_NODE_TYPE), INTENT(INOUT) :: COMPUTATIONAL_NODE !<The computational node to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    CALL ENTERS("COMPUTATIONAL_NODE_FINALISE",ERR,ERROR,*999)

    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=0
    COMPUTATIONAL_NODE%RANK=-1
    COMPUTATIONAL_NODE%NODE_NAME_LENGTH=0
    COMPUTATIONAL_NODE%NODE_NAME=""    

    CALL EXITS("COMPUTATIONAL_NODE_FINALISE")
    RETURN
999 CALL ERRORS("COMPUTATIONAL_NODE_FINALISE",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_NODE_FINALISE")
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the computational node data structures.
  SUBROUTINE COMPUTATIONAL_NODE_INITIALISE(COMPUTATIONAL_NODE,RANK,ERR,ERROR,*)
  
    !Argument Variables
    TYPE(COMPUTATIONAL_NODE_TYPE), INTENT(OUT) :: COMPUTATIONAL_NODE !<The computational node to initialise
    INTEGER(INTG), INTENT(IN) :: RANK !<The MPI rank of the computational node
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: MPI_IERROR

    CALL ENTERS("COMPUTATIONAL_NODE_INITIALISE",ERR,ERROR,*999)

!    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=COMP_DETECT_NUMBER_PROCESSORS(ERR)
!    IF(ERR/=0) GOTO 999
    COMPUTATIONAL_NODE%NUMBER_PROCESSORS=1
    COMPUTATIONAL_NODE%RANK=RANK
    CALL MPI_GET_PROCESSOR_NAME(COMPUTATIONAL_NODE%NODE_NAME,COMPUTATIONAL_NODE%NODE_NAME_LENGTH,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_PROCESSOR_NAME",MPI_IERROR,ERR,ERROR,*999)
    
    CALL EXITS("COMPUTATIONAL_NODE_INITIALISE")
    RETURN
999 CALL ERRORS("COMPUTATIONAL_NODE_INITIALISE",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_NODE_INITIALISE")
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalises the data structure containing the MPI type information for the COMPUTATIONAL_NODE_TYPE.
  SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(ERR,ERROR,*)
  
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,MPI_IERROR

    CALL ENTERS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE",ERR,ERROR,*999)

    DO i=1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES(i)=0
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS(i)=0
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(i)=0
    ENDDO !i
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS=0

    IF(MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE/=MPI_DATATYPE_NULL) THEN
      CALL MPI_TYPE_FREE(MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_TYPE_FREE",MPI_IERROR,ERR,ERROR,*999)
    ENDIF

    CALL EXITS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE")
    RETURN
999 CALL ERRORS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_NODE_MPI_TYPE_FINALISE")
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the data structure containing the MPI type information for the COMPUTATIONAL_NODE_TYPE.
  SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE(COMPUTATIONAL_NODE,ERR,ERROR,*)
  
    !Argument Variables
    TYPE(COMPUTATIONAL_NODE_TYPE), INTENT(IN) :: COMPUTATIONAL_NODE !<The computational node containing the MPI type to initialise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: I,MPI_IERROR

    CALL ENTERS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE",ERR,ERROR,*999)

    MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE=MPI_DATATYPE_NULL
    
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS=4
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES=(/MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_CHARACTER/)
    MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS=(/1,1,1,MPI_MAX_PROCESSOR_NAME/)
	
	
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NUMBER_PROCESSORS,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(1),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%RANK,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(2),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NODE_NAME_LENGTH,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(3),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)
    !CPB 19/02/07 AIX compiler complains about the type of the first parameter i.e., the previous 3 have been integers
    !and this one is not so cast the type.
    CALL MPI_GET_ADDRESS(COMPUTATIONAL_NODE%NODE_NAME,MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(4),MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_GET_ADDRESS",MPI_IERROR,ERR,ERROR,*999)

    DO i=4,1,-1
      MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(I)=MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(I)- &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS(1)
    ENDDO !i

    CALL MPI_TYPE_CREATE_STRUCT(MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS, &
      & MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS,MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES, &
      & MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_CREATE_STRUCT",MPI_IERROR,ERR,ERROR,*999)

    CALL MPI_TYPE_COMMIT(MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_TYPE_COMMIT",MPI_IERROR,ERR,ERROR,*999)
    
    IF(DIAGNOSTICS3) THEN
      CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"MPI Computational Node Type Data:",ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  MPI type = ",MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE,ERR,ERROR,*999)
      CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number blocks  = ",MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS, &
        & ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,4,4, &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%TYPES,'("  Block types =",4(X,I15))','(15X,4(X,I15))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,8,8, &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%BLOCK_LENGTHS,'("  Block lengths =",8(X,I5))','(17X,8(X,I5))',ERR,ERROR,*999)
      CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%NUM_BLOCKS,8,8, &
        & MPI_COMPUTATIONAL_NODE_TYPE_DATA%DISPLACEMENTS,'("  Displacements =",8(X,I5))','(17X,8(X,I5))',ERR,ERROR,*999)
    ENDIF

    CALL EXITS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE")
    RETURN
999 CALL COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(ERR,ERROR,*998)
998 CALL ERRORS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE")
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>
  !>Finalises the computational environment data structures and deallocates all memory.
  !>
  SUBROUTINE COMPUTATIONAL_ENVIRONMENT_FINALISE(ERR,ERROR,*)
  
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: COMPUTATIONAL_NODE,MPI_IERROR

    CALL ENTERS("COMPUTATIONAL_ENVIRONMENT_FINALISE",ERR,ERROR,*999)

    !Finalise PetSc
    !CALL PETSC_LOGPRINTSUMMARY(PETSC_COMM_WORLD,"OpenCMISSTest.petsc",ERR,ERROR,*999)
    CALL PETSC_FINALIZE(ERR,ERROR,*999)
    
    IF(ALLOCATED(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)) THEN
       DO COMPUTATIONAL_NODE=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
          CALL COMPUTATIONAL_NODE_FINALISE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(COMPUTATIONAL_NODE),ERR,ERROR,*999)
       ENDDO
       DEALLOCATE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)
    ENDIF
    COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES=0

    CALL COMPUTATIONAL_NODE_MPI_TYPE_FINALISE(ERR,ERROR,*999)

    CALL MPI_COMM_FREE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_FREE",MPI_IERROR,ERR,ERROR,*999)
    
    CALL MPI_FINALIZE(MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_FINALIZE",MPI_IERROR,ERR,ERROR,*999)

    CALL EXITS("COMPUTATIONAL_ENVIRONMENT_FINALISE")
    RETURN
999 CALL ERRORS("COMPUTATIONAL_ENVIRONMENT_FINALISE",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_ENVIRONMENT_FINALISE")
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_ENVIRONMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>
  !>Initialises the computational environment data structures.
  !>
  SUBROUTINE COMPUTATIONAL_ENVIRONMENT_INITIALISE(ERR,ERROR,*)
  
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: i,DUMMY_ERR,MPI_IERROR,RANK
    TYPE(VARYING_STRING) :: DUMMY_ERROR

    CALL ENTERS("COMPUTATIONAL_ENVIRONMENT_INITIALISE",ERR,ERROR,*999)

    !Initialise the MPI environment
    CALL MPI_INIT(MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_INIT",MPI_IERROR,ERR,ERROR,*999)

    !Create a (private) communicator for cmiss. For now just duplicate MPI_COMM_WORLD
    CALL MPI_COMM_DUP(MPI_COMM_WORLD,COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_DUP",MPI_IERROR,ERR,ERROR,*999)
    
    !Determine the number of ranks/computational nodes we have in our computational environment
    CALL MPI_COMM_SIZE(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_SIZE",MPI_IERROR,ERR,ERROR,*999)

    !Allocate the computational node data structures
    ALLOCATE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0:COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1),STAT=ERR)
    IF(ERR /=0) CALL FLAG_ERROR("Could not allocate computational nodes",ERR,ERROR,*999)

    !Determine my processes rank
    CALL MPI_COMM_RANK(COMPUTATIONAL_ENVIRONMENT%MPI_COMM,RANK,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_COMM_RANK",MPI_IERROR,ERR,ERROR,*999)
    COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER=RANK
    
    !Create the MPI type information for the COMPUTATIONAL_NODE_TYPE
    CALL COMPUTATIONAL_NODE_MPI_TYPE_INITIALISE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(RANK),ERR,ERROR,*999)
    !Fill in all the computational node data structures for this rank at the root position (will be changed later with an
    !allgather call)
    CALL COMPUTATIONAL_NODE_INITIALISE(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),RANK,ERR,ERROR,*999)

    !Now transfer all the computational node information to the other computational nodes so that each rank has all the
    !information.
    CALL MPI_ALLGATHER(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, &
      & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(0),1,MPI_COMPUTATIONAL_NODE_TYPE_DATA%MPI_TYPE, &
      & COMPUTATIONAL_ENVIRONMENT%MPI_COMM,MPI_IERROR)
    CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)

    !Initialise PETSc
    CALL PETSC_INITIALIZE(PETSC_NULL_CHARACTER,ERR,ERROR,*999)
    
    IF(DIAGNOSTICS1) THEN
      !Just let the master node write out this information
      IF(RANK==0) THEN
        CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"Computational environment:",ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  Number of computational nodes = ", &
          & COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  My computational node number = ", &
          & COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER,ERR,ERROR,*999)
        IF(DIAGNOSTICS2) THEN
          DO i=0,COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES-1
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"  Computational Node:",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Number of Processors = ", &
              & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(i)%NUMBER_PROCESSORS,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    MPI rank = ", &
              & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(i)%RANK,ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    Node Name = ", &
              & COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES(i)%NODE_NAME,ERR,ERROR,*999)
          ENDDO !i
        ENDIF
      ENDIF
    ENDIF

    CALL EXITS("COMPUTATIONAL_ENVIRONMENT_INITIALISE")
    RETURN
999 CALL COMPUTATIONAL_ENVIRONMENT_FINALISE(DUMMY_ERR,DUMMY_ERROR,*998)
998 CALL ERRORS("COMPUTATIONAL_ENVIRONMENT_INITIALISE",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_ENVIRONMENT_INITIALISE")
    RETURN 1
  END SUBROUTINE COMPUTATIONAL_ENVIRONMENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the number/rank of the computational nodes.
  FUNCTION COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
      
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: COMPUTATIONAL_NODE_NUMBER_GET !<On exit the computational node number/rank of the current process.
    !Local Variables

    CALL ENTERS("COMPUTATIONAL_NODE_NUMBER_GET",ERR,ERROR,*999)

    IF(ALLOCATED(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)) THEN
      COMPUTATIONAL_NODE_NUMBER_GET=COMPUTATIONAL_ENVIRONMENT%MY_COMPUTATIONAL_NODE_NUMBER
    ELSE
      CALL FLAG_ERROR("Computational environment not initialised",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COMPUTATIONAL_NODE_NUMBER_GET")
    RETURN
999 CALL ERRORS("COMPUTATIONAL_NODE_NUMBER_GET",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_NODE_NUMBER_GET")
    RETURN 
  END FUNCTION COMPUTATIONAL_NODE_NUMBER_GET

  !
  !================================================================================================================================
  !
  
  !>Returns the number of computational nodes.
  FUNCTION COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
     
    !Argument Variables
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Function variable
    INTEGER(INTG) :: COMPUTATIONAL_NODES_NUMBER_GET !<On exit, the number of computational nodes for the program.
    !Local Variables

    CALL ENTERS("COMPUTATIONAL_NODES_NUMBER_GET",ERR,ERROR,*999)

    IF(ALLOCATED(COMPUTATIONAL_ENVIRONMENT%COMPUTATIONAL_NODES)) THEN
      COMPUTATIONAL_NODES_NUMBER_GET=COMPUTATIONAL_ENVIRONMENT%NUMBER_COMPUTATIONAL_NODES
    ELSE
      CALL FLAG_ERROR("Computational environment not initialised",ERR,ERROR,*999)
    ENDIF
    
    CALL EXITS("COMPUTATIONAL_NODES_NUMBER_GET")
    RETURN
999 CALL ERRORS("COMPUTATIONAL_NODES_NUMBER_GET",ERR,ERROR)
    CALL EXITS("COMPUTATIONAL_NODES_NUMBER_GET")
    RETURN 
  END FUNCTION COMPUTATIONAL_NODES_NUMBER_GET

  !
  !================================================================================================================================
  !

END MODULE COMP_ENVIRONMENT
