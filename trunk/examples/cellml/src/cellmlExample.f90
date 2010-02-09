!> \file
!> $Id: cellml.f90 $
!> \author David Nickerson <nickerso@users.sourceforge.net>
!> \brief An example deomonstrating and testing the implementation of the openCMISS(cellml) API.
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

!> Doxygen comments get ignored when writing this out in the documentation

! Uses EXAMPLE_PATH from Doxyfile but putting path from examples folder to ensure uniqueness.
!> \example cellml/src/cellmlExample.f90
!! A complete example demonstrating and testing the methods defined in the openCMISS(cellml) API.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/cellml/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/cellml/build-gnu'>Linux GNU Build</a>
!<

! An example application demonstrating the complete openCMISS(cellml) API.
! This example is designed to both provide a guide to the use of CellML models 
! in openCMISS and to test (guide?) the openCMISS(cellml) API implementation.
! In this example we aim to define a regular mesh, define a field ....
PROGRAM CELLMLEXAMPLE

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE CMISS
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MESH_ROUTINES
  USE MPI
  USE PROBLEM_CONSTANTS
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

  IMPLICIT NONE

  !Program parameters
  !==================

  ! The mesh size (mm) that will be used.
  REAL(DP), PARAMETER :: HEIGHT = 1.0_DP
  REAL(DP), PARAMETER :: WIDTH = 1.0_DP
  REAL(DP), PARAMETER :: LENGTH = 1.0_DP

  !Program types
  !=============
  TYPE MY_MESH_DATA
     TYPE(BASIS_TYPE), POINTER :: BASIS
     TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
     TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
     TYPE(MESH_TYPE), POINTER :: MESH
     TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
     TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD
     TYPE(REGION_TYPE), POINTER :: REGION,WORLD_REGION
  END TYPE MY_MESH_DATA
  
  !Program variables
  !=================
  TYPE(MY_MESH_DATA), POINTER :: MESH_DATA
  INTEGER(INTG) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
       & NUMBER_GLOBAL_Z_ELEMENTS,NUMBER_OF_DOMAINS
  INTEGER(INTG), PARAMETER :: field_geometry_id=1
  INTEGER(INTG), PARAMETER :: field_geometry_number_of_variables=1
  INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
  INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
  INTEGER(INTG) :: MPI_IERROR
  LOGICAL :: EXPORT_FIELD
  TYPE(VARYING_STRING) :: FILE,METHOD
  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1), &
       & STOP_SYSTEM_TIME(1)
  INTEGER(INTG) i,count,status,uriL
  CHARACTER(256) BUFFER


  !Generic CMISS variables
  !=======================
  
  ! The error flag and string
  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR

  ! Intialise CMISS
  NULLIFY(MESH_DATA%WORLD_REGION)
  CALL CMISS_INITIALISE(MESH_DATA%WORLD_REGION,ERR,ERROR,*999)

  !Calculate the start times
  CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR,ERROR,*999)

  !Get the number of computational nodes
  NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  !Get my computational node number
  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  
  NUMBER_GLOBAL_X_ELEMENTS=3
  NUMBER_GLOBAL_Y_ELEMENTS=3
  NUMBER_GLOBAL_Z_ELEMENTS=3
  NUMBER_OF_DOMAINS=2

  IF (MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
     !GET THE PARAMETERS FROM THE COMMAND LINE
     IF (COMMAND_ARGUMENT_COUNT() == 4) THEN
        CALL GET_COMMAND_ARGUMENT(1,BUFFER)
        READ(BUFFER,*) NUMBER_OF_DOMAINS
        CALL GET_COMMAND_ARGUMENT(2,BUFFER)
        READ(BUFFER,*) NUMBER_GLOBAL_X_ELEMENTS
        CALL GET_COMMAND_ARGUMENT(3,BUFFER)
        READ(BUFFER,*) NUMBER_GLOBAL_Y_ELEMENTS
        CALL GET_COMMAND_ARGUMENT(4,BUFFER)
        READ(BUFFER,*) NUMBER_GLOBAL_Z_ELEMENTS
     ELSE
        CALL FLAG_ERROR("Invalid command line arguments.",ERR, &
             & ERROR,*999)
     END IF
  END IF

  !Broadcast the number of elements and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"COMPUTATIONAL ENVIRONMENT:",ERR,ERROR,*999)
  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE, &
       & "  Total number of computational nodes = ", &
       & NUMBER_COMPUTATIONAL_NODES,ERR,ERROR,*999)
  CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE, &
       & "  My computational node number = ", &
       & MY_COMPUTATIONAL_NODE_NUMBER,ERR,ERROR,*999)

  

  ! Finalise CMISS
  CALL CMISS_FINALISE(ERR,ERROR,*999)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
  STOP

END PROGRAM CELLMLEXAMPLE
