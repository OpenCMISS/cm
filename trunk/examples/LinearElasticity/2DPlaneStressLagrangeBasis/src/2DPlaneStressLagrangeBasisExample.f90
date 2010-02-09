!> \file
!> $Id: LinearElasticityExample.f90 20 2009-02-15 13:26:52Z cpb $
!> \author Chris Bradley
!> \brief This is an example program to solve a linear elasticity equation using openCMISS calls.
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

!> \example LinearElasticity/src/LinearElasticityExample.f90
!! Example program to solve a linear elasticity equation using openCMISS calls.
!<

!> Main program
PROGRAM LINEARELASTICITYEXAMPLE

  USE BASE_ROUTINES
  USE BASIS_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CMISS
  USE CMISS_MPI
  USE COMP_ENVIRONMENT
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DISTRIBUTED_MATRIX_VECTOR
  USE DOMAIN_MAPPINGS
  USE EQUATIONS_ROUTINES
  USE EQUATIONS_SET_CONSTANTS
  USE EQUATIONS_SET_ROUTINES
  USE FIELD_ROUTINES
  USE FIELD_IO_ROUTINES
  USE LINEAR_ELASTICITY_ROUTINES  
  USE GENERATED_MESH_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE LISTS
  USE MESH_ROUTINES
  USE MPI
  USE NODE_ROUTINES  
  USE PROBLEM_CONSTANTS
  USE PROBLEM_ROUTINES
  USE REGION_ROUTINES
  USE SOLVER_ROUTINES
  USE TIMER
  USE TYPES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters
!#### Index: ne
!###  Description:
!###    Index label for a element.
!#### Index: ng
!###  Description:
!###    Index label for a gauss point.
!#### Index: ni
!###  Description:
!###    Index label for a xi direction.
!#### Index: nk
!###  Description:
!###    Index label for a derivative with respect to the global directions.
!#### Index: nn
!###  Description:
!###    Index for a local node within an element.
!#### Index: np
!###  Description:
!###    Index for a node.
!#### Index: ns
!###  Description:
!###    Index for a element parameter within an element.
!#### Index: nu
!###  Description:
!###    Index for a partial derivative.

  !Program types
  TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
  TYPE(REGION_TYPE), POINTER :: REGION,WORLD_REGION
  TYPE(BASIS_PTR_TYPE), POINTER :: XYZ_BASES(:)
  TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
  TYPE(NODES_TYPE), POINTER :: NODES
  TYPE(MESH_TYPE), POINTER :: MESH
  TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
  TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,DEPENDENT_FIELD,MATERIAL_FIELD
  TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
  TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
  TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
  TYPE(SOLVER_TYPE), POINTER :: SOLVER
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
  TYPE(VARYING_STRING) :: FILE,METHOD

  !Element Node Types for prescribing Nodal positions
  !ELEMENTS(:) 					#TYPE(ELEMENT_TYPE)
  !  MESH_COMP(:) 				#TYPE(ELEMENT_NODES_TYPE)
  !    ELEMENT					#TYPE(MESH_ELEMENTS_TYPE)
  !    ELEM_NODE_NUMBERS(:)		#INTEGER(INTG),ALLOCATABLE
  !    DERIV(8)					#TYPE(NODE_COORDINATES_TYPE)
  !      NODE_COORDINATES(:)	#REAL(DP),ALLOCATABLE 
  TYPE NODE_COORDINATES_TYPE
    REAL(DP),ALLOCATABLE :: NODE_COORDINATES(:)
  END TYPE NODE_COORDINATES_TYPE
  TYPE ELEMENT_NODES_TYPE
    TYPE(MESH_ELEMENTS_TYPE),POINTER :: ELEMENT
    INTEGER(INTG),ALLOCATABLE :: ELEM_NODE_NUMBERS(:)
    TYPE(NODE_COORDINATES_TYPE) :: DERIV(8)
  END TYPE ELEMENT_NODES_TYPE
  TYPE ELEMENT_TYPE
    TYPE(ELEMENT_NODES_TYPE),POINTER :: MESH_COMP(:)
  END TYPE ELEMENT_TYPE
  TYPE(ELEMENT_TYPE), POINTER :: ELEMENTS(:)

  !Boundary Condition Types
  !DISP_BC 								#TYPE(BC_MESH_COMP)
  !  MESH_COMP(:) 						#TYPE(BC_ELEMENT_NODES_TYPE)
  !    NUMBER_OF_BC_NODES_IN_DERIV(8)	#INTEGER(INTG)
  !    DERIV(8)							#TYPE(BC_NODE_COORDINATES_TYPE)
  !      NODE_NUMBER(:)					#INTEGER(INTG),ALLOCATABLE
  !      NODE_COORDINATES(:)			#REAL(DP),ALLOCATABLE 
  TYPE BC_NODE_COORDINATES_TYPE
    INTEGER(INTG),ALLOCATABLE :: NODE_NUMBER(:)
    REAL(DP),ALLOCATABLE :: NODE_COORDINATES(:)
  END TYPE BC_NODE_COORDINATES_TYPE
  TYPE BC_ELEMENT_NODES_TYPE
    INTEGER(INTG) :: NUMBER_OF_BC_NODES_IN_DERIV(8)
    TYPE(BC_NODE_COORDINATES_TYPE) :: DERIV(8)
  END TYPE BC_ELEMENT_NODES_TYPE
  TYPE BC_MESH_COMP
    TYPE(BC_ELEMENT_NODES_TYPE),POINTER :: MESH_COMP(:)
  END TYPE BC_MESH_COMP
  TYPE(BC_MESH_COMP), POINTER :: DISP_BC,FORCE_BC

  !Program variables
  INTEGER(INTG) :: NUMBER_OF_GLOBAL_X_ELEMENTS,NUMBER_OF_GLOBAL_Y_ELEMENTS,NUMBER_OF_GLOBAL_Z_ELEMENTS
  INTEGER(INTG) :: NUMBER_OF_DOMAINS,NUMBER_COMPUTATIONAL_NODES,MY_COMPUTATIONAL_NODE_NUMBER,MPI_IERROR
  INTEGER(INTG) :: CS_USER_NUMBER, REGION_USER_NUMBER, BASIS_USER_NUMBER, MESH_USER_NUMBER, DECOMPOSITION_USER_NUMBER
  INTEGER(INTG) :: FIELD_USER_NUMBER,EQUATION_SET_USER_NUMBER, PROBLEM_USER_NUMBER
  INTEGER(INTG) :: NUMBER_OF_MESH_COMPS,NUMBER_OF_XI,NUMBER_OF_NODES,NUMBER_OF_ELEMENTS,NUMBER_OF_DERIV
  INTEGER(INTG) :: NUMBER_OF_FIELD_VARIABLES,NUMBER_OF_BC_NODES,NUMBER_OF_FIELD_COMPS
  INTEGER(INTG) :: MESH_COMP_IDX,FIELD_COMP_IDX,FIELD_DERIVATIVE_IDX,SOLVER_IDX,EQUATION_SET_IDX,np,xi,ne,nu
  INTEGER(INTG),ALLOCATABLE :: XI_INTERPOLATION(:,:),NUMBER_OF_GAUSS_POINTS(:),MESH_COMP_NUMBER_OF_ELEMENT_NODES(:)

  INTEGER(INTG) :: NUMBER_OF_GLOBAL_DEPENDENT_DOFS
  REAL(DP), POINTER :: FIELD_DATA(:)

  REAL(DP) :: l,w,h
  REAL(DP),ALLOCATABLE :: MATE_PARA(:)
  LOGICAL :: Export_Field

  !Timing Variables
  REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)

  !Generic CMISS variables
  INTEGER(INTG) :: ERR
  TYPE(VARYING_STRING) :: ERROR

  INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
  CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(3),TIMING_ROUTINE_LIST(1)

  !WIN32 Variables
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !=========================================1=======================================================================================

  !Intialise cmiss
  NULLIFY(WORLD_REGION)
  CALL CMISS_INITIALISE(WORLD_REGION,ERR,ERROR,*999)

  !Set all diganostic levels on for testing
  DIAG_LEVEL_LIST(1)=1
  DIAG_LEVEL_LIST(2)=2
  DIAG_LEVEL_LIST(3)=3
  DIAG_LEVEL_LIST(4)=4
  DIAG_LEVEL_LIST(5)=5
  !DIAG_ROUTINE_LIST(1)="MATRIX_VALUES_ADD_DP2"
  !DIAG_ROUTINE_LIST(2)="EQUATIONS_MATRIX_STRUCTURE_CALCULATE"
  !DIAG_ROUTINE_LIST(3)="EQUATIONS_MAPPING_CALCULATE"
  !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"LinearElasticity",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
  !TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
  !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)

  !Calculate the start times
  CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR,ERROR,*999)

  !=BROADCAST PARAMETERS TO COMPUTATIONAL NODES====================================================================================
  !Get the number of computational nodes
  NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  !Get my computational node number
  MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
  IF(ERR/=0) GOTO 999
  NUMBER_OF_GLOBAL_X_ELEMENTS = 1
  NUMBER_OF_GLOBAL_Y_ELEMENTS = 1
  NUMBER_OF_GLOBAL_Z_ELEMENTS = 1
  NUMBER_OF_DOMAINS = 1
  NUMBER_OF_XI = 2
  NUMBER_OF_ELEMENTS = 1
  NUMBER_OF_NODES = 8

  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_OF_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_OF_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_OF_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)

  !=CREATE COORDINATE SYSTEM=======================================================================================================
  !Start the creation of a new 3D RC coordinate system
  CS_USER_NUMBER = 1
  NULLIFY(COORDINATE_SYSTEM)
  CALL COORDINATE_SYSTEM_CREATE_START(CS_USER_NUMBER,COORDINATE_SYSTEM,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_TYPE_SET(COORDINATE_SYSTEM,COORDINATE_RECTANGULAR_CARTESIAN_TYPE,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,NUMBER_OF_XI,ERR,ERROR,*999)
  CALL COORDINATE_SYSTEM_ORIGIN_SET(COORDINATE_SYSTEM,(/0.0_DP,0.0_DP,0.0_DP/),ERR,ERROR,*999) 
  CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)

  !=CREATE REGION==================================================================================================================
  !Create Region and set CS to newly created 3D RC CS
  REGION_USER_NUMBER = 1
  NULLIFY(REGION)
  CALL REGION_CREATE_START(REGION_USER_NUMBER,WORLD_REGION,REGION,ERR,ERROR,*999)
  CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
  CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)

  !=CREATE BASIS ==================================================================================================================
  NULLIFY(XYZ_BASES)
  ALLOCATE(XYZ_BASES(NUMBER_OF_XI),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate XYZ_BASES",ERR,ERROR,*999)
  ALLOCATE(XI_INTERPOLATION(NUMBER_OF_XI,NUMBER_OF_XI),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate XI_INTERPOLATION",ERR,ERROR,*999)
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Interpolation in each direction
  !NOTE if you change interpolation you need to change Boundary Conditions
  XI_INTERPOLATION(1,:) = (/BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_LINEAR_LAGRANGE_INTERPOLATION/)
  XI_INTERPOLATION(2,:) = (/BASIS_LINEAR_LAGRANGE_INTERPOLATION,BASIS_LINEAR_LAGRANGE_INTERPOLATION/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ALLOCATE(NUMBER_OF_GAUSS_POINTS(NUMBER_OF_XI),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate NUMBER_OF_GAUSS_POINTS",ERR,ERROR,*999)
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Number of Gauss Points
  !NOTE:: Num of Gauss points must be the same across X,Y & Z coordinates and be sufficient to accurately integrate the hightest order interpolation being used
  NUMBER_OF_GAUSS_POINTS(:) = (/4,4,4/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    BASIS_USER_NUMBER = xi
    NULLIFY(XYZ_BASES(xi)%PTR)
    CALL BASIS_CREATE_START(BASIS_USER_NUMBER,XYZ_BASES(xi)%PTR,ERR,ERROR,*999) 
    CALL BASIS_TYPE_SET(XYZ_BASES(xi)%PTR,BASIS_LAGRANGE_HERMITE_TP_TYPE,ERR,ERROR,*999)
    CALL BASIS_NUMBER_OF_XI_SET(XYZ_BASES(xi)%PTR,NUMBER_OF_XI,ERR,ERROR,*999)
    CALL BASIS_INTERPOLATION_XI_SET(XYZ_BASES(xi)%PTR,XI_INTERPOLATION(xi,:),ERR,ERROR,*999)
    CALL BASIS_QUADRATURE_NUMBER_OF_GAUSS_XI_SET(XYZ_BASES(xi)%PTR,NUMBER_OF_GAUSS_POINTS,ERR,ERROR,*999)  
    CALL BASIS_CREATE_FINISH(XYZ_BASES(xi)%PTR,ERR,ERROR,*999)
  ENDDO !xi

  !=MANUAL MESH CREATION===========================================================================================================
  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE," *** USING MANUAL MESH CREATION ***",ERR,ERROR,*999)

  !=NODE CREATION==================================================================================================================
  !Create nodes in REGION and set initial coordinates to 0,0,0
  NULLIFY(NODES)
  CALL NODES_CREATE_START(REGION,NUMBER_OF_NODES,NODES,ERR,ERROR,*999)
  CALL NODES_CREATE_FINISH(NODES,ERR,ERROR,*999)
  
  !=CREATE MESH====================================================================================================================
  !Create a mesh with xi number of coordinates
  MESH_USER_NUMBER = 1
  NUMBER_OF_MESH_COMPS = NUMBER_OF_XI
  NULLIFY(MESH)

  CALL MESH_CREATE_START(MESH_USER_NUMBER,REGION,NUMBER_OF_XI,MESH,ERR,ERROR,*999)
  CALL MESH_NUMBER_OF_ELEMENTS_SET(MESH,NUMBER_OF_ELEMENTS,ERR,ERROR,*999)
  CALL MESH_NUMBER_OF_COMPONENTS_SET(MESH,NUMBER_OF_MESH_COMPS,ERR,ERROR,*999)

  ALLOCATE(ELEMENTS(NUMBER_OF_ELEMENTS),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENTS",ERR,ERROR,*999)
  DO ne=1,NUMBER_OF_ELEMENTS
    ALLOCATE(ELEMENTS(ne)%MESH_COMP(NUMBER_OF_MESH_COMPS),STAT=ERR)
    IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT(ne)%MESH_COMP",ERR,ERROR,*999)
  ENDDO !ne
  ALLOCATE(MESH_COMP_NUMBER_OF_ELEMENT_NODES(NUMBER_OF_MESH_COMPS),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate MESH_COMP_NUMBER_OF_ELEMENT_NODES",ERR,ERROR,*999)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    IF (XI_INTERPOLATION(1,1)<=3) THEN
      MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)=PRODUCT(XI_INTERPOLATION(xi,:)+1)
      NUMBER_OF_DERIV = 1
    ELSEIF (XI_INTERPOLATION(1,1)==4) THEN
      MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)=(2**NUMBER_OF_XI)
      NUMBER_OF_DERIV = 2**NUMBER_OF_XI
    ELSE
      CALL FLAG_ERROR("This example is only setup for Lagrange and cubic hermite bases.",ERR,ERROR,*999)
    ENDIF
  ENDDO !xi
  DO ne=1,NUMBER_OF_ELEMENTS
    DO xi=1,NUMBER_OF_XI
      MESH_COMP_IDX=xi
      ALLOCATE(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)),STAT=ERR)
      !!TODO:: Add Num to String
      IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS",ERR,ERROR,*999)
      ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS = 0
    ENDDO !xi
  ENDDO !ne
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Element node numbers
  ELEMENTS(1)%MESH_COMP(1)%ELEM_NODE_NUMBERS(:) = (/1,2,3,4/)
  ELEMENTS(1)%MESH_COMP(2)%ELEM_NODE_NUMBERS(:) = (/1,2,3,4/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !!TODO:: When diagnostics are on - MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET trys to output MESH%TOPOLOGY(1)%PTR%ELEMENTS%ELEMENTS(1)%MESH_ELEMENT_NODES which are only allocated when the MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH command is given
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    DO ne=1,NUMBER_OF_ELEMENTS
      NULLIFY(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT)
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_START(MESH,MESH_COMP_IDX,XYZ_BASES(xi)%PTR, & 
        & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT,ERR,ERROR,*999)
      CALL MESH_TOPOLOGY_ELEMENTS_ELEMENT_NODES_SET(ne,ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT, & 
        & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(:),ERR,ERROR,*999)
      CALL MESH_TOPOLOGY_ELEMENTS_CREATE_FINISH(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEMENT,ERR,ERROR,*999)
    ENDDO !ne
  ENDDO !xi
  CALL MESH_CREATE_FINISH(MESH,ERR,ERROR,*999) 

  !=CREATE DECOMPOSITION===========================================================================================================
  !Create mesh decomposition dividing mesh into number_of_domains for parallel solving
  DECOMPOSITION_USER_NUMBER = 1
  NULLIFY(DECOMPOSITION)
  CALL DECOMPOSITION_CREATE_START(DECOMPOSITION_USER_NUMBER,MESH,DECOMPOSITION,ERR,ERROR,*999)
  CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
  CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*999)
  CALL DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*999)

  !=CREATE GEOMETRIC FIELD=========================================================================================================
  !Start to create a default (geometric) field on the region
  NULLIFY(GEOMETRIC_FIELD)
  FIELD_USER_NUMBER = 1
  NUMBER_OF_FIELD_VARIABLES = 1 !Geometric Field Coordinates
  NUMBER_OF_FIELD_COMPS = NUMBER_OF_XI
  CALL FIELD_CREATE_START(FIELD_USER_NUMBER,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
  CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
  CALL FIELD_TYPE_SET(GEOMETRIC_FIELD,FIELD_GEOMETRIC_TYPE,ERR,ERROR,*999)  
  CALL FIELD_NUMBER_OF_VARIABLES_SET(GEOMETRIC_FIELD,NUMBER_OF_FIELD_VARIABLES,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_FIELD_COMPS,ERR,ERROR,*999) 
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    FIELD_COMP_IDX = xi
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_COMP_IDX,MESH_COMP_IDX, &
      & ERR,ERROR,*999)
  ENDDO !xi
  CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)
  !Initialize node coordinates
  DO ne=1,NUMBER_OF_ELEMENTS
    DO xi=1,NUMBER_OF_XI
      MESH_COMP_IDX=xi
      DO nu=1,NUMBER_OF_DERIV
        ALLOCATE(ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES&
          &(MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX)),STAT=ERR)
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate ELEMENT(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES", &
          & ERR,ERROR,*999)
        ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES = 0.0_DP
      ENDDO
    ENDDO
  ENDDO
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe node coordinates
  l = 120.0_DP
  w = 160.0_DP
  h = 10.0_DP
  ELEMENTS(1)%MESH_COMP(1)%DERIV(1)%NODE_COORDINATES = (/0.0_DP,l,0.0_DP,l/)
  ELEMENTS(1)%MESH_COMP(2)%DERIV(1)%NODE_COORDINATES = (/0.0_DP,0.0_DP,w,w/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    FIELD_COMP_IDX = xi
    DO ne=1,NUMBER_OF_ELEMENTS
      DO nu=1,NUMBER_OF_DERIV
        DO np=1,MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX) !Can also loop using nu
          CALL FIELD_PARAMETER_SET_UPDATE_NODE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,nu, &
            & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(np),FIELD_COMP_IDX, &
            & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(np),ERR,ERROR,*999)
        ENDDO !np
      ENDDO !nu
    ENDDO !ne
  ENDDO !xi

  !=CREATE DEPENDENT FIELD=========================================================================================================
  !Create a dependent field
  NULLIFY(DEPENDENT_FIELD)
  FIELD_USER_NUMBER = 2
  NUMBER_OF_FIELD_VARIABLES = 2 !Dependent Field Displacement Coordinates & Force 
  NUMBER_OF_FIELD_COMPS = NUMBER_OF_XI
  CALL FIELD_CREATE_START(FIELD_USER_NUMBER,REGION,DEPENDENT_FIELD,ERR,ERROR,*999)
  CALL FIELD_TYPE_SET(DEPENDENT_FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
  CALL FIELD_MESH_DECOMPOSITION_SET(DEPENDENT_FIELD,DECOMPOSITION,ERR,ERROR,*999)
  CALL FIELD_GEOMETRIC_FIELD_SET(DEPENDENT_FIELD,GEOMETRIC_FIELD,ERR,ERROR,*999)
  CALL FIELD_DEPENDENT_TYPE_SET(DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_VARIABLES_SET(DEPENDENT_FIELD,NUMBER_OF_FIELD_VARIABLES,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_FIELD_COMPS,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_FIELD_COMPS,ERR,ERROR,*999)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    FIELD_COMP_IDX = xi
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_COMP_IDX,MESH_COMP_IDX, &
      & ERR,ERROR,*999)
    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_COMP_IDX,MESH_COMP_IDX, &
      & ERR,ERROR,*999)
  ENDDO !xi
  CALL FIELD_CREATE_FINISH(DEPENDENT_FIELD,ERR,ERROR,*999) 

  !=CREATE MATERIAL FIELD==========================================================================================================
  !!TODO:: Set Material Field Interpolation to constant element based interpolation when field i/o and cmgui allows for this
  !Create a material field for a general 2D isotropic material
  NULLIFY(MATERIAL_FIELD)
  FIELD_USER_NUMBER = 3
  NUMBER_OF_FIELD_VARIABLES = 1
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe material properties E1,v13
  !!TODO::Bring 2D plane stress/strain element thickness in through a field - element constant when it can be exported by field i/o & cmgui allows for this. Currently brought in through material field (Temporary)
  !!Update ELASTICITY_EQUATIONS_SET_SETUP>>LINEAR_ELASTICITY_EQUATIONS_SET_SETUP>>FIELD_NUMBER_OF_COMPONENTS_CHECK once it is changed so for 2D max material components allowed is 2
  NUMBER_OF_FIELD_COMPS = 3 !Young's Modulus & Poisson's Ratio & Thickness
  ALLOCATE(MATE_PARA(NUMBER_OF_FIELD_COMPS),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate MATE_PARA",ERR,ERROR,*999)
  MATE_PARA = (/30E6_DP,0.25_DP,0.036_DP/)
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  CALL FIELD_CREATE_START(FIELD_USER_NUMBER,REGION,MATERIAL_FIELD,ERR,ERROR,*999)
  CALL FIELD_TYPE_SET(MATERIAL_FIELD,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
  CALL FIELD_MESH_DECOMPOSITION_SET(MATERIAL_FIELD,DECOMPOSITION,ERR,ERROR,*999)
  CALL FIELD_GEOMETRIC_FIELD_SET(MATERIAL_FIELD,GEOMETRIC_FIELD,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_VARIABLES_SET(MATERIAL_FIELD,NUMBER_OF_FIELD_VARIABLES,ERR,ERROR,*999)
  CALL FIELD_NUMBER_OF_COMPONENTS_SET(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_FIELD_COMPS,ERR,ERROR,*999)
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    DO FIELD_COMP_IDX=1,NUMBER_OF_FIELD_COMPS
      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_COMP_IDX,MESH_COMP_IDX, &
        & ERR,ERROR,*999)
    ENDDO !FIELD_COMP_IDX
  ENDDO !xi
  CALL FIELD_CREATE_FINISH(MATERIAL_FIELD,ERR,ERROR,*999)
  FIELD_DERIVATIVE_IDX = 1
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX = xi
    DO ne=1,NUMBER_OF_ELEMENTS
      DO np=1,MESH_COMP_NUMBER_OF_ELEMENT_NODES(MESH_COMP_IDX) !Can also loop using nu
        DO FIELD_COMP_IDX=1,NUMBER_OF_FIELD_COMPS
          CALL FIELD_PARAMETER_SET_UPDATE_NODE(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,FIELD_DERIVATIVE_IDX, &
            & ELEMENTS(ne)%MESH_COMP(MESH_COMP_IDX)%ELEM_NODE_NUMBERS(np),FIELD_COMP_IDX,MATE_PARA(FIELD_COMP_IDX), &
            & ERR,ERROR,*999)
        ENDDO !FIELD_COMP_IDX
      ENDDO !np
    ENDDO !ne
  ENDDO !xi

  !=CREATE EQUATION SET============================================================================================================
  !Create a Elasticity Class, Linear Elasticity type, no subtype, equations_set
  EQUATION_SET_USER_NUMBER = 1
  NULLIFY(EQUATIONS_SET)
  CALL EQUATIONS_SET_CREATE_START(EQUATION_SET_USER_NUMBER,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
  !Set the equations set to be a Elasticity Class, Linear Elasticity type, no subtype, equations_set
  CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_ELASTICITY_CLASS,EQUATIONS_SET_LINEAR_ELASTICITY_TYPE, &
        EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE,ERR,ERROR,*999)
  CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  !Create the equations set dependent field variables
  FIELD_USER_NUMBER = 2 !Dependent Field
  CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,FIELD_USER_NUMBER,DEPENDENT_FIELD,ERR,ERROR,*999)
  CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  !Create the equations set material field variables
  FIELD_USER_NUMBER = 3 !Material Field
  CALL EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,FIELD_USER_NUMBER,MATERIAL_FIELD,ERR,ERROR,*999)  
  CALL EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  !Create the equations set equations
  NULLIFY(EQUATIONS)
  CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
  CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
                                            !EQUATIONS_SPARSE_MATRICES=1 !<Use sparse matrices for the equations.
                                            !EQUATIONS_FULL_MATRICES=2 !<Use fully populated matrices for the equations. 
  CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_TIMING_OUTPUT,ERR,ERROR,*999)
                                          !EQUATIONS_ELEMENT_MATRIX_OUTPUT=3 !<All below and element matrices output.
                                          !EQUATIONS_MATRIX_OUTPUT=2 !<All below and equation matrices output.
                                          !EQUATIONS_TIMING_OUTPUT=1 !<Timing information output.
                                          !EQUATIONS_NO_OUTPUT=0 !<No output.
  CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999) 

  !=PRESCRIBE BOUNDARY CONDITIONS==================================================================================================
  NULLIFY(BOUNDARY_CONDITIONS)
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
  !Allocate Memory for Displacement & Force BC
  NULLIFY(DISP_BC,FORCE_BC)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DISP_BC",ERR,ERROR,*999)
  ALLOCATE(DISP_BC,STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate FORCE_BC",ERR,ERROR,*999)
  ALLOCATE(FORCE_BC,STAT=ERR)
  ALLOCATE(DISP_BC%MESH_COMP(MESH_COMP_IDX),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DISP_BC%MESH_COMP(MESH_COMP_IDX)",ERR,ERROR,*999)
  ALLOCATE(FORCE_BC%MESH_COMP(MESH_COMP_IDX),STAT=ERR)
  IF(ERR/=0) CALL FLAG_ERROR("Could not allocate FORCE_BC%MESH_COMP(MESH_COMP_IDX)",ERR,ERROR,*999)
  !Initialize Number of Displacement & Force BC
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    DISP_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV=0
    FORCE_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV=0
  ENDDO
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Number of BC & thier DERIV
  DISP_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV(1)=2
  DISP_BC%MESH_COMP(2)%NUMBER_OF_BC_NODES_IN_DERIV=DISP_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV
  FORCE_BC%MESH_COMP(1)%NUMBER_OF_BC_NODES_IN_DERIV(1)=2
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    DO nu=1,NUMBER_OF_DERIV
      NUMBER_OF_BC_NODES = DISP_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        ALLOCATE(DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(NUMBER_OF_BC_NODES))
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER", &
          & ERR,ERROR,*999) !!TODO:: Add Num to string
        ALLOCATE(DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(NUMBER_OF_BC_NODES))
        !Initialize Displacement BC Derivative_Node_Numbers & Node Coordinates
        DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER=0
        DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES=0.0_DP
      ENDIF
      NUMBER_OF_BC_NODES = FORCE_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        ALLOCATE(FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(NUMBER_OF_BC_NODES))
        IF(ERR/=0) CALL FLAG_ERROR("Could not allocate DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER", &
          & ERR,ERROR,*999) !!TODO:: Add Num to string
        ALLOCATE(FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(NUMBER_OF_BC_NODES))
        !Initialize Displacement BC Derivative_Node_Numbers & Node Coordinates
        FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER=0
        FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES=0.0_DP
      ENDIF
    ENDDO !nu
  ENDDO !xi
  !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Prescribe Displacement BC Derivative Node Numbers & Node Coordinates
  DISP_BC%MESH_COMP(1)%DERIV(1)%NODE_NUMBER(:)=(/1,3/)
  DISP_BC%MESH_COMP(2)=DISP_BC%MESH_COMP(1)
  !Prescribe Force BC Derivative Node Numbers & Node Coordinates
  FORCE_BC%MESH_COMP(1)%DERIV(1)%NODE_NUMBER(:)=(/2,4/)
  FORCE_BC%MESH_COMP(1)%DERIV(1)%NODE_COORDINATES(:)=800.0_DP
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  DO xi=1,NUMBER_OF_XI
    MESH_COMP_IDX=xi
    FIELD_COMP_IDX=xi
    DO nu=1,NUMBER_OF_DERIV
      !Displacement BC
      NUMBER_OF_BC_NODES = DISP_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        DO np=1,NUMBER_OF_BC_NODES
          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,nu, &
            & DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(np),FIELD_COMP_IDX,BOUNDARY_CONDITION_FIXED, &
            & DISP_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(np),ERR,ERROR,*999)
        ENDDO !np
      ENDIF
      !Force BC
      NUMBER_OF_BC_NODES = FORCE_BC%MESH_COMP(MESH_COMP_IDX)%NUMBER_OF_BC_NODES_IN_DERIV(nu)
      IF(NUMBER_OF_BC_NODES/=0) THEN
        DO np=1,NUMBER_OF_BC_NODES
          CALL BOUNDARY_CONDITIONS_SET_NODE(BOUNDARY_CONDITIONS,FIELD_DELUDELN_VARIABLE_TYPE,nu, &
            & FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_NUMBER(np),FIELD_COMP_IDX,BOUNDARY_CONDITION_FIXED, &
            & FORCE_BC%MESH_COMP(MESH_COMP_IDX)%DERIV(nu)%NODE_COORDINATES(np),ERR,ERROR,*999)
        ENDDO !np
      ENDIF
    ENDDO !nu
  ENDDO !xi
  CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
  
  !=CREATE PROBLEM=================================================================================================================
  !Create the problem
  NULLIFY(PROBLEM)
  PROBLEM_USER_NUMBER=1
  CALL PROBLEM_CREATE_START(PROBLEM_USER_NUMBER,PROBLEM,ERR,ERROR,*999)
  !Set the problem to be a elasticity class, linear elasticity type with no subtype.
  CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_ELASTICITY_CLASS,PROBLEM_LINEAR_ELASTICITY_TYPE, &
    & PROBLEM_NO_SUBTYPE,ERR,ERROR,*999)
  CALL PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=CREATE PROBLEM CONTROL LOOP====================================================================================================
  !Create the problem control loop
  CALL PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=CREATE PROBLEM SOLVER==========================================================================================================
  !Start the creation of the problem solvers
  NULLIFY(SOLVER)
  SOLVER_IDX = 1
  CALL PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,SOLVER_IDX,SOLVER,ERR,ERROR,*999)
  !CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
                                     !SOLVER_CMISS_LIBRARY=LIBRARY_CMISS_TYPE !<CMISS (internal) solver library 
                                     !SOLVER_PETSC_LIBRARY=LIBRARY_PETSC_TYPE !<PETSc solver library
  !CALL SOLVER_LINEAR_TYPE_SET(SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
                                    !SOLVER_LINEAR_DIRECT_SOLVE_TYPE=1 !<Direct linear solver type \NOT IMPLEMENTED
                                    !SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE=2 !<Iterative linear solver type
  !CALL SOLVER_LINEAR_DIRECT_TYPE_SET(SOLVER,SOLVER_DIRECT_LU,ERR,ERROR,*999)  
                                           !SOLVER_DIRECT_LU=1 !<LU direct linear solver 
                                           !SOLVER_DIRECT_CHOLESKY=2 !<Cholesky direct linear solver \NOT IMPLEMENTED
                                           !SOLVER_DIRECT_SVD=3 !<SVD direct linear solver \NOT IMPLEMENTED
  CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_MATRIX_OUTPUT,ERR,ERROR,*999)
                                    !SOLVER_MATRIX_OUTPUT=4 !<SolVER matrices output from the solver routines plus below
                                    !SOLVER_SOLVER_OUTPUT=3 !<Solver specific output from the solver routines plus below
                                    !SOLVER_TIMING_OUTPUT=2 !<Timing output from the solver routines plus below
                                    !SOLVER_PROGRESS_OUTPUT=1 !<Progress output from solver routines 
                                    !SOLVER_NO_OUTPUT=0 !<No output from the solver routines 
  CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=CREATE PROBLEM SOLVER EQUATIONS================================================================================================
  !Create the problem solver equations
  NULLIFY(SOLVER) !WHY NULLIFY HERE AFTER ALREADY NULLIFIED ABOVE????
  NULLIFY(SOLVER_EQUATIONS)
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
  CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,SOLVER_IDX,SOLVER,ERR,ERROR,*999)
  CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
  CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
                                                          !SOLVER_SPARSE_MATRICES=1 !<Use sparse solver matrices
                                                          !SOLVER_FULL_MATRICES=2 !<Use fully populated solver matrices
  EQUATION_SET_IDX = 0 !Initialize index of the equations set that has been added 
                        !(Variable is returned from PROBLEM_SOLVER_EQUATIONS_EQUATIONS_SET_ADD)
  CALL SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_SET,EQUATION_SET_IDX,ERR,ERROR,*999)
  !CALL SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SOLVER_MATRICES_LINEAR_ONLY,ERR,ERROR,*999)
                                         !SOLVER_MATRICES_ALL=1 !<Select all the solver matrices and vectors 
                                         !SOLVER_MATRICES_DYNAMIC_ONLY=2 !<Select only the dynamic solver matrices and vectors 
                                         !SOLVER_MATRICES_LINEAR_ONLY=3 !<Select only the linear solver matrices and vectors 
                                         !SOLVER_MATRICES_NONLINEAR_ONLY=4 !<Select only the nonlinear solver matrices and vectors
                                         !SOLVER_MATRICES_JACOBIAN_ONLY=5 !<Select only the Jacobian solver matrix
                                         !SOLVER_MATRICES_RESIDUAL_ONLY=6 !<Select only the residual solver vector
                                         !SOLVER_MATRICES_RHS_ONLY=7 !<Select only the RHS solver vector
                                         !SOLVER_MATRICES_RHS_RESIDUAL_ONLY=8 !<Select only the residual and RHS solver vectors
  CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)

  !=SOLVE PROBLEM==================================================================================================================
  !Solve the problem
  CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)

  !=OUTPUT SOLUTION================================================================================================================
  !!TODO:: Output reaction forces in ipnode files
  FILE="LinearElasticityExample"
  METHOD="FORTRAN"
  Export_Field=.TRUE.
  IF(Export_Field) THEN
    CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)  
    CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
  ENDIF
  NUMBER_OF_GLOBAL_DEPENDENT_DOFS=DEPENDENT_FIELD%VARIABLES(1)%NUMBER_OF_GLOBAL_DOFS
  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Displacments",ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,FIELD_DATA,ERR,ERROR,*999) 
  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,NUMBER_OF_GLOBAL_DEPENDENT_DOFS,1,1,FIELD_DATA,'(2x,f20.15)','(2x,f20.15)', &
    & ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,FIELD_DATA, &
    & ERR,ERROR,*999)
  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Reaction Forces",ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,FIELD_DATA,ERR,ERROR,*999) 
  CALL WRITE_STRING_VECTOR(GENERAL_OUTPUT_TYPE,1,1,NUMBER_OF_GLOBAL_DEPENDENT_DOFS,1,1,FIELD_DATA,'(2x,f20.15)','(2x,f20.15)', &
    & ERR,ERROR,*999)
  CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,FIELD_DATA, &
    & ERR,ERROR,*999)
  !================================================================================================================================

  !Calculate the stop times and write out the elapsed user and system times
  CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
  CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)

  CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
    & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)

  CALL CMISS_FINALISE(ERR,ERROR,*999)

  WRITE(*,'(A)') "Program successfully completed."

  STOP
999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
  STOP

END PROGRAM LINEARELASTICITYEXAMPLE
