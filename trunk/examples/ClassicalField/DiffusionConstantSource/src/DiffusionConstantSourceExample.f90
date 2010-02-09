!> \file
!> $Id: DiffusionConstantSourceExample.f90 20 2009-10-09 18:49:52Z anc $
!> \author Chris Bradley
!> \brief This is an example program to solve a diffusion equation using openCMISS calls.
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

!> \example ClassicalField/Diffusion/src/DiffusionExample.f90
!! Example program to solve a diffusion equation using openCMISS calls.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Diffusion/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Diffusion/build-gnu'>Linux GNU Build</a>
!<

!> Main program
PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE



  USE OPENCMISS
  USE MPI


#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: HEIGHT=1.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=3.0_CMISSDP
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=12
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR
  
    !CMISS variables

  TYPE(CMISSBasisType) :: Basis
  TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
  TYPE(CMISSCoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMISSDecompositionType) :: Decomposition
  TYPE(CMISSEquationsType) :: Equations
  TYPE(CMISSEquationsSetType) :: EquationsSet
  TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialsField,SourceField
  TYPE(CMISSFieldsType) :: Fields
  TYPE(CMISSGeneratedMeshType) :: GeneratedMesh  
  TYPE(CMISSMeshType) :: Mesh
  TYPE(CMISSProblemType) :: Problem
  TYPE(CMISSControlLoopType) :: ControlLoop
  TYPE(CMISSRegionType) :: Region,WorldRegion
  TYPE(CMISSSolverType) :: Solver, LinearSolver
  TYPE(CMISSSolverEquationsType) :: SolverEquations

  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
 

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: EquationsSetIndex
  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
  INTEGER(CMISSIntg) :: Err
  
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

  !Intialise OpenCMISS
  CALL CMISSInitialise(WorldCoordinateSystem,WorldRegion,Err)

  NUMBER_GLOBAL_X_ELEMENTS=10
  NUMBER_GLOBAL_Y_ELEMENTS=20
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=1


  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the coordinate system to be 2D
      CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,2,Err)
    ELSE
      !Set the coordinate system to be 3D
      CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,3,Err)
    ENDIF
    !Finish the creation of the coordinate system
    CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)


    !Start the creation of the region
    CALL CMISSRegionTypeInitialise(Region,Err)
    CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegion,Region,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL CMISSRegionCreateFinish(Region,Err)

    !Start the creation of a basis (default is trilinear lagrange)
    CALL CMISSBasisTypeInitialise(Basis,Err)
    CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear Lagrange basis
      CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
    ELSE
      !Set the basis to be a trilinear Lagrange basis
      CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(BASIS,Err)

    !Start the creation of a generated mesh in the region
    CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
    !Set up a regular x*y*z mesh
    CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,CMISSGeneratedMeshRegularMeshType,Err)
    !Set the default basis
    CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)   
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT/),Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
    ELSE
      CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,(/WIDTH,HEIGHT,LENGTH/),Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
        & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
    ENDIF    
    !Finish the creation of a generated mesh in the region
    CALL CMISSMeshTypeInitialise(Mesh,Err)
    CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

    !Create a decomposition
    CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
    CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
    CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    !Finish the decomposition
    CALL CMISSDecompositionCreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMISSFieldTypeInitialise(GeometricField,Err)
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(GeometricField,Err)

       
    !Update the geometric field parameters
    CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)
!  ENDIF

  !IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR
  
  !Create the equations_set
  CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,EquationsSet,Err)
  !Set the equations set to be a standard Laplace problem
  CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetDiffusionEquationType,CMISSEquationsSetConstantSourceDiffusionSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL CMISSFieldTypeInitialise(DependentField,Err)
  CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  CALL CMISSFieldTypeInitialise(MaterialsField,Err)
  CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

  !Create the equations set source field variables
  CALL CMISSFieldTypeInitialise(SourceField,Err)
  CALL CMISSEquationsSetSourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetSourceCreateFinish(EquationsSet,Err)

  !Create the equations set equations
  CALL CMISSEquationsTypeInitialise(Equations,Err)
  CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsNoOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)
  
  !Create the equations set boundary conditions
  CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,1,FirstNodeNumber,1, &
    & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldDeludelnVariableType,1,LastNodeNumber,1, &
    & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)


  !Create the problem
  CALL CMISSProblemTypeInitialise(Problem,Err)
  CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
  !Set the problem to be a No Source Diffusion problem
  CALL CMISSProblemSpecificationSet(Problem,CMISSProblemClassicalFieldClass,CMISSProblemDiffusionEquationType, &
    & CMISSProblemLinearSourceDiffusionSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(Problem,Err)


  !Create the problem control
  CALL CMISSProblemControlLoopCreateStart(Problem,Err)
  !Get the control loop
  !CALL CMISSProblemControlLoopGet(Problem,ControlLoopNode,ControlLoop,Err)
  !Set the times
  !CALL CMISSControlLoopTimesSet(ControlLoop,0.0_CMISSDP,1.0_CMISSDP,0.1_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(Problem,Err)


  !Start the creation of the problem solvers
!  
! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
! 

  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverTypeInitialise(LinearSolver,Err)
  CALL CMISSProblemSolversCreateStart(Problem,Err)
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverProgressOutput,Err)
  CALL CMISSSolverDynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL CMISSSolverLinearIterativeMaximumIterationsSet(LinearSolver,300,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(Problem,Err)


  !Create the problem solver equations
  CALL CMISSSolverTypeInitialise(Solver,Err)
  CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
  CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
  !Get the solve equations
  CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
  CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

  !Solve the problem
  CALL CMISSProblemSolve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL CMISSFieldsTypeInitialise(Fields,Err)
    CALL CMISSFieldsTypeCreate(Region,Fields,Err)
    CALL CMISSFieldIONodesExport(Fields,"DiffusionConstantSource","FORTRAN",Err)
    CALL CMISSFieldIOElementsExport(Fields,"DiffusionConstantSource","FORTRAN",Err)
    CALL CMISSFieldsTypeFinalise(Fields,Err)

  ENDIF

  !CALL CMISSFinalise(Err)
  WRITE(*,'(A)') "Program successfully completed."
  

  STOP


!   USE BASE_ROUTINES
!   USE BASIS_ROUTINES
!   USE BOUNDARY_CONDITIONS_ROUTINES
!   USE CMISS
!   USE CMISS_MPI
!   !USE CMISS_PETSC
!   USE COMP_ENVIRONMENT
!   USE CONSTANTS
!   USE CONTROL_LOOP_ROUTINES
!   USE COORDINATE_ROUTINES
!   USE DISTRIBUTED_MATRIX_VECTOR
!   USE DOMAIN_MAPPINGS
!   USE EQUATIONS_ROUTINES
!   USE EQUATIONS_SET_CONSTANTS
!   USE EQUATIONS_SET_ROUTINES
!   USE FIELD_ROUTINES
!   USE FIELD_IO_ROUTINES
!   USE GENERATED_MESH_ROUTINES
!   USE INPUT_OUTPUT
!   USE ISO_VARYING_STRING
!   USE KINDS
!   USE LISTS
!   USE MESH_ROUTINES
!   USE MPI
!   USE PROBLEM_CONSTANTS
!   USE PROBLEM_ROUTINES
!   USE REGION_ROUTINES
!   USE SOLVER_ROUTINES
!   USE TIMER
!   USE TYPES
! 
! #ifdef WIN32
!   USE IFQWIN
! #endif
! 
!   IMPLICIT NONE
! 
!   !Test program parameters
! 
!   REAL(DP), PARAMETER :: HEIGHT=1.0_DP
!   REAL(DP), PARAMETER :: WIDTH=2.0_DP
!   REAL(DP), PARAMETER :: LENGTH=3.0_DP
!   
!   !Program types
!   
!   !Program variables
! 
!   INTEGER(INTG) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
!   INTEGER(INTG) :: NUMBER_OF_DOMAINS
!   
!   INTEGER(INTG) :: NUMBER_COMPUTATIONAL_NODES
!   INTEGER(INTG) :: MY_COMPUTATIONAL_NODE_NUMBER
!   INTEGER(INTG) :: MPI_IERROR
! 
!   INTEGER(INTG) :: first_global_dof,first_local_dof,first_local_rank,last_global_dof,last_local_dof,last_local_rank,rank_idx
!   INTEGER(INTG) :: EQUATIONS_SET_INDEX
!   TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
!   
!   TYPE(BASIS_TYPE), POINTER :: BASIS
!   TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
!   TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: COORDINATE_SYSTEM
!   TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
!   TYPE(GENERATED_MESH_TYPE), POINTER :: GENERATED_MESH
!   TYPE(MESH_TYPE), POINTER :: MESH
!   TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
!   TYPE(EQUATIONS_TYPE), POINTER :: EQUATIONS
!   TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
!   TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD,DEPENDENT_FIELD,MATERIALS_FIELD,SOURCE_FIELD
!   TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
!   TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
!   TYPE(REGION_TYPE), POINTER :: REGION,WORLD_REGION
!   TYPE(SOLVER_TYPE), POINTER :: SOLVER,LINEAR_SOLVER
!   TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
!   
!   LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
!   TYPE(VARYING_STRING) :: FILE,METHOD
! 
!   REAL(SP) :: START_USER_TIME(1),STOP_USER_TIME(1),START_SYSTEM_TIME(1),STOP_SYSTEM_TIME(1)
! 
! #ifdef WIN32
!   !Quickwin type
!   LOGICAL :: QUICKWIN_STATUS=.FALSE.
!   TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
! #endif
!   
!   !Generic CMISS variables
!   
!   INTEGER(INTG) :: ERR
!   TYPE(VARYING_STRING) :: ERROR
! 
!   INTEGER(INTG) :: DIAG_LEVEL_LIST(5)
!   CHARACTER(LEN=MAXSTRLEN) :: DIAG_ROUTINE_LIST(1),TIMING_ROUTINE_LIST(1)
!   
! #ifdef WIN32
!   !Initialise QuickWin
!   QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
!   QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
!   QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
!   !Set the window parameters
!   QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
!   !If attempt fails set with system estimated values
!   IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
! #endif
! 
!   !Intialise cmiss
!   NULLIFY(WORLD_REGION)
!   CALL CMISS_INITIALISE(WORLD_REGION,ERR,ERROR,*999)
!   
!   !Set all diganostic levels on for testing
!   DIAG_LEVEL_LIST(1)=1
!   DIAG_LEVEL_LIST(2)=2
!   DIAG_LEVEL_LIST(3)=3
!   DIAG_LEVEL_LIST(4)=4
!   DIAG_LEVEL_LIST(5)=5
!   DIAG_ROUTINE_LIST(1)="SOLUTION_MAPPING_CALCULATE"
!   !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"DiffusionExample",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
!   !CALL DIAGNOSTICS_SET_ON(ALL_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
!   !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
!   !CALL DIAGNOSTICS_SET_ON(IN_DIAG_TYPE,DIAG_LEVEL_LIST,"",DIAG_ROUTINE_LIST,ERR,ERROR,*999)
! 
!   TIMING_ROUTINE_LIST(1)="PROBLEM_FINITE_ELEMENT_CALCULATE"
!   !CALL TIMING_SET_ON(IN_TIMING_TYPE,.TRUE.,"",TIMING_ROUTINE_LIST,ERR,ERROR,*999)
!   
!   !Calculate the start times
!   CALL CPU_TIMER(USER_CPU,START_USER_TIME,ERR,ERROR,*999)
!   CALL CPU_TIMER(SYSTEM_CPU,START_SYSTEM_TIME,ERR,ERROR,*999)
!   
!   !Get the number of computational nodes
!   NUMBER_COMPUTATIONAL_NODES=COMPUTATIONAL_NODES_NUMBER_GET(ERR,ERROR)
!   IF(ERR/=0) GOTO 999
!   !Get my computational node number
!   MY_COMPUTATIONAL_NODE_NUMBER=COMPUTATIONAL_NODE_NUMBER_GET(ERR,ERROR)
!   IF(ERR/=0) GOTO 999
!   
! 
!   NUMBER_GLOBAL_X_ELEMENTS=2
!   NUMBER_GLOBAL_Y_ELEMENTS=2
!   NUMBER_GLOBAL_Z_ELEMENTS=0
!   NUMBER_OF_DOMAINS=2
!   
!   !Read in the number of elements in the X & Y directions, and the number of partitions on the master node (number 0)
!   IF(MY_COMPUTATIONAL_NODE_NUMBER==0) THEN
!     WRITE(*,'("Enter the number of elements in the X direction :")')
!     READ(*,*) NUMBER_GLOBAL_X_ELEMENTS
!     WRITE(*,'("Enter the number of elements in the Y direction :")')
!     READ(*,*) NUMBER_GLOBAL_Y_ELEMENTS
!     WRITE(*,'("Enter the number of elements in the Z direction :")')
!     READ(*,*) NUMBER_GLOBAL_Z_ELEMENTS
!     WRITE(*,'("Enter the number of domains :")')
!     READ(*,*) NUMBER_OF_DOMAINS
!   ENDIF
!   !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
!   CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
!   CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
!   CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
!   CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
!   CALL MPI_ERROR_CHECK("MPI_BCAST",MPI_IERROR,ERR,ERROR,*999)
!   CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"COMPUTATIONAL ENVIRONMENT:",ERR,ERROR,*999)
!   CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  Total number of computational nodes = ",NUMBER_COMPUTATIONAL_NODES, &
!     & ERR,ERROR,*999)
!   CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"  My computational node number = ",MY_COMPUTATIONAL_NODE_NUMBER,ERR,ERROR,*999)
! 
!   !Start the creation of a new RC coordinate system
!   NULLIFY(COORDINATE_SYSTEM)
!   CALL COORDINATE_SYSTEM_CREATE_START(1,COORDINATE_SYSTEM,ERR,ERROR,*999)
!   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!     !Set the coordinate system to be 2D
!     CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,2,ERR,ERROR,*999)
!   ELSE
!     !Set the coordinate system to be 3D
!     CALL COORDINATE_SYSTEM_DIMENSION_SET(COORDINATE_SYSTEM,3,ERR,ERROR,*999)
!   ENDIF
!   !Finish the creation of the coordinate system
!   CALL COORDINATE_SYSTEM_CREATE_FINISH(COORDINATE_SYSTEM,ERR,ERROR,*999)
!   
!   !Start the creation of the region
!   NULLIFY(REGION)
!   CALL REGION_CREATE_START(1,WORLD_REGION,REGION,ERR,ERROR,*999)
!   !Set the regions coordinate system to the 2D RC coordinate system that we have created
!   CALL REGION_COORDINATE_SYSTEM_SET(REGION,COORDINATE_SYSTEM,ERR,ERROR,*999)
!   !Finish the creation of the region
!   CALL REGION_CREATE_FINISH(REGION,ERR,ERROR,*999)
!   
!   FILE="DiffusionConstantSourceExample"
!   METHOD="FORTRAN"
!   IMPORT_FIELD=.FALSE.
!   IF(IMPORT_FIELD) THEN
!      CALL FIELD_IO_FIELDS_IMPORT(FILE, METHOD, REGION, MESH, 1, DECOMPOSITION, 1,  &
!       & DECOMPOSITION_CALCULATED_TYPE, FIELD_U_VARIABLE_TYPE, FIELD_ARITHMETIC_MEAN_SCALING, ERR, ERROR, *999)
!   ELSE
!     !Start the creation of a basis (default is trilinear lagrange)
!     NULLIFY(BASIS)
!     CALL BASIS_CREATE_START(1,BASIS,ERR,ERROR,*999)  
!     IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!       !Set the basis to be a bilinear Lagrange basis
!       CALL BASIS_NUMBER_OF_XI_SET(BASIS,2,ERR,ERROR,*999)
!     ELSE
!       !Set the basis to be a trilinear Lagrange basis
!       CALL BASIS_NUMBER_OF_XI_SET(BASIS,3,ERR,ERROR,*999)
!     ENDIF
!     !Finish the creation of the basis
!     CALL BASIS_CREATE_FINISH(BASIS,ERR,ERROR,*999)
!     
!     !Start the creation of a generated mesh in the region
!     NULLIFY(GENERATED_MESH)
!     NULLIFY(MESH)
!     CALL GENERATED_MESH_CREATE_START(1,REGION,GENERATED_MESH,ERR,ERROR,*999)
!     !Set up a regular 100x100 mesh
!     CALL GENERATED_MESH_TYPE_SET(GENERATED_MESH,1, &
!         & ERR,ERROR,*999)
!     CALL GENERATED_MESH_BASIS_SET(GENERATED_MESH,BASIS,ERR,ERROR,*999)
!     
!     !Define the mesh on the region
!     IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
!       CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/WIDTH,HEIGHT/),ERR,ERROR,*999)
!       CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/), &
!         & ERR,ERROR,*999)
!     ELSE
!       CALL GENERATED_MESH_EXTENT_SET(GENERATED_MESH,(/WIDTH,HEIGHT,LENGTH/),ERR,ERROR,*999)
!       CALL GENERATED_MESH_NUMBER_OF_ELEMENTS_SET(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS, &
!         & NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/), ERR,ERROR,*999)
!     ENDIF  
!     
!     !Finish the creation of a generated mesh in the region
!     NULLIFY(MESH)
!     CALL GENERATED_MESH_CREATE_FINISH(GENERATED_MESH,1,MESH,ERR,ERROR,*999) 
! 
!     !Create a decomposition
!     NULLIFY(DECOMPOSITION)
!     CALL DECOMPOSITION_CREATE_START(1,MESH,DECOMPOSITION,ERR,ERROR,*999)
!     !Set the decomposition to be a general decomposition with the specified number of domains
!     CALL DECOMPOSITION_TYPE_SET(DECOMPOSITION,DECOMPOSITION_CALCULATED_TYPE,ERR,ERROR,*999)
!     CALL DECOMPOSITION_NUMBER_OF_DOMAINS_SET(DECOMPOSITION,NUMBER_OF_DOMAINS,ERR,ERROR,*999)
!     CALL DECOMPOSITION_CREATE_FINISH(DECOMPOSITION,ERR,ERROR,*999)
! 
! !!Start to create a default (geometric) field on the region
!     NULLIFY(GEOMETRIC_FIELD)
!     CALL FIELD_CREATE_START(1,REGION,GEOMETRIC_FIELD,ERR,ERROR,*999)
!     !Set the decomposition to use
!     CALL FIELD_MESH_DECOMPOSITION_SET(GEOMETRIC_FIELD,DECOMPOSITION,ERR,ERROR,*999)
!     !Set the domain to be used by the field components
!     !NB these are needed now as the default mesh component number is 1
!     CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1,1,ERR,ERROR,*999)
!     CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,2,1,ERR,ERROR,*999)
!     IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
!       CALL FIELD_COMPONENT_MESH_COMPONENT_SET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,3,1,ERR,ERROR,*999)
!     ENDIF
!     !Finish creating the field
!     CALL FIELD_CREATE_FINISH(GEOMETRIC_FIELD,ERR,ERROR,*999)
! 
!        
!     !Update the geometric field parameters
!     CALL GENERATED_MESH_GEOMETRIC_PARAMETERS_CALCULATE(GEOMETRIC_FIELD,GENERATED_MESH,ERR,ERROR,*999)
! 
!   ENDIF
! 
!   IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR
!   
!   !Create the equations_set
!   NULLIFY(EQUATIONS_SET)
!   CALL EQUATIONS_SET_CREATE_START(1,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,ERR,ERROR,*999)
!   !Set the equations set to be a standard Laplace problem
!   CALL EQUATIONS_SET_SPECIFICATION_SET(EQUATIONS_SET,EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_DIFFUSION_EQUATION_TYPE, &
!     & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE,ERR,ERROR,*999)
!   !Finish creating the equations set
!   CALL EQUATIONS_SET_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! 
!   !Create the equations set dependent field variables
!   NULLIFY(DEPENDENT_FIELD)
!   CALL EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,2,DEPENDENT_FIELD,ERR,ERROR,*999)
!   !Finish the equations set dependent field variables
!   CALL EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! 
!   !Create the equations set material field variables
!   NULLIFY(MATERIALS_FIELD)
!   CALL EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,3,MATERIALS_FIELD,ERR,ERROR,*999)
!   !Finish the equations set material field variables
!   CALL EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! 
!   !Create the equations set source field variables
!   NULLIFY(SOURCE_FIELD)
!   CALL EQUATIONS_SET_SOURCE_CREATE_START(EQUATIONS_SET,4,SOURCE_FIELD,ERR,ERROR,*999)
!   !Finish the equations set source field variables
!   CALL EQUATIONS_SET_SOURCE_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
! 
!   !Create the equations set equations
!   NULLIFY(EQUATIONS)
!   CALL EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
!   !Set matrix lumping
!   CALL EQUATIONS_LUMPING_TYPE_SET(EQUATIONS,EQUATIONS_UNLUMPED_MATRICES,ERR,ERROR,*999)
!   !CALL EQUATIONS_LUMPING_TYPE_SET(EQUATIONS,EQUATIONS_LUMPED_MATRICES,ERR,ERROR,*999)
!   !Set the equations matrices sparsity type
!   CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_SPARSE_MATRICES,ERR,ERROR,*999)
!   !CALL EQUATIONS_SPARSITY_TYPE_SET(EQUATIONS,EQUATIONS_FULL_MATRICES,ERR,ERROR,*999)
!   !Set the equations set output
!   !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_TIMING_OUTPUT,ERR,ERROR,*999)
!   !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_MATRIX_OUTPUT,ERR,ERROR,*999)
!   !CALL EQUATIONS_OUTPUT_TYPE_SET(EQUATIONS,EQUATIONS_ELEMENT_MATRIX_OUTPUT,ERR,ERROR,*999)
!   CALL EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999) 
! 
!   !Create the equations set boundary conditions
!   !Find the first and last dof numbers and ranks
!   NULLIFY(FIELD_VARIABLE)
!   CALL FIELD_VARIABLE_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
!   DEPENDENT_DOF_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
!   first_global_dof=1
!   first_local_dof=0
!   first_local_rank=0
!   DO rank_idx=1,DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%NUMBER_OF_DOMAINS
!     IF(DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
!       first_local_dof=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%LOCAL_NUMBER(rank_idx)
!       first_local_rank=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(first_global_dof)%DOMAIN_NUMBER(rank_idx)
!       EXIT
!     ENDIF
!   ENDDO !rank_idx  
!   NULLIFY(FIELD_VARIABLE)
!   CALL FIELD_VARIABLE_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VARIABLE,ERR,ERROR,*999)
!   DEPENDENT_DOF_MAPPING=>FIELD_VARIABLE%DOMAIN_MAPPING
!   last_global_dof=DEPENDENT_DOF_MAPPING%NUMBER_OF_GLOBAL
!   last_local_dof=0
!   last_local_rank=0
!   DO rank_idx=1,DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%NUMBER_OF_DOMAINS
!     IF(DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%LOCAL_TYPE(rank_idx)/=DOMAIN_LOCAL_GHOST) THEN
!       last_local_dof=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%LOCAL_NUMBER(rank_idx)
!       last_local_rank=DEPENDENT_DOF_MAPPING%GLOBAL_TO_LOCAL_MAP(last_global_dof)%DOMAIN_NUMBER(rank_idx)
!       EXIT
!     ENDIF
!   ENDDO !rank_idx
!   NULLIFY(BOUNDARY_CONDITIONS)
!   CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_START(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
!   IF(MY_COMPUTATIONAL_NODE_NUMBER==first_local_rank) &
!     & CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD_U_VARIABLE_TYPE,first_local_dof, &
!     & BOUNDARY_CONDITION_FIXED,1.0_DP,ERR,ERROR,*999)
!   IF(MY_COMPUTATIONAL_NODE_NUMBER==last_local_rank) &
!     & CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,FIELD_DELUDELN_VARIABLE_TYPE,last_local_dof, &
!     & BOUNDARY_CONDITION_FIXED,1.0_DP,ERR,ERROR,*999)
!   CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_CREATE_FINISH(EQUATIONS_SET,ERR,ERROR,*999)
!   
!   !Create the problem
!   NULLIFY(PROBLEM)
!   CALL PROBLEM_CREATE_START(1,PROBLEM,ERR,ERROR,*999)
!   !Set the problem to be a standard Laplace problem
!   CALL PROBLEM_SPECIFICATION_SET(PROBLEM,PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_DIFFUSION_EQUATION_TYPE, &
!     & PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE,ERR,ERROR,*999)
!   !Finish creating the problem
!   CALL PROBLEM_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! 
!   !Create the problem control
!   NULLIFY(CONTROL_LOOP)
!   CALL PROBLEM_CONTROL_LOOP_CREATE_START(PROBLEM,ERR,ERROR,*999)
!   !Get the control loop
!   CALL PROBLEM_CONTROL_LOOP_GET(PROBLEM,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
!   !Set the times
!   CALL CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,0.0_DP,1.0_DP,0.1_DP,ERR,ERROR,*999)
!   !Finish creating the problem control
!   CALL PROBLEM_CONTROL_LOOP_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! 
!   !Start the creation of the problem solvers
!   NULLIFY(SOLVER)
!   NULLIFY(LINEAR_SOLVER)
!   CALL PROBLEM_SOLVERS_CREATE_START(PROBLEM,ERR,ERROR,*999)
!   CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
!   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_NO_OUTPUT,ERR,ERROR,*999)
!   CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_PROGRESS_OUTPUT,ERR,ERROR,*999)
!   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_TIMING_OUTPUT,ERR,ERROR,*999)
!   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_SOLVER_OUTPUT,ERR,ERROR,*999)
!   !CALL SOLVER_OUTPUT_TYPE_SET(SOLVER,SOLVER_MATRIX_OUTPUT,ERR,ERROR,*999)
!   !CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_EULER_SCHEME,ERR,ERROR,*999)
!   !CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,ERR,ERROR,*999)
!   !CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE,ERR,ERROR,*999)
!   !CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME,ERR,ERROR,*999)
!   !Get the associated linear solver
!   CALL SOLVER_DYNAMIC_LINEAR_SOLVER_GET(SOLVER,LINEAR_SOLVER,ERR,ERROR,*999)
!   CALL SOLVER_LINEAR_ITERATIVE_MAXIMUM_ITERATIONS_SET(LINEAR_SOLVER,300,ERR,ERROR,*999)
!   !Finish the creation of the problem solvers
!   CALL PROBLEM_SOLVERS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! 
!   !Create the problem solver equations
!   NULLIFY(SOLVER)
!   NULLIFY(SOLVER_EQUATIONS)
!   CALL PROBLEM_SOLVER_EQUATIONS_CREATE_START(PROBLEM,ERR,ERROR,*999)
!   CALL PROBLEM_SOLVER_GET(PROBLEM,CONTROL_LOOP_NODE,1,SOLVER,ERR,ERROR,*999)
!   CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
!   CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
!   !CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_FULL_MATRICES,ERR,ERROR,*999)
!   !Add in the equations set
!   CALL SOLVER_EQUATIONS_EQUATIONS_SET_ADD(SOLVER_EQUATIONS,EQUATIONS_SET,EQUATIONS_SET_INDEX,ERR,ERROR,*999)
!   !Finish the problem solver equations
!   CALL PROBLEM_SOLVER_EQUATIONS_CREATE_FINISH(PROBLEM,ERR,ERROR,*999)
! 
!   !Turn of PETSc error handling
!   !CALL PETSC_ERRORHANDLING_SET_OFF(ERR,ERROR,*999)
!   
! 
!   !Solve the problem
!   CALL PROBLEM_SOLVE(PROBLEM,ERR,ERROR,*999)
! 
!   EXPORT_FIELD=.TRUE.
!   IF(EXPORT_FIELD) THEN
!     CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)  
!     CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, ERR,ERROR,*999)
!   ENDIF
! 
!   !Output timing summary
!   !CALL TIMING_SUMMARY_OUTPUT(ERR,ERROR,*999)
! 
!   !Calculate the stop times and write out the elapsed user and system times
!   CALL CPU_TIMER(USER_CPU,STOP_USER_TIME,ERR,ERROR,*999)
!   CALL CPU_TIMER(SYSTEM_CPU,STOP_SYSTEM_TIME,ERR,ERROR,*999)
! 
!   CALL WRITE_STRING_TWO_VALUE(GENERAL_OUTPUT_TYPE,"User time = ",STOP_USER_TIME(1)-START_USER_TIME(1),", System time = ", &
!     & STOP_SYSTEM_TIME(1)-START_SYSTEM_TIME(1),ERR,ERROR,*999)
!   
!   !CALL CMISS_FINALISE(ERR,ERROR,*999)
!     WRITE(*,'("problem solver finished")')  
!   WRITE(*,'(A)') "Program successfully completed."
!   
!   STOP
! 999 CALL CMISS_WRITE_ERROR(ERR,ERROR)
!   STOP
  
END PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE
