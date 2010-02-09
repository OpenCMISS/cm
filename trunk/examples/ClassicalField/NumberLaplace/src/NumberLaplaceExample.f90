!> \file
!> $Id: NumberLaplaceExample.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Chris Bradley
!> \brief This is an example program to solve a Laplace equation using OpenCMISS calls with objects accessed by user number.
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

!> \example ClassicalField/NumberLaplace/src/NumberLaplaceExample.f90
!! Example program to solve a Laplace equation using openCMISS calls with objects accessed by user number.
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Laplace/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/Laplace/build-gnu'>Linux GNU Build</a>
!!
!<

!> Main program
PROGRAM NUMBERLAPLACEEXAMPLE

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
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=10
 
  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
  INTEGER(CMISSIntg) :: MPI_IERROR
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: WorldCoordinateSystemUserNumber
  INTEGER(CMISSIntg) :: WorldRegionUserNumber
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

  !Intialise cmiss
  CALL CMISSInitialise(WorldCoordinateSystemUserNumber,WorldRegionUserNumber,Err)
 
  NUMBER_GLOBAL_X_ELEMENTS=2
  NUMBER_GLOBAL_Y_ELEMENTS=2
  NUMBER_GLOBAL_Z_ELEMENTS=0
  NUMBER_OF_DOMAINS=1
    
  !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the coordinate system to be 2D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystemUserNumber,2,Err)
  ELSE
    !Set the coordinate system to be 3D
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystemUserNumber,3,Err)
  ENDIF
  !Finish the creation of the coordinate system
  CALL CMISSCoordinateSystemCreateFinish(CoordinateSystemUserNumber,Err)

  !Start the creation of the region
  CALL CMISSRegionCreateStart(RegionUserNumber,WorldRegionUserNumber,Err)
  !Set the regions coordinate system to the 2D RC coordinate system that we have created
  CALL CMISSRegionCoordinateSystemSet(RegionUserNumber,CoordinateSystemUserNumber,Err)
  !Finish the creation of the region
  CALL CMISSRegionCreateFinish(RegionUserNumber,Err)

  !Start the creation of a basis (default is trilinear lagrange)
  CALL CMISSBasisCreateStart(BasisUserNumber,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    !Set the basis to be a bilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(BasisUserNumber,2,Err)
  ELSE
    !Set the basis to be a trilinear Lagrange basis
    CALL CMISSBasisNumberOfXiSet(BasisUserNumber,3,Err)
  ENDIF
  !Finish the creation of the basis
  CALL CMISSBasisCreateFinish(BasisUserNumber,Err)
    
  !Start the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshCreateStart(GeneratedMeshUserNumber,RegionUserNumber,Err)
  !Set up a regular x*y*z mesh
  CALL CMISSGeneratedMeshTypeSet(GeneratedMeshUserNumber,CMISSGeneratedMeshRegularMeshType,Err)
  !Set the default basis
  CALL CMISSGeneratedMeshBasisSet(GeneratedMeshUserNumber,BasisUserNumber,Err)
  !Define the mesh on the region
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    CALL CMISSGeneratedMeshExtentSet(GeneratedMeshUserNumber,(/WIDTH,HEIGHT/),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMeshUserNumber,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/),Err)
  ELSE
    CALL CMISSGeneratedMeshExtentSet(GeneratedMeshUserNumber,(/WIDTH,HEIGHT,LENGTH/),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMeshUserNumber,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
      & NUMBER_GLOBAL_Z_ELEMENTS/),Err)
  ENDIF
  !Finish the creation of a generated mesh in the region
  CALL CMISSGeneratedMeshCreateFinish(GeneratedMeshUserNumber,MeshUserNumber,Err)
  
  !Create a decomposition
  CALL CMISSDecompositionCreateStart(DecompositionUserNumber,RegionUserNumber,MeshUserNumber,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMISSDecompositionTypeSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,CMISSDecompositionCalculatedType,Err)
  CALL CMISSDecompositionNumberOfDomainsSet(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL CMISSDecompositionCreateFinish(RegionUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  
  !Start to create a default (geometric) field on the region
  CALL CMISSFieldCreateStart(GeometricFieldUserNumber,RegionUserNumber,Err)
  !Set the decomposition to use
  CALL CMISSFieldMeshDecompositionSet(RegionUserNumber,GeometricFieldUserNumber,MeshUserNumber,DecompositionUserNumber,Err)
  !Set the domain to be used by the field components.
  CALL CMISSFieldComponentMeshComponentSet(RegionUserNumber,GeometricFieldUserNumber,CMISSFieldUVariableType,1,1,Err)
  CALL CMISSFieldComponentMeshComponentSet(RegionUserNumber,GeometricFieldUserNumber,CMISSFieldUVariableType,2,1,Err)
  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
    CALL CMISSFieldComponentMeshComponentSet(RegionUserNumber,GeometricFieldUserNumber,CMISSFieldUVariableType,3,1,Err)
  ENDIF
  !Finish creating the field
  CALL CMISSFieldCreateFinish(RegionUserNumber,GeometricFieldUserNumber,Err)
       
  !Update the geometric field parameters
  CALL CMISSGeneratedMeshGeometricParametersCalculate(RegionUserNumber,GeometricFieldUserNumber,GeneratedMeshUserNumber,Err)
  
  !Create the equations_set
  CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,RegionUserNumber,GeometricFieldUserNumber,Err)
  !Set the equations set to be a standard Laplace problem
  CALL CMISSEquationsSetSpecificationSet(RegionUserNumber,EquationsSetUserNumber,CMISSEquationsSetClassicalFieldClass, &
    & CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetStandardLaplaceSubtype,Err)
  !Finish creating the equations set
  CALL CMISSEquationsSetCreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)

  !Create the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateStart(RegionUserNumber,EquationsSetUserNumber,DependentFieldUserNumber,Err)
  !Finish the equations set dependent field variables
  CALL CMISSEquationsSetDependentCreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)

  !Create the equations set equations
  CALL CMISSEquationsSetEquationsCreateStart(RegionUserNumber,EquationsSetUserNumber,Err)
  !Set the equations matrices sparsity type
  CALL CMISSEquationsSparsityTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISSEquationsSparseMatrices,Err)
  !Set the equations set output
  !CALL CMISSEquationsOutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISSEquationsNoOutput,Err)
  CALL CMISSEquationsOutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(RegionUserNumber,EquationsSetUserNumber,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMISSEquationsSetEquationsCreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)

  !Start the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateStart(RegionUserNumber,EquationsSetUserNumber,Err)
  !Set the first node to 0.0 and the last node to 1.0
  FirstNodeNumber=1
  IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
  ELSE
    LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  ENDIF
  CALL CMISSBoundaryConditionsSetNode(RegionUserNumber,EquationsSetUserNumber,CMISSFieldUVariableType,1,FirstNodeNumber,1, &
    & CMISSBoundaryConditionFixed,0.0_CMISSDP,Err)
  CALL CMISSBoundaryConditionsSetNode(RegionUserNumber,EquationsSetUserNumber,CMISSFieldUVariableType,1,LastNodeNumber,1, &
    & CMISSBoundaryConditionFixed,1.0_CMISSDP,Err)
  !Finish the creation of the equations set boundary conditions
  CALL CMISSEquationsSetBoundaryConditionsCreateFinish(RegionUserNumber,EquationsSetUserNumber,Err)
  
  !Start the creation of a problem.
  CALL CMISSProblemCreateStart(ProblemUserNumber,Err)
  !Set the problem to be a standard Laplace problem
  CALL CMISSProblemSpecificationSet(ProblemUserNumber,CMISSProblemClassicalFieldClass,CMISSProblemLaplaceEquationType, &
    & CMISSProblemStandardLaplaceSubtype,Err)
  !Finish the creation of a problem.
  CALL CMISSProblemCreateFinish(ProblemUserNumber,Err)

  !Start the creation of the problem control loop
  CALL CMISSProblemControlLoopCreateStart(ProblemUserNumber,Err)
  !Finish creating the problem control loop
  CALL CMISSProblemControlLoopCreateFinish(ProblemUserNumber,Err)
 
  !Start the creation of the problem solvers
  CALL CMISSProblemSolversCreateStart(ProblemUserNumber,Err)
  !CALL CMISSSolverOutputTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverNoOutput,Err)
  !CALL CMISSSolverOutputTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverProgressOutput,Err)
  !CALL CMISSSolverOutputTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverSolverOutput,Err)
  CALL CMISSSolverOutputTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverSolverMatrixOutput,Err)
  !Finish the creation of the problem solver
  CALL CMISSProblemSolversCreateFinish(ProblemUserNumber,Err)

  !Start the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateStart(ProblemUserNumber,Err)
  !Set the solver equations sparsity
  CALL CMISSSolverEquationsSparsityTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverEquationsSparseMatrices,Err)
  !CALL CMISSSolverEquationsSparsityTypeSet(ProblemUserNumber,CMISSControlLoopNode,1,CMISSSolverEquationsFullMatrices,Err)  
  !Add in the equations set
  CALL CMISSSolverEquationsEquationsSetAdd(ProblemUserNumber,CMISSControlLoopNode,1,RegionUserNumber,EquationsSetUserNumber, &
    & EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL CMISSProblemSolverEquationsCreateFinish(ProblemUserNumber,Err)

  !Solve the problem
  CALL CMISSProblemSolve(ProblemUserNumber,Err)

  !Finialise CMISS
  CALL CMISSFinalise(Err)
 
  WRITE(*,'(A)') "Program successfully completed."
  
  STOP
  
END PROGRAM NUMBERLAPLACEEXAMPLE
