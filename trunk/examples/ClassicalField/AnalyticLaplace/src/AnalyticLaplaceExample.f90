!> \file
!> $Id: AnalyticLaplaceExample.f90 20 2007-05-28 20:22:52Z cpb $
!> \author Chris Bradley
!> \brief This is an example program to solve an Analytic Laplace equation using openCMISS calls.
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

!> \example ClassicalField/AnalyticLaplace/src/AnalyticLaplaceExample.f90
!! Example illustrating the use of openCMISS to solve the Laplace problem and check with its Analytic Solution.
!! 
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/AnalyticLaplace/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/AnalyticLaplace/build-gnu'>Linux GNU Build</a>
!< 

!> Main program
PROGRAM ANALYTICLAPLACEEXAMPLE

  USE MPI
  USE OPENCMISS
  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: ORIGIN(2)=(/-3.141592653579_CMISSDP/2, -3.141592653579_CMISSDP/2/)
  REAL(CMISSDP), PARAMETER :: HEIGHT=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=2.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: LENGTH=2.0_CMISSDP

  !Program types

  !Program variables

  TYPE(CMISSRegionType) :: WORLD_REGION
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
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
  CALL CMISSInitialise(WorldCoordinateSystem,WORLD_REGION,Err)

  CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(2,6,2)
  CALL ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(2,10,2)
  CALL ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_EXPORT(2,2,0)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    TYPE(CMISSFieldType) :: FIELD

    CALL ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"AnalyticLaplaceBilinear",Err)
    
    CALL ANALYTICLAPLACE_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSDP,VALUE,0.5_CMISSDP,ERR)
    
    WRITE(*,'(A)') "Analytic Laplace Example Testcase1 - bilinear lagrange is successfully completed."
    
  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    CALL ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,3,X_VALUES,Y_VALUES)
    
   CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
   CALL TEST_FRAMEWORK_ASSERT_EQUALS(4.0_CMISSDP,VALUE,1.0_CMISSDP,Err)
   IF (Err/=0) THEN
     WRITE(*,'(A,F3.5)') "Analytic Laplace Example Testcase2 - bicubic Hermite failure: Convergence should be around 4.0" &
       & //", but it was ", VALUE
   ENDIF
   WRITE(*,'(A)') "Analytic Laplace Example Testcase2 - bicubic Hermite is successfully completed."

  END SUBROUTINE ANALYTICLAPLACE_TESTCASE_BICUBIC_HERMITE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of bilinear langrange interpolation is expected.
  SUBROUTINE ANALYTICLAPLACE_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
    & NUMBER_OF_ELEMENTS_XI_INTERVAL,INTERPOLATION_SPECIFICATIONS,X_VALUES,Y_VALUES)
  
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<interpolation specifications
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    !Local Variables
    REAL(CMISSDP) :: VALUE
    
    INTEGER(CMISSIntg) :: i
    TYPE(CMISSFieldType) :: FIELD
    
    ALLOCATE(X_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)
    ALLOCATE(Y_VALUES((NUMBER_OF_ELEMENTS_XI_END-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1),STAT=ERR)

    DO i = NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL
      
      CALL ANALYTICLAPLACE_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL CMISSAnalyticAnalysisAbsoluteErrorGetNode(FIELD,1,1,(i+1)**2/2+1,1,VALUE,Err)

      Y_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(VALUE)
      X_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(HEIGHT/i)
      CALL ANALYTICLAPLACE_GENERIC_CLEAN(1,1,1,1,1)
   
    ENDDO
  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CONVERGENCE
  
  
  !
  !================================================================================================================================
  !   
    
  SUBROUTINE ANALYTICLAPLACE_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_SPECIFICATIONS,DEPENDENT_FIELD)
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements on x axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements on y axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements on z axis
    INTEGER(CMISSIntg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<the interpolation specifications
    TYPE(CMISSFieldType) :: DEPENDENT_FIELD
    !Local Variables
    INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: ANALYTIC_FUNCTION
    INTEGER(CMISSIntg) :: EquationsSetIndex

    TYPE(CMISSBasisType) :: Basis
    TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
    TYPE(CMISSGeneratedMeshType) :: GENERATED_MESH
    TYPE(CMISSMeshType) :: MESH
    TYPE(CMISSDecompositionType) :: DECOMPOSITION
    TYPE(CMISSEquationsType) :: EQUATIONS
    TYPE(CMISSEquationsSetType) :: EQUATIONS_SET
    TYPE(CMISSFieldType) :: ANALYTIC_FIELD,GEOMETRIC_FIELD
    TYPE(CMISSProblemType) :: PROBLEM
    TYPE(CMISSRegionType) :: REGION
    TYPE(CMISSSolverType) :: SOLVER
    TYPE(CMISSSolverEquationsType) :: SOLVER_EQUATIONS
    
    NUMBER_OF_DOMAINS=1

    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !Start the creation of a new RC coordinate system
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(1,CoordinateSystem,Err)
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
    CALL CMISSRegionTypeInitialise(REGION,Err)
    CALL CMISSRegionCreateStart(1,WORLD_REGION,REGION,Err)
    !Set the regions coordinate system to the 2D RC coordinate system that we have created
    CALL CMISSRegionCoordinateSystemSet(REGION,CoordinateSystem,Err)
    !Finish the creation of the region
    CALL CMISSRegionCreateFinish(REGION,Err)

  
    !Start the creation of a basis (default is trilinear lagrange)
    CALL CMISSBasisTypeInitialise(Basis,Err)
    CALL CMISSBasisCreateStart(1,Basis,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      !Set the basis to be a bilinear basis
      CALL CMISSBasisNumberOfXiSet(Basis,2,Err)
      CALL CMISSBasisInterpolationXiSet(Basis,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS/),Err)
    ELSE
      !Set the basis to be a trilinear basis
      CALL CMISSBasisNumberOfXiSet(Basis,3,Err)
      CALL CMISSBasisInterpolationXiSet(Basis,(/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS, &
          & INTERPOLATION_SPECIFICATIONS/),Err)
    ENDIF
    !Finish the creation of the basis
    CALL CMISSBasisCreateFinish(Basis,Err)

    !Start the creation of a generated mesh in the region
    CALL CMISSGeneratedMeshTypeInitialise(GENERATED_MESH,Err)
    CALL CMISSGeneratedMeshCreateStart(1,REGION,GENERATED_MESH,Err)
    !Set up a regular 100x100 mesh
    CALL CMISSGeneratedMeshTypeSet(GENERATED_MESH,1,Err)
    CALL CMISSGeneratedMeshBasisSet(GENERATED_MESH,Basis,Err)
    !Define the mesh on the region
    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
      CALL CMISSGeneratedMeshExtentSet(GENERATED_MESH,(/WIDTH,HEIGHT/),Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS/), &
        & Err)
      CALL CMISSGeneratedMeshOriginSet(GENERATED_MESH,ORIGIN,Err)
    ELSE
      CALL CMISSGeneratedMeshExtentSet(GENERATED_MESH,(/WIDTH,HEIGHT,LENGTH/),Err)
      CALL CMISSGeneratedMeshNumberOfElementsSet(GENERATED_MESH,(/NUMBER_GLOBAL_X_ELEMENTS, &
        & NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/),Err)
    ENDIF
    !Finish the creation of a generated mesh in the region
    CALL CMISSMeshTypeInitialise(MESH,Err)
    CALL CMISSGeneratedMeshCreateFinish(GENERATED_MESH,1,MESH,Err)
    
    !Create a decomposition
    CALL CMISSDecompositionTypeInitialise(DECOMPOSITION,Err)
    CALL CMISSDecompositionCreateStart(1,MESH,DECOMPOSITION,Err)
    !Set the decomposition to be a general decomposition with the specified number of domains
    CALL CMISSDecompositionTypeSet(DECOMPOSITION,CMISSDecompositionCalculatedType,Err)
    CALL CMISSDecompositionNumberOfDomainsSet(DECOMPOSITION,NUMBER_OF_DOMAINS,Err)
    CALL CMISSDecompositionCreateFinish(DECOMPOSITION,Err)

    !Start to create a default (geometric) field on the region
    CALL CMISSFieldTypeInitialise(GEOMETRIC_FIELD,Err)
    CALL CMISSFieldCreateStart(1,REGION,GEOMETRIC_FIELD,Err)
    !Set the decomposition to use
    CALL CMISSFieldMeshDecompositionSet(GEOMETRIC_FIELD,DECOMPOSITION,Err)
    !Set the domain to be used by the field components
    !NB these are needed now as the default mesh component number is 1
    CALL CMISSFieldComponentMeshComponentSet(GEOMETRIC_FIELD,CMISSFieldUVariableType,1,1,Err)
    CALL CMISSFieldComponentMeshComponentSet(GEOMETRIC_FIELD,CMISSFieldUVariableType,2,1,Err)
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      CALL CMISSFieldComponentMeshComponentSet(GEOMETRIC_FIELD,CMISSFieldUVariableType,3,1,Err)
    ENDIF
    !Finish creating the field
    CALL CMISSFieldCreateFinish(GEOMETRIC_FIELD,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMeshGeometricParametersCalculate(GEOMETRIC_FIELD,GENERATED_MESH,Err)

    !Create the equations_set
    CALL CMISSEquationsSetTypeInitialise(EQUATIONS_SET,Err)
    CALL CMISSEquationsSetCreateStart(1,REGION,GEOMETRIC_FIELD,EQUATIONS_SET,Err)
    !Set the equations set to be a standard Laplace problem
    CALL CMISSEquationsSetSpecificationSet(EQUATIONS_SET,CMISSEquationsSetClassicalFieldClass, &
      & CMISSEquationsSetLaplaceEquationType,CMISSEquationsSetStandardLaplaceSubtype,Err)
    !Finish creating the equations set
    CALL CMISSEquationsSetCreateFinish(EQUATIONS_SET,Err)
  
    !Create the equations set analytic field variables
    CALL CMISSFieldTypeInitialise(DEPENDENT_FIELD,Err)
    CALL CMISSEquationsSetDependentCreateStart(EQUATIONS_SET,2,DEPENDENT_FIELD,Err)
    !Finish the equations set dependent field variables
    CALL CMISSEquationsSetDependentCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set analytic field variables
    IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
      ANALYTIC_FUNCTION=CMISSEquationsSetLaplaceEquationThreeDim2
    ELSE
      ANALYTIC_FUNCTION=CMISSEquationsSetLaplaceEquationTwoDim2
    ENDIF
    CALL CMISSFieldTypeInitialise(ANALYTIC_FIELD,Err)
    CALL CMISSEquationsSetAnalyticCreateStart(EQUATIONS_SET,ANALYTIC_FUNCTION,3,ANALYTIC_FIELD,Err)
    !Finish the equations set analtyic field variables
    CALL CMISSEquationsSetAnalyticCreateFinish(EQUATIONS_SET,Err)

    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(EQUATIONS,Err)
    CALL CMISSEquationsSetEquationsCreateStart(EQUATIONS_SET,EQUATIONS,Err)
    !Set the equations matrices sparsity type
    CALL CMISSEquationsSparsityTypeSet(EQUATIONS,CMISSEquationsSparseMatrices,Err)
    CALL CMISSEquationsSetEquationsCreateFinish(EQUATIONS_SET,Err)

    !Set up the boundary conditions as per the analytic solution
    CALL CMISSEquationsSetBoundaryConditionsAnalytic(EQUATIONS_SET,Err)
  
    !Create the problem
    CALL CMISSProblemTypeInitialise(PROBLEM,Err)
    CALL CMISSProblemCreateStart(1,PROBLEM,Err)
    !Set the problem to be a standard Laplace problem
    CALL CMISSProblemSpecificationSet(PROBLEM,CMISSProblemClassicalFieldClass,CMISSProblemLaplaceEquationType, &
      & CMISSProblemStandardLaplaceSubtype,Err)
    !Finish creating the problem
    CALL CMISSProblemCreateFinish(PROBLEM,Err)

    !Create the problem control loop
    CALL CMISSProblemControlLoopCreateStart(PROBLEM,Err)
    !Finish creating the problem control
    CALL CMISSProblemControlLoopCreateFinish(PROBLEM,Err)

    !Start the creation of the problem solvers
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSProblemSolversCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
    !Finish the creation of the problem solver
    CALL CMISSProblemSolversCreateFinish(Problem,Err)

    !Start the creation of the problem solver equations
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSSolverEquationsTypeInitialise(Solver_Equations,Err)
    CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
    !Get the solve equations
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,1,Solver,Err)
    CALL CMISSSolverSolverEquationsGet(Solver,Solver_Equations,Err)
    !Set the solver equations sparsity
    CALL CMISSSolverEquationsSparsityTypeSet(Solver_Equations,CMISSSolverEquationsSparseMatrices,Err)
    !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsFullMatrices,Err)
    !Add in the equations set
    CALL CMISSSolverEquationsEquationsSetAdd(Solver_Equations,Equations_Set,EquationsSetIndex,Err)
    !Finish the creation of the problem solver equations
    CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

    !Solve the problem
    CALL CMISSProblemSolve(PROBLEM,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC

  SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber,GeneratedMeshUserNumber, &
    & ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: GeneratedMeshUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL CMISSProblemDestroy(ProblemUserNumber,Err)
    CALL CMISSGeneratedMeshDestroy(GeneratedMeshUserNumber,Err)
    CALL CMISSBasisDestroy(BasisUserNumber,Err)
    CALL CMISSRegionDestroy(RegionUserNumber,Err)
    CALL CMISSCoordinateSystemDestroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE ANALYTICLAPLACE_GENERIC_CLEAN



END PROGRAM ANALYTICLAPLACEEXAMPLE 
