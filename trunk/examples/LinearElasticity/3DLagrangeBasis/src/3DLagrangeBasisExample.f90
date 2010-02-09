! \file
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
PROGRAM LinearElasticity3DLagrangeBasis

  USE MPI
  USE OPENCMISS

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: ORIGIN(3)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=1

  !Program types

  !Program variables

  TYPE(CMISSRegionType) :: World_Region
  TYPE(CMISSCoordinateSystemType) :: WorldCoordinateSystem
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif

  !Generic CMISS variables
  INTEGER(CMISSIntg) :: Err

  CHARACTER(LEN=255) :: ERROR

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
  CALL CMISSInitialise(WorldCoordinateSystem,World_Region,Err)


  CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,CMISSBasisLinearLagrangeInterpolation,CMISSSolverLinearDirectSolveType)
  CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,CMISSBasisLinearLagrangeInterpolation,CMISSSolverLinearIterativeSolveType)
  !CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,CMISSBasisQuadraticLagrangeInterpolation,CMISSSolverLinearDirectSolveType)
  !CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,CMISSBasisQuadraticLagrangeInterpolation,CMISSSolverLinearIterativeSolveType)
  !CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,CMISSBasisCubicLagrangeInterpolation,CMISSSolverLinearDirectSolveType)
  !CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,CMISSBasisCubicLagrangeInterpolation,CMISSSolverLinearIterativeSolveType)
  !CALL RUN_LINEAR_ELASTICITY_PROBLEM(1,1,1,MixedInterpolation,CMISSSolverLinearDirectSolveType)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE RUN_LINEAR_ELASTICITY_PROBLEM(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS,InterpolationType,SolverType)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<Number of global X elements
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<Number of global Y elements
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<Number of global Z elements
    INTEGER(CMISSIntg), INTENT(IN) :: InterpolationType !<Type of Interpolation
    INTEGER(CMISSIntg), INTENT(IN) :: SolverType !<Type of Solver

    CALL LINEAR_ELASTICITY_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,InterpolationType, & 
      & SolverType)
    
    CALL LINEAR_ELASTICITY_GENERIC_CLEAN(1,1,2,3,1,1)

  END SUBROUTINE RUN_LINEAR_ELASTICITY_PROBLEM

  !
  !================================================================================================================================
  !

  !>Sets the Error string specified by a character string and flags an Error 
  SUBROUTINE FLAG_ERROR(STRING,Err,ERROR,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: STRING !<The Error condition string
    INTEGER(CMISSIntg), INTENT(OUT) :: Err !<The Error code
    CHARACTER(LEN=255), INTENT(OUT) :: ERROR !<The Error string
    !Local variables
    INTEGER(CMISSIntg) :: STRING_LENGTH

    IF(Err==0) Err=1
    STRING_LENGTH=LEN_TRIM(STRING)
    ERROR=STRING(1:STRING_LENGTH)

    RETURN 1
  END SUBROUTINE FLAG_ERROR

  !
  !================================================================================================================================
  !
    
  SUBROUTINE LINEAR_ELASTICITY_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & InterpolationType,SolverType)
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements on x axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements on y axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements on z axis
    INTEGER(CMISSIntg), INTENT(IN) :: InterpolationType !<Type of Interpolation
    INTEGER(CMISSIntg), INTENT(IN) :: SolverType !<Type of Solver
    !Local Variables
    !Program types

    TYPE(CMISSBasisType) :: Basis(3)
    TYPE(CMISSBoundaryConditionsType) :: BoundaryConditions
    TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
    TYPE(CMISSDecompositionType) :: Decomposition
    TYPE(CMISSEquationsType) :: Equations
    TYPE(CMISSEquationsSetType) :: EquationsSet
    TYPE(CMISSFieldType) :: GeometricField,DependentField,MaterialField
    TYPE(CMISSFieldsType) :: Fields
    TYPE(CMISSMeshType) :: Mesh
    TYPE(CMISSNodesType) :: Nodes
    TYPE(CMISSProblemType) :: Problem
    TYPE(CMISSRegionType) :: Region
    TYPE(CMISSSolverType) :: Solver
    TYPE(CMISSSolverEquationsType) :: SolverEquations

    !Element Node Types for prescribing Nodal positions
    !Elements(:) 					#TYPE(ElementType)
    !  MeshComponent(:) 				#TYPE(ElementNodesType)
    !    Element					#TYPE(Mesh_Elements_TYPE)
    !    ElementNodeNumbers(:)		#INTEGER(CMISSIntg),ALLOCATABLE
    !    Derivative(8)					#TYPE(NodeCoordinatesType)
    !      NodeCoordinates(:)	#REAL(CMISSDP),ALLOCATABLE 
    TYPE NodeCoordinatesType
      REAL(CMISSDP),ALLOCATABLE :: NodeCoordinates(:)
    END TYPE NodeCoordinatesType
    TYPE ElementNodesType
      TYPE(CMISSMeshElementsType) :: Element
      INTEGER(CMISSIntg),ALLOCATABLE :: ElementNodeNumbers(:)
      TYPE(NodeCoordinatesType) :: Derivative(8)
    END TYPE ElementNodesType
    TYPE ElementType
      TYPE(ElementNodesType),POINTER :: MeshComponent(:)
    END TYPE ElementType
    TYPE(ElementType), POINTER :: Elements(:)

    !Boundary Condition Types
    !DisplacementBC 								#TYPE(BC_MeshComponent)
    !  MeshComponent(:) 						#TYPE(BC_ElementNodesType)
    !    NumberOfBCNodesForDerivative(8)	#INTEGER(CMISSIntg)
    !    Derivative(8)							#TYPE(BC_NodeCoordinatesType)
    !      NodeNumber(:)					#INTEGER(CMISSIntg),ALLOCATABLE
    !      NodeCoordinates(:)			#REAL(CMISSDP),ALLOCATABLE 
    TYPE BC_NodeCoordinatesType
      INTEGER(CMISSIntg),ALLOCATABLE :: NodeNumber(:)
      REAL(CMISSDP),ALLOCATABLE :: NodeCoordinates(:)
    END TYPE BC_NodeCoordinatesType
    TYPE BC_ElementNodesType
      INTEGER(CMISSIntg) :: NumberOfBCNodesForDerivative(8)
      TYPE(BC_NodeCoordinatesType) :: Derivative(8)
    END TYPE BC_ElementNodesType
    TYPE BC_MeshComponent
      TYPE(BC_ElementNodesType),POINTER :: MeshComponent(:)
    END TYPE BC_MeshComponent
    TYPE(BC_MeshComponent), POINTER :: DisplacementBC,ForceBC

    !Program variables

    INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: BasisUserNumber
    INTEGER(CMISSIntg) :: NumberOfMeshComponents,NumberOfXi,NumberOfNodes,NumberOfElements,NumberOfDerivatives
    INTEGER(CMISSIntg) :: NumberOfFieldVariables,NumberOfBCNodes,NumberOfFieldComponents
    INTEGER(CMISSIntg) :: MeshComponentIdx,FieldComponentIdx,FieldDerivativeIdx,SolverIdx,EquationSetIdx,np,xi,ne,nu
    INTEGER(CMISSIntg) :: Interpolation(3,3),NumberOfGaussPoints(3)
    INTEGER(CMISSIntg),ALLOCATABLE :: MeshComponentNumberOfElementNodes(:)
    
    REAL(CMISSDP) :: l,w,h,MaterialParameters(6)

    LOGICAL :: EXPORT_FIELD

    CHARACTER(LEN=255) :: ERROR

    !=BROADCAST PARAMETERS TO COMPUTATIONAL NODES====================================================================================
    NUMBER_OF_DOMAINS=1
    NumberOfXi = 3
    NumberOfElements = 1
    NumberOfNodes = 8
      
    !Broadcast the number of Elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

    !=CREATE COORDINATE SYSTEM=====================================================================================================
    !Start the creation of a new RC coordinate system
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfXi,Err)
    CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,ORIGIN,Err)
    CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

    !=CREATE REGION==================================================================================================================
    !Create Region and set CS to newly created 3D RC CS
    CALL CMISSRegionTypeInitialise(Region,Err)
    CALL CMISSRegionCreateStart(RegionUserNumber,World_Region,Region,Err)
    CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
    CALL CMISSRegionCreateFinish(Region,Err)

    !=CREATE BASIS ==================================================================================================================
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !Prescribe Interpolation in each direction
    !NOTE if you change interpolation you need to change Boundary Conditions
    Interpolation(1,:) = InterpolationType
    Interpolation(2,:) = InterpolationType
    Interpolation(3,:) = InterpolationType
    !Prescribe Number of Gauss Points
    !NOTE:: Num of Gauss points must be the same across X,Y & Z coordinates and be sufficient to accurately integrate the hightest order interpolation being used
    NumberOfGaussPoints = (/4,4,4/)
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    DO xi=1,NumberOfXi
      BasisUserNumber = xi
      CALL CMISSBasisTypeInitialise(Basis(xi),Err)
      CALL CMISSBasisCreateStart(BasisUserNumber,Basis(xi),Err)
      CALL CMISSBasisTypeSet(Basis(xi),CMISSBasisLagrangeHermiteTPType,Err)
      CALL CMISSBasisNumberOfXiSet(Basis(xi),NumberOfXi,Err)
      CALL CMISSBasisInterpolationXiSet(Basis(xi),Interpolation(xi,1:NumberOfXi),Err)
      CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis(xi),NumberOfGaussPoints(1:NumberOfXi),Err)
      CALL CMISSBasisCreateFinish(Basis(xi),Err)
    ENDDO !xi

    !=MANUAL Mesh CREATION===========================================================================================================

    !=NODE CREATION==================================================================================================================
    !Create nodes in REGION and set initial coordinates to 0,0,0
    CALL CMISSNodesTypeInitialise(NODES,Err)
    CALL CMISSNodesCreateStart(REGION,NumberOfNodes,NODES,Err)
    CALL CMISSNodesCreateFinish(NODES,Err)
    
    !=CREATE MESH====================================================================================================================
    !Create a Mesh with xi number of coordinates
    NumberOfMeshComponents = NumberOfXi
    CALL CMISSMeshTypeInitialise(Mesh,Err)
    CALL CMISSMeshCreateStart(MeshUserNumber,REGION,NumberOfXi,Mesh,Err)
    CALL CMISSMeshNumberOfElementsSet(Mesh,NumberOfElements,Err)
    CALL CMISSMeshNumberOfComponentsSet(Mesh,NumberOfMeshComponents,Err)

    ALLOCATE(Elements(NumberOfElements),STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate Elements",Err,ERROR,*999)
    DO ne=1,NumberOfElements
      ALLOCATE(Elements(ne)%MeshComponent(NumberOfMeshComponents),STAT=Err)
      IF(Err/=0) CALL FLAG_ERROR("Could not allocate Element(ne)%MeshComponent",Err,ERROR,*999)
    ENDDO !ne
    ALLOCATE(MeshComponentNumberOfElementNodes(NumberOfMeshComponents),STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate MeshComponentNumberOfElementNodes",Err,ERROR,*999)
    DO xi=1,NumberOfXi
      MeshComponentIdx=xi
      IF (Interpolation(1,1)<=3) THEN !Linear/quadratic/cubic Lagrange Interpolation
        MeshComponentNumberOfElementNodes(MeshComponentIdx)=PRODUCT(Interpolation(xi,:)+1)
        NumberOfDerivatives = 1
      ELSEIF (Interpolation(1,1)==4) THEN !CubicHermite Interpolation
        MeshComponentNumberOfElementNodes(MeshComponentIdx)=(2**NumberOfXi)
        NumberOfDerivatives = 2**NumberOfXi
      ELSE
        CALL FLAG_ERROR("This example is only setup for Lagrange and cubic hermite bases",Err,ERROR,*999)
      ENDIF
    ENDDO !xi
    DO ne=1,NumberOfElements
      DO xi=1,NumberOfXi
        MeshComponentIdx=xi
        ALLOCATE(Elements(ne)%MeshComponent(MeshComponentIdx)%ElementNodeNumbers(MeshComponentNumberOfElementNodes( &
                & MeshComponentIdx)),STAT=Err)
        IF(Err/=0) CALL FLAG_ERROR("Could not allocate Element(ne)%MeshComponent(MeshComponentIdx)%ElementNodeNumbers", &
            & Err,ERROR,*999)
        Elements(ne)%MeshComponent(MeshComponentIdx)%ElementNodeNumbers = 0
      ENDDO !xi
    ENDDO !ne
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !Prescribe Element node numbers
    SELECT CASE(InterpolationType)
    CASE(CMISSBasisLinearLagrangeInterpolation)
      Elements(1)%MeshComponent(1)%ElementNodeNumbers(:) = (/1,2,3,4,5,6,7,8/)
      Elements(1)%MeshComponent(2)%ElementNodeNumbers(:) = (/1,2,3,4,5,6,7,8/)
      Elements(1)%MeshComponent(3)%ElementNodeNumbers(:) = (/1,2,3,4,5,6,7,8/)
    CASE(CMISSBasisQuadraticLagrangeInterpolation)
      DO np=1,27
        Elements(1)%MeshComponent(1)%ElementNodeNumbers(np) = np
        Elements(1)%MeshComponent(2)%ElementNodeNumbers(np) = np
        Elements(1)%MeshComponent(3)%ElementNodeNumbers(np) = np
      ENDDO
    CASE(CMISSBasisCubicLagrangeInterpolation)
      DO np=1,64
        Elements(1)%MeshComponent(1)%ElementNodeNumbers(np) = np
        Elements(1)%MeshComponent(2)%ElementNodeNumbers(np) = np
        Elements(1)%MeshComponent(3)%ElementNodeNumbers(np) = np
      ENDDO
    END SELECT
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !!TODO:: When diagnostics are on - Mesh_TOPOLOGY_Elements_Element_NODES_SET trys to output Mesh%TOPOLOGY(1)%PTR%Elements%Elements(1)%Mesh_Element_NODES which are only allocated when the Mesh_TOPOLOGY_Elements_CREATE_FINISH command is given
    DO xi=1,NumberOfXi
      MeshComponentIdx=xi
      DO ne=1,NumberOfElements
        CALL CMISSMeshElementsCreateStart(Mesh,MeshComponentIdx,Basis(xi),Elements(ne)%MeshComponent(MeshComponentIdx)%Element,Err)
        CALL CMISSMeshElementsNodesSet(Elements(ne)%MeshComponent(MeshComponentIdx)%Element,ne, & 
          & Elements(ne)%MeshComponent(MeshComponentIdx)%ElementNodeNumbers(:),Err)
        CALL CMISSMeshElementsCreateFinish(Elements(ne)%MeshComponent(MeshComponentIdx)%Element,Err)
      ENDDO !ne
    ENDDO !xi
    CALL CMISSMeshCreateFinish(Mesh,Err)

    !=CREATE Decomposition===========================================================================================================
    !Create Mesh Decomposition dividing Mesh into number_of_domains for parallel solving
    CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
    CALL CMISSDecompositionCreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
    CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
    CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    CALL CMISSDecompositionCreateFinish(Decomposition,Err)

    !=CREATE GEOMETRIC FIELD=========================================================================================================
    !Start to create a default (geometric) field on the region
    NumberOfFieldVariables = 1 !Geometric Field Coordinates
    NumberOfFieldComponents = NumberOfXi
    CALL CMISSFieldTypeInitialise(GeometricField,Err)
    CALL CMISSFieldCreateStart(GeometricFieldUserNumber,REGION,GeometricField,Err)
    CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
    CALL CMISSFieldNumberOfVariablesSet(GeometricField,NumberOfFieldVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,NumberOfFieldComponents,Err)
    DO xi=1,NumberOfXi
      MeshComponentIdx = xi
      FieldComponentIdx = xi
      CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,FieldComponentIdx,MeshComponentIdx,Err)
    ENDDO !xi
    CALL CMISSFieldCreateFinish(GeometricField,Err)
    !Initialize node coordinates
    DO ne=1,NumberOfElements
      DO xi=1,NumberOfXi
        MeshComponentIdx=xi
        DO nu=1,NumberOfDerivatives
          ALLOCATE(Elements(ne)%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates&
            &(MeshComponentNumberOfElementNodes(MeshComponentIdx)),STAT=Err)
          IF(Err/=0) CALL FLAG_ERROR("Could not allocate Element(ne)%MeshComponent(MeshComponentIdx)%Derivative(nu)%&
              &NodeCoordinates",Err,ERROR,*999)
          Elements(ne)%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates = 0.0_CMISSDP
        ENDDO
      ENDDO
    ENDDO
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !Prescribe node coordinates
    l = 120.0_CMISSDP
    w = 160.0_CMISSDP
    h = 10.0_CMISSDP
    SELECT CASE(InterpolationType)
    CASE(CMISSBasisLinearLagrangeInterpolation)
      Elements(1)%MeshComponent(1)%Derivative(1)%NodeCoordinates = (/0.0_CMISSDP,l,0.0_CMISSDP,l,0.0_CMISSDP,l,0.0_CMISSDP,l/)
      Elements(1)%MeshComponent(2)%Derivative(1)%NodeCoordinates = (/0.0_CMISSDP,0.0_CMISSDP,w,w,0.0_CMISSDP,0.0_CMISSDP,w,w/)
      Elements(1)%MeshComponent(3)%Derivative(1)%NodeCoordinates = (/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,h,h,h,h/)
    CASE(CMISSBasisQuadraticLagrangeInterpolation)
      Elements(1)%MeshComponent(1)%Derivative(1)%NodeCoordinates = (/0.0_CMISSDP,l/2.0_CMISSDP,l,0.0_CMISSDP,l,0.0_CMISSDP,l, &
          & 0.0_CMISSDP,l/)
      Elements(1)%MeshComponent(2)%Derivative(1)%NodeCoordinates = (/0.0_CMISSDP,0.0_CMISSDP,w,w,0.0_CMISSDP,0.0_CMISSDP,w,w/)
      Elements(1)%MeshComponent(3)%Derivative(1)%NodeCoordinates = (/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,h,h,h,h/)
    CASE(CMISSBasisCubicLagrangeInterpolation)

    END SELECT
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=
    DO xi=1,NumberOfXi
      MeshComponentIdx = xi
      FieldComponentIdx = xi
      DO ne=1,NumberOfElements
        DO nu=1,NumberOfDerivatives
          DO np=1,MeshComponentNumberOfElementNodes(MeshComponentIdx) !Can also loop using nu
             CALL CMISSFieldParameterSetUpdateNode(GeometricField,CMISSFieldUVariableType,CMISSFieldValuesSetType,nu, &
              & Elements(ne)%MeshComponent(MeshComponentIdx)%ElementNodeNumbers(np),FieldComponentIdx, &
              & Elements(ne)%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates(np),Err)
          ENDDO !np
        ENDDO !nu
      ENDDO !ne
    ENDDO !xi

    !=CREATE DEPENDENT FIELD=========================================================================================================
    !Create a dependent field
    NumberOfFieldVariables = 2 !Dependent Field Displacement Coordinates & Force 
    NumberOfFieldComponents = NumberOfXi
    CALL CMISSFieldTypeInitialise(DependentField,Err)
    CALL CMISSFieldCreateStart(DependentFieldUserNumber,REGION,DependentField,Err)
    CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
    CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err)
    CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
    CALL CMISSFieldNumberOfVariablesSet(DependentField,NumberOfFieldVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,NumberOfFieldComponents,Err)
    DO xi=1,NumberOfXi
      MeshComponentIdx = xi
      FieldComponentIdx = xi
      CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,FieldComponentIdx,MeshComponentIdx,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,FieldComponentIdx,MeshComponentIdx,Err)
    ENDDO !xi
    CALL CMISSFieldCreateFinish(DependentField,Err)

    !=CREATE MATERIAL FIELD==========================================================================================================
    !!TODO:: Set Material Field Interpolation to constant Element based interpolation when field i/o and cmgui allows for this
    NumberOfFieldVariables = 1
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    SELECT CASE(NumberOfXi)
    CASE(1)
      !Prescribe material properties Area,E1
      NumberOfFieldComponents = 2 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/w*h,10.0E3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    CASE(2)
      !Prescribe material properties h,E1,v12
      NumberOfFieldComponents = 3 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/h,10.0E3_CMISSDP,0.3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    CASE(3)
      !Prescribe material properties E1,E2,E3 & v13,v23,v12
      NumberOfFieldComponents = 6 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/10.0E3_CMISSDP,10.0E3_CMISSDP,10.0E3_CMISSDP,0.3_CMISSDP,0.3_CMISSDP,0.3_CMISSDP/)
    END SELECT
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    CALL CMISSFieldTypeInitialise(MaterialField,Err)
    CALL CMISSFieldCreateStart(MaterialFieldUserNumber,REGION,MaterialField,Err)
    CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
    CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
    CALL CMISSFieldNumberOfVariablesSet(MaterialField,NumberOfFieldVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,NumberOfFieldComponents,Err)
    DO xi=1,NumberOfXi
      MeshComponentIdx = xi
      DO FieldComponentIdx=1,NumberOfFieldComponents
        CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,FieldComponentIdx,MeshComponentIdx,Err)
      ENDDO !FieldComponentIdx
    ENDDO !xi
    CALL CMISSFieldCreateFinish(MaterialField,Err)
    FieldDerivativeIdx = 1
    DO xi=1,NumberOfXi
      MeshComponentIdx = xi
      DO ne=1,NumberOfElements
        DO np=1,MeshComponentNumberOfElementNodes(MeshComponentIdx) !Can also loop using nu
          DO FieldComponentIdx=1,NumberOfFieldComponents
            CALL CMISSFieldParameterSetUpdateNode(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType, &
              & FieldDerivativeIdx,Elements(ne)%MeshComponent(MeshComponentIdx)%ElementNodeNumbers(np),FieldComponentIdx, &
              & MaterialParameters(FieldComponentIdx),Err)
          ENDDO !FieldComponentIdx
        ENDDO !np
      ENDDO !ne
    ENDDO !xi

    !=CREATE EQUATIONS SET===========================================================================================================
    !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,REGION,GeometricField,EquationsSet,Err)
    !Set the equations set to be a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass,CMISSEquationsSetLinearElasticityType, &
      & CMISSEquationsSetThreeDimensionalSubtype,Err)
    CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)
    !Create the equations set dependent field variables
    CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)
    !Create the equations set material field variables
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialFieldUserNumber,MaterialField,Err)
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)

    !=CREATE EQUATIONS SET EQUATION==================================================================================================
    !Create the equations set equations
    CALL CMISSEquationsTypeInitialise(Equations,Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
    CALL CMISSEquationsSparsityTypeSet(EQUATIONS,CMISSEquationsFullMatrices,Err)
                                                !CMISSEquationsSparseMatrices=1 !<Use sparse matrices for the equations.
                                                !CMISSEquationsFullMatrices=2 !<Use fully populated matrices for the equations. 
    CALL CMISSEquationsOutputTypeSet(EQUATIONS,CMISSEquationsTimingOutput,Err)
                                              !CMISSEquationsNoOutput !<No output from the equations.
                                              !CMISSEquationsTimingOutput !<Timing information output.
                                              !CMISSEquationsMatrixOutput !<All below and equation matrices output.
                                              !CMISSEquationsElementMatrixOutput !<All below and Element matrices output.
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

    !=PRESCRIBE BOUNDARY CONDITIONS==================================================================================================
    CALL CMISSBoundaryConditionsTypeInitialise(BoundaryConditions,Err)
    CALL CMISSEquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
    !Allocate Memory for Displacement & Force BC
    NULLIFY(DisplacementBC,ForceBC)
    ALLOCATE(DisplacementBC,STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate DisplacementBC",Err,ERROR,*999)
    ALLOCATE(ForceBC,STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate ForceBC",Err,ERROR,*999)
    ALLOCATE(DisplacementBC%MeshComponent(MeshComponentIdx),STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate DisplacementBC%MeshComponent",Err,ERROR,*999)
    ALLOCATE(ForceBC%MeshComponent(MeshComponentIdx),STAT=Err)
    IF(Err/=0) CALL FLAG_ERROR("Could not allocate ForceBC%MeshComponent",Err,ERROR,*999)
    !Initialize Number of Displacement & Force BC
    DO xi=1,NumberOfXi
      MeshComponentIdx=xi
      DisplacementBC%MeshComponent(MeshComponentIdx)%NumberOfBCNodesForDerivative=0
      ForceBC%MeshComponent(MeshComponentIdx)%NumberOfBCNodesForDerivative=0
    ENDDO
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !Prescribe Number of BC & thier Derivatives
    DisplacementBC%MeshComponent(1)%NumberOfBCNodesForDerivative(1)=4
    DisplacementBC%MeshComponent(2)%NumberOfBCNodesForDerivative=DisplacementBC%MeshComponent(1)%NumberOfBCNodesForDerivative
    DisplacementBC%MeshComponent(3)%NumberOfBCNodesForDerivative=DisplacementBC%MeshComponent(1)%NumberOfBCNodesForDerivative
    ForceBC%MeshComponent(1)%NumberOfBCNodesForDerivative(1)=4
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    DO xi=1,NumberOfXi
      MeshComponentIdx=xi
      DO nu=1,NumberOfDerivatives
        !Displacement BC
        NumberOfBCNodes = DisplacementBC%MeshComponent(MeshComponentIdx)%NumberOfBCNodesForDerivative(nu)
        IF(NumberOfBCNodes/=0) THEN
          ALLOCATE(DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber(NumberOfBCNodes),STAT=Err)
          IF(Err/=0) CALL FLAG_ERROR("Could not allocate DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber",&
              & Err,ERROR,*999)
          ALLOCATE(DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates(NumberOfBCNodes),STAT=Err)
          IF(Err/=0) CALL FLAG_ERROR("Could not allocate DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%&
              &NodeCoordinates",Err,ERROR,*999)
          !Initialize Displacement BC Derivativesative_NodeNumbers & Node Coordinates
          DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber=0
          DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates=0.0_CMISSDP
        ENDIF
        !Force BC
        NumberOfBCNodes = ForceBC%MeshComponent(MeshComponentIdx)%NumberOfBCNodesForDerivative(nu)
        IF(NumberOfBCNodes/=0) THEN
          ALLOCATE(ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber(NumberOfBCNodes),STAT=Err)
          IF(Err/=0) CALL FLAG_ERROR("Could not allocate ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber", &
              & Err,ERROR,*999)
          ALLOCATE(ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates(NumberOfBCNodes),STAT=Err)
          IF(Err/=0) CALL FLAG_ERROR("Could not allocate ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates", &
            & Err,ERROR,*999)
          !Initialize Displacement BC Derivativesative_NodeNumbers & Node Coordinates
          ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber=0
          ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates=0.0_CMISSDP
        ENDIF
      ENDDO !nu
    ENDDO !xi
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    !Prescribe Displacement BC Derivativesative Node Numbers & Node Coordinates
    DisplacementBC%MeshComponent(1)%Derivative(1)%NodeNumber(:)=(/1,3,5,7/)
    DisplacementBC%MeshComponent(2)=DisplacementBC%MeshComponent(1)
    DisplacementBC%MeshComponent(3)=DisplacementBC%MeshComponent(1)
    !Prescribe Force BC Derivativesative Node Numbers & Node Coordinates
    ForceBC%MeshComponent(1)%Derivative(1)%NodeNumber(:)=(/2,4,6,8/)
    ForceBC%MeshComponent(1)%Derivative(1)%NodeCoordinates(:)=400.0_CMISSDP
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    DO xi=1,NumberOfXi
      MeshComponentIdx=xi
      FieldComponentIdx=xi
      DO nu=1,NumberOfDerivatives
        !Displacement BC
        NumberOfBCNodes = DisplacementBC%MeshComponent(MeshComponentIdx)%NumberOfBCNodesForDerivative(nu)
        IF(NumberOfBCNodes/=0) THEN
          DO np=1,NumberOfBCNodes
            CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldUVariableType,nu, &
              & DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber(np),FieldComponentIdx, &
              & CMISSBoundaryConditionFixed,DisplacementBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates(np),Err)
          ENDDO !np
        ENDIF
        !Force BC
        NumberOfBCNodes = ForceBC%MeshComponent(MeshComponentIdx)%NumberOfBCNodesForDerivative(nu)
        IF(NumberOfBCNodes/=0) THEN
          DO np=1,NumberOfBCNodes
            CALL CMISSBoundaryConditionsSetNode(BoundaryConditions,CMISSFieldDelUDelNVariableType,nu, &
              & ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeNumber(np),FieldComponentIdx, &
              & CMISSBoundaryConditionFixed,ForceBC%MeshComponent(MeshComponentIdx)%Derivative(nu)%NodeCoordinates(np),Err)
          ENDDO !np
        ENDIF
      ENDDO !nu
    ENDDO !xi
    CALL CMISSEquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)
    
    !=CREATE Problem=================================================================================================================
    !Create the Problem
    CALL CMISSProblemTypeInitialise(Problem,Err)
    CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
    !Set the Problem to be a elasticity class, linear elasticity type with no subtype.
    CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemLinearElasticityType, &
      & CMISSProblemNoSubtype,Err)
    CALL CMISSProblemCreateFinish(Problem,Err)

    !=CREATE Problem CONTROL LOOP====================================================================================================
    !Create the Problem control loop
     CALL CMISSProblemControlLoopCreateStart(Problem,Err)
     CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

    !=CREATE Problem SOLVER==========================================================================================================
    !Start the creation of the Problem solvers
    SolverIdx = 1
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSProblemSolversCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverIdx,SOLVER,Err)
    CALL CMISSSolverLibraryTypeSet(SOLVER,CMISSSolverPETScLibrary,Err)
                                         !CMISSSolverCMISSLibrary     !<CMISS (internal) solver library.
                                         !CMISSSolverPETScLibrary     !<PETSc solver library.
                                         !CMISSSolverMUMPSLibrary     !<MUMPS solver library.
                                         !CMISSSolverSuperLULibrary   !<SuperLU solver library.
                                         !CMISSSolverSpoolesLULibrary !<SPOOLES solver library.
                                         !CMISSSolverUMFPACKLibrary   !<UMFPACK solver library.
                                         !CMISSSolverLUSOLLibrary     !<LUSOL solver library.
                                         !CMISSSolverESSLLibrary      !<ESSL solver library.
                                         !CMISSSolverLAPACKLibrary    !<LAPACK solver library.
                                         !CMISSSolverTAOLibrary       !<TAO solver library.
                                         !CMISSSolverHypreLibrary     !<Hypre solver library.
    CALL CMISSSolverLinearTypeSet(SOLVER,CMISSSolverLinearDirectSolveType,Err)
                                        !CMISSSolverLinearDirectSolveType    !<Direct linear solver type.
                                        !CMISSSolverLinearIterativeSolveType !<Iterative linear solver type.
    !CALL CMISSSolverLinearDirectTypeSet(SOLVER,CMISSSolverDirectLU,Err)
                                              !CMISSSolverDirectLU       !<LU direct linear solver.
                                              !CMISSSolverDirectCholesky !<Cholesky direct linear solver.
                                              !CMISSSolverDirectSVD      !<SVD direct linear solver.
    CALL CMISSSolverOutputTypeSet(SOLVER,CMISSSolverSolverMatrixOutput,Err)
                                        !CMISSSolverNoOutput !<No output from the solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                        !CMISSSolverProgressOutput !<Progress output from solver routines.
                                        !CMISSSolverTimingOutput !<Timing output from the solver routines plus below.
                                        !CMISSSolverSolverOutput !<Solver specific output from the solver routines plus below.
                                        !CMISSSolverSolverMatrixOutput !<Solver matrices output from the solver routines plus below.
    CALL CMISSProblemSolversCreateFinish(Problem,Err)

    !=CREATE Problem SOLVER EQUATIONS================================================================================================
    !Create the Problem solver equations
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
    CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverIdx,SOLVER,Err)
    CALL CMISSSolverSolverEquationsGet(SOLVER,SolverEquations,Err)
    CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
                                                            !CMISSSolverEquationsSparseMatrices !<Use sparse solver matrices.
                                                            !CMISSSolverEquationsFullMatrices !<Use fully populated solver matrices.
    EquationSetIdx = 1 !Initialize index of the equations set that has been added 
                          !(Variable is returned from Problem_SolverEquations_EquationsSet_ADD)
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationSetIdx,Err)
    CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

    !=SOLVE Problem==================================================================================================================
    !Solve the Problem
    CALL CMISSProblemSolve(Problem,Err)

    !=OUTPUT SOLUTION================================================================================================================
    !!TODO:: Output reaction forces in ipnode files
    EXPORT_FIELD=.TRUE.
    IF(EXPORT_FIELD) THEN
      CALL CMISSFieldsTypeInitialise(Fields,Err)
      CALL CMISSFieldsTypeCreate(Region,Fields,Err)
      CALL CMISSFieldIONodesExport(Fields,"LinearElasticity3DLinearLagrangeBasisExample","FORTRAN",Err)
      CALL CMISSFieldIOElementsExport(Fields,"LinearElasticity3DLinearLagrangeBasisExample","FORTRAN",Err)
      CALL CMISSFieldsTypeFinalise(Fields,Err)
    ENDIF

    RETURN
999 WRITE(*,'(A)') ERROR
    STOP 1

  END SUBROUTINE LINEAR_ELASTICITY_GENERIC

  !
  !================================================================================================================================
  !  

  SUBROUTINE LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber1, &
    & BasisUserNumber2,BasisUserNumber3,ProblemUserNumber)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: CoordinateSystemUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: RegionUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber1
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber2
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber3
    INTEGER(CMISSIntg), INTENT(IN) :: ProblemUserNumber

    CALL CMISSProblemDestroy(ProblemUserNumber,Err)
    CALL CMISSBasisDestroy(BasisUserNumber1,Err)
    CALL CMISSBasisDestroy(BasisUserNumber2,Err)
    CALL CMISSBasisDestroy(BasisUserNumber3,Err)
    CALL CMISSRegionDestroy(RegionUserNumber,Err)
    CALL CMISSCoordinateSystemDestroy(CoordinateSystemUserNumber,Err)

  END SUBROUTINE LINEAR_ELASTICITY_GENERIC_CLEAN

END PROGRAM LinearElasticity3DLagrangeBasis


