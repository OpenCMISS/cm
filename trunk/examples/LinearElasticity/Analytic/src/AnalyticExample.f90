!> \file
!> $Id: LinearElasticityExample.f90 20 2009-02-15 13:26:52Z cpb $
!> \author Chris Bradley
!> \brief This example illustrates the use of openCMISS to solve Linear Elasticity Problems in all dimensions and check with thier Analytic Solutions.
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

!> example LinearElasticity/Analytic/src/AnalyticExample.f90
!! Example illustrating the use of openCMISS to solve Linear Elasticity Problems in all dimensions and check with thier Analytic Solutions.
!! 
!! \par Latest Builds:
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/ANALYTIC_LINEAR_ELASTICITY/build-intel'>Linux Intel Build</a>
!! \li <a href='http://autotest.bioeng.auckland.ac.nz/opencmiss-build/logs_x86_64-linux/ClassicalField/ANALYTIC_LINEAR_ELASTICITY/build-gnu'>Linux GNU Build</a>
!< 

!> Main program
PROGRAM ANALYTIC_LINEAR_ELASTICITY_EXAMPLE

  USE MPI
  USE OPENCMISS
  USE TEST_FRAMEWORK_ROUTINES

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  !Test program parameters

  REAL(CMISSDP), PARAMETER :: ORIGIN(3)=(/0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
  REAL(CMISSDP), PARAMETER :: LENGTH=20.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: WIDTH=20.0_CMISSDP
  REAL(CMISSDP), PARAMETER :: HEIGHT=5.0_CMISSDP

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: MaterialFieldUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=4
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

  CALL ANALYTIC_TESTCASE_LINEAR_LAGRANGE_EXPORT(10,0,0)
  !CALL ANALYTIC_TESTCASE_CUBIC_HERMITE_EXPORT(1,1,1)

  CALL CMISSFinalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP

CONTAINS

  !
  !================================================================================================================================
  !  
    !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_TESTCASE_LINEAR_LAGRANGE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<Number of global X elements
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<Number of global Y elements
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<Number of global Z elements
    !Local Variables
    TYPE(CMISSFieldType) :: FIELD

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,1, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"ANALYTIC_LINEAR_ELASTICITYBilinear",Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTIC_TESTCASE_LINEAR_LAGRANGE_EXPORT
  
  !
  !================================================================================================================================
  !  
    !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_TESTCASE_CUBIC_HERMITE_EXPORT(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
    & NUMBER_GLOBAL_Z_ELEMENTS)

    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<increment interval number of elements per axis
    !Local Variables
    TYPE(CMISSFieldType) :: FIELD

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS,3, &
      & FIELD)

    CALL CMISSAnalyticAnalysisOutput(FIELD,"ANALYTIC_LINEAR_ELASTICITYBilinear",Err)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(1,1,1,1,1)

  END SUBROUTINE ANALYTIC_TESTCASE_CUBIC_HERMITE_EXPORT

  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)
    
    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,1,X_VALUES,Y_VALUES)
    
    CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)

    CALL TEST_FRAMEWORK_ASSERT_EQUALS(2.0_CMISSDP,VALUE,0.5_CMISSDP,ERR)
    
    WRITE(*,'(A)') "Analytic linear elasticity Example Testcase1 - bilinear lagrange is successfully completed."
    
  END SUBROUTINE ANALYTIC_TESTCASE_BILINEAR_LAGRANGE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_TESTCASE_BICUBIC_HERMITE_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START, &
    & NUMBER_OF_ELEMENTS_XI_END,NUMBER_OF_ELEMENTS_XI_INTERVAL)
  
    !Argument variables
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_START !<initial number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_END !<final number of elements per axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_OF_ELEMENTS_XI_INTERVAL !<increment interval number of elements per axis
    !Local Variables
    REAL(CMISSDP) :: VALUE
    REAL(CMISSDP), ALLOCATABLE :: X_VALUES(:),Y_VALUES(:)

    CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
      & NUMBER_OF_ELEMENTS_XI_INTERVAL,3,X_VALUES,Y_VALUES)
    
   CALL TEST_FRAMEWORK_GRADIENT_VALUE_GET(X_VALUES,Y_VALUES,VALUE)
   CALL TEST_FRAMEWORK_ASSERT_EQUALS(4.0_CMISSDP,VALUE,1.0_CMISSDP,Err)
   IF (Err/=0) THEN
     WRITE(*,'(A,F3.5)') "Analytic linear elasticity Example Testcase2 - bicubic Hermite failure: Convergence should be around " &
       & //"4.0, but it was ", VALUE
   ENDIF
   WRITE(*,'(A)') "Analytic linear elasticity Example Testcase2 - bicubic Hermite is successfully completed."

  END SUBROUTINE ANALYTIC_TESTCASE_BICUBIC_HERMITE_CONVERGENCE
  
  !
  !================================================================================================================================
  !   
  
  !>Check if the convergence of linear langrange interpolation is expected.
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CONVERGENCE(NUMBER_OF_ELEMENTS_XI_START,NUMBER_OF_ELEMENTS_XI_END, &
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
      
      CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC(i,i,0,INTERPOLATION_SPECIFICATIONS,FIELD)
      CALL CMISSAnalyticAnalysisNodeAbsoluteErrorGet(FIELD,1,(i+1)**2/2+1,1,1,VALUE,Err)

      Y_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(VALUE)
      X_VALUES((i-NUMBER_OF_ELEMENTS_XI_START)/NUMBER_OF_ELEMENTS_XI_INTERVAL+1)=log10(HEIGHT/i)
      CALL ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(1,1,1,1,1)
   
    ENDDO
  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CONVERGENCE
  
  
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
    
  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC(NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS, &
    & INTERPOLATION_SPECIFICATIONS,DependentField)
    !Argument variables 
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_X_ELEMENTS !<number of elements on x axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Y_ELEMENTS !<number of elements on y axis
    INTEGER(CMISSIntg), INTENT(IN) :: NUMBER_GLOBAL_Z_ELEMENTS !<number of elements on z axis
    INTEGER(CMISSIntg), INTENT(IN) :: INTERPOLATION_SPECIFICATIONS !<the interpolation specifications
    TYPE(CMISSFieldType) :: DependentField
    !Local Variables
    INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
    INTEGER(CMISSIntg) :: MPI_IERROR
    INTEGER(CMISSIntg) :: ANALYTIC_FUNCTION
    INTEGER(CMISSIntg) :: Interpolation(3)
    INTEGER(CMISSIntg) :: NumberOfElements(3)
    INTEGER(CMISSIntg) :: NumberOfGaussPoints(3)
    INTEGER(CMISSIntg) :: NumberOfXi
    INTEGER(CMISSIntg) :: NumberOfFieldVariables,NumberOfFieldComponents
    INTEGER(CMISSIntg) :: MeshComponentIdx,FieldComponentIdx,xi,SolverIdx,EquationsetIdx

    REAL(CMISSDP) :: MaterialParameters(6)
    REAL(CMISSDP) :: MeshDimensions(3)

    LOGICAL :: EXPORT_FIELD

    TYPE(CMISSBasisType) :: Basis
    TYPE(CMISSCoordinateSystemType) :: CoordinateSystem
    TYPE(CMISSGeneratedMeshType) :: GeneratedMesh
    TYPE(CMISSMeshType) :: Mesh
    TYPE(CMISSDecompositionType) :: Decomposition
    TYPE(CMISSEquationsType) :: Equations
    TYPE(CMISSEquationsSetType) :: EquationsSet
    TYPE(CMISSFieldType) :: AnalyticField,GeometricField,MaterialField
    TYPE(CMISSFieldsType) :: Fields
    TYPE(CMISSProblemType) :: Problem
    TYPE(CMISSRegionType) :: Region
    TYPE(CMISSSolverType) :: Solver
    TYPE(CMISSSolverEquationsType) :: SolverEquations
    
    NUMBER_OF_DOMAINS=1
    IF((NUMBER_GLOBAL_Y_ELEMENTS == 0) .AND. (NUMBER_GLOBAL_Z_ELEMENTS == 0)) THEN
      NumberOfXi = 1
    ELSEIF (NUMBER_GLOBAL_Z_ELEMENTS == 0) THEN
      NumberOfXi = 2
    ELSE
      NumberOfXi = 3
    ENDIF
    Interpolation = (/INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS,INTERPOLATION_SPECIFICATIONS/)
    NumberOfElements = (/NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS/)
    MeshDimensions = (/LENGTH,WIDTH,HEIGHT/)
    NumberOfGaussPoints = (/4,4,4/)

    !Broadcast the number of elements in the X & Y directions and the number of partitions to the other computational nodes
    CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_World,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_World,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_World,MPI_IERROR)
    CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_World,MPI_IERROR)
    CALL MPI_BCAST(INTERPOLATION_SPECIFICATIONS,1,MPI_INTEGER,0,MPI_COMM_World,MPI_IERROR)

    !=CREATE COORDINATE SYSTEM=====================================================================================================
    !Start the creation of a new RC coordinate system
    CALL CMISSCoordinateSystemTypeInitialise(CoordinateSystem,Err)
    CALL CMISSCoordinateSystemCreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
    CALL CMISSCoordinateSystemTypeSet(CoordinateSystem,CMISSCoordinateRectangularCartesianType,Err)
    CALL CMISSCoordinateSystemDimensionSet(CoordinateSystem,NumberOfXi,Err)
    CALL CMISSCoordinateSystemOriginSet(CoordinateSystem,ORIGIN,Err)
    CALL CMISSCoordinateSystemCreateFinish(CoordinateSystem,Err)

    !=CREATE Region================================================================================================================
    !Create Region and set CS to newly created 3D RC CS
    CALL CMISSRegionTypeInitialise(Region,Err)
    CALL CMISSRegionCreateStart(RegionUserNumber,World_Region,Region,Err)
    CALL CMISSRegionCoordinateSystemSet(Region,CoordinateSystem,Err)
    CALL CMISSRegionCreateFinish(Region,Err)

    !=CREATE BASIS ================================================================================================================
    CALL CMISSBasisTypeInitialise(Basis,Err)
    CALL CMISSBasisCreateStart(BasisUserNumber,Basis,Err)
    CALL CMISSBasisTypeSet(Basis,CMISSBasisLagrangeHermiteTPType,Err)
    CALL CMISSBasisNumberOfXiSet(Basis,NumberOfXi,Err)
    CALL CMISSBasisInterpolationXiSet(Basis,Interpolation(1:NumberOfXi),Err)
    CALL CMISSBasisQuadratureNumberOfGaussXiSet(Basis,NumberOfGaussPoints(1:NumberOfXi),Err)
    CALL CMISSBasisCreateFinish(Basis,Err)

    !=CREATE Mesh =================================================================================================================
    !Start the creation of a generated Mesh in the Region
    CALL CMISSGeneratedMeshTypeInitialise(GeneratedMesh,Err)
    CALL CMISSGeneratedMeshCreateStart(1,Region,GeneratedMesh,Err)
    CALL CMISSGeneratedMeshTypeSet(GeneratedMesh,1,Err)
    CALL CMISSGeneratedMeshBasisSet(GeneratedMesh,Basis,Err)
    !Define the Mesh on the Region
    CALL CMISSGeneratedMeshOriginSet(GeneratedMesh,ORIGIN(1:NumberOfXi),Err)
    CALL CMISSGeneratedMeshExtentSet(GeneratedMesh,MeshDimensions(1:NumberOfXi),Err)
    CALL CMISSGeneratedMeshNumberOfElementsSet(GeneratedMesh,NumberOfElements(1:NumberOfXi),Err)
    CALL CMISSMeshTypeInitialise(Mesh,Err)
    CALL CMISSGeneratedMeshCreateFinish(GeneratedMesh,1,Mesh,Err)
    
    !=CREATE Decomposition ========================================================================================================
    CALL CMISSDecompositionTypeInitialise(Decomposition,Err)
    CALL CMISSDecompositionCreateStart(1,Mesh,Decomposition,Err)
    !Set the Decomposition to be a general Decomposition with the specified number of domains
    CALL CMISSDecompositionTypeSet(Decomposition,CMISSDecompositionCalculatedType,Err)
    CALL CMISSDecompositionNumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
    CALL CMISSDecompositionCreateFinish(Decomposition,Err)

    !=CREATE GEOMETRIC FIELD=======================================================================================================
    !Start to create a default (geometric) field on the Region
    NumberOfFieldVariables = 1 !Geometric Field Coordinates
    NumberOfFieldComponents = NumberOfXi
    CALL CMISSFieldTypeInitialise(GeometricField,Err)
    CALL CMISSFieldCreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
    CALL CMISSFieldMeshDecompositionSet(GeometricField,Decomposition,Err)
    CALL CMISSFieldTypeSet(GeometricField,CMISSFieldGeometricType,Err)
    CALL CMISSFieldNumberOfVariablesSet(GeometricField,NumberOfFieldVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(GeometricField,CMISSFieldUVariableType,NumberOfFieldComponents,Err)
    DO xi=1,NumberOfXi
      MeshComponentIdx = 1
      FieldComponentIdx = xi
      CALL CMISSFieldComponentMeshComponentSet(GeometricField,CMISSFieldUVariableType,FieldComponentIdx,MeshComponentIdx,Err)
    ENDDO !xi
    CALL CMISSFieldCreateFinish(GeometricField,Err)

    !Update the geometric field parameters
    CALL CMISSGeneratedMeshGeometricParametersCalculate(GeometricField,GeneratedMesh,Err)

    !=CREATE DEPENDENT FIELD=======================================================================================================
    !Create a dependent field
    NumberOfFieldVariables = 2 !Dependent Field Displacement Coordinates & Force 
    NumberOfFieldComponents = NumberOfXi
    CALL CMISSFieldTypeInitialise(DependentField,Err)
    CALL CMISSFieldCreateStart(DependentFieldUserNumber,Region,DependentField,Err)
    CALL CMISSFieldTypeSet(DependentField,CMISSFieldGeneralType,Err)
    CALL CMISSFieldMeshDecompositionSet(DependentField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(DependentField,GeometricField,Err)
    CALL CMISSFieldDependentTypeSet(DependentField,CMISSFieldDependentType,Err)
    CALL CMISSFieldNumberOfVariablesSet(DependentField,NumberOfFieldVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(DependentField,CMISSFieldUVariableType,NumberOfFieldComponents,Err)
    DO xi=1,NumberOfXi
      MeshComponentIdx = 1
      FieldComponentIdx = xi
      CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldUVariableType,FieldComponentIdx,MeshComponentIdx,Err)
      CALL CMISSFieldComponentMeshComponentSet(DependentField,CMISSFieldDelUDelNVariableType,FieldComponentIdx,MeshComponentIdx,Err)
    ENDDO !xi
    CALL CMISSFieldCreateFinish(DependentField,Err)

    !=CREATE MATERIAL FIELD========================================================================================================
    !!TODO:: Set Material Field Interpolation to constant element based interpolation when field i/o and cmgui allows for this
    !Create a material field for a general 2D isotropic material
    NumberOfFieldVariables = 1
    !-=MODIFIY-HERE-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    SELECT CASE(NumberOfXi)
    CASE(1)
      !Prescribe material properties Area,E1
      NumberOfFieldComponents = 2 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/WIDTH*HEIGHT,10.0E3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    CASE(2)
      !Prescribe material properties h,E1,v12
      NumberOfFieldComponents = 3 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/HEIGHT,10.0E3_CMISSDP,0.3_CMISSDP,0.0_CMISSDP,0.0_CMISSDP,0.0_CMISSDP/)
    CASE(3)
      !Prescribe material properties E1,E2,E3 & v13,v23,v12
      NumberOfFieldComponents = 6 !Young's Modulus & Poisson's Ratio
      MaterialParameters = (/10.0E3_CMISSDP,10.0E3_CMISSDP,10.0E3_CMISSDP,0.3_CMISSDP,0.3_CMISSDP,0.3_CMISSDP/)
    END SELECT
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    CALL CMISSFieldTypeInitialise(MaterialField,Err)
    CALL CMISSFieldCreateStart(MaterialFieldUserNumber,Region,MaterialField,Err)
    CALL CMISSFieldTypeSet(MaterialField,CMISSFieldMaterialType,Err)
    CALL CMISSFieldMeshDecompositionSet(MaterialField,Decomposition,Err)
    CALL CMISSFieldGeometricFieldSet(MaterialField,GeometricField,Err)
    CALL CMISSFieldNumberOfVariablesSet(MaterialField,NumberOfFieldVariables,Err)
    CALL CMISSFieldNumberOfComponentsSet(MaterialField,CMISSFieldUVariableType,NumberOfFieldComponents,Err)
    DO xi=1,NumberOfXi
      MeshComponentIdx = 1
      DO FieldComponentIdx=1,NumberOfFieldComponents
        CALL CMISSFieldComponentMeshComponentSet(MaterialField,CMISSFieldUVariableType,FieldComponentIdx,MeshComponentIdx,Err)
      ENDDO !FieldComponentIdx
    ENDDO !xi
    CALL CMISSFieldCreateFinish(MaterialField,Err)
    DO FieldComponentIdx=1,NumberOfFieldComponents
      CALL CMISSFieldComponentValuesInitialise(MaterialField,CMISSFieldUVariableType,CMISSFieldValuesSetType,FieldComponentIdx, &
        & MaterialParameters(FieldComponentIdx),Err)
    ENDDO !FieldComponentIdx

    !=CREATE Equations SET=========================================================================================================
    !Create a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    CALL CMISSEquationsSetTypeInitialise(EquationsSet,Err)
    CALL CMISSEquationsSetCreateStart(EquationsSetUserNumber,Region,GeometricField,EquationsSet,Err)
    !Set the Equations set to be a Elasticity Class, Linear Elasticity type, no subtype, EquationsSet
    SELECT CASE(NumberOfXi)
    CASE(1)
      CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass,CMISSEquationsSetLinearElasticityType, &
        & CMISSEquationsSetOneDimensionalSubtype,Err)
    CASE(2)
      CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass,CMISSEquationsSetLinearElasticityType, &
        & CMISSEquationsSetPlaneStressSubtype,Err)
    CASE(3)
      CALL CMISSEquationsSetSpecificationSet(EquationsSet,CMISSEquationsSetElasticityClass,CMISSEquationsSetLinearElasticityType, &
        & CMISSEquationsSetThreeDimensionalSubtype,Err)
    END SELECT
    CALL CMISSEquationsSetCreateFinish(EquationsSet,Err)
    !Create the Equations set dependent field variables
    CALL CMISSEquationsSetDependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
    CALL CMISSEquationsSetDependentCreateFinish(EquationsSet,Err)
    !Create the Equations set material field variables
    CALL CMISSEquationsSetMaterialsCreateStart(EquationsSet,MaterialFieldUserNumber,MaterialField,Err)
    CALL CMISSEquationsSetMaterialsCreateFinish(EquationsSet,Err)
    !Create the Equations set analtyic field variables
    SELECT CASE(NumberOfXi)
    CASE(1)
      ANALYTIC_FUNCTION=CMISSEquationsSetLinearElasticityEquationOneDim1
    CASE(2)
      ANALYTIC_FUNCTION=CMISSEquationsSetLinearElasticityEquationTwoDim1
    CASE(3)
      ANALYTIC_FUNCTION=CMISSEquationsSetLinearElasticityEquationThreeDim1
    END SELECT
    CALL CMISSFieldTypeInitialise(AnalyticField,Err)
    CALL CMISSEquationsSetAnalyticCreateStart(EquationsSet,ANALYTIC_FUNCTION,AnalyticFieldUserNumber,AnalyticField,Err)
    CALL CMISSEquationsSetAnalyticCreateFinish(EquationsSet,Err)

    !=CREATE Equations SET EQUATION================================================================================================
    !Create the Equations set Equations
    CALL CMISSEquationsTypeInitialise(Equations,Err)
    CALL CMISSEquationsSetEquationsCreateStart(EquationsSet,Equations,Err)
    CALL CMISSEquationsSparsityTypeSet(Equations,CMISSEquationsFullMatrices,Err)
                                                !CMISSEquationsSparseMatrices=1 !<Use sparse matrices for the Equations.
                                                !CMISSEquationsFullMatrices=2 !<Use fully populated matrices for the Equations. 
    CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
                                              !CMISSEquationsNoOutput !<No output from the Equations.
                                              !CMISSEquationsTimingOutput !<Timing information output.
                                              !CMISSEquationsMatrixOutput !<All below and equation matrices output.
                                              !CMISSEquationsElementMatrixOutput !<All below and element matrices output.
    CALL CMISSEquationsSetEquationsCreateFinish(EquationsSet,Err)

    !Set up the boundary conditions as per the analytic solution
    CALL CMISSEquationsSetBoundaryConditionsAnalytic(EquationsSet,Err)
  
    !=CREATE Problem===============================================================================================================
    !Create the Problem
    CALL CMISSProblemTypeInitialise(Problem,Err)
    CALL CMISSProblemCreateStart(ProblemUserNumber,Problem,Err)
    !Set the Problem to be a elasticity class, linear elasticity type with no subtype.
    CALL CMISSProblemSpecificationSet(Problem,CMISSProblemElasticityClass,CMISSProblemLinearElasticityType, &
      & CMISSProblemNoSubtype,Err)
    CALL CMISSProblemCreateFinish(Problem,Err)

    !=CREATE Problem CONTROL LOOP==================================================================================================
    !Create the Problem control loop
     CALL CMISSProblemControlLoopCreateStart(Problem,Err)
     CALL CMISSProblemControlLoopCreateFinish(Problem,Err)

    !=CREATE Problem Solver========================================================================================================
    !Start the creation of the Problem Solvers
    SolverIdx = 1
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSProblemSolversCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverIdx,Solver,Err)
    CALL CMISSSolverLibraryTypeSet(Solver,CMISSSolverPETScLibrary,Err)
                                         !CMISSSolverCMISSLibrary     !<CMISS (internal) Solver library.
                                         !CMISSSolverPETScLibrary     !<PETSc Solver library.
                                         !CMISSSolverMUMPSLibrary     !<MUMPS Solver library.
                                         !CMISSSolverSuperLULibrary   !<SuperLU Solver library.
                                         !CMISSSolverSpoolesLULibrary !<SPOOLES Solver library.
                                         !CMISSSolverUMFPACKLibrary   !<UMFPACK Solver library.
                                         !CMISSSolverLUSOLLibrary     !<LUSOL Solver library.
                                         !CMISSSolverESSLLibrary      !<ESSL Solver library.
                                         !CMISSSolverLAPACKLibrary    !<LAPACK Solver library.
                                         !CMISSSolverTAOLibrary       !<TAO Solver library.
                                         !CMISSSolverHypreLibrary     !<Hypre Solver library.
    CALL CMISSSolverLinearTypeSet(Solver,CMISSSolverLinearDirectSolveType,Err)
                                        !CMISSSolverLinearDirectSolveType    !<Direct linear Solver type.
                                        !CMISSSolverLinearIterativeSolveType !<Iterative linear Solver type.
    !CALL CMISSSolverLinearDirectTypeSet(Solver,CMISSSolverDirectLU,Err)
                                               !CMISSSolverDirectLU       !<LU direct linear Solver.
                                               !CMISSSolverDirectCholesky !<Cholesky direct linear Solver.
                                               !CMISSSolverDirectSVD      !<SVD direct linear Solver.
    CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
                                        !CMISSSolverNoOutput !<No output from the Solver routines. \see OPENCMISS_SolverOutputTypes,OPENCMISS
                                        !CMISSSolverProgressOutput !<Progress output from Solver routines.
                                        !CMISSSolverTimingOutput !<Timing output from the Solver routines plus below.
                                        !CMISSSolverSolverOutput !<Solver specific output from the Solver routines plus below.
                                        !CMISSSolverSolverMatrixOutput !<Solver matrices output from the Solver routines plus below.
    CALL CMISSProblemSolversCreateFinish(Problem,Err)

    !=CREATE Problem Solver Equations==============================================================================================
    !Create the Problem Solver Equations
    CALL CMISSSolverTypeInitialise(Solver,Err)
    CALL CMISSSolverEquationsTypeInitialise(SolverEquations,Err)
    CALL CMISSProblemSolverEquationsCreateStart(Problem,Err)
    CALL CMISSProblemSolverGet(Problem,CMISSControlLoopNode,SolverIdx,Solver,Err)
    CALL CMISSSolverSolverEquationsGet(Solver,SolverEquations,Err)
    CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
                                                            !CMISSSolverEquationsSparseMatrices !<Use sparse Solver matrices.
                                                            !CMISSSolverEquationsFullMatrices !<Use fully populated Solver matrices.
    EquationsetIdx = 1 !Initialize index of the Equations set that has been added 
                          !(Variable is returned from Problem_SolverEquations_EquationsSet_ADD)
    CALL CMISSSolverEquationsEquationsSetAdd(SolverEquations,EquationsSet,EquationsetIdx,Err)
    CALL CMISSProblemSolverEquationsCreateFinish(Problem,Err)

    !=SOLVE Problem================================================================================================================
    !Solve the Problem
    CALL CMISSProblemSolve(Problem,Err)

    !=OUTPUT SOLUTION===============================================================================================================
    !!TODO:: Output reaction forces in ipnode files
    EXPORT_FIELD=.TRUE.
    IF(EXPORT_FIELD) THEN
      CALL CMISSFieldsTypeInitialise(Fields,Err)
      CALL CMISSFieldsTypeCreate(Region,Fields,Err)
      CALL CMISSFieldIONodesExport(Fields,"LinearElasticityAnalyticExample","FORTRAN",Err)
      CALL CMISSFieldIOElementsExport(Fields,"LinearElasticityAnalyticExample","FORTRAN",Err)
      CALL CMISSFieldsTypeFinalise(Fields,Err)
    ENDIF



  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC

  !
  !================================================================================================================================
  !  

  SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN(CoordinateSystemUserNumber,RegionUserNumber,BasisUserNumber, &
    & GeneratedMeshUserNumber,ProblemUserNumber)

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

  END SUBROUTINE ANALYTIC_LINEAR_ELASTICITY_GENERIC_CLEAN



END PROGRAM ANALYTIC_LINEAR_ELASTICITY_EXAMPLE
